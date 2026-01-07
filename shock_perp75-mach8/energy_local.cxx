// safe allocate
#define ALLOCATE(A,LEN,TYPE)                                             \
  if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) ) ERROR(("Cannot allocate."));

#define LOCAL_CELL_ID(x,y,z) \
  INDEX_FORTRAN_3(x,y,z,0,grid->nx+1,0,grid->ny+1,0,grid->nz+1)

#define SPEC_FILE_FORMAT "hydro/T.%d/spectrum-%s.%d.%d"

{
    // energy band diagnostics
    static int energy_diagnostics_first_time = 1;
    static float *dist;         // array to hold the distribution function
    static std::vector<edata *> edParams; // vector to hold species parameters, such as vth
    static std::size_t total_cells;

    // local energy spectrum diagnostics
    static float dloge, emin_spect_log;
    static float *local_spectrum;           // local particle energy spectrum
    static float *bx_avg, *by_avg, *bz_avg; // local averaged magnetic field
    static int nzones_x;    // number of zones along the x-direction
    static int nzones_y;    // number of zones along the y-direction
    static int nzones_z;    // number of zones along the z-direction
    static int nzones_tot;  // total number of zones in each MPI rank
    static int nbins;       // number of energy bins
    static float incells;   // inverse of total cells in each zone

    const static int nx = grid->nx;
    const static int ny = grid->ny;
    const static int nz = grid->nz;
    const static int nx2 = nx + 2;
    const static int ny2 = ny + 2;
    const static int nz2 = nz + 2;
    const interpolator_t * ALIGNED(16) fi;

    if (should_dump(ehydro)) {
        if (energy_diagnostics_first_time) {
            sim_log("initializing the energy diagnostics");

            total_cells = nx2 * ny2 * nz2;
            ALLOCATE(dist, global->nbands * total_cells, float);

            nzones_x = (grid->nx + global->nx_zone - 1) / global->nx_zone;
            nzones_y = (grid->ny + global->ny_zone - 1) / global->ny_zone;
            nzones_z = (grid->nz + global->nz_zone - 1) / global->nz_zone;
            nzones_tot = nzones_x * nzones_y * nzones_z;
            nbins = global->nbins;
            incells = 1.0 / static_cast<float>(global->nx_zone * global->ny_zone * global->nz_zone);

            ALLOCATE(local_spectrum, (nbins+3) * nzones_tot, float);
            ALLOCATE(bx_avg, nzones_tot, float);
            ALLOCATE(by_avg, nzones_tot, float);
            ALLOCATE(bz_avg, nzones_tot, float);

            emin_spect_log = log10(global->emin_spect);
            dloge = (log10(global->emax_spect) - emin_spect_log) / (nbins - 1);

            edParams.push_back(&global->ede);
            edParams.push_back(&global->edi);

            energy_diagnostics_first_time = 0;
        } // first time called

        // averaged magnetic field in each zone
        CLEAR(bx_avg, nzones_tot);
        CLEAR(by_avg, nzones_tot);
        CLEAR(bz_avg, nzones_tot);
        int zone_x, zone_y, zone_z, zone;
        for (int k = 1; k <= nz; ++k) {
            zone_z = (k - 1) / global->nz_zone;
            for (int j = 1; j <= ny; ++j) {
                zone_y = (j - 1) / global->ny_zone;
                fi = &interpolator(VOXEL(1, j, k, nx, ny, nz));
                for (int i = 1; i <= nx; ++i) {
                    zone_x = (i - 1) / global->nx_zone;
                    zone = zone_x + nzones_x * (zone_y + zone_z * nzones_y);
                    bx_avg[zone] += fi->cbx * incells;
                    by_avg[zone] += fi->cby * incells;
                    bz_avg[zone] += fi->cbz * incells;
                    fi++;
                }
            }
        }

        // Assuming electron and ion are the only two species in species_list
        // Modify if you have more species in the list
        for (int isp = 0; isp < 2; isp++) {   // loop over species
            species_t *sp = find_species_id(edParams.at(isp)->sp_id, species_list);
            sim_log("computing the distribution function for species "<< sp->name);
            double vth = edParams.at(isp)->vth;
            double dke_band = global->emax_band * (vth * vth / 2.0) / global->nbands;

            CLEAR(local_spectrum, (nbins + 3) * nzones_tot);
            for (int k = 0; k < nzones_tot; k++) {
                local_spectrum[k*(nbins+3)  ] = bx_avg[k];
                local_spectrum[k*(nbins+3)+1] = by_avg[k];
                local_spectrum[k*(nbins+3)+2] = bz_avg[k];
            }

            CLEAR(dist, total_cells * global->nbands);
            particle_t *p = sp->p;   // header of the particle array
            int zone_x, zone_y, zone_z, zone;
            int eband;               // energy band
            int ebin;                // energy bin for particle spectrum diagnostics
            for (std::size_t i = 0; i < sp->np; ++i ) { // loop over particles
                double ke = sqrt(1.0 + p->ux*p->ux + p->uy*p->uy + p->uz*p->uz) - 1.0;
                eband = int(ke / dke_band);

                // increment the corresponding bin for cell p->i
                if (eband < global->nbands && eband >= 0)
                    dist[eband * total_cells + p->i]++;

                ebin = floor((log10(ke) - emin_spect_log) / dloge) + 1;
                zone_z = (p->i / (nx2 * ny2) - 1) / global->nz_zone;
                zone_y = ((p->i % (nx2 * ny2)) / nx2 - 1) / global->ny_zone;
                zone_x = (p->i % nx2 - 1) / global->nx_zone;
                zone = zone_x + nzones_x * (zone_y + zone_z * nzones_y);
                if (ebin < nbins && ebin >= 0) {
                    local_spectrum[zone * (nbins + 3) + ebin + 3]++;
                }
                p++; // next particle
            }

#if 0
            // normalize the distribution function
            int ix, iy, iz, ixn, iyn, izn, gcell;

            double np;
            std::size_t icell;
            for (iz = 0; iz < nz2; iz++) {
                for (iy = 0; iy < ny2; iy++) {
                    for (ix = 0; ix < nx2; ix++) {
                        np = 0;  // particle counter

                        icell = LOCAL_CELL_ID(ix, iy, iz);

                        for (int k = 0; k < global->nbands; k++) {
                            np += dist[k * total_cells + icell];  // count the particles in the cell
                        }
                        if (np > 0) {
                            for ( int k = 0; k < global->nbands; k++) {
                                dist[k * total_cells + icell] /= np;  // normalize the distribution
                            }
                        }

                        // is this a ghost cell ?
                        gcell = (ix==0) || (ix==nx2-1) || (iy==0) || (iy==ny2-1) || (iz==0) || (iz==nz2-1);

                        // if yes, assign the value from a neighboring cell
                        if (gcell) {
                            // find the neighboring cell
                            ixn = ix;
                            iyn = iy;
                            izn = iz;

                            if (ix == 0) ixn++;
                            if (ix == nx2-1) ixn--;
                            if (iy == 0) iyn++;
                            if (iy == ny2-1) iyn--;
                            if (iz == 0) izn++;
                            if (iz == nz2-1) izn--;

                            std::size_t nid = LOCAL_CELL_ID (ixn,iyn,izn);

                            for (int k = 0; k < global->nbands; k++) {
                                dist[k * total_cells + icell] = dist[k * total_cells + nid];
                            }
                        } // if (gcell)
                    } // ix
                } // iy
            } // iz
#endif

            // dump all the energy bands
            char fname[256];
            sim_log(" writing the distribution function to file ");
            sprintf(fname, HYDRO_FILE_FORMAT, step(), edParams.at(isp)->fname, step(), (int)rank());
            sim_log("append data to "<<fname);

            FileIO fedump;
            FileIOStatus status = fedump.open(fname, io_append );
            if (status == fail) ERROR(("Could not open file."));

            fedump.write(dist, total_cells * global->nbands);
            fedump.close();

            // dump local energy spectra
            FileIO fespect;
            sim_log(" writing the local spectra to file ");
            sprintf(fname, SPEC_FILE_FORMAT, step(), edParams.at(isp)->fname, step(), (int)rank());
            status = fespect.open(fname, io_write);
            if (status == fail) ERROR(("Could not open spectrum file."));
            fespect.write(local_spectrum, (nbins + 3) * nzones_tot);
            fespect.close();

        } // end species loop

    } // if (should_dump(hydro))

} // end energy diagnostics
