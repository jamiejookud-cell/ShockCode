// 12/23-15: Do a low-Mach number shock similar to LAPD experiment
// NOTE THIS DOES NOT HAVE INJECTION!! - load particles, boost, stay away from boundary!
//
// COUNTERSTREAMING SIMULATIONS
// 12/12/14: Boosting particles in opposite directions depending on which half
// of grid located
// Version of shock.cxx that I was trying to run in 2D in early Nov., 2014
////////////////////////////////////////////////////////////////
//
//   Shock Example - First attempt to setup a relativistic shock.
//
//   Start with initially unmagnetized plasma, with specified temperature
//   and flow speed, interacting with a reflecting, conducting wall.
//
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////

struct edata {
  species_id sp_id;       /* species id */
  double     vth;         /* thermal energy */
  char fname[256];        /* file to save data */
};

// naming convention for the hydro dump files
#define HYDRO_FILE_FORMAT "hydro/T.%d/%s.%d.%d"
#define SPEC_FILE_FORMAT "hydro/T.%d/spectrum-%s.%d.%d"

begin_globals {

  int restart_interval;
  int energies_interval;
  int fields_interval;
  int ehydro_interval;
  int Hhydro_interval;
  int eparticle_interval;
  int Hparticle_interval;
  int quota_check_interval; //  How frequently to check if quote exceeded

  int rtoggle;              // enables save of last two restart dumps for safety
  double quota_sec;         // Run quota in seconds
  double b0;                // B0
  double va0;

  double topology_x;        // domain topology
  double topology_y;
  double topology_z;

  // Variables for new output format
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;

  // particle spectrum and energy band diagnostics
  edata ede;            // electron species information
  edata edi;            // ion species information
  double emax_band;     // maximum energy for energy band diagnostics
  double emin_band;     // minimum energy for energy band diagnostics
  int nbands;           // # of energy bands
  double emax_spect;    // maximum energy for energy spectrum diagnostics
  double emin_spect;    // minimum energy for energy spectrum diagnostics
  int nbins;            // # of energy bins for energy spectrum diagnostics
  int nx_zone;          // # of cells / zone along x for energy spectrum diagnostics
  int ny_zone;          // # of cells / zone along y for energy spectrum diagnostics
  int nz_zone;          // # of cells / zone along z for energy spectrum diagnostics

  // Vadim:  modified restart machinery
  int write_restart ;      // global flag for all to write restart files
  int write_end_restart ;  // global flag for all to write restart files
};

begin_initialization {

  // Use natural PIC units

  double ec   = 1;         // Charge normalization
  double me   = 1;         // Mass normalization
  double c    = 1;         // Speed of light
  double de   = 1;         // Length normalization (electron inertial length)
  double wpe  = 1;         // Time normalized by electron plasma frequency
  double eps0 = 1;         // Permittivity of space

  //  Note - these choices imply time is normalized to wpe

  //  Some basic numerical parameters

  double cfl_req   = 0.9;  // How close to Courant should we try to run
  double wpedt_max = 0.3;  // How big a timestep is allowed if Courant is not too restrictive
  double damp      = 0.0;  // Level of radiation damping
  int rng_seed     = 1;    // Random number seed increment

  // Set Physics parameters

  double mime   = 100.0;   // Ion mass / electron mass
  double gam = 1.02;       // Gamma of plasma flow
  double V1 = 0.045;       // speed of the flow
  double Te = 0.005;       // Electron temperature in units of me*c**2 in plasma rest frame
  double Ti = 0.005;       // Ion temperature in units of me*c**2 in plasma rest frame
  // Lambda_debye=0.316de for T=0.1
  // Relativistic => 0.1de - RESOLVE THIS - spatial grid

  double b0 = 0.05;        //  Initial magnetic field magnitude
  double theta = 80.0;     //  Shock angle

  double pi = 3.1415927;
  double cs = cos(theta/180.0*pi);
  double sn = sin(theta/180.0*pi);

  // Simulation time and quota

  double tau    = 600;            // Maximum simulation time in units of t*wpe
  double quota  = 15.5;           // run quota in hours
  double quota_sec = quota*3600;  // Run quota in seconds

  // Derived qunatities

  double mi = me*mime;        // Ion mass
  double vthe = sqrt(Te/me);  // Electron thermal velocity
  double vthi = sqrt(Ti/mi);  // Ion thermal velocity
  double va0 = b0/sqrt(mi);

  double Vflow = sqrt(1.0-1.0/(gam*gam));   // Flow Velocity
  double wpi   = wpe/sqrt(mime);            // ion plasma frequency
  double di    = c/wpi;                     // ion inertial length

  // Simulation  parameters

  double nppc = 100;        // Average number of macro particle per cell per species

  double Lx  = 400*de;      // size of box in x dimension
  double Ly  = 0.25*de;     // size of box in y dimension
  double Lz  = 200*de;      // size of box in z dimension

  double topology_x = 16;   // Number of domains in x, y, and z
  double topology_y = 1;
  double topology_z = 8;    // Note! top_x*top_y*top_z=total number cores!

  double nx = 3072;         // Number of cells in x-direction
  double ny = 1;            // Number of cells in y-direction
  double nz = 1536;         // Number of cells in z-direction

  double hx = Lx/nx;        // cell size in x
  double hy = Ly/ny;        // cell size in y
  double hz = Lz/nz;        // cell size in z

  // For 1D - set the transverse cell sizes to be same

  hy = hx;
  Ly = hx;

  double n0 = me*eps0*wpe*wpe/(ec*ec);  // Peak electron (ion) density
  double Ne = nppc*nx*ny*nz;            // total macro electrons in box
  double Np = n0*Lx*Ly*Lz;              // total number of physical electrons     
  Ne  = trunc_granular(Ne,nproc());     // Make it divisible by number of processors
  double weight = Np/Ne;
  double qe = -ec*Np/Ne;                // Charge per macro electron
  double qi = ec*Np/Ne;                 // Charge per macro ion       

  // Determine the time step

  double dg = courant_length(Lx,Ly,Lz,nx,ny,nz);  // Courant length
  double dt = cfl_req*dg/c;                       // Courant limited time step
  if( wpe*dt>wpedt_max) dt=wpedt_max/wpe;         // override timestep if plasma frequency limited

  // Intervals for various output and checks

  int restart_interval = 2000000;
  int energies_interval = 200;
  int interval = int(2.0/(wpi*dt));
  int fields_interval = interval;
  int ehydro_interval = interval;
  int Hhydro_interval = interval;
  int eparticle_interval = 1000*interval;
  int Hparticle_interval = 1000*interval;
  int quota_check_interval = 100;

  // particle spectrum and energy band diagnostics
  double emax_band = 120.0;     // maximum energy for energy band diagnostics
  double emin_band = 0.0;       // minimum energy for energy band diagnostics
  int nbands = 6;               // # of energy bands
  double emax_spect = 1.0E5;    // maximum energy for energy spectrum diagnostics
  double emin_spect = 1.0E-5;   // minimum energy for energy spectrum diagnostics
  int nbins = 1000;             // # of energy bins for energy spectrum diagnostics
  int nx_zone = 32;             // # of cells / zone along x for energy spectrum diagnostics
                                // MAKE SURE that (nx / topology_x) is divisible by nx_zone
  int ny_zone = 1;              // # of cells / zone along y for energy spectrum diagnostics
                                // MAKE SURE that (ny / topology_y) is divisible by ny_zone
  int nz_zone = 32;             // # of cells / zone along z for energy spectrum diagnostics
                                // MAKE SURE that (nz / topology_z) is divisible by nz_zone

  ///////////////////////////////////////////////
  // Setup high level simulation parameters
  num_step             = int(tau/(wpi*dt));
  status_interval      = 200;
  sync_shared_interval = status_interval/2;
  clean_div_e_interval = 2*status_interval;
  clean_div_b_interval = 2*status_interval;

  global->restart_interval   = restart_interval;
  global->energies_interval  = energies_interval;
  global->fields_interval    = fields_interval;
  global->ehydro_interval    = ehydro_interval;
  global->Hhydro_interval    = Hhydro_interval;
  global->eparticle_interval = eparticle_interval;
  global->Hparticle_interval = Hparticle_interval;
  global->quota_check_interval     = quota_check_interval;
  global->quota_sec          = quota_sec;

  global->rtoggle            = 0;
  global->b0  = b0;
  global->va0  = va0;

  global->topology_x  = topology_x;
  global->topology_y  = topology_y;
  global->topology_z  = topology_z;

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the grid

  // Setup basic grid parameters
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = c;
  grid->eps0 = eps0;

  // Define periodic the grid

  define_periodic_grid( 0, -0.5*Ly, -0.5*Lz,    // Low corner
                        Lx, 0.5*Ly, 0.5*Lz,     // High corner
                        nx, ny, nz,             // Resolution
                        topology_x, topology_y, topology_z); // Topology


 // From grid/partition.c: used to determine which domains are on edge
 // This will return the logical "brick" layout of the MPI domains (ix,iy,iz)

# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {                   \
    int _ix, _iy, _iz;                                                    \
    _ix  = (rank);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(global->topology_x);   /* iy = iy+gpy*iz */            \
    _ix -= _iy*int(global->topology_x);   /* ix = ix */                   \
    _iz  = _iy/int(global->topology_y);   /* iz = iz */                   \
    _iy -= _iz*int(global->topology_y);   /* iy = iy */                   \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \
  } END_PRIMITIVE
  int ix, iy, iz;
  RANK_TO_INDEX( int(rank()), ix, iy, iz );

  //  Over-ride the periodic boundary conditions as desired:

  // Set the boundary at x=0 to be reflecting conductor

  if ( ix==0 ) set_domain_field_bc( BOUNDARY( -1,0,0), pec_fields );
  if ( ix==0 ) set_domain_particle_bc( BOUNDARY( -1,0,0), reflect_particles );

 // Set the boundary at x=Lx to be absorbing for particles and fields

  if ( ix==topology_x-1 ) set_domain_field_bc( BOUNDARY(1,0,0), absorb_fields );
  if ( ix==topology_x-1 ) set_domain_particle_bc( BOUNDARY(1,0,0), absorb_particles );

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup the species

  sim_log("Setting up species. ");

  double electron_sort_interval = 25;
  double ion_sort_interval = 25;
  double nmax = 3.0*Ne/nproc();
  double nmovers = 0.1*nmax;
  double sort_method = 1;   //  0=in place and 1=out of place
  species_t *electron = define_species("electron",-ec, me, nmax, nmovers,
                                       electron_sort_interval, sort_method);
  species_t *ion = define_species("ion", ec, mi, nmax, nmovers,
                                  ion_sort_interval, sort_method);

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Setup materials

  sim_log("Setting up materials. ");

  define_material( "vacuum", 1 );

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "shapes" for how to define them and assign them to regions.
  // Also, space is initially filled with the first material defined.

  ////////////////////////////////////////////////////////////////////////////////////////////
  //  Finalize Field Advance

  sim_log("Finalizing Field Advance");

  define_field_array(NULL, damp);

  ///////////////////////////////////////////////////
  // Log diagnostic information about this simulation

  sim_log( "***********************************************" );
  sim_log("* Topology:                       "<<topology_x<<" "<<topology_y<<" "<<topology_z);
  sim_log ( "mi/me = " << mime );
  sim_log ( "Te/(me*c**2) = " << Te );
  sim_log ( "Ti/(me*c**2) = " << Ti );
  sim_log ( "Gamma = " << gam );
  sim_log ( "Vflow/c = " << Vflow );
  sim_log ( "tau*wpe = " << tau );
  sim_log ( "num_step = " << num_step );
  sim_log ( "Lx/di = " << Lx/di );
  sim_log ( "Lx/de = " << Lx/de );
  sim_log ( "Ly/di = " << Ly/di );
  sim_log ( "Ly/de = " << Ly/de );
  sim_log ( "Lz/di = " << Lz/di );
  sim_log ( "Lz/de = " << Lz/de );
  sim_log ( "nx = " << nx );
  sim_log ( "ny = " << ny );
  sim_log ( "nz = " << nz );
  sim_log ( "damp = " << damp );
  sim_log ( "courant = " << c*dt/dg );
  sim_log ( "nproc = " << nproc ()  );
  sim_log ( "nppc = " << nppc );
  sim_log ( " b0 = " << b0 );
  sim_log ( " va0 = " << va0 );
  sim_log ( " di = " << di );
  sim_log ( " Ne = " << Ne );
  sim_log ( "total # of particles = " << 2*Ne );
  sim_log ( "dt*wpe = " << wpe*dt );
  sim_log ( " energies_interval: " << energies_interval );
  sim_log ( "dx/de = " << Lx/(de*nx) );
  sim_log ( "dy/de = " << Ly/(de*ny) );
  sim_log ( "dz/de = " << Lz/(de*nz) );
  sim_log ( "dx/debye = " << (Lx/nx)/(vthe/wpe)  );
  sim_log ( "vthi/c = " << vthi/c );
  sim_log ( "vthe/c = " << vthe/c );

  // Dump simulation information to file "info"
  if (rank() == 0 ) {
    FILE *fp_info;
    if ( ! (fp_info=fopen("info", "w")) ) ERROR(("Cannot open file."));
    fprintf(fp_info, "           ***** Simulation parameters ***** \n");
    fprintf(fp_info, "       mi/me =        %e\n", mime );
    fprintf(fp_info, "       Te/(me*c**2) =     %e\n", Te );
    fprintf(fp_info, "       Ti/(me*c**2) =     %e\n", Ti );
    fprintf(fp_info, "       Gamma =        %e\n", gam );
    fprintf(fp_info, "       Vflow/c =      %e\n", Vflow );
    fprintf(fp_info, "       tau*wpe =      %e\n", tau );
    fprintf(fp_info, "       num_step =         %i\n", num_step );
    fprintf(fp_info, "       Lx/de =        %e\n", Lx/de );
    fprintf(fp_info, "       Ly/de =        %e\n", Ly/de );
    fprintf(fp_info, "       Lz/de =        %e\n", Lz/de );
    fprintf(fp_info, "       Lx/di =        %e\n", Lx/di );
    fprintf(fp_info, "       Ly/di =        %e\n", Ly/di );
    fprintf(fp_info, "       Lz/di =        %e\n", Lz/di );
    fprintf(fp_info, "       nx =           %e\n", nx );
    fprintf(fp_info, "       ny =           %e\n", ny );
    fprintf(fp_info, "       nz =           %e\n", nz );
    fprintf(fp_info, "       damp =         %e\n", damp );
    fprintf(fp_info, "       courant =      %e\n", c*dt/dg );
    fprintf(fp_info, "       nproc =        %e\n", nproc() );
    fprintf(fp_info, "       nppc =         %e\n", nppc );
    fprintf(fp_info, "       b0 =           %e\n", b0 );
    fprintf(fp_info, "       va0 =          %e\n", va0 );
    fprintf(fp_info, "       di =           %e\n", di );
    fprintf(fp_info, "       Ne =           %e\n", Ne );
    fprintf(fp_info, "       total # of particles = %e\n", 2*Ne );
    fprintf(fp_info, "       dt*wpe =       %e\n", wpe*dt );
    fprintf(fp_info, "       energies_interval:     %i\n", energies_interval);
    fprintf(fp_info, "       dx/de =        %e\n", Lx/(de*nx) );
    fprintf(fp_info, "       dy/de =        %e\n", Ly/(de*ny) );
    fprintf(fp_info, "       dz/de =        %e\n", Lz/(de*nz) );
    fprintf(fp_info, "       dx/debye =         %e\n", (Lx/nx)/(vthe/wpe) );
    fprintf(fp_info, "       vthi/c =       %e\n", vthi/c );
    fprintf(fp_info, "       vthe/c =       %e\n", vthe/c );
    fprintf(fp_info, "       ***************************\n");
    fclose(fp_info);
  }

  // for the fortran translation routine
  // write binary info file

  FileIO fp_info;
  if ( ! (fp_info.open("info.bin", io_write)==ok) ) ERROR(("Cannot open file."));
  fp_info.write(&topology_x, 1 );
  fp_info.write(&topology_y, 1 );
  fp_info.write(&topology_z, 1 );
  fp_info.write(&Lx, 1 );
  fp_info.write(&Ly, 1 );
  fp_info.write(&Lz, 1 );
  fp_info.write(&nx, 1 );
  fp_info.write(&ny, 1 );
  fp_info.write(&nz, 1 );
  fp_info.write(&dt, 1 );
  fp_info.write(&mime, 1 );
  fp_info.write(&vthe, 1 );
  fp_info.write(&vthi, 1 );
  fp_info.close();

  ////////////////////////////
  // Load fields

  sim_log( "Loading fields" );
  set_region_field( everywhere, 0, 0, 0,       // Electric field
                    cs*b0, 0.0 , sn*b0  );     // Magnetic field

  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specied as logical equations (i.e. x>0 && x+y<2)

  // LOAD PARTICLES

  sim_log( "Loading particles" );

  // Do a fast load of the particles

  seed_entropy( rank() );  //Generators desynchronized
  double xmin = grid->x0 , xmax = grid->x0+(grid->dx)*(grid->nx);
  double ymin = grid->y0 , ymax = grid->y0+(grid->dy)*(grid->ny);
  double zmin = grid->z0 , zmax = grid->z0+(grid->dz)*(grid->nz);

  double x, y, z, ux, uy, uz, upa1, upe1, uz1, gu1;
  double fs0, f, fs, u2, gg, ue_range, ui_range;

  repeat ( Ne/nproc() ) {

    // Now choose random position

    x = uniform(rng(0), xmin, xmax);
    y = uniform(rng(0), ymin, ymax);
    z = uniform(rng(0), zmin, zmax);

    //  Choose random momentum in rest frame of flow
    //  NOTE: In y and z directions, no boost, so u=v
    ux = normal(rng(0), 0, vthe);
    uy = normal(rng(0), 0, vthe);
    uz = normal(rng(0), 0, vthe);

    // Now boost particle momentum into flow frame. Note that the particle is boosted
    // the left if in the right half of the grid and to the right if in the left half of the
    // grid in order to create a counter stream.

    // ep = sqrt(1.0 + ux*ux + uy*uy + uz*uz);   // particle energy in rest frame of the flow
    // ux=gam*(ux - Vflow*ep);  // Boost x-component of momentum

    // Boost x component of momentum, direction depending on which half of grid particle located.
    // if (x < 0.5*xmax) {ux = gam*(ux + Vflow*ep);}
    // else {ux = gam*(ux - Vflow*ep);}

    // if (x < 0.5*Lx) { ux = gam*(ux + Vflow*ep); }
    // else  { ux = gam*(ux - Vflow*ep); }
    // if (x < 0.5*Lx) { ux = ux + V1; }
    // else { ux = ux - V1; }
    ux = ux - V1;


    // Inject electron into the simulation

    inject_particle( electron, x, y, z, ux, uy, uz, weight, 0, 0);

    // Pick a new momentum for the ion - but load at the same position as electron

    ux = normal(rng(0), 0, vthi);
    uy = normal(rng(0), 0, vthi);
    uz = normal(rng(0), 0, vthi);

    // Same boosting scheme but for ions
    // ep = sqrt(1.0 + ux*ux + uy*uy + uz*uz);   // particle energy in rest frame of the flow
    // ux=gam*(ux - Vflow*ep);  // Boost x-component of momentum

    // Boost x component of momentum, direction depending on which half of grid particle located.
    // if (x < 0.5*xmax) {ux = gam*(ux + Vflow*ep);}
    // else {ux = gam*(ux - Vflow*ep);}

    // if (x < 0.5*Lx) ux = ux + V1;
    // if (x >= 0.5*Lx) ux = ux - V1;
    ux = ux - V1;

    inject_particle( ion, x, y, z, ux, uy, uz, weight, 0, 0);

  }

  sim_log( "Finished loading particles" );

  /*--------------------------------------------------------------------------
    * New dump definition
    *------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------
   * Set data output format
   *
   * This option allows the user to specify the data format for an output
   * dump.  Legal settings are 'band' and 'band_interleave'.  Band-interleave
   * format is the native storage format for data in VPIC.  For field data,
   * this looks something like:
   *
   *   ex0 ey0 ez0 div_e_err0 cbx0 ... ex1 ey1 ez1 div_e_err1 cbx1 ...
   *
   * Banded data format stores all data of a particular state variable as a
   * contiguous array, and is easier for ParaView to process efficiently.
   * Banded data looks like:
   *
   *   ex0 ex1 ex2 ... exN ey0 ey1 ey2 ...
   *
   *------------------------------------------------------------------------*/

  global->fdParams.format = band;

  sim_log ( "Fields output format = band" );

  global->hedParams.format = band;

  sim_log ( "Electron species output format = band" );

  global->hHdParams.format = band;

  sim_log ( "Ion species output format = band" );

  /*--------------------------------------------------------------------------
   * Set stride
   *
   * This option allows data down-sampling at output.  Data are down-sampled
   * in each dimension by the stride specified for that dimension.  For
   * example, to down-sample the x-dimension of the field data by a factor
   * of 2, i.e., half as many data will be output, select:
   *
   *   global->fdParams.stride_x = 2;
   *
   * The following 2-D example shows down-sampling of a 7x7 grid (nx = 7,
   * ny = 7.  With ghost-cell padding the actual extents of the grid are 9x9.
   * Setting the strides in x and y to equal 2 results in an output grid of
   * nx = 4, ny = 4, with actual extents 6x6.
   *
   * G G G G G G G G G
   * G X X X X X X X G
   * G X X X X X X X G         G G G G G G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G   ==>   G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G G G G G G
   * G G G G G G G G G
   *
   * Note that grid extents in each dimension must be evenly divisible by
   * the stride for that dimension:
   *
   *   nx = 150;
   *   global->fdParams.stride_x = 10; // legal -> 150/10 = 15
   *
   *   global->fdParams.stride_x = 8; // illegal!!! -> 150/8 = 18.75
   *------------------------------------------------------------------------*/

  // relative path to fields data from global header
  sprintf(global->fdParams.baseDir, "fields");

  // base file name for fields output
  sprintf(global->fdParams.baseFileName, "fields");

  global->fdParams.stride_x = 1;
  global->fdParams.stride_y = 1;
  global->fdParams.stride_z = 1;

  // add field parameters to list
  global->outputParams.push_back(&global->fdParams);

  sim_log ( "Fields x-stride " << global->fdParams.stride_x );
  sim_log ( "Fields y-stride " << global->fdParams.stride_y );
  sim_log ( "Fields z-stride " << global->fdParams.stride_z );

  // relative path to electron species data from global header
  sprintf(global->hedParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(global->hedParams.baseFileName, "ehydro");

  global->hedParams.stride_x = 1;
  global->hedParams.stride_y = 1;
  global->hedParams.stride_z = 1;

  // add electron species parameters to list
  global->outputParams.push_back(&global->hedParams);

  sim_log ( "Electron species x-stride " << global->hedParams.stride_x );
  sim_log ( "Electron species y-stride " << global->hedParams.stride_y );
  sim_log ( "Electron species z-stride " << global->hedParams.stride_z );

  // relative path to electron species data from global header
  sprintf(global->hHdParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(global->hHdParams.baseFileName, "Hhydro");

  global->hHdParams.stride_x = 1;
  global->hHdParams.stride_y = 1;
  global->hHdParams.stride_z = 1;

  sim_log ( "Ion species x-stride " << global->hHdParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hHdParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hHdParams.stride_z );

  // add electron species parameters to list
  global->outputParams.push_back(&global->hHdParams);

  /*--------------------------------------------------------------------------
   * Set output fields
   *
   * It is now possible to select which state-variables are output on a
   * per-dump basis.  Variables are selected by passing an or-list of
   * state-variables by name.  For example, to only output the x-component
   * of the electric field and the y-component of the magnetic field, the
   * user would call output_variables like:
   *
   *   global->fdParams.output_variables( ex | cby );
   *
   * NOTE: OUTPUT VARIABLES ARE ONLY USED FOR THE BANDED FORMAT.  IF THE
   * FORMAT IS BAND-INTERLEAVE, ALL VARIABLES ARE OUTPUT AND CALLS TO
   * 'output_variables' WILL HAVE NO EFFECT.
   *
   * ALSO: DEFAULT OUTPUT IS NONE!  THIS IS DUE TO THE WAY THAT VPIC
   * HANDLES GLOBAL VARIABLES IN THE INPUT DECK AND IS UNAVOIDABLE.
   *
   * For convenience, the output variable 'all' is defined:
   *
   *   global->fdParams.output_variables( all );
   *------------------------------------------------------------------------*/
  /* CUT AND PASTE AS A STARTING POINT
   * REMEMBER TO ADD APPROPRIATE GLOBAL DUMPPARAMETERS VARIABLE

  output_variables( all );

  output_variables( electric | div_e_err | magnetic | div_b_err |
                    tca      | rhob      | current  | rhof |
                    emat     | nmat      | fmat     | cmat );

  output_variables( current_density  | charge_density |
                    momentum_density | ke_density     | stress_tensor );
  */

  //  These have just the most useful things turned on
  //   global->fdParams.output_variables( electric | magnetic | current | div_e_err );
  //    global->hedParams.output_variables( current_density | charge_density | ke_density | stress_tensor);
  //    global->hHdParams.output_variables( current_density | charge_density | ke_density | stress_tensor);


  // Here we are dumping everthing

  global->fdParams.output_variables( all );
  global->hedParams.output_variables( all );
  global->hHdParams.output_variables( all );

  /*--------------------------------------------------------------------------
   * Convenience functions for simlog output
   *------------------------------------------------------------------------*/

  char varlist[512];
  create_field_list(varlist, global->fdParams);

  sim_log ( "Fields variable list: " << varlist );

  create_hydro_list(varlist, global->hedParams);

  sim_log ( "Electron species variable list: " << varlist );

  create_hydro_list(varlist, global->hHdParams);

  sim_log ( "Ion species variable list: " << varlist );

  sim_log("*** Finished with user-specified initialization ***");

  global->ede.sp_id = electron->id;
  global->ede.vth = sqrt(2.0)*vthe;
  sprintf(global->ede.fname, global->hedParams.baseFileName);

  global->edi.sp_id = ion->id;
  global->edi.vth = sqrt(2.0)*vthi;
  sprintf(global->edi.fname, global->hHdParams.baseFileName);

  global->emax_band  = emax_band;
  global->emin_band  = emin_band;
  global->nbands     = nbands;
  global->emax_spect = emax_spect;
  global->emin_spect = emin_spect;
  global->nbins      = nbins;
  global->nx_zone    = nx_zone;
  global->ny_zone    = ny_zone;
  global->nz_zone    = nz_zone;

  // Upon completion of the initialization, the following occurs:
  // - The synchronization error (tang E, norm B) is computed between domains
  //   and tang E / norm B are synchronized by averaging where discrepancies
  //   are encountered.
  // - The initial divergence error of the magnetic field is computed and
  //   one pass of cleaning is done (for good measure)
  // - The bound charge density necessary to give the simulation an initially
  //   clean divergence e is computed.
  // - The particle momentum is uncentered from u_0 to u_{-1/2}
  // - The user diagnostics are called on the initial state
  // - The physics loop is started
  //
  // The physics loop consists of:
  // - Advance particles from x_0,u_{-1/2} to x_1,u_{1/2}
  // - User particle injection at x_{1-age}, u_{1/2} (use inject_particles)
  // - User current injection (adjust field(x,y,z).jfx, jfy, jfz)
  // - Advance B from B_0 to B_{1/2}
  // - Advance E from E_0 to E_1
  // - User field injection to E_1 (adjust field(x,y,z).ex,ey,ez,cbx,cby,cbz)
  // - Advance B from B_{1/2} to B_1
  // - (periodically) Divergence clean electric field
  // - (periodically) Divergence clean magnetic field
  // - (periodically) Synchronize shared tang e and norm b
  // - Increment the time step
  // - Call user diagnostics
  // - (periodically) Print a status message

} //begin_initialization

#define should_dump(x) \
    (global->x##_interval>0 && remainder(step(), global->x##_interval) == 0)

begin_diagnostics {

  if ( (step()%100)==0 ) sim_log( "Time step: " << step());

  /*--------------------------------------------------------------------------
   * Normal rundata dump
   *------------------------------------------------------------------------*/
  if(step()==0) {
      dump_mkdir("fields");
      dump_mkdir("hydro");
      dump_mkdir("rundata");
      dump_mkdir("data");
      dump_mkdir("restart0");
      dump_mkdir("restart1");  // 1st backup
      dump_mkdir("restart2");  // 2nd backup
      dump_mkdir("particles");

      dump_grid("rundata/grid");
      dump_materials("rundata/materials");
      dump_species("rundata/species");
      global_header("global", global->outputParams);
  } // if

  /*--------------------------------------------------------------------------
   * Normal rundata energies dump
   *------------------------------------------------------------------------*/
  if(should_dump(energies)) {
    dump_energies("rundata/energies", step() == 0 ? 0 : 1);
  } // if

  /*--------------------------------------------------------------------------
   * Field data output
   *------------------------------------------------------------------------*/

  if(step() == 1 || should_dump(fields)) field_dump(global->fdParams);

  /*--------------------------------------------------------------------------
   * Electron species output
   *------------------------------------------------------------------------*/

  if(should_dump(ehydro)) hydro_dump("electron", global->hedParams);

  /*--------------------------------------------------------------------------
   * Ion species output
   *------------------------------------------------------------------------*/

  if(should_dump(Hhydro)) hydro_dump("ion", global->hHdParams);

  /*--------------------------------------------------------------------------
   * Restart dump
   *------------------------------------------------------------------------*/

  // Vadim:
  if (step() && !(step()%global->restart_interval))
    global->write_restart = 1; // set restart flag. the actual restart files are written during the next step
  else if (global->write_restart) {
    global->write_restart = 0; // reset restart flag
    double dumpstart = uptime();
    if(!global->rtoggle) {
      global->rtoggle = 1;
      checkpt("restart1/restart", 0);
    }
    else {
      global->rtoggle = 0;
      checkpt("restart2/restart", 0);
    } // if

    //Vadim:
    mp_barrier(); // Just to be safe

    double dumpelapsed = uptime() - dumpstart;
    sim_log("Restart duration "<<dumpelapsed);

    //Vadim
    if (rank()==0) {
      FileIO fp_restart_info;
      if (!(fp_restart_info.open("latest_restart", io_write)==ok) ) ERROR(("Cannot open file."));
      if(!global->rtoggle) {
        fp_restart_info.print("restart restart2/restart");
      } else {
        fp_restart_info.print("restart restart1/restart");
      }
      fp_restart_info.close();
    }
  } // if

  // Dump particle data
  char subdir[36];
  if ( should_dump(eparticle) && step() !=0 &&
       step() > 0*(global->fields_interval)  ) {
    sprintf(subdir,"particle/T.%d",step());
    dump_mkdir(subdir);
    sprintf(subdir,"particle/T.%d/eparticle",step());
    dump_particles("electron",subdir);
    sprintf(subdir,"particle/T.%d/hparticle",step());
    dump_particles("ion",  subdir);
  }

  #include "energy_local.cxx"

  // Shut down simulation when wall clock time exceeds global->quota_sec.
  // Note that the uptime() is guaranteed to return the same value for all
  // processors (i.e., elapsed time on proc #0), and therefore the abort will
  // be synchronized across processors. Note that this is only checked every
  // few timesteps to eliminate the expensive uptime()all from every
  // timestep. uptime has an ALL_REDUCE in it!
  //
  //Vadim:
  if ((step()>0 && global->quota_check_interval>0 &&
      (step()&global->quota_check_interval)==0) || (global->write_end_restart) ) {
    if ((global->write_end_restart)) {
      global->write_end_restart = 0; // reset restart flag
      sim_log( "Allowed runtime exceeded for this job.  Terminating....\n");
      double dumpstart = uptime();
      checkpt("restart0/restart",0);
      mp_barrier(); // Just to be safe
      sim_log( "Restart dump restart completed." );
      double dumpelapsed = uptime() - dumpstart;
      sim_log("Restart duration "<<dumpelapsed);

      //Vadim:
      if (rank()==0) {
        FileIO fp_restart_info;
        if (!(fp_restart_info.open("latest_restart", io_write)==ok) )
            ERROR(("Cannot open file."));
        fp_restart_info.print("restart restart0/restart");
        fp_restart_info.close();
       }

      /* mp_finalize( grid->mp ); */
      exit(0); // Exit or abort?
    }
    if(uptime() > global->quota_sec) global->write_end_restart = 1;
  }
} // end diagnostics


begin_particle_injection {

  // No particle injection for this simulation
}

begin_current_injection {

  // No current injection for this simulation

}


begin_field_injection {

  // No field injection for this simulation

}


begin_particle_collisions {

  // No particle collisions in this simulation

}
