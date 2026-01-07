#!/bin/bash
#SBATCH -J shock_perp      # Job name
#SBATCH --output=%x-%j.out       # Name of stdout output file
#SBATCH --error=%x-%j.err     # Name of stderr error file
#SBATCH -q regular          # Queue (partition) name
#SBATCH -N 32               # Total # of nodes
#SBATCH -n 2000             # Total # of mpi tasks
#SBATCH -t 48:00:00      # Run time (hh:mm:ss)
#SBATCH --constraint=knl
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A m2407
#SBATCH --mail-user=xiaohang@lpl.arizona.edu




date
module swap craype-haswell craype-mic-knl
module unload craype-hugepages2M
module load cray-hdf5-parallel
module load dws
module list
RUN_DIR=/global/cscratch1/sd/xhchen/shock_perp/test_run
cd $RUN_DIR
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export PMI_MMAP_SYNC_WAIT_TIME=400
time srun -n 2000 -N 32 -c 4 --cpu_bind=cores ../shock_perp.Linux --tpp 4

date
echo 'Done'
