#!/bin/bash
#
#SBATCH --qos=premium
#SBATCH -N 16
#SBATCH -t 24:00:00
#SBATCH -C knl,quad,cache
#SBATCH -o vpic%j.out
#SBATCH -e vpic%j.err
#SBATCH -J reconnection
#SBATCH -A m3122
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=guofan@lanl.gov
#SBATCH -L SCRATCH,project

##DW persistentdw name=reconnection

##DW stage_in source=/global/cscratch1/sd/xiaocan/3D-Lx150-bg0.2-150ppc-2048KNL-tracking/input destination=$DW_PERSISTENT_STRIPED_reconnection/input type=directory
##DW stage_out source=$DW_PERSISTENT_STRIPED_reconnection/input/particle destination=/global/cscratch1/sd/xiaocan/3D-Lx150-bg0.2-150ppc-2048KNL-tracking/particle type=directory
##DW stage_out source=$DW_PERSISTENT_STRIPED_reconnection/input/spectrum destination=/global/cscratch1/sd/xiaocan/3D-Lx150-bg0.2-150ppc-2048KNL-tracking/spectrum type=directory
##DW stage_out source=$DW_PERSISTENT_STRIPED_reconnection/input/rundata destination=/global/cscratch1/sd/xiaocan/3D-Lx150-bg0.2-150ppc-2048KNL-tracking/rundata type=directory
##DW stage_out source=$DW_PERSISTENT_STRIPED_reconnection/input/restore0 destination=/global/cscratch1/sd/xiaocan/3D-Lx150-bg0.2-150ppc-2048KNL-tracking/restore0 type=directory
##DW stage_out source=$DW_PERSISTENT_STRIPED_reconnection/input/restore1 destination=/global/cscratch1/sd/xiaocan/3D-Lx150-bg0.2-150ppc-2048KNL-tracking/restore1 type=directory
##DW stage_out source=$DW_PERSISTENT_STRIPED_reconnection/input/restore2 destination=/global/cscratch1/sd/xiaocan/3D-Lx150-bg0.2-150ppc-2048KNL-tracking/restore2 type=directory
##DW stage_out source=$DW_PERSISTENT_STRIPED_reconnection/input/info destination=/global/cscratch1/sd/xiaocan/3D-Lx150-bg0.2-150ppc-2048KNL-tracking/info type=file
##DW stage_out source=$DW_PERSISTENT_STRIPED_reconnection/input/info.bin destination=/global/cscratch1/sd/xiaocan/3D-Lx150-bg0.2-150ppc-2048KNL-tracking/info.bin type=file


##### These are shell commands
date
module swap craype-haswell craype-mic-knl
module unload craype-hugepages2M
module load cray-hdf5-parallel
module load dws
module list
# module load lustre-default

#srun mkdir -p /tmp/xiaocanli
#for f in ./* ; do sbcast --compress --force --fanout=8 -t 240 $f /tmp/xiaocanli/$f; done
#wait
#echo "sbcast done!!!"
#export LD_LIBRARY_PATH=/tmp/xiaocanli/:$LD_LIBRARY_PATH
#ldd /tmp/xiaocanli/reconnection.Linux

#sbcast --compress --force --fanout=8 -t 240 ./untar.sh /tmp/untar.sh &  # Untar after stage in.
#sbcast --compress --force --fanout=8 -t 240 ./tar.sh /tmp/tar.sh &
#sbcast --compress --force --fanout=8 -t 480 ./reconnection.Linux /tmp/reconnection.Linux
# sbcast --compress --force --fanout=8 -t 480 /global/homes/x/xiaocan/vpic_test/reconnection.Linux /tmp/reconnection.Linux
# wait

#cd $DW_JOB_STRIPED
##chmod +x ./ff-coll.Linux
##cd $DW_PERSISTENT_STRIPED_FFCOLL
#mkdir KNL64-PROP
#cd KNL64-PROP

mkdir spectrum
mkdir particle
mkdir tracer
mkdir field_hdf5
mkdir hydro_hdf5
lfs setstripe -S 8388608 -c 8 spectrum
lfs setstripe -S 8388608 -c 8 particle
lfs setstripe -S 8388608 -c 8 tracer
lfs setstripe -S 8388608 -c 8 field_hdf5
lfs setstripe -S 8388608 -c 8 hydro_hdf5

# mkdir particle
# lfs setstripe -S 8388608 -c 72 particle
# mkdir spectrum
# lfs setstripe -S 8388608 -c 72 spectrum

# export MPICH_MPIIO_DVS_MAXNODES=276
# export DVS_MAXNODES=276
# export DVS_BLOCKSIZE=8388608

#mkdir $DW_PERSISTENT_STRIPED_reconnection/input
#cd $DW_PERSISTENT_STRIPED_reconnection/input
#cd $SLURM_SUBMIT_DIR 

export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#Quad Cache:
export PMI_MMAP_SYNC_WAIT_TIME=400
# srun -n 131072 -N 2048 -c 4 --cpu_bind=cores /tmp/reconnection.Linux --tpp 4
# srun -n 131072 -N 2048 -c 4 --cpu_bind=cores ./reconnection.Linux --tpp 4
# time srun -n 32768 -N 512 -c 4 --cpu_bind=cores /tmp/reconnection.Linux --tpp 4
# srun -n 131072 -N 2048 -c 4 --cpu_bind=cores /tmp/xiaocanli/reconnection.Linux --tpp 4
srun -n 1024 -N 16 -c 4 --cpu_bind=cores ./shock_perp.Linux --tpp 4
#srun -n 64 -N 1 -c 4 --cpu_bind=cores ./reconnection.Linux --tpp 4 --restore ./restore1/restore.0

date
echo 'Done'
