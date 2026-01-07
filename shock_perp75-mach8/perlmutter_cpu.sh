#!/bin/bash
#
#SBATCH -q regular
#SBATCH -N 128
#SBATCH -t 12:00:00
#SBATCH -o vpic%j.out
#SBATCH -e vpic%j.err
#SBATCH --constraint cpu
#SBATCH -J shock
#SBATCH -A m3737
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=guofan@lanl.gov
##SBATCH -L SCRATCH,project

##### These are shell commands
date
module swap PrgEnv-cray PrgEnv-gnu
module load cmake
module load cray-hdf5-parallel
# module load lustre-default
module list

#RUN_DIR=/global/cscratch1/sd/xiaocan/power_law_index/mime1836/mime1836_sigmaic1024_bg00_small

#cd $RUN_DIR/files_sbcast

#srun mkdir -p /tmp/xiaocanli
#for f in ./* ; do sbcast --compress --force --fanout=8 -t 240 $f /tmp/xiaocanli/$f; done
#wait
#echo "sbcast done!!!"
#export LD_LIBRARY_PATH=/tmp/xiaocanli/:$LD_LIBRARY_PATH
#ldd /tmp/xiaocanli/reconnection.Linux

###sbcast --compress --force --fanout=8 -t 240 ./untar.sh /tmp/untar.sh &  # Untar after stage in.
#sbcast --compress --force --fanout=8 -t 240 ./reconnection.Linux /tmp/reconnection.Linux &
###sbcast --compress --force --fanout=8 -t 240 ./tar.sh /tmp/tar.sh &
#wait


#mkdir spectrum
#mkdir particle
#mkdir tracer
#mkdir fields
#mkdir hydro
# # lfs setstripe -S 16777216 -c 32 field_hdf5
# # lfs setstripe -S 16777216 -c 32 hydro_hdf5
# # lfs setstripe -S 16777216 -c 32 particle
# # lfs setstripe -S 8388608 -c 32 tracer
# # stripe_medium field_hdf5
# # stripe_medium hydro_hdf5
# # stripe_medium particle
# # stripe_medium tracer
#lfs setstripe -S 8388608 -c 24 spectrum
#lfs setstripe -S 8388608 -c 24 particle
#lfs setstripe -S 8388608 -c 24 tracer
#lfs setstripe -S 8388608 -c 24 fields
#lfs setstripe -S 8388608 -c 24 hydro

#Quad Cache:
export PMI_MMAP_SYNC_WAIT_TIME=400
srun -n 16384 -N 128 -c 1 --cpu_bind=cores shock.Linux --tpp 1

date
echo 'Done'
