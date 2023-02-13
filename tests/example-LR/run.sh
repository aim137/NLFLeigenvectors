#!/bin/bash
#SBATCH --job-name=".shg"
#SBATCH --output=_out_%j
#SBATCH --error=_err_%j
#SBATCH --partition=main
#SBATCH --ntasks=36
#SBATCH --time=12:00:00

bins="/s/martin/SOFTWARE/6BUG_FIXES/yambo-20221026-RTeVecs/bin"
SAVEDIR="../"

mpirun -n 1 $bins/yambo_nl    -F i7-collisions      -J job0c -C out_job0c_collisions -I $SAVEDIR

#generate yambo_nl inputs
./command.sh

for i in 1 2 3 4 5 6 7 8 9 10 11 12; do 
mpirun -n $SLURM_NTASKS $bins/yambo_nl    -F i3-nl-shg_$i        -J job02_$i,job0c -C out_job02_shg_$i -I $SAVEDIR
sleep 3                                 
mpirun -n $SLURM_NTASKS $bins/ypp_nl      -F i4-nl-fourier2      -J job02_$i,job0c -C out_job02_shg_$i -I $SAVEDIR
sleep 3                                 
done
