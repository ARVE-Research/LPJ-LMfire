#!/usr/bin/env bash

module load Autotools pkg-config netCDF-Fortran

jobfile=${1}

tmp=${jobfile##*/}
job=NA5km_${tmp%%.*}

echo "running $job"

mpirun --mca pml ob1 --mca btl ^openib -np 31 src/lpj $jobfile -4350000/3345000/-3885000/3780000 output/$job.nc > output/$job.log
