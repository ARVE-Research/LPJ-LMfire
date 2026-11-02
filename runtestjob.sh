#!/usr/bin/env bash

jobfile=${1}

tmp=${jobfile##*/}
job=${tmp%%.*}_test

echo $job

mpirun --mca pml ob1 --mca btl ^openib -np 31 src/lpj $jobfile -130/-50/24/60 output/$job.nc > output/$job.log
