#!/usr/bin/env bash

jobfile=${1}

tmp=${jobfile##*/}
job=${tmp%%.*}

echo $job

mpirun --mca pml ob1 --mca btl ^openib -np 31 src/lpj $jobfile -180/180/-90/90 output/$job.nc > output/$job.log
