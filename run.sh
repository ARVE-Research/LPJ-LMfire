#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=40:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-T106778 # research project to submit under.
#SBATCH --ntasks=16 # specify number of processors.
#SBATCH --mail-type=END
#SBATCH --error=error
#SBATCH --mem=12288
#SBATCH --mail-user=rv237@exeter.ac.uk # email me at job completion

# Commands you wish to run must go here, after the PBS directives

mpirun -np 16  ~/LPJ_ox/src/lpj joboptions/standard_run.namelist -180/180/-90/90 /gpfs/ts0/projects/Research_Project-T106778/Rayanne/LPJLMfire/Outputs/all_full_38.nc
