#PBS -V # export all environment variables to the batch job.
#PBS -d . # set working directory to .
#PBS -q pq # submit to the parallel queue
#PBS -l walltime=36:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-T106778 # research project to submit under.
#PBS -l procs=10 # specify number of processors.
#PBS -m e -M rv237@exeter.ac.uk # email me at job completion

# Commands you wish to run must go here, after the PBS directives

mpirun -np 10  ~/LPJ-LMfire/src/lpj joboptions/standard_run.namelist -180/180/-90/90 /gpfs/ts0/projects/Research_Project-T106778/Rayanne/LPJLMfire/Outputs/standard.nc
