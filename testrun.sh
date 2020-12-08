#!/bin/bash

wdir=`pwd`

# namelist=$wdir/joboptions/NZ_testrun.namelist
# namelist=$wdir/joboptions/NZ_mid-east_natural.namelist
namelist=$wdir/joboptions/NZ_mid-east_people.namelist

executable=$wdir/src/lpj

# bounds=  #whole NZ domain
# bounds=1995500/3004500/5643500/6087500  #N-half South Island
# bounds=2401500/2404500/5702500/5703500  #8 pixels Canterbury
bounds=2452500/2527500/5684500/5743500  #Banks Peninsula
# bounds=2457500/5725500 #onepix Canterbury
# bounds=2457000/2459000/5724000/5726000 #fourpix Canterbury

pftpar=$wdir/pftpars_NZ.csv   #edited so that all temperate PFTs are evergreen

# outfile=test1pix
# outfile=test4pix
# outfile=test8pix
outfile=banks_pen

mpirun --prefix /usr/local/openmpi -hostfile $wdir/arvehosts.txt -np 102 $executable $namelist $bounds output/$outfile.nc $pftpar > $outfile.log

# mpirun --prefix /usr/local/openmpi -hostfile $wdir/arvehosts.txt -np 2 $executable $namelist $bounds output/$outfile.nc $pftpar > $outfile.log
