#!/bin/bash

wdir=`pwd`

# namelist=$wdir/joboptions/NZ_mid-east_natural.namelist
namelist=$wdir/joboptions/NZ_mid-east_people.namelist

executable=$wdir/src/lpj

# bounds=1995500/3004500/5305500/6759500  #whole NZ domain - run with the filtered landf file in order to get just a subset of NZ
bounds=2220000/2640000/5570000/6070000  #sub region of above: NE-half South Island

pftpar=$wdir/pftpars_NZ.csv   #edited so that all temperate PFTs are evergreen

outfile=NE_SIsland_people

mpirun --prefix /usr/local/openmpi -hostfile $wdir/arvehosts.txt -np 102 $executable $namelist $bounds output/$outfile.nc $pftpar > $outfile.log
