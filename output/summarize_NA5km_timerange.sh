#!/usr/bin/env bash

# -------------------
# load required modules

module load NCO netCDF netCDF-Fortran

# -------------------
# compile the landcover.f90 program (only really needs to be done once but ok to recompile)

netcdf_include=`nf-config --fflags`

gfortran -Wall -o landcover landcover.f90 $netcdf_include -lnetcdff

# -------------------
# specify the input file as the first argument on the command line

infile=${1}

# specify the start and end years

y0=2005
y1=2024

# list the variable names to summarize, as a comma separated list (no spaces allowed). the variables must be present in the input file!

vars="cover,burnedf,acflux_fire"

# should not need to modify anything below this line

# -------------------
# calculate time average for selected variables

timeavg=${infile%%.*}_avg_$y0"-"$y1.nc

echo "making time average over $y0-$y1"

ncra -O -v $vars -d time,$y0"-"01-01,$y1"-"12-31 $infile $timeavg

# -------------------
# aggregate land cover to tree, grass, and barren fractions

landcover=${timeavg%%.*}_landcover.nc

echo "make landcover" $timeavg $landcover

# work out the size of the output file

xdim=`ncdump -h $infile | grep "x = [0-9]"`
ydim=`ncdump -h $infile | grep "y = [0-9]"`

xlen=`echo $xdim | grep -Po "\d+"`
ylen=`echo $ydim | grep -Po "\d+"`

# create an empty output file based on the bounds above and the template .cdl file

sed -e "s/XLEN/$xlen/g" -e "s/YLEN/$ylen/g" landcover_NA5km.cdl | ncgen -4 -o $landcover

# aggregate PFT cover fractions to simple categories using the fortran program

./landcover $timeavg $landcover
