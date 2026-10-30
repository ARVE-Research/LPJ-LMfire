#!/usr/bin/env bash

module load NCO netCDF netCDF-Fortran

infile=${1}

# take the first tavg years of the simulation and use that for the summary

# ---

tavg=1801-1850

timeavg=${infile%%.*}_avg_$tavg.nc

echo "make time average" $tavg

ncra -O -v cover,burnedf,acflux_fire -d time,"1801-01-01 00:00:00","1851-01-01 00:00:00" $infile $timeavg

# ---

landcover=${timeavg%%.*}_landcover.nc

echo "make landcover" $timeavg $landcover

xdim=`ncdump -h $infile | grep "x = [0-9]"`
ydim=`ncdump -h $infile | grep "y = [0-9]"`

xlen=`echo $xdim | grep -Po "\d+"`
ylen=`echo $ydim | grep -Po "\d+"`

sed -e "s/XLEN/$xlen/g" -e "s/YLEN/$ylen/g" landcover_NA5km.cdl | ncgen -4 -o $landcover

./landcover $timeavg $landcover
