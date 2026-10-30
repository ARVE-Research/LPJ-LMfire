#!/usr/bin/env bash

cat << END > modellist.txt
MRI-CGCM3
MPI-ESM-P
MIROC-ESM
IPSL-CM5A-LR
HadGEM2-ES
GISS-E2-R
CNRM-CM5
MMM
END

while read model
do

  echo $model
  
  sed -e "s/MODEL/$model/g" midHolocene_TEMPLATE.namelist > midHolocene_$model.namelist
  
done < modellist.txt
