#!/usr/bin/env bash

cat << END > modellist.txt
AWI-ESM-1-1-LR
INM-CM4-8
MIROC-ES2L
MPI-ESM1-2-LR
END

while read model
do

  echo $model
  
  sed -e "s/MODEL/$model/g" LGM_TEMPLATE.namelist > LGM_$model.namelist
  
done < modellist.txt
