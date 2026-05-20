#!/usr/bin/env bash

template=midHolocene_template.namelist

climatedir=/home/terraces/projects/LPJ-LMfire_paleo2025/makeclimate/climate/NA5km

for tmp in `ls $climatedir/midHolocene*.nc`
do

  climate=${tmp##*/}
  
  model=${climate%%.*}
 
  echo $model
  
  sed "s/MODELNAME/$climate/g" midHolocene_template.namelist > $model.namelist
  
done
