#!/usr/bin/env bash

module load Autotools pkg-config netCDF-Fortran

./runjob.sh joboptions/CHELSA_spinup.namelist

./runjob.sh joboptions/midHolocene_nohumans.namelist
