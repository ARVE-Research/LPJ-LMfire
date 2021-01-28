#!/bin/bash

module purge
module load Autoconf Automake Autotools gettext libtool pkgconfig netCDF-Fortran/4.5.3-gmpich-2021.01

autoreconf -v
