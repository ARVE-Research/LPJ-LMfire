# LPJ-LMfire
The LPJ-LMfire Dynamic Global Vegetation Model

https://doi.org/10.5281/zenodo.1184588

The following software needs to be available in the path when building
- Autotools (including Autoconf, Automake, M4, and libtool)
- pkg-config
- a Fortran compiler
- an MPI implementation that supports the above Fortran compiler
- HDF5 compiled with support for the above MPI (only C library required)
- netCDF compiled against the above HDF5 and MPI, including C and Fortran libraries
