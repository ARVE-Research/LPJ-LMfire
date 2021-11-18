# makefile for LPJ-LMfire

#---------------------------------------------
# specify the fortran and C compilers here. normally this needs to be the MPI version of the compiler

CC = mpicc
FC = mpifort

# specify the compiler options. choose either optimized or debugging flags, or neither

# optimized flags for ifort

# FCFLAGS  = -xHost

# OR
# debugging flags for ifort

# FCFLAGS = -debug -fpe0 -traceback -g -check -check noarg_temp_created # -warn all

# OR
# debugging flags for gfortran

FCFLAGS  = -ffree-line-length-none -march='native' -funroll-loops -O3#-finit-local-zero -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -g -fbacktrace -Wall -pedantic 
#FCFLAGS  = -ffree-line-length-none -finit-local-zero -fcheck=all -ffpe-trap=invalid,zero -g -fbacktrace -Wall -pedantic	-fbounds-check
#---------------------------------------------
# specify the location of your netcdf install

# netcdf=/alcc/gpfs2/home/geo/georwa/local
netcdf=/home/public/easybuild/software/netCDF-Fortran/4.5.2-gompi-2020a

#---- should not need to change anything below this line ----

# pointers to netcdf directories

NC_LIB=$(netcdf)/lib
NC_INC=$(netcdf)/include

# preprocessor and linker flags

CPPFLAGS = -Difort -I$(NC_INC)
LDFLAGS  = -L$(NC_LIB)

# --- required libraries --- may not always need libcurl, depends on the HDF install

LIBS     = -lnetcdff -lnetcdf -lhdf5 -lhdf5_hl -lcurl

#note: files without the 'mod' suffix are f77 source that eventually need to be converted

COREOBJS = overprint.o       \
       parametersmod.o   \
       randomdistmod.o   \
       weathergenmod.o   \
       orbitmod.o        \
       geohashmod.o      \
       mpistatevarsmod.o \
       iovariablesmod.o  \
       errormod.o        \
       coordsmod.o       \
       initsoilmod.o     \
       initclimatemod.o  \
       getyrdatamod.o    \
       netcdfsetupmod.o  \
       initjobmod.o      \
       radiationmod.o    \
       bioclimmod.o      \
       alccmod.o         \
       hetrespmod.o      \
       snowmod.o         \
       lightmod.o        \
       waterbalancemod.o \
       pedotransfermod.o \
       simplesoilmod.o   \
       nppmod.o          \
       gppmod.o          \
       establishmentmod.o\
       mortalitymod.o    \
       allocationmod.o   \
       turnovermod.o     \
       killplantmod.o    \
       individualmod.o   \
       spitfiremod.o     \
       soiltemperaturemod.o   \
       summerphenology.o      \
       reproduction.o         \
       photosynthesis.o       \
       isotope.o              \
       pftparameters.o        \
       fire.o                 \
       netcdfoutputmod.o \
       landscape_geometrymod.o \
       foragersmod.o \
       lpjmod.o

DRIVER_PARALLEL = mpimod.o     \
                  drivermod.o  \
                  main.o

DRIVER_SERIAL = drivermod_serial.o \
                main_serial.o

.SUFFIXES: .o .f90 .f .mod

%.o : %.c
	$(CC) $(CFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

all::	parallel

parallel: $(COREOBJS) $(DRIVER_PARALLEL)
	$(FC) $(FCFLAGS) -o lpj $(COREOBJS) $(DRIVER_PARALLEL) $(LDFLAGS) $(LIBS)

serial: $(COREOBJS) $(DRIVER_SERIAL)
	$(FC) $(FCFLAGS) -o lpjserial $(COREOBJS) $(DRIVER_SERIAL) $(LDFLAGS) $(LIBS)

clean::	
	-rm lpj *.mod *.o
