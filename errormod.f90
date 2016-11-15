module errormod

!this module should be compiled with the compiler wrapper mpif90

use mpi, only : mpi_status_size

implicit none

public :: netcdf_err
public :: mpi_err

integer :: ncstat
integer :: ierr
integer, dimension(mpi_status_size) :: mpistat

contains

!Internal subroutines for error handline - checks error status after each call,
!prints out text message each time an error code is returned. 

!-------------------------

subroutine netcdf_err(ncstat)

use netcdf, only : nf90_strerror

implicit none

integer, intent(in) :: ncstat

write(0,'(a,i5,a,a)')' NetCDF error ',ncstat,' encountered: ',trim(nf90_strerror(ncstat))
stop

end subroutine netcdf_err

!-------------------------

subroutine mpi_err(ierr)

integer, intent(in) :: ierr

write(0,*)'MPI error',ierr
stop

end subroutine mpi_err

!-------------------------

end module errormod
