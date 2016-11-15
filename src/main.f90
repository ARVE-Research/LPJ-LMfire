program lpj3

use mpi
use drivermod,       only : driver
use mpimod,          only : mproc,initworker,worker

implicit none

integer  :: proc
integer  :: ierr

!-------------------------------------------------

call mpi_init(ierr)

call mpi_comm_rank(mpi_comm_world,proc,ierr)

if (proc == mproc) then  !I am the master node

  call driver()

else

  call initworker()
  call worker()
  
end if

call mpi_finalize(ierr)

end program lpj3
