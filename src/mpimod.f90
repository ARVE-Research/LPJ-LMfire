module mpimod

use mpi
use errormod,        only : mpistat,ierr
use parametersmod,   only : i1
use parametersmod,   only : stdout,stderr
use mpistatevarsmod, only : inputdata,statevars

implicit none

public  :: initmpi
public  :: master
public  :: initworker
public  :: worker
private :: checksum

!module parameters

integer, parameter :: mproc    = 0             !process ID of the master process
integer, parameter :: celltype = MPI_INTEGER1

!module variables

integer :: nprocs
integer :: tsize
integer :: bufsize

integer(i1), allocatable, dimension(:,:) :: iobuf    !the MPI transfer buffer in 1-byte integers

type(inputdata), allocatable, dimension(:) :: in  !data structures that get 'transfer'ed to/from the MPI buffer
type(statevars), allocatable, dimension(:) :: sv

integer :: a  !array indices for input data buffer
integer :: b

contains

!------------------------------------------------------------------------------------------------------------

subroutine initmpi(ncells,ntiles,nlayers)

use iovariablesmod, only : jobinfo,pftparsfile

implicit none

!arguments

integer, intent(in) :: ncells
integer, intent(in) :: ntiles
integer, intent(in) :: nlayers

!parameters

!integer, parameter :: minbufsize = 1000  !probably want to tune this to avoid too much message passing

!local variables

integer :: minbufsize

! integer, dimension(3) :: jobinfo

integer(i1), allocatable, dimension(:) :: bbuf

integer :: nelem1
integer :: nelem2
integer :: nelem

integer :: bs

!-------------------------------

!find out the number of processors

call mpi_comm_size(mpi_comm_world,nprocs,ierr)

if (nprocs < 2) then
  write(stdout,*)'this job requires more than one process, aborting!'
  stop
end if

write(stdout,*)'processes for this job: ',nprocs

!establish the size of the transfer buffer

minbufsize = min(100,max(10,ncells/(nprocs-1)))

bufsize = min(ncells,minbufsize)

write(stdout,*)'MPI buffer size: ',bufsize,' gridcells'

jobinfo%bufsize     = bufsize
jobinfo%ntiles      = ntiles
jobinfo%nlayers     = nlayers
jobinfo%pftparsfile = pftparsfile

bs = sizeof(jobinfo)

allocate(bbuf(bs))

bbuf = transfer(jobinfo,bbuf)

!broadcast it to all the nodes

write(stdout,*)'broadcasting pftparsfile: ',pftparsfile

call mpi_bcast(bbuf,bs,MPI_INTEGER1,mproc,MPI_COMM_WORLD,ierr)

!allocate the input and state variable structures on the master process

allocate(in(bufsize))
allocate(sv(bufsize))

nelem1 = sizeof(in(1))
nelem2 = sizeof(sv(1))

nelem = nelem1 + nelem2  !bytes

a = nelem1
b = a + 1

allocate(iobuf(nelem,bufsize))

tsize = size(iobuf)  !the total number of elements in the array that is passed as a message

write(stdout,'(a,f6.1,a)')' MPI message size: ',real(tsize)/1024.,' kB'

end subroutine initmpi

!------------------------------------------------------------------------------------------------------------

subroutine master(lastyear,ncells,in_master,sv_master)

use iovariablesmod, only : cellmask

implicit none

!arguments

logical, intent(in) :: lastyear
integer, intent(in) :: ncells

type(inputdata), dimension(:), intent(inout) :: in_master
type(statevars), dimension(:), intent(inout) :: sv_master

!local variables

integer :: proc    !the ID of any given process (node)
integer :: tflag

integer :: i,j,k  !counters
integer :: x,y

real, parameter :: missing = -32768.

!---------------------------------------
!loop through the cells of work

i = 1
j = 0

do
  
!   write(stdout,*)'MPIMOD MASTER',in_master(1)%soil%sand
  
  !wait for a message from any processor before sending them work
  !receive the gridcell data structure

  call mpi_recv(iobuf,tsize,celltype,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)

  do k = 1,bufsize
    in(k) = transfer(iobuf(:a,k),in(1))
    sv(k) = transfer(iobuf(b:,k),sv(1))
  end do
  
  !write(stdout,'(a,i5,5f8.2)')'master finished year ',in(1)%year,sv(1)%tile(1)%soil%sand
  !write(stdout,*)'master finished year ',in(1)%year,sv(1)%tile(1)%soil%sand(1),in(1)%climate%temp(1)

  !put the result into the master state variable structure
  
  do k = 1,bufsize

    j = in(k)%idx

    if (j > 0 .and. j <= ncells) then
    
      in_master(j) = in(k)
      sv_master(j) = sv(k)

    end if
    
  end do
  
  !send the next (bufsize) valid cells to the processor that just answered

  proc = mpistat(MPI_SOURCE)
  
  j = 1
  in%idx = -1
  
  do
    x = in_master(i)%xpos
    y = in_master(i)%ypos
    !if (in_master(i)%climate%temp(1) /= missing .and. in_master(i)%soil%sand(1) >= 0.) then
    if (cellmask(x,y)) then
      
      !load this valid pixel into transfer buffer

      in(j) = in_master(i)
      sv(j) = sv_master(i)
      j = j + 1
      
      if (j > bufsize) exit  !the buffer is full, so send it
    end if
    
    i = i + 1
    
    if (i > ncells) exit     !there is no more input data, so send the buffer
 
  end do
  
  do k = 1,bufsize
    iobuf(:a,k) = transfer(in(k),iobuf(:a,1))
    iobuf(b:,k) = transfer(sv(k),iobuf(b:,1))
  end do

 ! write(stdout,*)'master send job',checksum(iobuf)
  
  call mpi_send(iobuf,tsize,celltype,proc,1,MPI_COMM_WORLD,ierr)  !send the next gridcell structure

  i = i + 1
  
  if (i > ncells) exit

end do

!-------------------
!tell the workers there is no more work

do i = 1,nprocs-1  !always subtract the master process that has ID 0 anyway
    
  !collect the last results outstanding from all processors

  call mpi_recv(iobuf,tsize,celltype,i,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)

  do k = 1,bufsize
    in(k) = transfer(iobuf(:a,k),in(1))
    sv(k) = transfer(iobuf(b:,k),sv(1))
  end do
  
  !put the result into the master state variable structure
  
  do k = 1,bufsize

    j = in(k)%idx

    if (j > 0 .and. j <= ncells) then
      in_master(j) = in(k)
      sv_master(j) = sv(k)
    end if

  end do

  !send the termination message
  
  if (lastyear) then
    tflag = 3
  else
    tflag = 2
  end if
  
  proc = mpistat(MPI_SOURCE)

  call mpi_send(iobuf,1,MPI_INTEGER,proc,tflag,MPI_COMM_WORLD,ierr)
  
end do

end subroutine master

!------------------------------------------------------------------------------------------------------------

subroutine initworker()

use radiationmod,    only : initairmass
use parametersmod,   only : sp,bytes_dp,lm_sapl,sm_sapl,rm_sapl,hm_sapl,pftpar,pft, &
                            allom1,allom2,allom3,allom4,sla,latosa,wooddens,reinickerp
!uses pftparameters.f
!uses the common module variables: iobuf, mproc

use pftparametersmod, only : pftparameters
use mpistatevarsmod, only : initstatevars
use iovariablesmod,  only : jobinfo,pftparsfile

implicit none

! integer, dimension(3) :: jobinfo

integer :: ntiles
integer :: nlayers
real(sp), parameter, dimension(3) :: co2 = [ 273., -8., 0. ]

integer :: nelem1
integer :: nelem2
integer :: nelem

integer :: i

logical :: ismaster = .false.

integer(i1), allocatable, dimension(:) :: bbuf

integer :: rank
integer :: bs

!---
!get the size of the transfer buffer and allocate

bs = sizeof(jobinfo)

allocate(bbuf(bs))

call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)

! write(stderr,*)'init worker',rank,'waiting for broadcast'

call mpi_bcast(bbuf,bs,MPI_INTEGER1,mproc,MPI_COMM_WORLD,ierr)

jobinfo = transfer(bbuf,jobinfo)

bufsize     = jobinfo%bufsize
ntiles      = jobinfo%ntiles
nlayers     = jobinfo%nlayers
pftparsfile = jobinfo%pftparsfile

! write(stderr,*)'worker pftparsfile: ',pftparsfile

allocate(in(bufsize))
allocate(sv(bufsize))

do i = 1,bufsize
  call initstatevars(in(i),sv(i),ismaster,nlayers)
end do

nelem1 = sizeof(in(1))
nelem2 = sizeof(sv(1))

nelem = nelem1 + nelem2  !bytes

a = nelem1
b = a + 1

allocate(iobuf(nelem,bufsize))

tsize = size(iobuf)  !the total number of elements in the array that is passed as a message

! write(stdout,'(a,f6.1,a)')'worker MPI message size',real(tsize)/1024.,' kB'

in%idx = 0

!-------------------------------
!initialize the airmass parameters

call initairmass()

!--------------------------
!initialize the pft parameters

call pftparameters(pftparsfile,pftpar,sla,                                     &
                   pft%tree,pft%evergreen,pft%summergreen,                     &
                   pft%raingreen,pft%needle,pft%boreal,                        &
                   lm_sapl,sm_sapl,hm_sapl,rm_sapl,                            &
                   latosa,allom1,allom2,allom3,allom4,wooddens,co2)
                   
end subroutine initworker

!------------------------------------------------------------------------------------------------------------

subroutine worker()

use lpjmod, only : lpjcore

!uses the common module variables: tsize, bufsize, iobuf, mproc, in, sv

implicit none

integer :: i,k
integer :: tflag

!-------------------------------

! write(stdout,*)'worker checking in'

do  !until we receive a message from the master that there is no more work to do
  
  !send result to master (dummy call on first iteration to request work)

  !write(stdout,'(a,i5,5f8.2)')'worker sent year',in(1)%year,sv(1)%tile(1)%soil%sand

  do k = 1,bufsize
    iobuf(:a,k) = transfer(in(k),iobuf(:a,1))
    iobuf(b:,k) = transfer(sv(k),iobuf(b:,1))
  end do

  call mpi_send(iobuf,tsize,celltype,mproc,1,MPI_COMM_WORLD,ierr)
  
  !receive message from master

  call mpi_recv(iobuf,tsize,celltype,mproc,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)

!   write(stdout,*)'worker got something',bufsize

  do k = 1,bufsize
  
!     write(stdout,*)'k:',k

    in(k) = transfer(iobuf(:a,k),in(1))
    sv(k) = transfer(iobuf(b:,k),sv(1))

  end do

  tflag = mpistat(MPI_TAG)
  
  select case (tflag)
  case(3)
    !write(stdout,*)'case 3'
    exit              !all done so quit the subroutine

  case(2)
    !write(stdout,*)'case 2',i,bufsize
    
    in%idx = 0
    cycle             !go back to the top and wait for more work

  case default
    
!     write(stdout,*)'worker',checksum(iobuf)
    
    do i = 1,bufsize  !loop over grid cells in this sub-block

      if (in(i)%idx < 1) cycle

!     write(stdout,*)'worker starting year ',in(1)%year,in(1)%soil%sand,sv(1)%tile(1)%soil%sand   ! '(a,i5,4f8.2)'

      call lpjcore(in(i),sv(i))

!     write(stdout,'(a,i5,5f8.2)')'worker finished year ',in(1)%year,sv(1)%tile(1)%soil%sand

    end do

  end select
  
end do

end subroutine worker

!------------------------------------------------------------------------------------------------------------

integer(i8) function checksum(val)

use parametersmod, only : i1,i8

implicit none

integer(i1), dimension(:,:), intent(in) :: val

!---

checksum = sum(val)

end function checksum

!------------------------------------------------------------------------------------------------------------

end module mpimod
