program landcover

! gfortran -Wall -o landcover landcover.f90 -I/home/public/easybuild/software/netCDF-Fortran/4.6.1-gompi-2023a/include -lnetcdff

! summarize some LPJ output

use iso_fortran_env
use netcdf

implicit none

integer, parameter :: sp = real32
integer, parameter :: dp = real64
integer, parameter :: i2 = int16

character(200) :: infile
character(200) :: outfile

! character(50) :: ivarname
! character(50) :: ovarname

integer :: status
integer :: ncid
integer :: dimid
integer :: varid
integer :: xlen
integer :: ylen
integer :: npft
integer :: tlen

real(sp), dimension(2) :: lonrange
real(sp), dimension(2) :: latrange

real(sp), allocatable, dimension(:) :: x
real(sp), allocatable, dimension(:) :: y
real(dp), allocatable, dimension(:) :: time

real(sp), allocatable, dimension(:,:,:,:) :: cover
real(sp), allocatable, dimension(:,:,:) :: totalveg
real(sp), allocatable, dimension(:,:,:) :: treecover
real(sp), allocatable, dimension(:,:,:) :: grasscover
real(sp), allocatable, dimension(:,:,:) :: barren

real(sp), parameter :: missing = -9999.

! character(10) :: cstart
! character(10) :: cnyrs

integer :: srtt
integer :: nyrs

! ------------------------------------------------------------------------------------------------------------

! call getarg(3,cstart)
! call getarg(4,cnyrs)
! 
! read(cstart,*)srtt
! read(cnyrs,*)nyrs

srtt = 1
nyrs = 1

! ---------------------------------------

call getarg(1,infile)

status = nf90_open(infile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lon',dimid)
if (status == nf90_ebaddim) status = nf90_inq_dimid(ncid,'x',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lat',dimid)
if (status == nf90_ebaddim) status = nf90_inq_dimid(ncid,'y',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'pft',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=npft)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'time',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=tlen)
if (status /= nf90_noerr) call handle_err(status)

! write(0,*)xlen,ylen,tlen

allocate(x(xlen))
allocate(y(ylen))
allocate(time(nyrs))
allocate(cover(xlen,ylen,npft,nyrs))

! write(0,*)sizeof(cover) * 1.e-6

status = nf90_inq_varid(ncid,'lon',varid)
if (status == nf90_enotvar) status = nf90_inq_varid(ncid,'x',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,x)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'actual_range',lonrange)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'lat',varid)
if (status == nf90_enotvar) status = nf90_inq_varid(ncid,'y',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,y)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'actual_range',latrange)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'time',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,time,start=[srtt],count=[nyrs])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'cover',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,cover,start=[1,1,1,srtt],count=[xlen,ylen,npft,nyrs])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! -------------------------------------

allocate(totalveg(xlen,ylen,tlen))
allocate(treecover(xlen,ylen,tlen))
allocate(grasscover(xlen,ylen,tlen))
allocate(barren(xlen,ylen,tlen))

totalveg   = sum(cover,dim=3,mask=cover /= missing)
treecover  = sum(cover(:,:,1:7,:),dim=3,mask=cover /= missing)
grasscover = sum(cover(:,:,8:9,:),dim=3,mask=cover /= missing)
barren = 1. - totalveg

where (cover(:,:,1,:) == missing)
  totalveg   = missing
  treecover  = missing
  grasscover = missing
  barren = missing
end where

! -------------------------------------

call getarg(2,outfile)

status = nf90_open(outfile,nf90_write,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'lon',varid)
if (status == nf90_enotvar) status = nf90_inq_varid(ncid,'x',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,x)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'actual_range',lonrange)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'lat',varid)
if (status == nf90_enotvar) status = nf90_inq_varid(ncid,'y',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,y)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'actual_range',latrange)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'time',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,time)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'totalveg',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,totalveg)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'treecover',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,treecover)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'grasscover',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,grasscover)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'barren',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,barren)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

!-------------------------------------------------------

contains

subroutine handle_err(status)

!   Internal subroutine - checks error status after each netcdf call,
!   prints out text message each time an error code is returned. 

integer, intent (in) :: status

if(status /= nf90_noerr) then 
  write(0,*)'NetCDF error: ',trim(nf90_strerror(status))
  stop
end if

end subroutine handle_err

!-------------------------------------------------------

end program landcover