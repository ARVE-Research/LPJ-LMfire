module initclimatemod

use parametersmod, only : stdout,stderr

implicit none

public :: initclimate
public :: closeclimate

contains

!--------------------------------------------------------------------------------------------------------------------------

subroutine initclimate(climatefile,ncells)

use parametersmod,   only : sp,dp,i4,i8
use netcdf
use typesizes
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : cfid,cntx,cnty,srtx,srty,climatemonths,climateyears,ibuf,varinfo,nclimv, &
                           maxmem,input_sp,input_i2,timebuflen,geolon,geolat,projgrid
implicit none

!parameters

real(dp), parameter :: bmb  = 1048576.d0 !bytes in one Mb

character(4), dimension(nclimv) :: varname

!character(4), dimension(nclimv), parameter :: varname
!arguments

character(120), intent(in) :: climatefile
integer,        intent(in) :: ncells

!local variables

integer :: i,j
integer :: varid
integer :: dimid

integer :: xsize
integer :: ysize

integer :: timelen

integer, dimension(3) :: chunksizes

real(dp) :: memcheck

integer(i8) :: membytes

integer :: tchunk

integer :: mmos  !number of months of input data that can be held in memory

!-------------------------

write(stdout,'(a,a)')' using climate file: ',trim(climatefile)

!-------------------------
!open climate files 

ncstat = nf90_open(climatefile,nf90_nowrite,cfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!write(stdout,*)'cfid:',cfid

!-------------------------
!retrieve dimensions

ncstat = nf90_inq_dimid(cfid,'lon',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(cfid,dimid,len=xsize)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(cfid,'lat',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(cfid,dimid,len=ysize)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(cfid,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(cfid,dimid,len=climatemonths)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

climateyears = climatemonths / 12

write(stdout,'(a,i6,a)')'climate input data contains',climateyears,' years of data'

!-------------------------
!retrieve variable IDs, scale factor and add offset
varname(1) = 'tmp'
varname(2) = 'pre'
varname(3) = 'cld'
varname(4) = 'wet'
varname(5) = 'dtr'
varname(6) = 'wnd'
varname(7) = 'lght'


do i = 1,nclimv

  ncstat = nf90_inq_varid(cfid,varname(i),varinfo(i)%varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_att(cfid,varinfo(i)%varid,'scale_factor',varinfo(i)%scale_factor)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_att(cfid,varinfo(i)%varid,'add_offset',varinfo(i)%add_offset)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  write(stdout,'(a5,2f10.4)')trim(varname(i)),varinfo(i)%scale_factor,varinfo(i)%add_offset

end do

!-------------------------
!allocate the vector input buffer for climate data

timelen = 12               !this version of the LPJ driver works with only one year of climate data at a time.

if (.not.allocated(ibuf)) then


  allocate(ibuf(cntx,cnty))

  do j = 1,cnty
    do i = 1,cntx
      allocate(ibuf(i,j)%temp(timelen))
      allocate(ibuf(i,j)%prec(timelen))
      allocate(ibuf(i,j)%cldp(timelen))
      allocate(ibuf(i,j)%wetd(timelen))
      allocate(ibuf(i,j)%trng(timelen))
      allocate(ibuf(i,j)%temp0(timelen))
      allocate(ibuf(i,j)%lght(timelen))
      allocate(ibuf(i,j)%wind(timelen))
      allocate(ibuf(i,j)%cropfrac(2))
    end do
  end do

end if

!------------------------------
!read in the lon and lat arrays

if (projgrid.EQV..FALSE.) then

  ncstat = nf90_inq_varid(cfid,'lon',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_var(cfid,varid,ibuf(:,1)%lon,start=[srtx],count=[cntx])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inq_varid(cfid,'lat',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_var(cfid,varid,ibuf(1,:)%lat,start=[srty],count=[cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  !copy across
  !would not be necessary with grids that have explicit x and y for the dimension variables (e.g., projected grids).

  do j = 1,cnty
    ibuf(:,j)%lon = ibuf(:,1)%lon
  end do

  do i = 1,cntx
    ibuf(i,:)%lat = ibuf(1,:)%lat
  end do

else

  ibuf%lon = geolon
  ibuf%lat = geolat

end if

!-----
!memory checks and input buffer setup

write(stdout,'(a,f5.1,a)')'memory limit:',maxmem/1024.,' Gb'

!set timebuflen to the lesser of the number of months in the input climate data file and the 
!amount of memory available to hold these data
!if less than the whole dataset will be read in at once, choose a timebuflen that is an even
!multiple of the time dimension chunk size

mmos = 12 * int((bmb * maxmem) / (288 * ncells))  !6 climate variables x 12 months x 4 bytes = 288

if (mmos >= climatemonths) then
  
  timebuflen = climatemonths
  
  write(stdout,'(a,i5,a)')'all input data fits into memory, buflen:',timebuflen/12,' yrs'

else
  
  !get chunk sizes


  ncstat = nf90_inq_varid(cfid,'tmp',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inquire_variable(cfid,varid,chunksizes=chunksizes)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  tchunk = chunksizes(3)

  write(*,*) 'tchunk',tchunk
  write(*,*) 'mmos',mmos

  !calculate an even multiple of the months to hold in memory and the chunksizes
  
  timebuflen = mmos / tchunk * tchunk
  
  write(stdout,'(a,i5,a,i5,a)')'timebuflen:',mmos/tchunk,' chunks,',timebuflen/12,' yrs'
  
end if

write(stdout,'(a,3i8)')'allocate input structures:',cntx,cnty,timebuflen

if (allocated(input_i2)) deallocate(input_i2)
if (allocated(input_sp)) deallocate(input_sp)

allocate(input_i2(cntx,cnty,timebuflen))
allocate(input_sp(cntx,cnty,timebuflen,nclimv))

input_sp = -9999.

membytes = 4_i8 * (ncells * (timebuflen * 6_i8 + 10_i8))  !six transient climate variables, plus about 10 other time-invariant variables

memcheck = real(membytes,dp) / (bmb * 1024.)

write(stdout,'(a,f5.1,a)')      'total memory requirement of input data:',memcheck,' Gb'

end subroutine initclimate

!--------------------------------------------------------------------------------------------------------------------------

subroutine closeclimate()

use netcdf
use typesizes
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : cfid

implicit none

!---

ncstat = nf90_close(cfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---

end subroutine closeclimate

!--------------------------------------------------------------------------------------------------------------------------

end module initclimatemod
