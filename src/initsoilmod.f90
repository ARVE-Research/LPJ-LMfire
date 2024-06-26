module initsoilmod

use parametersmod, only : stdout,stderr

implicit none

public :: initsoil

contains

!-----------------------------------------------------------------------------------------------

subroutine initsoil(cal_year,layers)

use netcdf
use typesizes

use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : srtx,cntx,srty,cnty,inputlonlen,inputlatlen,lonvect,latvect,bounds, &
                           soilfile,soilfid,soil,topofile,projgrid,geolon,geolat,gridres
use parametersmod,  only : sp,i1

implicit none

integer, intent(in)  :: cal_year
integer, intent(out) :: layers

real(sp), allocatable, dimension(:) :: depth
real(sp), allocatable, dimension(:,:,:) :: sand
real(sp), allocatable, dimension(:,:,:) :: clay

integer :: i,j
integer :: ncid
integer :: varid
integer :: dimid
integer :: tlen
integer :: srt
integer :: ndims

integer, dimension(1) :: pos
integer, dimension(2) :: xpos,ypos

integer, allocatable, dimension(:) :: times

integer(i1) :: flag

character(10) :: projectedgrid

character(3), dimension(2) :: dimname

real(sp) :: scale_factor
real(sp) :: add_offset

real(sp), parameter, dimension(6) :: zpos_std = [ 2.5, 10., 22.5, 45., 80., 150. ]  !standard depths to soil midpoint (from soilgrids1km)

!-------------------------
!generate file names. Regardless of the path name the input files must always have these names.

write(stdout,'(a,a)')' soil file:    ',trim(soilfile)

!-------------------------
!open soil file

ncstat = nf90_open(soilfile,nf90_nowrite,soilfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! --------------------------
! check the dimension names and figure out if input is a projected or lon-lat grid

ncstat = nf90_inquire(soilfid,nDimensions=ndims)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do i = 1,2

  ncstat = nf90_inquire_dimension(soilfid,i,name=dimname(i))
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end do

ncstat = nf90_get_att(soilfid,nf90_global,'ProjectedGrid',projectedgrid)

if (trim(dimname(1)) == 'x') then

  projgrid = .true.
    
  write(stdout,*)'NB: input data is in a projected grid!'

else if (trim(dimname(1)) == 'lon') then

  projgrid = .false.

else

  write(0,*)'error the dimension names',dimname,'are in an unrecognized format'
  stop

end if

! ---

ncstat = nf90_inq_dimid(soilfid,dimname(1),dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(soilfid,dimid,len=inputlonlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(soilfid,dimname(2),dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(soilfid,dimid,len=inputlatlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(lonvect(inputlonlen))
allocate(latvect(inputlatlen))

ncstat = nf90_inq_varid(soilfid,dimname(1),varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(soilfid,varid,lonvect)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
ncstat = nf90_inq_varid(soilfid,dimname(2),varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)  

ncstat = nf90_get_var(soilfid,varid,latvect)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-------------------------
!calculate the starting indices and and number of  pixels to be calculated

gridres(1) = lonvect(2) - lonvect(1)
gridres(2) = latvect(2) - latvect(1)

write(stdout,*)'input grid resolution: ',gridres

pos = minloc(abs(lonvect - bounds(1)))
xpos(1) = pos(1)

pos = minloc(abs(lonvect - bounds(2)))
xpos(2) = pos(1)

pos = minloc(abs(latvect - bounds(3)))
ypos(1) = pos(1)

pos = minloc(abs(latvect - bounds(4)))
ypos(2) = pos(1)

srtx = minval(xpos)
srty = minval(ypos)

if (bounds(1) == bounds(2) .and. bounds(3) == bounds(4)) then  !special case, run just one nearest gridcell even if coords were ambiguous
  
  cntx = 1
  cnty = 1
  
else

  if (lonvect(srtx) < bounds(1)) srtx = srtx + 1
  cntx = 1 + abs(maxval(xpos) - srtx)

  if (latvect(srty) < bounds(3)) srty = srty + 1
  cnty = 1 + abs(maxval(ypos) - srty)

end if

!write(0,*)'init soil'
!write(0,*)bounds
!write(0,*)xpos,srtx
!write(0,*)ypos,srty
!write(0,*)cntx,cnty
!read(*,*)

!call calcpixels(bounds,srtx,srty,cntx,cnty,gridres,lboundlon,uboundlat)

! --------------------------
! if the grid is projected grid, read geodetic longitude and latitude into arrays and fix the cell area

if (projgrid) then
  
  write(stdout,*)'initsoil: NB input data is in a projected grid!'
  
  allocate(geolon(cntx,cnty))
  allocate(geolat(cntx,cnty))

  ncstat = nf90_inq_varid(soilfid,'lon',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_var(soilfid,varid,geolon,start=[srtx,srty],count=[cntx,cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_inq_varid(soilfid,'lat',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_var(soilfid,varid,geolat,start=[srtx,srty],count=[cntx,cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
    
  write(stdout,*)minval(geolon),maxval(geolon)
  write(stdout,*)minval(geolat),maxval(geolat)
  
  ! in_master(idx)%cellarea = 25000.

end if

!--------------------------
!get soil properties from file

write(stdout,*)'reading soil data'

ncstat = nf90_inq_dimid(soilfid,'layer',dimid)
if (ncstat == nf90_ebaddim) ncstat = nf90_inq_dimid(soilfid,'zpos',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(soilfid,dimid,len=layers)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)  

write(stdout,'(a,i5)')'soil layers input: ',layers
  
allocate(depth(layers))

allocate(sand(cntx,cnty,layers))
allocate(clay(cntx,cnty,layers))

ncstat = nf90_inq_varid(soilfid,'zpos',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(soilfid,varid,depth)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!---

ncstat = nf90_inq_varid(soilfid,'sand',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(soilfid,varid,sand,start=[srtx,srty,1],count=[cntx,cnty,layers])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(soilfid,varid,'scale_factor',scale_factor)
if (ncstat == -43) then  ! attribute not present
  scale_factor = 1.
else if (ncstat /= nf90_noerr) then
  call netcdf_err(ncstat)
end if

ncstat = nf90_get_att(soilfid,varid,'add_offset',add_offset)
if (ncstat == -43) then  ! attribute not present
  add_offset = 0.
else if (ncstat /= nf90_noerr) then
  call netcdf_err(ncstat)
end if

sand = sand * scale_factor + add_offset

!---

scale_factor = 1.
add_offset   = 0.

ncstat = nf90_inq_varid(soilfid,'clay',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(soilfid,varid,clay,start=[srtx,srty,1],count=[cntx,cnty,layers])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(soilfid,varid,'scale_factor',scale_factor)
if (ncstat == -43) then  ! attribute not present
  scale_factor = 1.
else if (ncstat /= nf90_noerr) then
  call netcdf_err(ncstat)
end if

ncstat = nf90_get_att(soilfid,varid,'add_offset',add_offset)
if (ncstat == -43) then  ! attribute not present
  add_offset = 0.
else if (ncstat /= nf90_noerr) then
  call netcdf_err(ncstat)
end if

clay = clay * scale_factor + add_offset

!------------

allocate(soil(cntx,cnty))

do j=1,cnty
  do i=1,cntx

    allocate(soil(i,j)%zpos(layers))
    allocate(soil(i,j)%whc(layers))
    allocate(soil(i,j)%cond(layers))
    allocate(soil(i,j)%sand(layers))
    allocate(soil(i,j)%clay(layers))
    allocate(soil(i,j)%orgm(layers))

    soil(i,j)%zpos = depth
    soil(i,j)%sand = sand(i,j,:)
    soil(i,j)%clay = clay(i,j,:)
    soil(i,j)%orgm = 0.

!     soil(i,j)%zpos(1:2) = depth(1)
!     soil(i,j)%zpos(3:5) = depth(2)
! 
!         soil(i,j)%sand(1:2) = sand(i,j,1)
!         soil(i,j)%sand(3:5) = sand(i,j,2)
! 
!         soil(i,j)%clay(1:2) = clay(i,j,1)
!         soil(i,j)%clay(3:5) = clay(i,j,2)

!     write(stdout,*)'DONE READING SOIL'
!     write(stdout,*)soil(i,j)%zpos
!     write(stdout,*)soil(i,j)%sand
!     write(stdout,*)soil(i,j)%clay

    end do
end do

!-------------------------

deallocate(sand)
deallocate(clay)

write(stdout,*)'finished reading soil'

end subroutine initsoil

!----------------------------------------------------------------------

end module initsoilmod
