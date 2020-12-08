module initjobmod

use parametersmod, only : stdout,stderr

implicit none

public :: initjob
public :: inithumans

contains

!-----------------------------------------------------------------------------------------------------------------

subroutine initjob(ncells,ntiles,nlayers,spinupyears,transientyears)

use parametersmod,   only : sp,dp,i8,area
use iovariablesmod,  only : cfile_spinup,cfile_transient,soilfile,                            &
                            dospinup,dotransient,co2file,topofile,topofid,elvid,slopeid,landfid,topotime,          &
                            fixedco2,ocean_uptake,cal_year,nspinyrsout,outputvar,    &
                            lu_turn_yrs,humanfile,maxmem,bounds,                     &
                            outputfile,srtx,srty,cntx,cnty,endx,endy,inputlonlen,inputlatlen, &
                            calchumans,cellindex,lonvect,latvect,co2vect,nolanduse,nclimv,projgrid,geolon,geolat,startyr_foragers
use coordsmod,       only : parsecoords
use initsoilmod,     only : initsoil
use initclimatemod,  only : initclimate
use getyrdatamod,    only : getco2
use mpistatevarsmod, only : in_master,sv_master,initstatevars
use netcdfsetupmod,  only : netcdf_create

use typesizes
use netcdf
use errormod,       only : ncstat,netcdf_err

implicit none

!arguments

integer, intent(out) :: ncells
integer, intent(out) :: ntiles
integer, intent(out) :: nlayers !number of soil layers present in the input data
integer, intent(out) :: spinupyears
integer, intent(out) :: transientyears

!local variables

integer :: dimid
integer :: varid

real(dp) :: minlon,maxlon
real(dp) :: minlat,maxlat

real(dp) :: ilon,ilat

integer  :: x,y
integer  :: i,j
integer  :: a,b
integer  :: tyears
logical  :: newsave

logical  :: ismaster = .true.

integer :: tlen

integer(i8) :: idx


character(45)  :: coords
character(200) :: jobfile

namelist /joboptions/ &
  cfile_spinup,       &
  cfile_transient,    &
  soilfile,           &
  topofile,           &
  co2file,            &
  humanfile,          &
  spinupyears,        &
  transientyears,     &
  dospinup,           &
  dotransient,        &
  fixedco2,           &
  ocean_uptake,       &
  cal_year,           &
  nspinyrsout,        &
  outputvar,          &
  lu_turn_yrs,        &
  maxmem,             &
  nolanduse,          &
  startyr_foragers

!-------------------------
!initialize variables with a default value if they are not specified in the namelist

spinupyears    = -9999
transientyears = 1
nspinyrsout    = -9999
nolanduse      = .false.
startyr_foragers = 1000

write(stdout,*)'==== Welcome to LPJ, the versatile DGVM ===='

!read the joboptions

call getarg(1,jobfile)

open(10,file=jobfile,status='old')

read(10,nml=joboptions)

close(10)

if (spinupyears <= 0) then
  write(stdout,*)'no years indicated for spinup!'
  stop
end if

if (nspinyrsout < 0) then
  nspinyrsout = spinupyears  !if not specified, write out all years of the spinup
end if

!-------------------------
!read the grid boundary coordinates for the run

call getarg(2,coords)

call parsecoords(coords,bounds)

!-------------------------
!get the name of the output file

call getarg(3,outputfile)

!-------------------------
!initialize the run

! !CO2
!
! if (co2file /= '') then
!   write(stdout,'(a,a)')'using co2file: ',trim(co2file)
!   call getco2(cal_year,transientyears)             !externally prescribed CO2 concentrations
! else
!   co2vect = fixedco2
! end if

!gridded input data
!open all of the specified input files and make sure the grids cover the same area
  
!climate

!  status = nf90_open(cfile_spinup,nf90_nowrite,cfid1)
!  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
!  status = nf90_inq_dimid(cfid1,'lon',dimid)
!  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
!  status = nf90_inquire_dimension(cfid1,dimid,len=xlen)
!
!  status = nf90_inq_varid(cfid1,'lon',varid)
!
!  status = nf90_get_var(cfid1,varid,minlon,start=[1])
!
!  status = nf90_get_var(cfid1,varid,maxlon,start=[xlen])
!
!  if (cfile_transient /= '') then
!    status = nf90_open(cfile_transient,nf90_nowrite,cfid2)
!    if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!  end if

!soil
!topography
!slope
!land fraction


!open the soil initial conditions files and allocate the soils input matrix (and lat and lon vect).
!allocates lonvect, latvect and soil%

call initsoil(cal_year,nlayers)

srtx = max(1,srtx)
srty = max(1,srty)

cntx = min(cntx,inputlonlen)
cnty = min(cnty,inputlatlen)

endx = srtx + cntx - 1
endy = srty + cnty - 1

ncells = cntx * cnty

!-------------------------------
!topo file

ncstat = nf90_open(topofile,nf90_nowrite,topofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!check for a time dimension. if there isn't one, set tlen to 1 and don't read the variable

ncstat = nf90_inq_dimid(topofid,'time',dimid)
if (ncstat /= nf90_noerr) then  !there isn't a time dimension
  
  tlen = 1
  
  allocate(topotime(tlen))

else

  ncstat = nf90_inquire_dimension(topofid,dimid,len=tlen)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  allocate(topotime(tlen))

  ncstat = nf90_inq_varid(topofid,'time',varid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_var(topofid,varid,topotime)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end if

ncstat = nf90_inq_varid(topofid,'elv',elvid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(topofid,'slope',slopeid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(topofid,'landf',landfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
write(stdout,*) 'Done reading topofile'  

!-------------------------------
!externally prescribed CO2 concentrations

if (co2file /= '') then
  write(stdout,'(a,a)')'using co2file: ',trim(co2file)
  write(stdout,*)cal_year,transientyears
  call getco2(cal_year,transientyears)
else
  co2vect = fixedco2
end if

!FLAGFLAGFLAG

!co2vect = 295.5638

write(stdout,*)'done with CO2'

!-------------------------------
!externally prescribed land use (ALCC)

if (humanfile /= '') then

  calchumans = .true.
  ntiles     = 3  ! never used, used (crop or pasture), abandoned/recovering

  call inithumans()        !open the humans file here

else  !world without people simulation

  calchumans = .false.
  ntiles     = 1

end if

write(stdout,'(a,i3)')'WARNING number of land use tiles used in this run: ',ntiles,calchumans

!-------------------

allocate(cellindex(ncells,3))  !this is a potential maximum. we will check for valid cells below

!-------------------

write(stdout,'(a,i5,a,i5,a,i5)')'input files have:  ',inputlonlen,' columns and',inputlatlen,' rows'
write(stdout,'(a,i5,a,i5,a,i8)')'cells to calculate:',cntx,' x',cnty,' =',ncells
write(stdout,'(a,2i5,2f12.1)')  'starting at:       ',srtx,srty,lonvect(srtx),latvect(srty)

write(stdout,'(a,2f8.4)')'landuse turnover: ',lu_turn_yrs

!-------------------
!allocate the input variable array (all cells for one year) and initialize with some data

allocate(in_master(ncells))

!allocate the soil layer elements

! write(0,*)'soil layers in_master',nlayers

! do i = 1,ncells
!   allocate(in_master(i)%soil%zpos(nlayers))
!   allocate(in_master(i)%soil%sand(nlayers))
!   allocate(in_master(i)%soil%clay(nlayers))
!   allocate(in_master(i)%soil%orgm(nlayers))
!   allocate(in_master(i)%soil%bulk(nlayers))
! end do

idx = 1

do y = 1,cnty
  b = srty + y - 1
  do x = 1,cntx
    
    a = srtx + x - 1

    in_master(idx)%idx  = idx
    in_master(idx)%xpos = x
    in_master(idx)%ypos = y
    
    !2015-12: new code to handle projected grids
    
    if (projgrid) then

      in_master(idx)%lon  = geolon(x,y)
      in_master(idx)%lat  = geolat(x,y)
      in_master(idx)%cellarea = 1.e8     !1.e8   !10 km grid. NB this should be flexible and handle grids of arbitrary cell size!
     
    else
	
      in_master(idx)%lon  = lonvect(a)
      in_master(idx)%lat  = latvect(b)
      in_master(idx)%cellarea = area(latvect(b),[30.,30.])  !NB this should be changed to allow geographic input grids of arbitrary resolution

    end if

    idx = idx + 1

  end do
end do

!---
!allocate the state variable array (all cells for one year) and initialize

allocate(sv_master(ncells))

do i = 1,ncells
  call initstatevars(in_master(i),sv_master(i),ismaster,nlayers)
end do

!---

call netcdf_create(ncells)

!---

end subroutine initjob

!-----------------------------------------------------------------------------------------------------------------

subroutine inithumans()

use netcdf
use typesizes
use iovariablesmod, only : humanfile,humanfid,popdtime,lu_sf,lu_ao
use errormod,       only : ncstat,netcdf_err

implicit none

integer :: dimid
integer :: varid
integer :: tlen

!-------
!open the consolidate land use file

ncstat = nf90_open(humanfile,nf90_nowrite,humanfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(humanfid,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(humanfid,dimid,len=tlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(popdtime(tlen))

ncstat = nf90_inq_varid(humanfid,'time',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(humanfid,varid,popdtime)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(humanfid,'lu_crop',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(humanfid,varid,'scale_factor',lu_sf)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(humanfid,varid,'add_offset',lu_ao)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
end subroutine inithumans

!-----------------------------------------------------------------------------------------------------------------

end module initjobmod
