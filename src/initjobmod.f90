module initjobmod

implicit none

public :: initjob
public :: initlucc

contains

!-----------------------------------------------------------------------------------------------------------------

subroutine initjob(ncells,ntiles,spinupyears,transientyears)

use parametersmod,   only : sp,dp,i8,area
use iovariablesmod,  only : cfile_spinup,cfile_transient,soilfile,                            &
                            dospinup,dotransient,co2file,topofile,topofid,elvid,slopeid,landfid,topotime,          &
                            fixedco2,ocean_uptake,cal_year,nspinyrsout,outputvar,    &
                            lu_turn_yrs,popdfile,poppfile,maxmem,bounds,                     &
                            outputfile,srtx,srty,cntx,cnty,endx,endy,inputlonlen,inputlatlen, &
                            lucc,cellindex,lonvect,latvect,co2vect,nolanduse,nclimv,calcforagers,projgrid,geolon,geolat
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
  popdfile,           &
  poppfile,           &
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
  nolanduse

!-------------------------
!initialize variables with a default value if they are not specified in the namelist

spinupyears    = -9999
transientyears = -9999
nspinyrsout    = -9999
nolanduse      = .false.

write(0,*)'==== Welcome to LPJ, the versatile DGVM ===='

!read the joboptions

call getarg(1,jobfile)

open(10,file=jobfile,status='old')

read(10,nml=joboptions)

close(10)

if (spinupyears <= 0) then
  write(0,*)'no years indicated for spinup!'
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
!   write(0,'(a,a)')'using co2file: ',trim(co2file)
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

call initsoil(cal_year)

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
  
write(0,*) 'Done reading topofile'  

!-------------------------------
!externally prescribed CO2 concentrations

if (co2file /= '') then
  write(0,'(a,a)')'using co2file: ',trim(co2file)
  call getco2(cal_year,transientyears)
else
  co2vect = fixedco2
end if

!FLAGFLAGFLAG

!co2vect = 295.5638

write(0,*)'done with CO2'

!-------------------------------
!externally prescribed land use (ALCC)

if (poppfile /= '' .and. popdfile /= '') then
  write(0,*)'you cannot specify both a pop_p and pop_d file!'
  stop
end if

if (poppfile /= '') then
  calcforagers = .true.
  lucc = .false.
  ntiles = 1
  call initpopp()

else if (popdfile /= '') then
  lucc = .true.
  calcforagers = .false.
  ntiles = 3
  call initlucc()        !open the LUCC file here

else
  calcforagers = .false.
  lucc = .false.
  ntiles = 1
end if

write(0,'(a,i3)')'WARNING number of land use tiles used in this run: ',ntiles,lucc

!-------------------

allocate(cellindex(ncells,3))  !this is a potential maximum. we will check for valid cells below

!-------------------

write(0,'(a,i5,a,i5,a,i5)')'input files have:  ',inputlonlen,' columns and',inputlatlen,' rows'
write(0,'(a,i5,a,i5,a,i8)')'cells to calculate:',cntx,' x',cnty,' =',ncells
write(*,*)   'starting at:       ',srtx,srty,lonvect(srtx),latvect(srty)

write(0,'(a,2f8.4)')'landuse turnover: ',lu_turn_yrs

!-------------------
!allocate the input variable array (all cells for one year) and initialize with some data

allocate(in_master(ncells))

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
      in_master(idx)%cellarea = 1.e8   !10 km grid. NB this should be flexible and handle grids of arbitrary cell size!
     
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
  call initstatevars(in_master(i),sv_master(i),ismaster)
end do

!---

call netcdf_create(ncells)

!---

end subroutine initjob

!-----------------------------------------------------------------------------------------------------------------

subroutine initlucc()

use netcdf
use typesizes
use iovariablesmod, only : lucctime,popdfile,popfid,popdtime
use errormod,       only : ncstat,netcdf_err

implicit none

integer :: dimid
integer :: varid
integer :: tlen

!-------
!open the population density file

ncstat = nf90_open(popdfile,nf90_nowrite,popfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(popfid,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(popfid,dimid,len=tlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(popdtime(tlen))

ncstat = nf90_inq_varid(popfid,'time',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(popfid,varid,popdtime)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
end subroutine initlucc

!-----------------------------------------------------------------------------------------------------------------

subroutine initpopp()

!since you can only have either a popd or a popp file but not both, we recycle variables

use netcdf
use typesizes
use iovariablesmod, only : lucctime,poppfile,popfid,popdtime
use errormod,       only : ncstat,netcdf_err

implicit none

integer :: dimid
integer :: varid
integer :: tlen

!-------
!open the forager potential population density file
  
ncstat = nf90_open(poppfile,nf90_nowrite,popfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(popfid,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(popfid,dimid,len=tlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(popdtime(tlen))

ncstat = nf90_inq_varid(popfid,'time',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(popfid,varid,popdtime)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

end subroutine initpopp

!-----------------------------------------------------------------------------------------------------------------

end module initjobmod
