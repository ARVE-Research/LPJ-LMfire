module initjobmod

use parametersmod, only : stdout,stderr

implicit none

public :: initjob
public :: initlucc

contains

! -----------------------------------------------------------------------------------------------------------------

subroutine initjob(ncells,ntiles,nlayers,spinupyears,transientyears)

use parametersmod,   only : sp,dp,i8
use iovariablesmod,  only : cfile_spinup,cfile_transient,soilfile,                                                        &
                            dospinup,dotransient,co2file,topofile,topofid,elvid,slopeid,landfid,topotime,                 &
                            fixedco2,ocean_uptake,cal_year,nspinyrsout,outputvar,                                         &
                            lu_turn_yrs,popdfile,poppfile,maxmem,bounds,                                                  &
                            outputfile,srtx,srty,cntx,cnty,endx,endy,inputlonlen,inputlatlen,                             &
                            lucc,cellindex,lonvect,latvect,co2vect,nolanduse,nclimv,calcforagers,projgrid,geolon,geolat,  &
                            startyr_foragers,pftparsfile,dosoilco2,timeunit_climate,timeunit_basedate
use coordsmod,       only : parsecoords
use initsoilmod,     only : initsoil
use initclimatemod,  only : initclimate
use getyrdatamod,    only : getco2
use mpistatevarsmod, only : in_master,sv_master,initstatevars
use netcdfsetupmod,  only : netcdf_create
use utilitiesmod,    only : tunit2year,area
use calendarmod,     only : timestruct,ymdt2jd,jd2ymdt

use typesizes
use netcdf
use errormod,       only : ncstat,netcdf_err

implicit none

! arguments

integer, intent(out) :: ncells
integer, intent(out) :: ntiles
integer, intent(out) :: nlayers ! number of soil layers present in the input data
integer, intent(out) :: spinupyears
integer, intent(out) :: transientyears

! local variables

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

integer  :: transmons

integer :: cfid1
integer :: cfid2

logical  :: ismaster = .true.

real(dp) :: day0

real(dp), allocatable, dimension(:) :: time

integer :: tlen

integer(i8) :: idx

real(sp) :: xres
real(sp) :: yres

character(45)  :: coords
character(200) :: jobfile

integer, dimension(8) :: ts

! character(8)  :: date
! character(10) :: time
! character(5)  :: zone

integer :: timeunit_baseyr

! type(timestruct) :: basedate
type(timestruct) :: spinstartd


namelist /joboptions/ &
  cfile_spinup,       &
  cfile_transient,    &
  soilfile,           &
  topofile,           &
  co2file,            &
  popdfile,           &
  poppfile,           &
  pftparsfile,        &
  spinupyears,        &
  transientyears,     &
  dospinup,           &
  dotransient,        &
  fixedco2,           &
  ocean_uptake,       &
  nspinyrsout,        &
  outputvar,          &
  lu_turn_yrs,        &
  maxmem,             &
  nolanduse,          &
  startyr_foragers

!   cal_year,           &  removed

! -------------------------

call date_and_time(values=ts)

write(stdout,*)'==== Welcome to LPJ, the versatile DGVM ===='
write(stdout,10)' Timestamp: ',ts(1),'-',ts(2),'-',ts(3),'T',ts(5),':',ts(6),':',ts(7)


write(stderr,*)'==== Welcome to LPJ, the versatile DGVM ===='
write(stderr,10)' Timestamp: ',ts(1),'-',ts(2),'-',ts(3),'T',ts(5),':',ts(6),':',ts(7)

10 format(a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)

! initialize variables with a default value if they are not specified in the namelist

spinupyears    = 10
transientyears = -1
nspinyrsout    = -9999
! nolanduse      = .false.  ! not used
startyr_foragers = 1000
dosoilco2 = .false.

! read the joboptions

call getarg(1,jobfile)

open(10,file=jobfile,status='old')

read(10,nml=joboptions)

close(10)

write(stdout,'(a,a)')' pftpars file: ',pftparsfile

if (spinupyears <= 0) then
  write(stdout,*)'no years indicated for spinup! '
  stop
end if

if (nspinyrsout < 0) then
  nspinyrsout = spinupyears  ! if not specified, write out all years of the spinup
end if

! -------------------------
! read the grid boundary coordinates for the run

call getarg(2,coords)

call parsecoords(coords,bounds)

! -------------------------
! get the name of the output file

call getarg(3,outputfile)

! -------------------------
! initialize the run

! gridded input data
! open all of the specified input files and make sure the grids cover the same area
  
! -------
! climate

ncstat = nf90_open(cfile_spinup,nf90_nowrite,cfid1)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(cfid1,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(cfid1,dimid,len=tlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(time(tlen))

ncstat = nf90_inq_varid(cfid1,'time',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(cfid1,varid,time)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(cfid1,varid,'units',timeunit_climate)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(cfid1)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

timeunit_baseyr = tunit2year(timeunit_climate)

! -------

if (dotransient .and. transientyears < 0) then

  write(stdout,*)'running all years in transient climate file'

  ncstat = nf90_open(cfile_transient,nf90_nowrite,cfid2)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_inq_dimid(cfid2,'time',dimid)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_inquire_dimension(cfid1,dimid,len=transmons)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  ncstat = nf90_close(cfid2)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  transientyears = transmons / 12

  write(stdout,*)'there are ',transientyears,' years of climate in the file'

end if

! -------
! return the Julian day of the basedate

timeunit_basedate = timestruct(timeunit_baseyr,1,1,0,0,0.)

call ymdt2jd(timeunit_basedate)

! calc Julian day of first day of spinup

! day0 = time(1 + tlen - 12 * spinupyears)
day0 = time(1)  ! days since 1950-01-01

spinstartd%jd = timeunit_basedate%jd + day0

call jd2ymdt(spinstartd)

! calculate the year CE and year BP of the first year of the spinup based on the number of years of spinup requested

cal_year = 1950 - spinstartd%y

if (spinstartd%y < 0) cal_year = cal_year - 1  ! adjust if the start year is in BCE time

write(stdout,'(a,i0,a,i0,a)')' spinup starts at ',spinstartd%y,' CE = ',cal_year,' BP'

! write(stdout,*)'transient ends at ',spinupyears + transientyears + spinstartd%y,' CE'

! open the soil initial conditions files and allocate the soils input matrix (and lat and lon vect).
! allocates lonvect, latvect and soil%

call initsoil(cal_year,nlayers)

srtx = max(1,srtx)
srty = max(1,srty)

cntx = min(cntx,inputlonlen)
cnty = min(cnty,inputlatlen)

endx = srtx + cntx - 1
endy = srty + cnty - 1

ncells = cntx * cnty

! -------------------------------
! topo file

write(stdout,*)'Read topofile'  

ncstat = nf90_open(topofile,nf90_nowrite,topofid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! check for a time dimension. if there isn't one, set tlen to 1 and don't read the variable

ncstat = nf90_inq_dimid(topofid,'time',dimid)
if (ncstat /= nf90_noerr) then  ! there isn't a time dimension
  
  tlen = 1
  
  allocate(topotime(tlen))
  
  topotime = 0

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
if (ncstat == nf90_enotvar) ncstat = nf90_inq_varid(topofid,'elev',elvid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(topofid,'slope',slopeid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(topofid,'landf',landfid)
if (ncstat == nf90_enotvar) ncstat = nf90_inq_varid(topofid,'areafrac',landfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
write(stdout,*) 'Done reading topofile'  

! -------------------------------
! externally prescribed CO2 concentrations

if (co2file /= '') then
  write(stdout,'(a,a)')'using co2file: ',trim(co2file)
  write(stdout,*)cal_year,transientyears
  call getco2(cal_year,spinupyears + max(transientyears,0))
else
  co2vect = fixedco2
end if

write(stdout,*)'done with CO2'

! -------------------------------
! externally prescribed land use (ALCC)

if (poppfile /= '' .and. popdfile /= '') then
  write(stdout,*)'you cannot specify both a pop_p and pop_d file! '
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
  call initlucc()        ! open the LUCC file here

else
  calcforagers = .false.
  lucc = .false.
  ntiles = 1
end if

write(stdout,'(a,i3)')'WARNING number of land use tiles used in this run: ',ntiles,lucc

! -------------------

allocate(cellindex(ncells,3))  ! this is a potential maximum. we will check for valid cells below

! -------------------

write(stdout,'(a,i5,a,i5,a,i5)')'input files have:  ',inputlonlen,' columns and',inputlatlen,' rows'
write(stdout,'(a,i5,a,i5,a,i8)')'cells to calculate:',cntx,' x',cnty,' =',ncells
write(stdout,'(a,i0,a,i0,a,f0.4,a,f0.4)')  'starting at:       ',srtx,' ',srty,' ',lonvect(srtx),' ',latvect(srty)

write(stdout,'(a,2f8.4)')'landuse turnover: ',lu_turn_yrs

! -------------------
! allocate the input variable array (all cells for one year) and initialize with some data

allocate(in_master(ncells))

! allocate the soil layer elements

! write(0,*)'soil layers in_master',nlayers

! do i = 1,ncells
!   allocate(in_master(i)%soil%zpos(nlayers))
!   allocate(in_master(i)%soil%sand(nlayers))
!   allocate(in_master(i)%soil%clay(nlayers))
!   allocate(in_master(i)%soil%orgm(nlayers))
!   allocate(in_master(i)%soil%bulk(nlayers))
! end do

idx = 1

if (projgrid) then

  xres = lonvect(2) - lonvect(1)
  yres = latvect(2) - latvect(1)
  
  write(stdout,'(a,f0.1,a)')'projected grid cell area: ',xres * yres * 1.e-6,' km2'

end if

do y = 1,cnty
  b = srty + y - 1
  do x = 1,cntx
    
    a = srtx + x - 1

    in_master(idx)%idx  = idx
    in_master(idx)%xpos = x
    in_master(idx)%ypos = y
    
    ! 2015-12: new code to handle equal-area projected grids
    
    if (projgrid) then
    
      in_master(idx)%lon  = geolon(x,y)
      in_master(idx)%lat  = geolat(x,y)
      in_master(idx)%cellarea = xres * yres
     
    else
  
      in_master(idx)%lon  = lonvect(a)
      in_master(idx)%lat  = latvect(b)
      in_master(idx)%cellarea = area(latvect(b),[30.,30.])  ! NB this should be changed to allow geographic input grids of arbitrary resolution

    end if

    idx = idx + 1

  end do
end do

! ---
! allocate the state variable array (all cells for one year) and initialize

allocate(sv_master(ncells))

do i = 1,ncells
  call initstatevars(in_master(i),sv_master(i),ismaster,nlayers)
end do

! ---

call netcdf_create(ncells)

! ---

end subroutine initjob

! -----------------------------------------------------------------------------------------------------------------

subroutine initlucc()

use netcdf
use typesizes
use iovariablesmod, only : popdfile,popfid,popdtime
use errormod,       only : ncstat,netcdf_err

implicit none

integer :: dimid
integer :: varid
integer :: tlen

! -------
! open the population density file

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

! -----------------------------------------------------------------------------------------------------------------

subroutine initpopp()

! since you can only have either a popd or a popp file but not both, we recycle variables

use netcdf
use typesizes
use iovariablesmod, only : poppfile,popfid,popdtime
use errormod,       only : ncstat,netcdf_err

implicit none

integer :: dimid
integer :: varid
integer :: tlen

! -------
! open the forager potential population density file
  
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

! -----------------------------------------------------------------------------------------------------------------

end module initjobmod
