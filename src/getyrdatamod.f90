module getyrdatamod

use parametersmod, only : stdout,stderr

! get all data that changes on a yearly basis

implicit none

public  :: getco2
public  :: getdata
private :: getclimate
private :: getlucc
private :: gettopo

contains

! ------------------------------------------------------------------------------------------------------------

subroutine getco2(cal_year,years)

use netcdf
use typesizes
use parametersmod,  only : dp
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : co2file,co2vect
use calendarmod,    only : timestruct,ymdt2jd

implicit none

! arguments
integer, intent(in) :: cal_year
integer, intent(in) :: years

! local variables

integer :: ncid
integer :: varid
integer :: dimid

integer :: tlen

integer, dimension(1) :: tloc
integer :: srtt
real(dp), allocatable, dimension(:) :: ghgtime

integer :: yearCE

type(timestruct) :: baseyear
type(timestruct) :: thisyear

character(80) :: basetimestring

real(dp) :: dt

integer :: i,t

! --------------------

ncstat = nf90_open(co2file,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ncid,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ncid,dimid,len=tlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(ghgtime(tlen))

ncstat = nf90_inq_varid(ncid,'time',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,ghgtime)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_att(ncid,varid,'units',basetimestring)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! ----
! calculate Julian day for base year 
! (this should read from the units in the file and be done only once in the init subroutine)

baseyear = timestruct(1950,1,1)

! calculate Julian day for base year

call ymdt2jd(baseyear)

! ----
! calculate Julian day for current year 

! convert yr BP (cal_year) to year CE

yearCE = 1950 - cal_year

if (yearCE <= 0) yearCE = yearCE - 1

thisyear = timestruct(yearCE,1,1)

call ymdt2jd(thisyear)

! get position in GHG file time dimension

dt = thisyear%jd - baseyear%jd

tloc = minloc(abs(ghgtime - dt))

srtt = tloc(1)

write(0,*)'getting CO2',cal_year,yearCE,srtt

! ! check to make sure the requested calendar year for the beginning of the run is available in the dataset
! 
! if (cal_year > times(1)) then
!   write(stdout,*)'ERROR: the requested starting year for the run is earlier than the first year of data in the CO2 file'
!   write(stdout,*)cal_year,times(1)
!   stop
! end if
! 
! ! scan the time vector to figure out where to start getting the co2 vector from  
! 
! srt = minloc(abs(times - cal_year))

allocate(co2vect(years))

ncstat = nf90_inq_varid(ncid,'co2',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
ncstat = nf90_get_var(ncid,varid,co2vect,start=[srtt],count=[years])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

! write(stdout,*)'co2',srtt,transientyears,cal_year

! i = srtt
! do t = 1,transientyears
!   
!   write(0,*)t,i,ghgtime(i),co2vect(t)
!   
!   i = i + 1
!   
! end do
! 
! stop

end subroutine getco2

! ------------------------------------------------------------------------------------------------------------

subroutine getdata(ncells,year,cal_year,firstyear,time0,in_master)

use mpistatevarsmod,only : inputdata
use iovariablesmod, only : ibuf,soil,lucc,climateyears,co2vect,calcforagers,startyr_foragers,dosoilco2
use orbitmod,       only : calcorbitpars,orbitpars
use parametersmod,  only : sp

implicit none

! arguments

integer, intent(in) :: ncells
integer, intent(in) :: year
integer, intent(in) :: cal_year
integer, intent(in) :: firstyear

integer, intent(inout) :: time0

type(inputdata), dimension(:), intent(inout) :: in_master

! local variables

integer :: i
integer :: x
integer :: y
integer :: lyear
integer :: itime
integer :: nl

type(orbitpars) :: orbit

! ------------------------------

lyear = mod(year,climateyears)

if (lyear == 0) lyear = climateyears

itime = 1 + 12 * (lyear - 1)

call calcorbitpars(cal_year,orbit)

! if (year == 1 .or..not.in_master(1)%spinup) then  ! we need to get annual topo data
  call gettopo(year,cal_year)

  if (calcforagers) then
    call getforg(cal_year)
    
    in_master%startyr_foragers = startyr_foragers

  else if (lucc) then
    call getlucc(cal_year)

  else  ! run without people

    do i = 1,3
      in_master%human%popd(i)  = 0.
    end do

    in_master%human%landuse(1) = 1.

    do i = 2,7
      in_master%human%landuse(i) = -1.
    end do

  end if

! end if

call getclimate(itime,time0)

! fill the values of in_master here

in_master%dosoilco2 = dosoilco2 ! FLAG: added by AK

if (year == 1) then
  do i = 1,ncells

    x = in_master(i)%xpos
    y = in_master(i)%ypos
    
    in_master(i)%elev      = soil(x,y)%elv
    in_master(i)%slope     = soil(x,y)%slopeangle  ! FLAG: added by MP, 13.12.2011
    in_master(i)%landf     = soil(x,y)%landf       ! FLAG: added by MP, 09.08.2012
    
    !  in_master expects two layers of soil. if soil has more layers, take an average of the top 30 cm (two layers) and the rest
    
    nl = size(soil(x,y)%zpos)
    
    if (nl > 2) then
    
      in_master(i)%soil%sand(1) = sum(soil(x,y)%sand(1:2)) / 2
      in_master(i)%soil%sand(2) = sum(soil(x,y)%sand(3:nl)) / (nl-2)

      in_master(i)%soil%clay(1) = sum(soil(x,y)%clay(1:2)) / 2
      in_master(i)%soil%clay(2) = sum(soil(x,y)%clay(3:nl)) / (nl-2)

      in_master(i)%soil%orgm(1) = sum(soil(x,y)%orgm(1:2)) / 2
      in_master(i)%soil%orgm(2) = sum(soil(x,y)%orgm(3:nl)) / (nl-2)

      in_master(i)%soil%zpos(1) = sum(soil(x,y)%zpos(1:2)) / 2
      in_master(i)%soil%zpos(2) = sum(soil(x,y)%zpos(3:nl)) / (nl-2)
    
    else
    
      in_master(i)%soil%sand = soil(x,y)%sand
      in_master(i)%soil%clay = soil(x,y)%clay
      in_master(i)%soil%orgm = soil(x,y)%orgm
      in_master(i)%soil%zpos = soil(x,y)%zpos

    end if
    
!      write(stdout,*)'SETTING OUTPUT SOIL' 
!      write(stdout,*)in_master(i)%soil%sand
!          write(stdout,*)in_master(i)%soil%clay
!          write(stdout,*)in_master(i)%soil%orgm
!          write(stdout,*)in_master(i)%soil%zpos

    ! FLAGFLAGFLAG
    
!     if(in_master(i)%slope /= in_master(i)%slope) in_master(i)%slope = 0.

  end do
end if

! write(0,*)'FLAG',firstyear,cal_year,1+firstyear-cal_year,co2vect(1+firstyear-cal_year)

do i = 1,ncells

  x = in_master(i)%xpos
  y = in_master(i)%ypos

  in_master(i)%co2          = co2vect(1+firstyear-cal_year)
  in_master(i)%climate%temp = ibuf(x,y)%temp
  in_master(i)%climate%prec = ibuf(x,y)%prec
  in_master(i)%climate%cldp = ibuf(x,y)%cldp
  in_master(i)%climate%wetd = ibuf(x,y)%wetd
  in_master(i)%climate%trng = ibuf(x,y)%trng
  in_master(i)%climate%wind = ibuf(x,y)%wind
  in_master(i)%climate%lght = ibuf(x,y)%lght
  in_master(i)%orbit%yrBP   = orbit%yrBP
  in_master(i)%orbit%ecc    = orbit%ecc
  in_master(i)%orbit%pre    = orbit%pre
  in_master(i)%orbit%perh   = orbit%perh
  in_master(i)%orbit%xob    = orbit%xob

  if (calcforagers) then
    in_master(i)%human%foragerPD = ibuf(x,y)%popd(1)
  else
    in_master(i)%human%foragerPD = 0.
  end if
  
  if (lucc) then
    in_master(i)%human%lu_turnover   = ibuf(x,y)%lu_turnover
    in_master(i)%human%popd          = ibuf(x,y)%popd          ! NB this is an array
    in_master(i)%human%landuse       =-1.                      ! initialize to missing
    in_master(i)%human%landuse(2:5)  = ibuf(x,y)%landuse(2:5)  ! cropland, pasture, rangeland, urban
    
    ! write(stdout,*)'input landuse',in_master(i)%human%landuse(1:2)
  else
    in_master(i)%human%popd        =  0          ! NB this is an array
    in_master(i)%human%landuse(1)  =  1  ! only natural landuse without lucc file
    in_master(i)%human%landuse(2:) = -1 
    
  end if

end do

end subroutine getdata

! ------------------------------------------------------------------------------------------------------------

subroutine getclimate(itime,time0)

use netcdf
use typesizes
use errormod,       only : netcdf_err,ncstat
use parametersmod,  only : sp,i2
use iovariablesmod, only : cfid,ibuf,timebuflen,climatemonths,inputclimlen,srtx,cntx,srty,cnty,varinfo,soil, &
                           cellmask,input_sp,input_i2,nclimv

implicit none

integer, intent(in)    :: itime
integer, intent(inout) :: time0

integer :: i
integer :: x
integer :: y
integer :: t0
integer :: t1

integer :: remainmon
integer :: tlen

integer(i2), parameter :: imissing = -32768
real(sp),    parameter :: rmissing =  -9999.

! -------------------------------------
! check if we need to read in data from the climate data file

t0 = mod(itime,timebuflen)

! write(stdout,*)'getyrdata',timebuflen,climatemonths,itime,time0,t0,t0+11

if (t0 == 1 .and. itime /= time0 .and. time0+11 /= climatemonths) then
  
  remainmon = 1 + climatemonths - itime
  tlen = min(remainmon,timebuflen)

  write(stdout,*)'read more data',itime,time0,remainmon,tlen
  
  ! ---
  ! read data from file
  
  do i = 1,nclimv
    
    ncstat = nf90_get_var(cfid,varinfo(i)%varid,input_i2,start=[srtx,srty,itime],count=[cntx,cnty,tlen])
    if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

    where (input_i2 /= imissing) input_sp(:,:,:,i) = real(input_i2) * varinfo(i)%scale_factor + varinfo(i)%add_offset

  end do

  ! ---
  
end if

t1 = t0 + 11

! transfer one year of climate data from ivals to ibuf

do y = 1,cnty
  do x = 1,cntx

    ibuf(x,y)%temp = input_sp(x,y,t0:t1,1)
    ibuf(x,y)%prec = input_sp(x,y,t0:t1,2)
    ibuf(x,y)%cldp = input_sp(x,y,t0:t1,3)
    ibuf(x,y)%wetd = input_sp(x,y,t0:t1,4)
    ibuf(x,y)%trng = input_sp(x,y,t0:t1,5)
    ibuf(x,y)%wind = input_sp(x,y,t0:t1,6)
    ibuf(x,y)%lght = input_sp(x,y,t0:t1,7)
    
    if (any(input_sp(x,y,t0:t1,:) == rmissing) .or. soil(x,y)%sand(1) < 0. .or. soil(x,y)%landf <= 0.) then
      cellmask(x,y) = .false.
    else
      cellmask(x,y) = .true.
    end if

!      write(stdout,*)x,y,soil(x,y)%landf,cellmask(x,y)
    
    ! if (ibuf(x,y)%temp(1) /= rmissing .and. soil(x,y)%sand(1) >= 0. .and. soil(x,y)%landf > 0.) cellmask(x,y) = .true.  
    
  end do
end do

time0 = itime

end subroutine getclimate

! ------------------------------------------------------------------------------------------------------------

subroutine getlucc(cal_year)

! this subroutine now modified to work only with integer(2) input files

use netcdf
use typesizes
use parametersmod,  only : i2,sp,dp
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : ibuf,srtx,cntx,srty,cnty,popfid,popdtime,nolanduse
use calendarmod,    only : timestruct,ymdt2jd
use utilitiesmod,   only : pos

implicit none

integer, intent(in) :: cal_year   ! should be calendar year in yr BP (1950)

integer :: x,y

real(sp),           dimension(cntx,cnty) :: rvals
integer(i2), allocatable, dimension(:,:) :: svals

real(sp) :: scale_factor
real(sp) :: add_offset

integer(i2) :: imissing

integer :: varid

integer :: srtt
integer, dimension(1) :: tloc

real :: lu_turn_yrs
real :: lu_turnover

integer :: xtype

integer :: yearCE

type(timestruct) :: baseyear
type(timestruct) :: thisyear

real(dp) :: dt

! -------------------------------------------------------------------------------------
! new version 02.2024, the time coordinate in the land use file is in days since a base time

! ----
! calculate Julian day for base year (this should be done only once in the init subroutine)

baseyear = timestruct(1950,1,1)

! calculate Julian day for base year
! subroutine populates the data structure with the Julian day corresponding to the YMD date

call ymdt2jd(baseyear)

! ----
! calculate Julian day for current year 

! convert yr BP (cal_year) to year CE

yearCE = 1950 - cal_year

if (yearCE <= 0) yearCE = yearCE - 1

thisyear = timestruct(yearCE,1,1)

call ymdt2jd(thisyear)

! get position in land use file time dimension

dt = thisyear%jd - baseyear%jd

! tloc = minloc(abs(popdtime - dt))
! 
! srtt = tloc(1)

srtt = pos(popdtime,dt)

! write(0,*)'landuse tpos: ',cal_year,yearCE,dt,srtt

! -----------------
! intensive land use: includes crops, pasture, rangeland, and urban

rvals = 0.

! ---------
! cropland

ncstat = nf90_inq_varid(popfid,'cropland',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_variable(popfid,varid,xtype=xtype)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

if (xtype == nf90_short) then
  
  ncstat = nf90_get_att(popfid,varid,'scale_factor',scale_factor)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_att(popfid,varid,'add_offset',add_offset)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_att(popfid,varid,'missing_value',imissing)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  allocate(svals(cntx,cnty))
  
  ncstat = nf90_get_var(popfid,varid,svals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  where (svals /= imissing) rvals = real(svals) * scale_factor + add_offset
  
  deallocate(svals)
  
else if (xtype == nf90_float) then
  
  ncstat = nf90_get_var(popfid,varid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

else  ! 

  write(stdout,*)'error, the landuse variable is in an invalid type!  ',xtype
  stop

end if

do y = 1,cnty
 do x = 1,cntx
   ibuf(x,y)%landuse(2) = rvals(x,y)
 end do
end do

! ---------
! pasture

ncstat = nf90_inq_varid(popfid,'pasture',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_variable(popfid,varid,xtype=xtype)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

if (xtype == nf90_short) then
  
  ncstat = nf90_get_att(popfid,varid,'scale_factor',scale_factor)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_att(popfid,varid,'add_offset',add_offset)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  allocate(svals(cntx,cnty))
  
  ncstat = nf90_get_var(popfid,varid,svals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
  where (svals /= imissing) rvals = rvals + real(svals) * scale_factor + add_offset
  
  deallocate(svals)
  
else if (xtype == nf90_float) then
  
  ncstat = nf90_get_var(popfid,varid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

else  ! 

  write(stdout,*)'error, the landuse variable is in an invalid type!  ',xtype
  stop

end if

do y = 1,cnty
 do x = 1,cntx
   ibuf(x,y)%landuse(3) = rvals(x,y)
 end do
end do

! ---------
! rangeland

ncstat = nf90_inq_varid(popfid,'rangeland',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_variable(popfid,varid,xtype=xtype)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

if (xtype == nf90_short) then
  
  ncstat = nf90_get_att(popfid,varid,'scale_factor',scale_factor)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_att(popfid,varid,'add_offset',add_offset)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  allocate(svals(cntx,cnty))
  
  ncstat = nf90_get_var(popfid,varid,svals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
    
  where (svals /= imissing) rvals = rvals + real(svals) * scale_factor + add_offset
  
  deallocate(svals)
  
else if (xtype == nf90_float) then
  
  ncstat = nf90_get_var(popfid,varid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

else  ! 

  write(stdout,*)'error, the landuse variable is in an invalid type!  ',xtype
  stop

end if

do y = 1,cnty
 do x = 1,cntx
   ibuf(x,y)%landuse(4) = rvals(x,y)
 end do
end do

! ---------
! urban

ncstat = nf90_inq_varid(popfid,'urban',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_variable(popfid,varid,xtype=xtype)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

if (xtype == nf90_short) then
  
  ncstat = nf90_get_att(popfid,varid,'scale_factor',scale_factor)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_get_att(popfid,varid,'add_offset',add_offset)
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

  allocate(svals(cntx,cnty))
  
  ncstat = nf90_get_var(popfid,varid,svals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
    
  where (svals /= imissing) rvals = rvals + real(svals) * scale_factor + add_offset
  
  deallocate(svals)
  
else if (xtype == nf90_float) then
  
  ncstat = nf90_get_var(popfid,varid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

else  ! 

  write(stdout,*)'error, the landuse variable is in an invalid type!  ',xtype
  stop

end if

do y = 1,cnty
 do x = 1,cntx
   ibuf(x,y)%landuse(5) = rvals(x,y)
 end do
end do

! -----------------
! land use turnover

lu_turn_yrs = 0.

if (lu_turn_yrs > 0.) then
  lu_turnover = 1. / lu_turn_yrs
else
  lu_turnover = 0.
end if

do y = 1,cnty
  do x = 1,cntx
    ibuf(x,y)%lu_turnover = lu_turnover
  end do
end do

! -----------------
! population density

ncstat = nf90_inq_varid(popfid,'hunter_gatherers',varid)

if (ncstat == -49) then  ! no population data in the input dataset, ignore the rest of this routine

  return

else if (ncstat /= nf90_noerr) then 

  call netcdf_err(ncstat)  ! write the error message and abort

end if

! otherwise, retrieve the data

ncstat = nf90_get_var(popfid,varid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do y = 1,cnty
  do x = 1,cntx
    ibuf(x,y)%popd(1) = rvals(x,y)
  end do
end do

ncstat = nf90_inq_varid(popfid,'farmers',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(popfid,varid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do y = 1,cnty
  do x = 1,cntx
    ibuf(x,y)%popd(2) = rvals(x,y)
  end do
end do

ncstat = nf90_inq_varid(popfid,'pastoralists',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(popfid,varid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
do y = 1,cnty
  do x = 1,cntx
    ibuf(x,y)%popd(3) = rvals(x,y)
  end do
end do

end subroutine getlucc

! ------------------------------------------------------------------------------------------------------------

subroutine getforg(cal_year)

use netcdf
use typesizes
use parametersmod,  only : sp
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : ibuf,srtx,cntx,srty,cnty,popfid,popdtime

implicit none

integer, intent(in) :: cal_year   ! should be calendar year in yr BP (1950)

integer :: varid

integer :: srtt
integer, dimension(1) :: tloc

real(sp), dimension(cntx,cnty) :: rvar

! -------------------------------------------------------------------------------------
! potential density of hunter-gatherers (based on archaeological site density)

tloc = minloc(abs(popdtime - cal_year))

srtt = tloc(1)
! write (*,*) "getforg(), srtt: ", srtt
ncstat = nf90_inq_varid(popfid,'foragerPD',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(popfid,varid,rvar,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ibuf%popd(1) = rvar
! write (*,*) "getforg() minmaxpopd1: ", minval(ibuf%popd(1)), maxval(ibuf%popd(1))
end subroutine getforg

! ------------------------------------------------------------------------------------------------------------

subroutine gettopo(year,cal_year)

use ieee_arithmetic, only : ieee_is_nan
use netcdf
use typesizes
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : ibuf,srtx,cntx,srty,cnty,topofid,elvid,slopeid,landfid,topotime,soil
use parametersmod,  only : i2,sp

implicit none

integer, intent(in) :: year       ! year counter for run years
integer, intent(in) :: cal_year   ! should be calendar year in yr BP (1950)

integer(i2), dimension(cntx,cnty) :: ivals
real(sp),    dimension(cntx,cnty) :: rvals

integer :: srtt
integer, dimension(1) :: tloc

!  -------------------------
!  get elevation, slope, and land fraction

! slope does not change from year to year, so only get it once
if (year == 1) then
  ncstat = nf90_get_var(topofid,slopeid,soil%slopeangle,start=[srtx,srty],count=[cntx,cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
end if

! check to make sure the requested calendar year for the beginning of the run is available in the dataset

if (cal_year > topotime(1)) then
  write(stdout,*)'WARNING: the requested starting year for the run is earlier than the first year of data in the topofile'
  write(stdout,*)'year cal BP, topotime: ',cal_year,topotime(1)
  write(stdout,*)'using topo data from nearest year available in data set'
end if  

! scan the time vector to figure out where to start getting the elevation- and landf- vector from 

tloc = minloc(abs(topotime - cal_year),1)

srtt = tloc(1)

! write(stdout,*)'reading topo data',srtx,srty,srtt,cntx,cnty

ncstat = nf90_get_var(topofid,elvid,ivals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

soil%elv = real(ivals)

ncstat = nf90_get_var(topofid,landfid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat) 

where (ieee_is_nan(rvals)) rvals = -9999.

soil%landf = rvals

!  write(stdout,*) 'end of gettopo',soil%landf

end subroutine gettopo

! ------------------------------------------------------------------------------------------------------------

end module getyrdatamod
