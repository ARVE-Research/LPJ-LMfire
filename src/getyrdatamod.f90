module getyrdatamod

!get all data that changes on a yearly basis

implicit none

public  :: getco2
public  :: getdata
private :: getclimate
private :: getlucc
private :: gettopo

contains

!------------------------------------------------------------------------------------------------------------

subroutine getco2(cal_year,transientyears)

use netcdf
use typesizes
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : co2file,co2vect

implicit none

!arguments
integer, intent(in) :: cal_year
integer, intent(in) :: transientyears

!local variables

integer :: ncid
integer :: varid
integer :: dimid

integer :: tlen

integer, dimension(1) :: srt
integer, allocatable, dimension(:) :: times

!--------------------

ncstat = nf90_open(co2file,nf90_nowrite,ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(ncid,'time',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(ncid,dimid,len=tlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(times(tlen))

ncstat = nf90_inq_varid(ncid,'time',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(ncid,varid,times)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!check to make sure the requested calendar year for the beginning of the run is available in the dataset

if (cal_year > times(1)) then
  write(0,*)'ERROR: the requested starting year for the run is earlier than the first year of data in the CO2 file'
  write(0,*)cal_year,times(1)
  stop
end if

!scan the time vector to figure out where to start getting the co2 vector from  

srt = minloc(abs(times - cal_year))

allocate(co2vect(transientyears))

ncstat = nf90_inq_varid(ncid,'co2',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
ncstat = nf90_get_var(ncid,varid,co2vect,start=[srt],count=[transientyears])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!write(0,*)'co2',srt,transientyears,cal_year

end subroutine getco2

!------------------------------------------------------------------------------------------------------------

subroutine getdata(ncells,year,cal_year,firstyear,time0,in_master)

use mpistatevarsmod,only : inputdata
use iovariablesmod, only : ibuf,soil,lucc,climateyears,co2vect,calcforagers
use orbitmod,       only : calcorbitpars,orbitpars
use parametersmod,  only : sp

implicit none

!arguments

integer, intent(in) :: ncells
integer, intent(in) :: year
integer, intent(in) :: cal_year
integer, intent(in) :: firstyear

integer, intent(inout) :: time0

type(inputdata), dimension(:), intent(inout) :: in_master

!local variables

integer :: i
integer :: x
integer :: y
integer :: lyear
integer :: itime

type(orbitpars) :: orbit

!------------------------------

lyear = mod(year,climateyears)

if (lyear == 0) lyear = climateyears

itime = 1 + 12 * (lyear - 1)

call calcorbitpars(cal_year,orbit)

call getclimate(itime,time0)

if (year == 1 .or..not.in_master(1)%spinup) then  !we need to get annual topo data
  call gettopo(year,cal_year)

  if (calcforagers) then
    call getforg(cal_year)

  else if (lucc) then
    call getlucc(cal_year)

  else  !run without people

    do i = 1,3
      in_master%human%popd(i)  = 0.
    end do

    in_master%human%landuse(1) = 1.

    do i = 2,7
      in_master%human%landuse(i) = -1.
    end do

  end if

end if

!fill the values of in_master here

if (year == 1) then
  do i = 1,ncells

    x = in_master(i)%xpos
    y = in_master(i)%ypos
    
    in_master(i)%elev      = soil(x,y)%elv
    in_master(i)%slope     = soil(x,y)%slopeangle  !FLAG: added by MP, 13.12.2011
    in_master(i)%landf     = soil(x,y)%landf       !FLAG: added by MP, 09.08.2012
    in_master(i)%soil%sand = soil(x,y)%sand
    in_master(i)%soil%clay = soil(x,y)%clay
    in_master(i)%soil%orgm = soil(x,y)%orgm
    in_master(i)%soil%zpos = soil(x,y)%zpos
    
    !FLAGFLAGFLAG
    
!    if(in_master(i)%slope /= in_master(i)%slope) in_master(i)%slope = 0.

  end do
end if

do i = 1,ncells

  x = in_master(i)%xpos
  y = in_master(i)%ypos

  in_master(i)%co2          = co2vect(1+firstyear-cal_year)
  in_master(i)%climate%temp = ibuf(x,y)%temp
  in_master(i)%climate%prec = ibuf(x,y)%prec
  in_master(i)%climate%cldp = ibuf(x,y)%cldp
  in_master(i)%climate%wetd = ibuf(x,y)%wetd
  in_master(i)%climate%trng = ibuf(x,y)%trng
  in_master(i)%climate%lght = 10**(ibuf(x,y)%lght) ! Pour la premier simulation on a mis: exp(ibuf(x,y)%lght)
  in_master(i)%climate%wind = ibuf(x,y)%wind
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
    in_master(i)%human%popd          = ibuf(x,y)%popd          !NB this is an array
    in_master(i)%human%landuse(1:2)  = ibuf(x,y)%cropfrac(1:2)
    in_master(i)%human%landuse(3:)   =-1.                      !all other types, set to missing for now
    
    !write(0,*)'input landuse',in_master(i)%human%landuse(1:2)
  else
    in_master(i)%human%popd        =  0          !NB this is an array
    in_master(i)%human%landuse(1)  =  1  !only natural landuse without lucc file
    in_master(i)%human%landuse(2:) = -1 
    
  end if

end do

end subroutine getdata

!------------------------------------------------------------------------------------------------------------

subroutine getclimate(itime,time0)

use netcdf
use typesizes
use errormod,       only : netcdf_err,ncstat
use parametersmod,  only : sp,i2
use iovariablesmod, only : cfid,ibuf,timebuflen,climatemonths,inputclimlen,srtx,cntx,srty,cnty,varinfo,soil,cellmask,input_sp,input_i2,nclimv

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

integer(i2), parameter :: missing  = -32768
real(sp),    parameter :: rmissing =  -9999.

!-------------------------------------
!check if we need to read in data from the climate data file

t0 = mod(itime,timebuflen)

!write(0,*)'getyrdata',timebuflen,climatemonths,itime,time0,t0,t0+11

if (t0 == 1 .and. itime /= time0 .and. time0+11 /= climatemonths) then
  
  remainmon = 1 + climatemonths - itime
  tlen = min(remainmon,timebuflen)

  write(0,*)'read more data',itime,time0,remainmon,tlen
  
  !---
  !read data from file
  
  do i = 1,nclimv

    ncstat = nf90_get_var(cfid,varinfo(i)%varid,input_i2,start=[srtx,srty,itime],count=[cntx,cnty,tlen])
    if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

    where (input_i2 /= missing) input_sp(:,:,:,i) = real(input_i2) * varinfo(i)%scale_factor + varinfo(i)%add_offset

  end do

  !---
  
end if

t1 = t0 + 11

!transfer one year of climate data from ivals to ibuf

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
    
    !if (ibuf(x,y)%temp(1) /= rmissing .and. soil(x,y)%sand(1) >= 0. .and. soil(x,y)%landf > 0.) cellmask(x,y) = .true.  
    
  end do
end do

time0 = itime


end subroutine getclimate

!------------------------------------------------------------------------------------------------------------

subroutine getlucc(cal_year)

!this subroutine now modified to work only with integer(2) input files

use netcdf
use typesizes
use parametersmod,  only : i2
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : ibuf,srtx,cntx,srty,cnty,lucctime,popfid,popdtime,nolanduse

implicit none

integer, intent(in) :: cal_year   !should be calendar year in yr BP (1950)

integer :: x,y

real,                     dimension(cntx,cnty) :: rvals
integer(i2), allocatable, dimension(:,:) :: svals

real :: scale_factor
real :: add_offset

integer :: varid

integer :: srtt
integer, dimension(1) :: tloc

real :: lu_turn_yrs
real :: lu_turnover

integer :: xtype

!-------------------------------------------------------------------------------------
!this part removed because now we use a consolidate population and land use file

!tloc = minloc(abs(lucctime - cal_year))

!srtt = tloc(1)

! if (nolanduse) then
!
!   do y = 1,cnty
!     do x = 1,cntx
!       ibuf(x,y)%cropfrac = 0.
!     end do
!   end do
!
! else
!
! !-----------------
! !unusable fraction
!
! ncstat = nf90_inq_varid(luccfid,'unusable',varid)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_var(luccfid,varid,svals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! ncstat = nf90_get_att(luccfid,varid,'scale_factor',scale_factor)
! if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
!
! do y = 1,cnty
!   do x = 1,cntx
!     ibuf(x,y)%cropfrac(1) = scale_factor * real(svals(x,y))
!   end do
! end do
!end if
!-------------------------------------------------------------------------------------

tloc = minloc(abs(popdtime - cal_year))

srtt = tloc(1)

!-----------------
!intensive land use

ncstat = nf90_inq_varid(popfid,'land_use',varid)
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
    
  rvals = real(svals) * scale_factor + add_offset
  
  deallocate(svals)
  
else if (xtype == nf90_float) then
  
  ncstat = nf90_get_var(popfid,varid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

else  !

  write(0,*)'error, the landuse variable is in an invalid type! ',xtype
  stop

end if

do y = 1,cnty
  do x = 1,cntx
    ibuf(x,y)%cropfrac(1) = 0.          !unusable fraction set to zero for now
    ibuf(x,y)%cropfrac(2) = rvals(x,y)  !intensive land use
  end do
end do

!write(0,'(a,2i5,5f10.4)')'getlucc',cal_year,srtt,ibuf(1,1)%cropfrac,ibuf(1,1)%popd
!write(0,*)'getlucc',cal_year,srtt

!-----------------
!land use turnover

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

!-----------------
!population density

ncstat = nf90_inq_varid(popfid,'hunter_gatherers',varid)

if (ncstat == -49) then  !no population data in the input dataset, ignore the rest of this routine

  return

else if (ncstat /= nf90_noerr) then 

  call netcdf_err(ncstat)  !write the error message and abort

end if

!otherwise, retrieve the data

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

!------------------------------------------------------------------------------------------------------------

subroutine getforg(cal_year)

use netcdf
use typesizes
use parametersmod,  only : sp
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : ibuf,srtx,cntx,srty,cnty,popfid,popdtime

implicit none

integer, intent(in) :: cal_year   !should be calendar year in yr BP (1950)

integer :: varid

integer :: srtt
integer, dimension(1) :: tloc

real(sp), dimension(cntx,cnty) :: rvar

!-------------------------------------------------------------------------------------
!potential density of hunter-gatherers (based on archaeological site density)

tloc = minloc(abs(popdtime - cal_year))

srtt = tloc(1)
!write (*,*) "getforg(), srtt: ", srtt
ncstat = nf90_inq_varid(popfid,'foragerPD',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(popfid,varid,rvar,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ibuf%popd(1) = rvar
!write (*,*) "getforg() minmaxpopd1: ", minval(ibuf%popd(1)), maxval(ibuf%popd(1))
end subroutine getforg

!------------------------------------------------------------------------------------------------------------

subroutine gettopo(year,cal_year)

use netcdf
use typesizes
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : ibuf,srtx,cntx,srty,cnty,topofid,elvid,slopeid,landfid,topotime,soil
use parametersmod,  only : i2,sp

implicit none

integer, intent(in) :: year       !year counter for run years
integer, intent(in) :: cal_year   !should be calendar year in yr BP (1950)

integer(i2), dimension(cntx,cnty) :: ivals
real(sp),    dimension(cntx,cnty) :: rvals

integer :: srtt
integer, dimension(1) :: tloc

!-------------------------
!get elevation, slope, and land fraction

!slope does not change from year to year, so only get it once
if (year == 1) then
  ncstat = nf90_get_var(topofid,slopeid,soil%slopeangle,start=[srtx,srty],count=[cntx,cnty])
  if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
end if

!check to make sure the requested calendar year for the beginning of the run is available in the dataset

if (cal_year > topotime(1)) then
  write(0,*)'WARNING: the requested starting year for the run is earlier than the first year of data in the topofile'
  write(0,*)cal_year,topotime(1)
  write(0,*)'using topo data from nearest year available in data set'
end if  

!scan the time vector to figure out where to start getting the elevation- and landf- vector from 

tloc = minloc(abs(topotime - cal_year),1)

srtt = tloc(1)

!write(0,*)'reading topo data',srtx,srty,srtt,cntx,cnty

ncstat = nf90_get_var(topofid,elvid,ivals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

soil%elv = real(ivals)

ncstat = nf90_get_var(topofid,landfid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat) 

soil%landf = rvals

!write(0,*) 'end of gettopo'

end subroutine gettopo

!------------------------------------------------------------------------------------------------------------

end module getyrdatamod
