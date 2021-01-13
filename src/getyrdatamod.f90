module getyrdatamod

use parametersmod, only : stdout,stderr,i2,sp

!get all data that changes on a yearly basis

implicit none

public  :: getco2
public  :: getdata
private :: getclimate
private :: gethumans
private :: gettopo

integer(i2), parameter :: missing  = -32768
real(sp),    parameter :: rmissing =  -9999.

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
  write(stdout,*)'ERROR: the requested starting year for the run is earlier than the first year of data in the CO2 file'
  write(stdout,*)cal_year,times(1)
  stop
end if

!scan the time vector to figure out where to start getting the co2 vector from  

srt = minloc(abs(times - cal_year))

allocate(co2vect(transientyears))

ncstat = nf90_inq_varid(ncid,'co2',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

write(stdout,*)'getco2',srt,transientyears
  
ncstat = nf90_get_var(ncid,varid,co2vect,start=[srt],count=[transientyears])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_close(ncid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!write(stdout,*)'co2',srt,transientyears,cal_year

end subroutine getco2

!------------------------------------------------------------------------------------------------------------

subroutine getdata(ncells,year,cal_year,firstyear,time0,in_master)

use mpistatevarsmod,only : inputdata
use iovariablesmod, only : ibuf,soil,calchumans,climateyears,co2vect,startyr_foragers
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

if (year == 1 .or..not.in_master(1)%spinup) then  !we need to get annual topo data
  call gettopo(year,cal_year)

  if (calchumans) then

    call gethumans(cal_year)

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

call getclimate(itime,time0)

!fill the values of in_master here

if (year == 1) then
  do i = 1,ncells

    x = in_master(i)%xpos
    y = in_master(i)%ypos
    
    in_master(i)%elev      = soil(x,y)%elv
    in_master(i)%slope     = soil(x,y)%slopeangle  !FLAG: added by MP, 13.12.2011
    in_master(i)%landf     = soil(x,y)%landf       !FLAG: added by MP, 09.08.2012
    in_master(i)%soil%sand = soil(x,y)%sand(1:2)	!ONLY FIRST 2 OF 6 ARE USED
    in_master(i)%soil%clay = soil(x,y)%clay(1:2)	!ONLY FIRST 2 OF 6 ARE USED
    in_master(i)%soil%orgm = soil(x,y)%orgm(1:2)	!ONLY FIRST 2 OF 6 ARE USED
    in_master(i)%soil%zpos = soil(x,y)%zpos(1:2)	!ONLY FIRST 2 OF 6 ARE USED
    
    
!         write(stdout,*)'SETTING OUTPUT SOIL' 
!     write(stdout,*)in_master(i)%soil%sand
! 				c
! 				write(stdout,*)in_master(i)%soil%orgm
! 				write(stdout,*)in_master(i)%soil%zpos

    


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
  in_master(i)%climate%lght = max(ibuf(x,y)%lght,0.)  ! also correct issue with incorrect unpacking, for canada simulations was = 10**(ibuf(x,y)%lght) ! Pour la premier simulation on a mis: exp(ibuf(x,y)%lght)
  in_master(i)%climate%wind = ibuf(x,y)%wind
  in_master(i)%orbit%ecc    = orbit%ecc
  in_master(i)%orbit%pre    = orbit%pre
  in_master(i)%orbit%perh   = orbit%perh
  in_master(i)%orbit%xob    = orbit%xob
  
  if (calchumans) then
  
    in_master(i)%human%hg_present   = ibuf(x,y)%hg_present
    in_master(i)%human%landuse(1:3) = ibuf(x,y)%landuse(1:3)
    in_master(i)%human%landuse(4:)  =-1.                      !set all other types to missing, for now
    
  else

    in_master(i)%human%hg_present   = .false.
    in_master(i)%human%landuse(1:2) = 0.
    in_master(i)%human%landuse(3:)  =-1.
    
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

!-------------------------------------
!check if we need to read in data from the climate data file

t0 = mod(itime,timebuflen)

!write(stdout,*)'getyrdata',timebuflen,climatemonths,itime,time0,t0,t0+11

if (t0 == 1 .and. itime /= time0 .and. time0+11 /= climatemonths) then
  
  remainmon = 1 + climatemonths - itime
  tlen = min(remainmon,timebuflen)

  write(stdout,*)'read more data',itime,time0,remainmon,tlen
  
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

!     write(stdout,*)x,y,soil(x,y)%landf,cellmask(x,y)
    
    !if (ibuf(x,y)%temp(1) /= rmissing .and. soil(x,y)%sand(1) >= 0. .and. soil(x,y)%landf > 0.) cellmask(x,y) = .true.  
    
  end do
end do

time0 = itime

end subroutine getclimate

!------------------------------------------------------------------------------------------------------------

subroutine gethumans(cal_year)

!this subroutine now modified to work only with integer(2) input files

! 2018.05: this subroutine now reads three quantities: crop and pasture land use fraction, and hunter-gatherer presence
! I believe that is all that is necessary for calculating managed fire and h-g dynamics

use netcdf
use typesizes
use parametersmod,  only : sp,i1,i2
use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : ibuf,srtx,cntx,srty,cnty,humanfid,popdtime,nolanduse,lu_sf,lu_ao

implicit none

integer, intent(in) :: cal_year   !should be calendar year in yr BP (1950)

integer :: x,y

integer(i1), dimension(cntx,cnty) :: bvals
integer(i2), dimension(cntx,cnty) :: svals
real(sp),    dimension(cntx,cnty) :: rvals

real :: scale_factor
real :: add_offset

integer :: varid

integer :: srtt
integer, dimension(1) :: tloc

real :: lu_turn_yrs
real :: lu_turnover

integer :: xtype

!-------------------------------------------------------------------------------------

tloc = minloc(abs(popdtime - cal_year))

srtt = tloc(1)

!-----------------
! crop land use fraction

ncstat = nf90_inq_varid(humanfid,'lu_crop',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(humanfid,varid,svals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do y = 1,cnty
  do x = 1,cntx
    if (svals(x,y) /= missing) then
      ibuf(x,y)%landuse(1) = real(svals(x,y)) * lu_sf + lu_ao  ! cropland land use
    else
      ibuf(x,y)%landuse(1) = 0.
    end if
  end do
end do

!-----------------
! pasture land use fraction

ncstat = nf90_inq_varid(humanfid,'lu_past',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(humanfid,varid,svals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do y = 1,cnty
  do x = 1,cntx
    if (svals(x,y) /= missing) then
      ibuf(x,y)%landuse(2) = real(svals(x,y)) * lu_sf + lu_ao  ! pasture land use
    else
      ibuf(x,y)%landuse(2) = 0.
    end if
  end do
end do

!-----------------
! hunter-gatherer presence

ncstat = nf90_inq_varid(humanfid,'hg_presence',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(humanfid,varid,bvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
    
do y = 1,cnty
  do x = 1,cntx
  
    if (bvals(x,y) > 0) then
      ibuf(x,y)%hg_present = .true.
    else
      ibuf(x,y)%hg_present = .false.
    end if

  end do
end do

!-----------------
! burned fraction on managed land

ncstat = nf90_inq_varid(humanfid,'ag_burn',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(humanfid,varid,svals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

do y = 1,cnty
  do x = 1,cntx
    if (svals(x,y) /= missing) then
      ibuf(x,y)%landuse(3) = real(svals(x,y)) * lu_sf + lu_ao  ! burned fraction on managed land
    else
      ibuf(x,y)%landuse(3) = 0.
    end if
  end do
end do

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
write(stdout,*)'lu_turnover',lu_turnover
end subroutine gethumans

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
  write(stdout,*)'WARNING: the requested starting year for the run is earlier than the first year of data in the topofile'
  write(stdout,*)cal_year,topotime(1)
  write(stdout,*)'using topo data from nearest year available in data set'
end if  

!scan the time vector to figure out where to start getting the elevation- and landf- vector from 

tloc = minloc(abs(topotime - cal_year),1)

srtt = tloc(1)

!write(stdout,*)'reading topo data',srtx,srty,srtt,cntx,cnty

ncstat = nf90_get_var(topofid,elvid,ivals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

soil%elv = real(ivals)

ncstat = nf90_get_var(topofid,landfid,rvals,start=[srtx,srty,srtt],count=[cntx,cnty,1])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat) 

soil%landf = rvals

!write(stdout,*) 'end of gettopo'

end subroutine gettopo

!------------------------------------------------------------------------------------------------------------

end module getyrdatamod
