module snowmod

implicit none

contains

!-------------------------------------------------------------------------------------------------------------

subroutine snow(dtemp,dprec,treefrac,snow0,snowpack,dmelt,idx)

!SUBROUTINE SNOW
!Adjust daily precipitation by dmelt(d) and accumulation in snowpack

!Ref: Haxeltine & Prentice 1996
!Ref: Semadini-Davies 1997
!Ref: Tarboton & Luce 1996 (Utah Energy Balance Snow Accumulation and Melt Model)

use parametersmod, only : i8

implicit none

integer(i8) :: idx

!PARAMETERS:

real, parameter :: utsnow =  3.0
real, parameter :: ltsnow = -1.0

!ARGUMENTS:      
real, intent(in) :: treefrac  !fractional coverage of trees in the gridcell
real, intent(in) :: snow0     !the amount of snowpack on the last day of the previous year
real, intent(in),  dimension(:) :: dtemp
real, intent(inout),  dimension(:) :: dprec

real, intent(out), dimension(:) :: snowpack
real, intent(out), dimension(:) :: dmelt

!LOCAL VARIABLES:

real :: kmelt  !the degree-day melt fraction in mm degC-1 d-1, proportional to tree cover
real :: dtsnow !the difference in the snow threshold temperatures

integer :: d   !day of the year
integer :: yd  !the day before

real :: newsnow

!------------------------------

dtsnow = utsnow - ltsnow

kmelt = 6. - 4. * treefrac

snowpack(365) = snow0

do d = 1,365
  
  if (d > 1) then
    yd = d-1
  else
    yd = 365
  end if
  
  !Calculate snow melt and new snow for today
  
  if (dtemp(d) >= utsnow) then
    newsnow = 0.
  else if (dtemp(d) <= ltsnow) then
    newsnow = dprec(d)
  else
    newsnow = dprec(d) - dprec(d) * (dtemp(d) - ltsnow) / dtsnow
  end if
  
  if (newsnow < 0.) then
    write(0,*)'invalid newsnow',newsnow
    stop
  end if
  
  dprec(d) = dprec(d) - newsnow

  snowpack(d) = snowpack(yd) + newsnow

  !for some reason the maximum snow pack is limited to 5000mm water equivalent, maybe otherwise you get a glacier!
  snowpack(d) = min(snowpack(d),5000.)

  if (dtemp(d) > 0.) then
    dmelt(d) = kmelt * dtemp(d)
    snowpack(d) = snowpack(yd) - dmelt(d)
  else
    dmelt(d) = 0.0
  end if
  
  if (snowpack(d) < 0.0) then
    snowpack(d) = 0.0
    dmelt(d) = snowpack(yd)
  end if
  
  !if (idx == 1) write(*,'(a,i5,7f8.2)')'snow',d,snow0,kmelt,dtemp(d),dprec(d),newsnow,snowpack(d),dmelt(d)
  !read(*,*)

end do

end subroutine snow

!-------------------------------------------------------------------------------------------------------------

end module snowmod
