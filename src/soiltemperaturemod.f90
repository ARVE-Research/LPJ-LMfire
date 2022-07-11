module soiltemperaturemod
  
implicit none

public :: soiltemp

contains

!----------------------------------------------------------------------------------------------

subroutine soiltemp(soilpar,tair,tair0,tsoil,mw1,idx)

use parametersmod, only : sp,pi,i8

implicit none

!parameter

real(sp), parameter :: ang2m = 6. / pi
real(sp), parameter :: sec2m = 2.628

!arguments

real(sp), dimension(:), intent(in)  :: soilpar
real(sp), dimension(:), intent(in)  :: tair
real(sp), dimension(:), intent(in)  :: tair0
real(sp), dimension(:), intent(in)  :: mw1

real(sp), dimension(:), intent(out) :: tsoil

integer(i8), intent(in) :: idx

!local variables

real(sp) :: meanw1
real(sp) :: diffus
real(sp) :: tempthismonth
real(sp) :: templastmonth
real(sp) :: avetemp
real(sp) :: alag
real(sp) :: amp
real(sp) :: lag
real(sp) :: lagtemp

integer :: m

real(sp), dimension(24) :: tair2

!--------
!Annual cycle of soil temperatures follows surface temperatures with damped oscillation about a common mean, 
!and temporal lag.  Soil temperature at depth z and time t given by (Carslaw & Jaeger 1959; Eqn 5.52, Jury et al 1991):

!T(z,t) = Tav + A*exp(-z/d)*sin(omega*t - z/d)

!Tav       = average (base) air/soil temperature (avetemp)
!A         = amplitude of air temp fluctuation
!exp(-z/d) = fractional amplitude of temp fluctuation at soil depth z, relative to surface temp fluctuation (amp)
!z/d       = oscillation lag in angular units at soil depth z (alag)
!z         = soil depth = 0.25 m
!d         = sqrt(2*K/omega), damping depth
!K         = soil thermal diffusivity, m2/month (diffus)
!omega     = 2*pi/tau, angular frequency of oscillation
!tau       = oscillation period = 12 months

!Assume a sinusoidal cycle, but estimate soil temperatures by
!linear interpolation between previous monthly air temperatures to
!implement a lag relative to air temperature, and damping of soil
!temperature amplitude relative to air temperature amplitude.

!Soil thermal diffusivities are calculated by linear interpolation
!between estimates for 0, 15% and 100% soil water content from
!van Duin (1963) and Jury et al (1991) Fig 5.11.

!--------

!if (idx == 1) then
!  write(*,*)soilpar
!  do m = 1,12
!    write(*,*)m,tair0(m),tair(m),mw1(m)
!  end do
!end if

!Calculate mean annual water content in soil layer 1

meanw1 = sum(mw1) / 12.

!In case of zero soil water, return with soil temp = air temp

if (meanw1 == 0.) then
  tsoil = tair
  return
end if

!Interpolate thermal diffusivity function against soil water content

if (meanw1 < 0.15) then
  diffus = (soilpar(6) - soilpar(5)) / 0.15 * meanw1 + soilpar(5)
else
  diffus = (soilpar(7) - soilpar(6)) / 0.85 * (meanw1 - 0.15) + soilpar(6)
end if

!Convert diffusivity from mm2/s to m2/month
!multiplication by 1e-6 (-> m2/s) * 2.628e6 (s/month) = 2.628

diffus = diffus * sec2m

!Record 24 months of air temperatures in tair2

tair2(:12) = tair0
tair2(13:) = tair

!Calculate amplitude fraction and lag at soil depth 0.25 m

alag = 0.25 / sqrt(12. * diffus / pi)
amp  = exp(-alag)
lag  = alag * ang2m  !convert lag from angular units to months

!Calculate monthly soil temperatures for this year.  For each month,
!calculate average air temp for preceding 11 months plus the current month

do m = 13,24  !the current year's air temperature is stored in position 13-24
  
  avetemp = sum(tair2(m-11:m)) / 12.

  !Estimate air temperature "lag" months ago by linear interpolation
  !between air temperatures for this and last month

  tempthismonth = tair2(m)
  templastmonth = tair2(m-1)

  lagtemp = (tempthismonth - templastmonth) * (1. - lag) + templastmonth

  !Adjust amplitude of lagged air temp to give estimated soil temp

  tsoil(m-12) = avetemp + amp * (lagtemp - avetemp)

end do

end subroutine soiltemp

!----------------------------------------------------------------------------------------------

end module soiltemperaturemod

