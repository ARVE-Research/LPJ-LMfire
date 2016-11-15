module weathergenmod

!This module includes subroutines to calculate daily maximum and minimum temperature and cloud cover fraction
!based on an annual timeseries of monthly values of these variables
!The weather generator is based on the WGEN model (Richardson, 1981) with extension to use monthly summary 
!data from Geng et al., 1986, and Geng and Auburn 1986.
!Additional statistical relationships for both temperature and cloudiness have been produced
!by J.O. Kaplan using global weather station datasets (GSOD and global synoptic cloud reports).

!Coded in 2007-2009 by Jed Kaplan and Joe Melton, ARVE Group, EPFL/UVic, jed.kaplan@epfl.ch

use parametersmod, only : sp,dp,i4
use randomdistmod, only : randomstate

implicit none

public  :: weathergen_driver
public  :: metvars_in
public  :: metvars_out
public  :: rmsmooth
public  :: daily

private :: weathergen
private :: daymetvars
private :: meancv
private :: esat

!-------------------------------

type metvars_in
  
  real(sp) :: prec    !monthly total precipitation amount (mm)
  real(sp) :: wetd    !number of days in month with precipitation
  real(sp) :: wetf    !fraction of days in month with precipitation

  real(sp) :: tmin    !minumum temperture (C)
  real(sp) :: tmax    !maximum temperture (C)
  real(sp) :: cldf    !cloud fraction (0=clear sky, 1=overcast) (fraction)

  real(sp)               :: NI      !previous day's Nesterov index (degC 2)
  logical                :: pday    !precipitation status: true if the previous day was a rain day
  type(randomstate)      :: rndst   !state of the random number generator
  real(sp), dimension(3) :: resid   !previous day's weather residuals

end type metvars_in

type metvars_out
  
  real(sp) :: prec    !24 hour total precipitation (mm)
  real(sp) :: tmin    !24 hour minimum temperature (C)
  real(sp) :: tmax    !24 hour maximum temperature (C)
  real(sp) :: tdew    !dewpoint temperature (C)
  real(sp) :: cldf    !24 hour mean cloud cover fraction 0=clear sky, 1=overcast (fraction)
  real(sp) :: lght    !lightning flashes (flashes ha-1 day-1)
  real(sp) :: wind    !wind speed (m s-1)
  real(sp) :: dayl    !daylength (h)
  real(sp) :: srad    !downwelling surface shortwave radiation (J m-2 d-1)
  real(sp) :: dpet    !total potential evapotranspiration (mm)
  real(sp) :: NI      !Nesterov index (degC 2)

  logical                :: pday    !precipitation status
  type(randomstate)      :: rndst   !state of the random number generator, 15 elements
  real(sp), dimension(3) :: resid   !previous day's weather residuals
  
end type metvars_out

type daymetvars

  real(sp) :: tmax_mn     !maximum temperature monthly mean (K)
  real(sp) :: tmin_mn     !minimum temperature mothly mean (K)
  real(sp) :: cldf_mn     !mean cloud fraction (fraction)

  real(sp) :: tmax_cv     !coefficients of variation of corresponding variable above
  real(sp) :: tmin_cv     ! "
  real(sp) :: cldf_cv     ! "

end type daymetvars

contains

!------------------------------------------------------------------------------------------------------------

subroutine weathergen_driver(dtmin,dtmax,dcldf,prec,wetd,lght,met_out)

use parametersmod, only : sp,ndaymonth
use randomdistmod, only : ranur

implicit none

!arguments

real(sp), dimension(:), intent(in) :: dtmin
real(sp), dimension(:), intent(in) :: dtmax
real(sp), dimension(:), intent(in) :: dcldf
real(sp), dimension(:),  intent(in) :: prec
real(sp), dimension(:),  intent(in) :: wetd
real(sp), dimension(:),  intent(in) :: lght  !average (flashes ha-1 day-1)

type(metvars_out), dimension(:), intent(inout) :: met_out

!local variables

integer :: a
integer :: b
integer :: m
integer :: dyr
integer :: d
integer :: yesterday

real(sp), dimension(365) :: prob  !probability of lightning on this day (0-1)

real(sp) :: mprec_sim

type(metvars_in) :: met_in

integer :: i
integer :: wd

!------------------------------------------------
!copy the random number state from the previous year's last value - passed to this subroutine in met_out(1)

met_in%rndst = met_out(1)%rndst

!initialize lightning to zero

met_out%lght = 0.

!calculate daily meteorology

do m = 1,12

  if (m == 1) then
    a = 1
  else
    a = 1 + sum(ndaymonth(1:m-1))
  end if

  b = a + ndaymonth(m) - 1

  !-----

  met_in%prec = prec(m)
  met_in%wetd = wetd(m)
  met_in%wetf = wetd(m) / ndaymonth(m)            !++TODO: this doesn't take into account leap years

  i = 1
  do
    dyr = a
    do d = 1,ndaymonth(m)

      yesterday = dyr - 1
      if (yesterday == 0) yesterday = 365

      met_in%tmin  = dtmin(dyr)
      met_in%tmax  = dtmax(dyr)
      met_in%cldf  = dcldf(dyr)
      met_in%NI    = met_out(yesterday)%NI
      met_in%pday  = met_out(yesterday)%pday
      met_in%resid = met_out(yesterday)%resid

      call weathergen(met_in,met_out(dyr))       ! MP: weather generator gets the smoothed pseudo-daily values from rmsmooth

      met_in%rndst = met_out(dyr)%rndst
      
      dyr = dyr + 1

    end do !days of month loop

    !enforce a total monthly precip that is within 5% of the input (or 100 iterations)

    mprec_sim = sum(met_out(a:b)%prec)

    if (prec(m) == 0. .or. abs((prec(m) - mprec_sim) / prec(m)) < 0.05) exit

    i = i + 1
    
    if (i > 100) exit

  end do !end of conditional precip loop
  
  !-----------
  !disaggregate lightning flashes randomly only on days with precip, only if there is ligtning in the input file in this month
    
  if (lght(m) > 0.) then
  
    do wd = a,b
      if (met_out(wd)%prec > 0.) then
        prob(wd) = ranur(met_in%rndst)  !random real value from [0,1]
      else
        prob(wd) = 0.
      end if
    end do

    !disaggregate
    do wd = a,b
      if(sum(prob(a:b)) /= 0.) then     ! there is a likelihood for lightning in this month
         met_out(wd)%lght = lght(m) * ndaymonth(m) * prob(wd) / sum(prob(a:b)) !total flashes * fraction of total monthly strikes on this day
      else
         met_out(wd)%lght = 0.
      end if
    end do

  end if

  !-----------

end do  !end of month loop

end subroutine weathergen_driver

!------------------------------------------------------------------------------------------------------------

subroutine weathergen(met_in,met_out)

use parametersmod, only : sp,dp,i4,ndaymonth,tfreeze
use randomdistmod, only : ranur,ran_normal,ran_gamma

implicit none

!arguments

type(metvars_in),  intent(in)    :: met_in
type(metvars_out), intent(inout) :: met_out
  
!parameters

real(sp), parameter :: pmin = 2.16d0 / 0.83d0 !minimum value for pbar when using Geng linear relationship,
                                              !below this value we use a 1:1 line, see below
real(sp), parameter :: small = 5.e-5

real(sp), dimension(9), parameter :: corva = [ 0.567, 0.086,-0.002, &   !lag day-1 correlation coefficients 
                                               0.253, 0.504,-0.050, &
                                              -0.006,-0.039, 0.244 ]

real(sp), dimension(3,3), parameter :: cor_a = reshape(corva,[3,3])

real(sp), dimension(9), parameter :: corvb = [ 0.781, 0.000, 0.000, &   !current day correlation coefficients
                                               0.328, 0.637, 0.000, &
                                               0.238,-0.341, 0.873 ]

real(sp), dimension(3,3), parameter :: cor_b = reshape(corvb,[3,3])

!local variables

integer  :: i

real(sp) :: pre    !monthly total precipitation amount (mm)
real(sp) :: wetd   !number of days in month with precipitation (fraction)
real(sp) :: wetf   !fraction of days in month with precipitation (fraction)
real(sp) :: tmn    !minumum temperture (C)
real(sp) :: tmx    !maximum temperture (C)
real(sp) :: cld    !cloud fraction (0=clear sky, 1=overcast) (fraction)

real(sp), pointer :: tmax_mn
real(sp), pointer :: tmin_mn
real(sp), pointer :: cldf_mn
real(sp), pointer :: tmax_cv
real(sp), pointer :: tmin_cv
real(sp), pointer :: cldf_cv

type(randomstate) :: rndst               !integer state of the random number generator
logical       :: pday
real(sp), dimension(3) :: resid     !previous day's weather residuals

real(sp) :: tdiff

real(sp) :: tmean  !rough approximation of daily mean temperature (max + min / 2)
real(sp) :: temp   !rough approximation of daily mean temperature (max + min / 2)
real(sp) :: prec
real(sp) :: tmin
real(sp) :: tmax
real(sp) :: cldf
real(sp) :: es
real(sp) :: tdew
real(sp) :: NI

real(sp) :: pbar     !mean amount of precipitation per wet day (mm)
real(sp) :: pwd      !transition probability of a wet day following a dry day (fraction)
real(sp) :: pww      !transition probability of a wet day following a wet day (fraction)
real(sp) :: alpha    !shape parameter for the precipitation amount function
real(sp) :: beta     !shape parameter for the precipitation amount function
real(sp) :: u        !uniformly distributed random number (0-1)

type(daymetvars), target :: dmetvars

real(sp), dimension(3) :: unorm             !vector of uniformly distributed random numbers (0-1)

!---------------------------------------------------------
!input

pre   = met_in%prec
wetd  = met_in%wetd
wetf  = met_in%wetf
tmn   = met_in%tmin
tmx   = met_in%tmax
cld   = met_in%cldf
rndst = met_in%rndst
pday  = met_in%pday
resid = met_in%resid
NI    = met_in%NI

!shorthand to mean and CV structure

tmin_mn => dmetvars%tmin_mn
tmax_mn => dmetvars%tmax_mn
cldf_mn => dmetvars%cldf_mn
tmin_cv => dmetvars%tmin_cv
tmax_cv => dmetvars%tmax_cv
cldf_cv => dmetvars%cldf_cv

!---------------------------
!1) Precipitation occurrence

!if there is precipitation this month, calculate the precipitation state for today

if (wetf > 0. .and. pre > 0.) then

  !calculate transitional probabilities for dry to wet and wet to wet days
  !Relationships from Geng & Auburn, 1986, Weather simulation models based on summaries of long-term data  

  pwd = 0.75 * wetf
  pww = 0.25 + pwd

  !determine the precipitation state of the current day using the Markov chain approach

  u = ranur(rndst)  !random value from [0,1]

  !the precip status of the current day is conditioned on the status of the previous day    
  if (pday) then   !previous day's precip saved from last call to this subroutine

    if (u - pww > 0.) then
      pday = .false.
    else
      pday = .true.
    end if

  else if (u - pwd > 0.) then
    pday = .false.
  else
    pday = .true.
  end if

  !-----
  
  !2) precipitation amount
  
  if (pday) then  !today is a wet day, calculate the rain amount

    !calculate parameters for the distribution function of precipitation amount

    pbar = pre / wetd

    if (pbar > pmin) then
      beta = -2.16 + 1.83 * pbar
    else
      beta = pbar
    end if

    beta  = max(beta,small)   !put here to avoid infinity values of alpha
    alpha = pbar / beta 

    !today's precipitation

    call ran_gamma(rndst,alpha,beta,.true.,prec)

  else
      
    prec = 0.

  end if

else

  pday = .false.
  prec = 0.

end if


!-----

!3) temperature min and max, cloud fraction

!calculate a baseline mean and cv for today's weather dependent on precip status
!need temp in K here so convert from C.

call meancv(pday,tmn+tfreeze,tmx+tfreeze,cld,dmetvars)

!use random number generator for the normal distribution

do i = 1,3
  call ran_normal(rndst,unorm(i))
end do

!calculate today's residuals for weather variables

resid = matmul(cor_a,resid) + matmul(cor_b,unorm)  !Richardson 1981, eqn 5; WGEN tech report eqn. 3

tmin = tmin_mn * (resid(2) * tmin_cv + 1.)
tmax = tmax_mn * (resid(1) * tmax_cv + 1.)      !WGEN tech report eqn. 13
cldf = cldf_mn * (resid(3) * cldf_cv + 1.)

if (tmin > tmax) then             !set the values equal to the mean between the two  FLAG - this happens too frequently - check with improved parameters
  tdiff = 0.5* (tmin - tmax)
  tmin = tmin - tdiff
  tmax = tmax + tdiff
end if

!ensure cldf results in a value inside its valid range

cldf = min(max(cldf,0.),1.)

!---
!calculate dewpoint
!To estimate dewpoint temperature we use the day's minimum temperature
!this makes the asumption that there is a close correlation between Tmin and dewpoint
!see, e.g., Glassy & Running, Ecological Applications, 1994

if (tmin < 0.) then
  write(0,*) 'tmin weathergenmod ', tmn,tmin_mn,tmin
  stop
end if

es = 0.01 * esat(tmin) !saturation vapor pressure (mbar)  

tdew = 34.07 + 4157. / log(2.1718e8 / es) !Josey et al., Eqn. 10 (K)

!---
!convert calculated temperatures from K to degC

tmin = tmin - Tfreeze
tmax = tmax - Tfreeze
tdew = tdew - Tfreeze

temp = 0.5 * (tmx + tmn)

!---
!Nesterov index, Thonicke et al. (2010) eqn. 5

if (prec <= 3. .and. temp > 0.) then
  NI = NI + tmx * (tmx - (tmn - 4.))
else
  NI = 0.
end if

!write(*,'(a,7f12.2)')'weathergen',wetf,pre,wetd,prec,tmx,tmn,NI

!Nesterov index, Venevsky et al. (2002) eqn. 3

tmean = 0.5 * (tmax + tmin)

!if (prec <= 3. .and. tmean > 0.) then
!  NI = NI + tmax * (tmax - (tmin - 4.))
!else
!  NI = 0.
!end if

!write(*,'(6f12.2)')wetf,pre/wetd,prec,tmean,tmn,NI

!---

met_out%prec  = prec
met_out%tmin  = tmin
met_out%tmax  = tmax
met_out%tdew  = tdew
met_out%cldf  = cldf
met_out%pday  = pday
met_out%rndst = rndst
met_out%resid = resid
met_out%NI    = NI

!10 format(2i4,l4,f8.3,8f9.2,3f9.5)

end subroutine weathergen

!------------------------------------------------------------------------------------------------------------

subroutine meancv(pday,tmn,tmx,cld,dmetvars)

!calculate the mean and CV for a single day value of tmax, tmin, and cloud fraction
!requires temperatures in K

implicit none

!arguments
logical,  intent(in) :: pday  !precipitation status
real(sp), intent(in) :: tmn   !smooth interpolation of monthly minimum temperature (K)
real(sp), intent(in) :: tmx   !smooth interpolation of monthly maximum temperature (K)
real(sp), intent(in) :: cld   !fraction (0-1)
type(daymetvars), intent(out) :: dmetvars

!parameters
!coefficients for the temperature wet day:dry day mean split (based on temperatures in K)

real(sp), dimension(3), parameter :: tmnc = [ -7.424e+01,  4.573e-01, -6.944e-04 ]
real(sp), dimension(3), parameter :: tmxc = [ -3.794e+02,  2.604e+00, -4.440e-03 ]

!coefficients for the temperature mean:CV

real(sp), dimension(3), parameter :: cvtmnc = [ 2.699e-01, -1.362e-03, 1.596e-06 ]
real(sp), dimension(3), parameter :: cvtmxc = [ 2.399e-01, -1.220e-03, 1.517e-06 ] 

real(sp), parameter :: a  =  2.414532    !coefficient for the cloud dry day mean split
real(sp), parameter :: b  =  2.436341    !coefficient for the cloud wet day mean split
real(sp), parameter :: cc = -1.0479262   !slope of the regression line through the cloud wet day mean:CV relationship

!coefficients for the regression line through the cloud dry day:CV relationship

real(sp), dimension(3), parameter :: cd = [ 4.084e+00, 4.165e-03, 1.694e-01 ]

!local variables
real(sp) :: tmaxdiff
real(sp) :: tmindiff

real(sp) :: tmax_mn
real(sp) :: tmin_mn
real(sp) :: cldf_mn
real(sp) :: tmax_cv
real(sp) :: tmin_cv
real(sp) :: cldf_cv

!---

tmax_mn = dmetvars%tmax_mn
tmin_mn = dmetvars%tmin_mn
cldf_mn = dmetvars%tmin_mn
tmax_cv = dmetvars%tmax_cv
tmin_cv = dmetvars%tmin_cv
cldf_cv = dmetvars%cldf_cv

!----------------------------------

tmaxdiff = tmxc(1) + tmxc(2) * tmx + tmxc(3) * tmx**2
tmindiff = tmnc(1) + tmnc(2) * tmn + tmnc(3) * tmn**2

if (pday) then   !wet day
  
  !---tmax---
  tmax_mn = tmx - tmaxdiff                                              !mean
  tmax_cv = cvtmxc(1) + cvtmxc(2) * tmx + cvtmxc(3) * tmx**2            !CV

  !---tmin---
  tmin_mn = tmn - tmindiff                                              !mean
  tmin_cv = cvtmnc(1) + cvtmnc(2) * tmn + cvtmnc(3) * tmn**2            !CV
  
  !---cloud---
  cldf_mn = 1. + (cld - 1.) / (b * cld + 1.)                            !mean
  cldf_cv = cc * cld                                                    !CV

else   !dry day
  
  !---tmax---
  tmax_mn = tmx + tmaxdiff                                              !mean
  tmax_cv = cvtmxc(1) + cvtmxc(2) * tmx + cvtmxc(3) * tmx**2            !CV

  !---tmin---
  tmin_mn = tmn + tmindiff                                              !mean
  tmin_cv = cvtmnc(1) + cvtmnc(2) * tmn + cvtmnc(3) * tmn**2            !CV
  
  !---cloud---
  cldf_mn = - (a * cld) / (cld - a - 1.)                                !mean
  cldf_cv = 1. /(cd(1) * (cld + cd(2))) - cd(3)                         !CV

end if

dmetvars%tmax_mn = tmax_mn
dmetvars%tmin_mn = tmin_mn
dmetvars%cldf_mn = cldf_mn

dmetvars%tmax_cv = tmax_cv
dmetvars%tmin_cv = tmin_cv
dmetvars%cldf_cv = cldf_cv

!write(0,'(a,l4,4f9.3)')'meancv',pday,tmn,tmin_mn,tmin_cv,tmin_mn*tmin_cv
!write(0,'(a,l4,4f9.3)')'meancv',pday,tmx,tmax_mn,tmax_cv,tmax_mn*tmax_cv

end subroutine meancv

!----------------------------------------------------------------------------------------------------------------

real(sp) function esat(temp)
  
  !Function to calculate saturation vapor pressure (Pa) in water and ice
  !From CLM formulation, table 5.2, after Flatau et al. 1992
  
  use parametersmod, only : tfreeze

  implicit none
  
  real(sp), intent(in) :: temp !temperature in K

  real(sp) :: T        !temperature (degC)
  
  real(sp), dimension(9)   :: al !coefficients for liquid water
  real(sp), dimension(9)   :: ai !coefficients for ice
  real(sp), dimension(0:8) :: a  !coefficients
  
  integer :: i
  
  !--------------
  
  al(1) = 6.11213476
  al(2) = 4.44007856e-1
  al(3) = 1.43064234e-2
  al(4) = 2.64461437e-4
  al(5) = 3.05903558e-6
  al(6) = 1.96237241e-8
  al(7) = 8.92344772e-11
  al(8) =-3.73208410e-13
  al(9) = 2.09339997e-16

  ai(1) = 6.11123516
  ai(2) = 5.03109514e-1
  ai(3) = 1.88369801e-2
  ai(4) = 4.20547422e-4
  ai(5) = 6.14396778e-6
  ai(6) = 6.02780717e-8
  ai(7) = 3.87940929e-10
  ai(8) = 1.49436277e-12
  ai(9) = 2.62655803e-15
    
  if (temp <= tfreeze) then   !these coefficients are for temperature values in Celcius
    a(0:8) = ai
  else
    a(0:8) = al
  end if
  
  T = temp - tfreeze
  
  esat = a(0)
  
  do i = 1,8
    esat = esat + a(i) * T**i
  end do
  
  esat = 100. * esat
   
end function esat

!------------------------------------------------------------------------------------------------------------

subroutine rmsmooth(m,dmonth,bcond,r)

!Iterative, mean preserving method to smoothly interpolate mean data to pseudo-sub-timestep values
!From Rymes, M.D. and D.R. Myers, 2001. Solar Energy (71) 4, 225-231

use parametersmod, only : sp,dp

implicit none

!arguments
real(sp), dimension(:), intent(in)  :: m      !vector of mean values at super-time step (e.g., monthly), minimum three values
integer,  dimension(:), intent(in)  :: dmonth !vector of number of intervals for the time step (e.g., days per month)
real(sp), dimension(2), intent(in)  :: bcond  !boundary conditions for the result vector (1=left side, 2=right side)
real(sp), dimension(:), intent(out) :: r      !result vector of values at chosen time step

!parameters
real(sp), parameter :: ot = 1. / 3

!local variables
integer :: n
integer :: ni
integer :: a
integer :: b
integer :: i
integer :: j
integer :: k
integer :: l
integer, dimension(size(r)) :: g
real(sp) :: ck

real(sp), dimension(2) :: bc

!----------

n  = size(m)
ni = size(r)

bc = bcond

!initialize the result vector
i = 1
do a = 1,n
  j = i
  do b = 1,dmonth(a)
    r(i) = m(a)
    g(i) = j
    i = i + 1
  end do
end do

!iteratively smooth and correct the result to preserve the mean

!iteration loop
do i = 1,ni

  do j = 2,ni-1                           
    r(j) = ot * (r(j-1) + r(j) + r(j+1))   !Eqn. 1
  end do

  r(1)  = ot * (bc(1)   + r(1)  +  r(2))   !Eqns. 2
  r(ni) = ot * (r(ni-1) + r(ni) + bc(2))

  j = 1
  do k = 1,n                               !calculate one correction factor per super-timestep
    
    a = g(j)                               !index of the first timestep value of the super-timestep
    b = g(j) + dmonth(k) - 1               !index of the last timestep value of the super-timestep
        
    ck = sum(m(k) - r(a:b)) / ni           !Eqn. 4
    
    do l = 1,dmonth(k)                     !apply the correction to all timestep values in the super-timestep
      r(j) = r(j) + ck
      j = j + 1
    end do

    !correction for circular conditions when using climatology (do not use for transient simulations)
    bc(1) = r(ni)
    bc(2) = r(1)
  
  end do
end do

end subroutine rmsmooth

!------------------------------------------------------------------------------------------------------------

subroutine daily(mval,dval,means)

!linear interpolation of monthly to pseudo-daily values

implicit none

integer, parameter, dimension(12) :: ndaymo = (/  31,28,31,30,31,30,31,31,30,31,30,31 /)
integer, parameter, dimension(14) :: midday = (/ -15,16,44,75,105,136,166,197,228,258,289,319,350,381 /) !middle day of each month

real,    intent(in),  dimension(:)  :: mval
logical, intent(in) :: means

real, intent(out), dimension(:) :: dval

real, dimension(14) :: emval
real, dimension(13) :: slope

integer :: day
integer :: d
integer :: m

!-------------------------------------------------------
!interpolate the monthly values to daily ones in a cyclical format

!copy last month's data to first and vice versa

emval(2:13) = mval
emval(1) = mval(12)
emval(14) = mval(1)

if (means) then

  !calculate slopes

  forall (m = 1:13)
    slope(m) = (emval(m+1) - emval(m)) / (midday(m+1) - midday(m))
  end forall

  !calculate daily values based on monthly means

  m = 1
  do day = 1,365
    if (day > midday(m+1)) m = m + 1
    d = day - midday(m+1)
    dval(day) = slope(m) * d + emval(m+1)
  end do

else
  
  !distribute the total evenly among the days of the month
  
  day = 1
  do m = 1,12
    do d = 1,ndaymo(m)
      dval(day) = mval(m) / ndaymo(m)
      day = day + 1
    end do
  end do

end if

end subroutine daily

!------------------------------------------------------------------------------------------------------------

end module weathergenmod
