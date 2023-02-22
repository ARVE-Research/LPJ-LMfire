module photosynthesismod

implicit none

contains

! SUBROUTINE PHOTOSYNTHESIS
! Adapted from Farquhar (1982) photosynthesis model, as simplified by
! Collatz et al 1991, Collatz et al 1992, and Haxeltine & Prentice 1996

subroutine photosynthesis(ca,temp,fpar,par,dayl,c4,lambda,rd,agd,adtmm,x1,x2,x3,x4,pfti)

use parametersmod, only : sp,dp,npft,npftpar,ncvar,pftpar

implicit none

! parameters
  
real(sp), parameter :: alphaa    =    0.5     ! fraction of PAR assimilated at ecosystem level
real(sp), parameter :: alphac3   =    0.08    ! intrinsic quantum efficiency of CO2 uptake in C3 plants
real(sp), parameter :: alphac4   =    0.053   ! C4 intrinsic quantum efficiency
real(sp), parameter :: bc3       =    0.015   ! leaf respiration as fraction of Vmax for C3 plants
real(sp), parameter :: bc4       =    0.02    ! leaf respiration as fraction of vmax for C4 plants
real(sp), parameter :: cmass     =   12.      ! atomic mass of carbon
real(sp), parameter :: cq        =    4.6e-6  ! conversion factor for solar radiation at 550 nm from J/m2 to E/m2 (E=mol quanta)
real(sp), parameter :: e0        =  308.56    ! parameter in Arrhenius temp response function
real(sp), parameter :: kc25      =   30.      ! value of kc at 25 deg C
real(sp), parameter :: ko25      =    3.e4    ! value of ko at 25 deg C
real(sp), parameter :: lambdamc3 =    0.8     ! optimal (maximum?) ci:ca ratio for C3 plants
real(sp), parameter :: lambdamc4 =    0.4     ! optimal ci:ca ratio for c4 plants
real(sp), parameter :: m         =   25.      ! corresponds to parameter p in Eqn 28, Haxeltine & Prentice 1996
real(sp), parameter :: n0        =    7.15    ! leaf N concentration (mg/g) not involved in photosynthesis
real(sp), parameter :: p         =    1.e5    ! atmospheric pressure (Pa)
real(sp), parameter :: po2       =   20.9e3   ! O2 partial pressure (Pa)
real(sp), parameter :: q10kc     =    2.1     ! q10 for temperature-sensitive parameter kc
real(sp), parameter :: q10ko     =    1.2     ! q10 for temperature-sensitive parameter ko
real(sp), parameter :: q10tau    =    0.57    ! q10 for temperature-sensitive parameter tau
real(sp), parameter :: t0c3      =  250.      ! base temperature (K) in Arrhenius temperature response function for C3 plants
real(sp), parameter :: t0c4      =  260.      ! base temperature in Arrhenius func for C4 plants
real(sp), parameter :: tau25     = 2600.      ! value of tau at 25 deg C
real(sp), parameter :: theta     =    0.7     ! colimitation (shape) parameter
real(sp), parameter :: tk25      =  298.15    ! 25 deg C in Kelvin
real(sp), parameter :: tmc3      =   45.      ! maximum temperature for C3 photosynthesis
real(sp), parameter :: tmc4      =   55.      ! maximum temperature for C4 photosynthesis

! arguments

logical,  intent(in)  :: c4
integer,  intent(in)  :: pfti
real(sp), intent(in)  :: ca
real(sp), intent(in)  :: dayl
real(sp), intent(in)  :: fpar
real(sp), intent(in)  :: lambda
real(sp), intent(in)  :: par
real(sp), intent(in)  :: temp
real(sp), intent(in)  :: x1
real(sp), intent(in)  :: x2
real(sp), intent(in)  :: x3
real(sp), intent(in)  :: x4
real(sp), intent(out) :: adtmm
real(sp), intent(out) :: agd
real(sp), intent(out) :: rd

! local variables

real(sp) :: adt
real(sp) :: and
real(sp) :: apar
real(sp) :: b
real(sp) :: c1
real(sp) :: c2
real(sp) :: gammastar
real(sp) :: high
real(sp) :: jc
real(sp) :: je
real(sp) :: k1
real(sp) :: k2
real(sp) :: k3
real(sp) :: kc
real(sp) :: ko
real(sp) :: lambdam
real(sp) :: low
real(sp) :: pa
real(sp) :: phipi
real(sp) :: pi
real(sp) :: s
real(sp) :: sigma
real(sp) :: t0
real(sp) :: tau
real(sp) :: tstress
real(sp) :: vm

! ----------------------------------------------------------------------------------------------------

lambdam = pftpar(pfti,26) 

! Return without performing calculations if daylength too short

if (dayl < 0.01) then

  agd   = 0.
  adtmm = 0.
  rd    = 0.

  return

end if

! APAR in J/m2/day
! alphaa = scaling factor for absorbed PAR at ecosystem, versus leaf, scale
! See Eqn 4, Haxeltine & Prentice 1996

apar = par * fpar * alphaa

! calculate temperate inhibition function

if (temp < x4) then

   k1  = 2. * alog((1. / 0.99) - 1.) / (x1 - x2)
   k2  = (x1 + x2) / 2.
   low = 1./(1.+ exp(k1*(k2-temp)))

   k3   = alog(0.99 / 0.01) / (x4 - x3)
   high = 1. - 0.01 * exp(k3 * (temp - x3))

   tstress = low * high

else

   tstress = 0.

end if

!  write(0,*)'tstress',temp,x4,tstress
if (tstress < 1.e-2) tstress = 0.

! First calculate catalytic capacity of rubisco, Vm, assuming optimal (non-water-stressed) value for lambda, i.e. lambdamc3

if (c4) then  !C4 photosynthesis

  ! Specify C1C4, C2C4; Eqns 14,15, Haxeltine & Prentice 1996
  ! Notes:
  ! - alphaa, the upscaling parameter accounting for the reduction in PAR utilisation in ecosystems compared with leaf level,
  !   appears in the calculation of APAR instead of here
  ! - Cmass, the atomic weight of carbon, used in unit conversion from molC to g appears in the calculation of Vm instead of here
  ! - parameter phipi is not needed for calculation of optimal Vm which assumes optimal intercellular CO2 concentration (lambdamc4)

  c1 = tstress * alphac4

  ! High-temperature inhibition modelled conservatively as a step function
  ! prohibiting photosynthesis above 55 deg C (Table 3.7, Larcher 1983)

  if (temp > tmc4) c1 = 0.

  c2 = 1.

  b  = bc4    !Choose C4 value of b for Eqn 10, Haxeltine & Prentice 1996
  t0 = t0c4   !base temperature for temperature response of rubisco

else  ! C3 photosynthesis

  ! Temperature-adjusted values of kinetic parameters; Eqn 22, Haxeltine & Prentice 1996a

  ko  = ko25  *  q10ko**((temp - 25.) / 10.)  !Michaelis constant of rubisco for O2
  kc  = kc25  *  q10kc**((temp - 25.) / 10.)  !Michaelis constant for CO2
  tau = tau25 * q10tau**((temp - 25.) / 10.)  !CO2/O2 specificity ratio

  ! CO2 compensation point (CO2 partial pressure, Pa) Eqn 8, Haxeltine & Prentice 1996

  gammastar = po2 / (2. * tau)
 
  ! Convert ambient CO2 level, ca, from mole fraction to partial pressure (Pa)

  !write(0,*) 'CO2 in photosynthesis: ', ca * 1e6
  !stop

  pa = ca * p

  ! Non-water-stressed intercellular CO2 partial pressure (Pa); Eqn 7, Haxeltine & Prentice 1996

  pi = lambdam * pa

  ! Calculation of C1C3, Eqn 4, Haxeltine & Prentice 1996
  ! Notes: - there is an error in this equation in the above paper (missing 2.0* in denominator),
  ! which is fixed here (see Eqn A2, Collatz et al 1991)

  ! - There is no longer an explicit temperature inhibition function 
  !   (low-temperature inhibition is now done mechanistically by imposing a temperature-dependent upper limit on Vm, see below)
  ! - There is no longer any reduction in maximum photosynthesis due to leaf age (phic)
  ! - alphaa, the upscaling parameter accounting for the reduction in PAR utilisation in ecosystems compared with leaf level, 
  !   appears in the calculation of APAR instead of here
  ! - Cmass, the atomic weight of carbon, used in unit conversion from molC to g appears in the calculation of Vm instead of here
 
  c1 = tstress * alphac3 * ((pi - gammastar)/(pi + 2. * gammastar))

  ! High temperature inhibition modelled primarily by suppression of LUE by decreased relative affinity of rubisco for CO2
  ! relative to O2 with increasing temperature, but we also implement a step function to prohibit any C3 photosynthesis 
  ! above 45 degrees (Table 3.7, Larcher 1983)

  if (temp  > tmc3) c1 = 0.

  ! Calculation of C2C3, Eqn 6, Haxeltine & Prentice 1996

  c2 = (pi - gammastar) / (pi + kc * (1. + po2 / ko))

  b  = bc3   !Choose C3 value of b for Eqn 10, Haxeltine & Prentice 1996
  t0 = t0c3  !base temperature for temperature response of rubisco

end if

! Eqn 13, Haxeltine & Prentice 1996

s = (24. / dayl) * b

! Eqn 12, Haxeltine & Prentice 1996

sigma = sqrt(max(0., 1. - (c2 - s)/(c2 - theta * s)))

! Calculation of optimal rubisco capacity, Vm, in gC/m2/day; Eqn 11, Haxeltine & Prentice 1996

vm = (1. / b) * (c1 / c2) * ((2. * theta - 1.) * s - (2. * theta * s - c2) * sigma) * apar * cmass * cq

! Now use this Vm value to calculate actual photosynthesis

if (.not.c4) then  !C3 photosynthesis

  ! Intercellular CO2 partial pressure in Pa
  ! Eqn 7, Haxeltine & Prentice 1996

  pi = lambda * pa

  ! Recalculation of C1C3, C2C3 with actual pi

  c1 = tstress * alphac3 * ((pi - gammastar) / (pi + 2. * gammastar))

  if (temp > tmc3) c1 = 0.  !high-temperature inhibition

  c2 = (pi - gammastar) / (pi + kc * (1. + po2 / ko))

else  !C4 photosynthesis

  ! Parameter accounting for effect of reduced intercellular CO2 concentration on photosynthesis, Phipi
  ! Eqn 14,16, Haxeltine & Prentice 1996; Fig 1b, Collatz et al 1992
 
  phipi = min(lambda / lambdam, 1.)

  c1 = tstress * phipi * alphac4

  if (temp> tmc4) c1 = 0.  !high-temperature inhibition

end if

! je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3; Convert je from daytime to hourly basis

! Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h ; Eqn 3, Haxeltine & Prentice 1996

je = c1 * apar * cmass * cq / dayl

! Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h; Eqn 5, Haxeltine & Prentice 1996

jc = c2 * vm / 24.

!write(0,'(i5,4f12.2,4f10.4)')pfti,temp,par,fpar,dayl,je,jc,tstress,c1

if (je < 1.e-10 .or. jc <= 1.e-10) then

  agd=0.0

else

  !  write(0,*)'daily',je,jc,(je+jc)**2.0,4.0*theta*je*jc

  ! Calculation of daily gross photosynthesis, Agd, gC/m2/day; Eqn 2, Haxeltine & Prentice 1996
  ! Note: - there is an error in this equation in the above paper (missing theta in 4*theta*je*jc term) which is fixed here

  agd = (je + jc - sqrt((je + jc)**2. - 4. * theta * je * jc)) / (2. * theta) * dayl

end if

! Daily leaf respiration, Rd, gC/m2/day; Eqn 10, Haxeltine & Prentice 1996

rd = b * vm

! Daily net photosynthesis (at leaf level), And, gC/m2/day

and = agd - rd

! Total daytime net photosynthesis, Adt, gC/m2/day; Eqn 19, Haxeltine & Prentice 1996

adt = and + (1. - dayl / 24.) * rd

! Convert adt from gC/m2/day to mm/m2/day using ideal gas equation

adtmm = adt / cmass * 8.314 * (temp + 273.3) / p * 1000.

end subroutine photosynthesis

! ----------------------------------------------------------------------------------------------------
! scrap for unused variable declarations and other code from f77 version

!subroutine photosynthesis(ca,temp,fpar,par,dayl,c4,sla,nmax,lambda,rd,agd,adtmm,x1,x2,x3,x4,pfti)

! real(sp), intent(in) :: nmax
! real(sp), intent(in) :: sla

! integer  :: pft
! integer :: i
! real(sp) :: cn
! real(sp) :: tk
! real(sp) :: vmmax

end module photosynthesismod
