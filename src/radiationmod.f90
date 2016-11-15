module radiationmod

!consolidated module for calculation airmass, surface downwelling shortwave, net longwave, net radiation and PET.

use parametersmod, only : sp,dp,pi,pir

implicit none

!module subroutines and functions

public  :: initairmass
public  :: calcPjj
public  :: elev_corr
public  :: radpet

private :: airmass
private :: surf_sw
private :: surf_lw
private :: netrad_pet

private :: m
private :: F
private :: esat
private :: desdT

!module parameters

real(sp), parameter :: w   = 15.  !solar angular velocity (degrees hr-1)
real(sp), parameter :: rw  = pir * w !solar angular velocity (radians hr-1)

real(sp), parameter :: m0  =  1.   !air mass at 0 degree solar zenith angle
real(sp), parameter :: m80 =  5.6  !air mass at 80 degree solar zenith angle
real(sp), parameter :: m90 = 39.7  !air mass at 90 degree solar zenith angle

real(sp), parameter :: cos80 = cos(80. * pir)  !(degrees)

real(sp), parameter :: albedo = 0.17  !surface shortwave albedo (fraction)

!module shared variables

real(sp), dimension(3) :: c00  !air mass coefficients for solar zenith angle <=80 degrees
real(sp), dimension(3) :: c80  !air mass coefficients for solar zenith angle  >80 degrees

type airmasspars
  real(sp) :: Ratm    !relative atmospheric pressure 1=sea level
  real(sp) :: mbar    !daytime mean optical air mass (unitless, 1 at equatorial noon)
  real(sp) :: mo      !air mass at cosine zenith angle maximum
  real(sp) :: mc      !air mass at cosine zenith angle medium
  real(sp) :: ml      !air mass at cosine zenith angle bottom quarter range point
end type airmasspars

real(sp) :: lvap    !Latent heat of vaporization of water (temperature dependent) (kJ kg-1)
real(sp) :: gamma   !psychrometer constant (Pa K-1)
real(sp) :: ss      !rate of increase of saturated vapor pressure with temperature (desdT) (Pa K-1)

contains   !the following subroutines and functions

!----------------------------------------------------------------------------------------------------------------

subroutine initairmass()

!calculate parameters used in the airmass calculations

implicit none

c00(1) = 0.008307
c00(2) = (m0 - m80) * (c00(1) + 1.) * (c00(1) + cos80) / (cos80 - 1.)
c00(3) = m0 - c00(2) / (c00(1) + 1.)

c80(1) = 0.037160
c80(2) = (m90 - m80) * c80(1) * (c80(1) + cos80) / cos80
c80(3) = m90 - c80(2) / c80(1)

end subroutine initairmass

!----------------------------------------------------------------------------------------------------------------

subroutine calcPjj(temp,prec,Pjj)

implicit none

!arguments
real(sp), dimension(:), intent(in) :: temp  !annual time series of monthly temperature
real(sp), dimension(:), intent(in) :: prec  !annual time series of monthly precipitation

real(sp), intent(out) :: Pjj

!local variables

integer, dimension(1) :: wm  !index position of the warmest month
integer, dimension(1) :: cm  !index position of the warmest month

real(sp) :: p_wm
real(sp) :: p_cm

!-----------------------------------------------------

!calculation of the "precipitation equitability index"

wm   = maxloc(temp)  !temperature of the warmest month (should be carried as the last warmest month, or running average)
cm   = minloc(temp)  !temperature of the coldest month (should be carried as the last coldest month, or running average)
p_wm = prec(wm(1))   !total precipitation in the warmest month
p_cm = prec(cm(1))   !total precipitation in the coldest month

if (p_wm + p_cm > 0.) then
  Pjj = 2. * (p_wm - p_cm) / (p_wm + p_cm)
  Pjj = max(Pjj,0.)
else
  Pjj = 0.
end if

end subroutine calcPjj

!----------------------------------------------------------------------------------------------------------------

real(sp) function elev_corr(elevation)

implicit none

real(sp), intent(in)  :: elevation 

real(sp), parameter :: z0 = 1. / 8000.

!----

elev_corr = exp(-elevation * z0)

end function elev_corr

!----------------------------------------------------------------------------------------------------------------

subroutine radpet(orbit,lat,tcm,Pjj,dyr,Ratm,met)

use parametersmod, only : midday
use orbitmod,      only : orbitpars,toa_insolation
use weathergenmod, only : metvars_out

implicit none

!arguments

type(orbitpars), intent(in) :: orbit
real(dp),        intent(in) :: lat
real(sp),        intent(in) :: tcm   !temperature of the coldest month
real(sp),        intent(in) :: Pjj
integer,         intent(in) :: dyr
real(sp),        intent(in) :: Ratm  !relative atmospheric pressure (based on elevation)

type(metvars_out), intent(inout) :: met

!local variables

real(sp) :: temp  !daytime mean temperature (C)
real(sp) :: prec  !total precipitation
real(sp) :: cldf  !cloud cover fraction
real(sp) :: tdew  !dewpoint temperature (C)

type(airmasspars) :: air

real(sp) :: toa_sw  !top of the atmosphere downwelling shortwave rad (kJ m-2 d-1)
real(sp) :: delta   !solar declination (degrees)
real(sp) :: pet0    !previous value for PET (mm d-1)
real(sp) :: direct  !direct beam surface downwelling shortwave (kJ m-2 d-1)
real(sp) :: diffuse !diffuse surface downwelling shortwave (kJ m-2 d-1)
real(sp) :: lw_rad  !net longwave (kJ m-2 d-1)

real(sp) :: dayl    !day length (h)
real(sp) :: sw_rad  !total surface downwelling shortwave (kJ m-2 d-1)
real(sp) :: pet     !day potential evapotranspiraton (mm)

!counters

integer :: i

!----------------------------------------------------------------------------------

temp = 0.5 * (met%tmax + met%tmin)
prec = met%prec
cldf = met%cldf

call toa_insolation(orbit,dyr,lat,toa_sw,dayl,delta)

call airmass(lat,delta,dayl,Ratm,air)

call surf_lw(temp,met%tmin,cldf,dayl,lw_rad,tdew)

i = 1

pet  = 0.
pet0 = 0.

do !because of the weak dependence of surface shortwave on PET, we equilibrate PET and surf_sw

  call surf_sw(Pjj,Ratm,toa_sw,cldf,dayl,air,prec,tcm,pet,direct,diffuse)

  sw_rad = direct + diffuse

  call netrad_pet(sw_rad,lw_rad,pet)

  if (abs(pet - pet0) < 0.01 .or. i > 100) exit

  pet0 = pet

  i = i + 1

end do

!write(0,'(7f10.3)')temp,prec,cldf,toa_sw,sw_rad,lw_rad,pet

met%tdew = tdew
met%dayl = dayl
met%srad = sw_rad
met%dpet = pet

!write(0,'(a,3f12.4)')'radpet',sw_rad,lw_rad,pet

end subroutine radpet

!----------------------------------------------------------------------------------------------------------------

subroutine airmass(lat,delta,dayl,Ratm,air)

!This code is based on the paper:
!X. Yin (1997) Optical air mass: Daily integration and its applications, Meteorol. Atmos. Phys. 63, 227-233
!Jed Kaplan, EPFL, 2008

implicit none

!arguments

real(dp), intent(in) :: lat    !latitude (degrees)
real(sp), intent(in) :: delta  !solar declination (degrees)
real(sp), intent(in) :: dayl   !day length (hours)
real(sp), intent(in) :: Ratm   !relative atmospheric pressure

type(airmasspars), intent(out) :: air

!local variables

real(sp) :: mbar    !daytime mean optical air mass (unitless, 1 at equatorial noon)
real(sp) :: mo      !air mass at cosine zenith angle maximum
real(sp) :: mc      !air mass at cosine zenith angle medium
real(sp) :: ml      !air mass at cosine zenith angle bottom quarter range point

real(sp) :: rlat      !latitude (radians)
real(sp) :: rdelta    !solar declination (radians)

real(sp) :: t1        !number of hours between sunrise/sunset and solar noon (hr)
real(sp) :: t80       !solar hour corresponding to the 80 degree zenith angle

real(sp) :: t    !solar hour (hr)

real(sp) :: Z    !solar zenith angle (degrees)
real(sp) :: Zn   !lesser of solar zenith angle at sunset or at midnight (degrees)
real(sp) :: Z0   !zenith angle at solar noon (degrees)
real(sp) :: cosZ !cosine solar zenith angle (fraction), used in calculation of instantaneous air mass

real(sp) :: l

integer :: steps  !integer number of time steps
integer :: i      !counter

real(sp) :: sinlat
real(sp) :: coslat
real(sp) :: sindel
real(sp) :: cosdel

real(sp)               :: a   !values in equation 2.6b
real(sp)               :: b
real(sp), dimension(3) :: c

real(sp) :: tmp1
real(sp) :: tmp2
real(sp) :: tmp3

real(sp) :: tinv

real(sp) :: rZ0
real(sp) :: rZn

real(sp), parameter :: mindayl = 2. * tiny(0._sp)

!-------------------------------------
!calculate daily mean air mass (mbar)

if (dayl == 0.) then

  mbar = m90
  mc   = m90
  ml   = m90
  mo   = m90

else
  
  !basic setup

  rlat   = pir * lat
  rdelta = pir * delta
  
  sinlat = sin(rlat)
  sindel = sin(rdelta)
  coslat = cos(rlat)
  cosdel = cos(rdelta)

  !------

  !Eqn. 2.5 -- commented out because of floating point invalid in rare cases, plus we already have the day length calculated elsewhere (JOK 06.2016)
  
!   if (abs(lat - delta) < 90. .and. abs(lat + delta) >= 90.) then
!    t1 = 12.
!   else
!    t1 = (12. / pi) * acos(-tan(rlat) * tan(rdelta))
!   end if
  
  if (dayl > mindayl) then
    t1 = 0.5 * dayl
  else
    t1 = dayl
  end if
  
  !Eqn. 2.9
  if (abs(lat + delta) >= 90.) then
    Zn = acos(sinlat * sindel - coslat * cosdel) / pir
  else
    Zn = 90.
  end if

  !Eqn. 2.10
  if (abs(lat - delta) >= 90.) then
    Z0 = 90.
  else
    Z0 = lat - delta
  end if
  
  rZ0 = Z0 * pir  !convert to radians
  rZn = Zn * pir
  
  !--------------------------

  b = coslat * cosdel
  
  if (t1 == 0.) then

    mbar = m90
    
  else if (abs(Zn) <= 80.) then
  
    tinv = 1. / t1

    c = c00
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else if (abs(Z0) >= 80.) then
  
    tinv = 1. / t1

    c = c80
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else
    
    t80 = 1. / w * acos((cos80 - sinlat * sindel) / (coslat * cosdel)) / pir  !Eqn. 2.8

    c = c00
    a = c(1) + sinlat * sindel
    
    !write(*,*)'crash',t80,a,b,c
    
    tmp1 = F(t80,a,b,c)

    c = c80
    a = c(1) + sinlat * sindel
    tmp2 = F(t1,a,b,c)

    c = c80
    a = c(1) + sinlat * sindel
    tmp3 = F(t80,a,b,c)

    tinv = 1. / t1

    
    mbar = tinv * (tmp1 + tmp2 - tmp3)
    
  end if

  !--------------------------
  !calculate instantaneous air mass at max, mid, and bottom quarter solar zenith angle (m0, mc, ml)
  
  Z = Z0

  cosZ = cos(Z * pir)

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if
  
  mo = m(cosZ,c)

  !--

  Z = (Z0 + Zn) / 2.
  
  cosz = (cos(rZ0) + cos(rZn)) / 2.

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if

  mc = m(cosZ,c)

  !--

  Z = (Z0 + 3. * Zn) / 4.
  
  cosz = (cos(rZ0) + 3. * cos(rZn)) / 4.

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if

  ml = m(cosZ,c)

end if

!-------------------------------------
!correct calculated air mass for elevation

air%mbar = Ratm * mbar
air%mo   = Ratm * mo
air%mc   = Ratm * mc
air%ml   = Ratm * ml

end subroutine airmass

!----------------------------------------------------------------------------------------------------------------

subroutine surf_sw(Pjj,Ratm,r0,cldf,dayl,air,prec,tcm,pet,direct,diffuse)

!This code is based on the paper:
!X. Yin (1998) Temporally-aggregated atmospheric optical properties as a function of common climatic information:
!Systems development and application, Meteorol. Atmos. Phys. 68, 99-113
!Jed Kaplan, EPFL, 2008

implicit none

!arguments

real(sp), intent(in)  :: Pjj
real(sp), intent(in)  :: Ratm
real(sp), intent(in)  :: r0    !top-of-atmospere insolation (kJ m-2 d-1)
real(sp), intent(in)  :: cldf  !bright sunshine duration fraction, n/N (percent)
real(sp), intent(in)  :: dayl  !daylength (hr)

type(airmasspars), intent(in) :: air
  
real(sp), intent(in)  :: prec  !precipitation mm/day
real(sp), intent(in)  :: tcm   !temperature of the coldest month (used as tropics indicator)
real(sp), intent(in)  :: pet   !potential evapotranspiration mm/day

real(sp), intent(out) :: direct  !direct-beam downwelling shortwave (kJ m-2 d-1)
real(sp), intent(out) :: diffuse !diffuse downwelling shortwave (kJ m-2 d-1)

!local variables

real(sp) :: mbar    !daytime mean optical air mass (unitless, 1 at equatorial noon)
real(sp) :: mo      !air mass at cosine zenith angle maximum
real(sp) :: mc      !air mass at cosine zenith angle medium
real(sp) :: ml      !air mass at cosine zenith angle bottom quarter range point

real(sp) :: tau   !direct insolation atmospheric turbidity factor
real(sp) :: zeta0 !diffuse insolation atmospheric turbidity factor
real(sp) :: x     !tropics indicator (tropical = 1, else 0)
real(sp) :: fm    !atmospheric transmittance function

real(sp) :: j2w
real(sp) :: fdif
real(sp) :: stmp
real(sp) :: sunf   !bright sunshine duration fraction, n/N (fraction)

!-----------------------------------
!parameters

real(sp), parameter :: kp  = 0.500 !links absorption coeff. to trans. coeff.
real(sp), parameter :: kag = 3.300
real(sp), parameter :: kan = 2.320
real(sp), parameter :: kn  = 0.686 !cloud parameter

!----------------------------------------------------------------------------

mbar = air%mbar
mo   = air%mo
mc   = air%mc
ml   = air%ml

!------

sunf = 1. - cldf

if (tcm < 10.) then
  x = 0.
else if (tcm > 20.) then
  x = 1.
else
  x = sin(pi / 2. * (tcm / 10. - 1.))
end if

!Yin Eqn. 4.1
tau = exp(-0.115 * Ratm * ((2.15 - 0.713 * x + exp(-6.74 / (prec + 1.))) * exp(0.0971 * pet) - 0.650 * (1. - x) * Pjj))

fm = 0.01452 * (mbar + ml) * exp(1.403 * tau) - 0.1528 * mo + mc + 0.48700 * (mc - ml) + 0.2323   !Eqn. 2.4 2nd term

direct = sunf * tau**kp * r0 * tau**fm   !Eqn. 2.4

!Yin Eqn. 4.2
zeta0 = 0.503 * exp(-1.20 * Ratm * exp(-0.633 / (prec + 1.) - 0.226 * pet)) * kag**albedo * kan**(1. - sunf) * (1. - kn * (1. - sunf))

diffuse = zeta0 * kag**albedo * kan**(1. - sunf) * (1 - kn * (1. - sunf)) * (tau**kp * r0 - direct)   !Eqn. 2.5

!write(0,'(a,6f12.4)')'shortwave',r0,Ratm,prec,Pjj,sunf,direct+diffuse

end subroutine surf_sw

!----------------------------------------------------------------------------------------------------------------

subroutine surf_lw(temp,tmin,cldf,dayl,lw_rad,tdew)

!This code is based on the paper:
!A. Haxeltine and Prentice, I.C., BIOME3..., Glob. Biogeochem. Cycles, 10, 693-709
!With a new calculations of:
!downwelling longwave   (Josey et al., 2003. J. Geophys. Res., 108(C4), 3108, doi:10.1029/2002JC001418)
!dEsat/dT               (Oleson et al., 2004, CLM 3.0 technical note)
!lvap                   (Henderson-Sellers, 1984. Quart. J. R. Met. Soc., 110, 1186-1190)
!Other references:
!Linacre (1968) Agr. Meteorol., 5, 49-63
!Prentice et al. (1993) Ecological Modelling, 65, 51-70.
!Jed Kaplan, EPFL, 2008, 2011

use parametersmod, only : tfreeze

implicit none

!arguments

real(sp), intent(in)  :: temp    !surface air (2m) temperature (C)
real(sp), intent(in)  :: tmin    !surface air (2m) temperature (C)
real(sp), intent(in)  :: cldf    !cloud cover fraction 
real(sp), intent(in)  :: dayl    !daylength (h)

real(sp), intent(out) :: lw_rad  !daytime net longwave radiation (kJ m-2 d-1)
real(sp), intent(out) :: tdew    !dew point temperature (based on input temperature) (C)

!parameters

real(sp), parameter :: sb = 5.6704e-8  !Stefan-Bolzmann constant (W m-2 K-4)
real(sp), parameter :: e  = 0.98     !emissivity ()
real(sp), parameter :: al = 0.045    !longwave reflectivity (lw albedo), Josey et al., pg 5-9

real(sp), parameter :: a  =  10.77   !parameters in Josey et al.
real(sp), parameter :: b  =   2.34
real(sp), parameter :: c  = -18.44

real(sp), parameter :: cs = 1.5 !shape parameter for the curve relating fractional cloud cover to fractional sunshine duration

!local variables

real(sp) :: Tk     !surface air temperature (K)
real(sp) :: Ts     !ground surface temperature (K)
real(sp) :: TdewK  !dewpoint temperature (K)
real(sp) :: D      !dew point depression (K)
real(sp) :: es     !saturation vapor pressure

real(sp) :: f      !Linacre parameter (function of sunshine fraction)

real(sp) :: Ql     !net longwave radiation (W m-2)
real(sp) :: Ql_up  !upwelling longwave radiation (W m-2)
real(sp) :: Ql_dn  !downwelling longwave radiation (W m-2)

real(sp) :: sunf   !bright sunshine duration fraction, n/N (fraction)

!-------------------------------------------------

sunf = 1. - cldf

Tk = temp + Tfreeze

!calculate gamma, lvap

gamma = 65.05 + temp * 0.064  !psychrometer constant

lvap = 0.001 * 1.91846e6 * (Tk / (Tk - 33.91))**2  !(kJ kg-1) Eqn. from Henderson-Sellers (1984)

ss = desdT(Tk)

f = 0.2 + 0.8 * sunf  !Linacre Eqn. 7

!-------------------------------------------------
!calculate longwave radiation

Ts = Tk !approximation that mean daily surface temperature equals air temp.

!black body upwelling longwave (W m-2)  !various sources e.g., Oleson et al.

Ql_up = e * sb * Ts**4

!--
!Josey formulation for downwelling longwave

!To estimate dewpoint temperature we use the day's minimum temperature
!this makes the asumption that there is a close correlation between Tmin and dewpoint
!see, e.g., Glassy & Running, Ecological Applications, 1994

!es = 0.01 * esat(tmin+tfreeze) !saturation vapor pressure (mbar)  
es = 0.01 * esat(Tk) !saturation vapor pressure (mbar)  

TdewK = 34.07 + 4157. / log(2.1718e8 / es)  !Josey et al., Eqn. 10

D = TdewK - Tk

Ql_dn = sb * (Tk + a*cldf**2 + b*cldf + c + 0.84 * (D + 4.01))**4  !downwelling longwave (W m-2) Josey et al. Eqn. 14,J2

Ql = Ql_up - (1. - al) * Ql_dn   !Josey et al., Eqn 1

!----

lw_rad = 0.001 * 3600. * dayl * Ql  !daytime net longwave (kJ m-2 d-1)

tdew = TdewK - Tfreeze

!write(0,*)'longwave',cldf,Tdewk

end subroutine surf_lw

!----------------------------------------------------------------------------------------------------------------

subroutine netrad_pet(sw_rad,lw_rad,pet)

implicit none

real(sp), intent(in)  :: sw_rad   !downwelling shortwave radiation (kJ m-2 d-1)
real(sp), intent(in)  :: lw_rad   !net longwave radiation (kJ m-2 d-1)
real(sp), intent(out) :: pet      !potential evapotranspiration (mm d-1)

!local variable

real(sp) :: netrad !net radiation (kJ m-2 d-1)

!----

netrad = (1. - albedo) * sw_rad - lw_rad             !(kJ m-2 d-1)

pet = max((ss / (ss + gamma)) * netrad / lvap, 0.)   !(mm d-1)

end subroutine netrad_pet

!----------------------------------------------------------------------------------------------------------------

function esat(temp)
  
  !Function to calculate saturation vapor pressure in water and ice
  !From CLM formulation, table 5.2, after Flatau et al. 1992
  
  use parametersmod, only : dp,tfreeze

  implicit none
  
  real(sp) :: esat  !saturation vapor pressure (Pa)
  real(sp), intent(in) :: temp !temperature in K
  
  real(sp), dimension(9) :: al !coefficients for liquid water
  real(sp), dimension(9) :: ai !coefficients for ice

  real(sp), dimension(0:8) :: a !coefficients
  
  real(sp) :: T
  
  integer :: i
  
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

!----------------------------------------------------------------------------------------------------------------

function desdT(temp)

  !Function to calculate the first derivative of saturation vapor pressure in water and ice vs. temperature
  !From CLM formulation, table 5.3, after Flatau et al. 1992
  
  use parametersmod, only : dp,tfreeze
  
  implicit none
  
  real(sp) :: desdT    !derivative of saturation vapor pressure
  real(sp), intent(in) :: temp !temperature in K
  
  real(sp), dimension(9) :: bl !coefficients for liquid water
  real(sp), dimension(9) :: bi !coefficients for ice

  real(sp), dimension(0:8) :: b !coefficients
  
  real(sp) :: tmp

  real(sp) :: T
  
  integer :: i
  
  bl(1) = 4.44017302e-1
  bl(2) = 2.86064092e-2
  bl(3) = 7.94683137e-4
  bl(4) = 1.21211669e-5
  bl(5) = 1.03354611e-7
  bl(6) = 4.04125005e-10
  bl(7) =-7.88037859e-13
  bl(8) =-1.14596802e-14
  bl(9) = 3.81294516e-17

  bi(1) = 5.03277922e-1
  bi(2) = 3.77289173e-2
  bi(3) = 1.26801703e-3
  bi(4) = 2.49468427e-5
  bi(5) = 3.13703411e-7
  bi(6) = 2.57180651e-9
  bi(7) = 1.32268878e-11
  bi(8) = 3.94116744e-14
  bi(9) = 4.98070196e-17
  
  if (temp <= tfreeze) then
    b(0:8) = bi
  else
    b(0:8) = bl
  end if

  T = temp - tfreeze  !these coefficients are for temperature values in Celcius

  desdT = b(0)

  do i = 1,8
    desdT = desdT + b(i) * T**i
  end do
  
  desdT = 100. * desdT

end function desdT

!----------------------------------------------------------------------------------------------------------------

real(sp) function m(cosZ,c)

!Instantaneous air mass m, equation 2.1 in Yin, 1997

implicit none

real(sp),               intent(in) :: cosZ
real(sp), dimension(:), intent(in) :: c

m = c(2) / (c(1) + cosZ) + c(3)

end function m

!----------------------------------------------------------------------------------------------------------------

real(sp) function F(t1,a,b,c)

!integral air mass function F, equation 2.6b in Yin, 1997
!section inside curly braces only - multiply result by 1/t1 to get mbar

implicit none

real(sp),               intent(in) :: t1
real(sp),               intent(in) :: a
real(sp),               intent(in) :: b
real(sp), dimension(:), intent(in) :: c

real(sp) :: wt1
real(sp) :: wpi

real(sp) :: e1
real(sp) :: e2

wpi  = 180. / (pi * w)
wt1  = rw * t1

if (a > b) then
  
  F = wpi * c(2) / sqrt(a**2 - b**2) * acos((b + a * cos(wt1)) / (a + b * cos(wt1))) + c(3) * t1

else if (a < b) then
  
  e1 = sqrt((b + a) * (1. + cos(wt1))) + sqrt((b - a) * (1. - cos(wt1)))
  e2 = sqrt((b + a) * (1. + cos(wt1))) - sqrt((b - a) * (1. - cos(wt1)))
  
  F = wpi * c(2) / sqrt(b**2 - a**2) * log(e1 / e2) + c(3) * t1

else

  F = wpi * c(2) / a * tan(wt1 / 2.) + c(3) * t1

end if

!write(0,*)'F ab ',a,b
!write(0,*)'F X  ',wpi * c(2) / sqrt(b**2 - a**2) * log(e1 / e2)
!write(0,*)'Fc3t1',c(3) * t1
  
end function F

!----------------------------------------------------------------------------------------------------------------

end module radiationmod
