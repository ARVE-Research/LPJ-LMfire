module orbitmod

use parametersmod,   only : dp,pir

implicit none

!module subprograms

public :: calcorbitpars
public :: toa_insolation

!module variables

type orbitpars
  real(dp) :: ecc   !eccentricity parameter
  real(dp) :: pre   !precession parameter
  real(dp) :: perh  !longitude of perhelion
  real(dp) :: xob   !obliquity (tilt) (degrees)
end type orbitpars

type insolation
  real(dp) :: ww      !total daily top of the atmosphere insolation (kJ m-2 d-1)
  real(dp) :: dayl    !day length (hours)
  real(dp) :: delta   !solar declination (degrees)
end type insolation

!---------------------------------------------------------------------------------
!parameters

integer, parameter :: nef = 19
integer, parameter :: nob = 47
integer, parameter :: nop = 78

real(dp), parameter :: pirr = pir / 3600.d0

!   1.earth orbital elements : eccentricity           ecc   table 1
!***************************   precessional parameter pre
!                              obliquity              xob   table 2     
!                              general precession     prg
!                              longitude perihelion   perh  table 3

real(dp), dimension(nef), parameter :: ae = [ 0.01860798, &
   0.01627522, -0.01300660, 0.00988829, -0.00336700,  0.00333077, -0.00235400, & 
   0.00140015,  0.00100700, 0.00085700,  0.00064990,  0.00059900,  0.00037800, &
  -0.00033700,  0.00027600, 0.00018200, -0.00017400, -0.00012400,  0.00001250 ]

real(dp), dimension(nef), parameter :: bi = [ 4.2072050, &
   7.3460910,  17.8572630, 17.2205460,  16.8467330,   5.1990790,  18.2310760, &
  26.2167580,   6.3591690, 16.2100160,   3.0651810,  16.5838290,  18.4939800, &
   6.1909530,  18.8677930, 17.4255670,   6.1860010,  18.4174410,   0.6678630 ]

real(dp), dimension(nef), parameter :: ci = [ 28.620089, &
 193.788772,  308.307024, 320.199637,  279.376984,   87.195000,  349.129677, &
 128.443387,  154.143880, 291.269597,  114.860583,  332.092251,  296.414411, &
 145.769910,  337.237063, 152.092288,  126.839891,  210.667199,   72.108838 ]

real(dp), dimension(nef), parameter :: be = bi * pirr
real(dp), dimension(nef), parameter :: ce = ci * pir

real(dp), dimension(nob), parameter :: aob = [ &
  -2462.2214466, -857.3232075, -629.3231835, -414.2804924, -311.7632587, & 
    308.9408604, -162.5533601, -116.1077911,  101.1189923,  -67.6856209, &
     24.9079067,   22.5811241,  -21.1648355,  -15.6549876,   15.3936813, &
     14.6660938,  -11.7273029,   10.2742696,    6.4914588,    5.8539148, &
     -5.4872205,   -5.4290191,    5.1609570,    5.0786314,   -4.0735782, &
      3.7227167,    3.3971932,   -2.8347004,   -2.6550721,   -2.5717867, &
     -2.4712188,    2.4625410,    2.2464112,   -2.0755511,   -1.9713669, &
     -1.8813061,   -1.8468785,    1.8186742,    1.7601888,   -1.5428851, &
      1.4738838,   -1.4593669,    1.4192259,   -1.1818980,    1.1756474, &
     -1.1316126,    1.0896928 ]

real(dp), dimension(nob), parameter :: bib = [ &
  31.609974, 32.620504, 24.172203, 31.983787, 44.828336, 30.973257, &
  43.668246, 32.246691, 30.599444, 42.681324, 43.836462, 47.439436, &
  63.219948, 64.230478,  1.010530,  7.437771, 55.782177,  0.373813, &
  13.218362, 62.583231, 63.593761, 76.438310, 45.815258,  8.448301, &
  56.792707, 49.747842, 12.058272, 75.278220, 65.241008, 64.604291, &
   1.647247,  7.811584, 12.207832, 63.856665, 56.155990, 77.448840, &
   6.801054, 62.209418, 20.656133, 48.344406, 55.145460, 69.000539, &
  11.071350, 74.291298, 11.047742,  0.636717, 12.844549 ]

real(dp), dimension(nob), parameter :: cib = [ &
  251.9025, 280.8325, 128.3057, 292.7252,  15.3747, 263.7951, 308.4258, &
  240.0099, 222.9725, 268.7809, 316.7998, 319.6024, 143.8050, 172.7351, &
   28.9300, 123.5968,  20.2082,  40.8226, 123.4722, 155.6977, 184.6277, &
  267.2772,  55.0196, 152.5268,  49.1382, 204.6609,  56.5233, 200.3284, &
  201.6651, 213.5577,  17.0374, 164.4194,  94.5422, 131.9124,  61.0309, &
  296.2073, 135.4894, 114.8750, 247.0691, 256.6114,  32.1008, 143.6804, &
   16.8784, 160.6835,  27.5932, 348.1074,  82.6496 ]
 
real(dp), dimension(nob), parameter :: bob = bib * pirr
real(dp), dimension(nob), parameter :: cob = cib * pir

real(dp), dimension(nop), parameter :: aop = [ &
  7391.0225890, 2555.1526947,  2022.7629188, -1973.6517951,  1240.2321818,   953.8679112,  -931.7537108, &
   872.3795383,  606.3544732,  -496.0274038,   456.9608039,   346.9462320,  -305.8412902,   249.6173246, &
  -199.1027200,  191.0560889,  -175.2936572,   165.9068833,   161.1285917,   139.7878093,  -133.5228399, &
   117.0673811,  104.6907281,    95.3227476,    86.7824524,    86.0857729,    70.5893698,   -69.9719343, &
   -62.5817473,   61.5450059,   -57.9364011,    57.1899832,   -57.0236109,   -54.2119253,    53.2834147, &
    52.1223575,  -49.0059908,   -48.3118757,   -45.4191685,   -42.2357920,   -34.7971099,    34.4623613, &
   -33.8356643,   33.6689362,   -31.2521586,   -30.8798701,    28.4640769,   -27.1960802,    27.0860736, &
   -26.3437456,   24.7253740,    24.6732126,    24.4272733,    24.0127327,    21.7150294,   -21.5375347, &
    18.1148363,  -16.9603104,   -16.1765215,    15.5567653,    15.4846529,    15.2150632,    14.5047426, &
   -14.3873316,   13.1351419,    12.8776311,    11.9867234,    11.9385578,    11.7030822,    11.6018181, &
   -11.2617293,  -10.4664199,    10.4333970,   -10.2377466,    10.1934446,   -10.1280191,    10.0289441, &
   -10.0034259 ]

real(dp), dimension(nop), parameter :: bip = [ &
  31.609974,  32.620504,  24.172203,   0.636717,  31.983787,   3.138886,  30.973257, &
  44.828336,   0.991874,   0.373813,  43.668246,  32.246691,  30.599444,   2.147012,  10.511172,  42.681324, &
  13.650058,   0.986922,   9.874455,  13.013341,   0.262904,   0.004952,   1.142024,  63.219948,   0.205021, &
   2.151964,  64.230478,  43.836462,  47.439436,   1.384343,   7.437771,  18.829299,   9.500642,   0.431696, &
   1.160090,  55.782177,  12.639528,   1.155138,   0.168216,   1.647247,  10.884985,   5.610937,  12.658184, &
   1.010530,   1.983748,  14.023871,   0.560178,   1.273434,  12.021467,  62.583231,  63.593761,  76.438310, &
   4.280910,  13.218362,  17.818769,   8.359495,  56.792707,   8.448301,   1.978796,   8.863925,   0.186365, &
   8.996212,   6.771027,  45.815258,  12.002811,  75.278220,  65.241008,  18.870667,  22.009553,  64.604291, &
  11.498094,   0.578834,   9.237738,  49.747842,   2.147012,   1.196895,   2.133898,   0.173168 ]

real(dp), dimension(nop), parameter :: cip = [ &
 251.9025, 280.8325, 128.3057, 348.1074, 292.7252, 165.1686, 263.7951,  15.3747, &
  58.5749,  40.8226, 308.4258, 240.0099, 222.9725, 106.5937, 114.5182, 268.7809, 279.6869,  39.6448, 126.4108, &
 291.5795, 307.2848,  18.9300, 273.7596, 143.8050, 191.8927, 125.5237, 172.7351, 316.7998, 319.6024,  69.7526, &
 123.5968, 217.6432,  85.5882, 156.2147,  66.9489,  20.2082, 250.7568,  48.0188,   8.3739,  17.0374, 155.3409, &
  94.1709, 221.1120,  28.9300, 117.1498, 320.5095, 262.3602, 336.2148, 233.0046, 155.6977, 184.6277, 267.2772, &
  78.9281, 123.4722, 188.7132, 180.1364,  49.1382, 152.5268,  98.2198,  97.4808, 221.5376, 168.2438, 161.1199, &
  55.0196, 262.6495, 200.3284, 201.6651, 294.6547,  99.8233, 213.5577, 154.1631, 232.7153, 138.3034, 204.6609, &
 106.5938, 250.4676, 332.3345,  27.3039 ]

real(dp), dimension(nop), parameter :: bop = bip * pirr
real(dp), dimension(nop), parameter :: cop = cip * pir

!------------------------------------

contains

!--------------------------------------------------------------------------------------------


subroutine calcorbitpars(cal_year,orbit)

!-----------------------------------------------------------------------------
!This routine uses the orbital solutions of Berger (1978) and is valid only
!for calculations within +- 1.000.000 yr centered on 1950 AD. 
!For longer periods the Berger (1990) solution should be used.
!(Contact Berger for this 1990 solution).
!
!Recoded by J.O. Kaplan in 2002.
!
!Please refer to :
!  Berger A., 1978. A simple algorithm to compute long term
!                   variations of daily or monthly insolation.
!                   Contr. 18  Inst. of Astronomy and Geophysics,
!                   Universite Catholique de Louvain,
!                   Louvain-la-Neuve, Belgium
!
!  Berger A., 1978. Long term variations of daily insolation and
!                   Quaternary climatic changes.
!                   J. of Atmospheric Sciences 35, 2362-2367
!
!The function value returned by atan is assumed to be a real(dp)
!ranging from -pi/2 to pi/2
!
!Input parameters for the orbital solution are provided in a separate file.
!
!The read and write statements might have to be changed.
!-----------------------------------------------------------------------------

use parametersmod, only : pi,pir,piri

implicit none

!arguments

integer,         intent(in)  :: cal_year
type(orbitpars), intent(out) :: orbit

!parameters

real(dp), parameter :: step = 360.d0 / 365.25d0

real(dp), parameter :: xod = 23.320556d0
real(dp), parameter :: xop =  3.392506d0
real(dp), parameter :: prm = 50.439273d0                                       

!variables

integer :: i

real(dp) :: t
real(dp) :: xes
real(dp) :: xec
real(dp) :: arg
real(dp) :: tra
real(dp) :: rp
real(dp) :: prg

!final calculated variables

real(dp) :: ecc   !eccentricity parameter
real(dp) :: pre   !precession parameter
real(dp) :: perh  !longitude of perhelion
real(dp) :: xob   !obliquity (tilt) (degrees)


!*******************************************                            
!   daily insolation - long term variation *                            
!*******************************************                            
!
!
!This program computes the total daily irradiation received at the top
!of the atmosphere for a given latitude and time in the ka (in kj m-2).

!   3.numerical value for ecc pre xob
!************************************
!   t is negative for the past

t   = real(-cal_year)
xes = 0.d0
xec = 0.d0

do i = 1,nef
  arg = be(i) * t + ce(i)
  xes = xes + ae(i) * sin(arg)
  xec = xec + ae(i) * cos(arg)
end do

ecc = sqrt(xes**2 + xec**2)
tra = abs(xec)                                                     

if (tra > 1.d-8) then

  rp = atan(xes / xec)

  if(xec > 0.d0) then !line 12

    if (xes > 0.d0) then !line 13

      perh = rp * piri

    else if (xes < 0.d0) then !line 14

      rp   = rp + 2.d0 * pi
      perh = rp * piri

    else !line 13

      perh = rp * piri

    end if

  else if (xec < 0.d0) then !line 11

    rp   = rp + pi
    perh = rp * piri

  else !line 10

    if (xes > 0.d0) then !line 17

      rp   = pi / 2.d0
      perh = rp * piri

    else if (xes < 0.d0) then !line 15

      rp   = 1.5d0 * pi
      perh = rp * piri

    else !line 16

      rp   = 0.d0
      perh = rp * piri

    end if
  end if

else

  if (xes > 0.d0) then !line 17

    rp   = pi / 2.d0
    perh = rp * piri

  else if (xes < 0.d0) then !line 15

    rp   = 1.5d0 * pi
    perh = rp * piri

  else !line 16

    rp   = 0.d0
    perh = rp * piri

  end if
end if

prg=prm*t

do i=1,nop

  arg = bop(i) * t + cop(i)
  prg = prg + aop(i) * sin(arg)

end do

prg = prg / 3600.d0 + xop
perh = perh + prg

if (perh > 0.d0) then !line 53
  if (perh > 360.d0) then
    perh = perh - 360.d0
  end if
else if (perh < 0.d0) then
  perh = perh + 360.d0
end if

pre = ecc * sin(perh * pir)

xob = xod

do i = 1,nob
  arg = bob(i) * t + cob(i)
  xob = xob + aob(i) / 3600.d0 * cos(arg)
end do

orbit%ecc  = ecc
orbit%pre  = pre
orbit%perh = perh
orbit%xob  = xob

end subroutine calcorbitpars

!--------------------------------------------------------------------------------------------

subroutine toa_insolation(orbit,nd,phi,ww,dayl,delta)

!calculates top-of-the-atmosphere insolation given day of the year and latitude
!orbital parameters should be precalculated

use parametersmod, only : sp,dp,pi,pir

implicit none

!arguments

type(orbitpars), intent(in)  :: orbit  !orbital parameters for this year, calculated with calcorbitpars
integer,         intent(in)  :: nd     !day of year
real(dp),        intent(in)  :: phi    !latitude (degrees N)

real(sp),        intent(out) :: ww     !top of the atmosphere insolation (kJ m-2 d-1)
real(sp),        intent(out) :: dayl   !day length (h)
real(sp),        intent(out) :: delta  !solar declination (degrees)

!parameters

real(dp), parameter :: ss   = 1366.5d0    !solar constant (W m-2), updated with grab from NASA web site
real(dp), parameter :: tau  =   86.4d0
real(dp), parameter :: test =    1.e-8
real(dp), parameter :: step =  360.d0/365.25d0

!variables

real(dp) :: ecc   !eccentricity parameter
real(dp) :: pre   !precession parameter
real(dp) :: perh  !longitude of perhelion
real(dp) :: xob   !obliquity (tilt) (degrees)

real(dp) :: xec
real(dp) :: dlam

real(dp) :: sf
real(dp) :: so
real(dp) :: xl
real(dp) :: xllp
real(dp) :: xee
real(dp) :: xse
real(dp) :: xlam
real(dp) :: dlamm
real(dp) :: anm
real(dp) :: ranm
real(dp) :: ranv
real(dp) :: anv
real(dp) :: tls

real(dp) :: rphi
real(dp) :: rau
real(dp) :: s
real(dp) :: rlam
real(dp) :: sd
real(dp) :: cd
real(dp) :: rdelta  !declination in radians
real(dp) :: spa
real(dp) :: cp
real(dp) :: aphi    !absolute value of the latitude
real(dp) :: adelta  !absolute value of the solar zenith angle
real(dp) :: tt

real(dp) :: at
real(dp) :: spd
real(dp) :: tp
real(dp) :: stp
real(dp) :: rdayl

!----------------------------------------------------------------------

if (phi < -90.d0 .or. phi > 90.d0) then
  
  write(*,*)'invalid latitude',phi
  stop
  
else if (nd < 1 .or. nd > 366) then
    
  write(*,*)'invalid day of year',nd
  stop
    
end if


ecc  = orbit%ecc
pre  = orbit%pre
perh = orbit%perh
xob  = orbit%xob

!---

sf  = tau * ss / pi
so  = sin(xob * pir)
xl  = perh + 180.d0

xllp = xl * pir
xee  = ecc * ecc
xse  = sqrt(1.d0 - xee)
xlam = (ecc / 2.d0 + ecc * xee / 8.d0) * (1.d0 + xse) * sin(xllp) - xee / 4.d0 * (0.5d0 + xse) * &
        sin(2.d0 * xllp) + ecc * xee / 8.d0 * (1.d0 / 3.d0 + xse) * sin(3.d0 * xllp)

xlam  = 2.d0 * xlam / pir
dlamm = xlam + (nd - 80) * step
anm   = dlamm - xl

ranm = anm * pir
xec  = xee * ecc

ranv = ranm + (2.d0 * ecc - xec / 4.d0) * sin(ranm) + 5.d0 / 4.d0 * ecc**2 * sin(2.d0 * ranm) + 13.d0 / 12.d0 * xec * sin(3.d0 * ranm)

anv  = ranv / pir
tls  = anv + xl

dlam = tls

rphi    =  phi * pir
ranv    =  (dlam - xl) * pir
rau     =  (1.d0 - ecc * ecc) / (1.d0 + ecc * cos(ranv))

s       =  sf / rau / rau
rlam    =  dlam * pir
sd      =  so * sin(rlam)
cd      =  sqrt(1.d0 - sd * sd)

rdelta  =  atan(sd / cd)
delta   =  rdelta / pir
spa     =  sd * sin(rphi)

cp      =  cd * cos(rphi)
aphi    =  abs(phi)
adelta  =  abs(delta)

!singularity for aphi = 90 and delta = 0
!particular cases for lat = 0 or delta = 0

tt = abs(aphi - 90.d0)

if (tt <= test .and. adelta <= test) then

  dayl = 0.00d0
  ww = 0.00d0

else if (adelta <= test) then

  dayl = 12.d0
  ww = s * cos(rphi)

else if (aphi <= test) then

  dayl = 12.d0
  ww = s * cos(rdelta)

else

  at = 90.d0 - adelta
  spd = phi * delta

  if (aphi < at) then

    tp = -spa / cp
    stp = sqrt(1.d0 - tp * tp)
    rdayl = acos(tp)
    dayl = 24.d0 * rdayl / pi
    ww = s * (rdayl * spa + cp * stp)

  else if (spd > 0.d0) then

    dayl = 24.d0
    ww = s * spa * pi

  else if (spd < 0.d0) then

    dayl = 0.00d0
    ww = 0.00d0

  else

    tp =  - spa / cp
    stp = sqrt(1.d0 - tp * tp)
    rdayl = acos(tp)
    dayl = 24.d0 * rdayl / pi
    ww = s * (rdayl * spa + cp * stp)

  end if
end if

end subroutine toa_insolation

!------------------------------------

end module orbitmod
