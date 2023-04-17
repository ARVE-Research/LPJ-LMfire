module parametersmod

use iso_fortran_env, only : int8,int16,int32,int64,real32,real64,input_unit,output_unit,error_unit

implicit none

public :: area

!-----------------------------------------------------------------------------------------------

!type

integer, parameter :: i1 = int8    ! 1 byte integer
integer, parameter :: i2 = int16   ! 2 byte integer
integer, parameter :: i4 = int32   ! 4 byte integer
integer, parameter :: i8 = int64   ! 8 byte integer
integer, parameter :: sp = real32  ! 4 byte real
integer, parameter :: dp = real64  ! 8 byte real

!size

integer, parameter :: bytes_sp = sizeof(0._sp)
integer, parameter :: bytes_dp = sizeof(0._dp)

!parameters
 
integer, parameter :: stdin  = input_unit
integer, parameter :: stdout = output_unit
integer, parameter :: stderr = error_unit

integer, parameter :: npft     =  9 
integer, parameter :: npftpar  = 51
integer, parameter :: nsoilpar =  7
integer, parameter :: ncvar    =  3
integer, parameter :: climbuf  = 20
integer, parameter :: nspec    =  6   !trace gas species from biomass burning
integer, parameter :: nhclass  =  7   !number of height classes for discretization

real(sp), parameter :: allom1 = 100.
real(sp), parameter :: allom2 =  40.
real(sp), parameter :: allom3 =   1.  ! 0.5
real(sp), parameter :: allom4 =   0.3

real(sp), parameter :: reinickerp = 1.6
real(sp), parameter :: latosa     = 8.e3
real(sp), parameter :: wooddens   = 2.e5  !(g C m-2)  NB this value would actually be typical for total dry mass in some hardwoods, should probably PFT specific

integer, parameter :: maxoutvars = 100   !maximum number of variables for output

!code for land use type: (1) unmanaged, (2) rainfed crop, (3) regrowing natural vegetation 

integer, parameter,  dimension(3) :: lutype = [ 1,2,3 ]  !could eventually be more categories

!variables

integer :: ntiles

type pftflag
  logical :: evergreen   !whether PFT is evergreen
  logical :: needle      !whether PFT is needleleaved (alternative: broadleaved)
  logical :: boreal      !whether PFT is boreal
  logical :: raingreen   !whether PFT is raingreen
  logical :: summergreen !whether PFT is summergreen
  logical :: tree        !whether PFT is a tree (alternative: grass)
end type pftflag

type(pftflag), dimension(npft) :: pft

real(sp), dimension(npft) :: sla          !PFT specific leaf area (m2/gC)

!type sapling_pars
  real(sp), dimension(npft,ncvar) :: lm_sapl  !initial (sapling) leaf mass (gC/m2)
  real(sp), dimension(npft,ncvar) :: hm_sapl  !initial (sapling) heartwood mass (gC/m2)
  real(sp), dimension(npft,ncvar) :: sm_sapl  !initial (sapling) sapwood mass (gC/m2)
  real(sp), dimension(npft,ncvar) :: rm_sapl  !initial (sapling) fine root mass (gC/m2)
!end type sapling_pars

!type(sapling_pars), save, target, dimension(npft) :: sapl

real(sp), save, target, dimension(npft,npftpar) :: pftpar   !PFT parameters

integer, parameter, dimension(12) :: firstday  = [ 1, 32,60, 91,121,152,182,213,244,274,305,335 ]  !day number of the first day of each month
integer, parameter, dimension(12) :: midday    = [ 16,44,75,105,136,166,197,228,258,289,319,350 ]  !day number of mid-month day
integer, parameter, dimension(12) :: ndaymonth = [ 31,28,31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]  !number of days in each month

!-----------------------------------------------------------------------------------------------
!parameters used in the soil physics and surface energy balance routines

integer, parameter :: nl     = 11         !number of soil layers
integer, parameter :: snomax = 5          !maximum number of snow layers
integer, parameter :: ns     = 1 - snomax !index of top snow layer

real(sp), parameter :: dt    = 86400.d0  !number of seconds in 24hrs

!fixed CLM parameters

real(dp), parameter :: pi    = 3.14159265358979323846d0 !26433 83279 50288 41971 69399 37510 (unitless)
real(dp), parameter :: grav  = 9.80616d0     !Gravitational acceleration  (m s-2)
real(dp), parameter :: Pstd  = 101325.d0     !Standard pressure           (Pa)
real(dp), parameter :: sigsb = 5.67e-8       !Stefan-Boltzmann constant   (W m-2 K-4)
real(dp), parameter :: kbz   = 1.38065e-23   !Boltzmann constant          (J K-1 molecule-1)
real(dp), parameter :: avo   = 6.0221415e26  !Avogadro's number           (molecule kmol-1)
real(dp), parameter :: Rgas  = avo * kbz     !Universal gas constant      (J K-1 kmol-1)
real(dp), parameter :: MWda  = 28.966d0      !Molecular weight of dry air (kg kmol-1)
real(dp), parameter :: Rda   = Rgas / MWda   !Dry air gas constant        (J K-1 kg-1)
real(dp), parameter :: MWwv  = 18.016d0      !Mol. weight of water vapor  (kg kmol-1)
real(dp), parameter :: Rwv   = Rgas / MWwv   !Water vapor gas constant    (J K-1 kg-1)
real(dp), parameter :: vkarm = 0.4d0         !von Karman constant (unitless)

real(sp), parameter :: Tfreeze = 273.15      !freezing temperature of freshwater (K)

real(dp), parameter :: pliq  = 1000.d0     !density of water (kg m-3)
real(dp), parameter :: pice  =  917.d0     !density of ice (kg m-3)

real(dp), parameter :: Cair   =     1.00464e3  !Heat capacity of dry air (J kg-1 K-1) (CLM parameterization)
real(dp), parameter :: Cliq   =     4.18800e3  !Heat capacity of liquid water  (J kg-1 K-1)
real(dp), parameter :: Cice   =     2.11727e3  !Heat capacity of ice (typical) (J kg-1 K-1)

real(dp), parameter :: lvap  = 2.501e6       !Latent heat of vaporization (J kg-1)
real(dp), parameter :: Lf    = 3.337e5       !Water latent heat of fusion (J kg-1)
real(dp), parameter :: lsub  = lvap + Lf     !Latent heat of sublimation  (J kg-1)

real(dp), parameter :: Kliq   =     0.57d0   !Thermal conductivity of liquid water (typical) (W m-1 K-1)
real(dp), parameter :: Kice   =     2.29d0   !Thermal conductivity of ice (typical) (W m-1 K-1)
real(dp), parameter :: Kair   =     0.0243d0 !Thermal conductivity of dry air (W m-1 K-1)

!derived CLM parameters

real(dp), parameter :: Timp = 0.05d0  !volumetric liquid water content at which the soil is impermeable (fraction)
real(dp), parameter :: Pmax =-1.e8    !maximum allowed value for soil matric potential

real(dp), parameter :: Esoil  = 0.96d0 !thermal emissivity of bare soil (fraction)
real(dp), parameter :: Ewater = 0.96d0 !thermal emissivity of water (fraction)
real(dp), parameter :: Esnow  = 0.97d0 !thermal emissivity of snow (fraction)

real(dp), parameter :: Tstune = 0.34d0    !Tuning factor to turn first layer T into surface T (CLM parameterization, pg 88).
real(dp), parameter :: CNfac  = 0.5d0     !Crank Nicholson factor between 0 and 1

real(dp), parameter :: pdrysnow =  50.00000d0 !density of dry snow falling at under -15 deg C (kg m-3)
real(dp), parameter :: pwetsnow = 169.15775d0 !density of wet snow falling at over 2 deg C (kg m-3) (50.d0 + 1.7d0 * 17.d0**1.5, CLM eqn 7.18)
real(dp), parameter :: Tc = 2.5d0 !critical threshold temperature separating rain from snow (deg C)

real(dp), parameter :: lb = 1.e-5  !baseflow parameter (mm s-1) CLM eqn. 7.117
real(dp), parameter :: kd = 0.04d0 !saturated soil hydraulic cond. contributing to baseflow (mm s-1) CLM eqn. 7.118

real(dp), parameter :: wpond = 10.d0 !max. quantity of water for ponding (kg m-2)

real(dp), parameter :: Smin  = 0.033d0  !irreducible water saturation of snow (fraction)

real(dp), parameter :: z0mg  = 0.01d0   !momentum roughness length for soil (m, CLM eqn 3.49)

!other parameters

real(dp), parameter :: porg   =  1300.d0  !density of soil organic matter (typical) (kg m-3)
real(dp), parameter :: pmine  =  2660.d0  !density of mineral soil (typical) (kg m-3)

real(dp), parameter :: Korg   =     0.25d0   !Thermal conductivity of soil organic matter (typical) (W m-1 K-1)
real(dp), parameter :: Corg   =     2.496e6  !Heat capacity of soil organic matter (typical) (J m-3 K-1)

real(dp), parameter :: Kmine  =     5.85d0   !Thermal conductivity of soil minerals (typical) (W m-1 K-1) NB average bet. quartz and other minerals, from Hillel
real(dp), parameter :: Cmine  =     2.000e6  !Heat capacity of soil minerals (typical) (J m-3 K-1)

real(dp), parameter :: pir     = pi / 180.0d0
real(dp), parameter :: piri    = 180.d0 / pi

contains

!-----------------------------------------------------------------------------------------------

real(sp) function area(lat,minutes)

  !this function returns the size of a regular grid cell in square meters.

  implicit none

  real(dp), intent(in) :: lat
  real(sp), intent(in), dimension(2) :: minutes

  real(dp), parameter :: pi      =    3.14159265359d0
  real(dp), parameter :: radius  = 6378.137d0 !km, WGS-84 spherical approximation
  real(dp), parameter :: deg2rad = pi / 180.d0

  real(dp) :: cellarea
  real(dp) :: deltalat
  real(dp) :: deltalon
  real(dp) :: elevation
  real(dp), dimension(2) :: resolution

  resolution = minutes / 60.d0

  elevation = deg2rad * (lat + 0.5d0 * resolution(2))

  deltalat = deg2rad * resolution(1)
  deltalon = deg2rad * resolution(2)

  cellarea = 2.d0 * radius**2 * deltalon * cos(elevation) * sin(0.5d0 * deltalat)

  area = real(cellarea * 1.e6)

end function area

!-----------------------------------------------------------------------------------------------

!note from wikipedia: air conductivity 0.024, snow (typical) 0.11, soil (typical) 0.17-1.13

!1.9e6 J m-3 K-1 = 0.45 cal cm-3 degC-1
!3.2217e-6 W m-1 K-1 = 7.694902e-9 cal cm-1 sec-1 degC-1

end module parametersmod
