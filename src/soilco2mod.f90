module soilco2mod

use parametersmod, only : stdout,stderr

implicit none

public :: soilco2

contains

!-------------------------------------------------------------------------------------------------------------
subroutine soilco2(co2,soil,soilprop,mtemp_air,mtemp_soil,mwet_soil,hetresp_mon,soilcconc_dec,soilco2conc)

use parametersmod,   only : sp,stdin,stdout !,npftpar,nsoilpar,ncvar,i8
use mpistatevarsmod, only : soildata

implicit none

! arguments

real(sp),                       intent(in)  :: co2
type(soildata), intent(in)                  :: soil          ! needed for soil layer depth to midpoint zpos
real(sp),       dimension(:,:), intent(in)  :: soilprop      ! field capacity and total porosity for each layer (nl,value)
real(sp),       dimension(:),   intent(in)  :: mtemp_air     ! monthly air temperature (degC)
real(sp),       dimension(:),   intent(in)  :: mtemp_soil    ! monthly soil temperature (degC)
real(sp),       dimension(:),   intent(in)  :: mwet_soil     ! monthly soil moisture as a fraction of field capacity (unitless)
real(sp),       dimension(:,:), intent(in)  :: hetresp_mon   ! monthly pool-specific heterotrophic respiration (gC/m2); cflux_litter_atmos; cflux_surf_atmos; cflux_fast_atmos; cflux_slow_atmos
real(sp),       dimension(:),   intent(out) :: soilco2conc   ! monthly soil column CO2 concentration (ppm)
real(sp),       dimension(:),   intent(inout) :: soilcconc_dec   ! December CO2 concentration (mgCO2 m^-3)



! parameters
real(sp), parameter :: B_SHAPE = 4.547 ! shape parameter for porosity 0.9
real(sp), parameter :: M_SHAPE = 2.3 ! shape parameter for porosity, depends on B
real(sp), parameter :: Dg0st = 0.0000139 ! Gas diffusivity standard temperature and pressure (m2 s-1)
real(sp), parameter :: T0 = 273.15
real(sp), parameter :: dt = 24. ! time steps (hours)
real(sp), parameter :: dz = 0.025 ! depth increment (m)

! local variables
real(sp) :: surfco2   !  (g m-3)
real(sp) :: Dg0_max

integer :: nz  !number of soil layers
integer :: nl  !number of discretized soil layers
integer :: nl1  !number of discretized soil layers in top soil layer
integer :: nl2  !number of discretized soil layers in bottom soil layer
integer :: m  ! month counter
integer :: l  ! soil layer counter
integer :: Ndt ! number of Euler time steps
integer :: tt
real(sp) :: dt_i ! dt discretized to Ndt


real(sp) :: Tsat
real(sp), dimension(2) :: T33
real(sp) :: logDg0 ! log of diffusivity at surface

real(sp), allocatable, dimension(:) :: logDgs ! log of diffusivity
real(sp), allocatable, dimension(:) :: Dgs ! Diffusivity
real(sp), allocatable, dimension(:) :: porespace ! porosity
real(sp), allocatable, dimension(:) :: z ! layer thickness (m)
real(sp), allocatable, dimension(:,:) :: hetresp_mon_layers ! hetresp_mon discretized into layers
real(sp), allocatable, dimension(:) :: soilcconc     ! CO2 concentration per layer (mg CO2 m^3)
real(sp), allocatable, dimension(:) :: soilcconc_old ! soilcconc at previous timestep
real(sp), allocatable, dimension(:) :: soilcconc_nz ! soilcconc aggregated to nz


! ----------------------------------------------------------------------
Dg0_max = 0.000025 ! from Dg0st at 50C and 755 hPa
Dg0_max = Dg0_max * 3600. ! Convert gas diffusivity (m2 s-1) to (m2 hr-1)

nz = size(soil%zpos,dim=1)

allocate(z(nz+1)) ! layers incl. atmosphere/surface layer
allocate(porespace(nz))
allocate(soilcconc_nz(nz+1))


z(1) = 0. ! layer thickness (m)
z(2) = (soil%zpos(1) * 2.) / 100.
z(3) = ((soil%zpos(2) - soil%zpos(1)) * 2.) / 100.

! number of layers incl. atmosphere/surface layer excl. zero-bottom layer
nl1 = ceiling(z(2) / dz)
nl2 = ceiling(z(3) / dz)
nl = 1 + nl1 + nl2

allocate(hetresp_mon_layers(12, nl+1))
allocate(logDgs(nl+1))
allocate(Dgs(nl+1))
allocate(soilcconc(nl+1))
allocate(soilcconc_old(nl+1))

! Compute timesteps Ndt = (dt*Dgs)/((dz^2)*0.5); 0.05 as safety buffer
Ndt = ceiling((dt * Dg0_max) / ((dz**2) * (0.5 - 0.05)))
dt_i = dt / Ndt ! n hours discretized in nx timesteps

! monthly loop starts here
do m = 1,12
  ! calculate CO2 concentration in mg CO2 m-3 from atmospheric concentration in ppm
  surfco2 = (co2 * 44.01 / 1000.) * 101325. / (8.3143 * (T0 + mtemp_air(m)))

  ! separate monthly respiration into top and bottom layer
  ! surface flux into top layer; fast and slow flux into bottom layer
  ! convert monthly respiration (gC m^-2 month^-1) to volumetric respiration (mg CO2 m^-3 h^-1)
  hetresp_mon_layers(m,1) = 0.
  do l=2, nl, 1
      if (l .le. nl2) then
          hetresp_mon_layers(m,l) = (((hetresp_mon(m,1) + hetresp_mon(m,2)) * 1000) / (30. * 24.)) / nl1
      else
          hetresp_mon_layers(m,l) = (((hetresp_mon(m,1) + hetresp_mon(m,2)) * 1000) / (30. * 24.)) / nl2
      end if
  end do

  ! Diffusion coefficient (m2 s-1) to (m2 hr-1)
  ! Dgs(1) = Dg0st * ((T0 + mtemp_soil(m)) / T0)**1.75
  logDg0 = log(Dg0st) + (1.75 * log((T0 + mtemp_soil(m)) / T0))
  Dgs(1) = exp(logDg0) * 3600

  ! Porosity
  do l = 1, nz, 1
      ! calculate air-filled soil pore space as a function of texture and soil moisture (m3 m-3)
      Tsat = soilprop(l, 1)
      T33(l) = soilprop(l, 2)

      porespace(l) = T33(l) * (1. - mwet_soil(m)) + (Tsat - T33(l))
  end do

  ! Moldrup et al (2000) 
  do l=2, nl+1, 1
      if (l .le. nl2) then
          !Dgs(l) = Dgs(1) * B_SHAPE * porespace(1)**M_SHAPE
          logDgs(l) = logDg0 + log((2. * T33(1)**3.) + (0.04 * T33(1))) + ((2. + (3. / B_SHAPE)) * log(porespace(1) / T33(1)))
      else
          ! Dgs(l) = Dgs(1) * B_SHAPE * porespace(2)**M_SHAPE
          logDgs(l) = logDg0 + log((2. * T33(2)**3.) + (0.04 * T33(2))) + ((2. + (3. / B_SHAPE)) * log(porespace(2) / T33(2)))
      end if
	  Dgs(l) = exp(logDgs(l)) * 3600.
      ! Dgs(l) = 0.0137
  end do

  ! Compute soil CO2 concentrations following Ryan et al 2015
  ! https://doi.org/10.5194/gmd-11-1909-2018
  ! Initial conditions of soilcconc at January year 0 need to be set in lpjmod.f90;
  ! after that take December from previous year
  soilcconc(1) = surfco2

  if (m .eq. 1) then
      !soilcconc(2:(nl1 + 1)) = soilcconc_dec(2) * nl1
      !soilcconc((nl1 + 2):nl) = soilcconc_dec(3) * nl2
      soilcconc = soilcconc_dec
  end if

  do tt = 1, Ndt, 1
      soilcconc_old = soilcconc
      do l = 2, nl, 1
          if (l .eq. nl) then
              soilcconc(l) = soilcconc_old(l) + dt_i * ((Dgs(l) / dz**2) * (soilcconc_old(l-1) - soilcconc_old(l)) + hetresp_mon_layers(m, l))
          else
              soilcconc(l) = soilcconc_old(l) + dt_i * ((Dgs(l) / dz**2) * (soilcconc_old(l+1) - 2 * soilcconc_old(l) + soilcconc_old(l-1)) + ((Dgs(l+1) - Dgs(l-1)) * (soilcconc_old(l+1) - soilcconc_old(l-1))) / (4 * dz**2) + hetresp_mon_layers(m, l))
          end if
      end do
  end do

  ! aggregate soilcconc to soilcco2onc layers
  !soilcconc_nz(1) = soilcconc(1)
  !soilcconc_nz(2) = sum(soilcconc(2:(nl1 + 1))) / nl1
  !soilcconc_nz(3) = sum(soilcconc((nl1 + 2):nl)) / nl2

  ! convert mgCO2 m^-3 to ppm CO2
  soilco2conc(m) = (sum(soilcconc(2:nl)) * 8.3143 * (T0 + mtemp_soil(m)) / (44.01 / 1000. * 101325.)) / nl
end do

! December soilcconc for next year
soilcconc_dec = soilcconc

! write(stdout,*)soilco2conc

end subroutine soilco2
!-------------------------------------------------------------------------------------------------------------

end module soilco2mod
