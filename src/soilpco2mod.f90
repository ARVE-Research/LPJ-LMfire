subroutine soilpco2(co2,zpos,soilprop,mtemp_air,mtemp_soil,mwet_soil,hetresp_mon,soilcconc)

use parametersmod,   only : sp,stdin,stdout
! use mpistatevarsmod, only : soildata

implicit none

! arguments

real(sp),                       intent(in)  :: co2
real(sp), dimension(:)                      :: zpos          ! needed for soil layer depth to midpoint zpos
real(sp),       dimension(:,:),  intent(in)  :: soilprop      ! field capacity and total porosity for each layer (nl,value)
real(sp),       dimension(:),   intent(in)  :: mtemp_air     ! monthly air temperatuer (degC)
real(sp),       dimension(:),   intent(in)  :: mtemp_soil    ! monthly soil temperature (c or K?)
real(sp),       dimension(:),   intent(in)  :: mwet_soil     ! monthly soil moisture as a fraction of field capacity (unitless)
real(sp),       dimension(:,:), intent(in)  :: hetresp_mon   ! monthly pool-specific heterotrophic respiration by layer (gC/m2)
real(sp),       dimension(:),  intent(out) :: soilcconc     ! column averaged CO2 concentration (ppm)
real(sp), dimension(12) :: mw1         !monthly soil layer 1 water content (fraction of available water holding capacity)


! parameters



! local variables
real(sp), allocatable, dimension(:) :: k !  gas conductivity (units)
real(sp) :: surfco2   !  (g m-3)

integer :: nl  !number of soil layers
integer :: m  ! month counter
integer :: l  ! soil layer counter

real(sp), allocatable, dimension(:) :: porespace

real(sp) :: Tsat
real(sp) :: T33

real(sp), allocatable, dimension(:) :: z ! depth (m) midpoint to midpoint Z(1) = atmosphere - is not layer thickness

real(sp), allocatable, dimension(:)    :: a,b,c,d
real(sp)                  :: FG !
real(sp)                  :: B_SHAPE ! shape parameter for porosity
real(sp)                  :: M_SHAPE ! shape parameter for porosity, depends on B_SHAPE

real(sp), allocatable, dimension(:)   :: uptake ! uptake rate (g m^-3 s^-1) at level i we have /z/sec per month
real(sp), allocatable, dimension(:)   :: diffus ! gas diffusivity at layer i
real(sp), allocatable, dimension(:,:)  :: hetresp_mon_sec ! hetresp_mon convertedd to sec----------------------------------------------------------------------
real(sp), allocatable, dimension(:)   :: soilcconc_l ! gas diffusivity at layer i

nl = size(soilprop,dim=1)

allocate(soilcconc_l(nl))

allocate(a(nl))
allocate(b(nl))
allocate(c(nl))
allocate(d(nl))

allocate(k(nl))

allocate(uptake(nl))
allocate(diffus(nl))
allocate(hetresp_mon_sec(12,nl))

allocate(porespace(nl))

allocate(z(nl+2))
z(1:2) = 0. ! atmosphere z(1) and soil surface z(2)
z(3:) = zpos(:) * 2.

k(1) = 0.01 ! typical soil surface-air gas conductivity (units) (typical value, actual value would depend, e.g., on turbulence)

! monthly loop starts here

do m = 1,12

  ! calculate CO2 concentration in g m-3 from atmospheric concentration in ppm

  surfco2 = (co2 * 1e-6 * 1.013e3) * 44.01 / (8.3143 * (273.15 + mtemp_air(m)))
  soilcconc_l(1) = surfco2 ! soil surface co2

  do l = 1,nl

    ! calculate air-filled soil pore space as a function of texture and soil moisture (m3 m-3)

    Tsat = soilprop(l,1)
    T33 = soilprop(l,2)

    porespace(l) = T33 * (1. - mw1(m)) + (Tsat - T33)

    hetresp_mon_sec(m,l) = hetresp_mon(m,l) / z(l) * 30 * 24 * 60 * 60 ! gC m^-2 month^-1 to gC m^-3 s^-1

  end do

  uptake(1) = hetresp_mon_sec(m,1) ! uptake at surface

  do l = 2, nl, 1
          diffus(l) = porespace(l) * B_SHAPE * FG**M_SHAPE
          uptake(l) = hetresp_mon_sec(m,l) * EXP(-z(l) / 0.3) * (z(l+1) - z(l-1)) / 2 ! uptake rate decreases with depth

          if (l < nl) then
              k(l) = diffus(l) / (z(l + 1) - z(l))
          else
              k(l) = 0 ! set to zero at the bottom of soil profile to avoid flow out of bottom
          end if
          a(l+1) = -k(l)
          b(l) = k(l-1) + k(l)
          c(l) = -k(l)
          d(l) = uptake(l)
  end do

  d(2) = d(2) + K(1) * soilcconc_l(1)

  do l = 2, nl-1, 1
      c(l) = c(l) / b(l)
      d(l) = d(l) / b(l)
      b(l+1) = b(l+1) - a(l+1) * c(l)
      d(l+1) = d(l+1) - a(l+1) * d(l)
  end do

  soilcconc_l(nl) = d(nl) / b(nl)

  do l = nl, 2, -1
      soilcconc_l(l) =  d(l) - c(l) * soilcconc_l(l+1)
  end do

  soilcconc_l = sum(soilcconc_l) / nl

  do l = 1, nl, 1
      soilcconc(l) = (soilcconc(l) * 8.3143 * (273.15 + mtemp_soil(m))) / 44.01
  end do


end do

end subroutine soilpco2
