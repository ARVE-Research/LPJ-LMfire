subroutine soilpco2(co2,soilpar,mtemp_air,mtemp_soil,mwet_soil,hetresp_mon,soilcconc)

use parametersmod,   only : sp,stdin,stdout

implicit none

! arguments

real(sp),                       intent(in)  :: co2
real(sp),       dimension(:,:)  intent(in)  :: soilprop      ! field capacity and total porosity for each layer (nl,value)
real(sp),       dimension(:),   intent(in)  :: mtemp_air     ! monthly air temperatuer (degC)
real(sp),       dimension(:),   intent(in)  :: mtemp_soil    ! monthly soil temperature (c or K?)
real(sp),       dimension(:),   intent(in)  :: mwet_soil     ! monthly soil moisture as a fraction of field capacity (unitless)
real(sp),       dimension(:,:), intent(in)  :: hetresp_mon   ! monthly, and by pool
real(sp),                       intent(out) :: soilcconc     ! column averaged CO2 concentration (ppm)

! parameters

real(sp) :: ksoil_air =  XX ! typical soil surface-air gas conductivity (units) (typical value, actual value would depend, e.g., on turbulence)

! local variables

real(sp) :: surfco2   !  (g m-3)

integer :: nl  !number of soil layers
integer :: m  ! month counter
integer :: l  ! soil layer counter

real(sp), allocatable, dimension(:) :: porespace

real(sp) :: Tsat
real(sp) :: T33



! ----------------------------------------------------------------------

nl = size(soilprop,dim=1)

allocate(porespace(nl))

! monthly loop starts here

do m = 1,12

  ! calculate CO2 concentration in g m-3 from atmospheric concentration in ppm

  surfco2 = (co2 * 1e-6 * 1.013e3) * 44.01 / (8.3143 * (273.15 + mtemp_air(m))

  do l = 1,nl

    ! calculate air-filled soil pore space as a function of texture and soil moisture (m3 m-3)
  
    Tsat = soilprop(l,1)
    T33 = soilprop(l,2)
  
    porespace(l) = T33 * (1. - mw1(m)) + (Tsat - T33)
    
  end do


  
end do







! z(2) = 0. ! soil surface
! K(1) = XXX ! TODO calc conductance from hetres output
! nlayers = nl
! up_surf = 1 - XXX ! TODO calc surface uptake as inverse of respiration from hetres output

! theta = XXX !TODO calc from total porosity; Theta - T1500 wilting point; T33 field capacity

! do i = 2, nlayers, 1
!         z(i+1) = z(i) + dz
!         DF(i) = theta*B_SHAPE*FG**M_SHAPE
!         U(i) = up_surf*EXP(-Z(i)/.3)*(z(i+1)-z(i-1))/2
!         if ( i < nlayers ) then
!             K(i) = DF(i)/(z(i+1)-z(i))
!         else
!             K(i) = 0 ! set to zero at the bottom of soil profile to avoid flow out of bottom
!         end if
!         A(i+1) = -K(i)
!         b(i) = K(i-1) + K(i)
!         c(i) = -K(i)
!         d(i) = uptake(i)
! end do
!
! d(2) = d(2) + K(1) * conc(1)
!
! do i = 2, nlayers-1, 1
!     c(i) = c(i) / b(i)
!     d(i) = d(i) / b(i)
!     b(i+1) = b(i+1) - A(i+1) * c(i)
!     d(i+1) = d(i+1) - A(i+1) * d(i)
! end do
!
! conc(nlayers) = d(nlayers) / b(nlayers)
!
! do i = nlayers, 2, -1
!     conc(i) =  d(i) - c(i) * conc(i+1)
! end do
!
! do i = 1, nlayers, 1
!     pp(i) = (conc(i) * 8.3143 * (273.15+tsoil(i))) / 44.01
! end do


end subroutine soilpco2
