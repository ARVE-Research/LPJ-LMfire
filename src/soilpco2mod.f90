subroutine soilpco2(litter_ag_fast,litter_ag_slow,litter_bg,mw1,tsoil,cpool_surf,cpool_fast,cpool_slow,  &
                     arh,mrh,year,k_fast_ave,k_slow_ave,litter_decom_ave,clay,bulk,spinup,idx, &
                     cflux_surf_atmos,cflux_fast_atmos,cflux_slow_atmos,soil,soilpar,cflux_mean)
! Tsat,T33,T1500,temp,co2,
use parametersmod, only : sp,i8,stdin,stdout
use mpistatevarsmod, only : soildata
use hetrespmod, only : hetresp
use simplesoilmod, only : simplesoil

implicit none
! real(sp) :: co2
! real(sp), dimension(12) :: temp   !mean monthly temperature (degC)
! real(sp) :: dz ! soil layer thickness (m); TODO add from lpjmod
! integer :: nlayers
! integer, parameter :: nl = 3!size(clay)+1 ! incl. atmosphere

! hetrespmod
integer,  intent(in) :: year
logical,  intent(in) :: spinup

real(sp), dimension(:) :: mw1    !soil moisture (relative wetness, volumetric)
real(sp), dimension(:) :: tsoil  !soil temperature
real(sp), dimension(:) :: clay   !soil clay fraction
real(sp), dimension(:) :: bulk   !soil bulk density (kg m-3)

real(sp) :: k_fast_ave
real(sp) :: k_slow_ave

real(sp), dimension(:,:) :: litter_ag_fast
real(sp), dimension(:,:) :: litter_ag_slow
real(sp), dimension(:,:) :: litter_bg

real(sp), dimension(:) :: cpool_surf   !surface SOM pool (2yr turnover time)
real(sp), dimension(:) :: cpool_fast
real(sp), dimension(:) :: cpool_slow
real(sp), dimension(:) :: litter_decom_ave

real(sp), dimension(:) :: arh
real(sp), dimension(:,:) :: mrh
integer(i8), intent(in) :: idx

real(sp) :: cflux_surf_atmos    !soil fast pool decomposition flux to atmosphere
real(sp) :: cflux_fast_atmos    !soil fast pool decomposition flux to atmosphere
real(sp) :: cflux_slow_atmos    !soil slow pool decomposition flux to atmosphere

! simplesoil
type(soildata) :: soil  !state variables sent back out with MPI
real(sp), dimension(:) :: soilpar
real(sp) :: Tsat
real(sp) :: T33
real(sp) :: T1500

! new variables
! integer :: i
real(sp) :: cflux_mean
! real(sp), dimension(nl)   :: b,c,d
! real(sp), dimension(nl)   :: K ! conductance ()
! real(sp), dimension(nl)   :: conc ! atmospheric concentration (ppm -> g m^3)
! real(sp), dimension(nl)   :: uptake ! uptake rate (g m^-3 s^-1) at level i
! real(sp), dimension(nl)   :: DF
! real(sp), dimension(nl)   :: theta ! air filled porosity (m^3 m^-3)
! real(sp), dimension(nl)   :: pp ! Partial pressure ()
! real(sp), dimension(nl+1) :: A ! area (m^2) ?
! real(sp), dimension(nl+1) :: z ! depth (m) midpoint to midpoint Z(1) = atmosphere
! real(sp)                  :: up_surf ! uptake rate at surface based on inverse of CO2 emissions (gC m2 month -> gC m2/X/s) - respiration; TODO can this be uptake(1)?

! real(sp)                  :: FG !
! real(sp)                  :: B_SHAPE ! shape parameter for porosity
! real(sp)                  :: M_SHAPE ! shape parameter for porosity, depends on B_SHAPE

! Start calculations
call hetresp(litter_ag_fast,litter_ag_slow,litter_bg,mw1,tsoil,cpool_surf,cpool_fast,cpool_slow,  &
                     arh,mrh,year,k_fast_ave,k_slow_ave,litter_decom_ave,clay,bulk,spinup,idx, &
                     cflux_surf_atmos,cflux_fast_atmos,cflux_slow_atmos)

call simplesoil(soil,soilpar) !,Tsat,T33,T1500

write(stdout,*)size(clay)
! (volume fraction * Air pressure) * (Molar mass / (Gas Constant * Kelvin T + temp))
! conc(1) = (co2 * 1e-6 * 1.013e3) * 44.01 / (8.3143 * (273.15+temp) !ppm to g m^3; temp needs index?

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
