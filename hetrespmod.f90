module hetrespmod

implicit none

public :: littersom2
public :: hetresp

contains

!-------------------------------------------------------------------------------------------------------------

subroutine littersom2(litter_ag_fast,litter_ag_slow,litter_bg,mw1,tsoil,cpool_fast,cpool_slow,  &
                     arh,mrh,year,k_fast_ave,k_slow_ave,litter_decom_ave,clay,bulk,spinup)

use parametersmod, only : sp,npft,npftpar,nsoilpar,ncvar

implicit none

!arguments

integer,  intent(in) :: year
logical,  intent(in) :: spinup

real(sp), dimension(:), intent(in)    :: mw1    !soil moisture (relative wetness, volumetric)
real(sp), dimension(:), intent(in)    :: tsoil  !soil temperature
real(sp), dimension(:), intent(inout) :: clay   !soil clay fraction
real(sp), dimension(:), intent(inout) :: bulk   !soil bulk density (kg m-3)

real(sp), intent(inout) :: k_fast_ave
real(sp), intent(inout) :: k_slow_ave

real(sp), dimension(:,:), intent(inout) :: litter_ag_fast
real(sp), dimension(:,:), intent(inout) :: litter_ag_slow
real(sp), dimension(:,:), intent(inout) :: litter_bg

real(sp), dimension(:), intent(inout) :: cpool_fast
real(sp), dimension(:), intent(inout) :: cpool_slow
real(sp), dimension(:), intent(inout) :: litter_decom_ave

real(sp), dimension(:),   intent(out) :: arh
real(sp), dimension(:,:), intent(out) :: mrh

!parameters

integer,  parameter :: soil_equil_year = 750      !number of years until pool sizes for soil decomposition solved analytically

real(sp), parameter :: syr = real(soil_equil_year)

real(sp), parameter :: k_litter_fast10 = 1./2.     !fast litter decomposition rate at 10 deg C (yr-1)
real(sp), parameter :: k_litter_slow10 = 1./20.    !slow litter decomposition rate at 10 deg C (yr-1)
real(sp), parameter :: k_soil_fast10   = 1./20.    !fast soil decomposition rate at 10 deg C (yr-1)
real(sp), parameter :: k_soil_slow10   = 1./1000.  !slow soil decomposition rate at 10 deg C (yr-1)

real(sp), parameter :: atmfrac         = 0.7      !fraction of litter decomposition going directly into the atmosphere
real(sp), parameter :: soilfrac        = 1. - atmfrac   !fraction of litter decomposition going to soil C pools

real(sp), dimension(7), parameter :: dz = [ 0.2, 0.2, 0.2, 0.2, 0.2, 1.0, 1.0 ]  !soil layer thicknesses (only for this subroutine)

real(sp), parameter :: sb  = 0.18  !parameters for the moisture effect on respiration equation below
real(sp), parameter :: sfc = 0.7   !units relative soil wetness S

real(sp), parameter :: t0  = 0.15  !maximum SOM partitioning coefficient between fast and slow soil pools

!local variables

integer  :: pft
integer  :: m
integer  :: c

real(sp), parameter :: dt = 1./12.  !timestep

real(sp), dimension(12) :: temp_resp      !monthly temperature response of decomposition
real(sp), dimension(12) :: moist_resp     !monthly moisture response of decomposition

real(sp), dimension(12) :: k_litter_fast  !monthly fast litter decomposition rate (/month)
real(sp), dimension(12) :: k_litter_slow  !monthly slow litter decomposition rate (/month)
real(sp), dimension(12) :: k_soil_fast    !monthly fast pool decomposition rate (/month)
real(sp), dimension(12) :: k_soil_slow    !monthly slow pool decomposition rate (/month)

real(sp), dimension(12) :: litter_decom   !total litter decomposition (sum fast, slow and bg)

real(sp), dimension(12) :: ek_lf
real(sp), dimension(12) :: ek_ls
real(sp), dimension(12) :: ek_sf
real(sp), dimension(12) :: ek_ss

real(sp), dimension(3,12) :: respflux

real(sp), dimension(npft) :: litterdag_fast  !above-ground component of fast litter decomp
real(sp), dimension(npft) :: litterdag_slow  !above-ground component of slow litter decomp
real(sp), dimension(npft) :: litterdag_bg    !below-ground component of litter decomposition

real(sp) :: cflux_litter_soil   !litter decomposition flux to soil
real(sp) :: cflux_litter_atmos  !litter decomposition flux to atmosphere
real(sp) :: cflux_fast_atmos    !soil fast pool decomposition flux to atmosphere
real(sp) :: cflux_slow_atmos    !soil slow pool decomposition flux to atmosphere

real(sp) :: fastfrac  !fraction of litter entering fast soil decomposition pool (opposed to slowC)
real(sp) :: slowfrac  !fraction of litter entering slow soil decomposition pool

real(sp) :: claytot   !total mass of clay in the soil column (to 3m, kg)
real(sp) :: soilmass  !total mass of the soil column (depth * bulk density) (to 3m, kg)

!---------------------------------------------------------------
!calculation of organic matter decomposition in litter and soil,
!including transfer or organic matter from litter to soil

!original LPJ derivation of formula
!  (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate

!from (1),
!  (2) c = c0*exp(-kt) where c0=initial pool size

!from (2), decomposition in any month given by
!  (3) delta_c = c0 - c0*exp(-k)

!from (3),
! (4) delta_c = c0*(1.-exp(-k))

!Temperature response function is a modified Q10 relationship (Lloyd & Taylor 1994)

where (tsoil > -25.)
  temp_resp = exp(308.56 * ((1. / 56.02) - (1. / (tsoil + 273.15 - 227.13))))
elsewhere
  temp_resp = 0.
end where

!moisture response based on soil layer 1 moisture content (eqn. 7, Manzoni & Porporato, Soil Biol. Biochem. 39, 2007)
!this does not work in lpj because mw1 represents only the range between wilting point and field capacity, not total soil wetness

!where (mw1 <= sb)
!  moist_resp = 0.
!elsewhere (mw1 <= sfc)
!  moist_resp = (mw1 - sb) / (sfc - sb)
!elsewhere
!  moist_resp = sfc / mw1
!end where

!original LPJ function after Foley et al. (1995)

moist_resp = 0.25 + (0.75 * mw1)

!write(0,*)mw1
!read(*,*)

!Calculate decomposition rates as a function of annual rate, temperature, and moisture

k_litter_fast = k_litter_fast10 * dt * temp_resp * moist_resp
k_litter_slow = k_litter_slow10 * dt * temp_resp * moist_resp
k_soil_fast   = k_soil_fast10   * dt * temp_resp * moist_resp
k_soil_slow   = k_soil_slow10   * dt !* temp_resp * moist_resp  !experiment to take away climate response

!write(0,*)'hetresp',k_soil_fast,temp_resp,moist_resp

ek_lf = 1. - exp(-k_litter_fast)  !monthly vectors
ek_ls = 1. - exp(-k_litter_slow)
ek_sf = 1. - exp(-k_soil_fast)
ek_ss = 1. - exp(-k_soil_slow)

!calculate partitioning fractions for fast and slow SOM based on soil clay content

!calculate total mass of clay in soil
!this values is summed over the layers present
!if properties for the 5th soil layer are provided (80-100 cm layer)
!clay mass (and therefore carbon storage potential) is calculated down to 3m (following Jobbagy & Jackson)

bulk = max(0.,bulk)
clay = max(0.,clay)

claytot  = sum(0.2 * clay * bulk) + 1. * clay(5) * bulk(5) + 1. * clay(5) * bulk(5)  !final units kg
soilmass = sum(0.2 * bulk) + 2. * bulk(5)

!the partitioning between fast and slow SOM pools is a function of total soil clay mass
!after Jobbagy & Jackson (Ecol. App. 10, 2000)

slowfrac = t0 * claytot / 300. !300 kg m-2 is ca. the maximum value in the ISRIC WISE data, extrapolating layer 5 to 3m
fastfrac = 1. - slowfrac

!write(0,*)bulk
!write(0,*)clay
!write(0,*)claytot,soilmass
!write(0,*)fastfrac,slowfrac
!read(*,*)

mrh = 0.

do m = 1,12

  !litter decomposition: pft vectors

  do pft = 1,npft

    litterdag_fast(pft) = ek_lf(m) * litter_ag_fast(pft,1)
    litterdag_slow(pft) = ek_ls(m) * litter_ag_slow(pft,1)
    litterdag_bg(pft)   = ek_lf(m) * litter_bg(pft,1)        !belowground litter has the same turnover rate as fast litter

    !reduction of litter pool sizes

    litter_ag_fast(pft,1) = max(litter_ag_fast(pft,1) - litterdag_fast(pft),0.)
    litter_ag_slow(pft,1) = max(litter_ag_slow(pft,1) - litterdag_slow(pft),0.)
    litter_bg(pft,1)      = max(litter_bg(pft,1)      - litterdag_bg(pft),0.)

  end do

  !sum up decomposition over all pfts and all pools
  
  litter_decom(m) = sum(litterdag_fast + litterdag_slow + litterdag_bg)

  !write(0,*)'hetresp',mw1
  !write(0,*)tsoil
  write(0,'(9f8.2)')litter_ag_fast(:,1) !,litterdag_fast
  !write(0,'(9f8.2)')litter_ag_slow(:,1) !,litterdag_slow
  !write(0,'(9f8.2)')litter_bg(:,1) !,litterdag_bg
  !write(0,'(i5,f8.2)')m,litter_decom(m)
  !write(0,*)'end hetresp'
  !read(*,*)

  !partition decomposing litter flux between atmosphere and SOM components

  cflux_litter_atmos = atmfrac  * litter_decom(m)
  cflux_litter_soil  = soilfrac * litter_decom(m)

  !further subdivide soil fraction between fast and slow soil pools

  cpool_fast(1) = cpool_fast(1) + cflux_litter_soil * fastfrac
  cpool_slow(1) = cpool_slow(1) + cflux_litter_soil * slowfrac

  !calculate SOM decomposition
  
  cflux_fast_atmos = cpool_fast(1) * ek_sf(m) !eqn 4
  cflux_slow_atmos = cpool_slow(1) * ek_ss(m) !eqn 4

  mrh(m,1) = cflux_litter_atmos + cflux_fast_atmos + cflux_slow_atmos
  
  respflux(1,m) = cflux_litter_atmos
  respflux(2,m) = cflux_fast_atmos
  respflux(3,m) = cflux_slow_atmos

  cpool_fast(1) = max(cpool_fast(1) - cflux_fast_atmos,0.)
  cpool_slow(1) = max(cpool_slow(1) - cflux_slow_atmos,0.)
  
end do

arh = sum(mrh,dim=1)

!write(0,'(4f10.4)')arh(1),sum(respflux(1,:)),sum(respflux(2,:)),sum(respflux(3,:))
!read(*,*)


!write(0,*)sum(litter_decom),arh(1)
!read(*,*)

!SOIL DECOMPOSITION EQUILIBRIUM CALCULATION
!Analytical solution of differential flux equations for fast and slow
!soil carbon pools.  Implemented after (soil_equil_year) simulation
!years, when annual litter inputs should be close to equilibrium.  Assumes
!average climate (temperature and soil moisture) from all years up to
!soil_equil_year.

if (spinup .and. year == soil_equil_year + 1 .and. k_fast_ave > 0. .and. k_slow_ave > 0.) then

  !Analytically calculate pool sizes this year only

  !Rate of change of soil pool size = litter input - decomposition
  !  (5) dc/dt = litter_decom - kc

  !At equilibrium,
  !  (6) dc/dt = 0

  !From (5) & (6),
  !  (7) c = litter_decom / k

  cpool_fast(1) = (soilfrac * fastfrac * litter_decom_ave(1)) / k_fast_ave   !eqn 7

  cpool_slow(1) = (soilfrac * slowfrac * litter_decom_ave(1)) / k_slow_ave   !eqn 7
     
else if (year <= soil_equil_year) then

  !Update running average respiration rates and litter input

  k_fast_ave = k_fast_ave + sum(k_soil_fast) / syr
  k_slow_ave = k_slow_ave + sum(k_soil_slow) / syr
  
  litter_decom_ave(1) = litter_decom_ave(1) + sum(litter_decom) / syr

end if

end subroutine littersom2

!in standard LPJ:
!litter has three types (pools) per pft: aboveground fast, aboveground slow, belowground
!aboveground fast and belowground have the same turnover time
!soil has two pools: fast, slow

!-------------------------------------------------------------------------------------------------------------

subroutine hetresp(litter_ag_fast,litter_ag_slow,litter_bg,mw1,tsoil,cpool_surf,cpool_fast,cpool_slow,  &
                     arh,mrh,year,k_fast_ave,k_slow_ave,litter_decom_ave,clay,bulk,spinup,idx)

use parametersmod, only : sp,npft,npftpar,nsoilpar,ncvar,i8

implicit none

!arguments

integer,  intent(in) :: year
logical,  intent(in) :: spinup

real(sp), dimension(:), intent(in)    :: mw1    !soil moisture (relative wetness, volumetric)
real(sp), dimension(:), intent(in)    :: tsoil  !soil temperature
real(sp), dimension(:), intent(inout) :: clay   !soil clay fraction
real(sp), dimension(:), intent(inout) :: bulk   !soil bulk density (kg m-3)

real(sp), intent(inout) :: k_fast_ave
real(sp), intent(inout) :: k_slow_ave

real(sp), dimension(:,:), intent(inout) :: litter_ag_fast
real(sp), dimension(:,:), intent(inout) :: litter_ag_slow
real(sp), dimension(:,:), intent(inout) :: litter_bg

real(sp), dimension(:), intent(inout) :: cpool_surf
real(sp), dimension(:), intent(inout) :: cpool_fast
real(sp), dimension(:), intent(inout) :: cpool_slow
real(sp), dimension(:), intent(inout) :: litter_decom_ave

real(sp), dimension(:),   intent(out) :: arh
real(sp), dimension(:,:), intent(out) :: mrh

integer(i8), intent(in) :: idx

!parameters

integer,  parameter :: soil_equil_year = 750      !number of years until pool sizes for soil decomposition solved analytically

real(sp), parameter :: syr = real(soil_equil_year)

real(sp), parameter :: k_litter_fast10 = 1./2.     !fast litter decomposition rate at 10 deg C (yr-1)
real(sp), parameter :: k_litter_slow10 = 1./20.    !slow litter decomposition rate at 10 deg C (yr-1)
real(sp), parameter :: k_soil_fast10   = 1./20.    !fast soil decomposition rate at 10 deg C (yr-1)
real(sp), parameter :: k_soil_slow10   = 1./1000.  !slow soil decomposition rate at 10 deg C (yr-1)

real(sp), parameter :: atmfrac         = 0.7      !fraction of litter decomposition going directly into the atmosphere
real(sp), parameter :: soilfrac        = 1. - atmfrac   !fraction of litter decomposition going to soil C pools

real(sp), dimension(7), parameter :: dz = [ 0.2, 0.2, 0.2, 0.2, 0.2, 1.0, 1.0 ]  !soil layer thicknesses (only for this subroutine)

real(sp), parameter :: sb  = 0.18  !parameters for the moisture effect on respiration equation below
real(sp), parameter :: sfc = 0.7   !units relative soil wetness S

real(sp), parameter :: t0  = 0.15  !maximum SOM partitioning coefficient between fast and slow soil pools

real(sp), parameter :: klit2som = 1. - exp(-2.)  !transfer rate for fast litter to surface SOM pool

!local variables

integer  :: pft
integer  :: m
integer  :: c

real(sp), parameter :: dt = 1./12.  !timestep

real(sp), dimension(12) :: temp_resp      !monthly temperature response of decomposition
real(sp), dimension(12) :: moist_resp     !monthly moisture response of decomposition

real(sp), dimension(12) :: k_litter_fast  !monthly fast litter decomposition rate (/month)
real(sp), dimension(12) :: k_litter_slow  !monthly slow litter decomposition rate (/month)
real(sp), dimension(12) :: k_soil_fast    !monthly fast pool decomposition rate (/month)
real(sp), dimension(12) :: k_soil_slow    !monthly slow pool decomposition rate (/month)

real(sp), dimension(12) :: litter_decom   !total litter decomposition (sum fast, slow and bg)

real(sp), dimension(12) :: ek_lf
real(sp), dimension(12) :: ek_ls
real(sp), dimension(12) :: ek_sf
real(sp), dimension(12) :: ek_ss

real(sp), dimension(3,12) :: respflux

real(sp), dimension(npft) :: litterdag_fast  !above-ground component of fast litter decomp
real(sp), dimension(npft) :: litterdag_slow  !above-ground component of slow litter decomp
real(sp), dimension(npft) :: litterdag_bg    !below-ground component of litter decomposition
real(sp), dimension(npft) :: lit2som         !amount of C transferred from fast litter to surface SOM pool

real(sp) :: cflux_litter_soil   !litter decomposition flux to soil
real(sp) :: cflux_litter_atmos  !litter decomposition flux to atmosphere
real(sp) :: cflux_surf_atmos    !soil fast pool decomposition flux to atmosphere
real(sp) :: cflux_fast_atmos    !soil fast pool decomposition flux to atmosphere
real(sp) :: cflux_slow_atmos    !soil slow pool decomposition flux to atmosphere

real(sp) :: fastfrac  !fraction of litter entering fast soil decomposition pool (opposed to slowC)
real(sp) :: slowfrac  !fraction of litter entering slow soil decomposition pool

real(sp) :: claytot   !total mass of clay in the soil column (to 3m, kg)
real(sp) :: soilmass  !total mass of the soil column (depth * bulk density) (to 3m, kg)

!---------------------------------------------------------------
!calculation of organic matter decomposition in litter and soil,
!including transfer or organic matter from litter to soil

!original LPJ derivation of formula
!  (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate

!from (1),
!  (2) c = c0*exp(-kt) where c0=initial pool size

!from (2), decomposition in any month given by
!  (3) delta_c = c0 - c0*exp(-k)

!from (3),
! (4) delta_c = c0*(1.-exp(-k))

!Temperature response function is a modified Q10 relationship (Lloyd & Taylor 1994)

where (tsoil > -25.)
  temp_resp = exp(308.56 * ((1. / 56.02) - (1. / (tsoil + 273.15 - 227.13))))
elsewhere
  temp_resp = 0.
end where

!moisture response based on soil layer 1 moisture content (eqn. 7, Manzoni & Porporato, Soil Biol. Biochem. 39, 2007)
!this does not work in lpj because mw1 represents only the range between wilting point and field capacity, not total soil wetness

!where (mw1 <= sb)
!  moist_resp = 0.
!elsewhere (mw1 <= sfc)
!  moist_resp = (mw1 - sb) / (sfc - sb)
!elsewhere
!  moist_resp = sfc / mw1
!end where

!original LPJ function after Foley et al. (1995)

moist_resp = 0.25 + (0.75 * mw1)

!write(0,*)mw1
!read(*,*)

!Calculate decomposition rates as a function of annual rate, temperature, and moisture

k_litter_fast = k_litter_fast10 * dt * temp_resp * moist_resp
k_litter_slow = k_litter_slow10 * dt * temp_resp * moist_resp
k_soil_fast   = k_soil_fast10   * dt * temp_resp * moist_resp
k_soil_slow   = k_soil_slow10   * dt !* temp_resp * moist_resp  !experiment to take away climate response

!if (idx == 1) then
!  do m = 1,12
!    write(*,*)m,temp_resp(m),moist_resp(m)
!  end do
!  write(*,*)
!end if


ek_lf = 1. - exp(-k_litter_fast)  !monthly vectors
ek_ls = 1. - exp(-k_litter_slow)
ek_sf = 1. - exp(-k_soil_fast)
ek_ss = 1. - exp(-k_soil_slow)

!calculate partitioning fractions for fast and slow SOM based on soil clay content

!calculate total mass of clay in soil
!this values is summed over the layers present
!if properties for the 5th soil layer are provided (80-100 cm layer)
!clay mass (and therefore carbon storage potential) is calculated down to 3m (following Jobbagy & Jackson)

bulk = max(0.,bulk)
clay = max(0.,clay)

claytot  = sum(0.2 * clay * bulk) + 1. * clay(5) * bulk(5) + 1. * clay(5) * bulk(5)  !final units kg
soilmass = sum(0.2 * bulk) + 2. * bulk(5)

!the partitioning between fast and slow SOM pools is a function of total soil clay mass
!after Jobbagy & Jackson (Ecol. App. 10, 2000)

slowfrac = t0 * claytot / 300. !300 kg m-2 is ca. the maximum value in the ISRIC WISE data, extrapolating layer 5 to 3m
fastfrac = 1. - slowfrac

!write(0,*)bulk
!write(0,*)clay
!write(0,*)claytot,soilmass
!write(0,*)fastfrac,slowfrac
!read(*,*)

mrh = 0.

do m = 1,12

  !litter decomposition: pft vectors

  do pft = 1,npft

    litterdag_fast(pft) = ek_lf(m) * litter_ag_fast(pft,1)
    litterdag_slow(pft) = ek_ls(m) * litter_ag_slow(pft,1)
    litterdag_bg(pft)   = ek_lf(m) * litter_bg(pft,1)        !belowground litter has the same turnover rate as fast litter

    !reduction of litter pool sizes

    litter_ag_fast(pft,1) = max(litter_ag_fast(pft,1) - litterdag_fast(pft),0.)
    litter_ag_slow(pft,1) = max(litter_ag_slow(pft,1) - litterdag_slow(pft),0.)
    litter_bg(pft,1)      = max(litter_bg(pft,1)      - litterdag_bg(pft),0.)

  end do

  !sum up decomposition over all pfts and all pools
  
  litter_decom(m) = sum(litterdag_fast + litterdag_slow + litterdag_bg)

  !write(0,*)'hetresp',mw1
  !write(0,*)tsoil
  !write(0,'(9f8.2)')litter_ag_fast(:,1) !,litterdag_fast
  !write(0,'(9f8.2)')litter_ag_slow(:,1) !,litterdag_slow
  !write(0,'(9f8.2)')litter_bg(:,1) !,litterdag_bg
  !write(0,'(i5,f8.2)')m,litter_decom(m)
  !write(0,*)'end hetresp'
  !read(*,*)

  !partition decomposing litter flux between atmosphere and SOM components

  cflux_litter_atmos = atmfrac  * litter_decom(m)
  cflux_litter_soil  = soilfrac * litter_decom(m)

  !further subdivide soil fraction between fast and slow soil pools

  cpool_fast(1) = cpool_fast(1) + cflux_litter_soil * fastfrac
  cpool_slow(1) = cpool_slow(1) + cflux_litter_soil * slowfrac
  
  !transfer remaining fast litter to surface SOM pool with turnover time of klit2som
  
  lit2som = dt * klit2som * litter_ag_fast(:,1)  !amount of litter to transfer to soil (vector over PFTs)

  litter_ag_fast(:,1) = max(litter_ag_fast(:,1) - lit2som,0.)

  cpool_surf(1) = cpool_surf(1) + sum(lit2som)  !sum across all pfts

  !calculate SOM decomposition
  
  cflux_surf_atmos = cpool_surf(1) * ek_lf(m) !eqn 4
  cflux_fast_atmos = cpool_fast(1) * ek_sf(m) !eqn 4
  cflux_slow_atmos = cpool_slow(1) * ek_ss(m) !eqn 4

  mrh(m,1) = cflux_litter_atmos + cflux_surf_atmos + cflux_fast_atmos + cflux_slow_atmos

  cpool_surf(1) = max(cpool_surf(1) - cflux_surf_atmos,0.)
  cpool_fast(1) = max(cpool_fast(1) - cflux_fast_atmos,0.)
  cpool_slow(1) = max(cpool_slow(1) - cflux_slow_atmos,0.)
  
end do

arh = sum(mrh,dim=1)

!write(0,'(4f10.4)')arh(1),sum(respflux(1,:)),sum(respflux(2,:)),sum(respflux(3,:))
!read(*,*)


!write(0,*)sum(litter_decom),arh(1)
!read(*,*)

!SOIL DECOMPOSITION EQUILIBRIUM CALCULATION
!Analytical solution of differential flux equations for fast and slow
!soil carbon pools.  Implemented after (soil_equil_year) simulation
!years, when annual litter inputs should be close to equilibrium.  Assumes
!average climate (temperature and soil moisture) from all years up to
!soil_equil_year.

if (spinup .and. year == soil_equil_year + 1 .and. k_fast_ave > 0. .and. k_slow_ave > 0.) then

  !Analytically calculate pool sizes this year only

  !Rate of change of soil pool size = litter input - decomposition
  !  (5) dc/dt = litter_decom - kc

  !At equilibrium,
  !  (6) dc/dt = 0

  !From (5) & (6),
  !  (7) c = litter_decom / k

  cpool_fast(1) = (soilfrac * fastfrac * litter_decom_ave(1)) / k_fast_ave   !eqn 7

  cpool_slow(1) = (soilfrac * slowfrac * litter_decom_ave(1)) / k_slow_ave   !eqn 7
     
else if (year <= soil_equil_year) then

  !Update running average respiration rates and litter input

  k_fast_ave = k_fast_ave + sum(k_soil_fast) / syr
  k_slow_ave = k_slow_ave + sum(k_soil_slow) / syr
  
  litter_decom_ave(1) = litter_decom_ave(1) + sum(litter_decom) / syr

end if
                  
end subroutine hetresp

!-------------------------------------------------------------------------------------------------------------

end module hetrespmod
