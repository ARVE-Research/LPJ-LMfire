module mpistatevarsmod

use mpi
use parametersmod, only : i2,i4,i8,sp,dp,npft,ncvar,climbuf,nspec,nistage,nisex,nirep
use weathergenmod, only : metvars_out
use orbitmod,      only : orbitpars

implicit none

public :: initstatevars

!-------------------------------------------
!input data

type climatedata
  real(sp), dimension(12) :: temp0   !mean monthly temperature of the previous year (degC)
  real(sp), dimension(12) :: temp    !mean monthly temperature (degC)
  real(sp), dimension(12) :: prec    !total monthly precipitation (mm)
  real(sp), dimension(12) :: cldp    !mean monthly bright cloudiness (percent)
  real(sp), dimension(12) :: wetd    !total monthly days with precipitation (days)
  real(sp), dimension(12) :: trng    !mean monthly diurnal temperature range (degC)
  real(sp), dimension(12) :: wind    !climatological mean monthly wind speed (m s-1)
  real(sp), dimension(12) :: lght    !climatological mean monthly lightning flashes (flashes ha-1)
end type climatedata  !84 elements

!---

type soildata
  real(sp), dimension(5)  :: zpos    !depth of layer midpoint from soil surface
  real(sp), dimension(5)  :: sand    !mass fraction
  real(sp), dimension(5)  :: clay    !mass fraction
  real(sp), dimension(5)  :: orgm    !organic matter mass (g m-2 ?) FLAG
  real(sp), dimension(5)  :: bulk    !bulk density (units?)
end type soildata  !26 elements

!---

type humandata
  real(sp)               :: foragerPD    !potential density of hunter-gatherers
  real(sp), dimension(3) :: popd         !human population density (persons km-2) (three different types)
  real(sp), dimension(8) :: landuse      !land use fractions (up to 8 different types of land use)
  real(sp)               :: lu_turnover  !fallow period before re-cultivation (yrs)
end type humandata  !13 elements

!---

type inputdata

  !values input to LPJ
  
  integer(i8)        :: idx      !gridcell index
  integer(i4)        :: xpos     !gridcell position in x
  integer(i4)        :: ypos     !gridcell position in y
  real(dp)           :: lon      !longitude of gridcell center
  real(dp)           :: lat      !latitude of gridcell center
  real(sp)           :: elev     !elevation of gridcell center (m.a.s.l.)
  real(sp)           :: slope    !mean gridcell slope (degrees)
  real(sp)           :: cellarea !area of the gridcell (m2)
  real(sp)           :: landf    !fraction of grid cell that is land
  logical            :: spinup   !are we in the model spinup
  real(sp)           :: co2      !co2 concentration
  integer            :: year     !simulation year (not calendar year), starts at 1
  type(orbitpars)    :: orbit
  type(climatedata)  :: climate
  type(soildata)     :: soil
  type(humandata)    :: human

end type inputdata   !150 elements

!integer, parameter :: nelem_input = 150

type(inputdata), allocatable, dimension(:) :: in_master

!---------------
!state variables

!---

type subgrid
  
  !scalars: 14 elements
  
  real(sp) :: afire_frac    !fraction of gridcell burnt this year
  real(sp) :: k_fast_ave    !running average k_fast for subroutine littersom
  real(sp) :: k_slow_ave    !running average k_slow for subroutine littersom
  real(sp) :: litterC_bg
  real(sp) :: litterC_fast
  real(sp) :: litterC_slow
  real(sp) :: livebiomass   !sum of living plant carbon across all PFTs
  real(sp) :: albiomass	    ! total carbon in aboveground living biomass (g m-2 yr-1) Chaste 03.03.2016
  real(sp) :: snow0         !Dec 31 snow pack
  real(sp) :: coverfrac     !land use cover fraction
  real(sp) :: tilecarbon    !total carbon by land use tile (living + littler + SOM)
  real(sp) :: soilerosion   !total annual soilerosion by tile (g m-2)
  real(sp) :: erflux        !soil erosion carbon flux (g m-2)
  real(sp) :: forager_pd    !forager potential population density (persons 100km-2)
  integer  :: cumfires      !cumulative number of burning fires
  real(sp) :: grasscover    !decadal running mean grass cover fraction
  real(sp) :: dgrassdt      !rate of change of decadal running mean grass cover fraction

  !soil layer: 2 elements

  !real(sp), dimension(2) :: zpos   !depth to soil layer midpoint
  !real(sp), dimension(2) :: dz     !soil layer thickness (m)
  !real(sp), dimension(2) :: sand   !mass %
  !real(sp), dimension(2) :: silt   !mass %
  !real(sp), dimension(2) :: clay   !mass %
  !real(sp), dimension(2) :: OM     !mass %
  !real(sp), dimension(2) :: OrgM   !organic matter (g m-2)
  !real(sp), dimension(2) :: bulk   !bulk density (g cm-3)
  !real(sp), dimension(2) :: Tsat   !saturation (volume fraction)
  !real(sp), dimension(2) :: T33    !33 kPa (volume fraction)
  !real(sp), dimension(2) :: T1500  !1500 kPa (volume fraction)
  !real(sp), dimension(2) :: whc    !water holding capacity T33 - T1500 (volume fraction)
  !real(sp), dimension(2) :: Ksat   !saturated conductivity (mm hr-1)
  real(sp), dimension(2) :: w      !instantaneous soil water content

  !npft 12*npft = 216 elements (ASSUMING 9 PFTS)

  logical,  dimension(npft) :: present       !whether PFT present in gridcell
  logical,  dimension(npft) :: leafon
  integer,  dimension(npft) :: leafondays
  integer,  dimension(npft) :: leafoffdays
  real(sp), dimension(npft) :: crownarea     !crown area (m2)
  real(sp), dimension(npft) :: dwscal365     !daily water scalar for day 365
  real(sp), dimension(npft) :: fpc_grid      !gridcell foliar projective cover (FPC)
  real(sp), dimension(npft) :: fpc_inc       !increment (if +ve) in FPC since last year
  real(sp), dimension(npft) :: height        !tree height (m)
  real(sp), dimension(npft) :: pftalbiomass  !Biomasse par PFT, E.C. 21.10.2016
  real(sp), dimension(npft) :: lai_ind       !individual leaf area index
  real(sp), dimension(npft) :: nind          !gridcell individual density (indiv/m2)
  real(sp), dimension(npft) :: plant_carbon
  real(sp), dimension(npft) :: above_carbon
  real(sp), dimension(npft,12) :: mlai	     !monthly mean LAI, per pft, per tile
  real(sp), dimension(npft,12) :: mBBpft     !biomass burned (kg dry matter, per pft, per month)

  !ncvar 9*ncvar = 27 elements (3 NCVARS)
  
  real(sp), dimension(ncvar) :: acflux_conv       !C flux to atmosphere due to anthropogenic deforestation (gC/m2)
  real(sp), dimension(ncvar) :: acflux_estab      !annual biomass increment due to establishment (gC/m2)
  real(sp), dimension(ncvar) :: acflux_fire       !C flux to atmosphere due to fire (gC/m2)
  real(sp), dimension(ncvar) :: arh               !annual heterotrophic respiration (gC/m2)
  real(sp), dimension(ncvar) :: cpool_surf        !surface SOM pool (gC/m2)
  real(sp), dimension(ncvar) :: cpool_fast        !fast-decomposing soil C pool (gC/m2)
  real(sp), dimension(ncvar) :: cpool_slow        !slow-decomposing soil C pool (gC/m2)
  real(sp), dimension(ncvar) :: grid_npp          !gridcell total npp (gC/m2)
  real(sp), dimension(ncvar) :: grid_gpp          !gridcell total gpp (gC/m2) ! E.C. avril 2016
  real(sp), dimension(ncvar) :: litter_decom_ave  !running average litter_decom for subroutine littersom
  
  !npft,ncvar 9*npft*ncvar = 243 (9 PFTS * 3 NCVARS)
  
  real(sp), dimension(npft,ncvar) :: anpp           !annual gridcell NPP (gC/m2)
  real(sp), dimension(npft,ncvar) :: agpp           !annual gridcell GPP (gC/m2) ! E.C. avril 2016
  real(sp), dimension(npft,ncvar) :: lm_ind         !individual leaf mass (gC)
  real(sp), dimension(npft,ncvar) :: sm_ind         !individual sapwood mass (gC)
  real(sp), dimension(npft,ncvar) :: hm_ind         !individual heartwood mass (gC)
  real(sp), dimension(npft,ncvar) :: rm_ind         !individual fine root mass (gC)
  real(sp), dimension(npft,ncvar) :: cstore         !year-to-year labile carbon store (gC)
  real(sp), dimension(npft,ncvar) :: litter_ag_fast !gridcell above-ground litter (gC/m2)
  real(sp), dimension(npft,ncvar) :: litter_ag_slow !gridcell above-ground litter (gC/m2)
  real(sp), dimension(npft,ncvar) :: litter_bg      !gridcell below-ground litter (gC/m2)
  
  !annual trace gas emissions nspec (6 SPECIES)
  real(sp), dimension(nspec) :: aMx
  
  !monthly burned fraction of gridcell area
  real(sp), dimension(12) :: mburnedf
  
  !state of the soil, 5 layers * 5 elements = 25
  type(soildata) :: soil
  
  !historical burned fraction (20 years) = 20 elements
  real(sp), dimension(climbuf) :: burnedf_buf
  real(sp), dimension(climbuf) :: forager_pd_buf
  
  !insect mass
  real(sp), dimension(nistage,nisex,nirep) :: insectstate   !fraction of life stage completed (0-1)
  real(sp), dimension(nistage,nisex,nirep) :: insectmass    !biomass of insects in each life stage (mg/m2)
  real(sp), dimension(nistage,nisex,nirep) :: insectenergy  !fractional energetic state of the insect (0-1)

end type subgrid !14 + 2 + 216 + 27 + 243 + 6 + 25 + 20 + 20 = 573

!---

type carbonstock
  real(sp), dimension(npft) :: crop_harvest
  real(sp) :: wood_fast
  real(sp) :: wood_slow
  real(sp), dimension(2) :: prod_flux
end type carbonstock  !5 elements

!---

type statevars
  
  type(metvars_out) :: met     !30 
  type(carbonstock) :: carbon  !5
  real(sp) :: mtemp_min20 !20-year average minimum monthly temperature (deg C)
  real(sp) :: gdd20       !20-year average growing degree days
  real(sp) :: mat20       !20-year average mean annual temperature
  real(sp), dimension(climbuf) :: mat_buf        !buffer to store 'climbuf' years of mean annual temperature (20 elements each)
  real(sp), dimension(climbuf) :: gdd_buf        !buffer to store 'climbuf' years of GDD totals (20 elements each)
  real(sp), dimension(climbuf) :: mtemp_min_buf  !buffer to store 'climbuf' years of coldest month temperatures
  type(subgrid), dimension(3) :: tile

  real(sp), dimension(3) :: annburntarget  !desired fraction of the gridcell to burn: foragers, farmers, pastoralists
  !type(subgrid), allocatable, dimension(:) :: tile

end type statevars

!integer, parameter :: nelem_sv_scalar = 30+5+1+1+1+20+20+20   != 98
!integer, parameter :: nelem_sv_tile   = 14+2+216+27+243+6+25+20+20  != 573

type(statevars), allocatable, dimension(:) :: sv_master

!example of how this is indexed

!sv(1:ncells)%tile(1:ntiles)%fpc_grid(1:npft)

contains

!-------------------------------------------------------------------------------------------------------------------------------------

subroutine initstatevars(in,sv,ismaster)

use geohashmod,    only : geohash
use randomdistmod, only : ran_seed

implicit none

type(inputdata), intent(inout)  :: in
type(statevars), intent(out) :: sv

logical, intent(in) :: ismaster

integer :: i
integer :: j
integer :: k

!---------
!initialize the random number state based on random coordinates

if (ismaster) then

  !write(*,*)'seeding with',in%lon,in%lat,geohash(in%lon,in%lat)

  call ran_seed(geohash(in%lon,in%lat),sv%met%rndst)
  
end if
  
!write(*,*)'random seed',sv%met%rndst
!write(*,*)


in%soil%bulk = 0.

!initialize arbitrary for first day of first spinup year

sv%met%NI    = 0.
sv%met%pday  = .false.
sv%met%resid = 0.1

sv%met%tmin  = 15.
sv%met%tmax  = 15.
sv%met%cldf  = 0.5
sv%met%wind  = 0.
sv%met%lght  = 0.
sv%met%resid = 0.

!---

sv%annburntarget = [ 0.5, 0.05, 0.2 ]  !start out letting foragers try to burn everything

!---

sv%mat20       = 0.
sv%gdd20       = 0.
sv%mtemp_min20 = 0.

sv%carbon%wood_fast    = 0.
sv%carbon%wood_slow    = 0.
sv%carbon%prod_flux(1) = 0.
sv%carbon%prod_flux(2) = 0.

forall (i=1:climbuf)
  sv%mat_buf(i) = -9999.
  sv%gdd_buf(i) = -9999.
  sv%mtemp_min_buf(i) = -9999.
end forall

sv%tile%w(1) = 1.
sv%tile%w(2) = 1.
sv%tile%snow0 = 0.
sv%tile%cumfires = 0
sv%tile%grasscover = 0.
sv%tile%dgrassdt = 0.

sv%tile%afire_frac  = 0.
sv%tile%k_fast_ave  = 0.
sv%tile%k_slow_ave  = 0.

sv%tile%forager_pd = 7.e-3    !Initialize every grid cell with roughly two people per 3000 km-2

forall (i=1:climbuf)
  sv%tile%burnedf_buf(i) = 0.
  sv%tile%forager_pd_buf(i) = 0. 
end forall

forall (i=1:npft)
  sv%tile%present(i)      = .false.
  sv%tile%dwscal365(i)    = 1.
  sv%tile%crownarea(i)    = 0.
  sv%tile%fpc_grid(i)     = 0.
  sv%tile%fpc_inc(i)      = 0.
  sv%tile%height(i)       = 2.
  sv%tile%lai_ind(i)      = 0.
  sv%tile%nind(i)         = 0.
  sv%tile%plant_carbon(i) = 0.
  sv%tile%above_carbon(i) = 0.
  
  forall (j = 1:12)
    sv%tile%mBBpft(i,j) = 0.
  end forall

  forall (j=1:ncvar)
    sv%tile%hm_ind(i,j) = 0.
    sv%tile%lm_ind(i,j) = 0.
    sv%tile%sm_ind(i,j) = 0.
    sv%tile%rm_ind(i,j) = 0.
    sv%tile%cstore(i,j) = 0.

    sv%tile%litter_ag_fast(i,j) = 0.
    sv%tile%litter_ag_slow(i,j) = 0.
    sv%tile%litter_bg(i,j) = 0.

    sv%tile%anpp(i,j) = 0.
  end forall
end forall

forall (j=1:ncvar)
  sv%tile%grid_npp(j)     = 0.
  sv%tile%grid_gpp(j)     = 0. ! E.C. avril 2016
  sv%tile%acflux_estab(j) = 0.
  sv%tile%acflux_fire(j)  = 0.
  sv%tile%acflux_conv(j)  = 0.
  sv%tile%arh(j)          = 0.

  sv%tile%cpool_surf(j)       = 0.   
  sv%tile%cpool_fast(j)       = 0.   
  sv%tile%cpool_slow(j)       = 0.
  sv%tile%litter_decom_ave(j) = 0.
end forall

forall (k = 1:nirep)
  forall (j = 1:nisex)
    forall (i = 1:nistage)
      sv%tile%insectstate(i,j,k) = 0.
      sv%tile%insectmass(i,j,k)  = 0.
      sv%tile%insectenergy(i,j,k)  = 0.
    end forall

    !initialization mass of 1 l2o spruce budworm per m2 of leaf area
    !assume: initial LAI of 3 and L2o budworm mass (0.042 mg per worm)
    !0.042 * 3 = 0.126 / 2 for sex = 0.063 (starting mass)
    
    !now not done here because we will initialize in LPJmod under appropriate conditions for colonization

!     sv%tile%insectmass(1,j,k)   = 0.063
!     sv%tile%insectenergy(1,j,k) = 2.1571  !starting calories for L2o

  end forall
end forall

sv%tile%livebiomass  = 0.
sv%tile%albiomass  = 0.
sv%tile%litterC_fast = 0.
sv%tile%litterC_slow = 0.
sv%tile%litterC_bg   = 0.

sv%tile%coverfrac    = 0.   !initialize all categories to zero
sv%tile(1)%coverfrac = 1.   !set natural vegetation to 100% cover

end subroutine initstatevars

!-------------------------------------------------------------------------------------------------------------------------------------

end module mpistatevarsmod
