module spitfiremod

use parametersmod,    only : sp,dp,npft

implicit none

public  :: spitfire
public  :: burnedbiomass
public  :: managedburn 
private :: calcROS
private :: firemortality

real(sp), parameter :: me_lf  =  0.2       !moisture of extinction for live grass fuels: fractional moisture above which fuel does not burn (fraction) (20%)
real(sp), parameter :: h      = 18.        !heat content of fuel (kJ g-1)  5kwh/kgs
real(sp), parameter :: c2om   = 1. / 0.45  !conversion factor between Carbon and total mass of vegetation
real(sp), parameter :: ST     = 0.055      !fraction of total vegetation mass that is mineral (non-flammable)
real(sp), parameter :: ot     = 1. / 3.    !one third
real(sp), parameter :: orgf   = 1. - ST    !organic fraction of biomass

real(sp), parameter :: walkdist = 10000.  !distance one arsonist can cover in one day

real(sp), dimension(4) :: me_fc = [ 0.404, 0.487, 0.525, 0.544 ]

!PFT-specific parameters for fire damage (only for trees)
!NB parameters that are commented out are those given in Thonicke et al. (2010)
!alternate parameter sets come from the SPITFIRE f77 code - lpjmain_spitfire270106.f
  

real(sp), dimension(npft) :: ieffpft !ignition efficiency parameter

real(sp), dimension(npft) :: m_ex    !moisture of extinction

real(sp), dimension(npft) :: emCO2   !emission factor for CO2
real(sp), dimension(npft) :: emCO    !emission factor for CO
real(sp), dimension(npft) :: emCH4   !emission factor for CH4
real(sp), dimension(npft) :: emVOC   !emission factor for VOC
real(sp), dimension(npft) :: emTPM   !emission factor for TPM
real(sp), dimension(npft) :: emNOx   !emission factor for NOx from F77 code

!-----------------
!parameters for characterization of fuel (bulk density and surface area to volume ratio)

!fuel bulk density (kg m-3)
real(sp), dimension(npft) :: rhobPFT !these values averages from different USFS fuel models

real(sp), parameter :: rho_soilC = 75.

real(sp), parameter :: cm2ft  =  30.48     !cm in one foot

!Surface-area-to-volume ratios of fuel classes
!real(sp), parameter, dimension(3) :: surf2vol = [ 2750., 109., 30. ]  !average values from Andrews, 1986 (ft2 ft-3)
real(sp), parameter, dimension(3) :: surf2vol = [ 2012., 109., 30. ]  !average values from Andrews, 1986 (ft2 ft-3)
real(sp), parameter, dimension(3) :: sigma_i  = surf2vol / cm2ft       !here converted to (cm2 cm-3)

!-----------------

integer, parameter :: nspec = 6

!emission factors in units g kg-1 (dry matter)

!real(sp), parameter, dimension(npft,nspec) :: emfact = [emCO2, emCO, emCH4, emVOC, emTPM, emNOx ]  !FLAG - this assignment not compatible with gfortran

real(sp) :: precsum = 0.

real(sp) :: sumfdi    = 0.
real(sp) :: totburn   = 0.
integer  :: anumfires = 0
integer  :: lastyear  = 0
integer  :: burndays  = 0
integer  :: lasttile  = 1
real(sp) :: unburneda  !unburned area (ha)

real(sp), dimension(npft,5) :: annBBdead    !annual total biomass burned from dead fuel by fuel type and PFT (g m-2)
real(sp), dimension(npft,3) :: annBBlive    !annual biomass burned from live fuel by fuel type and PFT (g m-2)
real(sp), dimension(npft)   :: ann_kill     !annual total probability of mortality

contains

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine managedburn(i,j,acflux_fire,afire_frac,litter_ag_fast,pftCflux)

use parametersmod,   only : npft

implicit none

integer, intent(in) :: i  !tile index
integer, intent(in) :: j  !gridcell index

real(sp),               intent(out)   :: acflux_fire     !carbon flux from biomass burning (gC m-2 d-1)
real(sp),               intent(out)   :: afire_frac      !fraction of the gridcell burned
real(sp), dimension(:), intent(inout) :: litter_ag_fast  !aboveground litter, fast turnover (g C m-2) per PFT

real(sp), dimension(npft), intent(inout) :: pftCflux

!-------

!implementation of managed burning, 20% of used land area, with C emissions based on 20% of aboveground biomass

acflux_fire = 0.

afire_frac = 0.2   !20% of agricultural tile

pftCflux   = litter_ag_fast * afire_frac    !100% of the crop remains are burned on the burned fraction of managed land

litter_ag_fast = litter_ag_fast - pftCflux   !reduction in aboveground litter

acflux_fire = sum(pftCflux)

end subroutine managedburn

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine spitfire(year,i,j,d,input,met,soilwater,snowpack,dphen,wscal,osv,spinup,avg_cont_area,burnedf20,forager_pd20,FDI,omega_o0,omega0,BBpft,Ab,ind)

use parametersmod,   only : pir
use parametersmod,   only : npft,pi,pft,pftpar
use weathergenmod,   only : metvars_out
use mpistatevarsmod, only : inputdata,statevars
use randomdistmod,   only : randomstate,ranu,rng1,half
use individualmod,   only : sizeind

implicit none

!arguments

integer, intent(in) :: year  !year number, for diagnostic output, to be removed later
integer, intent(in) :: i  !tile index
integer, intent(in) :: j  !gridcell index
integer, intent(in) :: d  !day of the year
type(inputdata),   intent(in) :: input
type(metvars_out), intent(inout) :: met
real(sp),          intent(in) :: soilwater
real(sp),          intent(in) :: snowpack
real(sp),          intent(in) :: burnedf20           !20-year running mean of burned fraction
real(sp), 	   intent(in) :: forager_pd20		!20-year running mean forager population density km-1 as calculated by foragersmod
real(sp),          intent(inout) :: avg_cont_area        !average contiguous area size of non-used part of the gridcell (ha)
real(sp),	   intent(inout) :: FDI
real(sp), 	   intent(inout) :: omega_o0
real(sp), dimension(4), intent(inout) :: omega0     !moisture content of each fuel class on the previous day
real(sp), dimension(:), intent(in) :: dphen                !phenological state (1=full leaf; 0=no leaves) per pft, on this day
real(sp), dimension(:), intent(in) :: wscal          !water scalar (supply/demand <= 1) per pft, on this day
logical, intent(in) :: spinup

type(sizeind), dimension(:,:) :: ind   !size of the individuals by PFT and height class

real(sp), dimension(:), intent(out) :: BBpft  !daily total biomass burned from each PFT, sum of live and dead consumed (kg dry matter m-2)

type(statevars), target, intent(inout) :: osv

!parameters

real(sp), parameter :: b  = 10. / 9.
real(sp), parameter :: c  =  1. / 9.

real(sp), parameter, dimension(4) :: alpha = 1.5 * [1.e-3, 5.424e-5, 1.485e-5, 1.e-6]  !drying parameter for 1- 10- 100-, and 1000-h fuel classes (degC-2)

real(sp), parameter :: min2sec = 1. / 60.

!real(sp), parameter, dimension(3) :: annburntarget = [ 0.25, 0.05, 0.2 ]  !desired fraction of the gridcell to burn: foragers, farmers, pastoralists

!state variables (scalars)

real(sp), dimension(3) :: PD           !population density (individuals km-2)  (1=hunter gatherers, 2=farmers, 3=pastoralists)
real(sp) :: NI           !Nesterov fuel dryness Index (degC^2)
real(sp) :: Ustar        !day mean wind speed (m s-1)
real(sp) :: area         !gridcell area (m2)
real(sp) :: light        !frequency of total lighting flashes (ha-1 d-1)
real(sp) :: dwind        !wind speed

real(sp), pointer :: acflux_fire  !carbon flux from biomass burning (gC m-2 d-1)
real(sp), pointer :: afire_frac   !fraction of the gridcell burned

!pointers to PFT state variables (vectors)

logical,  pointer, dimension(:) :: present         !PFT is present
real(sp), pointer, dimension(:) :: fpc_grid        !PFT cover fraction
real(sp), pointer, dimension(:) :: height          !tree height (m)
real(sp), pointer, dimension(:) :: litter_ag_fast  !aboveground litter, fast turnover (g C m-2)
real(sp), pointer, dimension(:) :: litter_ag_slow  !aboveground litter, slow turnover (g C m-2)
real(sp), pointer, dimension(:) :: litter_bg       !belowground litter (g C m-2)
real(sp), pointer, dimension(:) :: cpool_surf      !surface SOM pool (2yr turnover time)
real(sp), pointer, dimension(:) :: nind            !gridcell individual density (indiv/m2)
real(sp), pointer, dimension(:) :: lm_ind          !leaf mass, average individual (g C)
real(sp), pointer, dimension(:) :: sm_ind          !sapwood mass, average individual (g C)
real(sp), pointer, dimension(:) :: hm_ind          !heartwood mass, average individual (g C)
real(sp), pointer, dimension(:) :: rm_ind          !root mass, average individual (g C)
real(sp), pointer, dimension(:) :: aMx             !annual trace gas emissions (g m-2; per species)

integer, pointer :: cumfires  !number of fires currently burning


!local variables

real(sp), dimension(3) :: annburntarget   !desired fraction of the gridcell to burn: foragers, farmers, pastoralists

real(sp), dimension(4) :: dry
real(sp) :: wet
real(sp) :: dry_o
real(sp) :: prob

integer  :: l

real(sp), intent(inout) :: Ab             !area burned (ha d-1)
real(sp) :: Abfrac         !fractional area burned on the gridcell (fraction)
real(sp) :: DT             !length of the major axis of the fire (total distance traveled) (m)


real(sp) :: Isurface       !surface fire intensity (kW m-1)
real(sp) :: LB
real(sp) :: LBgrass
real(sp) :: LBtree
real(sp) :: ROSfsurface    !overall forward rate of spread of fire (m min-1)
real(sp) :: ROSfsurface_g  !forward rate of spread of fire in herbaceous fuel (m min-1)
real(sp) :: ROSfsurface_w  !forward rate of spread of fire in woody fuels (m min-1)
real(sp) :: ROSfcrown      !forward rate of spread of fire in crown (m min-1)
real(sp) :: ROSbsurface    !backward rate of spread of fire (m min-1)
real(sp) :: Uforward       !mean wind speed (m min-1)
real(sp) :: abarf
real(sp) :: area_ha        !gridcell area (ha)
real(sp) :: kPD

real(sp) :: omega_lg       !relative moisture content of live grass
real(sp) :: omega_n        !relative moisture content of the 1-h fuel class
real(sp) :: omega_nl       !mean relative moisture contents of 1-h fuel class and live gras
real(sp) :: omega_o        !relative daily litter moisture
real(sp) :: omega_s1       !relative moisture content of top soil layer
real(sp) :: rm             !moisture relative to moisture of extinction (omega_o / m_e)
real(sp) :: tau_l          !fire residence time (min)
real(sp) :: tfire          !fire duration (min)
real(sp) :: wlivegrass     !mass of live grass (g C m-2)
real(sp) :: wo             !total mass of dead fuel summed across all fuel classes (kg m-2)
real(sp) :: BBtot          !total biomass burned (g C m-2)

integer, save  :: nhig     !number of human-caused ignition events (per gridcell d-1)
real(sp) :: nlig           !frequency of lightning-caused ignition events (ha-1 d-1)

real(sp) :: rdf
real(sp) :: rlf
real(sp) :: caf
real(sp) :: wtot
real(sp) :: wfinefuel

real(sp) :: canopywater  !integrated mean wscal

real(sp) :: m_e     !mass weighted average moisture of extinction for dead fuel
real(sp) :: me_avg  !mass weighted average moisture of extinction for all fuels (dead + live grass)
real(sp) :: me_nl   !mass weighted average moisture of extinction for live grass and 1hr fuel

real(sp) :: alpha_df
real(sp) :: alpha_lg
real(sp) :: CFlg      !consumed fraction of live grass

real(sp) :: totfuel   !total carbon in all live biomass, litter, and surface soil C
real(sp) :: netfuel   !total organic part of the 1- 10- and 100-h fuel
real(sp) :: Ukmh      !wind speed in km h-1

real(sp) :: treecover
real(sp) :: grascover
real(sp) :: totvcover

real(sp) :: Z_ratio                         !ratio between cloud-cloud and cloud-ground lightning, depending on latitude (Prentice and Mackerras, 1977)
real(sp) :: cg_light                         !number of cloud-ground lightning strikes
real(sp) :: lat                                 !latitude in radians

real(sp), dimension(4)      :: woi       !1-, 10-, 100- and 1000-h dead fuel mass summed across all PFTs (g m-2)
real(sp), dimension(3)      :: woi_g     !dead fuel mass of grass PFTs (g m-2), 1-, 10- and 100-h class
real(sp), dimension(3)      :: woi_w     !dead fuel mass of woody PFTs (g m-2), 1-, 10- and 100-h class
real(sp), dimension(4)      :: omega     !moisture content of each fuel class
real(sp), dimension(4)      :: CF        !fractional consumption of dead fuel per fuel class
real(sp), dimension(4)      :: FC        !amount of dead fuel consumed, per fuel class (g m-2)
real(sp), dimension(npft)   :: P_m       !total probability of mortality
real(sp), dimension(npft)   :: CK        !fraction of crown scorch (fraction)
real(sp), dimension(npft)   :: nind_kill !fraction of PFT killed by fire (per pft)
real(sp), dimension(npft,4) :: livefuel  !live fuel load in 1- 10- 100- and 1000- hour fuel classes (g DM m-2)
real(sp), dimension(npft,4) :: deadfuel  !dead fuel load in 1- 10- 100- and 1000- hour fuel classes (g DM m-2)
real(sp), dimension(npft,5) :: BBdead    !biomass burned from dead fuel by fuel type and PFT (g m-2)
real(sp), dimension(npft,3) :: BBlive    !biomass burned from live fuel by fuel type and PFT (g m-2)
  
real(sp), dimension(npft,4) :: lfuelg  !live fuel load in 1- 10- 100- and 1000- hour fuel classes that is grass (g DM m-2)
real(sp), dimension(npft,4) :: dfuelg  !dead fuel load in 1- 10- 100- and 1000- hour fuel classes that is grass (g DM m-2)

real(sp), dimension(npft,4) :: lfuelw  !live fuel load in 1- 10- 100- and 1000- hour fuel classes that is wood (g DM m-2)
real(sp), dimension(npft,4) :: dfuelw  !dead fuel load in 1- 10- 100- and 1000- hour fuel classes that is wood (g DM m-2)

real(sp), dimension(nspec) :: Mx     !trace gas emissions (g x m-2)

integer        :: numfires  !number of fires started on the grid cell today
integer               :: numfires_hum
integer               :: numfires_nat
integer, save  :: peopfire_account
real(sp) :: riskfact         !FDI dependent factor influencing people's fire behavior; people become more careful when fire danger is high

!real(sp), parameter, dimension(3) :: max_ig = [15., 2. , 8.]                 !maximum number of human possible fire ignitions (Mha-1 day-1) 

real(sp) :: burnedf
real(sp) :: ieff
real(sp) :: rho_livegrass
real(sp) :: wind_multiplier
real(sp) :: ieff_avg
real(sp) :: wf
real(sp) :: latscale
real(sp) :: slopefact
real(sp) :: dayburntarget
real(sp) :: targetfact

real(sp) :: SOM_surf   !mass of the total organic matter in the surface soil C pool

real(sp) :: gscale

!---
!variables below copied from original ROS subroutine

real(sp) :: wn
real(sp) :: livemass
real(sp) :: deadmass
real(sp) :: totlfuelw
real(sp) :: totdfuelw
real(sp) :: totlfuelg
real(sp) :: totdfuelg
real(sp) :: rho_b          !weighted average of fuel bulk density (kg m-3)
real(sp) :: sigma          !surface area to volume ratio of the fuel (cm2 cm-3)
real(sp) :: relmoist       !relative moisture content of the fuel relative to its moisture of extinction (unitless; omega_o / m_e)

logical, dimension(npft) :: fuelpft      !true if there is nonzero fuel for this PFT
real(sp),dimension(npft) :: pftlivefuel  !mass of live herbaceous fuel per pft summed over 1- 10- and 100-h fuel classes (g m-2)
real(sp),dimension(npft) :: pftdeadfuel  !mass of dead fuel per pft summed over 1- 10- and 100-h fuel classes (g m-2)
real(sp),dimension(npft) :: rho_pft      !bulk density of dead fuel per pft mass weighted average over 1- 10- and 100-h fuel classes (kg m-3)

logical :: crownfire
logical :: bavard = .false.

integer :: people   !total number of people on the gridcell
logical :: calchumanfire

integer  :: uniquefires  !number of unique fires one person can start in one day (without burning into one another)
real(sp) :: firewidth    !width of the burning ellipse (m)
real(sp) :: humfire1     !area one person can burn in one day (ha)

integer, save :: arsonists
integer :: group

!------------------
!assignment

! Enregistrement des parametres 
ieffpft  = pftpar(:,44)
emCO2    = pftpar(:,45)
emCO     = pftpar(:,46)
emCH4    = pftpar(:,47)
emVOC    = pftpar(:,48)
emTPM    = pftpar(:,49)
emNOx    = pftpar(:,50) 
rhobPFT  = pftpar(:,51) 


!write(*,'(2i5,2f7.2,f7.1,f7.2)')year,d,met%tmin,met%tmax,met%prec,met%dayl

Ab = 0.
abarf = 0.
BBpft = 0.

calchumanfire = .false.

!PD    = input%human%popd 

PD = 0.

!PD(1) = forager_pd20

!write(0,*) 'Forager PD spitfire: ', year, PD 

!PD(2) = 0.

area  = input%cellarea   !m2

area_ha = 1.e-4 * area    !convert m2 to ha

if(input%slope >= 1.72) then   !0.03 for radians
!   slopefact = 1. / (100. * input%slope - 2.)                !slope > 0.03, for slope coming in as radians
   slopefact = 1. / (5. / 9. * pi * input%slope - 2)	     ! this one for slope coming in in degrees
else 
   slopefact = 1.
end if 

light = met%lght * 0.01 !convert from km-2 to ha-1
Ustar = met%wind
NI    = met%NI

!if( .not. spinup .and. year>=114) write(0,'(2i4,2f14.4)') year, d, light, light*area_ha

omega_s1 = soilwater   !top layer soil water content

!take only the C-12 content for the vegetation state variables below

present => osv%tile(i)%present

litter_ag_fast => osv%tile(i)%litter_ag_fast(:,1)  !dead biomass
litter_ag_slow => osv%tile(i)%litter_ag_slow(:,1)
litter_bg      => osv%tile(i)%litter_bg(:,1)

cpool_surf     => osv%tile(i)%cpool_surf
SOM_surf = cpool_surf(1) * c2om

nind           => osv%tile(i)%nind                 !individual density
lm_ind         => osv%tile(i)%lm_ind(:,1)          !live biomass state variables
sm_ind         => osv%tile(i)%sm_ind(:,1)
hm_ind         => osv%tile(i)%hm_ind(:,1)
rm_ind         => osv%tile(i)%rm_ind(:,1)

fpc_grid       => osv%tile(i)%fpc_grid             !PFT cover fraction
acflux_fire    => osv%tile(i)%acflux_fire(1)       !fire flux (only total C not isotopes for now)
afire_frac     => osv%tile(i)%afire_frac           !annual burned fraction - summed from each daily timestep

height         => osv%tile(i)%height               !vectors across all PFTs

aMx            => osv%tile(i)%aMx
cumfires       => osv%tile(i)%cumfires

!------------------
!annual stats

!write(0,'(a,9f12.2)') 'height:', height

crownfire = .false.

if (d == 1) then
  !write(0,*)'writing annual stats'
  sumfdi = aMx(5)  
    
  anumfires = 0                        !FLAG
  sumfdi   = 0.
  totburn  = 0.
  burndays = 0 
  lastyear = year
  lasttile = i
  acflux_fire = 0.
  aMx = 0.
  peopfire_account = 365
  
  annBBdead = 0.
  annBBlive = 0.
  ann_kill  = 0.
  
  precsum = 0.
  
  omega0 = 1.  !set fuel moisture to fully wet on first day (should handle this differently and carry through as state variable)

  !write(0,*)'done writing annual stats'

  !unburneda = area_ha

  !if (cumfires > 0) then
  !  write(0,*) 'BEGIN SPITFIRE: ', year, PD, cumfires !light * area_ha
  !end if

end if

if (snowpack > 0.) then  !no fire on days with snow on the ground
  if(bavard)  write(*,'(a16,4i6,6f14.7)') 'snowpack ', year,i,d,cumfires, met%prec, NI, afire_frac, PD
  cumfires = 0 		! extinguish all smouldering or burning fires  
  return
end if  

if (met%prec == 0.) then
  precsum = 0.
else
  precsum = precsum + met%prec
end if

if ((precsum >= 10. .and. sum(fpc_grid(8:9)) <= 0.6) .or. (precsum >= 3. .and. sum(fpc_grid(8:9)) > 0.6)) then  !extinguish all currently burning or smoldering fires
  cumfires = 0 
end if

  

!---
!bulk density of standing grass biomass (function of average GDD)
!ranges from about 12 at tundra GDDs to 1 for tropical GDDs
  
rho_livegrass = 2.e4 / (osv%gdd20 + 1000.) - 1.                

!---

unburneda = area_ha - totburn

!write(0,'(a,18f10.2)')'BEGIN spitfire',litter_ag_fast,litter_ag_slow

if (any(litter_ag_fast < 0.) .or. any(litter_ag_slow < 0.)) then
  write(0,*)'negative litter spitfire 1!',i,j
  write(0,'(18f12.5)')litter_ag_fast,litter_ag_slow
  stop
end if

!calculations start here
!write(0,) 'entering SPITFIRE ',year,i,d,PD,NI

BBlive = 0.
BBdead = 0.

!------------------
!calculate fuel by PFT and by fuel type
!litter_ag_slow includes all woody litter, whereas fast is leaves only

deadfuel(:,1) = c2om * (0.045 * litter_ag_slow + litter_ag_fast)  !1-h fuel class 
deadfuel(:,2) = c2om * 0.075 * litter_ag_slow                   !10-h fuel class 
deadfuel(:,3) = c2om * 0.21  * litter_ag_slow                   !100-h fuel class 
deadfuel(:,4) = c2om * 0.67  * litter_ag_slow                   !1000-h fuel class 

livefuel(:,1) = c2om * nind * (0.045 * (hm_ind + sm_ind) + lm_ind)! * dphen)  !1-h fuel class
livefuel(:,2) = c2om * nind * (0.075 * (hm_ind + sm_ind))           !10-h fuel class
livefuel(:,3) = c2om * nind * (0.21  * (hm_ind + sm_ind))           !100-h fuel class
livefuel(:,4) = c2om * nind * (0.67  * (hm_ind + sm_ind))           !1000-h fuel class

!===============================    This distinction only for calculation of ROSfsurface, separately calculate it in grassland and woodland, then weight it a posteriori using fpc_grid

lfuelg = 0.
lfuelw = 0.
dfuelg = 0.
dfuelw = 0.

do l = 1,4
  where (pft%tree)
    lfuelw(:,l)= livefuel(:,l)		! FLAG: changed the loop variable from "i" to "l", have no idea how or if that worked at all before (MP, 28.02.12)
    dfuelw(:,l)= deadfuel(:,l)
  elsewhere
    lfuelg(:,l)= livefuel(:,l)
    dfuelg(:,l)= deadfuel(:,l)
  end where
end do

!===============================

if (i == 0 .and. d == 1) then

  write(0,*)'spitfire initial fuel loads',year
  
  write(0,*)'dead'
  write(0,'(a,9f10.1)')'1',deadfuel(:,1)
  write(0,'(a,9f10.1)')'2',deadfuel(:,2)
  write(0,'(a,9f10.1)')'3',deadfuel(:,3)
  write(0,'(a,9f10.1)')'4',deadfuel(:,4)

  write(0,*)'live'
  write(0,'(a,9f10.1)')'1',livefuel(:,1)
  write(0,'(a,9f10.1)')'2',livefuel(:,2)
  write(0,'(a,9f10.1)')'3',livefuel(:,3)
  write(0,'(a,9f10.1)')'4',livefuel(:,4)

end if

!dead fuel load

woi(1) = sum(deadfuel(:,1))  !1-h fuel class summed across all PFTs (g dry biomass m-2)
woi(2) = sum(deadfuel(:,2))  !10-h fuel class summed across all PFTs (g dry biomass m-2)
woi(3) = sum(deadfuel(:,3))  !100-h fuel class summed across all PFTs (g dry biomass m-2)
woi(4) = sum(deadfuel(:,4))  !1000-h fuel class summed across all PFTs (g dry biomass m-2)

!===============================    This distinction only for calculation of ROSfsurface, separately calculate it in grassland and woodland, then weight it a posteriori using fpc_grid

woi_g = sum(deadfuel(8:9,1:3),dim=1)
woi_w = sum(deadfuel(1:7,1:3),dim=1)

!===============================

!weighted average reduction factor to account for woody and non woody parts of the gridcell, also convert m s-1 to m min-1

treecover = sum(fpc_grid,mask=pft%tree)
grascover = sum(fpc_grid,mask=.not.pft%tree)
totvcover = sum(fpc_grid)

!write(0,'(a,4i5,9f10.2)')'enter spf',year,d,i,j,fpc_grid  !lm_ind(9),dphen(9)

netfuel = (1. - ST) * sum(woi(1:3))  !net organic part of the fuel (g m-2)

totfuel = sum(deadfuel) + sum(livefuel) + cpool_surf(1)

if (totfuel < 1000. .or. totvcover < 0.5) then
  if(bavard)  write(*,'(a16,4i6,13f14.4)')'no_fuel',year,i,d, cumfires, Ab, abarf, afire_frac, light*area_ha, 0., PD, met%prec, NI, PD
  cumfires = 0  
  return  !no fuel
end if

!---------------------------

!Uforward = 60. * (1. - totvcover) * Ustar + totvcover * (0.4 * Ustar * treecover + 0.6 * Ustar * grascover)

!from f77 spitfire version - only on the vegetated part of the gridcell

Ustar = max(Ustar, 0.)

Uforward = 60. * Ustar !(0.4 * Ustar * treecover + 0.6 * Ustar * grascover)

!Uforward = 3. * Uforward !FLAG remove!!!!

!if (totvcover > 0.) then
!  Uforward = 60. * (0.4 * Ustar * treecover + 0.6 * Ustar * grascover) / totvcover
!else
  !write(0,*)'no live vegcover'
!  Uforward = 60. * Ustar
!end if

!---------------------------
!part 2.2.1, ignition events

!lightning

!efficiency of lightning ignition decreases with increasing burned area

!calculate weighted average ieff

ieff_avg = sum(fpc_grid * ieffpft) / sum(fpc_grid)
!ieff_avg = 1.

burnedf = totburn / area_ha
!ieff = 0.2 * (1. - burnedf) / (1. + 25. * burnedf) * ieff_avg   !1. hyperbolic function that results in a steep decline in ignition efficiency with increasing area burned
!ieff = 0.8 * (1. - burnedf) / (1. + 100. * burnedf) * ieff_avg  !eq2
!ieff = 0.5 * (1. - burnedf) / (1. + 100. * burnedf) * ieff_avg  !eq3
!ieff = 0.4 * (1. - burnedf) / (1. + 100. * burnedf) * ieff_avg  !eq4
!ieff = 0.5 * (1. - burnedf) / (1. + 150. * burnedf) * ieff_avg  !eq5
!ieff = 0.4 * (1. - burnedf) / (1. + 150. * burnedf) * ieff_avg  !eq6
!ieff = 0.3 * (1. - burnedf) / (1. + 50. * burnedf) * ieff_avg  !eq7
 


!if(abs(input%lat) <= 60.d0) then
!   latscale = 1. / (4.16 + 2.16 * cos(3.d0 * input%lat * pir))
!else
!   latscale = 0.5
!end if

latscale = 1. 	!FLAG: we are now already using a climate file that holds the transformed strikes from the flashes (pre-processed)

!=========================================================================================================================

!if(.not. input%spinup .and. year > 115) then
!light = light  * 1.2 				! assuming an 80% detection efficiency for the ground sensors (Alaska groundstrike data)  
!end if    

!=========================================================================================================================

!if (met%prec < 10. .and. cumfires == 0 .and. burnedf == 0. .and. light*area_ha >= 1.) then
!if(met%prec < 10.) then 
  !nlig = light * latscale * ieff  ! 0.1  !total flashes * 20% cloud-to-ground * 10% ignition efficiency (ha-1 d-1)
  !nlig = 1./area_ha
!else
!  nlig = 0.
!end if

!--------------------------------------------------------------------------------
!human ignitions

if(FDI > (0.25)) then
   riskfact = 1. / (1.2172 * pi * FDI) * exp(-1. * (log(FDI) + 1.2963) ** 2 / 0.18)	!peak at FDI 0.25 
else
  riskfact = 1.
end if 

!number of people on the gridcell who are most active in maintaining the fire regime

!FLAG=====================

!PD(1) = 0.01

!FLAG=====================

if (PD(1) > 0.) then
  people = int(PD(1) * area * 1.e-6)  !hunter-gatherers
  group  = 1
else if (PD(2) > 0.) then
  people = int(PD(2) * area * 1.e-6)  !farmers
  group  = 2
else
  people = int(PD(3) * area * 1.e-6)  !pastoralists
  group  = 3
end if

if (people > 0) people = max(people / 10, 1)   !only every 10th person lights fire unless there are less than 10 people 

if (people > 0) then
  calchumanfire = .true.

  annburntarget = osv%annburntarget
    
else
  calchumanfire = .false.
end if

if(input%spinup .and. year < 800) then
  calchumanfire = .false.
end if

!human ignitions ends here

!--------------------------------------------------------------------------------
!part 2.2.2, fuel moisture content

!Nesterov index is now calculated in the weather generator and passed to this subroutine in met

!---
!moisture content of each fuel class individually

!omega = exp(-alpha * NI)

!------------------------------------------------
!=====new water balance approach=====

wet =  min(met%prec / 50.,1.)

dry = met%tmax * (met%tmax - met%tdew) * alpha * omega0

omega = min(max(omega0 - dry + wet,0.),1.)

omega0 = omega

!------------------------------------------------


!write(0,*)'omega',omega

!---
!mass of live grass

wlivegrass = sum(livefuel(:,1) * fpc_grid * dphen,mask=.not.pft%tree)

!---
!live grass fuel moisture content

omega_lg = max(0., b * omega_s1 - c)

!back calculate alpha_lg from omega_lg (don't really know why but this is the way it is done in orignal SPITFIRE)

!write(0,*)omega_lg,NI

if (omega_lg > 0. .and. NI > 0.) then 
  alpha_lg = -log(omega_lg) / NI
else
  alpha_lg = 0.
end if

!---
!moisture content weighted among live grass and 1-h fuel

wfinefuel = woi(1) + wlivegrass

!omega_nl = exp(-(alpha(1) * woi(1) + alpha_lg * wlivegrass) / wfinefuel * NI)

omega_nl = (omega(1) * woi(1) + omega_lg * (wlivegrass + SOM_surf)) / (wfinefuel + SOM_surf)  !FLAG in original spitfire code it is: omega(1) + omega_lg * wlivegrass / wfinefuel

!omega(1) = omega_nl

!---
!weighted average estimate of fuel relative moisture content

wo   = sum(woi(1:3))
wtot = wo + wlivegrass 

!write(0,*)'w mass',wtot,wo

rdf = wo / wtot             !dead fuel : total fuel
rlf = (wlivegrass) / wtot     !live fuel : total fuel

alpha_df = sum(alpha(1:3) * woi(1:3)) / wo    !mass weighted average alpha for dead fuel

caf = alpha_df * rdf + alpha_lg * rlf         !mass weighted average alpha for live and dead fuels combined




!------------------------------------------------
!=====new water balance approach=====

!omega_o = exp(-caf * NI)                      !integrated fuel moisture

!------------------------------------------------

dry_o = met%tmax * (met%tmax - met%tdew) * caf * omega_o0

omega_o = min(max(omega_o0 - dry_o + wet,0.),1.)

omega_o0 = omega_o

!if(.not. spinup .and. year==127) then
!  write(0,'(i3,6f14.7)')d, omega_o0,dry_o,wet,omega_o, met%tmax, met%tdew
!end if




!if (omega_o < 1.) write(0,'(3f10.4)')soilwater,snowpack,omega_o

!---
!moisture of extinction over dead fuel and live grass

!pftdeadfuel = sum(litter_ag_slow + litter_ag_fast)  !across all pfts

!m_e = sum(m_ex * (litter_ag_slow + litter_ag_fast) / pftdeadfuel)
!m_e = sum(m_ex * fpc_grid / sum(fpc_grid) * (litter_ag_slow + litter_ag_fast) / pftdeadfuel)        !FLAG MP: fpc-weight pft-dependent m_ex as well???
m_e = sum(woi * me_fc) / sum(woi)

me_avg = m_e * rdf + me_lf * rlf

me_nl = (me_fc(1) * woi(1) + me_lf * (wlivegrass + SOM_surf)) / (wfinefuel + SOM_surf)

if (me_avg == 0.) then
  write(*,*)'me_problem',netfuel,m_e , rdf , me_lf , rlf
  stop
end if

!---------------------------
!part 2.2.3, fire danger

if (grascover >= 0.6) then        !BURNING BEHAVIOR DOMINATED BY GRASS CHARACTERISTICS
  FDI = max(0.,(1. - omega_nl / me_nl))  !eqn. 8
else
  FDI = max(0.,(1. - omega_o / me_avg))  !eqn. 8
end if

!write(0,*)'FDI',d,omega_o,me_avg,FDI,NI,caf

!if (FDI == 0.) then
!  cumfires = 0        !extinguish all burning fires
!  write(*,'(a16,4i6,10f14.4)')'FDI_zero',year,i,d,cumfires,Ab,abarf,afire_frac,light*area_ha,nlig*area_ha, PD,omega_o/me_avg,NI
!   write(*,'(a16,4i6,29f14.4)') 'FDI_zero', year,i,d,cumfires,Ab,abarf,afire_frac,light*area_ha, nlig*area_ha,PD, area_ha, treecover, grascover, DT/1000., tfire/60., ROSfsurface*60./1000., ROSbsurface*60./1000., Uforward*60./1000.,  &
!                                 LB,LBtree,LBgrass,woi,omega_o,omega_o/me_avg,FDI,Isurface,slopefact, input%slope
!  return !no fire danger, so no burning
!end if






!---------------------------
!lightning ignitions

if (light * area_ha > 0.) then

!  ieff = FDI * (1. - burnedf) * 0.5  !constant 0.8 for the fact that not all of any landscape is flammable
  
  ieff = FDI * 1.0 * (1. - burnedf) / (1. + 25. * burnedf) * ieff_avg
  
  prob = real(ranu(met%rndst)) * rng1 + half  !random value from (0,1)

  if (ieff > prob) then
    nlig = 1. 
  else
    nlig = 0.
  end if
  
else

  nlig = 0.

end if

!---------------------------
!part 2.2.4, mean fire area (rate of spread)

!write(0,*)'go calcros'

! write(0,*)'me_avg',me_avg
! write(0,*)'omega_o',omega_o
! write(0,*)'uforward',uforward
! write(0,'(a,3f10.3)')'woi',woi(1:3)
! write(0,'(a,9f8.2)')'fpc',fpc_grid
!
! write(0,'(a,9f8.2)')'lf1',livefuel(:,1)
! write(0,'(a,9f8.2)')'lf2',livefuel(:,2)
! write(0,'(a,9f8.2)')'lf3',livefuel(:,3)
! write(0,'(a,9f8.2)')'lf4',livefuel(:,4)
!
! write(0,'(a,9f8.2)')'df1',deadfuel(:,1)
! write(0,'(a,9f8.2)')'df2',deadfuel(:,2)
! write(0,'(a,9f8.2)')'df3',deadfuel(:,3)
! write(0,'(a,9f8.2)')'df4',deadfuel(:,4)

!---------------------------
!characterize fuel that is carrying the fire in terms of mass, bulk density, surface area to volume ratio, and moisture content

!---
!fuel mass

totlfuelw = sum(lfuelw)
totdfuelw = sum(dfuelw)
totlfuelg = sum(lfuelg)
totdfuelg = sum(dfuelg)

pftlivefuel = sum(livefuel(:,1:3),dim=2)
pftdeadfuel = sum(deadfuel(:,1:3),dim=2)

where (pftdeadfuel > 0.)
  fuelpft = .true.
elsewhere
  fuelpft = .false.
end where

livemass = sum(pftlivefuel,mask=.not.pft%tree)
deadmass = sum(pftdeadfuel,mask=fuelpft)

!mean bulk density for each PFT across fuel classes

rhobPFT(8:9) = rho_livegrass !because the bulk density of grass litter should not be less than the bulk density of the live grass

where (fuelpft)
  rho_pft = rhobPFT * (deadfuel(:,1) + 1.2 * deadfuel(:,2) + 1.4 * deadfuel(:,3)) / sum(deadfuel(:,1:3),dim=2)
end where

!------------------------------------------------------------------------------------------------------------------
!Rate of spread calculations

ROSfcrown = 0.

wn       = totlfuelg + totdfuelg
rho_b    = rho_livegrass
sigma    = sigma_i(1) !typical value for the fire-carrying fuelbed (cf. Scott & Burgan, 2005)
relmoist = omega_nl / me_nl

if (relmoist < 1.) then
  
  gscale = -0.0848 * min(rho_livegrass,12.) + 1.0848    !scale factor for bulk density of grass fuel, reduces rate of spread in high bulk density grasses (e.g., tundra)

  ROSfsurface_g = (0.165 + 0.534 * Uforward / 60.) * exp(-0.108 * relmoist * 100.) * gscale * 60.     ! from Mell_etal2008, eqn. 2
  
!  write(*,'(a,2i5,6f14.4)') 'grass_ROS :', year, d, Uforward, relmoist, ROSfsurface_g, ROSfsurface_g/Uforward,osv%gdd20,rho_livegrass

else

  ROSfsurface_g = 0.

end if

!---
!surface ROS in woody litter

wn = livemass + deadmass

rho_b = (rho_livegrass * livemass + sum(rho_pft * pftdeadfuel)) / wn

!mass weighted average surface area to volume ratio across all fuel (sigma)

sigma = 5. !(sum(sigma_i * woi(1:3)) + sigma_i(1) * (livemass + SOM_surf)) / wn

relmoist = omega_o / me_avg

if (relmoist < 1.) then       !FDI not zero for this landscape component

  !write(0,*)'calcros',orgf*wn,rho_b,omega_o,relmoist

  call calcROS(orgf*wn,rho_b,sigma,omega_o,relmoist,Uforward,ROSfsurface_w)

else
  ROSfsurface_w = 0.
end if
  
!  write(*,'(a,2i4,6f14.7)') 'tree-dominated_ROS ', year,d, ROSfsurface, Uforward, rho_b, sigma, relmoist, wn
  
  if (treecover >= 0.6) then
  
    !---
    !crown ROS

    canopywater = sum(wscal * fpc_grid,mask=pft%tree) / treecover

    if (canopywater < 0.01) then   !calculate rate of spread in crown

  !    write(*,*) 'low canopywater! '

      wn = sum(livefuel(1:7,1))  !mass of 1-hr fuel in living biomass  !FLAG maximum leaf dry mass set to 8 kg m-2!

     ! write(0,*)'leafmass',wn

      wn = min(8000.,wn)

      rho_b = 0.1
      sigma = sigma_i(1)
      relmoist = 0.99  !at moisture content just below extinction fire can spread

      call calcROS(orgf*wn,rho_b,sigma,0.3,relmoist,Uforward,ROSfcrown)

    end if

  end if

  !---
  !calculate weighted average surface rate of spread on tree and grass fractions
  
  ROSfsurface = (ROSfsurface_w * treecover + ROSfsurface_g * grascover) / (treecover + grascover)
  
  !choose max from above calculations
  
!  if (ROSfcrown > ROSfsurface) then        !set definition of logical "crownfire" here.......
!    crownfire = .true.
!    ROSfsurface = ROSfcrown
    
!    write(*,'(a,2i4,2f7.2)') 'crownfire!', year, d, input%lon, input%lat    
!  else
!    crownfire = .false.
!  end if
  


!---------------------------
!backward rate of spread - decreases with stronger wind

ROSbsurface = ROSfsurface * exp(-0.012 * Uforward)

!write(0,*)'done ros,',ROSfsurface,ROSbsurface

!---
!length-to-breadth ratio of burn ellipse

if (Uforward >= 16.67) then        ! windspeed exceeds 1 km/h

  Ukmh = 0.06 * Uforward  !convert to km h-1 for the following equations

  LBtree  = 1.  + 8.729 * (1. - exp(-0.03 * Ukmh))**2.155  !eqn. 12

  LBgrass = 1.1 + Ukmh**0.0464    !eqn. 13 as in the f77 code

  !LBgrass = 1.1 + Ukmh**0.464     !eqn. 13 as published in the paper NB difference in the exponent

  if (totvcover > 0.) then
    !LB = (LBtree * treecover + LBgrass * grascover) / totvcover !weighted average
    LB = (LBtree * treecover + LBgrass * grascover) !weighted average
  else
    LB = LBgrass
  end if

else

  LB = 1.

end if

LB = min(LB,8.)  !limit LB to a maximum of 8 as in f77 code

!---
!fire duration

tfire = 241. / (1. + 240. * exp(-11.06 * FDI))  !eqn. 14 (minutes)

!---
!total distance traveled

DT = tfire * (ROSfsurface + ROSbsurface)

!---
!mean fire area

avg_cont_area = max(avg_cont_area,10.)                ! assumption is that the smallest possible size of a kernel is 1 ha: nothing will be fractionated to a size less than 1 ha

abarf = (pi / (4. * LB) * DT**2) * 0.0001 * slopefact  !average size of an individual fire (eqn. 11) (ha)

!write(0,*) 'abarf', abarf, LB, DT, slopefact, avg_cont_area

abarf = min(abarf,avg_cont_area)                ! the size of an individual fire is not allowed to be greater than the average contiguous patch size

!---------------------------
!human-caused fires

!before recalculating uniquefires (caused by humans), reduce the cumulative number of currently burning human-caused fires by
!merging all fires caused by one single person down into one single fire (since the number of unique fires caused in one day 
!is considered to be spaced at the minimum distance between...

cumfires = max(cumfires - nhig + arsonists,0)

!total area and fraction of gridcell burned

if (calchumanfire .and. abarf > 0. .and. abarf < 100.) then  !avoid starting very large fires

  !human burning
  !calculate number of human ignitions based on mean fire size and total number of people on the gridcell

  !potential number of unique fires that can be started by one person in one day given the conditions

  firewidth = DT / LB

  uniquefires = walkdist / firewidth
  
  !potential burned area of one person on this day

  humfire1 = uniquefires * abarf
  
  !target area of burning for this day (ha)

  dayburntarget = area_ha * max(annburntarget(group) -  burnedf20 - afire_frac,0.)
  
  !number of human ignitions is the lesser of the burn target and the number of fires people can cause
  
  if (people * humfire1 >= dayburntarget) then

    !there are more than enough people around to start enough fires large enough to reach the burn target

    nhig = int(dayburntarget / abarf)
    
    arsonists = max(nhig / uniquefires,1)  !the number of people who actually caused fires

  else
     
    nhig = uniquefires * people   !per gridcell on this day

    arsonists = people
 
  end if
  
  if (nhig == 0) arsonists = 0
  
  !write(0,'(3i12,3f12.4)')people,uniquefires,nhig,abarf,burnedf20,afire_frac
  
else

  nhig = 0
  arsonists = 0
  uniquefires = 0

end if

nhig = nhig * riskfact			  ! reduce number of fires caused when FDI gets above 0.25

!numfires_nat = int(FDI * nlig * area_ha)  !lightning fires started on this day
numfires_nat = int(nlig)! * area_ha)  !lightning fires started on this day
numfires_hum = nhig                       !human fires started on this day


numfires = numfires_nat + numfires_hum

!write(*,'(a,3i5,4f14.8)')'potential burnday',year,i,d,FDI,NI,light*area_ha,nlig*area_ha

!cumfires = cumfires + numfires             !accumulate all fire since last time FDI was not zero
cumfires = cumfires + numfires - nint(burnedf * real(cumfires + numfires))          ! allow naturally-caused fires to carry over from one day to the next

!if (calchumanfire .and. abarf > 0. .and. abarf < 100.) then  !avoid starting very large fires
!  write(*,'(i6,f10.4,6i10,2f10.6,f10.1,f10.6)')d,firewidth, uniquefires, numfires_nat, numfires_hum, cumfires, arsonists,people, PD(1),afire_frac,totburn,abarf
!end if



unburneda = area_ha - totburn

unburneda = max(0.,unburneda)

!Ab = max(0.,min(unburneda-0.05*area_ha,cumfires * abarf))  !ha
Ab = max(0.,min(unburneda,cumfires * abarf))  !ha



!reduce cumfires by nhig (human caused fires only last for one day)

!cumfires = max(0,cumfires - nhig)




!Ab = numfires * abarf  !ha

!write(0,*)year,i,d,area_ha,Ab,unburneda

Abfrac = Ab / area_ha

!if (d == 365) then
!  write(*,'(a13,4i6,2f10.1,f12.6,5f10.1,f10.5)')'burn_day',year,i,d,numfires,abarf,Ab,Abfrac,totburn,unburneda,PD,nhig*1.e6,ROSfsurface,FDI
!end if

!if (year == 200) then
  !write(0,'(a,4i6,2f10.1,f12.6,2f10.1,2f10.5)')'flag2',year,i,d,numfires,abarf,Ab,Abfrac,totburn,unburneda,dphen(9),omega_s1
  !write(*,'(a,3i6,6f10.1)')'flag2',year,i,d,abarf,LB,DT,tfire,ROSfsurface,ROSbsurface
!end if

if (Abfrac > 1. .or. Abfrac < 0. .or. unburneda < 0.) then
  write(0,*)'ABfrac problem'
  write(0,*)FDI,abarf,Abfrac,Ab,unburneda,numfires,totburn,area_ha,input%lon,input%lat
  stop
end if

!total area burned (ha)

!Ab = Abfrac * area_ha        !eqn. 2

if (Ab == 0.) then 
  if(bavard) write(*,'(a16,4i6,16f14.4)') 'Area_burned_zero ', year, i,d,cumfires,Ab,abarf,afire_frac,light*area_ha,nlig,FDI, met%prec, NI, grascover, omega_o, me_avg, omega_nl, me_nl, PD
!  write(*,'(a16,4i6,29f14.4)') 'Area_burned_zero', year,i,d,cumfires,Ab,abarf,afire_frac,light*area_ha, nlig*area_ha,PD, area_ha, treecover, grascover, DT/1000., tfire/60., ROSfsurface*60./1000., ROSbsurface*60./1000., Uforward*60./1000.,  &
!                                LB,LBtree,LBgrass,woi,omega_o,omega_o/me_avg,FDI,Isurface,slopefact, input%slope   
  return
end if 

!write(0,*)'area burned',abfrac,ab

!----------------------------------------------
!part 2.2.4, fractional combustion of dead fuel

!---
!fraction of live grass consumed in surface fire

rm = omega_lg / me_lf   !FLAG in spitfire code the division by moisture of extinction is not included

if (rm <= 0.18) then
  CFlg = 1.
else if (rm > 0.73) then
  CFlg = 2.45 - 2.45 * rm
else
  CFlg = 1.10 - 0.62 * rm   !coefficients from code
end if

!---
!fraction of 1-h fuel consumed in surface fire

!rm = omega_nl / me_avg  !0.22  !m_e  !rm: ratio of moisture content to moisture of extinction

rm = omega(1) / me_fc(1)

if (rm <= 0.18) then
  CF(1) = 1.
else if (rm > 0.73) then
  CF(1) = 2.45 - 2.45 * rm
else
  !CF(1) = 1.20 - 0.62 * rm  !coefficients from paper
  CF(1) = 1.10 - 0.62 * rm   !coefficients from code
end if

!---
!fraction of 10-h fuel consumed in surface fire

!rm = omega(2) / 0.43  !m_e
rm = omega_o / me_avg
rm = omega(2) / me_fc(2)

if (rm <= 0.12) then
  CF(2) = 1.
else if (rm > 0.51) then
  CF(2) = 1.47 - 1.47 * rm
else
  CF(2) = 1.09 - 0.72 * rm
end if

!---
!fraction of 100-h fuel consumed in surface fire

!rm = omega(3) / 0.52  !m_e
rm = omega(3) / me_fc(3)

if (rm <= 0.38) then
  CF(3) = 0.98 - 0.85 * rm
else
  CF(3) = 1.06 - 1.06 * rm
end if

!---
!fraction of 1000-h fuel consumed in surface fire

!rm = omega(4) / 0.52  !m_e
rm = omega(4) / me_fc(4)

CF(4) = -0.8 * rm + 0.8

!---
!correct for spurious values of consumed fraction

CF = min(max(CF,0.),1.)

!---
!total fuel consumed (g m-2), mineral fraction subtracted

FC = CF * woi * (1. - ST)

!---
!surface fire intensity (kW m-1)

Isurface = h * sum(FC(1:3)) * ROSfsurface * min2sec  !eqn. 15

!write(0,*)'Isurf',isurface

if (Isurface < 50.) then  !ignitions are extinguished and there is no fire on this day
  if(bavard)  write(*,'(a16,4i6,19f14.4)')'Isurface_low',year,i,d,cumfires,Ab,abarf,afire_frac,light*area_ha,nlig,FDI,Isurface, met%prec, NI, grascover, omega_o, me_avg, omega_nl, me_nl, PD
! write(*,'(a16,4i6,29f14.4)') 'Isurface_low', year,i,d,cumfires,Ab,abarf,afire_frac,light*area_ha, nlig*area_ha,PD, area_ha, treecover, grascover, DT/1000., tfire/60., ROSfsurface*60./1000., ROSbsurface*60./1000., Uforward*60./1000.,  &
!                               LB,LBtree,LBgrass,woi,omega_o,omega_o/me_avg,FDI,Isurface,slopefact, input%slope
  Ab     = 0.
  Abfrac = 0.
  return       !leave the subroutine
end if

!---
!if there is fire today, update the fractional area burned and the live and dead biomass pools

afire_frac = afire_frac + Abfrac

!write(*,'(a16,4i6,29f14.4)') 'BURNDAY', year,i,d,cumfires,Ab,abarf,afire_frac,light*area_ha, nlig*area_ha,PD, area_ha, treecover, grascover, DT/1000., tfire/60., ROSfsurface*60./1000., ROSbsurface*60./1000., Uforward*60./1000.,  &
!                                LB,LBtree,LBgrass,woi,omega_o,omega_o/me_avg,FDI,Isurface,slopefact,input%slope
if(bavard) write(*,'(a16,4i10,23f14.3)') 'BURNDAY', year,i,d,cumfires,AB,abarf,afire_frac,light*area_ha,nlig,FDI,woi,omega_o,omega_o/me_avg,Isurface, met%prec, NI, grascover, omega_o, me_avg, omega_nl, me_nl, PD

!-------------------------------------------
!update litter pools to remove burned litter

!write(0,'(a,8f10.3)')'consumed fraction: ',CF,omega

!BBdead(:,1) = CF(1) * litter_ag_fast          !1-hr fast litter (leaves and grass)

!BBdead(:,2) = CF(1) * litter_ag_slow * 0.045  !1-hr slow litter (twigs)
!BBdead(:,3) = CF(2) * litter_ag_slow * 0.075  !10-hr
!BBdead(:,4) = CF(3) * litter_ag_slow * 0.21   !100-hr
!BBdead(:,5) = CF(4) * litter_ag_slow * 0.67   !1000-hr

BBdead(:,1) = Abfrac * CF(1) * litter_ag_fast          !1-hr fast litter (leaves and grass)

BBdead(:,2) = Abfrac * CF(1) * litter_ag_slow * 0.045  !1-hr slow litter (twigs)
BBdead(:,3) = Abfrac * CF(2) * litter_ag_slow * 0.075  !10-hr
BBdead(:,4) = Abfrac * CF(3) * litter_ag_slow * 0.21   !100-hr
BBdead(:,5) = Abfrac * CF(4) * litter_ag_slow * 0.67   !1000-hr

annBBdead = annBBdead + BBdead

if (d == 0) then
  write(0,*)'last day spitfire'

  write(0,'(a,9f10.1)')'1-hr grass ',annBBdead(:,1)
  write(0,'(a,9f10.1)')'1-hr wood  ',annBBdead(:,2)
  write(0,'(a,9f10.1)')'10-hr w    ',annBBdead(:,3)
  write(0,'(a,9f10.1)')'100-hr w   ',annBBdead(:,4)
  write(0,'(a,9f10.1)')'1000-hr w  ',annBBdead(:,5)

end if
 

!-----------------------
!2.2.6 fire damage to living plants

!---
!fire residence time (min), Peterson & Ryan 1986, eqn. 8, used conversion factor because original eqn. needs fuel load in g cm-2

!write(0,'(a,i5,f12.5,f9.1)')'Abfrac',d,Abfrac,litter_ag_fast(9)
!write(0,'(a,4f10.2)')'dead fuel   ',woi
!write(0,'(a,4f10.2)')'dead burnedf',CF

tau_l = 39.4 * sum(woi * (1. - sqrt(1. - CF))) * 0.0001

!write(0,'(a,f10.2)')'res time',tau_l

!write(0,*)'flag 3',Isurface

!----------------------------------------------
!tree mortality due to crown and cambial damage

CK  = 0.
!P_m = 0.

do l = 1,npft
  if (present(l) .and. pft(l)%tree) then
    call firemortality(l,ind(:,l),Isurface,tau_l,CK(l),P_m(l),dphen(l))
  end if
end do

!CK  = 0.
!P_m = 0.


!write(*,'(a,3i5,8f10.3)')'crown CK ',year, d, i, CK(1:7), ROSfsurface
!write(*,'(a,3i5,8f10.3)')'Pm total ',year, d, i, P_m(1:7), ROSfsurface

!---
!calculate combusted live biomass (determined by CK)
!the equations below reflect removal of 100% of leaves and 1-h wood and 5% of 10-h wood as written on the next line
!sm_ind = sm_ind - (CK * 0.045 * sm_ind + CK * 0.05 * 0.075 * sm_ind)
!NB this relationship is from Appendix B, and does not seem to be the same as what is written on
!pg 1997 after eqn. 17


if (crownfire) ann_kill = ann_kill + Abfrac   !FLAG! test

where (present)  !over PFTs
  where (pft%tree)

    BBlive(:,1) = Abfrac * CK * lm_ind * nind
    BBlive(:,2) = Abfrac * CK * sm_ind * nind * 0.04875
    BBlive(:,3) = Abfrac * CK * hm_ind * nind * 0.04875

    !---
    !calculate running sum of mortality probability times burned fraction
        
    ann_kill = P_m * Abfrac + ann_kill

  elsewhere !grass

    !update live grass biomass, burning all grass at the same rate as 1-h dead fuel 
    !the below equation works because the nind of grass is always 1, so to remove biomass we reduce the lm_ind

    BBlive(:,1) = Abfrac * CF(1) * lm_ind

  end where
end where

annBBlive = BBlive + annBBlive

!---
!carbon flux accounting

BBtot = sum(BBdead) + sum(BBlive)  !total C emissions from burning, across all PFTs

acflux_fire = acflux_fire + BBtot

!trace gas emissions

BBpft = 0.001 * c2om * (sum(BBlive,dim=2) + sum(BBdead,dim=2))  !units kg (dry matter) m-2

!do l = 1,nspec
!  Mx(l) = sum(emfact(:,l) * BBpft)  !units g (species) m-2
!end do

!Mx = 0.

!-----------------------

sumfdi   = sumfdi   + FDI        !FDI summed up over the course of a year
burndays = burndays + 1
totburn  = totburn  + Ab

!aMx(1) = aMx(1) + Mx(1) ! CO2
!aMx(2) = aMx(2) + Mx(2) ! CO
!aMx(3) = aMx(3) + Mx(3) ! CH4
!aMx(4) = aMx(4) + Mx(4) ! VOC
!aMx(5) = aMx(5) + Mx(5) ! TPM
!aMx(6) = aMx(6) + Mx(6) ! NOx               

aMx(1) = aMx(1) + (sum(emCO2 * BBpft)) ! CO2
aMx(2) = aMx(2) + (sum(emCO * BBpft)) ! CO
aMx(3) = aMx(3) + (sum(emCH4 * BBpft)) ! CH4
aMx(4) = aMx(4) + (sum(emVOC * BBpft)) ! VOC
aMx(5) = aMx(5) + (sum(emTPM * BBpft)) ! TPM
aMx(6) = aMx(6) + (sum(emNOx * BBpft)) ! NOx           



!aMx(4) = anumfires / (area * 1e-6)    ! annual number of fires within a gridcell, divided by the area of the gridcell in km2 => annual number of fires per km2
!aMx(5) = sumfdi !/ burndays
!aMx(6) = totburn                     ! total annual area burned in the gridcell

!if (year == 919 .and. i == 1) then
  !write(0,'(a,3i6,3f10.1,f12.6,4f10.1)')'flag2',year,i,d,abarf,Ab,numfires,Abfrac,totburn,unburneda,PD,nhig
  !write(0,'(i10,f10.1,4f10.1,4f10.4,f10.4,f10.4,f10.1,f10.1,4f10.4,f10.1)')d+365*(year-1),NI,woi,omega,omega_o,FDI,ROSfsurface,Isurface,CF,Abarf
!end if

991 format(a,3i6,3f10.1,f12.6,4f10.1)
992 format(2i6,l4,9f6.2)
993 format(3i6,3f10.1,f12.6,2f10.1,l4,f6.2)

!write(*,993)year,i,d,abarf,Ab,numfires,Abfrac,totburn,unburneda,crownfire,dphen(2)

end subroutine spitfire

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine calcROS(wn,rho_b,sigma,omega_o,relmoist,Uforward,rateofspread)

!arguments

real(sp), intent(in)  :: wn            !total fuel mass of the organic part of the fuel (mineral fraction subtracted) (g DM m-2)
real(sp), intent(in)  :: rho_b         !bulk density of the fuel (kg m-3)
real(sp), intent(in)  :: sigma         !surface area to volume ratio of the fuel (cm2 / cm3)
real(sp), intent(in)  :: omega_o       !relative moisture content of the fuel (unitless; 0=completely dry fuel)
real(sp), intent(in)  :: relmoist      !relative moisture content of the fuel relative to its moisture of extinction (unitless; omega_o / m_e)
real(sp), intent(in)  :: Uforward      !windspeed (m min-1)

real(sp), intent(out) :: rateofspread  !fire rate of spread (m min-1)

!parameters

real(sp), parameter :: rho_p = 513.       !oven-dry particle density (kg m-3)
real(sp), parameter :: nu_s  =   0.41739  !mineral dampening coefficient (unitless)

!local variables

real(sp) :: A                   ! empirical function of sigma, between approx. 0.2 to 1.0, reduces sensitivity of the ratio gammaprime / gammaprimemax for large values of sigma
real(sp) :: B
real(sp) :: C
real(sp) :: E
real(sp) :: Beta           !packing ratio (unitless) (fuel bulk density / oven dry particle density)
real(sp) :: Beta_op        !optimum packing ratio (unitless)
real(sp) :: Gammaprime     !potential reaction velocity; reaction velocity that would occur if fuel were free of moisture and mineral content; = Gamma / nu_M / nu_s 
real(sp) :: Gammaprimemax  !maximum reaction velocity; would occur if fuel elements in the fuelbed were arranged in the most efficient way possible (min-1)
real(sp) :: IR             !reaction intensity (kJ m-2 min-1)
real(sp) :: Phiw           !wind coefficient (?)
real(sp) :: pratio         !ratio of packing ratio to optimum packing ratio (unitless)
real(sp) :: Qig            !heat of pre-ignition (kJ kg-1)
real(sp) :: epsilon        !effective heating number; proportion of fuel that is heated before ignition occurs (dimensionless), between 0 and 1 (1 at infinitive surface to volume ratio))
real(sp) :: nu_M           !moisture dampening coefficient (unitless)
real(sp) :: xi             !propagating flux ratio; proportion of IR transferred to unburned fuels (dimensionless); in reality influenced by convection, radiation, flame contact and ignition-point transfer
real(sp) :: windfact       !high wind multiplier for rate of spread (unitless)
real(sp) :: Ums            !wind speed (m s-1)

!-------------------------------
!Terms used in several equations

!---
!packing ratio (unitless)

Beta = rho_b / rho_p

!---
!optimum packing ratio (unitless)

Beta_op = 0.200395 * sigma**(-0.8189)

pratio = Beta / Beta_op

!-------------------------------
!Term 1: reaction intensity 

!---
!maximum reaction velocity (min-1)

Gammaprimemax = 1. / (0.0591 + 2.926 * sigma**(-1.5))

!---
!optimum reaction velocity (min-1)

A = 8.9033 * sigma**(-0.7913)

Gammaprime = Gammaprimemax * (pratio)**A * exp(A * (1. - pratio))

!---
!moisture dampening coefficient

nu_M = 1. - 2.59 * relmoist + 5.11 * relmoist**2 - 3.52 * relmoist**3

!---
!reaction intensity

IR = Gammaprime * wn * h * nu_M * nu_s  !eqn. A1 (kJ m-2 min-1)

!-------------------------------
!Term 2: 'ratio of propagating flux to reaction intensity' (fraction)

xi = exp(0.792 + 3.7597 * sqrt(sigma) * (Beta + 0.1)) / (192. + 7.9095 * sigma)

!-------------------------------
!Term 3: wind coefficient (unitless?)

C = 7.47 * exp(-0.8711 * sigma**0.55)

B = 0.15988 * sigma**0.54

E = 0.715 * exp(-0.01094 * sigma)

Phiw = C * (3.281 * Uforward)**B * pratio**(-E)

!-------------------------------
!Term 4: effective heating number (fraction)

epsilon = exp(-4.528 / sigma)

!-------------------------------
!Term 5. heat of pre-ignition (kJ kg-1)

Qig = 581. + 2594. * omega_o

!-------------------------------
!Term 6. wind multiplier for high wind conditions
!at a windspeed of 10 m s-2 and above, the calculated ROS will be doubled,
!as the BEHAVE-based ROS is increasingly too low at higher wind speeds (see Morvan et al., 2008, Figure 13)

Ums = Uforward / 60.

if(Ums <= 10.) then
  windfact = 1. + exp(2. * Ums - 20.)
else
  windfact = 2.
end if    

!-------------------------------
!Rate of spread (m min-1)

rateofspread = IR * xi * (1. + Phiw) / (rho_b * epsilon * Qig) * windfact

if(rateofspread < 0.) then
  write(0,'(a,7f14.7)') 'end calcROS: ', IR, xi, 1. + Phiw, rho_b, epsilon, Qig, windfact
  write(0,'(a,5f14.7)') 'IR components: ', Gammaprime, wn, h, nu_M, relmoist
end if

end subroutine calcROS

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine firemortality(pft,ind,Isurface,tau_l,CKbar,Pmbar,dphen)

use parametersmod, only : npft,nhclass,pftpar
use individualmod, only : sizeind

implicit none

!calculate mortality by fire due to crown scorch and cambial damage
!PFTs are divided up into height classes in the individualmod module

!arguments

integer,  intent(in) :: pft         !pft number
type(sizeind), dimension(:), intent(in) :: ind       !typical allometry across the height classes (vector structure)

real(sp), intent(in) :: Isurface   !surface fire intensity (kW m-1)
real(sp), intent(in) :: tau_l      !fire residence time (min)
real(sp), intent(in) :: dphen      !phenological state (1=full leaf; 0=no leaves)

real(sp), intent(out) :: CKbar     !mean fraction of crown scorch (fraction)
real(sp), intent(out) :: Pmbar     !mean total probability of mortality

!parameters


real(sp), dimension(npft) :: F  !scorch height parameters

real(sp), dimension(npft) :: RCK 
real(sp), dimension(npft) :: p   

!local variables

real(sp) :: SH         !scorch height

real(sp), dimension(nhclass) :: CK        !fraction of crown scorch (fraction)                                        
real(sp), dimension(nhclass) :: tau_c     !critical time for cambial damage (min)                                     
real(sp), dimension(nhclass) :: tau_r     !ratio of fire residence time to critical time for cambial damage (ratio)   
real(sp), dimension(nhclass) :: P_mtau    !probability of mortality due to cambial damage                             
real(sp), dimension(nhclass) :: P_mCK     !probability of mortality due to crown damage                               
real(sp), dimension(nhclass) :: P_m       !total probability of mortality                                             

!--------------------------
! Enregistrement des parametres
F   = pftpar(:,41)
RCK = pftpar(:,42)
p   = pftpar(:,43)

!crown kill by height class

SH = F(pft) * Isurface**0.667       !scorch height (m)

CK = (SH - ind%height + ind%lcrown) / ind%lcrown      !proportion of the crown affected by fire (combusted)

CK = max(min(CK,1.),0.)            !keep the proportion between 0 and 1

P_mCK = RCK(pft) * CK**p(pft) * dphen     !Eqn. 22 probability of mortality due to crown damage

!--------------------------
!cambial kill by height class

tau_c = 2.9 * ind%barkt**2                !critical time for cambial damage (min)

tau_r = tau_l / tau_c              !tau ratio: fire residence time / critical time for cambial damage

where (tau_r >= 2.)
  P_mtau = 1.
elsewhere
  P_mtau = 0.563 * tau_r - 0.125
end where

P_mtau = max(0.,P_mtau)  !do not allow P_mtau to be less than zero

!--------------------------
!Eqn. 18 total probability of mortality (per woody PFT and height class)

P_m = P_mtau + P_mCK - P_mtau * P_mCK

!---
!calculate scalar averages of CK and P_m over all height classes

CKbar = sum(CK  * ind%massf)
Pmbar = sum(P_m * ind%massf)

end subroutine firemortality

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine burnedbiomass(i,j,osv,Abfrac,year)

use parametersmod,   only : pft
use mpistatevarsmod, only : statevars

!several things happen to biomass as a result of fire
!1. dead litter is combusted
!2. live biomass is killed and added to dead litter
!3. live biomass is directly combusted

implicit none

!arguments

integer,  intent(in) :: i  !tile index
integer,  intent(in) :: j  !gridcell index
integer,  intent(in) :: year
real(sp), intent(in) :: Abfrac

type(statevars), target, intent(inout) :: osv

!pointers

logical,  pointer, dimension(:) :: present         !PFT is present

real(sp), pointer, dimension(:) :: litter_ag_fast  !aboveground litter, fast turnover (g C m-2)
real(sp), pointer, dimension(:) :: litter_ag_slow  !aboveground litter, slow turnover (g C m-2)
real(sp), pointer, dimension(:) :: litter_bg       !belowground litter (g C m-2)

real(sp), pointer, dimension(:) :: nind            !gridcell individual density (indiv/m2)
real(sp), pointer, dimension(:) :: lm_ind          !leaf mass, average individual (g C)
real(sp), pointer, dimension(:) :: sm_ind          !sapwood mass, average individual (g C)
real(sp), pointer, dimension(:) :: hm_ind          !heartwood mass, average individual (g C)
real(sp), pointer, dimension(:) :: rm_ind          !root mass, average individual (g C)

!local variables

real(sp), dimension(npft)   :: nind_kill !fraction of PFT killed by fire (per pft)
integer :: a

!-----------------------------------

present => osv%tile(i)%present

litter_ag_fast => osv%tile(i)%litter_ag_fast(:,1)  !dead biomass
litter_ag_slow => osv%tile(i)%litter_ag_slow(:,1)
litter_bg      => osv%tile(i)%litter_bg(:,1)

nind           => osv%tile(i)%nind                 !individual density
lm_ind         => osv%tile(i)%lm_ind(:,1)          !live biomass state variables
sm_ind         => osv%tile(i)%sm_ind(:,1)
hm_ind         => osv%tile(i)%hm_ind(:,1)
rm_ind         => osv%tile(i)%rm_ind(:,1)

!-----------------------------------

!write(0,*)'===end of year biomass removal:', year, '==='

!write(0,*)'step 0a'
!
!write(0,'(a,9f10.1)')'1-hr grass ',annBBdead(:,1)
!write(0,'(a,9f10.1)')'1-hr wood  ',annBBdead(:,2)
!write(0,'(a,9f10.1)')'10-hr w    ',annBBdead(:,3)
!write(0,'(a,9f10.1)')'100-hr w   ',annBBdead(:,4)
!write(0,'(a,9f10.1)')'1000-hr w  ',annBBdead(:,5)
!
!!annBBdead = annBBdead * Abfrac
!
!write(0,*)'step 0b'
!
!write(0,'(a,9f10.1)')'1-hr grass ',annBBdead(:,1)
!write(0,'(a,9f10.1)')'1-hr wood  ',annBBdead(:,2)
!write(0,'(a,9f10.1)')'10-hr w    ',annBBdead(:,3)
!write(0,'(a,9f10.1)')'100-hr w   ',annBBdead(:,4)
!write(0,'(a,9f10.1)')'1000-hr w  ',annBBdead(:,5)
!
!write(0,*)
!write(0,*)'step 1'
!write(0,'(a,9f10.1)')'fast',litter_ag_fast
!write(0,'(a,9f10.1)')'slow',litter_ag_slow

!write(0,'(9f10.1)')annBBlive(:,1)
!write(0,'(a,9f10.4)')'kill',ann_kill(:)

!1. consumption of dead biomass
!calculations across PFT vector

litter_ag_fast = litter_ag_fast - annBBdead(:,1)

litter_ag_slow = litter_ag_slow - annBBdead(:,2)
litter_ag_slow = litter_ag_slow - annBBdead(:,3)
litter_ag_slow = litter_ag_slow - annBBdead(:,4)
litter_ag_slow = litter_ag_slow - annBBdead(:,5)

litter_ag_fast = max(litter_ag_fast,0.)
litter_ag_slow = max(litter_ag_slow,0.)

!2. consumption of live biomass

!direct combustion of live biomass is not handled here for the moment
!in the main spitfire subroutine in the BBlive term which goes into the
!calculation of acflux_fire

!write(0,*)'step 2'
!write(0,'(a,9f10.1)')'fast',litter_ag_fast
!write(0,'(a,9f10.1)')'slow',litter_ag_slow


!3. transfer of killed but not consumed live biomass to litter

!ann_kill = 0.

where (present)  !over PFTs
  where (pft%tree)

    nind_kill = ann_kill * nind  !ann kill is the annual sum P_m * Abfrac

    litter_ag_fast = litter_ag_fast + nind_kill *  lm_ind
    litter_ag_slow = litter_ag_slow + nind_kill * (sm_ind + hm_ind)
    litter_bg      = litter_bg      + nind_kill *  rm_ind

    nind = max(0.,nind - nind_kill)

  elsewhere !grass
    
    !note we don't remove the live grass because the time required to rebuild grass
    !biomass in this version of LPJ is such that it will be unrealistically low
    !in years subsequent to a fire and there is no way to manipulate nind with grass

    !lm_ind = lm_ind - annBBlive(:,1)

  end where
end where

do a = 1, npft
  if(nind(a) == 0.) then
    present(a) = .false.
    lm_ind(a) = 0.
    sm_ind(a) = 0.
    hm_ind(a) = 0.
    rm_ind(a) = 0.
  end if
end do    

!write(0,*)'step 3'
!write(0,'(a,9f10.1)')'fast',litter_ag_fast
!write(0,'(a,9f10.1)')'slow',litter_ag_slow


!correct for mathematical slop

!if (any(litter_ag_fast < 0.) .or. any(litter_ag_slow < 0.)) then
!  write(0,*)'negative litter spitfire 2!'
!  write(0,'(18f10.2)')litter_ag_fast,litter_ag_slow
!  stop
!end if

!litter_ag_fast = max(litter_ag_fast,0.)
!litter_ag_slow = max(litter_ag_slow,0.)

end subroutine burnedbiomass

!-----------------------------------------------------------------------------------------------------------------------------------------------------

end module spitfiremod
