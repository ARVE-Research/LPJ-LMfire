module lpjmod

use parametersmod, only : stdout,stderr

implicit none

public :: lpjcore

contains

!------------------------------------------------------------------------------------------------------------

subroutine lpjcore(in,osv)

use parametersmod,    only : sp,dp,npft,ncvar,ndaymonth,midday,pftpar,pft, &
                             lm_sapl,sm_sapl,rm_sapl,hm_sapl,sla,          &
                             allom1,allom2,allom3,latosa,wooddens,reinickerp,lutype,climbuf,nhclass
use mpistatevarsmod,  only : inputdata,statevars
use weathergenmod,    only : metvars_in,metvars_out,rmsmooth,weathergen_driver,daily
use radiationmod,     only : elev_corr,calcPjj,radpet
use alccmod,          only : harvest,alcc,tile_landuse
use bioclimmod,       only : climate20,bioclim
use spitfiremod,      only : spitfire,burnedbiomass,managedburn
use snowmod,          only : snow
use hetrespmod,       only : littersom2,hetresp
use lightmod,         only : light
use simplesoilmod,    only : simplesoil
use nppmod,           only : calcnpp
use gppmod,           only : calcgpp
use mortalitymod,     only : mortality
use establishmentmod, only : establishment
use allocationmod,    only : allocation
use turnovermod,      only : turnover
use killplantmod,     only : killplant
use soiltemperaturemod, only : soiltemp
use foragersmod,      only : foragers,popgrowth,simpleforagers
use individualmod,    only : sizeind,allomind

use landscape_geometrymod, only : landscape_fractality

!use lpjstatevarsmod,  only : gsv,sv,ov  !TEMPORARY

implicit none

!argument

type(inputdata), intent(in)    :: in
type(statevars), target, intent(inout) :: osv  !state variables sent back out with MPI

integer :: ntiles
logical :: spinup
integer :: year
  
!---------------------------------------------------------------
!array location

integer :: i,j
!integer :: x,y
integer :: a ,b

!counters

integer :: m
integer :: d
integer :: dm
integer :: dyr
!integer :: yesterday

!----------------------------
!local variables
!integer :: dry
!integer :: maxdry

!type(metvars_in)  :: met_in
type(metvars_out), dimension(365) :: met_out

real(sp) :: tcm
real(sp) :: Pjj
real(sp) :: Ratm
real(sp) :: co2
real(sp) :: burnedf20
real(sp) :: forager_pd20

real(sp), dimension(12) :: temp0  !mean monthly temperature of the previous year (degC)
real(sp), dimension(12) :: prec   !total monthly precip (mm)
real(sp), dimension(12) :: temp   !mean monthly temperature (degC)
real(sp), dimension(12) :: tmin   !mean minimum monthly temperature (degC)
real(sp), dimension(12) :: tmax   !mean maximum monthly temperature (degC)
real(sp), dimension(12) :: cldf   !cloud cover fraction
real(sp), dimension(12) :: wetd   !
real(sp), dimension(12) :: lght   !
real(sp), dimension(12) :: wind   !

real(sp), dimension(12) :: mpet
real(sp), dimension(12) :: mdayl
real(sp), dimension(12) :: mpar_day

real(sp), dimension(365) :: dtemps  !smooth interpolation daily temperature
real(sp), dimension(365) :: dtemp
real(sp), dimension(365) :: dprec
real(sp), dimension(365) :: dtmn
real(sp), dimension(365) :: dtmx
real(sp), dimension(365) :: dcld
real(sp), dimension(365) :: dlght
real(sp), dimension(365) :: dwind

!real(sp), dimension(365) :: ddayl
real(sp), dimension(365) :: dpet
!real(sp), dimension(365) :: dswrad

real(sp), parameter :: fpc_tree_max = 0.95

logical, parameter :: dospitfire = .true.

!---------------------------------------------------------------
!TEMPORARY DECLARATIONS WE SHOULD TRY TO organize LATER

real(sp) :: recoverf  !fraction of the gridcell that is recovering from land use (difference)
real(sp) :: treefrac

real(sp), dimension(7) :: soilpar

!real, dimension(npft) :: gpp_temp
!real, dimension(npft) :: npp_temp

logical :: lucc

!annual gridcell state variables

real(sp) :: arunoff       !total annual runoff (mm)
real(sp) :: arunoff_drain !annual drainage (layer 2 runoff) (mm)
real(sp) :: arunoff_surf  !annual surface (layer 1) runoff (mm)

real(sp) :: aaet          !annual actual evapotranspiration (mm/year)
real(sp) :: apet          !annual potential evapotranspiration (mm/year)

!real(sp), dimension(npft) :: aresp  !annual maintenance + growth respiration (g m-2)

real(sp) :: cropfrac
real(sp) :: unusable

!instantaneous gridcell state variables

real(sp) :: gdd         !current-year growing degree days
real(sp) :: mtemp_max   !temperature of the warmest month (deg C)  

!local state variables

!annual pft state variables

real(sp), dimension(npft) :: estab_pft

logical,  dimension(npft) :: estab_lim   !whether PFT within bioclimatic limits for establishment
logical,  dimension(npft) :: survive     !whether PFT within bioclimatic limits for survival
logical,  dimension(npft) :: estab       !whether actually establishes (subject to land use)

real(sp), dimension(npft) :: turnover_ind  !total turnover of living biomass per individual (gC)
!real(sp), dimension(npft) :: heatstress    !reduction in individual density (& establishment) due to heat induced mortality  (indiv/m2)
real(sp), dimension(npft) :: pftCflux	   !total annual carbon flux from agricultural biomass burning, per PFT (g m-2 yr-1)

real(sp), dimension(npft,ncvar) :: bm_inc  !annual biomass increment (gC/m2)
real(sp), dimension(npft,ncvar) :: alresp  !annual gridcell leaf respiration (gC/m2)

!height class variables

type(sizeind), dimension(nhclass,npft) :: hclass

!daily pft state variables

real(sp), dimension(npft) :: wscal        !mean daily water scalar (among leaf-on days) (0-1 scale) 

real(sp), dimension(365,npft) :: dphen    !net daily leaf-on fraction
real(sp), dimension(365,npft) :: dphen_t  !daily leaf-on fraction due to temperature phenology
real(sp), dimension(365,npft) :: dphen_w  !daily leaf-on fraction due to drought phenology

real(sp), dimension(365,npft) :: wscal_v  !daily supply/demand ratio

!monthly gridcell state variables

real(sp), dimension(12) :: mrunoff     !total monthly runoff (mm)
real(sp), dimension(12) :: mtemp_soil  !monthly soil temperature (deg C)
real(sp), dimension(12) :: mw1         !monthly soil layer 1 water content (fraction of available water holding capacity)

real(sp), dimension(12,ncvar) :: mrh   !monthly heterotrophic respiration (gC/m2)

!monthly pft state variables

real(sp), dimension(12,npft,ncvar) :: mgpp   !monthly grid cell GPP (gC/m2)
real(sp), dimension(12,npft,ncvar) :: mlresp !monthly grid cell leaf respiration (gC/m2)
real(sp), dimension(12,npft,ncvar) :: mnpp   !monthly gridcell NPP (gC/m2)

real(sp), dimension(npft,12)	   :: mpftCflux		!monthly carbon flux from agricultural biomass burning, per pft (g C m-2 month-1)

real(sp), dimension(4) :: treecarbon

!----------------------------
!pointers to shared variables

!scalars

real(sp), pointer :: afire_frac
real(sp), pointer :: gdd20        !20-year average growing degree days
real(sp), pointer :: k_fast_ave   !running average k_fast for subroutine littersom
real(sp), pointer :: k_slow_ave   !running average k_slow for subroutine littersom
real(sp), pointer :: litterC_bg  
real(sp), pointer :: litterC_fast
real(sp), pointer :: litterC_slow
real(sp), pointer :: livebiomass
real(sp), pointer :: albiomass	  ! total carbon in aboveground living biomass (g m-2 yr-1) Chaste 03.03.2016
real(sp), pointer :: mtemp_min20  !20-year average minimum monthly temperature (deg C)
real(sp), pointer :: mat20        !20-year average mean annual temperature (deg C)
real(sp), pointer :: snow0
real(sp), pointer :: tilecarbon
real(sp), pointer :: forager_pd  !forager population density (persons km-2)
real(sp), pointer :: grasscover  !grass cover
real(sp), pointer :: dgrassdt    !grass annual change

!ntiles

real(sp), pointer, dimension(:) :: coverfrac

!daily meteorology gridcell state variables

real(sp), dimension(365) :: dmelt       !daily snowmelt (mm)
real(sp), dimension(365) :: snowpack    !daily snow water equivalent (mm)
real(sp), dimension(365) :: dtemp_soil  !daily soil temperature (deg C)
real(sp), dimension(365) :: dw1         !daily soil layer 1 water content

!climbuf

real(sp), pointer, dimension(:) :: mat_buf
real(sp), pointer, dimension(:) :: gdd_buf
real(sp), pointer, dimension(:) :: mtemp_min_buf

!soil layer

!real(sp), pointer, dimension(:) :: cond
!real(sp), pointer, dimension(:) :: whc
real(sp), pointer, dimension(:) :: w

!npft

logical,  pointer, dimension(:) :: present
logical,  pointer, dimension(:) :: leafon
integer,  pointer, dimension(:) :: leafondays
integer,  pointer, dimension(:) :: leafoffdays
real(sp), pointer, dimension(:) :: crownarea
real(sp), pointer, dimension(:) :: dwscal365
real(sp), pointer, dimension(:) :: fpc_grid
real(sp), pointer, dimension(:) :: fpc_inc
real(sp), pointer, dimension(:) :: height
real(sp), pointer, dimension(:) :: pftalbiomass
real(sp), pointer, dimension(:) :: lai_ind
real(sp), pointer, dimension(:) :: nind
real(sp), pointer, dimension(:) :: plant_carbon
real(sp), pointer, dimension(:) :: above_carbon

!ncvar

real(sp), pointer, dimension(:) :: acflux_fire
real(sp), pointer, dimension(:) :: arh
real(sp), pointer, dimension(:) :: cpool_surf   !surface SOM pool (2yr turnover time)
real(sp), pointer, dimension(:) :: cpool_fast
real(sp), pointer, dimension(:) :: cpool_slow
real(sp), pointer, dimension(:) :: grid_npp
real(sp), pointer, dimension(:) :: grid_gpp ! E.C. avril 2016
real(sp), pointer, dimension(:) :: litter_decom_ave
real(sp), pointer, dimension(:) :: acflux_estab

!npft,ncvar

real(sp), pointer, dimension(:,:) :: anpp
real(sp), pointer, dimension(:,:) :: agpp    !annual gridcell GPP (gC/m2)
real(sp), pointer, dimension(:,:) :: lm_ind
real(sp), pointer, dimension(:,:) :: hm_ind
real(sp), pointer, dimension(:,:) :: sm_ind
real(sp), pointer, dimension(:,:) :: rm_ind
real(sp), pointer, dimension(:,:) :: cstore
real(sp), pointer, dimension(:,:) :: litter_ag_fast
real(sp), pointer, dimension(:,:) :: litter_ag_slow
real(sp), pointer, dimension(:,:) :: litter_bg

!npft,month

real(sp), pointer, dimension(:,:) :: mBBpft
real(sp), pointer, dimension(:,:) :: mLAI

!month
real(sp), pointer, dimension(:) :: mburnedf    !monthly burned area fraction of gridcell
integer, pointer, dimension(:)  :: mnfire      !monthly number of fires
real(sp), pointer, dimension(:) :: mieff       !monthly average ignition efficieny
integer, dimension(12) :: nosnowdays           !number of days in a month with temp > 0. without snowcover

!additional local variables
real(sp) :: avg_cont_area	! average size of a natural patch at a given landuse fraction of the gridcell, in m2
!real(sp) :: nolanduse_frac	! 1. - coverfrac(2)
real(sp) :: totnat_area		! total natural area of a gridcell at a given landuse fraction, in m2
real(sp) :: avg_patch_number 	! average number of natural patches per gridcell at a given landuse fraction
real(sp) :: median_distance	! auxiliary for median distance to the edge of a natural patch, in (m)
real(sp) :: nbl			! normalized boundary length; for boundary between natural and used part; normalized to 
                                ! the max. possible boundary length when having a chessboard-type distribution of kernels
integer(sp) :: allnosnowdays                                  

real(sp) :: forager_ppd
real(sp) :: forager_fin
real(sp) :: forager_fout
real(sp) :: FDI
real(sp) :: omega_o0

real(sp), dimension(4) :: omega0

real(sp), dimension(npft)  :: BBpft  !biomass burned from each PFT, sum of live and dead consumed (g dry matter m-2)

real(sp) :: Ab	!area burned, ha per day
integer :: numfires_nat !number of natural fires per day
real (sp) :: ieff !ignition efficiency per day

real(sp), allocatable, dimension(:,:) :: help_me
integer,  allocatable, dimension(:,:) :: help_me2
real(sp), allocatable, dimension(:,:) :: help_me3

integer :: startyr_foragers

!declarations end here
!---------------------------------------------------------------

!write(*,*)

! write(stdout,*) 'working on gridcell ', in%year, in%lon, in%lat, in%human%foragerPD, in%slope !0,'(a,i5,3f14.2)'

ntiles = count(in%human%landuse >= 0.)

spinup = in%spinup
year   = in%year

! write(stdout,*)'start',in%lon,in%lat,year

co2 = in%co2

!co2 = 324.
 
if (ntiles > 1) then
  lucc = .true.
else
  lucc = .false.
end if

!import handle incoming input and state variables

met_out(1) = osv%met   !FLAG FLAG FLAG FLAG FLAG - correct?

tmin = in%climate%temp - 0.5 * in%climate%trng      ! MP: temp is average temperature, trng is temperature range
tmax = in%climate%temp + 0.5 * in%climate%trng
cldf = in%climate%cldp * 0.01
lght = in%climate%lght
wind = in%climate%wind
wetd = in%climate%wetd
prec = in%climate%prec

startyr_foragers = in%startyr_foragers

!write(0,*)'startyr_foragers',startyr_foragers

!if (in%idx == 1 .and. in%year <= 10) then
!   write(*,*)'year', year, osv%tile(1)%litter_ag_fast(:,1)
!  write(*,*)
!end if
 
!  if (in%idx == 1 .and. in%year >= 1) then
!   write(stdout,*) 'year: ', year,in%lon,in%lat
!   write(stdout,'(12f9.2)')in%climate%temp
!   write(stdout,'(12f9.2)')in%climate%temp0
!   write(stdout,'(12f9.2)')in%climate%prec
!   write(stdout,'(12f9.2)')in%climate%cldp
!   write(stdout,'(12f9.2)')in%climate%wetd
!   write(stdout,'(12f9.2)')in%climate%trng
!   write(stdout,'(12f9.2)')in%climate%wind
!   write(stdout,'(12f9.5)')in%climate%lght
!   write(stdout, *)
!   read(*,*)
!   ! stop
!  end if

!  write(stdout,'(a,12f9.2)')'WIND',in%climate%wind


! write(stdout,*)'initializing soil state'
! write(stdout,*)in%soil%sand
! write(stdout,*)in%soil%clay
! write(stdout,*)in%soil%orgm
! write(stdout,*)in%soil%zpos

if (year == 1) then  !initialize the soil state
  do i = 1,ntiles
    osv%tile(i)%soil%sand = in%soil%sand !* 100.
    osv%tile(i)%soil%clay = in%soil%clay !* 100.
    osv%tile(i)%soil%orgm = in%soil%orgm
    osv%tile(i)%soil%zpos = in%soil%zpos
  end do
end if

    !osv%tile(1)%coverfrac = 1. - in%human%landuse(i)  !should not be necessary for correct initialization

!--------------

!smooth interpolation to pseudo-daily values for temperature, cloud, lightning and windspeed

!call rmsmooth(tmin,ndaymonth,[tmin(12),tmin(1)],dtmn)
!call rmsmooth(tmax,ndaymonth,[tmax(12),tmax(1)],dtmx)
!call rmsmooth(cldf,ndaymonth,[cldf(12),cldf(1)],dcld)
!call rmsmooth(lght,ndaymonth,[lght(12),lght(1)],dlght)
!call rmsmooth(wind,ndaymonth,[wind(12),wind(1)],dwind)

call daily(tmin,dtmn,.true.)
call daily(tmax,dtmx,.true.)
call daily(cldf,dcld,.true.)
call daily(lght,dlght,.true.)
call daily(wind,dwind,.true.)

call daily(in%climate%prec,dprec,.false.)

dtemps = 0.5 * (dtmn + dtmx)
  
!write(stdout,*)'flag A'

!correct for potential out of bounds interpolation

dcld  = max(min(dcld,1.),0.)
dprec = max(dprec,0.)
met_out%lght = max(dlght,0.)
met_out%wind = max(dwind,0.)

!---------------------------------------------------
!initialize annual values for radiation calculations

Ratm = elev_corr(in%elev)

tcm = minval(in%climate%temp)  !temperature of the coldest month

!day loop for weather generator and radiation calculations

call weathergen_driver(dtmn,dtmx,dcld,prec,wetd,lght,met_out)

!we need to set a value for temp in order to calculate Pjj

temp = in%climate%temp

!write(stdout,*)'flag A1',temp

call calcPjj(temp,prec,Pjj)  !precipitation equitability index

!write(stdout,*)'flag A2'

do dyr = 1,365  !calculate radiation budget and PET

 call radpet(in%orbit,in%lat,tcm,Pjj,dyr,Ratm,met_out(dyr))
 !write(*,*)'after radpet', dyr,met_out(dyr)%prec,met_out(dyr)%lght
 
end do

!write(stdout,*)'flag A3'

!met_out%wind = max(dwind,0.)

!---------------------------------------------------
!write(stdout,*)'flag B'

!dry = 0
!maxdry = 0.
!maxprec = 0.

!do d = 1, 365
!  if(dprec(d) == 0.) then
!    dry = dry + 1
!  else if(dprec(d) /= 0.) then
!    dry = 0
!  end if
  
!  if(dry > maxdry) maxdry = dry
!  if(dprec(d) > maxprec) maxprec = dprec(d)  
!end do

!write(*,'(i5,2f10.4,i5,2f10.2)') year, minval(met_out%tmin), maxval(met_out%tmax), maxdry, maxprec, sum(dprec)      

osv%met = met_out(365)  !store the last day of the year for the next year's simulation

dtmn = met_out%tmin
dtmx = met_out%tmax
dpet = met_out%dpet

apet = sum(dpet)

!-----------------------------------------------------
!FLAG TEMPORARY SOLUTIONS

dtemp = 0.5 * (dtmn + dtmx)

call daily(in%climate%temp,dtemp,.true.)
!call daily(in%climate%temp,dsunp,.true.)

!-----------------------------------------------------

!make midmonth vectors
!since we now have daily variable climate, we take an average over 5 days around the middle of the month

do i = 1,12
  !a = midday(i)-2
  !b = midday(i)+2
  !mpet(i)     = 1./(b-a) * sum(met_out(a:b)%dpet)
  !mdayl(i)    = 1./(b-a) * sum(met_out(a:b)%dayl)
  !mpar_day(i) = 1./(b-a) * sum(met_out(a:b)%srad) * 0.5 * 1000. !divide by half to get to PAR (J m-2 d-1)

  mpet(i)     = met_out(midday(i))%dpet
  mdayl(i)    = met_out(midday(i))%dayl
  mpar_day(i) = met_out(midday(i))%srad * 0.5 * 1000. !divide by half to get to PAR (J m-2 d-1)

end do

!--------------------------------------------------------------------------------------

mat_buf       => osv%mat_buf
mat20         => osv%mat20
mtemp_min20   => osv%mtemp_min20
gdd20         => osv%gdd20
mtemp_min_buf => osv%mtemp_min_buf
gdd_buf       => osv%gdd_buf
coverfrac     => osv%tile%coverfrac

if (year == 1) then
  gdd_buf       = -9999.
  mtemp_min_buf = -9999.
end if

!END HANDLE OLD-STYLE STATEVARS ARRAYS
!--------------------------------------------------------------------------------------

!write(stdout,*) 'start lpjcore: ', present

!write(stdout,*)'flag C'

!summergreen phenology

dphen   = 0.
dphen_t = 0.

temp0 = in%climate%temp0
temp  = in%climate%temp
prec  = in%climate%prec

!write(stdout,'(a,12f8.2)')'temp0',temp0
!write(stdout,'(a,12f8.2)')'temp ',temp

call summerphenology(pftpar,temp,dtemp,gdd,dphen_t,pft%summergreen,pft%tree)  !FLAG: subroutine does not exist???

!if (year==225) then
!do i = 1,365
!  write(stdout,*)i,dtemp(i),dprec(i),dphen_t(i,3)
!end do
!end if

!--------
!20-year average climate variables

call climate20(temp,dtemps,gdd,mtemp_min_buf,gdd_buf,mtemp_min20,gdd20,mtemp_max,mat20,mat_buf)

!--------------------------------------------------------------------------------------
!set survival and establishment based on PFT bioclimatic limits

call bioclim(mtemp_min20,gdd,mtemp_max,survive,estab_lim)

!--------------------------------------------------------------------------------------
if (lucc) then

  unusable = in%human%landuse(1)
  cropfrac = in%human%landuse(2)
  
  call alcc(j,in,osv,cropfrac,unusable,coverfrac,recoverf)
  
  !write(stdout,*)'alcc',unusable,cropfrac,coverfrac

end if

!----------------------------------------------------------------------------------------------------------------------
!inner loop for grid-cell level tiles

!write(stdout,'(a,i5,5f10.4)')'go tile loop',ntiles,in%human%landuse(1:2),coverfrac

!write(stdout,*)'flag D'

do i = 1,3 !ntiles

  if (coverfrac(i) <= 0.) cycle

  !point pointers

  acflux_estab     => osv%tile(i)%acflux_estab
  acflux_fire      => osv%tile(i)%acflux_fire
  afire_frac       => osv%tile(i)%afire_frac
  anpp             => osv%tile(i)%anpp
  agpp             => osv%tile(i)%agpp
  arh              => osv%tile(i)%arh
  cpool_surf       => osv%tile(i)%cpool_surf
  cpool_fast       => osv%tile(i)%cpool_fast
  cpool_slow       => osv%tile(i)%cpool_slow
  crownarea        => osv%tile(i)%crownarea
  dwscal365        => osv%tile(i)%dwscal365
  fpc_grid         => osv%tile(i)%fpc_grid
  fpc_inc          => osv%tile(i)%fpc_inc
  grid_npp         => osv%tile(i)%grid_npp
  grid_gpp         => osv%tile(i)%grid_gpp ! E.C. avril 2016
  height           => osv%tile(i)%height
  pftalbiomass     => osv%tile(i)%pftalbiomass  ! E.C. octobre 2016
  hm_ind           => osv%tile(i)%hm_ind
  k_fast_ave       => osv%tile(i)%k_fast_ave
  k_slow_ave       => osv%tile(i)%k_slow_ave
  lai_ind          => osv%tile(i)%lai_ind
  leafoffdays      => osv%tile(i)%leafoffdays
  leafon           => osv%tile(i)%leafon
  leafondays       => osv%tile(i)%leafondays
  litterC_bg       => osv%tile(i)%litterC_bg
  litterC_fast     => osv%tile(i)%litterC_fast
  litterC_slow     => osv%tile(i)%litterC_slow
  litter_ag_fast   => osv%tile(i)%litter_ag_fast
  litter_ag_slow   => osv%tile(i)%litter_ag_slow
  litter_bg        => osv%tile(i)%litter_bg
  litter_decom_ave => osv%tile(i)%litter_decom_ave
  livebiomass      => osv%tile(i)%livebiomass
  albiomass        => osv%tile(i)%albiomass
  lm_ind           => osv%tile(i)%lm_ind
  cstore           => osv%tile(i)%cstore
  nind             => osv%tile(i)%nind
  plant_carbon     => osv%tile(i)%plant_carbon
  above_carbon     => osv%tile(i)%above_carbon
  present          => osv%tile(i)%present
  rm_ind           => osv%tile(i)%rm_ind
  sm_ind           => osv%tile(i)%sm_ind
  snow0            => osv%tile(i)%snow0
  tilecarbon       => osv%tile(i)%tilecarbon
  w                => osv%tile(i)%w
  forager_pd       => osv%tile(i)%forager_pd
  grasscover       => osv%tile(i)%grasscover
  dgrassdt         => osv%tile(i)%dgrassdt
  mLAI             => osv%tile(i)%mLAI
  mBBpft           => osv%tile(i)%mBBpft
  mburnedf         => osv%tile(i)%mburnedf
  mnfire           => osv%tile(i)%mnfire
  mieff            => osv%tile(i)%mieff

  !--------------------------------------------------------------------------------------
  !initializations (needed?)
    
  pftCflux = 0.
  mpftCflux = 0.
  mBBpft = 0. 
  mburnedf = 0.
  mnfire = 0.
  mieff = 0.
  
  !set up soil parameters
  
  call simplesoil(osv%tile(i)%soil,soilpar)
  
  !write(stdout,*) 'after simplesoil'
  
  !write(stdout,*)'soilpar',soilpar

  !write(stdout,*)osv%tile(i)%soil
  
 !soilpar(1) =  19.060502
 !soilpar(2) =   5.789572
 !soilpar(3) =  30.
 !soilpar(4) = 280.
 !soilpar(5) = 0.2
 !soilpar(6) = 0.65
 !soilpar(7) = 0.4
  
  !write(stdout,*)soilpar(1),soilpar(3)
  !write(stdout,*)soilpar(2),soilpar(4)
  
  !--------------------------------------------------------------------------------------
  !control establishment of woody vegetation based on alcc and total density for recovering vegetation

  call tile_landuse(lutype(i),recoverf,estab_lim,estab,nind,fpc_grid)
    
  do a = 1, npft
   if ((nind(a) <= 0.) .and. ((lm_ind(a,1) .ne. 0.) .or. (sm_ind(a,1) .ne. 0.) .or. (hm_ind(a,1) .ne. 0.) .or. (rm_ind(a,1) .ne. 0.))) then
     write(stdout,'(a,i4,7f14.7)') 'contradiction nind - C-pools, lu: ',a, in%lon, in%lat, nind(a), lm_ind(a,1), sm_ind(a,1), hm_ind(a,1), rm_ind(a,1)
   end if
  end do

  !Establishment of new individuals (saplings) of woody PFTs, grass establishment, 
  !removal of PFTs not adapted to current climate, update of individual structure and FPC

  !if (i == 2) then
  !  write(stdout,*)'ag litter -1',litter_ag_fast(8,1),litter_ag_slow(8,1)
  !end if

! ====== NB special conditions for Canada version ONLY =======
!
! set up limits to etablishment for special conditions
! clay_mean = (osv%tile(i)%soil%clay(1) + osv%tile(i)%soil%clay(3))/2
! 
! if (clay_mean >= 20.) estab(1) = .false.
! if (clay_mean >= 13.) estab(3) = .false.
! if (clay_mean >= 18.) estab(4) = .false.
! if (clay_mean >= 23.) estab(8) = .false.
! 
! if (afire_frac == 0.) estab(4) = .false.

! ====== NB special conditions for Canada version ONLY =======

! if (year > 1) then
! 
!   burnedf20 = sum(osv%tile(i)%burnedf_buf) / real(climbuf)
! 
! 		if (burnedf20 > 0.) then
! 				FRI20 = 1. / burnedf20
! 		else
! 				FRI20 = 100000.  !just choose some arbitrary big number
! 		end if
! else
!   FRI20 = 100000.
! end if

! if (FRI20 < 50.) estab(1) = .false.
! if (FRI20 < 30.) estab(4) = .false.

! ===== end special conditions for Canada version ONLY =====
  
  call establishment(pftpar,present,survive,estab,nind,lm_ind,sm_ind,rm_ind,hm_ind,lm_sapl,sm_sapl,rm_sapl,hm_sapl,pft%tree, &
                     crownarea,fpc_grid,lai_ind,height,sla,wooddens,latosa,prec,reinickerp,litter_ag_fast,litter_ag_slow,litter_bg,  &
                     allom1,allom2,allom3,acflux_estab,leafondays,leafoffdays,leafon,estab_pft,afire_frac,osv%tile(i)%soil%clay, osv%tile(i)%soil%sand)
                     
!  if(i==2) write(stdout,'(a,i3,9f14.4)') 'after establishment',i, litter_ag_fast(:,1)               
                     
  do a = 1, npft
   if ((nind(a) <= 0.) .and. ((lm_ind(a,1) .ne. 0.) .or. (sm_ind(a,1) .ne. 0.) .or. (hm_ind(a,1) .ne. 0.) .or. (rm_ind(a,1) .ne. 0.))) then
     write(stdout,'(a,i4,7f14.7)') 'contradiction nind - C-pools, estab: ', a, in%lon, in%lat, nind(a), lm_ind(a,1), sm_ind(a,1), hm_ind(a,1), rm_ind(a,1)
   end if
  end do                     
                     
                     
  !if (i == 2) then
  !  write(stdout,*)'ag litter 0',litter_ag_fast(8,1),litter_ag_slow(8,1)
  !end if

!  write(stdout,'(a,9f14.7)') 'after establishment', nind

                     
  where (lm_ind(:,1) <= 0.) present = .false. 
    
!  write(stdout,*) 'reset after establishment', present    

  !light competition among trees and between trees and grasses

  call light(present,pft%tree,lm_ind,sm_ind,hm_ind,rm_ind,crownarea,fpc_grid,fpc_inc,nind,litter_ag_fast,litter_ag_slow,litter_bg,sla,fpc_tree_max) 
  
!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after light1',i, litter_ag_fast(:,1) 

!  write(stdout,'(a,10f8.3)')'after light',fpc_grid,sum(fpc_grid)

  !Adjust daily precipitation by snowmelt and accumulation in snowpack

  treefrac = sum(fpc_grid,mask = pft%tree)

  !if (in%idx == 1) write(*,*)'incoming',snow0,treefrac

  call snow(dtemp,dprec,treefrac,snow0,snowpack,dmelt,in%idx)  !new f90 routine

  snow0 = snowpack(365)
  
  nosnowdays = 0.

  a = 1
  do m = 1,12

    b = a + ndaymonth(m) - 1
    
    if(temp(m) > 0.) nosnowdays(m) = count(snowpack(a:b) <= 0.)

    a = b + 1

  end do
  
  allnosnowdays = sum(nosnowdays)

  !if (in%idx == 1) write(*,*)'outgoing',snow0,treefrac
  
!  write(stdout,*) 'after snow'

  !--------------------------------------------------------------------------------------
  !growth, respiration and vegetation dynamics

  !gross primary productivity

  !write(stdout,'(12f10.2)')temp
  !write(stdout,'(12f10.2)')in%climate%prec
  !write(stdout,'(12f10.3)')mpet
  !write(stdout,'(12f10.3)')mdayl
  !write(stdout,'(12f10.1)')mpar_day*0.001
  !write(stdout,'(12f10.1)')dcld(midday)*100.
  
  !do a = 1,npft
  !  write(stdout,*)a,fpc_grid(a),lai_ind(a)
!end do

  call calcgpp(present,[co2,-8.,0.],soilpar,pftpar,lai_ind,fpc_grid,mdayl,temp,mpar_day,dphen_t,w,dpet,dprec,dmelt,   &
               sla,agpp,alresp,arunoff_surf,arunoff_drain,arunoff,mrunoff,dwscal365,dphen_w,dphen,wscal,mgpp,mlresp,  &
               mw1,dw1,aaet,leafondays,leafoffdays,leafon,pft%tree,pft%raingreen,year,mat20,wscal_v,in%idx)           
  
!  write(stdout,*)'after calcgpp ',present

  !write(stdout,'(13f10.1)')agpp(8,1),mgpp(:,8,1)

!  write(stdout,*)'flag D2',lm_ind(7:8,1),present

  !calculate mid-month soil temperatures

  call soiltemp(soilpar,temp,temp0,mtemp_soil,mw1,in%idx)

  !interpolate monthly soil temperature to daily values
  
!  write(stdout,*)mtemp_soil
  
  !call rmsmooth(mtemp_soil,ndaymonth,[mtemp_soil(12),mtemp_soil(1)],dtemp_soil)
  call daily(mtemp_soil,dtemp_soil,.true.)

  !autotrophic respiration and NPP
  
!  write(stdout,*)'flag D2a',lm_ind(4,1) !dtemp(1),dtemp_soil(1)

  call calcnpp(dtemp,dtemp_soil,dphen,present,nind,lm_ind(:,1),sm_ind(:,1),hm_ind(:,1),rm_ind(:,1),cstore(:,1),mgpp(:,:,1),mnpp(:,:,1),anpp(:,1))
    
  bm_inc = anpp
  anpp = max(0.,anpp)
  
!  write(stdout,*)'flag D2b',lm_ind(7:8,1)

  !call npp(pftpar,dtemp,dtemp_soil,pft%tree,dphen,nind,year,lm_ind,sm_ind,rm_ind,mgpp,anpp,mnpp,bm_inc,present,agpp,co2,aresp)

  !write(stdout,'(i5,14f10.1)')year,anpp(8,1),bm_inc(8,1),mnpp(:,8,1)
    
  !write(stdout,'(a,9f11.4)')'GPP    ',agpp(:,1)
  !write(stdout,'(a,9f11.4)')'NPP    ',anpp(:,1)
  !write(stdout,'(a,9f11.4)')'respire',aresp(:)

!  write(stdout,*)'flag D2b1'

  !allocation to reproduction

  !if (i == 2) then
  !  write(stdout,*)'ag litter 0.5',litter_ag_fast(8,1),litter_ag_slow(8,1)
  !end if

  call reproduction(bm_inc,lm_sapl,sm_sapl,hm_sapl,rm_sapl,litter_ag_fast,litter_ag_slow,present,pft%tree,[co2,-8.,0.])
 
!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after reproduction',i, litter_ag_fast(:,1) 
!  write(stdout,*)'flag D2d',lm_ind(:,1)

  !leaf, sapwood, and fine-root turnover

  call turnover(pftpar,present,lm_ind,sm_ind,hm_ind,rm_ind,litter_ag_fast,litter_ag_slow,litter_bg,nind,turnover_ind,year)
  
!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after turnover',i, litter_ag_fast(:,1) 

  !if (i == 2) then
  !  write(stdout,*)'ag litter 2',litter_ag_fast(8,1),litter_ag_slow(8,1)
  !end if

!  write(stdout,*)'flag D2a',lm_ind(7:8,1),present(7:8),bm_inc(7:8,1)

  !-------------------------------------------------------------------------
  !removal of PFTs with negative C increment this year

  !write(stdout,*)

  call killplant(bm_inc,present,pft%tree,lm_ind,rm_ind,hm_ind,sm_ind,nind,litter_ag_fast,litter_ag_slow,litter_bg)
  
!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after killplant',i, litter_ag_fast(:,1) 

  !allocation of annual carbon increment to leaf, stem and fine root compartments
  
!  write(stdout,*)'flag D3a',present
  !write(stdout,*)'flag D3b',lm_ind(:,1)
  !write(stdout,*)'flag D3a',hm_ind(:,1)
  !write(stdout,*)'flag D2b',lm_ind(:,1)
  
  call allocation(pftpar,allom1,allom2,allom3,latosa,wooddens,reinickerp,pft%tree,sla,wscal,nind,bm_inc,lm_ind,sm_ind,     &
                  hm_ind,rm_ind,crownarea,fpc_grid,lai_ind,height,litter_ag_fast,litter_ag_slow,litter_bg,fpc_inc,present)
                  
!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after allocation',i, litter_ag_fast(:,1)                               
                  
  !check validity of allocation and correct
  !heartwood can be zero, but all other pools have to be positive to have valid allometry
  
  do a = 1,npft
    if (.not.pft(a)%tree) cycle

    treecarbon(1) = lm_ind(a,1)
    treecarbon(2) = sm_ind(a,1)
    treecarbon(3) = rm_ind(a,1)
    treecarbon(4) = hm_ind(a,1)
    
    if (any(treecarbon(1:3) <= 0.) .and. (sum(treecarbon) > 0. .or. nind(a) > 0.)) then
      
      write(stdout,*)'invalid allometry, resetting',year, present(a)
      write(stdout,'(2f10.2,i4,5f14.7)')in%lon,in%lat,a,nind(a),lm_ind(a,1),sm_ind(a,1),hm_ind(a,1),rm_ind(a,1)
      
      !read(*,*) 
      
      !stop

      litter_ag_fast(a,1) = litter_ag_fast(a,1) + nind(a) * lm_ind(a,1)
      litter_ag_slow(a,1) = litter_ag_slow(a,1) + nind(a) *(sm_ind(a,1) + hm_ind(a,1))
      litter_bg(a,1)      = litter_bg(a,1)      + nind(a) * rm_ind(a,1)

      present(a)  = .false.
      nind(a)     = 0.
      lm_ind(a,1) = 0.
      rm_ind(a,1) = 0.
      sm_ind(a,1) = 0.
      hm_ind(a,1) = 0.

    end if
  end do
  
  !Implement light competition and background mortality among tree PFTs 
  !(including heat damage and due to lower limit of npp for boreal trees) 

  !write(stdout,*)'after alloc'
  !do a = 1,npft
  !  write(stdout,*)a,fpc_grid(a),lai_ind(a)
  !end do
  !write(stdout,*)
  
  !if (any(lm_ind(:,1) < 0.)) then
!    write(stdout,*)'flag D3a',in%lon,in%lat,lm_ind(:,1)
  !end if
  
  call mortality(pftpar,present,pft%tree,pft%boreal,bm_inc,turnover_ind,sla,lm_ind,sm_ind,hm_ind,rm_ind,nind,year,litter_ag_fast,litter_ag_slow, &
                 litter_bg,dtemp,anpp,mtemp_max) 
                 
  do a = 1, npft
   if ((nind(a) <= 0.) .and. ((lm_ind(a,1) .ne. 0.) .or. (sm_ind(a,1) .ne. 0.) .or. (hm_ind(a,1) .ne. 0.) .or. (rm_ind(a,1) .ne. 0.))) then
     write(stdout,'(a,i4,7f14.7)') 'contradiction nind - C-pools, mortality: ',a, in%lon, in%lat, nind(a), lm_ind(a,1), sm_ind(a,1), hm_ind(a,1), rm_ind(a,1)
   end if  
  end do                 
                 
!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after mortality',i, litter_ag_fast(:,1)                           


  !light competition between trees and grasses

  call light(present,pft%tree,lm_ind,sm_ind,hm_ind,rm_ind,crownarea,fpc_grid,fpc_inc,nind,litter_ag_fast,litter_ag_slow,litter_bg,sla,fpc_tree_max)
  
!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after light2',i, litter_ag_fast(:,1) 

  !write(stdout,'(a,i4,10f8.3)')'before fire',year,fpc_grid,sum(fpc_grid)

!  write(stdout,*)'flag D3b'

  !----------------------------------------------------------------------------
  !disaggregate average individual to allometrically consistent height/age classes

  do a = 1,npft
    if (present(a) .and. pft(a)%tree) then
     
    
     call allomind(a,height(a),hclass(:,a))

     !write(stderr,*)'call height',a,height(a)   !,hclass(:,a)


    end if
  end do

  !----------------------------------------------------------------------------
  !section: harvest and burning on managed land - wildfire on unmanaged land

  osv%carbon%crop_harvest = 0.
  
  !----------------------------------------------------------------------------
  !landscape fractioning, based on assumption that subgrid-kernels will be randomely distributed; calculations based on a testgrid of 10000 sub-kernels

  call landscape_fractality(coverfrac, in%cellarea, avg_cont_area, totnat_area, avg_patch_number, median_distance, nbl) 
  
  !---------------------------------------------------------------------------- 
  
!  goto 20
  
  if (lutype(i) == 2) then   !harvest of agricultural crops and reset of agricultural biomass on used land
    
!    write(stdout,'(a,i3,9f14.4)') 'before harvest',i, litter_ag_fast(:,1)    

    call harvest(i,j,osv)
    
!    write(stdout,'(a,i3,9f14.4)') 'after harvest',i, litter_ag_fast(:,1)
        
    call managedburn(i,j,acflux_fire(1),afire_frac,litter_ag_fast(:,1),pftCflux)
    
!    write(stdout,'(a,13i5,f14.7)') 'after managedburn', nosnowdays, allnosnowdays, afire_frac
    
    do m = 1, 12
    
         mBBpft(:,m) = pftCflux(:) * real(nosnowdays(m))/real(allnosnowdays)  * 0.001 * 1./0.45	!distribute the pftCflux equally on all days with a monthly temperature > 0. degrees Celsius
      									!convert from g C to kg dry matter, hence the multiplication factors
         mburnedf(m) = afire_frac * real(nosnowdays(m))/real(allnosnowdays)       
                    
    end do     
      
!    write(stdout,'(9f14.7)')  mpftCflux
!    write(stdout,'(12f14.7)')  temp

  else     !biomass destruction by wildfire

    afire_frac = 0.
    
!    goto 20

    if (dospitfire .and. ((spinup .and. year > 0) .or. .not. spinup)) then
      
      burnedf20 = sum(osv%tile(i)%burnedf_buf) / real(climbuf)
      
      forager_pd20 = sum(osv%tile(i)%forager_pd_buf) / real(climbuf) 
    
!      write(stdout,*) 'Forager density: ', year, forager_pd20

      !calculate annual burn target

!      write(stdout,*)osv%annburntarget(1),grasscover,dgrassdt
      
      if (year > 800) then
!        if (grasscover < 0.5 .and. abs(dgrassdt) > 0.01) then
!          osv%annburntarget(1) = min(osv%annburntarget(1) + 0.05,1.)
!        else
!          osv%annburntarget(1) = max(osv%annburntarget(1) - 0.05,0.)
!       end if
!        osv%annburntarget(1)= max(min(1.*(1-grasscover)*max(-dgrassdt,0.)/0.01,1.),0.)
        osv%annburntarget(1)= max(min(1.*(1-grasscover)*max(-dgrassdt,0.)/0.05,1.),0.) 
      end if

!      write(stdout,*)osv%annburntarget(1),grasscover,dgrassdt

      omega_o0 = 0.    !initialize omega_o0 here, this should carry over from one year to the next but ok for Boreal simulations
      d = 1
      do m = 1,12
        
        mburnedf(m) = 0.
        mBBpft(:,m) = 0.
        mnfire(m)   = 0.
        mieff(m)    = 0.
        
        do dm = 1,ndaymonth(m)

          call spitfire(year,i,j,d,in,met_out(d),dw1(d),snowpack(d),dphen(d,:),wscal_v(d,:),osv,spinup,avg_cont_area,burnedf20,forager_pd20,FDI,omega_o0,omega0,BBpft,Ab,hclass,numfires_nat,ieff)
          mBBpft(:,m) = mBBpft(:,m) + BBpft  !accumulate biomass burned totals
          
          mburnedf(m) = mburnedf(m) + Ab/(in%cellarea * 1e-4) !convert cell area to ha, as Ab is in ha

          if (numfires_nat <0) then
           numfires_nat = 0
          end if  
          mnfire(m)   = mnfire(m) + numfires_nat

          if (ieff <0) then
           ieff = 0
          end if
          mieff(m)    = mieff(m) + (ieff/ndaymonth(m))

          d = d + 1

        end do
      end do
      
      !add the agricultural burned biomass to the total burned biomass, by PFT

      !write(stdout,*)'enter biomassburned',afire_frac
      if (afire_frac > 0.) call burnedbiomass(i,j,osv,afire_frac,year)  !account for annual burned biomass
        
      !fill burnedf buffer with burned fraction for this year
      
      osv%tile(i)%burnedf_buf = eoshift(osv%tile(i)%burnedf_buf,-1,afire_frac)

      !if (i==1) write(*,'(f7.4)')burnedf20 !osv%tile(i)%burnedf_buf

    else         !option to use old LPJ fire routine

      call fire(pftpar,dtemp,litter_ag_fast,litter_ag_slow,acflux_fire,afire_frac,lm_ind,rm_ind,sm_ind,hm_ind,nind,dw1,present,pft%tree,year)

    end if
  end if
  
!  write(stdout,*) 'tile', i
!  write(stdout,'(12f12.7)') mburnedf
!  write(stdout,*)
  
  20 continue

!  if (i==1) write(*,'(a,2i4,6f14.7)') 'ANNUAL BURNEDF: ', year, i, afire_frac, in%human%popd
!  write(*,*)

  do a = 1, npft
   if ((nind(a) <= 0.) .and. ((lm_ind(a,1) .ne. 0.) .or. (sm_ind(a,1) .ne. 0.) .or. (hm_ind(a,1) .ne. 0.) .or. (rm_ind(a,1) .ne. 0.))) then
     write(stdout,'(a,i4,7f14.7)') 'contradiction nind - C-pools, fire: ',a, in%lon, in%lat, nind(a), lm_ind(a,1), sm_ind(a,1), hm_ind(a,1), rm_ind(a,1)
   end if
  end do

!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after fire',i, litter_ag_fast(:,1)   
!  write(stdout,'(a,i4,10f8.3)') 'after fire: ', year, fpc_grid,sum(fpc_grid)

  !----------------------------------------------------------------------------
  !section: litter and soil decomposition calculations

  !This is done after fire, so that fire probability is calculated on litter remaining
  !before year's decomposition. Include here agricultural litter etc.

  call hetresp(litter_ag_fast,litter_ag_slow,litter_bg,mw1,mtemp_soil,cpool_surf,cpool_fast,cpool_slow,arh,mrh,year, & 
                  k_fast_ave,k_slow_ave,litter_decom_ave,osv%tile(i)%soil%clay,osv%tile(i)%soil%bulk,spinup,in%idx)
                  
!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after hetresp',i, litter_ag_fast(:,1)                   

  !call littersom2(litter_ag_fast,litter_ag_slow,litter_bg,mw1,mtemp_soil,cpool_fast,cpool_slow,arh,mrh,year, & 
  !                k_fast_ave,k_slow_ave,litter_decom_ave,osv%tile(i)%soil%clay,osv%tile(i)%soil%bulk,spinup)


  !light competition between trees and grasses

  if (any(fpc_grid < 0.)) then
    write(stdout,*)'light1',in%lon,in%lat,fpc_grid
    stop
  end if

  call light(present,pft%tree,lm_ind,sm_ind,hm_ind,rm_ind,crownarea,fpc_grid,fpc_inc,nind,litter_ag_fast,litter_ag_slow,litter_bg,sla,fpc_tree_max)
  
!  if(i==2)   write(stdout,'(a,i3,9f14.4)') 'after light3',i, litter_ag_fast(:,1) 
  
  if (any(fpc_grid < 0.)) then
    write(stdout,*)'light2',in%lon,in%lat,fpc_grid  !'(a,2f10.2,9f14.7)'
    stop
  end if
  
!  write(stdout,'(a,i4,10f8.3)')'after light',year,fpc_grid,sum(fpc_grid)
!  write(stdout,*) 




  !radioactive decay of 14C

  !call decay(present,lm_ind,sm_ind,hm_ind,rm_ind,litter_ag_fast,litter_ag_slow,litter_bg,cpool_fast,cpool_slow)

  !calculate soil erosion

  !call rusle(i,j)

  !update soils for addition or loss of organic matter

  !call updatesoil(i,j)

!  write(stdout,*)'flag D4'

  !---------------
  !summary section
  
  !plant_carbon = nind * (lm_ind(:,1) + sm_ind(:,1) + hm_ind(:,1))  !only aboveground part for now
  plant_carbon = nind * (lm_ind(:,1) + sm_ind(:,1) + hm_ind(:,1) + rm_ind(:,1))  !pft vector
  above_carbon = nind * (lm_ind(:,1) + sm_ind(:,1) + hm_ind(:,1))  !pft vector
 
!  plant_carbon = nind * (lm_ind(:,1) + sm_ind(:,1) + hm_ind(:,1))
  
  !sum across pfts

  livebiomass  = sum(plant_carbon)
  albiomass    = sum(above_carbon)
  litterC_fast = sum(litter_ag_fast(:,1))
  litterC_slow = sum(litter_ag_slow(:,1))
  litterC_bg   = sum(litter_bg(:,1))
  grid_npp(1)  = sum(anpp(:,1))
  grid_gpp(1)  = sum(agpp(:,1)) ! E.C. avril 2016

  tilecarbon = livebiomass + litterC_fast + litterC_slow + litterC_bg + cpool_surf(1) + cpool_fast(1) + cpool_slow(1)
  
!  write(stdout,*) 'before forager routines'

!  goto 40

  if (i /= 2 .and. ((spinup .and. year > 100) .or. .not. spinup) .and. in%human%foragerPD > 0.) then 
    
!    call foragers(apet,aaet,in%elev,in%lat,grid_npp(1),livebiomass,soilpar(3),temp,prec,mw1,forager_ppd)
    
!    write(stdout,*) 'after call to foragers: ', forager_ppd
    
!    call popgrowth(forager_ppd,forager_pd,forager_fin,forager_fout)
    
!    write(stdout,*) 'after call to popgrowth:', forager_pd

    !write(stdout,*)'before call to simpleforagers',year,forager_pd

    call simpleforagers(in%human%foragerPD,pft%tree,anpp(:,1),fpc_grid,forager_pd)

    !write(*,*)year,forager_pd

    !forager_pd = 0.025
    
    osv%tile(i)%forager_pd_buf = eoshift(osv%tile(i)%forager_pd_buf,-1,forager_pd)
      
  else if(in%human%foragerPD == 0.) then
    
    forager_pd = 0. 
    
    osv%tile(i)%forager_pd_buf = 0.   

  end if
  
!  40 continue

  a = 1
  do m = 1,12

    b = a + ndaymonth(m) - 1
    
    mLAI(:,m) = lai_ind * sum(dphen(a:b,:),dim=1) / ndaymonth(m)

    a = b + 1

  end do

  !---------------------------------------------------------------------------- 
  !10-year running mean grass cover and long-term rate of change in grass cover

  dgrassdt   = grasscover - (0.9 * grasscover + 0.1 * sum(fpc_grid(8:9)))
  grasscover = 0.9 * grasscover + 0.1 * sum(fpc_grid(8:9))

  !---------------------------------------------------------------------------- 
    
end do  !sub-grid tile loop

!write(stdout,*)'flag E'

!write(*,*)osv
!write(*,*)
!write(stdout,*)'end',year,in%lon,in%lat,osv%tile(1)%fpc_grid(8) !sum(met_out%NI)

99 format(a,i5,2f8.2,f10.6)

!write(stdout,*)'outgoing',osv%tile(1)%w

!if (in%idx == 1) then
!  write(*,*)'outgoing W',osv%tile(1)%w
!end if

allocate(help_me(ntiles,12))
allocate(help_me2(ntiles,12))
allocate(help_me3(ntiles,12))

do i = 1, ntiles
!  write(stdout,'(i8,13f14.7)') i, osv%tile(i)%mburnedf, osv%tile(i)%coverfrac
  help_me(i,:) = osv%tile(i)%mburnedf * osv%tile(i)%coverfrac
  help_me2(i,:) = osv%tile(i)%mnfire
  help_me3(i,:) = osv%tile(i)%mieff
end do

!to make it easier in netcdfoutputmod, and since we are only interested in the total, not tilewise, already do the tile-integration here and then put it 
!back in osv%tile(1)%mburnedf, and then put that out in netcdfoutputmod

osv%tile(1)%mburnedf = sum(help_me(:,:), dim=1)
osv%tile(1)%mnfire   = sum(help_me2(:,:), dim=1)
osv%tile(1)%mieff    = sum(help_me3(:,:), dim=1)

!write(stdout,'(a,12f14.7)') 'integrated', osv%tile(1)%mburnedf

deallocate(help_me)
deallocate(help_me2)
deallocate(help_me3)

end subroutine lpjcore

!------------------------------------------------------------------------------------------------------------

end module lpjmod
