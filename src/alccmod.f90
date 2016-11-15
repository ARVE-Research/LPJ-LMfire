module alccmod

implicit none

public :: alcc
public :: tile_landuse
public :: harvest

contains

!------------------------------------------------------------------------------------------------------------

subroutine alcc(j,input,osv,cropfrac,unusable,coverfrac,recoverf)

use parametersmod,   only : sp,npft
use mpistatevarsmod, only : inputdata,statevars

implicit none

!parameters

real(sp), dimension(2), parameter :: wd    = [ exp(-1./2.), exp(-1./20.) ]  !inverse of the wood pool turnover times (yrs-1)
integer,  dimension(2), parameter :: lu_id = [ 1,3 ]

!arguments

integer,  intent(in)    :: j
type(inputdata), intent(in) :: input
real(sp), intent(out)   :: recoverf
real(sp), intent(inout) :: unusable   !fraction of the gridcell that cannot be used
real(sp), intent(inout) :: cropfrac
type(statevars), target, intent(inout) :: osv

real(sp), dimension(:), intent(inout) :: coverfrac

!pointers

logical,  pointer, dimension(:) :: present
real(sp), pointer, dimension(:) :: lm_ind
real(sp), pointer, dimension(:) :: sm_ind
real(sp), pointer, dimension(:) :: hm_ind
real(sp), pointer, dimension(:) :: rm_ind

real(sp), pointer, dimension(:) :: litter_ag_fast
real(sp), pointer, dimension(:) :: litter_ag_slow
real(sp), pointer, dimension(:) :: litter_bg

real(sp) :: tf  !fraction of the used land that is cycled through every year

!local variables

integer  :: i
integer  :: l

real(sp) :: lndturn
real(sp) :: usable

real(sp) :: convf
real(sp) :: excess

real(sp), dimension(2)    :: convft
real(sp), dimension(npft) :: nindcf

!--------------------------------------

!write(0,*)'flag 1'

tf = input%human%lu_turnover

convf  = 0.
convft = convf

cropfrac = min(max(0.,cropfrac),1.) !make sure cropfrac is between 0 and 1

!calculate the change in anthropogenic land use fraction
!positive means increased land use, negative means abandonment

!write(0,*)'flag 1a',tf

lndturn = tf * cropfrac
usable  = 1. - unusable

if (cropfrac > usable) then   !this should only happen during the interpolation between 1850 and present-day
  usable   = cropfrac
  unusable = 1. - usable
end if

!write(0,*)'flag 2'

!if the deforestation + turnover is greater than the usable amount, limit clearance to the usable fraction

convf = min(cropfrac - coverfrac(2) + lndturn,usable)  

coverfrac(2) = cropfrac

if (convf > 0.) then        !conversion to ag land use, there will always be some, because we have turnover of ag land

  !first take from the natural vegetation
  
  convft(1) = min(convf,coverfrac(1) - unusable)

  coverfrac(1) = coverfrac(1) - convft(1)  !includes the turnover fraction
    
  !if all natural vegetation is gone, take from the recovering fraction
  
  if (convf > convft(1)) then
    
    excess = convf - convft(1)

    coverfrac(3) = coverfrac(3) - excess
    convft(2)    = excess

    coverfrac(1) = unusable

  end if

end if

!recovery occurs all the time when there is land turnover, even if convf is positive (deforestation)

recoverf     = lndturn - min(convf,0.)         !convf will be negative under afforesting conditions
coverfrac(3) = coverfrac(3) + recoverf 
coverfrac(3) = max(coverfrac(3),0.)            !reset if below zero - this shouldn't happen, but there seems to be rounding error

!------------------
!transfer carbon stock of converted fraction to the C pools of the anthro land use fraction

litter_ag_fast => osv%tile(2)%litter_ag_fast(:,1)
litter_ag_slow => osv%tile(2)%litter_ag_slow(:,1)
litter_bg      => osv%tile(2)%litter_bg(:,1)

!if (convf /= 0.) then
!  write(0,*)'convf',convf,convft,coverfrac,unusable
!  read(*,*)
!end if

do i = 1,2  !once for each case - conversion from natural veg or recovering veg

  if (convft(i) > 0.) then

    l = lu_id(i)

    !reduce the individual density of living biomass by the fraction converted to human land use

    nindcf = osv%tile(l)%nind * convft(i)

    !pft vectors

    lm_ind => osv%tile(l)%lm_ind(:,1)
    sm_ind => osv%tile(l)%sm_ind(:,1)
    hm_ind => osv%tile(l)%hm_ind(:,1)
    rm_ind => osv%tile(l)%rm_ind(:,1)
    
    present => osv%tile(l)%present
    
    !carbon accounting of conversion fractions following Strassman et al., 2008 more or less

    !leaf mass: 100% to litter-ag-fast
       
    litter_ag_fast = litter_ag_fast + nindcf * lm_ind

    !root mass: 100% to litter-bg

    litter_bg = litter_bg + nindcf * rm_ind

    !stem mass: 25% burned immediately, 75% to two product pools (50% each)

    osv%tile(l)%acflux_conv(1) = 0.25 * sum(nindcf * (sm_ind + hm_ind),mask=present)
    
    osv%carbon%wood_fast = osv%carbon%wood_fast + 0.375 * sum(nindcf * (sm_ind + hm_ind), mask=present)
    osv%carbon%wood_slow = osv%carbon%wood_slow + 0.375 * sum(nindcf * (sm_ind + hm_ind), mask=present)
    
    !soil organic matter is the weighted average of whatever is currently on the used tile and the state of the land being converted
    
    !write(0,*)convft(i),sv(2,j)%cpool_fast,sv(l,j)%cpool_fast
    !read(*,*)
    
    osv%tile(2)%cpool_fast = (1. - convft(i)) * osv%tile(2)%cpool_fast + convft(i) * osv%tile(l)%cpool_fast
    osv%tile(2)%cpool_slow = (1. - convft(i)) * osv%tile(2)%cpool_slow + convft(i) * osv%tile(l)%cpool_slow
    
  end if
end do

!reduce the size of the fast and slow wood pools according to exponential decay with fixed turnover time

osv%carbon%prod_flux(1) = (1. - wd(1)) * osv%carbon%wood_fast
osv%carbon%prod_flux(2) = (1. - wd(2)) * osv%carbon%wood_slow

osv%carbon%wood_fast = osv%carbon%wood_fast - osv%carbon%prod_flux(1) !* wd(1)
osv%carbon%wood_slow = osv%carbon%wood_slow - osv%carbon%prod_flux(2) !* wd(2)

!---------------------------------------------

end subroutine alcc

!------------------------------------------------------------------------------------------------------------

subroutine tile_landuse(lutype,recoverf,estab_lim,estab,nind,fpc_grid)

use parametersmod, only : sp,pft

implicit none

integer,  intent(in) :: lutype
real(sp), intent(in) :: recoverf

logical,  dimension(:), intent(in)  :: estab_lim  !whether PFT is within bioclimatic limits for establishment
logical,  dimension(:), intent(out) :: estab      !whether PFT actually establishes

real(sp), dimension(:), intent(inout) :: nind
real(sp), dimension(:), intent(inout) :: fpc_grid

real(sp), dimension(size(nind)) :: nind0

!---------------

estab = estab_lim

if (lutype == 2) then  ! FLAG
  
  !set human land use, do not allow woody vegetation to establish
  !cover the used part of the gridcell with herbaceous vegetation

  where (pft%tree)
    estab = .false.
!  elsewhere
!    fpc_grid = 1.
!    lm_ind(:,1) = lm_sapl
  end where
  
else if (lutype == 3 .and. recoverf > 0.) then  !land abandonment
  
  !decrease the density across woody PFTs to allow more establishment of vegetation
  
  where (nind > 0. .and. pft%tree) 

    nind0 = nind
    nind  = nind - recoverf * nind

    fpc_grid = fpc_grid * nind / nind0

  end where

end if

end subroutine tile_landuse

!------------------

subroutine harvest(i,j,osv)

!Herbaceous PFTs have only leaf and root biomass (no stem).
!The harvest fraction = the leaf biomass. All roots are turned into belowground litter.
!A refinement of this scheme would be to separate some fraction of the harvest amount as actually edible,
!and set the rest to litter or to a burned pool. For now, we assume that all aboveground biomass
!is respired within the current year.

use parametersmod,   only : sp, npft
use mpistatevarsmod, only : statevars

implicit none

integer, intent(in) :: i
integer, intent(in) :: j
type(statevars), target, intent(inout) :: osv

logical,  pointer, dimension(:) :: present
real(sp), pointer, dimension(:) :: lm_ind
real(sp), pointer, dimension(:) :: sm_ind
real(sp), pointer, dimension(:) :: hm_ind
real(sp), pointer, dimension(:) :: rm_ind
real(sp), pointer, dimension(:) :: litter_bg
real(sp), pointer, dimension(:) :: litter_ag_fast

real(sp), pointer :: cropfrac

real(sp), dimension(npft) :: crop_biomass
real(sp), dimension(npft) :: crop_harvest

!----

present        => osv%tile(i)%present        
lm_ind         => osv%tile(i)%lm_ind(:,1)    
sm_ind         => osv%tile(i)%sm_ind(:,1)    
hm_ind         => osv%tile(i)%hm_ind(:,1)    
rm_ind         => osv%tile(i)%rm_ind(:,1)    
cropfrac       => osv%tile(2)%coverfrac
litter_bg      => osv%tile(i)%litter_bg(:,1) 
litter_ag_fast => osv%tile(i)%litter_ag_fast(:,1) 

!-------

!set cflux harvest amount (sum of leafmass across all present PFTs - should be only herbaceous PFTs)

!crop_biomass = sum(lm_ind,mask=present)

where(present) crop_biomass = lm_ind	!do it pft-wise, for the PFTs that are present
  
crop_harvest   = 0.8 * crop_biomass   !80% of the crop biomass is removed from the field

where(present) litter_ag_fast = litter_ag_fast + 0.2 * crop_biomass   !20% is left on the field a stubble, of which 20% is burned in the managedburn subroutine

!remove aboveground living biomass from the cultivated tiles (this way does not go into soil)

!where (present)  !should only be herbaceous pfts

  !lm_ind = 0.

  !belowground biomass is transferred to litter bg

  !litter_bg = litter_bg + rm_ind

  !rm_ind = 0.
  
  !present = .false.

!end where

osv%carbon%crop_harvest = crop_harvest

end subroutine harvest

!------------------

end module alccmod
