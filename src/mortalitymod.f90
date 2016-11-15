module mortalitymod

use parametersmod, only : sp

implicit none

public :: mortality

!parameters

integer, parameter :: npft = 9

real(sp), parameter :: mort_max   =   0.01  !asymptotic maximum mortality rate (/year)
real(sp), parameter :: k_mort     =   0.3   !coefficient of growth efficiency in mortality equation
real(sp), parameter :: ramp_gddtw = 300.    !ramp for heat damage function

!note: above 200 growing degree days above the upper limit tw, establishment is zero and mortality 100%  

contains

!---------------------------------------------------------

subroutine mortality(pftpar,present,tree,boreal,bm_inc,              &
       turnover_ind,sla,lm_ind,sm_ind,hm_ind,rm_ind,nind,year,       &
       litter_ag_fast,litter_ag_slow,litter_bg,dtemp,anpp,mtemp_max)

!arguments

integer,                  intent(in)    :: year

logical,  dimension(:),   intent(in)    :: boreal
logical,  dimension(:),   intent(in)    :: tree
real(sp), dimension(:),   intent(in)    :: sla
real(sp), dimension(:),   intent(in)    :: dtemp
real(sp), dimension(:),   intent(in)    :: turnover_ind

real(sp), dimension(:,:), intent(in)    :: anpp
real(sp), dimension(:,:), intent(in)    :: bm_inc
real(sp), dimension(:,:), intent(in)    :: pftpar

logical,  dimension(:),   intent(inout) :: present
real(sp), dimension(:),   intent(inout) :: nind

real(sp), dimension(:,:), intent(inout) :: lm_ind
real(sp), dimension(:,:), intent(inout) :: sm_ind
real(sp), dimension(:,:), intent(inout) :: hm_ind
real(sp), dimension(:,:), intent(inout) :: rm_ind
real(sp), dimension(:,:), intent(inout) :: litter_ag_fast
real(sp), dimension(:,:), intent(inout) :: litter_ag_slow
real(sp), dimension(:,:), intent(inout) :: litter_bg

!local variables

integer  :: d
integer  :: pft

real(sp) :: gddtw
real(sp) :: mtemp_max
real(sp) :: twmax
real(sp) :: greffic
real(sp) :: bm_delta   !net individual living biomass increment (incorporating loss through leaf, root and sapwood turnover) (gC)
real(sp) :: mort       !tree mortality rate
real(sp) :: nind_kill  !reduction in individual density due to mortality (indiv/m2)
real(sp) :: litter_inc
real(sp) :: heatstress  !reduction in individual density (& establishment) due to heat induced mortality  (indiv/m2)

!---------------------------------------------------------

do pft = 1,npft

  twmax=pftpar(pft,30)  !PFT-specific upper limit of warmest-month temperature

  if (present(pft) .and. tree(pft) .and. nind(pft) > 0.) then
    
    !Calculate net individual living biomass increment

    bm_delta = max(0.,bm_inc(pft,1) / nind(pft) - turnover_ind(pft))

    !Calculate growth efficiency (net biomass increment per unit leaf area)

    greffic = bm_delta / lm_ind(pft,1) / sla(pft)

    !Mortality rate inversely related to growth efficiency (Prentice et al 1993)

    mort = mort_max / (1. + k_mort * greffic)

    !heat damage mortality in boreal trees

    if (mtemp_max > twmax) then  ! heat damage

      !calculate growing degree days above twmax

      gddtw = sum(dtemp - twmax,mask = dtemp > twmax)

      heatstress = min(1.,gddtw / ramp_gddtw)

    else

      heatstress = 0.

    end if

    !write(0,'(i5,4f10.1,f10.6)')pft,nind(pft)*(lm_ind(pft,1)+sm_ind(pft,1) + hm_ind(pft,1)+rm_ind(pft,1)),  &
    !                                  bm_inc(pft,1) / nind(pft),turnover_ind(pft),greffic,mort

    !Reduce individual density (and thereby gridcell-level biomass) by mortality rate

    mort      = min(1.,mort + heatstress)
    nind_kill = nind(pft) * mort
    nind(pft) = max(0., nind(pft) - nind_kill)   

    !Transfer lost biomass to litter

    litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + nind_kill * lm_ind(pft,1)
    litter_ag_slow(pft,1) = litter_ag_slow(pft,1) + nind_kill *(sm_ind(pft,1) + hm_ind(pft,1))
    litter_bg(pft,1)      = litter_bg(pft,1)      + nind_kill * rm_ind(pft,1)
    
  end if
  
  if(nind(pft) == 0.) then
    present(pft)   = .false.
    lm_ind(pft,1)  = 0.
    sm_ind(pft,1)  = 0.
    hm_ind(pft,1)  = 0.
    rm_ind(pft,1)  = 0.
  end if      
     
end do

end subroutine mortality

!---------------------------------------------------------

end module mortalitymod
