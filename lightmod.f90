module lightmod

implicit none

public :: light

contains

!------------------------------------------------------------------------------------------

subroutine light(present,tree,lm_ind,sm_ind,hm_ind,rm_ind,crownarea,fpc_grid,fpc_inc,  &
                 nind,litter_ag_fast,litter_ag_slow,litter_bg,sla,fpc_tree_max)

!recoded in f90 by Jed Kaplan, 04/2010. Should give the same result as the original code if 
!tree competition section is not commented out.

use parametersmod, only : sp,npft

implicit none

!arguments

real(sp), intent(in) :: fpc_tree_max

logical,  dimension(:),   intent(in)    :: tree
real(sp), dimension(:),   intent(in)    :: sla
real(sp), dimension(:),   intent(in)    :: fpc_inc
logical,  dimension(:),   intent(inout) :: present
real(sp), dimension(:),   intent(inout) :: nind
real(sp), dimension(:),   intent(inout) :: fpc_grid
real(sp), dimension(:),   intent(inout) :: crownarea
real(sp), dimension(:,:), intent(inout) :: lm_ind
real(sp), dimension(:,:), intent(inout) :: sm_ind
real(sp), dimension(:,:), intent(inout) :: hm_ind
real(sp), dimension(:,:), intent(inout) :: rm_ind
real(sp), dimension(:,:), intent(inout) :: litter_ag_fast
real(sp), dimension(:,:), intent(inout) :: litter_ag_slow
real(sp), dimension(:,:), intent(inout) :: litter_bg

!local variables

integer :: pft
integer :: ntree         !no. of tree PFTs currently present
integer :: ngrass

real(sp) :: fpc_inc_tree     !this years total FPC increment for tree PFTs
real(sp) :: fpc_tree_total   !total grid FPC for tree PFTs      
real(sp) :: fpc_grass_max    !max allowed grass FPC given the current tree cover
real(sp) :: fpc_grass_total  !total grid FPC for grass PFTs
real(sp) :: grasscover       !grass PFT proportional cover ("crown area")
real(sp) :: excess           !total tree FPC or grass cover to be reduced

real(sp) :: nind_kill        !reduction in individual density to reduce tree FPC to permitted maximum (indiv/m2)
real(sp) :: rm_kill          !reduction in grass PFT root mass to reduce grass cover to permitted maximum (gC)  
real(sp) :: lm_kill          !reduction in grass PFT leaf mass to reduce grass cover to permitted maximum (gC)
real(sp) :: lm_old

real(sp), dimension(npft) :: fpc_ind
real(sp), dimension(npft) :: lai_ind
real(sp), dimension(npft) :: pft_excess

integer,  dimension(npft) :: pftsorted

integer :: i

!----------------------------------------------------------------------------------
!calculate total woody FPC, FPC increment and grass cover (= crown area)

where (crownarea > 0.) 

  lai_ind  = lm_ind(:,1) * sla / crownarea
  fpc_ind  = 1. - exp(-0.5 * lai_ind)
  fpc_grid = fpc_ind * nind * crownarea

elsewhere

  lai_ind  = 0.
  fpc_ind  = 0.
  fpc_grid = 0.

end where

ntree           = count(present .and. tree)
fpc_inc_tree    = sum(fpc_inc,  mask = present .and. tree)
fpc_tree_total  = sum(fpc_grid, mask = present .and. tree)

ngrass          = count(present .and. .not. tree)
grasscover      = sum(crownarea, mask = present .and. .not. tree)
fpc_grass_total = sum(fpc_grid,  mask = present .and. .not. tree)

!-------------------------------
!light competition, woody plants

pft_excess = 0.

!write(0,*)'light',fpc_tree_total,fpc_tree_max

if (fpc_tree_total > fpc_tree_max) then  !reduce tree cover

  excess = fpc_tree_total - fpc_tree_max
  
  do pft = 1,npft
    if (present(pft) .and. tree(pft) .and. fpc_grid(pft) > 0.) then

      !this formulation ensures equal competition (precludes total dominance by one PFT)

      pft_excess(pft) = min(fpc_grid(pft),excess *  fpc_grid(pft) / fpc_tree_total)

      !original LPJ formulation allows one PFT to become dominant if it has no fpc_inc (so the others are reduced)

      !if (fpc_inc_tree > 0.) then
      !  pft_excess(pft) = min(fpc_grid(pft),excess * (fpc_inc(pft) / fpc_inc_tree))
      !else
      !  pft_excess(pft) = min(fpc_grid(pft),excess / real(ntree))
      !end if

    else
      pft_excess(pft) = 0.
    end if
    
  end do

  do pft = 1,npft


    if (pft_excess(pft) > 0.) then
      
      !Reduce individual density (and thereby gridcell-level biomass) so that total tree FPC reduced to 'fpc_tree_max'

      nind_kill = nind(pft) * pft_excess(pft) / fpc_grid(pft)

      nind(pft) = nind(pft) - nind_kill

      !Transfer lost biomass to litter

      litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + nind_kill * lm_ind(pft,1)                  !leaves
      litter_ag_slow(pft,1) = litter_ag_slow(pft,1) + nind_kill *(sm_ind(pft,1) + hm_ind(pft,1)) !stems
      litter_bg(pft,1)      = litter_bg(pft,1)      + nind_kill * rm_ind(pft,1)                  !roots
      
      !zero out any pft that has been reduced to zero nind
      
      if (nind(pft) <= 0.) then
        
        present(pft)  =.false.
        nind(pft)     = 0.
        lm_ind(pft,1) = 0.
        sm_ind(pft,1) = 0.
        hm_ind(pft,1) = 0.
        rm_ind(pft,1) = 0.

      end if
      
      !update isotopes      
      !ignored for now

    end if
  
  end do
    
end if

!-----------------
!grass competition

fpc_grass_max = 1. - min(fpc_tree_total,fpc_tree_max)

if (fpc_grass_total > fpc_grass_max) then  !reduce grass cover

  excess = fpc_grass_total - fpc_grass_max

  do pft = 1,npft
    
    if (present(pft) .and. .not. tree(pft) .and. lm_ind(pft,1) > 0.) then
        
      lm_old = lm_ind(pft,1)

      lm_ind(pft,1) = max(0., -2. * log(1. - (fpc_grid(pft) - excess)) / sla(pft))

      lm_kill = lm_old - lm_ind(pft,1)
      rm_kill = rm_ind(pft,1) * lm_kill / lm_old

      rm_ind(pft,1) = rm_ind(pft,1) - rm_kill

      !Transfer lost biomass to litter
      litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + lm_kill
      litter_bg(pft,1)      = litter_bg(pft,1)      + rm_kill

      !update isotopes
      !ignored for now

    end if

  end do
  
end if

!-----------------

where (crownarea > 0.) 
  lai_ind  = lm_ind(:,1) * sla / crownarea
  fpc_ind  = 1. - exp(-0.5 * lai_ind)
  fpc_grid = fpc_ind * nind * crownarea
elsewhere
  lai_ind  = 0.
  fpc_ind  = 0.
  fpc_grid = 0.
end where

!correct for mathematical overshoot
      
litter_ag_fast = max(litter_ag_fast,0.)
litter_ag_slow = max(litter_ag_slow,0.)
litter_bg      = max(litter_bg,0.)

!do pft=1,npft
!  if (fpc_grid(pft) < 0. .or. fpc_grid(pft) > 1.) write(0,*)'resetting pft ',pft,fpc_grid(pft)
!end do

!fpc_grid = max(fpc_grid,0.)
!fpc_grid = min(fpc_grid,1.)

end subroutine light

!----------------------------------

end module lightmod
