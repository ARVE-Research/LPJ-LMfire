module killplantmod

use parametersmod, only : sp

implicit none

public :: killplant

integer,  parameter :: npft = 9

contains

!---------------------------------------------------------

subroutine killplant(bm_inc,present,tree,lm_ind,rm_ind,hm_ind,sm_ind,nind,litter_ag_fast,litter_ag_slow,litter_bg)

!Removal of PFTs with negative annual C increment so that there are not problems in allocation
!Note: Killing of PFTs newly beyond their bioclimatic limits is done in subroutine establishment

implicit none

!arguments

logical,  dimension(:),   intent(in)    :: tree
real(sp), dimension(:,:), intent(in)    :: bm_inc

logical,  dimension(:),   intent(inout) :: present
real(sp), dimension(:),   intent(inout) :: nind
real(sp), dimension(:,:), intent(inout) :: lm_ind
real(sp), dimension(:,:), intent(inout) :: sm_ind
real(sp), dimension(:,:), intent(inout) :: hm_ind
real(sp), dimension(:,:), intent(inout) :: rm_ind
real(sp), dimension(:,:), intent(inout) :: litter_ag_fast
real(sp), dimension(:,:), intent(inout) :: litter_ag_slow
real(sp), dimension(:,:), intent(inout) :: litter_bg

!parameter

real(sp), parameter :: bminc_ind_min = 1e-6 !minimum annual productivity per individual has to be more than this value (gC)

!local variables

integer  :: pft
real(sp) :: bm_inc_ind

!--------------------------------------

do pft = 1,npft

  if (present(pft)) then

    bm_inc_ind = bm_inc(pft,1) / nind(pft)
    
    if (bm_inc_ind < bminc_ind_min) then

      !not enough C increment this year, kill PFT and transfer carbon to litter

      !all PFTs leaf and root mass

      litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + lm_ind(pft,1) * nind(pft)
      litter_bg(pft,1)      = litter_bg(pft,1)      + rm_ind(pft,1) * nind(pft)

      lm_ind(pft,1) = 0.
      rm_ind(pft,1) = 0.

      if (tree(pft)) then  !also remove the stem mass

        litter_ag_slow(pft,1) = litter_ag_slow(pft,1) + (sm_ind(pft,1) + hm_ind(pft,1)) * nind(pft)

        sm_ind(pft,1) = 0.
        hm_ind(pft,1) = 0.

      end if

      !reset nind and present

      nind(pft)    = 0.
      present(pft) =.false.
      lm_ind(pft,1) = 0.
      sm_ind(pft,1) = 0.
      hm_ind(pft,1) = 0.
      rm_ind(pft,1) = 0.

    end if
  end if

end do

end subroutine killplant

!---------------------------------------------------------

end module killplantmod
