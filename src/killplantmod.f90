module killplantmod

use parametersmod, only : sp

implicit none

public :: killplant

integer,  parameter :: npft = 9

contains

!---------------------------------------------------------

subroutine killplant(bm_inc,present,tree,crownarea,lm_ind,rm_ind,hm_ind,sm_ind,nind,litter_ag_fast,litter_ag_slow,litter_bg)

!Removal of PFTs with negative annual C increment so that there are not problems in allocation
!Note: Killing of PFTs newly beyond their bioclimatic limits is done in subroutine establishment

implicit none

!arguments

logical,  dimension(:),   intent(in)    :: tree
real(sp), dimension(:,:), intent(in)    :: bm_inc
real(sp), dimension(:), intent(in)    :: crownarea

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

real(sp), parameter :: bminc_ind_min = 1.e-6 ! minimum annual productivity per individual has to be more than this value (gC)

!local variables

integer  :: pft
real(sp) :: bm_inc_ind

! --------------------------------------
! we cannot use this subroutine for herbaceous PFTs because the long recovery time in the model
! means that the recovery in leafmass is too slow (several years) compared to within one year
! in reality

do pft = 1,npft

  if (present(pft) .and. tree(pft)) then  

    ! bm_inc_ind = bm_inc(pft,1) / nind(pft)
    bm_inc_ind = bm_inc(pft,1) * crownarea(pft)
    
    ! because of leaf turnover, it is possible that lm_ind = 0 at this stage, so we cannot use a leafmass threshold for killing the plant

    ! if (bm_inc_ind < bminc_ind_min .or. lm_ind(pft,1) < 1.) then  ! does not work

    if (bm_inc_ind < bminc_ind_min) then   ! original LPJ formulation

      !not enough C increment this year, kill PFT and transfer carbon to litter
      
      ! write(0,*)'killing pft',pft,bm_inc_ind,crownarea(pft),lm_ind(pft,1)

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
