module turnovermod

use parametersmod, only : sp

implicit none

public :: turnover

contains

!---------------------------------------------------------

subroutine turnover(pftpar,present,lm_ind,sm_ind, hm_ind,rm_ind,litter_ag_fast,litter_ag_slow,litter_bg,nind,turnover_ind,year)

!Turnover of PFT-specific fraction from each living C pool, leaf and root C transferred to litter, sapwood C to heartwood

use parametersmod, only : npft

implicit none

integer, intent(in) :: year

real(sp), dimension(:,:), intent(in)    :: pftpar

real(sp), dimension(:), intent(out) :: turnover_ind

logical,  dimension(:),   intent(inout) :: present

real(sp), dimension(:,:), intent(inout) :: lm_ind
real(sp), dimension(:,:), intent(inout) :: sm_ind
real(sp), dimension(:,:), intent(inout) :: hm_ind
real(sp), dimension(:,:), intent(inout) :: rm_ind
real(sp), dimension(:,:), intent(inout) :: litter_ag_fast
real(sp), dimension(:,:), intent(inout) :: litter_ag_slow
real(sp), dimension(:,:), intent(inout) :: litter_bg
real(sp), dimension(:),   intent(inout) :: nind

!parameters

real(sp), parameter :: small = 1.

!local variables

integer  :: pft
real(sp) :: l_torate
real(sp) :: s_torate
real(sp) :: r_torate
real(sp) :: lm_turn
real(sp) :: sm_turn
real(sp) :: rm_turn

!--------------------------

do pft=1,npft

  if (.not.present(pft)) cycle
       
  !Turnover rates are reciprocals of tissue longevity

  l_torate = 1. / pftpar(pft,8)
  s_torate = 1. / pftpar(pft,9)
  r_torate = 1. / pftpar(pft,10)

  !Calculate the biomass turnover in this year

  lm_turn = lm_ind(pft,1) * l_torate
  sm_turn = sm_ind(pft,1) * s_torate
  rm_turn = rm_ind(pft,1) * r_torate

  !Update the pools

  lm_ind(pft,1) = lm_ind(pft,1) - lm_turn
  sm_ind(pft,1) = sm_ind(pft,1) - sm_turn
  rm_ind(pft,1) = rm_ind(pft,1) - rm_turn 

  !Convert sapwood to heartwood

  hm_ind(pft,1) = hm_ind(pft,1) + sm_turn

  !Transfer to litter pools

  litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + lm_turn * nind(pft)

  litter_bg(pft,1) = litter_bg(pft,1) + rm_turn * nind(pft)

  !Record total turnover

  turnover_ind(pft) = lm_turn + sm_turn + rm_turn

end do

end subroutine turnover

!---------------------------------------------------------

end module turnovermod
