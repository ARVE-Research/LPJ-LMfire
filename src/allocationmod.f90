module allocationmod

use parametersmod, only : sp,dp

implicit none

public  :: allocation
private :: root

integer,  parameter :: npft = 9
integer,  parameter :: nseg = 20
real(sp), parameter :: pi   =  3.14159265
real(sp), parameter :: xacc =  0.1     !x-axis precision threshold for the allocation solution
real(dp), parameter :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution

!NB The final result is pretty sensitive to the accuracy values xacc and yacc
!NB See extensive comments and derivation at end of this code section

contains

!---------------------------------------------------------

subroutine allocation(pftpar,allom1,allom2,allom3,latosa,wooddens,                  &
                      reinickerp,tree,sla,wscal,nind,bm_inc,lm_ind,sm_ind,hm_ind,   &
                      rm_ind,crownarea,fpc_grid,lai_ind,height,litter_ag_fast,      &
                      litter_ag_slow,litter_bg,fpc_inc,present)

implicit none

!arguments

real(sp), intent(in) :: allom1
real(sp), intent(in) :: allom2
real(sp), intent(in) :: allom3
real(sp), intent(in) :: latosa
real(sp), intent(in) :: wooddens
real(sp), intent(in) :: reinickerp

logical,  dimension(:),   intent(in)    :: present
logical,  dimension(:),   intent(in)    :: tree
real(sp), dimension(:),   intent(in)    :: sla
real(sp), dimension(:),   intent(in)    :: wscal

real(sp), dimension(:,:), intent(in)    :: pftpar
real(sp), dimension(:,:), intent(in)    :: bm_inc

real(sp), dimension(:),   intent(out)   :: fpc_inc

real(sp), dimension(:),   intent(inout) :: crownarea
real(sp), dimension(:),   intent(inout) :: fpc_grid
real(sp), dimension(:),   intent(inout) :: height
real(sp), dimension(:),   intent(inout) :: lai_ind
real(sp), dimension(:),   intent(inout) :: nind

real(sp), dimension(:,:), intent(inout) :: lm_ind
real(sp), dimension(:,:), intent(inout) :: sm_ind
real(sp), dimension(:,:), intent(inout) :: hm_ind
real(sp), dimension(:,:), intent(inout) :: rm_ind
real(sp), dimension(:,:), intent(inout) :: litter_ag_fast
real(sp), dimension(:,:), intent(inout) :: litter_ag_slow
real(sp), dimension(:,:), intent(inout) :: litter_bg

!local variables

logical  :: normal
integer  :: pft

real(sp) :: bm_inc_ind     !individual total biomass increment this year
real(sp) :: crownarea_max  !maximum crown area (m2)
real(sp) :: fpc_grid_old   !previous year's FPC
real(sp) :: fpc_ind        !individual FPC
real(sp) :: lm2rm          !ratio of leafmass to fine rootmass
real(sp) :: lminc_ind      !individual leafmass increment this year
real(sp) :: rminc_ind      !individual fineroot mass increment this year
real(sp) :: lminc_ind_min  !min leafmass increment to maintain current sapwood
real(sp) :: rminc_ind_min  !min rootmass increment to support new leafmass
real(sp) :: sap_xsa        !cross sectional area of sapwood  
real(sp) :: sminc_ind      !individual sapmass increment this year
real(sp) :: stemdiam       !stem diameter 

real(sp) :: lm
real(sp) :: sm
real(sp) :: hm
real(sp) :: rm

real(sp) :: x1             !working vars in bisection
real(sp) :: x2
real(sp) :: rtbis
real(sp) :: dx
real(sp) :: xmid
real(sp) :: sign

real(dp) :: fx1
real(dp) :: fmid

real(sp) :: lm1     !allometric leafmass requirement (leafmass req'd to keep sapwood alive; gC ind-1)

integer :: i

!--------------------------------------
!The goal of the allocation routine is to partition the annual biomass increment (bm_inc_ind)
!over the three living carbon pools (leaf, root, sapwood) in a way such that the basic allometric
!relationships are satisfied. This is accomplished by calculation of a leaf mass increment
!(lminc_ind) that satisfies Eqn. 22 (see bottom of code file for mathematical derivation).

!The values for leafmass and rootmass going into this subroutine can be zero because they will have been
!turned over (e.g., deciduous leaves), so the annual biomass increment will first have to be used to put
!on new leaf and root mass to support the sapwood that is standing. Additional carbon gain will be allocated
!to sapwood, which increases tree height.

do pft = 1,npft

  if (.not.present(pft)) cycle
  
  lm = lm_ind(pft,1)
  sm = sm_ind(pft,1)
  hm = hm_ind(pft,1)
  rm = rm_ind(pft,1)
  
  bm_inc_ind = bm_inc(pft,1) / nind(pft)

  !calculate this year's leaf to fine root mass ratio from mean annual water scalar and pft specific parameter
  !minumum value set at 10%, values lower than this are not generally found in nature

  lm2rm = max(pftpar(pft,16) * wscal(pft),0.1)

  if (tree(pft)) then

    !TREE ALLOCATION

    lm1 = latosa * sm / (wooddens * height(pft) * sla(pft))  !allometric leaf mass requirement

    lminc_ind_min = lm1 - lm  !eqn (27)

    !calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
    !May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.

    rminc_ind_min = lm1 / lm2rm - rm      !eqn (30)

    if (rminc_ind_min > 0. .and. lminc_ind_min > 0. .and. rminc_ind_min + lminc_ind_min <= bm_inc_ind) then

      !Normal allocation (positive increment to all living C compartments)

      normal = .true.

      !Calculation of leaf mass increment (lminc_ind) that satisfies Eqn (22)
      !Since this is normal allocation, we set the lower bound for the leafmass allocation (x1)
      !to its allometric minimum, because it should be able to be fulfilled, i.e.:

      x1 = lminc_ind_min
      x2 = (bm_inc_ind - (lm / lm2rm - rm)) / (1. + 1. / lm2rm)
      
      dx = x2 - x1
      
      if (dx < 0.01) then

        !there seems to be rare cases where lminc_ind_min (x1) is almost equal to x2. In this case,
        !assume that the leafmass increment is equal to the midpoint between the values and skip 
        !the root finding procedure

        lminc_ind = x1 + 0.5 * dx

      else

        !Find a root for non-negative lminc_ind, rminc_ind and sminc_ind using Bisection Method (Press et al 1986, p 346)
        !There should be exactly one solution (no proof presented, but Steve has managed one).
              
        dx = dx / real(nseg)

        !evaluate f(x1) = LHS of eqn (22) at x1

        fx1 = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla(pft),x1)

        !Find approximate location of leftmost root on the interval (x1,x2).
        !Subdivide (x1,x2) into nseg equal segments seeking change in sign of f(xmid) relative to f(x1).

        fmid = fx1
        xmid = x1

        i = 1

        do

          xmid = xmid + dx

          fmid = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla(pft),xmid)

          if (fmid * fx1 <= 0. .or. xmid >= x2) exit  !sign has changed or we are over the upper bound

          if (i > 20) write(0,*)'first alloc loop flag',i,pft,fmid*fx1,xmid,x1,x2,dx,bm_inc_ind
          if (i > 50) stop 'Too many iterations allocmod'

          i = i + 1

        end do

        !the interval that brackets zero in f(x) becomes the new bounds for the root search

        x1 = xmid - dx
        x2 = xmid

        !Apply bisection method to find root on the new interval (x1,x2)

        fx1 = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla(pft),x1)

        if (fx1 >= 0.) then
          sign = -1.
        else
          sign =  1.
        end if

        rtbis = x1
        dx    = x2 - x1

        !Bisection loop: search iterates on value of xmid until xmid lies within xacc of the root,
        !i.e. until |xmid-x| < xacc where f(x) = 0. the final value of xmid with be the leafmass increment

        i = 1

        do 

          dx   = 0.5 * dx
          xmid = rtbis + dx

          !calculate fmid = f(xmid) [eqn (22)]

          fmid = root(lm,sm,hm,rm,bm_inc_ind,lm2rm,sla(pft),xmid)

          if (fmid * sign <= 0.) rtbis = xmid

          if (dx < xacc .or. abs(fmid) <= yacc) exit

          if (i > 20) write(0,*)'second alloc loop flag',i,pft,dx,abs(fmid)
          if (i > 50) stop 'Too many iterations allocmod'

          i = i + 1

        end do

        !Now rtbis contains numerical solution for lminc_ind given eqn (22)

        lminc_ind = rtbis
      
      end if  !x2-x1 block
      
      !Calculate increments in other compartments using allometry relationships

      rminc_ind = (lm + lminc_ind) / lm2rm - rm       !eqn (9)

      sminc_ind = bm_inc_ind - lminc_ind - rminc_ind  !eqn (1)

    else
      
      !Abnormal allocation: reduction in some C compartment(s) to satisfy allometry

      normal = .false.

      !Attempt to distribute this year's production among leaves and roots only

      lminc_ind = (bm_inc_ind - lm / lm2rm + rm) / (1. + 1. / lm2rm)  !eqn (33)

      if (lminc_ind > 0.) then

        !Positive allocation to leafmass

        rminc_ind = bm_inc_ind - lminc_ind  !eqn (31)

        !Add killed roots (if any) to below-ground litter

        if (rminc_ind < 0.) then

          lminc_ind = bm_inc_ind
          rminc_ind = (lm + lminc_ind) / lm2rm - rm

          litter_bg(pft,1) = litter_bg(pft,1) + abs(rminc_ind) * nind(pft)

        end if
        
        i = 1

      else

        !Negative allocation to leaf mass

        rminc_ind = bm_inc_ind
        lminc_ind = (rm + rminc_ind) * lm2rm - lm  !from eqn (9)

        !Add killed leaves to litter

        litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + abs(lminc_ind) * nind(pft)
        
        i = 2

      end if

      !Calculate sminc_ind (must be negative)
      
      sminc_ind = (lm + lminc_ind) * sla(pft) / latosa * wooddens * height(pft) - sm  !eqn (35)

      !Convert killed sapwood to heartwood

      hm = hm + abs(sminc_ind)
      
      !write(0,*)'abnormal case',i,lminc_ind,rminc_ind,lminc_ind+rminc_ind,bm_inc_ind,sminc_ind

    end if  !normal/abnormal allocation

    !Increment C compartments

    lm_ind(pft,1) = lm + lminc_ind
    rm_ind(pft,1) = rm + rminc_ind
    sm_ind(pft,1) = sm + sminc_ind
    hm_ind(pft,1) = hm

    !Calculate new height, diameter and crown area

    if (lm_ind(pft,1) > 0.) then

      sap_xsa = lm_ind(pft,1) * sla(pft) / latosa  !eqn (5)

      height(pft) = sm_ind(pft,1) / sap_xsa / wooddens

      stemdiam    = (height(pft) / allom2)**(1./allom3)                  !eqn (C)

      crownarea_max  = pftpar(pft,18)

      crownarea(pft) = min(allom1 * stemdiam**reinickerp,crownarea_max)  !eqn (D)

    end if

  else

    !GRASS ALLOCATION

    !Distribute this year's production among leaves and fine roots according to leaf to rootmass ratio [eqn (33)] (see below)
    !Relocation of C from one compartment to the other not allowed: negative increment in either compartment transferred to litter
    !but the total negative amount cannot be more than the existing pool plus the increment

    lminc_ind = (bm_inc_ind - lm / lm2rm + rm) / (1. + 1. / lm2rm)

    rminc_ind = bm_inc_ind - lminc_ind
    
    if (lminc_ind > 0.) then
 
      if (rminc_ind < 0.) then  !negative allocation to root mass

        if (rminc_ind + rm < 0.) rminc_ind = -rm  !cannot be more than the root mass that is actually present

        !Add killed roots to below-ground litter

        litter_bg(pft,1) = litter_bg(pft,1) + abs(rminc_ind) * nind(pft)

      end if

    else

      !Negative allocation to leaf mass

      rminc_ind = bm_inc_ind
      lminc_ind = lm2rm * (rm + rminc_ind) - lm  !from eqn (9)
      
      if (lminc_ind > 0.) lminc_ind = -lm  !cannot be more than the leaf mass that is actually present

      !Add killed leaf mass to litter

      litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + abs(lminc_ind) * nind(pft)

    end if

    !Increment C compartments

    lm_ind(pft,1) = lm + lminc_ind
    rm_ind(pft,1) = rm + rminc_ind

  end if  !tree/grass

  !Update LAI and FPC

  if (crownarea(pft) > 0.) then
    lai_ind(pft) = (lm_ind(pft,1) * sla(pft)) / crownarea(pft)
  else
    lai_ind(pft) = 0.
  end if

  fpc_grid_old  = fpc_grid(pft)
  fpc_ind       = 1. - exp(-0.5 * lai_ind(pft))
  fpc_grid(pft) = crownarea(pft) * nind(pft) * fpc_ind
  fpc_inc(pft)  = max(fpc_grid(pft) - fpc_grid_old,0.)

end do  !pft loop

end subroutine allocation

!---------------------------------------------------------

real(dp) function root(lm,sm,hm,rm,inc,lm2rm,sla,x)

use parametersmod, only : pi,allom2,allom3,latosa,wooddens

implicit none

!parameters

real(sp), parameter :: pi4 = pi / 4.
real(sp), parameter :: a1  = 2. / allom3
real(sp), parameter :: a2  = 1. + a1
real(sp), parameter :: a3  = allom2**a1

!arguments

real(sp), intent(in) :: sla    !specific leaf area
real(sp), intent(in) :: lm2rm  !leaf mass to root mass ratio
real(sp), intent(in) :: lm     !individual leaf mass
real(sp), intent(in) :: sm     !individual sapwood mass
real(sp), intent(in) :: hm     !individual heartwood mass
real(sp), intent(in) :: rm     !individual root mass
real(sp), intent(in) :: inc    !individual biomass increment
real(sp), intent(in) :: x      !leafmass allocation amount as input

!---

root = a3 * ((sm + inc - x - ((lm + x) / lm2rm) + rm + hm) / wooddens) / pi4 - &
            ((sm + inc - x - ((lm + x) / lm2rm) + rm) / ((lm + x) * sla * wooddens / latosa))**a2

end function root  

!---------------------------------------------------------
!NOTES

!(A) (leaf area) = latosa * (sapwood xs area)
!      (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)

!(B) (leaf mass) = lm2rm * (root mass)

!(C) height = allom2 * (stem diameter)**allom3
!      (source?)

!(D) (crown area) = min (allom1 * (stem diameter)**reinickerp,crownarea_max)

!Mathematical derivation:   

!  (1) bm_inc_ind = lminc_ind + sminc_ind + rminc_ind
!  (2) leaf_area_new = latosa * sap_xsa_new   [from (A)]
!  (3) leaf_area_new = (lm_ind + lminc_ind) * sla

!from (2) & (3),           
!  (4) (lm_ind + lminc_ind) * sla = latosa * sap_xsa_new 

!from (4),
!  (5) sap_xsa_new = (lm_ind + lminc_ind) * sla / latosa
!  (6) (lm_ind + lminc_ind) = lm2rm * (rm_ind + rminc_ind) [from (B)]
!  (7) height_new = allom2 * stemdiam_new**allom3  [from (C)]

!from (1),  
!  (8) sminc_ind = bm_inc_ind - lminc_ind - rminc_ind 

!from (6), 
!  (9) rminc_ind=((lm_ind + lminc_ind) / lm2rm) - rm_ind

!from (8) & (9),
! (10) sminc_ind = bm_inc_ind - lminc_ind - ((lm_ind + lminc_ind)  / lm2rm) + rm_ind  
! (11) wooddens = (sm_ind + sminc_ind + hm_ind) / stemvolume_new 
! (12) stemvolume_new = height_new * pi * stemdiam_new**2 / 4  

!from (10), (11) & (12)                   
! (13) stemdiam_new = [ ((sm_ind + bm_inc_ind - lminc_ind
!        - ((lm_ind + lminc_ind) / lm2rm) + rm_ind + hm_ind)
!        / wooddens) / (height_new * pi / 4) ]**(1/2)

!combining (7) and (13),
! (14) height_new = allom2 * [ ((sm_ind + bm_inc_ind - lminc_ind
!        - ((lm_ind + lminc_ind) / lm2rm) + rm_ind + hm_ind) 
!        / wooddens) / (height_new * pi / 4) ]**(1/2 * allom3)

!from (14),
! (15) height_new**(1 + 2 / allom3) = allom2**(2 / allom3)
!        * ((sm_ind + bm_inc_ind - lminc_ind - ((lm_ind + lminc_ind) 
!        / lm2rm) + rm_ind + hm_ind) / wooddens) / (pi / 4)
! (16) wooddens = (sm_ind + sminc_ind) / sapvolume_new 

!from (10) and (16),
! (17) wooddens = (sm_ind + bm_inc_ind - lminc_ind 
!        - ((lm_ind + lminc_ind) / lm2rm) + rm_ind) / sapvolume_new 
! (18) sapvolume_new = height_new * sap_xsa_new

!from (17) and (18),
! (19) sap_xsa_new = (sm_ind + bm_inc_ind - lminc_ind
!        - ((lm_ind + lminc_ind) / lm2rm) + rm_ind)
!        / (height_new * wooddens) 

!from (19),
! (20) height_new = (sm_ind + bm_inc_ind - lminc_ind
!        - ((lm_ind + lminc_ind) / lm2rm) + rm_ind )
!        / (sap_xsa_new * wooddens)

!from (5) and (20),
! (21) height_new**(1 + 2 / allom3) = [ (sm_ind + bm_inc_ind 
!        - lminc_ind - ((lm_ind + lminc_ind) / lm2rm) + rm_ind )
!        / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
!        **(1 + 2 / allom3)

!-------------------------------------------------------------------

! (15) and (21) are two alternative expressions for
!      height_new**(1 + 2 / allom3). Combining these,

! (22) allom2**(2 / allom3) * ((sm_ind + bm_inc_ind - lminc_ind
!        - ((lm_ind + lminc_ind) / lm2rm) + rm_ind + hm_ind)
!        / wooddens) / (pi / 4) - [ (sm_ind + bm_inc_ind - lminc_ind
!        - ((lm_ind + lminc_ind) / lm2rm) + rm_ind )
!        / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
!        **(1 + 2 / allom3)
!        = 0

!Equation (22) can be expressed in the form f(lminc_ind) = 0
!now contained in the function root

!Numerical methods are used to solve the equation for the unknown lminc_ind.

!---------------------------------------------------------

!calculate minimum leaf production to maintain current sapwood mass

! (23) sap_xsa = sm_ind / wooddens / height

!from (A) and (23),
! (24) leaf_mass * sla = latosa * sap_mass / wooddens / height

!from (24),
! (25) leaf_mass = latosa * sap_mass / (wooddens * height * sla)

!from (25), assuming sminc_ind=0,
! (26) lm_ind + lminc_ind_min = latosa * sm_ind / (wooddens * height * sla)

!from (26),
! (27) lminc_ind_min = latosa * sm_ind / (wooddens * height * sla) - lm_ind

!from (B) and (25),
! (28) root_mass = latosa * sap_mass / (wooddens * height * sla) / lm2rm

!from (28), assuming sminc_ind=0,
! (29) rm_ind + rminc_ind_min = latosa * sm_ind / (wooddens * height * sla * lm2rm)

!from (29),
! (30) rminc_ind_min = latosa * sm_ind / (wooddens * height * sla * lm2rm) - rm_ind

!---------------------------------------------------------

!Abnormal allocation: reduction in some C compartment(s) to satisfy allometry
!Attempt to distribute this year's production among leaves and roots only

! (31) bm_inc_ind = lminc_ind + rminc_ind

!from (31) and (9),
! (32) bm_inc_ind = lminc_ind + ((lm_ind + lminc_ind) / lm2rm) - rm_ind

!from (32)
! (33) lminc_ind = (bm_inc_ind - lmind / lm2rm + rm_ind) / (1 + 1 / lm2rm)

!from (25),
! (34) lm_ind + lminc_ind = latosa * (sm_ind + sminc_ind) / (wooddens * height * sla)

!from (34),
! (35) sminc_ind = (lm_ind + lminc_ind) * wooddens * height * sla / latosa - sm_ind

end module allocationmod
