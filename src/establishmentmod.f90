module establishmentmod

use parametersmod, only : sp

implicit none

public :: establishment

integer, parameter :: npft = 9

real(sp), parameter :: pi               =   3.14159265

contains

!---------------------------------------------------------

subroutine establishment(pftpar,present,survive,estab,nind,lm_ind,sm_ind,rm_ind,hm_ind,lm_sapl,sm_sapl,rm_sapl,hm_sapl, &
                         tree,crownarea,fpc_grid,lai_ind,height,sla,wooddens,latosa,mprec,reinickerp, &
                         litter_ag_fast,litter_ag_slow,litter_bg,allom1,allom2,allom3,acflux_estab,   &
                         leafondays,leafoffdays,leafon,estab_pft,burnedf, clay, sand)

!Establishment of new individuals (saplings) of woody PFTs, grass establishment,
!removal of PFTs not adapted to current climate, update of individual structure and FPC.

use parametersmod, only : sp

implicit none

!arguments

real(sp), dimension(:), intent(inout) :: clay   !soil clay fraction
real(sp), dimension(:), intent(inout) :: sand   !soil sand fraction

real(sp), intent(in) :: allom1
real(sp), intent(in) :: allom2
real(sp), intent(in) :: allom3
real(sp), intent(in) :: latosa
real(sp), intent(in) :: wooddens
real(sp), intent(in) :: reinickerp
real(sp), intent(in) :: burnedf

logical,  dimension(:),   intent(in)    :: survive
logical,  dimension(:),   intent(in)    :: estab
logical,  dimension(:),   intent(in)    :: tree
real(sp), dimension(:),   intent(in)    :: sla
real(sp), dimension(:),   intent(in)    :: mprec

real(sp), dimension(:,:), intent(in)    :: pftpar
real(sp), dimension(:,:), intent(in)    :: lm_sapl
real(sp), dimension(:,:), intent(in)    :: sm_sapl
real(sp), dimension(:,:), intent(in)    :: hm_sapl
real(sp), dimension(:,:), intent(in)    :: rm_sapl

logical,  dimension(:),   intent(inout) :: present
logical,  dimension(:),   intent(inout) :: leafon
integer,  dimension(:),   intent(inout) :: leafondays
integer,  dimension(:),   intent(inout) :: leafoffdays
real(sp), dimension(:),   intent(inout) :: crownarea
real(sp), dimension(:),   intent(inout) :: estab_pft
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

real(sp), dimension(:),   intent(out)   :: acflux_estab

!parameters
real(sp), parameter :: aprec_min_estab  = 100.          !minimum annual precipitation for establishment (mm)
!real(sp), dimension(npft), parameter :: aprec_min_estab  = [ 100.,100.,100.,100.,100.,100.,100.,0.,0.]          !minimum annual precipitation for establishment (mm)
real(sp), parameter :: estab_max        =   0.15        !maximum sapling establishment rate (indiv/m2)
real(sp), parameter :: nind_min         =   1.e-10      !minimum individual density for persistence of PFT (indiv/m2)
real(sp), parameter :: eps              =   1.e-6       !epsilon parameter (?)

!local variables

integer  :: pft              !counter
integer  :: npft_estab       !number of regenerating tree PFTs
integer  :: ngrass           !number of grasses present
real(sp) :: aprec            !annual precipitation (mm)
real(sp) :: bare             !gridcell bare ground fraction
real(sp) :: crownarea_max    !maximum crown area (m2) 
real(sp) :: estab_grid       !grid-level establishment rate (indiv/m2)
real(sp) :: estab_mass       !mass of establishing PFT (gC ind-1)
real(sp) :: estab_rate       !sapling establishment rate over area available for establishment (indiv/m2)
real(sp) :: fpc_ind          !individual FPC
real(sp) :: fpc_total        !total grid FPC
real(sp) :: fpc_tree_total   !total grid FPC for tree PFTs
real(sp) :: nind0            !number of individuals /m2 before establishment
real(sp) :: sm_ind_tmp       !preliminary sapwood mass (g/m2)
real(sp) :: stemdiam         !stem diameter (m)

real(sp) :: clay_mean ! Moyenne du clay dans les differents tiles
real(sp) :: sand_mean ! Moyenne du clay dans les differents tiles

!--------------------------------------------
!write(0,*)'burnedf',burnedf

! Faire la moyenne du pourcentage de clay et de sand pour les deux couches 
clay_mean = (clay(1) + clay(3))/2
sand_mean = (sand(1) + sand(3))/2

!Kill PFTs not adapted to current climate, introduce newly "adapted" PFTs

aprec = sum(mprec)  !Calculate annual precipitation

do pft = 1,npft
        
  if (present(pft) .and. (.not.survive(pft) .or. nind(pft) < nind_min)) then

    !kill PFT

    present(pft) = .false.

    !Add killed biomass to litter

    litter_ag_fast(pft,1) = litter_ag_fast(pft,1) + nind(pft) * lm_ind(pft,1)
    litter_ag_slow(pft,1) = litter_ag_slow(pft,1) + nind(pft) *(sm_ind(pft,1) + hm_ind(pft,1))
    litter_bg(pft,1)      = litter_bg(pft,1)      + nind(pft) * rm_ind(pft,1)

    nind(pft) = 0.

    lm_ind(pft,:) = 0.
    sm_ind(pft,:) = 0.
    hm_ind(pft,:) = 0.
    rm_ind(pft,:) = 0.

  else if (.not.present(pft) .and. survive(pft) .and. estab(pft) .and. aprec > aprec_min_estab) then

    !Introduce PFT if conditions suitable for establishment

    present(pft) = .true.

    if (tree(pft)) then
      nind(pft) = 0.
    else
      nind(pft) = 1.    !each grass PFT = 1 "individual"
    end if

    fpc_grid(pft) = 0.

    if (.not.tree(pft)) crownarea(pft) = 1.

    leafon(pft)      =.true.
    leafondays(pft)  = 0
    leafoffdays(pft) = 0

  end if

end do

!sapling and grass establishment

!Calculate total woody FPC and number of woody PFTs present and able to establish
 
ngrass         = count(present .and. .not.tree)
npft_estab     = count(present .and. estab .and. tree)
fpc_total      = sum(fpc_grid)
fpc_tree_total = sum(fpc_grid,mask = tree)

acflux_estab = 0.

! Rajouter maximum clay content for establishment 
! Creer un vecteur d'establishment ! E.C. 21.10.16

if (aprec >= aprec_min_estab .and. npft_estab > 0) then

  !Calculate establishment rate over available space, per tree PFT
  !Maximum establishment rate reduced by shading as tree FPC approaches 1
  !Total establishment rate partitioned equally among regenerating woody PFTs

  estab_rate = estab_max * (1. - exp(5. * (fpc_tree_total - 1.))) / real(npft_estab) 

  !Calculate grid-level establishment rate per woody PFT
  !Space available for woody PFT establishment is proportion of grid cell
  !not currently occupied by woody PFTs and not affected by fire in this year

  estab_grid = max(estab_rate * (1. - fpc_tree_total),0.)
  
  !write(0,*)'estab rate',estab_rate,fpc_tree_total, burnedf,estab_grid,npft_estab

else  !unsuitable climate for establishment

  estab_grid = 0.

end if

do pft = 1,npft
  
  if (present(pft) .and. estab(pft)) then
    if (tree(pft)) then
      if (estab_grid > 0.) then
	  
	  ! Ici on va limiter letablissement que si le pourcentage de clay pour les 4 PFTs est inferieur a un seuil 
	  ! et si il y a pas eu de feux alors le pin ne peut pas setablir
	  
! 		if ((pft == 1. .and. clay_mean < 20.0) .OR. (pft == 3. .and. clay_mean < 13.0) .OR. (pft == 4. .and. clay_mean < 18.0) .OR. (pft == 8. .and. clay_mean < 23.0)) then 
! 		
! 			if ((pft == 1. .and. burnedf >= 0.) .OR. (pft == 3. .and. burnedf >= 0.) .OR. (pft == 4. .and. burnedf > 0.) .OR. (pft == 8. .and. burnedf >= 0.)) then
			
				crownarea_max = pftpar(pft,18)

		  !      write(0,*)'estab_grid',pft,estab_grid,npft_estab,crownarea_max

				!Add new saplings to current population

				nind0 = nind(pft)
				nind(pft) = nind0 + estab_grid   

				sm_ind_tmp    = (sm_ind(pft,1) * nind0 + sm_sapl(pft,1) * estab_grid) / nind(pft) ! Sapwood mass
				hm_ind(pft,1) = (hm_ind(pft,1) * nind0 + hm_sapl(pft,1) * estab_grid) / nind(pft) ! Heartwood mass
				lm_ind(pft,1) = (lm_ind(pft,1) * nind0 + lm_sapl(pft,1) * estab_grid) / nind(pft) ! Leaf mass
				rm_ind(pft,1) = (rm_ind(pft,1) * nind0 + rm_sapl(pft,1) * estab_grid) / nind(pft) ! Root mass
				
				!Accumulate biomass increment due to sapling establishment

				estab_mass = lm_sapl(pft,1) + sm_sapl(pft,1) + hm_sapl(pft,1) + rm_sapl(pft,1)

				estab_pft(pft) = estab_mass * estab_grid

				if (estab_mass * estab_grid > eps) acflux_estab(1) = acflux_estab(1) + estab_mass * estab_grid

				stemdiam = (4. * (sm_ind_tmp + hm_ind(pft,1)) / wooddens / pi / allom2)**(1./(2. + allom3)) !Eqn 9

				height(pft) = allom2 * stemdiam**allom3                           !Eqn C

				crownarea(pft) = min(crownarea_max,allom1 * stemdiam**reinickerp) !Eqn D

				!Recalculate sapwood mass, transferring excess sapwood to heartwood compartment, if necessary to satisfy Eqn A

				sm_ind(pft,1) = lm_ind(pft,1) * height(pft) * wooddens * sla(pft) / latosa 

				hm_ind(pft,1) = max(hm_ind(pft,1) + (sm_ind_tmp - sm_ind(pft,1)),0.)

				!if (pft == 5) write(0,*)'establishment',estab_grid,crownarea(pft),sm_ind(pft,1),sm_ind_tmp,sm_ind_tmp - sm_ind(pft,1)
			  
! 			 end if 
! 		
! 		end if 
		
		end if  !estab grid
 
    else !grass

      !Grasses can establish in non-vegetated areas

      bare = (1. - fpc_total) / real(ngrass)

      lm_ind(pft,1) = lm_ind(pft,1) + bare * lm_sapl(pft,1)
      rm_ind(pft,1) = rm_ind(pft,1) + bare * rm_sapl(pft,1)

      !Accumulate biomass increment due to grass establishment

      estab_mass    = bare * (lm_sapl(pft,1) + rm_sapl(pft,1))

      estab_pft(pft) = estab_mass * crownarea(pft)

      if (estab_mass * crownarea(pft) > eps) acflux_estab(1) = acflux_estab(1) + estab_mass * crownarea(pft)

    end if !tree grass

  end if  !present and estab true

  if (present(pft)) then

    !Update LAI and FPC

    if (crownarea(pft) > 0.) then
       lai_ind(pft) = (lm_ind(pft,1) * sla(pft)) / crownarea(pft)
    else
       lai_ind(pft) = 0.
    end if

    fpc_ind       = (1. - exp(-0.5 * lai_ind(pft)))
    fpc_grid(pft) = crownarea(pft) * nind(pft) * fpc_ind
  
  end if
  
end do  !PFT loop

end subroutine establishment

!---------------------------------------------------------

end module establishmentmod
