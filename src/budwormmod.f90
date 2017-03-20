module budwormmod

use parametersmod,    only : sp,dp,npft

implicit none

public  :: RateDevelopment
public  :: DevelopmentStatus

contains

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine RateDevelopment(year,i,j,d,met,osv,rate)

use parametersmod,   only : pft
use weathergenmod,   only : metvars_out
use mpistatevarsmod, only : statevars
use iovariablesmod, only : budwormParamfile

implicit none

integer, intent(in) :: year  !year number, for diagnostic output, to be removed later
integer, intent(in) :: i  !tile index
integer, intent(in) :: j  !gridcell index
integer, intent(in) :: d  !day of the year
type(metvars_out), intent(in) :: met
type(statevars), target, intent(inout) :: osv  !state variables; contains biomass
real(sp), dimension(12), intent(out):: rate				! Variable 1 pour calcul du taux de developpement


!pointers

logical,  pointer, dimension(:) :: present         !PFT is present

real(sp), pointer, dimension(:) :: nind            !gridcell individual density (indiv/m2)
real(sp), pointer, dimension(:) :: lm_ind          !leaf mass, average individual (g C)
 
! Variables temperature 
real(sp) :: tmin     ! Temperature journaliere minimum
real(sp) :: tmax     ! Temperature journaliere minimum
real(sp) :: tmean    ! Temperature journaliere minimum

!     LOCAL VARIABLES
!real budwormpar(1:nparam_budworm,1:nstage)	! Pour enregistrement dans une table des parametres du taux developpement
!real table(1:nparam_budworm,1:nstage)     	! Pour enregistrement dans une table des parametres du taux developpement
integer :: a, b, g, s, genre					! Variables des boucles 
real, dimension(12) :: tau				! Variable 2 pour calcul du taux de developpement
integer nstage,nsexe,ngroupe,nparam_budworm,ierror
        parameter (nstage=12,nsexe=2,ngroupe=5,nparam_budworm=6)
 
! Parameters for development 
real(sp), dimension(12) :: b1 = [ 0.194, 0.910, 0.430, 1.210, 0.269,  0.28,  0.317, 0.259, 0.205,  57.8, 0.228, 0.277]
real(sp), dimension(12) :: b2 = [   3.0,  2.91,  3.06,   3.8,  3.02,  2.67,   3.06,  2.75,  2.85, -3.08,  3.12, 32.14]
real(sp), dimension(12) :: b3 = [  5.84,  5.32,  6.85,  7.55,  8.57,  5.03,   4.66,  4.66,  6.28, 0.045,  5.94, 11.63]
real(sp), dimension(12) :: b4 = [ 0.034, 0.061, 0.061, 0.148, 0.005, 0.151,  0.136, 0.053, 0.044,   0.0, 0.073,   0.0]
real(sp), dimension(12) :: Tb = [   2.5,   4.4,   4.4,   4.4,   4.4,   4.4,    4.4,   4.4,   4.4,   8.0,   6.0,   6.2]
real(sp), dimension(12) :: Tm = [  35.0,  38.0,  38.0,  38.0,  38.0,  38.0,   38.0,  35.0,  35.0,  35.0,  35.0,  40.0]

!--------------------------------------------------------------------
!Enregistrement des temperetures MIN & MAX
tmax = met%tmax
tmin = met%tmin
tmean = tmin + (tmax - tmin)

! Read table of development parameters for each stages ! 
!call getarg(5,budwormParamfile)  
	  
!open (unit = 16, file=budwormParamfile, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
	  
!if (ierror .ne. 0) then
!	write(*,*) 'File budworm parameters cannot be open'
!else 
!    read (16, *)
!    do a = 1, nparam_budworm
!		read (16, *) table(a,:)
!		do b=1,nstage
!			budwormpar(a,b)=table(a,b)
!		end do  
!    end do
!end if
!close (unit = 16)
!write(*,*) 'jai correctement lu les parametres du budworm module'
	
!--------------------------------------------------------------------------
! CALCULATION RATE DEVELOPMENT BY STAGE AND DEPENDANT TO TEMPERATURE

do s = 1,nstage ! L2o, L2, L3, L4, L5, L6, Pupa, Adult, Egg, L1

	! (1) Calculate rate development: depend of temperature and stage
	! 	  Regniere et al. (2012) : Equation 8 for L1, Equation 9 for Adult, and the others Equation 7
	
	! L1
	if (s == 12. .and. tmean >= Tb(s)) then  
		tau(s) = 0.
		rate(s) = b1(s) * exp (-0.5 * (((tmean-b2(s)) / b3(s))**2)) 
	
	! Adult 
	else if (s == 10.) then 
		tau(s) = min(max(tmean,Tb(s)),Tm(s))
		rate(s) = 1 / (b1(s) + (b2(s) * tau(s)) + (b3(s) * tau(s)**2)) 
	
	! Egg, L20 to L5
	else if (s == 1. .or. s == 2. .or. s == 3. .or. s == 4. .or. s == 5. .or. s == 6. .or. s == 7. .or. s == 8. .or. s == 9. .or. s == 11.) then
		if(tmean >= Tb(s) .and. tmean <= Tm(s)) then 
			tau(s) = (tmean - Tb(s)) / (Tm(s) - Tb(s))
			rate(s) = b1(s) * ((1/ (1 + exp((b2(s) - (b3(s) * tau(s)))))) - (exp (((tau(s) - 1)/b4(s)))))
		else 
			rate(s) = 0.
		end if
	else 
	    tau(s) = 0.  
        rate(s) = 0.
	end if 
	
end do ! End stage

select case(d)
case(150)
write(0,*)'tau',year,d,tmin,tmax,tmean,tau
write(0,*)'rate',year,d,tmin,tmax,tmean,rate 
end select

end subroutine RateDevelopment

!-----------------------------------------------------------------------------------------------------------------------------------------------------
subroutine DevelopmentStatus(year,i,j,d,rate)

use parametersmod,   only : pft
use mpistatevarsmod, only : statevars

implicit none

integer, intent(in) :: year  !year number, for diagnostic output, to be removed later
integer, intent(in) :: i  !tile index
integer, intent(in) :: j  !gridcell index
integer, intent(in) :: d  !day of the year
real(sp), dimension(12), intent(in):: rate				! Variable 1 pour calcul du taux de developpement

!pointers
integer nstage,nsexe,ngroupe,nparam_budworm,ierror
        parameter (nstage=12,nsexe=2,ngroupe=5,nparam_budworm=6)
real(sp), parameter :: seuil   =  0.15    		   ! Minimum rate of development between groups
integer :: g, s						           	   ! Variables des boucles 
real state(1:5,1:12)				   			   ! Pour enregistrement dans une table des parametres du taux developpement

!-----------------------------------------------------------------
!Remember: Situation for year 1 and others are different because year 1 we need to create population of sprucebudworm

do s = 1,nstage ! Egg, L1, L2o, L2, L3, L4, L5, L6, Pupa, Adult 

	! We work with five groups because we want to put variability in our model: all individus in one stage dont develop at the same time
	do g = 1,ngroupe 
	
		if (d == 1.) then 
			state(g,s) = rate(s)
		
		!-------------------------
		! L2o
		!-------------------------
		
		! Groupe 1
		else if (d > 1. .and. g == 1. .and. s == 1.) then 
			state(g,s) = state(g,s) + rate(s)  ! Ici on ajoute l'ancien state avec le rate 
		
		! Autres Groupes / 2 cas en fonction de si le seuil est atteint ou non 
		else if (d > 1. .and. g > 1. .and. s == 1. .and. state(g-1,s) < seuil) then 
			state(g,s) = 0.
		else if (d > 1. .and. g > 1. .and. s == 1. .and. state(g-1,s) >= seuil) then 
			state(g,s) = state(g,s) + rate(s)
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HERE
		
		!-------------------------
		! L2 to L5, L6 male, Adult, Egg, L1 (car depend de la colonne d'avant)
		!------------------------- 
		  
		! 3 cas possibles: 	(A) si le state du stade d'avant n'est pas 1 alors pas de developpement au stade d'après 
		!					(B) si state au state d'avant >= 1 && que c'est le premier pas de temps de dvmpt pour ce stade alors dvmpt mais pondération en fonction du temps de passage: 
		!							On détermine le pourcentage restant de la journée de 24h qu'il a pour évolué dans le stade d'après, 
        ! 							afin de ne pas prendre 100% de ce développement car c'est faux puisqu'il y a une partie de la journée ou il était au stade d'avant. 
		! 					(C) si même chose que (2) mais que c'est pas le premier pas de temps de dvmpt alors calcul du dvmpt normal 
		
		else if (s == 2. .or. s == 3. .or. s == 4. .or. s == 5. .or. s == 6. .or. s == 10. .or. s == 11. .or. s == 12.) then 
		 
			if (d > 1. .and. state(g,d,s-1) < 1.) then 
				state(g,d,s) = 0.
			else if (d > 1. .and. state(g,d,s-1) >= 1. .and. state(g,d-1,s) == 0.) then 
			    state(g,d,s) = ((state(g,d-1,s) + rate(s)) * (((state(g,d,s-1)-1)*24)/(state(g,d,s-1)-state(g,d-1,s-1))))/100
			else if (d > 1. .and. state(g,d,s-1) >= 1. .and. state(g,d-1,s) > 0.) then 
			    state(g,d,s) = state(g,d-1,s) + rate(s)
			end if 
			
		!-------------------------
		! L6 female and Pupa (male and female) : ici different car on ne fait pas reference a la colonne juste avant mais deux colonnes avant car on a les sexes qui sont pris en compte
		!------------------------- 	
		
		! Female du L6 and Pupa
		else if (s == 7. .or. s == 8. .or. s == 9.) then 
		
			if (d > 1. .and. state(g,d,s-2) < 1.) then 
				state(g,d,s) = 0. 	
			else if (d > 1. .and. state(g,d,s-2) >= 1. .and. state(g,d-1,s) == 0.) then 
			    state(g,d,s) = ((state(g,d-1,s) + rate(s)) * (((state(g,d,s-2)-1)*24)/(state(g,d,s-2)-state(g,d-1,s-2))))/100
			else if (d > 1. .and. state(g,d,s-2) >= 1. .and. state(g,d-1,s) > 0.) then 
			    state(g,d,s) = state(g,d-1,s) + rate(s)
			end if
		
		else 

			state(g,d,s) = 0.
			
		end if ! Fin des conditions pour le calcul du state
		
	end do ! Boucle groupe
	
end do ! Boucle stage

end subroutine DevelopmentStatus

!-----------------------------------------------------------------------------------------------------------------------------------------------------
end module budwormmod
