module budwormmod

use parametersmod,    only : sp,dp,npft,nistage,nisex,nirep

implicit none

public  :: RateDevelopment
public  :: DevelopmentStatus
public  :: Oviposition

contains

!-----------------------------------------------------------------------------------------------------------------------------------------------------

!!!!!!!!! SUMMARY OF THIS MODULE 
! (1) Calculate development rate function of temperature 
! (2) Calculate state development to determine transit time between two stages 
! (3) Number of individuals in each stages : NOT SURE
! (4) Oviposition depend of temperature 
! (5) Mortality by freezing < -10 degree
! (5) Loss of energy by temperature in winter : mortality overwintering 
! (6) Also mortality by attrition: depend of the survival rate 
! (6) Gain of energy for L2 to L6 if they eat foliage
! (7) Mortality by attrition if Energy = 0 after diapause, , mortality by age if they cant develop to an other stage

!s == 1  : L2o    / s == 2  : L2       / s == 3  : L3       / s == 4  : L4           / s == 5  : L5
!s == 6  : L6Male / s == 7  : L6Female / s == 8  : PupaMale / s == 9  : PupaFemale   / s == 10 : Adult
!s == 11 : Egg    / s == 12 : L1

subroutine RateDevelopment(year,i,j,d,met,osv,rate,tmean)

use parametersmod,   only : pft
use weathergenmod,   only : metvars_out
use mpistatevarsmod, only : statevars
! use iovariablesmod, only : budwormParamfile

implicit none

!arguments

integer, intent(in) :: year  !year number, for diagnostic output, to be removed later
integer, intent(in) :: i  !tile index
integer, intent(in) :: j  !gridcell index
integer, intent(in) :: d  !day of the year

type(metvars_out), intent(in) :: met
type(statevars), target, intent(inout) :: osv  !state variables; contains biomass
real(sp), intent(out):: tmean		           ! Temperature journaliere mean

real(sp), dimension(:,:), intent(out) :: rate	   ! Variable 1 pour calcul du taux de developpement
	
!pointers
! logical,  pointer, dimension(:) :: present         !PFT is present
! 
! real(sp), pointer, dimension(:) :: nind            !gridcell individual density (indiv/m2)
! real(sp), pointer, dimension(:) :: lm_ind          !leaf mass, average individual (g C)
 
! Variables temperature 
real(sp) :: tmin     ! Temperature journaliere minimum
real(sp) :: tmax     ! Temperature journaliere maximum

!     LOCAL VARIABLES
integer :: s
integer :: g

real(sp) :: tau				! Variable 2 pour calcul du taux de developpement
 
!development parameters
 
real(sp), dimension(nistage,nisex) :: b1
real(sp), dimension(nistage,nisex) :: b2
real(sp), dimension(nistage,nisex) :: b3
real(sp), dimension(nistage,nisex) :: b4
real(sp), dimension(nistage,nisex) :: Tb
real(sp), dimension(nistage,nisex) :: Tm

!--------------------------------------------------------------------
!could move this variable assignment somewhere else so it is only called once at the beginning of the run
 
! Parameters for development FEMALES     
!             L2o     L2     L3     L4     L5    L6f  pupaf  adult    egg     L1 
b1(:,1) = [ 0.194, 0.910, 0.430, 1.210, 0.269, 0.317, 0.205,  57.8, 0.228, 0.277]
b2(:,1) = [   3.0,  2.91,  3.06,   3.8,  3.02,  3.06,  2.85, -3.08,  3.12, 32.14]
b3(:,1) = [  5.84,  5.32,  6.85,  7.55,  8.57,  4.66,  6.28, 0.045,  5.94, 11.63]
b4(:,1) = [ 0.034, 0.061, 0.061, 0.148, 0.005, 0.136, 0.044,   0.0, 0.073,   0.0]
Tb(:,1) = [   2.5,   4.4,   4.4,   4.4,   4.4,   4.4,   4.4,   8.0,   6.0,   6.2]
Tm(:,1) = [  35.0,  38.0,  38.0,  38.0,  38.0,  38.0,  35.0,  35.0,  35.0,  40.0]

! Parameters for development MALES       
!             L2o     L2     L3     L4     L5    L6m  pupam  adult    egg     L1 
b1(:,2) = [ 0.194, 0.910, 0.430, 1.210, 0.269,  0.28, 0.259,  57.8, 0.228, 0.277]
b2(:,2) = [   3.0,  2.91,  3.06,   3.8,  3.02,  2.67,  2.75, -3.08,  3.12, 32.14]
b3(:,2) = [  5.84,  5.32,  6.85,  7.55,  8.57,  5.03,  4.66, 0.045,  5.94, 11.63]
b4(:,2) = [ 0.034, 0.061, 0.061, 0.148, 0.005, 0.151, 0.053,   0.0, 0.073,   0.0]
Tb(:,2) = [   2.5,   4.4,   4.4,   4.4,   4.4,   4.4,   4.4,   8.0,   6.0,   6.2]
Tm(:,2) = [  35.0,  38.0,  38.0,  38.0,  38.0,  38.0,  35.0,  35.0,  35.0,  40.0]

!--------------------------------------------------------------------
!Enregistrement des temperetures MIN & MAX
tmax = met%tmax
tmin = met%tmin
tmean = tmin + (tmax - tmin)

!--------------------------------------------------------------------------
! CALCULATION RATE DEVELOPMENT BY STAGE AND DEPENDANT TO TEMPERATURE

!initialization values (if not otherwise set)

tau  = 0.
rate = 0.

do g = 1,nisex
		do s = 1,nistage

				select case(s)
				case(10)  !L1

						if (tmean >= Tb(s,g)) rate(s,g) = b1(s,g) * exp (-0.5 * (((tmean-b2(s,g)) / b3(s,g))**2)) 

				case(8)  !Adult

						tau = min(max(tmean,Tb(s,g)),Tm(s,g))
						rate(s,g) = 1 / (b1(s,g) + (b2(s,g) * tau) + (b3(s,g) * tau**2)) 

				case default  !for all other stages

						if(tmean >= Tb(s,g) .and. tmean <= Tm(s,g)) then

								tau = (tmean - Tb(s,g)) / (Tm(s,g) - Tb(s,g))
								rate(s,g) = b1(s,g) * ((1/ (1 + exp((b2(s,g) - (b3(s,g) * tau))))) - (exp (((tau - 1)/b4(s,g)))))
				
						end if
				
				end select
		end do
end do

select case(d)
case(150)
  write(0,*)'tau',year,d,tmin,tmax,tmean,tau
  write(0,*)'rate',year,d,tmin,tmax,tmean,rate 
end select

end subroutine RateDevelopment

!-----------------------------------------------------------------------------------------------------------------------------------------------------
subroutine DevelopmentStatus(year,i,j,d,rate,state,mass)

implicit none

!arguments

integer, intent(in) :: year  !year number, for diagnostic output, to be removed later
integer, intent(in) :: i     !tile index
integer, intent(in) :: j     !gridcell index
integer, intent(in) :: d     !day of the year

real(sp), dimension(:,:),   intent(inout) :: rate
real(sp), dimension(:,:,:), intent(inout) :: state
real(sp), dimension(:,:,:), intent(inout) :: mass

!parameters
integer, parameter  :: nisex = 2
real(sp), parameter :: seuil = 0.15    		   ! Minimum rate of development between groups

! Local variables 

integer :: r   !replicate
integer :: g   !gender
integer :: s   !stage
integer :: ns

real(sp) :: excess

!-----------------------------------------------------------------
! We work with five groups because we want to put variability in our model: all individus in one stage dont develop at the same time

do r = 1,nirep
  do g = 1,nisex

    do s = 1,nistage
  
						if (mass(s,g,r) == 0.) cycle  !skip calculations if no mass in this class

      !figure out index of next stage. if at last stage, loop back to 1

      ns = s+1  
      
      if (ns > nistage) ns = 1

						!account for increase in development stage

						state(s,g,r) = state(s,g,r) + rate(s,g)
				
						!check if development stage is completed, if yes, roll biomass up to next stage
				
						if (state(s,g,r) >= 1.) then  !stage completed

								!reset stage and roll mass up to next stage
				
								state(s,g,r) = 0.
								mass(ns,g,r) = mass(s,g,r)
								mass(s,g,r)  = 0.

								!adjust rate of the next stage for this day
				
								excess = max(state(s,g,r) + rate(s,g) - 1.,0.)
								
								rate(ns,g) = excess * rate(ns,g)
				
						end if

				end do
		end do
end do

end subroutine DevelopmentStatus

!-----------------------------------------------------------------------------------------------------------------------------------------------------
subroutine Oviposition() !(year,i,j,d,presentTBE)
! 
! implicit none
! 
! integer, intent(in) :: year  !year number, for diagnostic output, to be removed later
! integer, intent(in) :: i  !tile index
! integer, intent(in) :: j  !gridcell index
! integer, intent(in) :: d  !day of the year
! 
! !pointers
! integer nstage,nisex,ngroupe,nparam_budworm,ierror
!         parameter (nstage=12,nisex=2,ngroupe=5,nparam_budworm=6)
! integer :: g, s						           	   				! Variables des boucles 
! 
! logical, dimension(1:5,1:12), intent(in):: presentTBE			
! 
! !-------------------------------
! 
! do g = 1,ngroupe 
! 
! 	if (presentTBE(g,10) == .TRUE.) then 
! 
! 	! When we have females they can lay eggs but it depends of: 
! 	! (1) their age: more females are old more they are not efficace to lay eggs
! 	! (2) Lay eggs depend of the temperature : between 10 and 25 degree celsius and not for the first 24 hours 
! 	! (3) They can lay a maximum of 200 eggs 
! 	! (4) We hypothezise that the sex ratio of their eggs is 0.5
! 
!     !(1) Determine the age of individus: Egg = 0; L1= 1; L20 = 2; L2 = 3; L3 = 4; L4 = 5; L5 = 6; L6 = 7; Pupa = 8; Adult = 9
! 	
! 	end if 
! 
! end do
! 
! 			!if (Adult == TRUE) else 
! 			
! 				! Fecondity potential of each group
! 				!potFecon = 200 * (number of Adult /2) ! Sex ratio 0.5
! 			
! 				! if the first time to lay eggs, we determine age of female depend of temperature 
! 				!if (state(g,s) > 0 .and. stateDminus1Stage >= 0 .and. tmean >=8 .and. tmean <=35) then 
! 		
! 			        !aging[(d-1),1,gr] = 1 / (57.8 - 3.08 * tmean[(d-1),3] + 0.0451 * (tmean[(d-1),3])^2) ! Maybe here do a ponderation for the remaining time of the day
! 					!age[d,1,gr] = age[(d-1),1,gr] + aging[(d-1),1,gr]
! 
! 					!# Oviposition : Equation [A6] de Régnière (83) et code de Rémi#
! 					!if (age[t,1,gr] >= 0.1 .and. tmp[(t-1),3] >= 10 .and. tmp[(t-1),3] <= 25){
! 					  
! 					!  oogenesis[t,1,gr] = (pot_fec_gr[a,t,gr] * ((0.035 * tmp[(t-1),3]) - 0.32))
! 					  
! 					!  dev[a,t,(s+1),gr,ge] = floor(oogenesis[t,1,gr]) # Pour les femelles et les males
! 
! 					!} else {  
! 					  
! 					! oogenesis[t,1,gr] = 0 
! 					!  dev[a,t,(s+1),gr,ge]  =  oogenesis[t,1,gr] 
! 
! 					!} # fin de la pondaison
! 					
! 			!end if ! Fin de la condition Adulte 
! 
end subroutine Oviposition

!-----------------------------------------------------------------------------------------------------------------------------------------------------
end module budwormmod
