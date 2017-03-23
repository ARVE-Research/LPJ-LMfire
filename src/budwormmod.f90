module budwormmod

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

!-----------------------------------------------------------------------------------------------------------------------------------------------------

use parametersmod,    only : sp,dp,npft,nistage,nisex,nirep

implicit none

public :: insect_init
public :: insect_develrate
public :: insect_growth
public :: insect_updatestate
public :: insect_layeggs

!parameters used in the different subroutines below
!development parameters
 
real(sp), dimension(nistage,nisex) :: b1
real(sp), dimension(nistage,nisex) :: b2
real(sp), dimension(nistage,nisex) :: b3
real(sp), dimension(nistage,nisex) :: b4
real(sp), dimension(nistage,nisex) :: Tb
real(sp), dimension(nistage,nisex) :: Tm

real(sp), dimension(nistage) :: survi

real(sp), dimension(nistage,nisex) :: scons
real(sp), dimension(nistage,nisex) :: indmass
real(sp), dimension(nistage,nisex) :: b2l

integer, parameter :: diagfid = 73  !unit number for diagnostic output file, can be any number but don't choose <= 10

logical :: bavard

contains

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine insect_init(dobavard)

implicit none

logical, intent(in) :: dobavard

bavard = dobavard

!subroutine to initialize parameters

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

! stage-specific consumption FEMALES (mg foliage mass / individual, for the entire stage)
!             L2o     L2     L3     L4     L5     L6f  pupaf  adult    egg     L1 
scons(:,1) = [ 0.,  1.53,  6.11, 17.22, 31.84, 271.83,    0.,    0.,     0.,   0. ]

! stage-specific consumption MALES (mg foliage mass / individual, for the entire stage)
!             L2o     L2     L3     L4     L5    L6f  pupaf  adult    egg   L1 
scons(:,2) = [ 0.,  1.71,  6.85, 11.02, 27.68, 150.42,   0.,    0.,    0.,   0. ]

! mean individual mass FEMALES (mg)
!                 L2o    L2    L3    L4    L5    L6f  pupaf  adult    egg    L1 
indmass(:,1) = [ 0.04, 0.04, 0.04,  0.3, 0.93,  3.92, 29.98, 20.48,  0.04, 0.04 ]

! mean individual mass MALES (mg)
!                 L2o    L2    L3    L4    L5  L6f    pupaf  adult    egg    L1 
indmass(:,2) = [ 0.04, 0.04, 0.04, 0.32, 0.72, 3.5,   16.26,  9.95,  0.04, 0.04 ]

! bodymass gain per leafmass consumed FEMALES (mg / mg) (empirical estimate)
!             L2o     L2     L3     L4     L5    L6f  pupaf  adult    egg   L1 
b2l(:,1) = [ 0., 0.026, 0.036, 0.037, 0.094, 0.096, 0., 0., 0., 0. ]

! bodymass gain per leafmass consumed FEMALES (mg / mg) (empirical estimate)
!             L2o     L2     L3     L4     L5    L6f  pupaf  adult    egg   L1 
b2l(:,2) = [ 0., 0.023, 0.035, 0.036, 0.1, 0.085, 0., 0., 0., 0. ]

! Survival rate modifier for attrition mode 
!                  L2o     L2     L3    L4     L5    L6    pupa  adult    egg   L1 
survi = [ 1.0,  0.79,  0.73, 0.62,   0.4, 0.66,   1.0,   0.39,  1.0,   1.0 ]

!----

if (bavard) open(diagfid,file='budworm_diagnostics.txt',status='unknown')

end subroutine insect_init

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine insect_develrate(year,i,j,d,met,osv,rate,tmean)

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
 
!--------------------------------------------------------------------
!Enregistrement des temperetures MIN & MAX

tmax  = met%tmax
tmin  = met%tmin
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

! if (bavard) then
!   write(diagfid,*)'tau',year,d,tmin,tmax,tmean,tau
!   write(diagfid,*)'rate',year,d,tmin,tmax,tmean,rate 
! end if

end subroutine insect_develrate

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine insect_updatestate(year,i,j,d,rate,state,mass,energy)

implicit none

!arguments

integer, intent(in) :: year  !year number, for diagnostic output, to be removed later
integer, intent(in) :: i     !tile index
integer, intent(in) :: j     !gridcell index
integer, intent(in) :: d     !day of the year

real(sp), dimension(:,:),   intent(inout) :: rate
real(sp), dimension(:,:,:), intent(inout) :: state
real(sp), dimension(:,:,:), intent(inout) :: mass
real(sp), dimension(:,:,:), intent(inout) :: energy

!parameters
integer, parameter  :: nisex = 2
real(sp), parameter :: seuil = 0.15    		   ! Minimum rate of development between groups

! Local variables 

integer :: r   !replicate
integer :: g   !gender
integer :: s   !stage
integer :: ns

real(sp) :: rrate !realized rate, set to 0 for L2o cohorts that are not yet developing

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
						
						if (s == 1 .and. r > 1) then
						  if (sum(state(:,g,r-1)) < 0.15 .and. state(1,g,r) == 0.) then
						    rrate = 0.
  				  else
								  rrate = rate(s,g)
						  end if
				  else
						  rrate = rate(s,g)
				  end if

						!account for increase in development stage

						state(s,g,r) = state(s,g,r) + rrate
				
						!check if development stage is completed, if yes, roll biomass up to next stage
				
						if (state(s,g,r) >= 1.) then  !stage completed

								!reset stage and roll mass up to next stage
				
								state(s,g,r) = 0.
								energy(s,g,r) = 0.
								
								if (s /= 8) then  !do in all cases except the adult --> egg transition (stage 8-9)
								  mass(ns,g,r)   = mass(s,g,r)
								  energy(ns,g,r) = energy(s,g,r)
						  end if

								mass(s,g,r)    = 0.
								energy(s,g,r)  = 0.

								!adjust rate of the next stage for this day
				
								excess = max(state(s,g,r) + rate(s,g) - 1.,0.)
								
								rate(ns,g) = excess * rate(ns,g)
				
						end if

				end do
		end do
end do

if (bavard) then
  if (any(mass > 0.)) then
    write(diagfid,*)year,d
    write(diagfid,10)'1st',mass(:,1,1)
    write(diagfid,10)'2nd',mass(:,1,2)
  end if
  
   10 format(a,10f15.4)
  
end if

end subroutine insect_updatestate

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine insect_layeggs(year,i,j,d,tmean,famass,massegg,neggi)

implicit none

!arguments

integer, intent(in) :: year  !year number, for diagnostic output, to be removed later
integer, intent(in) :: i  !tile index
integer, intent(in) :: j  !gridcell index
integer, intent(in) :: d  !day of the year

real(sp), intent(in) :: tmean
real(sp), intent(in) :: famass
real(sp), intent(inout) :: neggi

real(sp), dimension(:),   intent(inout) :: massegg  !sex

!local variables

real(sp), parameter :: unitmass = 22.882 !(mg/individual)
real(sp), parameter :: eggmass  =  0.042 !(mg/egg)

real(sp) :: nfemales
real(sp) :: neggs
real(sp) :: oogenesis

!------------------------

nfemales = famass / unitmass

neggs = nfemales * (200. - neggi)

if (tmean > 10. .and. tmean <= 25.) then

		oogenesis = neggs * (0.035 * tmean - 0.32)
		
else

		oogenesis = 0.

end if

massegg(:) = massegg(:) + 0.5 * oogenesis * eggmass

neggi = neggi + oogenesis / nfemales
  
end subroutine insect_layeggs

!-----------------------------------------------------------------------------------------------------------------------------------------------------

subroutine insect_growth(tmean,rate,state,mass,energy,leafmass)

implicit none

!arguments
real(sp), intent(in) :: tmean
real(sp), intent(in) :: leafmass
real(sp), dimension(:,:),   intent(inout) :: rate
real(sp), dimension(:,:,:), intent(inout) :: state
real(sp), dimension(:,:,:), intent(inout) :: mass
real(sp), dimension(:,:,:), intent(inout) :: energy

!parameters

! real(sp), parameter :: a0 = 2.1571
real(sp), parameter :: a1 = 8.6623e-11
real(sp), parameter :: a2 = 6.5241

real(sp), parameter :: specenergy = 5.  !cal per mg foliage mass

!local variables

integer :: r
integer :: g
integer :: s

real(sp) :: fcons
real(sp) :: eloss

real(sp), dimension(nistage,nisex,nirep) :: nindiv
real(sp), dimension(nistage,nisex,nirep) :: potcons
real(sp), dimension(nistage,nisex,nirep) :: actcons
real(sp), dimension(nistage,nisex,nirep) :: egain

!-------------------------------------------------------------------------------------------------
!temperature mortality - kill all mass in all stages except L2o when temperature is less than -10C

if (tmean < -10.) then 

  mass(2:,:,:) = 0.
  energy(2:,:,:) = 0.

  return

end if

!----
!energy loss and survival mode for low and high temperatures

if (tmean > 0.) then
  eloss = a1 * tmean**a2 !(cal)
else
  eloss = 0.
end if

if (tmean <= 2.5 .or. tmean >= 32) then

  !calculate energy loss and quit the subroutine
  
  energy = energy - eloss

  return  !no growth can happen on this day if temperatures are outside of range development and eaten foliage
  
end if

!----
!normal growth

!estimate number of individuals

do g = 1,nirep

  nindiv(:,:,g) = mass(:,:,g) / indmass

  potcons(:,:,g) = nindiv(:,:,g) * scons * rate

end do

!calculate actual consumption as a function of available leaf mass
!requires leafmass in terms of total dry mass leaf mass per m2

if (sum(potcons) > 0.) then
  fcons = min(leafmass,sum(potcons)) / sum(potcons) 
else
  fcons = 0.
end if  

actcons = fcons * potcons  !(mg m-2)

!increment insect mass and energy as a function of consumption

do g = 1,nirep
  mass(:,:,g) = mass(:,:,g) + actcons(:,:,g) * b2l
end do

where (nindiv > 0.)
  egain = actcons / nindiv * specenergy  !(mg m-2 / ind m-2 * cal mg-1 = cal ind-1)
elsewhere
  egain = 0.
end where

energy = max(energy + egain - eloss,0.)  !(cal)

!assume everything dies if the energy supply runs out

where (energy == 0.) mass = 0.

!adjust the development rate when the foliage supply does not meet demand

do s = 1,nisex
  rate(:,s) = rate(:,s) * (1. - fcons) * survi
end do

end subroutine insect_growth

!-----------------------------------------------------------------------------------------------------------------------------------------------------

end module budwormmod
