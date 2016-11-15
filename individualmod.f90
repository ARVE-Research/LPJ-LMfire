module individualmod

use parametersmod,   only : sp,npft,npftpar,pftpar

implicit none

public :: allomind

type sizeind
  real(sp) :: massf
  real(sp) :: height
  real(sp) :: lcrown
  real(sp) :: stemd
  real(sp) :: barkt
end type sizeind

!-------------------------------------------------------------------

contains

subroutine allomind(pft,height,ind) !hclass,cheight,clength,stemd,barkt)

use parametersmod, only : sp,npft,npftpar,pftpar,nhclass

implicit none

!arguments

real(sp), intent(in) :: height

type(sizeind), dimension(:), intent(out) :: ind

!parameters

real(sp), parameter :: ot = 1. / 3.      !one third

real(sp), parameter :: hsapl = 1.

real(sp) :: crownarea_max    !maximum crown area (m2) 

real(sp) , dimension(npft) :: CLf   !crown length
real(sp) , dimension(npft) :: par1  !bark thickness 1
real(sp) , dimension(npft) :: par2  !bark thickness 2

!slope and intercept are the result of a comparison between LPJ-modeled maximum PFT-height and empirically observed maximum PFT-height

real(sp) , dimension(npft) :: hs   !height slope (hmax_r-hsap_r)/(hmax_m-hsap_m)
real(sp) , dimension(npft) :: hi   !height intercept (hsap_r - hs * hsap_m)

real(sp) , dimension(npft) :: ds   !diameter slope (cm diameter gained per m height)
real(sp) , dimension(npft) :: di   !diameter intercept

real(sp), parameter, dimension(nhclass) :: hsd    = [ -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0 ]

!local variables

integer :: i
integer  :: pft

real(sp) :: hbar
real(sp) :: hmax
real(sp) :: hsl

!------

!initialize height class biomass fraction (based on skewed normal distribution)

ind%massf = [ 0.0185, 0.0557, 0.1312, 0.2419, 0.3339, 0.1974, 0.0154 ]

! Enregistrement des valeurs de pfts ! E.C. 
CLf  = pftpar(:,32)
par1 = pftpar(:,33)
par2 = pftpar(:,34)
hs   = pftpar(:,35)
hi   = pftpar(:,36)
ds   = pftpar(:,37)
di   = pftpar(:,38) 
  
!calculate mean of the height class distribution based on the average individual height

hbar = hs(pft) * height + hi(pft)

!maximum individual height can be 25% taller than the current mean height of the height class distribution

hmax = 1.25 * hbar  

!calculate difference

hsl = hmax - hbar 	!well, that just equals 0.25 * hbar

!calculate height of each height class

ind%height = hsd * hsl + hbar

ind%lcrown = ind%height * CLf(pft)

!calculate fraction of total biomass in each height class
!check if one of the lower height classes is shorter than a sapling, if so, roll the distribution up to taller height classes


do i = 1,6	! we have 7 height classes, so go up to 6 to roll up to class 7 max and not step over the array boundary (M.P. 12.07.2013)
  if (ind(i)%height < hsapl) then
    
    ind(i+1)%massf = ind(i+1)%massf + ind(i)%massf
    
   ! write(0,*) i, ind(i)%height, hsapl, ind(i)%massf, ind(i+1)%massf
    
    ind(i)%massf = 0.
    
  end if
end do


!calculate stem diameter

ind%stemd = max(1.,ds(pft) * ind%height + di(pft)) !constrained to be not less than 1 cm

!calculate bark thickness

ind%barkt = par1(pft) * ind%stemd + par2(pft)                              !bark thickness (cm) eqn. 21

end subroutine allomind

end module individualmod
