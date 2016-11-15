module pedotransfermod

implicit none

public  :: fbulk
public  :: calctheta
public  :: fKsat
private :: fTsat

!type
integer, parameter :: sp = selected_real_kind(4)

!parameters
real(sp), parameter :: ombd = 0.224  !bulk density of organic matter (g cm-3)
real(sp), parameter :: omcf = 1.724  !conversion factor from organic carbon to organic matter

contains

!---------------------------------------------------------------------------

subroutine calctheta(sand,clay,OM,bulk,Tsat,T33,T1500,soiltype)

!equations from JAP Pollacco, 2008. Can. J. Soil Sci. 88: 761-774

implicit none

!arguments

real(sp), intent(in)  :: sand               !mass fraction of sand (g g-1) (0-1)
real(sp), intent(in)  :: clay               !mass fraction of clay (g g-1) (0-1)
real(sp), intent(in)  :: OM                 !mass fraction of organic matter (g g-1) (0-1) (multiply by 1.724 if input data is in terms of Carbon)
real(sp), intent(in)  :: bulk               !total soil bulk density (g cm-3)
integer,  intent(in), optional :: soiltype  !soil type code 1=typical, 2=tropical, 3=humic, 4=vitric

real(sp), intent(out) :: Tsat               !volumetric saturated water content (m3 m-3)
real(sp), intent(out) :: T33                !volumetric water content at field capacity (Psi = -33 kPa)   (m3 m-3)
real(sp), intent(out) :: T1500              !volumetric water content at wilting point  (Psi = -1500 kPa) (m3 m-3)

!parameters

type fitpars
  real(sp), dimension(2) :: Pmax
  real(sp), dimension(2) :: Pmin
  real(sp), dimension(2) :: Psand
  real(sp), dimension(2) :: Pclay
  real(sp), dimension(2) :: Pp
end type fitpars

type(fitpars) :: typical
type(fitpars) :: tropical
type(fitpars) :: humic
type(fitpars) :: vitric

type(fitpars) :: pars

!local variables

real(sp) :: sand_om   !sand fraction corrected for organic matter (0-1) (sand_om = (1-OM)*sand)
real(sp) :: clay_om

real(sp), dimension(2) :: Pmin
real(sp), dimension(2) :: Pmax
real(sp), dimension(2) :: Pclay
real(sp), dimension(2) :: Psand
real(sp), dimension(2) :: Pp

real(sp) :: Wsat    !saturated gravimetric soil moisture content (g g-1)

real(sp), dimension(2) :: W  !gravimetric soil moisture content (g g-1)

!---------------------------------------------------------------------------
!optimal fitting parameters for W33 and W1500, from Table 7 in Pollacco 2008

!                    FC     WP
typical%Pmax   = [ 0.953, 0.895 ]
typical%Pmin   = [ 0.608, 0.165 ]
typical%Psand  = [ 0.215, 0.000 ]
typical%Pclay  = [ 0.914, 0.759 ]
typical%Pp     = [-0.102, 1.468 ]

!                    FC     WP
tropical%Pmax  = [ 1.000, 0.891 ]
tropical%Pmin  = [ 0.781, 0.197 ]
tropical%Psand = [ 0.338, 0.000 ]
tropical%Pclay = [ 2.104, 0.521 ]
tropical%Pp    = [-2.009, 0.767 ]

!                    FC     WP   
humic%Pmax     = [ 0.685, 0.551 ]
humic%Pmin     = [ 0.654, 0.190 ]
humic%Psand    = [ 0.217, 0.000 ]
humic%Pclay    = [ 3.010, 0.372 ]
humic%Pp       = [-1.810, 0.273 ]

!                    FC     WP   
vitric%Pmax    = [ 1.000, 1.000 ]
vitric%Pmin    = [ 0.371, 0.094 ]
vitric%Psand   = [ 0.187, 0.000 ]
vitric%Pclay   = [ 0.563, 0.757 ]
vitric%Pp      = [-0.030, 0.616 ]

!---------------------------------------------------------------------------

if (present(soiltype)) then
  select case(soiltype)
  case default
    pars = typical
  case(2)
    pars = tropical
  case(3)
    pars = humic
  case(4)
    pars = vitric
  end select
else
  pars = typical
end if

Pmin  = pars%Pmin
Pmax  = pars%Pmax
Psand = pars%Psand
Pclay = pars%Pclay
Pp    = pars%Pp

!-----

sand_om = (1. - OM) * sand
clay_om = (1. - OM) * clay

Tsat = fTsat(bulk,OM,sand_om)

if (bulk > 0.) then

  Wsat = Tsat / bulk

  !Pollacco PTF Model 4, eqn. 7a
  W = Wsat * (Pmin + (Pmax - Pmin) * clay_om**(Pclay + Pp * Wsat**2)) * exp(-(Psand * sand_om**3) / Wsat)
  
else
  
  Wsat = 0.
  W    = 0.

end if

T33   = W(1) * bulk   !eqn. 1, NB pw (water density assumed = 1)
T1500 = W(2) * bulk

end subroutine calctheta

!---------------------------------------------------------------------------

function fTsat(bulk,OM,sand_om)

!equations from JAP Pollacco, 2008. Can. J. Soil Sci. 88: 761-774

implicit none

real(sp) :: fTsat

real(sp), intent(in) :: bulk     !bulk density (g cm-3)
real(sp), intent(in) :: OM       !organic matter mass fraction (0-1)
real(sp), intent(in) :: sand_om  !sand fraction corrected for organic matter (0-1) (sand_om = (1-OM)*sand)

real(sp) :: pmin   !mineral density (g cm-3)
real(sp) :: pom    !particle density of organic matter (g cm-3)
real(sp) :: pp     !soil particle density (g cm-3)

!-----

pmin  = 2.69 - 0.03 * sand_om  !pg 763

pom   = 1.127 + 0.373 * OM     !pg. 763

pp    = 1. / ((OM / pom) + ((1. - OM) / pmin))  !eqn. 3

fTsat = 1. - bulk / pp   !eqn. 2a

end function fTsat

!---------------------------------------------------------------------------

function fbulk(orgC,W15,clay,depth,silt)

!equations from Heuscher et al., (2005) Soil Sci. Soc. Am. J. 69

real(sp) :: fbulk  !bulk density (g cm-3)

real(sp), intent(in) :: orgC   !mass %
real(sp), intent(in) :: W15    !water content at -15 bar (-1500 kPa) (mass %)
real(sp), intent(in) :: clay   !mass %
real(sp), intent(in) :: silt   !mass %
real(sp), intent(in) :: depth  !horizon depth (cm)

real(sp), parameter :: incp  =  1.685    !intercept coefficient
real(sp), parameter :: pOC   = -0.198    !organic carbon coefficient 
real(sp), parameter :: pwc   = -0.0133   !wilting point water content coefficient
real(sp), parameter :: pclay =  0.0079   !clay coef
real(sp), parameter :: pdep  =  0.00014  !depth coef
real(sp), parameter :: psilt = -0.0007   !silt coef

fbulk = incp + pOC * sqrt(orgC) + pwc * W15 + pclay * clay + pdep * depth + psilt * silt

end function fbulk

!---------------------------------------------------------------------------

function fKsat(Tsat,T33,T1500)

!equations from Saxton & Rawls (2006) Soil Sci. Soc. Am. J. 70

real(sp) :: fKsat  !(mm h-1)

real(sp), intent(in) :: Tsat
real(sp), intent(in) :: T33
real(sp), intent(in) :: T1500

real(sp), parameter :: l1500 = log(1500.)
real(sp), parameter :: l33   = log(33.)
real(sp), parameter :: num   = l1500 - l33

real(sp) :: B
real(sp) :: lambda

B = num / (log(T33) - log(T1500))

lambda = 1. / B

fKsat = 1930. * (Tsat - T33)**(3.-lambda)

end function fKsat

!---------------------------------------------------------------------------

end module pedotransfermod
