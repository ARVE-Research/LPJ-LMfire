module waterbalancemod

implicit none

public  :: waterbalance
private :: soilwater
private :: hydraulics

contains

! ----------------------------------------------------------------------------------------------------------------------

! subroutine waterbalance
! Calculations of soil hydrology, plant water balance, daily canopy
! conductance and pft daily water stress factor

subroutine waterbalance(d,present,rootprop,w,dgp,dpet,dphen,      &
                        dgc,dmelt,dprec,ksat,awc,drunoff_drain,   &
                        drunoff_surf,dwscal,daet,fpc_grid,mat20,  &
                        dtemp,soilprop,lm_ind,sm_ind,rm_ind,fpc_ind,sla,latosa,tree)   

use parametersmod, only : sp,npft,i8,pliq,grav

implicit none 

! parameters

real(sp), parameter :: emax   = 5.   ! maximum daily transpiration rate (mm/day)
real(sp), parameter :: alpham = 1.4  ! maximum Priestly-Taylor coefficient
real(sp), parameter :: gm     = 5.   ! scaling conductance (mm s-1)?

! arguments

integer,  intent(in) :: d
real(sp), intent(in) :: mat20  ! 20 year running mean annual temperature

real(sp), dimension(:), intent(in) :: dpet
real(sp), dimension(:), intent(in) :: dmelt
real(sp), dimension(:), intent(in) :: dprec   ! daily total precipitation (mm)
real(sp), dimension(:), intent(in) :: ksat    ! soil saturated hydraulic conductivity (mm hr-1)
real(sp), dimension(:), intent(in) :: awc     ! soil water holding capacity (field capacity - wilting point) (mm)
real(sp), dimension(:), intent(in) :: fpc_grid
logical,  dimension(:), intent(in) :: present

real(sp), dimension(:,:), intent(in) :: rootprop
real(sp), dimension(:,:), intent(in) :: dphen  ! phenological state (0=no leaves, 1=full canopy)
real(sp), dimension(:,:), intent(in) :: dgp    ! non water stressed potential canopy conductance (mm d-1)

real(sp), dimension(:), intent(inout) :: w

real(sp), intent(out) :: daet
real(sp), intent(out) :: drunoff_surf
real(sp), intent(out) :: drunoff_drain

real(sp), intent(out), dimension(:) :: dwscal
real(sp), intent(out), dimension(:,:) :: dgc    ! actual canopy conductance (mm d-1)

! passthrough arguments for hydraulics

real(sp), dimension(:),   intent(in) :: dtemp
real(sp), dimension(:,:), intent(in) :: soilprop
real(sp), dimension(:),   intent(in) :: lm_ind
real(sp), dimension(:),   intent(in) :: sm_ind
real(sp), dimension(:),   intent(in) :: rm_ind
real(sp), dimension(:),   intent(in) :: fpc_ind
real(sp), dimension(:),   intent(in) :: sla
logical,  dimension(:),   intent(in) :: tree

real, intent(in) :: latosa

! local variables

integer :: pft

real(sp) :: supply
real(sp) :: demand
real(sp) :: demandpot
real(sp) :: wr

real(sp), dimension(2) :: aettotal
real(sp), dimension(2) :: beta

real(sp), dimension(npft) :: aet

! -----------------------------------------

daet = 0.
aettotal = 0.
drunoff_surf = 0.
drunoff_drain = 0.

do pft = 1,npft

  ! Calculate actual canopy conductance, potential water scalar and actual evapotranspiration (aet) for each pft

  if (present(pft)) then

    ! 2023.04 need to calculate demand function before supply, in order to determine leaf water potential     

    ! Calculate actual evapotranspiration demand function and potential demand assuming full leaf cover
    ! Eqn 23, Haxeltine & Prentice 1996

    if (dphen(d,pft) > 0.) then 

      demand = dpet(d) * alpham * (1. - exp(-dgp(d,pft) * dphen(d,pft) / gm))
      
      ! alternative formulation used by Gerten et al. 2004 (results in slightly lower demand for any given dgp)
      ! demand = (1. - wet) * dpet(d) * alpham / (1. + gm / (dgp(d,pft) * dphen(d,pft))) 
      
    else
      demand = 0.
    end if

    if (dphen(d,pft) == 1.) then
      demandpot = demand 
    else
      demandpot = dpet(d) * alpham * (1. - exp(-dgp(d,pft) / gm))
    end if

    ! --------------------------
    ! root-zone weighted soil moisture

    wr = rootprop(1,pft) * w(1) + rootprop(2,pft) * w(2)

    ! maximum potential water supply based on root-zone soil moisture and a maximum transpiration rate
    ! Eqn 24, Haxeltine & Prentice 1996

    supply = emax * wr  ! supply at maximum transpiration rate
    
    ! --------------------------
    ! if max supply is less than demand, then calculate actual supply based on hydraulic properties

    if (supply < demand .and. tree(pft)) then
    
      call hydraulics(w,rootprop(:,pft),soilprop,dtemp(d),pft,tree(pft),   &
                      lm_ind(pft),rm_ind(pft),sm_ind(pft),fpc_ind(pft),sla(pft),latosa,supply)
    
    end if
    
    ! --------------------------

    ! Calculate actual canopy conductance (gc), according to balance between supply and demand
    ! Eqn 25, Haxeltine & Prentice 1996

    if (supply >= demand) then

      dgc(d,pft) = dgp(d,pft) * dphen(d,pft) 

    else if (dpet(d) > 0.) then

      dgc(d,pft) = max(-gm * log(1. - supply / (dpet(d) * alpham)),0.)

    else

      dgc(d,pft) = 0.

    end if

    ! AET is smaller of supply and demand
    ! Eqn 22, Haxeltine & Prentice 1996

    aet(pft) = min(supply,demand) 

    daet = daet + aet(pft)
    
    ! Accumulate total AET
    ! Eqns 28,29, Haxeltine & Prentice 1996

    if (wr == 0.) then
      beta(1) = 0.
      beta(2) = 0.
    else 
      beta(1) = rootprop(1,pft) * w(1) / wr 
      beta(2) = rootprop(2,pft) * w(2) / wr 
    end if

    aettotal(1) = aettotal(1) + beta(1) * aet(pft) 
    aettotal(2) = aettotal(2) + beta(2) * aet(pft)

    ! --------------------------
    ! Calculate daily potential water scalar, used in calculating the raingreen phenology                      

    if (demandpot > 1.e-10) then                 ! FLAG changed here from demand to demandpot to avoid flucutating leafout situation
      dwscal(pft) = min(supply / demandpot,1.) 
    else 
      dwscal(pft) = 1.
    end if

  end if ! present
end do  ! pft loop

! put this statement in to approximate some evaporation from the soil surface, even if there are no plants
! fpc_grid is used to estimate the bare ground fraction

aettotal(1) = aettotal(1) + w(1) * dpet(d) * (1. - sum(fpc_grid))

call soilwater(awc,ksat,dmelt(d),dprec(d),aettotal,w,drunoff_surf,drunoff_drain,mat20)

! write(90,*)year,d,dpet(d),sum(aettotal),dprec(d)

end subroutine waterbalance                                                                        

! ----------------------------------------------------------------------------------------------------------------------

subroutine soilwater(awc,ksat,melt,prec,aettotal,swf,runoff,drainage,mat20)

use parametersmod, only : sp,i8

implicit none

! NB whc is defined as the difference in theta between theta(Psi-33) and theta(Psi-1500)

! arguments

real(sp), intent(in),    dimension(2) :: awc       ! soil available water content (mm or kg)
real(sp), intent(in),    dimension(2) :: ksat      ! saturated conductivity (mm h-1)

real(sp), intent(in)                  :: melt      ! daily total snowmelt (mm)
real(sp), intent(in)                  :: prec      ! daily total precipitation (mm)

real(sp), intent(in),    dimension(2) :: aettotal  ! total water removed through evapotranspiration from each soil layer (mm)

real(sp), intent(inout), dimension(2) :: swf       ! water status in each soil layer (fraction of awc)

real(sp), intent(out)                 :: runoff    ! surface runoff (mm)
real(sp), intent(out)                 :: drainage  ! groundwater drainage (mm)

real(sp), intent(in) :: mat20  ! 20 year running mean annual temperature

! parameters

integer,  parameter :: ts = 24
real(sp), parameter :: dt = 1. / real(ts)
real(sp), parameter :: a  = 1.e-5  ! minimum unsaturated conductivity (mm hr-1)

! local variables

integer  :: t         ! counters
real(sp) :: qsurf     ! surface flux of water (mm h-1)
real(sp) :: perc      ! percolation of soil water (mm) between layers
real(sp) :: reversef  ! reverse flow: if lower layer is saturated, water backs up into the upper layer

real(sp) :: infil

real(sp), dimension(2) :: whcmm   ! water holding capacity of the soil in each layer (mm)
real(sp), dimension(2) :: kmmh    ! saturated conductivity of the soil in each layer (mm hr-1)
real(sp), dimension(2) :: swater  ! soil water content (mm)
real(sp), dimension(2) :: aet     ! evapotranspiration (mm h-1)

! ------------------------------

infil    = 0.
perc     = 0.
drainage = 0.
runoff   = 0.
reversef = 0.

whcmm = awc   ! (/  84.56579590,  300.7981567   /) 
kmmh  = ksat  ! (/   7.46948719,    5.689514160 /)

! instantaneous hydraulic conductivity is based on the power4 scaling to wetness fraction (soil at FC behaves as saturated flow)

qsurf = dt * (melt + prec)   ! mm hr-1
aet   = dt * aettotal

do t = 1,ts  ! hourly loop

  where (swf < 1.e-5) swf = 0.
  
  kmmh = ksat * swf**4 - a * swf**4 + a  ! unsaturated conductivity (mm hr-1)
  
  swater = swf * whcmm  ! total water mass (mm)
  
  ! drainage out the bottom
    
  if (mat20 <= 0.) then    ! very simple parameterization of permafrost impedance of drainage
    drainage = 0.
  else
    drainage = min(swater(2),kmmh(2))
  end if
  
  ! percolation between layers

  perc = min(swater(1),kmmh(1))
  
  ! infiltration is a function of the available pore volume and the conductivity
  
  infil = qsurf

  swater(1) = max(swater(1) + infil - aet(1) - perc, 0.)
  swater(2) = max(swater(2) + perc  - aet(2) - drainage, 0.)

  if (swater(2) > whcmm(2)) then
    reversef  = swater(2) - whcmm(2)
    swater(1) = swater(1) + reversef
    swater(2) = whcmm(2)
  end if

  if (swater(1) > whcmm(1)) then
    runoff = runoff + swater(1) - whcmm(1)
    swater(1) = whcmm(1)
  end if

  ! update swf

  where (whcmm > 0.)
    swf = swater / whcmm
  elsewhere
    swf = 0.
  end where

end do

if (swf(2) <= 0.) then

  ! write(stdout,'(5f10.4)')prec,melt,qsurf,infil,runoff
  ! write(stdout,'(5f10.4)')whcmm(1),kmmh(1),swater(1),swf(1),perc
  ! write(stdout,'(6f10.4)')whcmm(2),kmmh(2),swater(2),swf(2),drainage,reversef

!  stop
end if

end subroutine soilwater

! ----------------------------------------------------------------------------------------------------------------------

subroutine hydraulics(w,rootprop,soilprop,dtemp,pft,tree,lm_ind,rm_ind,sm_ind,fpc_ind,sla,latosa,supply)

! Hickler et al. (2006) hydraulics scheme for calculating supply

use parametersmod, only : sp,wooddens,grav,pliq

implicit none

! arguments

real(sp), dimension(:),   intent(in)  :: w         ! soil moisture as a proportion of field capacity (per layer, fraction)
real(sp), dimension(:,:), intent(in)  :: soilprop  ! per layer and property
real(sp), dimension(:),   intent(in)  :: rootprop  ! per layer
real(sp),                 intent(in)  :: dtemp     ! mean daily temperature (degC)
integer,                  intent(in)  :: pft
logical,                  intent(in)  :: tree
real(sp),                 intent(in)  :: lm_ind
real(sp),                 intent(in)  :: sm_ind
real(sp),                 intent(in)  :: rm_ind
real(sp),                 intent(in)  :: fpc_ind
real(sp),                 intent(in)  :: sla
real(sp),                 intent(in)  :: latosa    ! leaf area to sapwood area ratio (m2 m-2)
real(sp),                 intent(out) :: supply

! parameters

real(sp), parameter :: a1 = 0.556   ! water viscosity scaling parameters
real(sp), parameter :: a2 = 0.022

real(sp), parameter :: Kl = 1.5 ! leaf specific conductivity per unit leaf area (10^-7 m s-1 MPa-1)

! minimum leaf water potential (MPa)
real(sp), dimension(9), parameter :: Plmin = [ -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2 ]

! maximum sapwood specific conductivity  (10^-4 m s-1 MPa-1)
real(sp), dimension(9), parameter :: Ks_max = [ 50., 10., 8., 5., 30., 8., 20., -9999., -9999.]

! shape parameter for sapwood specific conductivity curve
real(sp), dimension(9), parameter :: c = [ 5., 3., 3., 3., 3., 3., 3., -9999., -9999. ]

! soil water potential that causes a 50% lost of sapwood conductance
real(sp), parameter, dimension(9) :: P50 = [ -0.3, -1.3, -2., -1.7, -1., -2., -1., -9999., -9999. ]

real(sp), parameter :: Kr = 4.  ! fine root specific conductivity (10^-7 m3 kg-1 s-1 MPa-1) only for trees

real(sp), parameter :: mm2MPa = grav / 1.e6   ! convert pressure in mm H2O to MPa

! local variables

integer  :: l        ! counter
integer  :: nl       ! number of soil layers

real(sp) :: height   ! plant height (m)
real(sp) :: sap_xsa  ! sapwood cross-sectional area (m2)
real(sp) :: dP       ! water potential gradient (MPa)
real(sp) :: Ks       ! sapwood specific conductivity (10^-4 m s-1 MPa-1) calculated dynamically
real(sp) :: Psat     ! texture-dependent soil water potential at saturation (from CLM, mm)
real(sp) :: Bexp     ! Brooks-Corey B exponent
real(sp) :: Pleaf    ! leaf water potential (MPa)
real(sp) :: Psoil    ! soil water potential (MPa)
real(sp) :: Tf       ! temperature scaling factor (dimensionless)

! hydraulic resistances (MPa)
real(sp) :: Rr       ! fine roots
real(sp) :: Rs       ! stem
real(sp) :: Rl       ! leaves
real(sp) :: R        ! total hydraulic resistance MPa mm-1

real(sp) :: Tsat     ! total soil porosity
real(sp) :: T33      ! soil porosity at field capacity
real(sp) :: T1500    ! soil porosity at wilting point
real(sp) :: Theta    ! soil porosity at current water content

real(sp) :: Trz      ! root-zone weighted Theta / Tsat 
real(sp) :: awc      ! available water holding capacity (fraction) Theta(Psi-33) - Theta(Psi-1500) 

! ------------------

nl = size(w)

Trz = 0.

do l = 1,nl
  
  Tsat  = soilprop(l,1)
  T33   = soilprop(l,2)
  T1500 = soilprop(l,3)
  
  awc = T33 - T1500
  
  Theta = awc * w(l) + T1500
  
  Trz = Trz + Theta / Tsat * rootprop(l)
  
end do

! tree allometry

sap_xsa = lm_ind * sla / latosa * 1.e4  ! unit conversion from m2 to cm2 because Ksmax requires cm2

height = sm_ind / sap_xsa / wooddens

! write(0,*)'hydraulics',lm_ind,sm_ind,sap_xsa,height

! root zone weighted average soil properties

Psat = rootprop(1) * soilprop(1,4) + rootprop(2) * soilprop(2,4)  

Bexp = rootprop(1) * soilprop(1,5) + rootprop(2) * soilprop(2,5)  

! soil water potential

Psoil = Psat * Trz**(-Bexp) * mm2MPa    ! this equation returns pressure in MPa

! leaf water potential
! since we are in a supply-limited situation, leaf water potential is set to a typical minumum value

Pleaf = Plmin(pft) ! (MPa)

! water potential gradient

dP = Pleaf - Psoil + (height * pliq * grav) * 1.e-6 ! scale to result in MPa

! compartment resistances

! write(0,*)'hydraulics',pft,Psoil,Psat,wr,Bexp

Ks = Ks_max(pft) * exp(-(Psoil / P50(pft))**c(pft))

! write(0,*)'Ks',Ks_max(pft),exp(-(Psoil / P50(pft))**c(pft)),sap_xsa,1. / (Ks * sap_xsa)

if (Ks < 1.e-6) then  ! the soil is too dry and there is no sapwood conductivity
  
  supply = 0.
  return

end if

Rr = 1. / (Kr * rm_ind / 1000.)  ! root mass needs to be in units Kg m-2
Rs = 1. / (Ks * sap_xsa)
Rl = 1. / (Kl * fpc_ind)         ! the rationale for scaling leaf resistance by fpc is not really clear

! write(0,*)'resistance',rm_ind,Rr,fpc_ind,Rl

! temperature scalar for increase in viscosity at low temps

Tf = a1 + a2 * dtemp ! could separate into air and soil components and scale temperature differentially

! total path resistance

R = (Rr + Rs + Rl) / Tf

! supply

supply = max(-dP / R,0.)   ! Hickler et al. eqn 10

! write(0,*)'hydraulics ',w,Trz,Psoil,Pleaf,(height * pliq * grav) * 1.e-6,dP,Rr,Rs,Rl,supply
! 10 format (a,3f8.4,4f12.6,4f12.2)


! write(0,*)'hydraulics',pft,dP,R,supply,height

end subroutine hydraulics

! ----------------------------------------------------------------------------------------------------------------------

end module waterbalancemod

