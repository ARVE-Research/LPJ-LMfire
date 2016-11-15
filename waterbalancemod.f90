module waterbalancemod

implicit none

public  :: waterbalance
private :: soilwater

contains

!----------------------------------------------------------------------------------------------------------------------

!subroutine waterbalance
!Calculations of soil hydrology, plant water balance, daily canopy
!conductance and ipft daily water stress factor

subroutine waterbalance(d,present,rootprop,w,dgp,dpet,dphen,  &
                        dgc,dmelt,dprec,ksat,awc,drunoff_drain,  &
                        drunoff_surf,dwscal,daet,fpc_grid,mat20,idx)   

use parametersmod, only : sp,npft,i8

implicit none 

integer(i8) :: idx

!parameters

real(sp), parameter :: emax   = 5.   !maximum daily transpiration rate (mm/day)
real(sp), parameter :: alpham = 1.4  !maximum Priestly-Taylor coefficient
real(sp), parameter :: gm     = 5.   !scaling conductance (mm s-1)?

!arguments

integer,  intent(in) :: d
real(sp), intent(in) :: mat20  !20 year running mean annual temperature

real(sp), intent(in), dimension(:) :: dpet
real(sp), intent(in), dimension(:) :: dmelt
real(sp), intent(in), dimension(:) :: dprec
real(sp), intent(in), dimension(:) :: ksat  !(mm hr-1)
real(sp), intent(in), dimension(:) :: awc
real(sp), intent(in), dimension(:) :: fpc_grid
logical,  intent(in), dimension(:) :: present

real(sp), intent(in), dimension(:,:) :: rootprop
real(sp), intent(in), dimension(:,:) :: dphen  !phenological state (0=no leaves, 1=full canopy)
real(sp), intent(in), dimension(:,:) :: dgp    !non water stressed canopy conductance (?)

real(sp), intent(inout), dimension(:) :: w

real(sp), intent(out) :: daet
real(sp), intent(out) :: drunoff_surf
real(sp), intent(out) :: drunoff_drain

real(sp), intent(out), dimension(:) :: dwscal
real(sp), intent(out), dimension(:,:) :: dgc    !actual canopy conductance (?)


!local variables

integer :: ipft

real(sp) :: supply
real(sp) :: demand
real(sp) :: demandpot
real(sp) :: wr

real(sp), dimension(2) :: aettotal
real(sp), dimension(2) :: beta

real(sp), dimension(npft) :: aet

!-----------------------------------------

daet = 0.
aettotal = 0.
drunoff_surf = 0.
drunoff_drain = 0.

do ipft = 1,npft

  !Calculate actual canopy conductance, potential water scalar and actual evapotranspiration (aet) for each ipft

  if (present(ipft)) then

    !Calculate effective supply function in root zone today
    !Eqn 24, Haxeltine & Prentice 1996

    wr = rootprop(1,ipft) * w(1) + rootprop(2,ipft) * w(2)

    supply = emax * wr

    !Calculate actual evapotranspiration demand function and potential demand assuming full leaf cover
    !Eqn 23, Haxeltine & Prentice 1996

    if (dphen(d,ipft) > 0.) then 

      demand = dpet(d) * alpham * (1. - exp(-dgp(d,ipft) * dphen(d,ipft) / gm))
      
      !alternative formulation used by Gerten et al. 2004 (results in slightly lower demand for any given dgp)
      !demand = (1. - wet) * dpet(d) * alpham / (1. + gm / (dgp(d,pft) * dphen(d,pft))) 
      
    else
      demand = 0.
    end if

    if (dphen(d,ipft) == 1.) then
      demandpot = demand 
    else
      demandpot = dpet(d) * alpham * (1. - exp(-dgp(d,ipft) / gm))
    end if

    !Calculate daily potential water scalar                        

    if (demandpot > 1.e-10) then                 !FLAG changed here from demand to demandpot to avoid flucutating leafout situation
      dwscal(ipft) = min(supply/demandpot,1.) 
    else 
      dwscal(ipft) = 1.
    end if
    
    !--------------------------
    !if (ipft == 9) then
    !  write(0,'(a,4f9.4)')'waterbalance dwscal out',supply,demandpot,supply/demandpot,dwscal(ipft)
    !end if
    !--------------------------

    !Calculate actual canopy conductance (gc), according to balance between supply and demand
    !Eqn 25, Haxeltine & Prentice 1996

    if (supply >= demand) then
      dgc(d,ipft) = dgp(d,ipft) * dphen(d,ipft) 
    else if (dpet(d) > 0.) then
      dgc(d,ipft) = max(-gm * log(1. - supply / (dpet(d) * alpham)),0.)
    else
      dgc(d,ipft) = 0.
    end if

    !AET is smaller of supply and demand
    !Eqn 22, Haxeltine & Prentice 1996

    aet(ipft) = min(supply,demand) 

    daet = daet + aet(ipft)
    
    !if (idx == 1) then
    !  write(*,*)ipft,aet(ipft),supply,demand
    !end if
    

    !Accumulate total AET
    !Eqns 28,29, Haxeltine & Prentice 1996

    if (wr == 0.) then
      beta(1) = 0.
      beta(2) = 0.
    else 
      beta(1) = rootprop(1,ipft) * w(1) / wr 
      beta(2) = rootprop(2,ipft) * w(2) / wr 
    end if

    aettotal(1) = aettotal(1) + beta(1) * aet(ipft) 
    aettotal(2) = aettotal(2) + beta(2) * aet(ipft)
    
  end if

end do

!put this statement in to approximate some evaporation from the soil surface, even if there are no plants
!fpc_grid is used to estimate the bare ground fraction

aettotal(1) = aettotal(1) + w(1) * dpet(d) * (1. - sum(fpc_grid))

call soilwater(awc,ksat,dmelt(d),dprec(d),aettotal,w,drunoff_surf,drunoff_drain,mat20,idx)

!write(90,*)year,d,dpet(d),sum(aettotal),dprec(d)

end subroutine waterbalance                                                                        

!----------------------------------------------------------------------------------------------------------------------

subroutine soilwater(awc,ksat,melt,prec,aettotal,swf,runoff,drainage,mat20,idx)

use parametersmod, only : sp,i8

implicit none

integer(i8) :: idx

!arguments

real(sp), intent(in),    dimension(2) :: awc       !soil available water content (mm or kg)
real(sp), intent(in),    dimension(2) :: ksat      !saturated conductivity (mm h-1)

real(sp), intent(in)                  :: melt      !daily total snowmelt (mm)
real(sp), intent(in)                  :: prec      !daily total precipitation (mm)

real(sp), intent(in),    dimension(2) :: aettotal  !total water removed through evapotranspiration from each soil layer (mm)

real(sp), intent(inout), dimension(2) :: swf       !water status in each soil layer (fraction of awc)

real(sp), intent(out)                 :: runoff    !surface runoff (mm)
real(sp), intent(out)                 :: drainage  !groundwater drainage (mm)

real(sp), intent(in) :: mat20  !20 year running mean annual temperature

!parameters

integer,  parameter :: ts = 24
real(sp), parameter :: dt = 1. / real(ts)
real(sp), parameter :: a  = 1.e-5  !minimum unsaturated conductivity (mm hr-1)

!local variables

integer  :: t         !counters
real(sp) :: qsurf     !surface flux of water (mm h-1)
real(sp) :: perc      !percolation of soil water (mm) between layers
real(sp) :: reversef  !reverse flow: if lower layer is saturated, water backs up into the upper layer

real(sp) :: infil


real(sp), dimension(2) :: whcmm   !water holding capacity of the soil in each layer (mm)
real(sp), dimension(2) :: kmmh    !saturated conductivity of the soil in each layer (mm hr-1)
real(sp), dimension(2) :: swater  !soil water content (mm)
real(sp), dimension(2) :: aet     !evapotranspiration (mm h-1)

!------------------------------

infil    = 0.
perc     = 0.
drainage = 0.
runoff   = 0.
reversef = 0.

whcmm = awc   !(/  84.56579590,  300.7981567   /) 
kmmh  = ksat  !(/   7.46948719,    5.689514160 /)

!instantaneous hydraulic conductivity is based on the power4 scaling to wetness fraction (soil at FC behaves as saturated flow)

qsurf = dt * (melt + prec)   !mm hr-1
aet   = dt * aettotal

do t = 1,ts  !hourly loop
  
  kmmh = ksat * swf**4 - a*swf**4 + a  !unsaturated conductivity (mm hr-1)
  
  swater = swf * whcmm  !total water mass (mm)
  
  !drainage out the bottom
    
  if (mat20 <= 0.) then    !very simple parameterization of permafrost impedance of drainage
    drainage = 0.
  else
    drainage = min(swater(2),kmmh(2))
  end if
  
  !percolation between layers

  perc = min(swater(1),kmmh(1))
  
  !infiltration is a function of the available pore volume and the conductivity
  
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

  !update swf

  where (whcmm > 0.)
    swf = swater / whcmm
  elsewhere
    swf = 0.
  end where

end do

!if (idx == 1) then
!  write(*,*)melt,prec,runoff,drainage,aettotal
!end if


if (swf(2) <= 0.) then

  !write(0,'(5f10.4)')prec,melt,qsurf,infil,runoff
  !write(0,'(5f10.4)')whcmm(1),kmmh(1),swater(1),swf(1),perc
  !write(0,'(6f10.4)')whcmm(2),kmmh(2),swater(2),swf(2),drainage,reversef

!  stop
end if

end subroutine soilwater

!----------------------------------------------------------------------------------------------------------------------

end module waterbalancemod

