module simplesoilmod

implicit none

public :: simplesoil

contains

! ---------------------------------------------------------------------------

subroutine simplesoil(soil,soilpar,soilprop)

use parametersmod,   only : sp,stdout
use mpistatevarsmod, only : soildata
use pedotransfermod, only : fbulk,calctheta,fKsat,fBexp,ombd,omcf

implicit none

type(soildata),           intent(inout) :: soil  ! state variables sent back out with MPI
real(sp), dimension(:),   intent(out)   :: soilpar
real(sp), dimension(:,:), intent(out)   :: soilprop

integer :: nl
integer :: l
integer :: it

integer :: soiltype

real(sp) :: sand
real(sp) :: clay
real(sp) :: OM    ! organic matter (mass %)
real(sp) :: OMf   ! organic matter (mass fraction)

real(sp) :: silt
real(sp) :: bulk
real(sp) :: Tsat
real(sp) :: T33
real(sp) :: T1500

real(sp) :: Bexp_min
real(sp) :: Bexp_org
real(sp) :: Bexp

real(sp) :: Psat_min
real(sp) :: Psat_org
real(sp) :: Psat

real(sp), allocatable, dimension(:) :: dz
real(sp), allocatable, dimension(:) :: zpos
real(sp), allocatable, dimension(:) :: OrgM  ! (g m-2)
real(sp), allocatable, dimension(:) :: whc
real(sp), allocatable, dimension(:) :: Ksat

real(sp) :: Csoil     ! (g m-2)
real(sp) :: orgC      ! organic carbon (mass %)
real(sp) :: soilmass  ! (g m-3)
real(sp) :: blk0      ! (g m-3)
real(sp) :: dOM       ! change in organic matter (g m-2)
real(sp) :: dzOM      ! interlayer transport of SOM (g m-2)
real(sp) :: dzx       ! excess change in top layer thickness

logical, allocatable, dimension(:) :: valid

! ----------

nl = size(soil%sand)

! write(0,*)'SOIL NLAYERS',nl

allocate(dz(nl))
allocate(zpos(nl))
allocate(OrgM(nl))
allocate(whc(nl))
allocate(Ksat(nl))
allocate(valid(nl))

valid = .true.

zpos  = soil%zpos
dz = 0.2

OrgM  = 0.

do l = 1,nl
  sand  = soil%sand(l)
  clay  = soil%clay(l)
  OM    = soil%orgm(l)
  
  silt = 100. - (sand + clay)
  
!   write(0,*)'SOIL',l,sand,silt,clay,OM
  
  if (sand < 0.) then
    valid(l) = .false.
    cycle
  end if

  Csoil = OrgM(l) / omcf  ! layer
  
  orgC = Csoil
  
  ! write(stdout,*)'lyr ',l,' tile',i
  ! write(stdout,*)'dz  ',dz(l)
  ! write(stdout,*)'Csol',Csoil

  ! because bulk density depends strongly on organic matter content and
  ! weakly on wilting point water content, we guess an initial value and
  ! iterate to a stable solution
  
  T1500 = 0.1
  
  bulk = fbulk(orgC,T1500*100.,clay,zpos(l),silt)
  
  it = 1

  do
    
    ! convert SOM to mass fraction and calculate the difference in OM content

    soilmass = bulk * 1.e6 * dz(l)     ! (g) (m-2)  note: g cm-3 * (m*100=cm) 100cm * 100cm = g

    orgC = max(100. * Csoil / soilmass,0.)  ! (organic carbon, mass %)

    ! recalculate bulk

    blk0 = fbulk(orgC,T1500*100.,clay,zpos(l),silt)  ! units (g cm-3)

    ! calculate wilting point, field capacity, and saturation, needs input in fractions not percent

    OM = orgC * omcf             ! organic matter (mass %)
    
    OMf = OM / 100.              ! organic matter (mass fraction)

    if (OM >= 30.) soiltype = 3  ! humic soil

    call calctheta(sand/100.,clay/100.,OM/100.,blk0,Tsat,T33,T1500,soiltype)

    ! recalculate bulk

    bulk = fbulk(orgC,T1500*100.,clay,zpos(l),silt)  ! units (g cm-3)

    if (abs(bulk - blk0) < 0.001 .or. it > 50) exit

    blk0 = bulk
    
    it = it + 1

  end do

  ! with the final value for bulk density, recalculate porosity

  call calctheta(sand/100.,clay/100.,OM/100.,bulk,Tsat,T33,T1500)

  ! update layer-integrated WHC
      
  whc(l) = 1000. * dz(l) * (T33 - T1500)  ! mm

  ! calculate saturated conductivity

  Ksat(l) = fKsat(Tsat,T33,T1500)
  
  ! calculate B exponent for water retention curve (Psi-Theta relationship)
  
  ! Bexp = fBexp(clay,sand,0.,OM)  ! Bloemen formulation, setting rock fraction to zero temporarily
  
  Bexp_min = 2.91 + 0.159 * clay     ! CLM4 eqn. 7.84
  
  Bexp_org = 2.7
  
  Bexp = (1. - OMf) * Bexp_min + OMf * Bexp_org
  
  ! calculate Psi_sat (soil matric potential at saturation)
  
  Psat_min = -10. * 10**(1.88 - 0.0131 * sand)  ! CLM4 eqn. 7.87
  
  Psat_org = -10.3
  
  Psat = (1. - OMf) * Psat_min + OMf * Psat_org  ! CLM4 eqn. 7.86 (mm)
    
  ! save derived soil properties
  
  soilprop(l,1) = Tsat
  soilprop(l,2) = T33
  soilprop(l,3) = T1500
  soilprop(l,4) = Psat
  soilprop(l,5) = Bexp
  
!   write(stdout,*)'layer',l
!   write(stdout,*)'sand',sand
!   write(stdout,*)'silt',silt
!   write(stdout,*)'clay',clay
!   write(stdout,*)'OM  ',OM
!   write(stdout,*)'bulk',bulk

!   write(stdout,*)'Tsat',Tsat
!   write(stdout,*)'T33 ',T33
!   write(stdout,*)'Twp ',T1500
!   write(stdout,*)'Ksat',Ksat(l)
!   write(stdout,*)'Bexp',Bexp
!   write(stdout,*)'Psat',Psat
  
  soil%bulk(l) = bulk
  
end do  ! layers

where (.not. valid)  ! assign a typical rock value for non-soil layers
  whc  = 0.03
  Ksat = 1.e-3
  valid = .true.
end where

! saturated conductivity (mm h-1)
soilpar(1) = Ksat(1)
soilpar(2) = sum(Ksat(2:nl))/count(valid(2:nl))

! water holding capacity (mm)
soilpar(3) = whc(1)
soilpar(4) = sum(whc(2:nl),mask=valid(2:nl)) + whc(nl) / dz(nl) * 2. ! mm/m layer add 2 more meters like LPJ2

soilpar(5) = 0.2    ! thermal diffusivity (mm2/s) at wilting point (0% WHC)
soilpar(6) = 0.650  ! thermal diffusivity (mm2/s) at 15% WHC
soilpar(7) = 0.4    ! thermal diffusivity at field capacity (100% WHC)

!  write(stdout,*)'input soilpars'
!  write(stdout,*)'Ksat',soilpar(1:2)
!  write(stdout,*)'whc ',soilpar(3:4)

end subroutine simplesoil 

! ---------------------------------------------------------------------------

end module simplesoilmod
