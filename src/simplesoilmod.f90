module simplesoilmod

implicit none

public :: simplesoil

contains

!---------------------------------------------------------------------------

subroutine simplesoil(soil,soilpar)

use parametersmod,   only : sp
use mpistatevarsmod, only : soildata
use pedotransfermod, only : fbulk,calctheta,fKsat,ombd,omcf

implicit none

type(soildata), intent(inout) :: soil  !state variables sent back out with MPI
real(sp), dimension(:), intent(out) :: soilpar

integer, parameter :: nl = size(soil%sand)
integer :: l
integer :: it

integer :: soiltype

real(sp) :: sand
real(sp) :: clay
real(sp) :: OM    !(mass %)

real(sp) :: silt
real(sp) :: bulk
real(sp) :: Tsat
real(sp) :: T33
real(sp) :: T1500

real(sp), dimension(nl) :: dz
real(sp), dimension(nl) :: zpos
real(sp), dimension(nl) :: OrgM  !(g m-2)
real(sp), dimension(nl) :: whc
real(sp), dimension(nl) :: Ksat

real(sp) :: Csoil     !(g m-2)
real(sp) :: orgC      !(mass %)
real(sp) :: soilmass  !(g m-3)
real(sp) :: blk0      !(g m-3)
real(sp) :: dOM       !change in organic matter (g m-2)
real(sp) :: dzOM      !interlayer transport of SOM (g m-2)
real(sp) :: dzx       !excess change in top layer thickness

logical, dimension(nl) :: valid

!----------

valid = .true.

zpos  = soil%zpos
dz = 0.2

OrgM  = 0.

do l = 1,nl
  sand  = soil%sand(l)
  clay  = soil%clay(l)
  OM    = soil%orgm(l)
  
  silt = 100. - (sand + clay)
  
  if (sand < 0.) then
    valid(l) = .false.
    cycle
  end if

  Csoil = OrgM(l) / omcf  !layer
  
  orgC = Csoil
  
  !write(0,*)'lyr ',l,' tile',i
  !write(0,*)'dz  ',dz(l)
  !write(0,*)'Csol',Csoil

  !because bulk density depends strongly on organic matter content and
  !weakly on wilting point water content, we guess an initial value and
  !iterate to a stable solution
  
  T1500 = 0.1
  
  bulk = fbulk(orgC,T1500*100.,clay,zpos(l),silt)
  
  it = 1

  do
    
    !convert SOM to mass fraction and calculate the difference in OM content

    soilmass = bulk * 1.e6 * dz(l)     !(g) (m-2)  note: g cm-3 * (m*100=cm) 100cm * 100cm = g

    orgC = max(100. * Csoil / soilmass,0.)  !(mass %)

    !recalculate bulk

    blk0 = fbulk(orgC,T1500*100.,clay,zpos(l),silt)  !units (g cm-3)

    !calculate wilting point, field capacity, and saturation, needs input in fractions not percent

    OM = orgC * omcf             !(mass %)

    if (OM >= 30.) soiltype = 3  !humic soil

    call calctheta(sand/100.,clay/100.,OM/100.,blk0,Tsat,T33,T1500,soiltype)

    !recalculate bulk

    bulk = fbulk(orgC,T1500*100.,clay,zpos(l),silt)  !units (g cm-3)

    if (abs(bulk - blk0) < 0.001 .or. it > 50) exit

    blk0 = bulk
    
    it = it + 1

  end do

  !with the final value for bulk density, recalculate porosity

  call calctheta(sand/100.,clay/100.,OM/100.,bulk,Tsat,T33,T1500)

  !update layer-integrated WHC
      
  whc(l) = 1000. * dz(l) * (T33 - T1500)  !mm

  !calculate saturated conductivity

  Ksat(l) = fKsat(Tsat,T33,T1500)

  !write(0,*)'layer',l
  !write(0,*)'sand',sand
  !write(0,*)'silt',silt
  !write(0,*)'clay',clay
  !write(0,*)'OM  ',OM
  !write(0,*)'bulk',bulk

  !write(0,*)'Tsat',Tsat
  !write(0,*)'T33 ',T33
  !write(0,*)'Twp ',T1500
  !write(0,*)'Ksat',Ksat(l)
  
  soil%bulk(l) = bulk
  
end do  !layers

where (.not. valid)  !assign a typical rock value for non-soil layers
  whc  = 0.03
  Ksat = 1.e-3
  valid = .true.
end where

!saturated conductivity (mm h-1)
soilpar(1) = Ksat(1)
soilpar(2) = sum(Ksat(2:nl))/count(valid(2:nl))

!water holding capacity (mm)
soilpar(3) = whc(1)
soilpar(4) = sum(whc(2:nl),mask=valid(2:nl)) + whc(nl) / dz(nl) * 2. !mm/m layer add 2 more meters like LPJ2

soilpar(5) = 0.2    !thermal diffusivity (mm2/s) at wilting point (0% WHC)
soilpar(6) = 0.650  !thermal diffusivity (mm2/s) at 15% WHC
soilpar(7) = 0.4    !thermal diffusivity at field capacity (100% WHC)

!  write(0,*)'input soilpars'
!  write(0,*)'Ksat',soilpar(1:2)
!  write(0,*)'whc ',soilpar(3:4)

end subroutine simplesoil 

!---------------------------------------------------------------------------

end module simplesoilmod
