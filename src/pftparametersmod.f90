module pftparametersmod

use parametersmod, only : sp

implicit none

public :: pftparameters

! ----------------------------

contains

subroutine pftparameters(pftparsfile,pftpar,sla,tree,evergreen,summergreen,raingreen,needle,boreal,        &
                         lm_sapl,sm_sapl,hm_sapl,rm_sapl,latosa,allom1,allom2,allom3,allom4,wooddens,co2)

! Assignment of PFT-specific parameters and bioclimatic limits
! Definition of initial sapling and grass mass structure
! Calculation of maximum crown area for woody PFTs

use parametersmod, only : sp,npft,npftpar,ncvar,pi,reinickerp
   
implicit none

! parameters

! arguments

character(*) :: pftparsfile

real(sp), dimension(npft,npftpar), intent(out) :: pftpar

real(sp), intent(in) :: latosa
real(sp), intent(in) :: wooddens
real(sp), intent(in) :: allom1
real(sp), intent(in) :: allom2
real(sp), intent(in) :: allom3
real(sp), intent(in) :: allom4

real(sp), dimension(npft,ncvar), intent(out) :: lm_sapl
real(sp), dimension(npft,ncvar), intent(out) :: sm_sapl
real(sp), dimension(npft,ncvar), intent(out) :: rm_sapl
real(sp), dimension(npft,ncvar), intent(out) :: hm_sapl

real(sp), dimension(ncvar), intent(in) :: co2

real(sp), dimension(npft), intent(out) :: sla

logical, dimension(npft), intent(out) :: tree
logical, dimension(npft), intent(out) :: evergreen
logical, dimension(npft), intent(out) :: summergreen
logical, dimension(npft), intent(out) :: raingreen
logical, dimension(npft), intent(out) :: needle
logical, dimension(npft), intent(out) :: boreal

! parameters

real(sp), parameter :: sla1 = 2.e-4 * exp(6.15)

! local variables

integer :: ierror
integer :: n
integer :: pft

real(sp) :: lai_sapl       ! sapling or initial grass LAI
real(sp) :: x              ! sapwood:heartwood diameter ratio
real(sp) :: lmtorm         ! non-water-stressed leafmass to rootmass ratio
real(sp) :: stemdiam       ! sapling stem diameter (m)
real(sp) :: height_sapl    ! sapling height (m)
real(sp) :: crownarea      ! sapling crown area (m2)
real(sp) :: longevity      ! leaf longevity (months)

! -----------------------------------------------------------------------------
! List of the PFT parameters
!
! ==== ecophysiological parameters ====
!
!  1  fraction of roots in upper soil layer                             ! used in gppmod
!  2  plants with C4 (1) or C3 (0) photosynthetic pathway               ! used in gppmod & pftparameter pour 13C valeurs initiales pour les jeunes arbres et l'herbe
!  3  water scalar value at which leaves shed by drought deciduous PFT  ! used in gppmod, les feuillus perdent leurs feuilles si la valeur d'eau scalaire est inférieure a celle-ci
!  4  canopy conductance component (gmin, mm/s) not associated with     ! used in gppmod, minimum canopy conductance réduit par la fraction de couverture projective du feuillage en supposant la couverture complète de la feuille photosynthesis (Haxeltine & Prentice 1996, Table 4)
!  5  maintenance respiration coefficient                               ! used in nppmod, coefficient de respiration
!  6  maximum foliar N content (mg/g)                                   ! used in gppmod (Haxeltine & Prentice 1996a, Fig 4)
!  7  leaf longevity (years)                                            ! used in gppmod et pftparameter pour calculer le specific leaf area
!  8  leaf turnover period (years)                                      ! used in turnovermod sachant que les taux de rotation sont inverses de la longévité des tissus
!  9  sapwood turnover period (sapwood converted to heartwood) (years)  ! même chose que précédent
! 10  root turnover period (years)                                      ! même chose que précédent
! 11  leaf C:N mass ratio                                               ! used in nppmod pour calculer la NPP
! 12  sapwood C:N mass ratio                                            ! même chose que précédent
! 13  root C:N mass ratio                                               ! même chose que précédent
! 14  leaf type: broadleaved (1), needleleaved (2) or grass (3)
! 15  phenology type: evergreen (1), summergreen (2), raingreen (3), any type (4) 
! 16  leaf to root ratio under non-water stressed conditions            ! used in pftparameter: Calculate sapling or initial grass rootmass et est égale a (leafmass) / (rootmass)
!                                                                       ! et call dans allocationmod: prend le max entre ce ratio et la valeur water scalar per pft on this day. Valeur minimale fixée à 10%, les valeurs inférieures à cette ne sont généralement pas trouvé dans la nature
! 17  summergreen phenology ramp, GDD5 requirement to grow full leaf canopy
! 18  tree maximum crown area (m2)
! 19  sapling (or grass on initialisation) LAI
! 20  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
! 21  boreal pft (1), non-boreal pft (0)     
! 22  low temperature limit for CO2 uptake
! 23  lower range of temperature optimum for photosynthesis
! 24  upper range of temperature optimum for photosynthesis
! 25  high temperature limit for CO2 unptake
! 26  optimal Ci/Ca ratio
!
! ==== bioclimatic limits ====
! 
! 27 minimum coldest monthly mean temperature
! 28 maximum coldest monthly mean temperature
! 29 minimum growing degree days (at or above 5 deg C)
! 30 upper limit of temperature of the warmest month 
! 31 lower limit of growth efficiency (g/m2) REMARQUE: Ce parametre n'est pas pris en compte
!
! ==== individual parameters ====
!
! 32 crown length
! 33 bark thickness 1
! 34 bark thickness 2
! 35 height slope (hmax_r-hsap_r)/(hmax_m-hsap_m)
! 36 height intercept (hsap_r - hs * hsap_m)
! 37 diameter slope (cm diameter gained per m height)
! 38 diameter intercept
!
! ==== fire parameters ====
!
! 39 flammability threshold                          ! used in fire.f, pour calculer un facteur de pondération de l'humidité de la litière
! 40 fire resistance index                           ! used in fire.f, permet de calculer la fraction of individuals in grid cell which die
! 41 scorch height parameters = F
! 42 crown damage 1 = RCK
! 43 crown damage 2 = p
! 44 ignition efficiency parameter
! 45 emission factor for CO2
! 46 emission factor for CO
! 47 emission factor for CH4
! 48 emission factor for VOC
! 49 emission factor for TPM
! 50 emission factor for NOx from F77 code
! 51 fuel bulk density (kg m-3)
!
! -----------------------------------------------------------------------------
    
open(15,file=pftparsfile,status='old')

read(15,*)  ! header line

pft = 1

do pft = 1,npft

  read(15,*)pftpar(pft,:)
  
end do

close(15)

do pft = 1,npft
        
  ! Assign leaf and phenology logicals

  if (pftpar(pft,14) == 3.) then

    tree(pft) = .false.

  else

    tree(pft) = .true.

  end if

  if (pftpar(pft,15) <= 2.) then

    ! tree(pft)=.true.

    if (pftpar(pft,14) == 2.) then

      needle(pft) = .true.

    else

      needle(pft) = .false.

    end if

  else

    ! tree(pft) = .false.
    needle(pft) = .false.

  end if

  if (pftpar(pft,15) == 1.) then

    evergreen(pft)   = .true.
    summergreen(pft) = .false.
    raingreen(pft)   = .false.

  else if (pftpar(pft,15) == 2.) then

    evergreen(pft)   = .false.
    summergreen(pft) = .true.
    raingreen(pft)   = .false.

  else if (pftpar(pft,15) == 3.) then

    evergreen(pft)   = .false.
    summergreen(pft) = .false.
    raingreen(pft)   = .true.

  else

    evergreen(pft)   = .true.
    summergreen(pft) = .true.
    raingreen(pft)   = .true.

  end if

  if (pftpar(pft,21) == 1.) then

    boreal(pft) = .true.

  else

    boreal(pft) = .false.

  end if            

  ! Calculate specific leaf area (SLA) for each PFT from leaf longevity
  ! Include conversion (multiplier of 2.) from m2/g(dry wt) to m2/gC
  ! Equation based on Reich et al 1997, Fig 1f.
  
  longevity = 12. * pftpar(pft,7)  ! convert years to months

  ! SLA = 2e-4 * exp(6.15) / (12 * leaf_longevity)^0.46) ! Sitch et al. 2003 Eqn. 6
  
  ! the term sla1 (2.e-4 * exp(6.15)) is precalculated above
  
  sla(pft) = sla1 / longevity**0.46

  ! Define initial mass structure
  
  ! In the original LPJ, the sapling mass structure is solved by defining an initial LAI from which
  ! the sapling leaf mass is solved to (I think) assume the minimum allometrically consistent mass (leaf, sapwood, etc.)
  ! that can support the prescribed LAI. This results in saplings of highly unrealistic size and shape.
  ! As an alternative (Mar 2023), I suggest prescribing a sapling leaf mass of 25gC and calculating the
  ! rest of the allometry based on this. I have not made any change for grasses. 

  lai_sapl = pftpar(pft,19)

  if (tree(pft)) then ! woody PFTs

    ! Calculate leafmass for a sapling individual
    !  (1) lai = leafmass * sla / (crown area)
    !  (2) (leaf area) = latosa * (sapwood xs area)
    !         (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
    !  (3) (crown area) = allom1 * (stem diameter) ** reinickerp
    !         (Reinickes theory)
    ! From (1),
    !  (4) leafmass = lai * (crown area) / sla
    ! From (1) & (3),
    !  (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
    ! From (2),
    !  (6) leafmass = latosa * (sapwood xs area) / sla
    !  (7) (sapwood xs area) = pi * (sapwood diameter)**2 / 4
    ! From (6) and (7),
    !  (8) leafmass = latosa * pi * (sapwood diameter)**2 / 4 / sla
    ! From (8),
    !  (9) (sapwood diameter) = [ 4 * leafmass * sla / pi / latosa ]**0.5
    ! (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
    ! Define x,
    ! (11) x = [ (sapwood diameter)+(heartwood diameter) ] / (sapwood diameter)
    ! From (10) & (11),
    ! (12) (stem diameter) = x * (sapwood diameter)
    ! From (5), (9) & (12),
    ! (13) leafmass = lai * allom1 * x**reinickerp * (4*leafmass*sla/pi/latosa)**(reinickerp*0.5) / sla
    ! From (13),
    ! (14) leafmass = [ lai * allom1 * x**reinickerp * (4*sla/pi/latosa)**(reinickerp*0.5) / sla ]**(1-1/reinickerp)

    x = pftpar(pft,20)

!     lm_sapl(pft,1) = (lai_sapl * allom1 * x**reinickerp * (4. * sla(pft) / pi / latosa)**(reinickerp * 0.5) / & 
!                       sla(pft))**(1. - 1. / reinickerp)  !eqn 14
                      
    lm_sapl(pft,1) = 25. ! try setting to a standard 25 gC

    lm_sapl(pft,2) = 17.8 - co2(2)   ! initial 13C value from lloyd & farquhar, 1994
    lm_sapl(pft,3) = 0.

    ! Calculate sapling stem diameter
    ! From (9) & (12),
    ! (15) (stem diameter) = x * [ 4 * leafmass * sla / pi / latosa ]**0.5

    stemdiam = x * (4. * lm_sapl(pft,1) * sla(pft) / pi / latosa)**0.5  ! Eqn 15

    ! Calculate sapling height
    ! (16) height = allom2 * (stem diameter)**allom3 (source?)

    ! height_sapl = allom2 * stemdiam**allom3   !Eqn 16

    height_sapl = allom2 * stemdiam   !recent studies show that the exponent is close to 1 for saplings

    ! Calculate sapling sapwood mass
    ! (17) (sapwood volume) = height * (sapwood xs area)
    ! (18) (sapwood xs area) = leafmass * sla / latosa

    ! From (17) & (18),
    ! (19) (sapwood volume) = height * leafmass * sla / latosa
    ! (20) (sapwood mass) = (wood density) * (sapwood volume)

    ! From (19) & (20),
    ! (21) (sapwood mass) = (wood density) * height * leafmass * sla / latosa

    sm_sapl(pft,1) = wooddens * height_sapl * lm_sapl(pft,1) * sla(pft) / latosa   ! Eqn 21
    sm_sapl(pft,2) = lm_sapl(pft,2)   ! 13C value in permille
    sm_sapl(pft,3) = lm_sapl(pft,3)

    ! Calculate sapling heartwood mass
    ! From (11),
    ! (22) (heartwood mass) = (x-1) * (sapwood mass)

    hm_sapl(pft,1) = (x - 1.) * sm_sapl(pft,1)  ! Eqn 22
    hm_sapl(pft,2) = sm_sapl(pft,2)             ! 13C value in permille
    hm_sapl(pft,3) = sm_sapl(pft,3)
    
    crownarea = allom1 * stemdiam**reinickerp

    lai_sapl = lm_sapl(pft,1) * sla(pft) / crownarea

    ! write(0,*)'sapling ',pft,sla(pft),lai_sapl,lm_sapl(pft,1),height_sapl,stemdiam*100.,crownarea

  else ! grass PFT

    lm_sapl(pft,1 )= lai_sapl / sla(pft)

    ! Set initial 13C values for saplings, grass

    if (pftpar(pft,2) == 1.) then ! C4 plants

      lm_sapl(pft,2) = 3.6 - co2(2)   ! lloyd & farquhar,1994           
      lm_sapl(pft,3) = 0.

    else                          ! C3 plpants

      lm_sapl(pft,2) = 17.8 - co2(2)  ! lloyd & farquhar,1994          
      lm_sapl(pft,3) = 0.

    end if
    
    ! no sapwood and heartwood for grass PFTs

    sm_sapl(pft,2) = 0.
    hm_sapl(pft,2) = 0.
    sm_sapl(pft,3) = 0.
    hm_sapl(pft,3) = 0.

  end if

  ! Calculate sapling or initial grass rootmass
  ! (23) lmtorm = (leafmass) / (rootmass)

  lmtorm = pftpar(pft,16)
  rm_sapl(pft,1) = (1. / lmtorm) * lm_sapl(pft,1)  ! From Eqn 23
  rm_sapl(pft,2) = lm_sapl(pft,2)                  ! 13C value in permille
 

end do ! pft loop

end subroutine pftparameters

! ----------------------------

end module pftparametersmod
