!     SUBROUTINE PFTPARAMETERS
!     Assignment of PFT-specific parameters and bioclimatic limits
!     Definition of initial sapling and grass mass structure
!     Calculation of maximum crown area for woody PFTs

! Vegetation is represented by a combination of the following plant functional
! types (PFTs)

!       1. Picea
!       2. No-boreal
!       3. Abies
!       4. Pinus
!       5. No-boreal
!       6. No-boreal
!       7. No-boreal
!       8. Populus
! 	9. No-boreal

subroutine pftparameters(pftpar,sla,tree,evergreen, &
     &  summergreen,raingreen,needle,boreal,lm_sapl,sm_sapl,hm_sapl, &
     &  rm_sapl,latosa,allom1,allom2,allom3, &
     &  allom4,wooddens,co2)

use iovariablesmod, only : pftfile
	 
implicit none

!     PARAMETERS:
      integer npft,npftpar,nsoilpar,ncvar, ierror
        parameter (npft=9,npftpar=51,nsoilpar=7,ncvar=3)
      real pi
        parameter (pi=3.14159265)
      real reinickerp
        parameter (reinickerp=1.6)

!     ARGUMENTS:
      real pftpar(1:npft,1:npftpar),sla(1:npft),co2(1:3)
      logical tree(1:npft),evergreen(1:npft)
      logical summergreen(1:npft),raingreen(1:npft),needle(1:npft)
      logical boreal(1:npft)
      real lm_sapl(1:npft,1:ncvar),sm_sapl(1:npft,1:ncvar)
      real rm_sapl(1:npft,1:ncvar),hm_sapl(1:npft,1:ncvar)
      real latosa,wooddens
      real allom1,allom2,allom3,allom4


!     LOCAL VARIABLES:
      integer n,pft
      real table(1:npft,1:npftpar)
      real lai_sapl       !sapling or initial grass LAI
      real x
      real lmtorm         !non-waterstressed leafmass to rootmass ratio
      real stemdiam       !sapling stem diameter
      real height_sapl    !sapling height

!-----------------------------------------------------------------------------
 
!     PFT PARAMETERS
 
!      1  fraction of roots in upper soil layer 									! call in gppmod
!      2  plants with C4 (1) or C3 (0) photosynthetic pathway                       ! call in gppmod & pftparameter pour 13C valeurs initiales pour les jeunes arbres et l'herbe
!      3  water scalar value at which leaves shed by drought deciduous PFT 			! call in gppmod, les feuillus perdent leurs feuilles si la valeur d'eau scalaire est inférieure a celle-ci
!      4  canopy conductance component (gmin, mm/s) not associated with				! call in gppmod, minimum canopy conductance réduit par la fraction de couverture projective du feuillage en supposant la couverture complète de la feuille
!         photosynthesis (Haxeltine & Prentice 1996, Table 4)
!      5  maintenance respiration coefficient										! call in nppmod, coefficient de respiration
!      6  maximum foliar N content (mg/g)											! call in gppmod
!         (Haxeltine & Prentice 1996a, Fig 4)
!      7  leaf longevity (years) 													! call in gppmod et pftparameter pour calculer le specific leaf area
!      8  leaf turnover period (years)												! call in turnovermod sachant que les taux de rotation sont inverses de la longévité des tissus
!      9  sapwood turnover period (sapwood converted to heartwood) (years)			! même chose que précédent
!     10  root turnover period (years)												! même chose que précédent
!     11  leaf C:N mass ratio														! call in nppmod pour calculer la NPP
!     12  sapwood C:N mass ratio													! même chose que précédent
!     13  root C:N mass ratio														! même chose que précédent
!     14  leaf type: broadleaved (1), needleleaved (2) or grass (3)
!     15  phenology type: evergreen (1), summergreen (2), raingreen (3),
!         any type (4) 
!     16  leaf to root ratio under non-water stressed conditions 					! call in pftparameter: Calculate sapling or initial grass rootmass et est égale a (leafmass) / (rootmass)
!																					! et call dans allocationmod: prend le max entre ce ratio et la valeur water scalar per pft on this day. Valeur minimale fixée à 10%, les valeurs inférieures à cette ne sont généralement pas trouvé dans la nature
!     17  summergreen phenology ramp, GDD5 requirement to grow full leaf canopy
!     18  tree maximum crown area (m2)
!     19  sapling (or grass on initialisation) LAI
!     20  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
!     21  boreal pft (1), non-boreal pft (0)     
!     22  low temperature limit for CO2 uptake
!     23  lower range of temperature optimum for photosynthesis
!     24  upper range of temperature optimum for photosynthesis
!     25  high temperature limit for CO2 unptake
!	  26  optimal Ci/Ca ratio

!     BIOCLIMATIC LIMITS
 
!     27 minimum coldest monthly mean temperature
!     28 maximum coldest monthly mean temperature
!     29 minimum growing degree days (at or above 5 deg C)
!     30 upper limit of temperature of the warmest month 
!     31 lower limit of growth efficiency (g/m2) REMARQUE: Ce parametre n'est pas pris en compte

!     INDIVIDUAL PARAMETERS

!     32 crown length
!     33 bark thickness 1
!     34 bark thickness 2
!     35 height slope (hmax_r-hsap_r)/(hmax_m-hsap_m)
!     36 height intercept (hsap_r - hs * hsap_m)
!     37 diameter slope (cm diameter gained per m height)
!     38 diameter intercept

!     FIRE PARAMETERS

!     39 flammability threshold													! call in fire.f, pour calculer un facteur de pondération de l'humidité de la litière
!     40 fire resistance index 													! call in fire.f, permet de calculer la fraction of individuals in grid cell which die
!     41 scorch height parameters = F
!     42 crown damage 1 = RCK
!     43 crown damage 2 = p
!     44 ignition efficiency parameter
!     45 emission factor for CO2
!     46 emission factor for CO
!     47 emission factor for CH4
!     48 emission factor for VOC
!     49 emission factor for TPM
!     50 emission factor for NOx from F77 code
!     51 fuel bulk density (kg m-3)

! Lecture de la table des parametres a partir d'un fichier csv
! Emeline Chaste, avril 2016
      
	call getarg(4,pftfile)  
	  
    open (unit = 15, file=pftfile, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
	  
	!write(*,*) 'jai ouvert la table'

	if (ierror .ne. 0) then
		write(*,*) 'File pft parameters cannot be open'
	else 
        read (15, *)
        do pft = 1, npft
		!write(*,*) pft 

				read (15, *) table(pft,:)
				!write(*,*) table(pft,:)
				
				do n=1,npftpar
					pftpar(pft,n)=table(pft,n)
					!write(*,*) table(pft,n)
				end do
        
        end do
    end if
		  		 
    close (unit = 15)
	
	!write(*,*) 'jai correctement lu les parametres des PFTs'

! Déchiffrage des parametres des PFTs

      do pft=1,npft

!       Transfer parameter values to array pftpar

!       do n=1,npftpar
!          pftpar(pft,n)=table(pft,n)
!        enddo
        
!       Assign leaf and phenology logicals

        if (pftpar(pft,15).le.2.0) then
          tree(pft)=.true.
          if (pftpar(pft,14).eq.2.0) then
            needle(pft)=.true.
          else
            needle(pft)=.false.
          endif
        else
          tree(pft)=.false.
          needle(pft)=.false.
        endif

        if (pftpar(pft,15).eq.1.0) then
          evergreen(pft)=.true.
          summergreen(pft)=.false.
          raingreen(pft)=.false.
        elseif (pftpar(pft,15).eq.2.0) then
          evergreen(pft)=.false.
          summergreen(pft)=.true.
          raingreen(pft)=.false.
        elseif (pftpar(pft,15).eq.3.0) then
          evergreen(pft)=.false.
          summergreen(pft)=.false.
          raingreen(pft)=.true.
        else
          evergreen(pft)=.true.
          summergreen(pft)=.true.
          raingreen(pft)=.true.
        endif
        
        if (pftpar(pft,21).eq.1.0) then
          boreal(pft)=.true.
        else
          boreal(pft)=.false.
        endif            

!       Calculate specific leaf area (SLA) for each PFT from leaf longevity
!       Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
!       Equation based on Reich et al 1997, Fig 1f:

!       SLA = 2e-4 * exp(6.15 - 0.46 ln (leaf_longevity * 12))

!       SLA in m2/gC, leaf_longevity in years

        sla(pft)=2.e-4*exp(6.15-0.46*log(pftpar(pft,7)*12.))

!       Define initial mass structure

        lai_sapl=pftpar(pft,19)

        if (tree(pft)) then  !woody PFTs

!         Calculate leafmass for a sapling individual
!          (1) lai = leafmass * sla / (crown area)
!          (2) (leaf area) = latosa * (sapwood xs area)
!                 (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
!          (3) (crown area) = allom1 * (stem diameter) ** reinickerp
!                 (Reinickes theory)
!         From (1),
!          (4) leafmass = lai * (crown area) / sla
!         From (1) & (3),
!          (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
!         From (2),
!          (6) leafmass = latosa * (sapwood xs area) / sla
!          (7) (sapwood xs area) = pi * (sapwood diameter)**2 / 4
!         From (6) and (7),
!          (8) leafmass = latosa * pi * (sapwood diameter)**2 / 4 / sla
!         From (8),
!          (9) (sapwood diameter) = [ 4 * leafmass * sla / pi / latosa ]**0.5
!         (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
!         Define x,
!         (11) x = [ (sapwood diameter)+(heartwood diameter) ] / 
!                  (sapwood diameter)
!         From (10) & (11),
!         (12) (stem diameter) = x * (sapwood diameter)
!         From (5), (9) & (12),
!         (13) leafmass = lai * allom1 * x**reinickerp * 
!                       (4*leafmass*sla/pi/latosa)**(reinickerp*0.5) / sla
!         From (13),
!         (14) leafmass = [ lai * allom1 * x**reinickerp *
!                (4*sla/pi/latosa)**(reinickerp*0.5) / sla ]**(1-1/reinickerp)

          x=pftpar(pft,20)

          lm_sapl(pft,1)=(lai_sapl*allom1*x**reinickerp*(4.0*sla(pft)/pi/latosa)**(reinickerp*0.5) / sla(pft))**(1.0-1.0/reinickerp)  !eqn 14
          lm_sapl(pft,2)=17.8-co2(2)   ! initial 13C value from llyod & farquhar, 1994
          lm_sapl(pft,3)=0.

!         Calculate sapling stem diameter
!         From (9) & (12),
!         (15) (stem diameter) = x * [ 4 * leafmass * sla / pi / latosa ]**0.5

          stemdiam=x*(4.0*lm_sapl(pft,1)*sla(pft)/pi/latosa)**0.5  !Eqn 15

!         Calculate sapling height
!         (16) height = allom2 * (stem diameter)**allom3 (source?)

          height_sapl=allom2*stemdiam**allom3   !Eqn 16

!		  Calculate sapling sapwood mass
!         (17) (sapwood volume) = height * (sapwood xs area)
!         (18) (sapwood xs area) = leafmass * sla / latosa
!         From (17) & (18),

!     (19) (sapwood volume) = height * leafmass * sla / latosa
!         (20) (sapwood mass) = (wood density) * (sapwood volume)
!         From (19) & (20),
!         (21) (sapwood mass) = (wood density) * height * leafmass * sla /
!                latosa

         sm_sapl(pft,1)=wooddens*height_sapl*lm_sapl(pft,1)*sla(pft)/ &
     &      latosa   !Eqn 21
          sm_sapl(pft,2)=lm_sapl(pft,2)   ! 13C value in permille
          sm_sapl(pft,3)=lm_sapl(pft,3)

!         Calculate sapling heartwood mass
!         From (11),
!         (22) (heartwood mass) = (x-1) * (sapwood mass)

          hm_sapl(pft,1)=(x-1.0)*sm_sapl(pft,1)  !Eqn 22
          hm_sapl(pft,2)=sm_sapl(pft,2)   ! 13C value in permille
          hm_sapl(pft,3)=sm_sapl(pft,3)

        else ! grass PFT

          lm_sapl(pft,1)=lai_sapl/sla(pft)

!         Set initial 13C values for saplings, grass

          if (pftpar(pft,2).eq.1) then   !C4 plants
            lm_sapl(pft,2)=3.6-co2(2)  !from lloyd & farquhar,1994           
            lm_sapl(pft,3)=0.
          else                           !C3 plpants
            lm_sapl(pft,2)=17.8-co2(2)  !from lloyd & farquhar,1994          
            lm_sapl(pft,3)=0.
          endif

          sm_sapl(pft,2)=0.0             ! no sapwood and hartwood
          hm_sapl(pft,2)=0.0             ! for grass PFT
          sm_sapl(pft,3)=0.0             ! no sapwood and hartwood
          hm_sapl(pft,3)=0.0             ! for grass PFT

        endif

!       Calculate sapling or initial grass rootmass
!       (23) lmtorm = (leafmass) / (rootmass)

        lmtorm=pftpar(pft,16) 
        rm_sapl(pft,1)=(1.0/lmtorm)*lm_sapl(pft,1)  !From Eqn 23
        rm_sapl(pft,2)=lm_sapl(pft,2)       ! 13C value in permille

      !write(0,*)pft,1000.*stemdiam,height_sapl,lai_sapl,lm_sapl(pft,1)

      enddo ! pft loop

      !write(0,*)'done pftinitassign'
      !read(*,*)

      return
      end
