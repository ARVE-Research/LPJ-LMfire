module foragersmod
  
public :: foragers
public :: popgrowth

contains

subroutine foragers(pet,aet,elev,lat,anpp,livebiomass,whc,mtemp,mprec,soilm,pd)

!minimalist model of forager population density developed by Binford (XXXX)
!Kristen Krumhardt & Jed Kaplan 04.2012

use parametersmod, only : sp,dp

implicit none

!arguments

real(sp),                intent(in) :: pet             !annual potential evapotranspiration (mm)
real(sp),                intent(in) :: aet             !annual actual evapotranspiration (mm)
real(sp),                intent(in) :: elev            !elevation (m)
real(dp),                intent(in) :: lat             !latitude (degrees)
real(sp),                intent(in) :: anpp            !annual NPP (g m-2)
real(sp),                intent(in) :: livebiomass     !primary biomass used in binford equations
real(sp),                intent(in) :: whc             !soil water holding capacity, top layer (mm)

!real(sp),                intent(in) :: wstorage        !ordinal scaling of soil type in terms of water holding capacity (not used in this version)

real(sp), dimension(12), intent(in) :: mtemp           !mean monthly temperature (deg C)
real(sp), dimension(12), intent(in) :: mprec           !mean monthly precipitation (mm)
real(sp), dimension(12), intent(in) :: soilm           !mean monthly soil water content (fraction of field capacity?) (called mw1 in lpjcore)

real(sp), intent(out)               :: pd              !potential population density

!parameters

!constants needed for the calculation of EXWGT, the weight of hunter-gatherers

real(sp), parameter :: a = 39.154974
real(sp), parameter :: b = 16.276180
real(sp), parameter :: c =  0.146003

!local variables

real(sp) :: elef   !elevation (ft)
real(sp) :: llat   !log10 of latitude
real(sp) :: lnpp   !log10 of net aboveground productivity
real(sp) :: watd   !annual water deficit, aet - pet
real(sp) :: lwatd  !log10 of watd

integer :: watrgrc   !number of months during the growing season in which water was retained in the soil at the end of the month
real(sp) :: lwatrgrc  !log10 of watrgrc

integer, dimension(1) :: m  !variable used for finding location in monthly data

integer :: warmest_month
integer :: coldest_month
integer :: driest_month
integer :: wettest_month

real(sp) :: rrcorr          !measure of seasonality of rainfall
real(sp) :: rrcorr2         !modified version of rrcorr, making the values positive (values 0-12)

real(sp) ::wret             !yearly total of water stored in the soil (mm)

real(sp) :: exprey          !expected prey biomass (km/km2)
real(sp) :: expreya

real(sp) :: trange          !temperature range
real(sp) :: ltrange         !log10 of trange        

real(sp) :: mrain           !measure of evenness of rainfall over the year

real(sp) :: exwgt           !expected body weight of hunter-gatherer person

real(sp) :: bar5            !biomass accumulation ratio

real(sp) :: exprim1         !expressions for calculating the number of people supported by plant food alone per 100km2
real(sp) :: exprim2
real(sp) :: exprim3

real(sp) :: growc           !effective growing season, number of consecutive months in which mean temperature exceeds 8 deg C

real(sp) :: et              !effective temperature

real(sp) :: termh2          !people supported per 100km2 by hunting
real(sp) :: termg2          !people supported per 100km2 by plant foods
real(sp) :: termd2          !potential population density
real(sp) :: plant_percent   !percent dependence on plants
real(sp) :: animal_percent  !percent dependence on animals

!end of var declarations

!-----------------------------------------
!calculations start here

!write(0,*)'enter foragers'
!write(0,*)pet,aet,elev,lat,anpp,livebiomass,whc,mtemp,mprec

!write(0,*) 'entering forager-routine'

if (anpp == 0.) then
  pd = 0.
  return
end if

!convert elevation to feet
elef = elev * 3.28084

!log of latitude
llat = log10(abs(lat))

!log of net aboveground production

lnpp = log10(anpp)

!WATD, total annual water deficit (page 109)
watd = max(pet - aet,0.)

!write(0,*) 'watd:', watd

lwatd = log10(watd + 1.)  !log of above 

!WATRGRC, count number of months during the growing season in which water was retained in the soil (page 109),
!months when mtemp is > 8. is considering "effective growing season"

watrgrc = count(soilm > 0. .and. mtemp > 8.)

if (watrgrc > 0) then
  lwatrgrc = log10(real(watrgrc))
else
  lwatrgrc = 0.
end if

!write(0,*) 'lwatrgrc', lwatrgrc

!for subsequent calculations using monthly extremes in temperature and precipitation

m = maxloc(mtemp)

warmest_month = m(1)

m = minloc(mtemp)

coldest_month = m(1)

m = maxloc(mprec)

wettest_month = m(1)

m = minloc(mprec)

driest_month = m(1)

!for calculating RRCORR2 (described on pg 71) ----- not sure about this...

rrcorr = real(wettest_month - warmest_month)

rrcorr2 = rrcorr + 4.5

if (rrcorr2 < 0.) rrcorr2 = 12. - rrcorr2
  
!WRET, yearly total of water stored in the soil (page 75, 109) ----- not sure about this...

wret = sum(soilm * whc)

!write(0,*) 'wret:', wret

!calculation of EXPREY (kg/km2)

exprey = 10.**(elef * 5.3081e-5 + llat * -0.300235 + lnpp * 1.200771 + lwatd * -0.11661 + lwatrgrc * 0.216493 + anpp * -4.26495e-4 + rrcorr2 * -0.028577 + wret * -0.008066) !+(wstorage * 0.005171))

expreya = 100 * (exprey + 0.01)

!write(0,*) 'exprey:', exprey

!------------------------------------------------------------------

!TRANGE = temperature range (page 59)

trange = mtemp(warmest_month) - mtemp(coldest_month)       

ltrange = log10(trange)

!write(0,*) 'ltrange:', ltrange, mprec(wettest_month)

!MRAIN, a measure of rainfall evenness (page 72)

mrain = mprec(driest_month)/max(mprec(wettest_month) * 100.,1e-6)		!FLAG MP: this led to a division-by-zero in some desert places with no rain at all, threfore the max

!write(0,*) 'mrain:', mrain

!calculation of EXWGT, anticipated weight of individual EXWGT (kg), see page 182

exwgt = a + (b * ltrange) + (c * mrain)

!write(0,*) 'exwgt:', exwgt

!variables for calculating TERMG2

!BAR5, biomass accumulation ratio, see page 85
bar5 = livebiomass/anpp

exprim1 = ((anpp / 1000.) * 1.e8) * (1. - (bar5/85.)) * (1. - (livebiomass / 61000.)**2) * (1. - (livebiomass / 61000.))

exprim2 = exprim1 * (1. - (exprey / 20000.)**2) * (1. - (anpp / 63000.))

!write(0,*) 'exprim2:', exprim2

!GROWC, effective growing season = number of consecutive months in which mean temp exceeds 8 deg C
growc =  count(mtemp > 8.) 
 
!ET, effective temperature
et = ((18. * mtemp(warmest_month)) - (10. * mtemp(coldest_month))) / (mtemp(warmest_month) - mtemp(coldest_month) + 8.)

exprim3 = exprim2 - exprim2 * (1. - (growc / 12.)**2) * 1. - ((et - 7.) / 23.) !unbalanced parentheses in text, pg 180

!write(0,*) 'exprim3:', exprim3

!MINIMALIST TERRESTRIAL MODEL

!TERMH2, pop dens (person/100km2) of people supported by ungulates alone

termh2 = (expreya * (1 - (exprey / 20000.)) * 0.026142) / (exwgt / 0.045)

!TERMG2, pop dens (person/100km2) of people supported by plant foods alone

termg2 = (exprim3 * 0.00006) / (exwgt / 0.043748)

!TERMD2, overall population density (persons/100km2)

termd2 = termh2 + termg2

pd = max((termd2 * 0.01),0.)  !convert to persons km-2	!FLAG MP: would go below zero in some low-productivity places, therefore constrained it

!write(0,*) 'pd:', pd

!Percent plant dependence

plant_percent = termg2/termd2 * 100.

!Percent animal dependence

animal_percent =  termh2/termd2 * 100.

!write(0,*) 'animal_percent:', animal_percent

!write(0,*)termh2,termg2,termd2,plant_percent,animal_percent

end subroutine foragers

!-----------------------------------------------------

subroutine popgrowth(ppd,pop,fout,fin)

use parametersmod, only : sp

implicit none

real(sp), intent(in)    :: ppd     !potential popoulation density at gridcell carrying capacity (persons km-2)
real(sp), intent(inout) :: pop     !actual popoulation density (persons km-2)
real(sp), intent(inout) :: fin     !flux of people into the gridcell (persons km-2 yr-1)
real(sp), intent(inout) :: fout

real(sp), parameter :: pgr = 8.e-5  !Population growth rate (yr-1) Livi-Bacci (2007) table 1.2

real(sp) :: pchange  !dP/dt change in population per time

!first calculate endogenous growth

!rP * (1-P/K)

!pchange = pgr * pop * (1. - pop / ppd)  !Verhulst equation FLAG MP: THIS IS NOT USED => purpose?

if(pop==0. .and. ppd==0.) then
  
  pop = 0.

else  

  pop = ppd * pop * exp(pgr) / (ppd + pop *(exp(pgr) - 1.)) 	!FLAG MP: caused divide-by-zero for cases when both ppd and pop equal zero 
    
end if    

!write(0,*)pop,ppd

!if below potential there is no out-migration

if (pop <= ppd) then
  
  fout = 0.

else
  !actual population > carrying capacity
  
  fout = pop - ppd
  pop = ppd
  fin = 0.

end if

end subroutine popgrowth 

!-----------------------------------------------------

subroutine simpleforagers(PD0,tree,anpp,fpc_grid,forager_pd)

!very simple subroutine to calculate potential density of foragers given herbaceous npp

use parametersmod, only : sp

implicit none

real(sp),               intent(in)    :: PD0    !baseline population density from kernel interpolation of archaeological sites
logical,  dimension(:), intent(in)    :: tree
real(sp), dimension(:), intent(in)    :: anpp
real(sp), dimension(:), intent(in)    :: fpc_grid
real(sp),               intent(out) :: forager_pd

!---

real(sp), parameter :: npp0  = 600.  !baseline herbaceous NPP for baseline PD
!real(sp), parameter :: PDmax = 1.646 !maximum forager PD under LGM conditions

real(sp) :: grassnppfrac
real(sp) :: grassnpp

!---
!calculate the total annual NPP of the herbaceous vegetation

grassnpp = sum(anpp * fpc_grid,mask=.not.tree)  !total NPP of the herbaceous vegetation in the cell

grassnppfrac = grassnpp / (npp0 * 0.50)

forager_pd = PD0 * grassnppfrac

end subroutine simpleforagers

!-----------------------------------------------------

end module foragersmod
