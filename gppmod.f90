module gppmod

implicit none

public :: calcgpp

contains 

subroutine calcgpp(present,co2,soilpar,pftpar,lai_ind,fpc_grid,mdayl,mtemp,mpar_day,dphen_t,w,dpet,dprec,dmelt,sla,   &
                   agpp,alresp,arunoff_surf,arunoff_drain,arunoff,mrunoff,dwscal365,dphen_w,dphen,wscal,mgpp,mlresp,  &
                   mw1,dw1,aaet,leafondays,leafoffdays,leafon,tree,raingreen,year,mat20,wscal_v,idx)

!Calculation of GPP, explicitly linking photosynthesis and water balance through canopy conductance feedback

use parametersmod,   only : sp,npft,npftpar,nsoilpar,ncvar,ndaymonth,i8
use waterbalancemod, only : waterbalance
use weathergenmod,   only : rmsmooth,daily

implicit none

integer(i8) :: idx

!parameters

real(sp), parameter :: epsilon = 0.05   !min precision of solution in bisection method
real(sp), dimension(npft) :: lambdam ! optimal Ci/Ca ratio

!arguments

integer,  intent(in) :: year
real(sp), intent(in) :: mat20  !20 year running mean annual temperature

logical,  dimension(:), intent(in) :: present
logical,  dimension(:), intent(in) :: tree
logical,  dimension(:), intent(in) :: raingreen
real(sp), dimension(:), intent(in) :: co2
real(sp), dimension(:), intent(in) :: soilpar
real(sp), dimension(:), intent(in) :: lai_ind
real(sp), dimension(:), intent(in) :: fpc_grid
real(sp), dimension(:), intent(in) :: mdayl
real(sp), dimension(:), intent(in) :: mtemp
real(sp), dimension(:), intent(in) :: mpar_day
real(sp), dimension(:), intent(in) :: dpet
real(sp), dimension(:), intent(in) :: dmelt
real(sp), dimension(:), intent(in) :: dprec
real(sp), dimension(:), intent(in) :: sla

real(sp), dimension(:,:), intent(in) :: pftpar

real(sp), intent(out) :: aaet
real(sp), intent(out) :: arunoff
real(sp), intent(out) :: arunoff_surf
real(sp), intent(out) :: arunoff_drain

real(sp), dimension(:), intent(out) :: mw1
real(sp), dimension(:), intent(out) :: mrunoff
real(sp), dimension(:), intent(out) :: dw1
integer,  dimension(:), intent(out) :: leafondays
integer,  dimension(:), intent(out) :: leafoffdays

real(sp), dimension(:,:),   intent(out) :: agpp
real(sp), dimension(:,:),   intent(out) :: alresp
real(sp), dimension(:,:,:), intent(out) :: mgpp
real(sp), dimension(:,:,:), intent(out) :: mlresp

real(sp), dimension(:),   intent(inout) :: w
logical,  dimension(:),   intent(inout) :: leafon
real(sp), dimension(:),   intent(inout) :: wscal
real(sp), dimension(:),   intent(inout) :: dwscal365
real(sp), dimension(:,:), intent(inout) :: dphen_w
real(sp), dimension(:,:), intent(inout) :: dphen_t
real(sp), dimension(:,:), intent(inout) :: dphen
real(sp), dimension(:,:), intent(inout) :: wscal_v

!local variables

integer  :: m
integer  :: b
integer  :: d
integer  :: dm
integer  :: i
integer  :: c
integer  :: pft
integer  :: dyr
real(sp) :: ca
real(sp) :: fpar
real(sp) :: rd
real(sp) :: agd
real(sp) :: drunoff_surf
real(sp) :: drunoff_drain
real(sp) :: gpd
real(sp) :: x1
real(sp) :: x2
real(sp) :: xmid
real(sp) :: fmid
real(sp) :: rtbis
real(sp) :: dx
real(sp) :: adt1
real(sp) :: adt2
real(sp) :: adtmm
real(sp) :: daet

real(sp), dimension(2)    :: ksat
real(sp), dimension(2)    :: awc
real(sp), dimension(2)    :: bound
real(sp), dimension(12)   :: mrunoff_surf
real(sp), dimension(12)   :: mrunoff_drain
real(sp), dimension(12)   :: gp
real(sp), dimension(12)   :: tsecs
real(sp), dimension(365)  :: dval

logical,  dimension(npft) :: C4
integer,  dimension(npft) :: aleafdays
real(sp), dimension(npft) :: awscal
real(sp), dimension(npft) :: dwscal
real(sp), dimension(npft) :: minwscal
real(sp), dimension(npft) :: gminp
real(sp), dimension(npft) :: nmax
real(sp), dimension(npft) :: longevity
real(sp), dimension(npft) :: inhibx1
real(sp), dimension(npft) :: inhibx2
real(sp), dimension(npft) :: inhibx3
real(sp), dimension(npft) :: inhibx4

real(sp), dimension(2,npft)   :: rootprop
real(sp), dimension(12,npft)  :: meanfpc
real(sp), dimension(12,npft)  :: meangc
real(sp), dimension(12,npft)  :: meangmin
real(sp), dimension(12,npft)  :: Cratio
real(sp), dimension(365,npft) :: dgp
real(sp), dimension(365,npft) :: dgc
      
!------------------------------------------
!initializations

ksat(1) = soilpar(1)
ksat(2) = soilpar(2)
awc(1)  = soilpar(3)
awc(2)  = soilpar(4)

!write(0,*) 'CO2 in gppmod: ', co2(1)

ca = co2(1) * 1.e-6  !from ppmv to mole fraction

mrunoff = 0.
mrunoff_surf  = 0.
mrunoff_drain = 0.

aaet = 0.
arunoff_surf  = 0.
arunoff_drain = 0.

!------------------------------------------

do pft=1,npft

  lambdam(pft) = pftpar(pft,26)
      
  longevity(pft)=pftpar(pft,7)

  !define pft inhibition function parameters

  inhibx1(pft)=pftpar(pft,22)
  inhibx2(pft)=pftpar(pft,23)
  inhibx3(pft)=pftpar(pft,24)
  inhibx4(pft)=pftpar(pft,25)

  if (present(pft)) then

    aleafdays(pft) = 0

    awscal(pft) = 0.
    minwscal(pft) = pftpar(pft,3)

    if (pftpar(pft,2) == 1.) then
      C4(pft) = .true.
    else
      C4(pft) = .false.
    end if

    !gminp = PFT-specific min canopy conductance scaled by fpc assuming full leaf cover

    gminp(pft) = pftpar(pft,4) * fpc_grid(pft)

    rootprop(1,pft) = pftpar(pft,1)
    rootprop(2,pft) = 1. - pftpar(pft,1)

    agpp   = 0.
    alresp = 0.

    nmax(pft) = pftpar(pft,6)

  end if

end do  !pft

!------------------------------------------

do pft = 1,npft

  if (present(pft)) then

    !find the potential canopy conductance realisable under non-water-stressed conditions

    do m = 1,12

      !Initialisations

      meanfpc(m,pft)  = 0.
      meangc(m,pft)   = 0.
      meangmin(m,pft) = 0.

      tsecs(m) = 3600. * mdayl(m)  !number of daylight seconds/day

      !Calculate non-water-stressed net daytime photosynthesis assuming full leaf cover

      fpar = fpc_grid(pft)

      call photosynthesis(ca,mtemp(m),fpar,mpar_day(m),mdayl(m),c4(pft),sla(pft),nmax(pft),lambdam(pft), &
                          rd,agd,adtmm,inhibx1(pft),inhibx2(pft),inhibx3(pft),inhibx4(pft),pft)

      if (tsecs(m) > 0.) then

        !Calculate non-water-stressed canopy conductance (gp) mm/sec basis averaged over entire grid cell
        !Eqn 21 Haxeltine & Prentice 1996

        gp(m)=(((1.6 * adtmm) / (ca * (1. - lambdam(pft)))) / tsecs(m)) + gminp(pft)

      else

        gp(m) = 0.

      end if

    end do  !month

    !Linearly interpolate mid-monthly gp to daily values

    bound(1) = gp(12)
    bound(2) = gp(1)

    call daily(gp,dval,.true.)
    !call rmsmooth(gp,ndaymonth,bound,dval)

    dgp(:,pft) = max(dval,0.)

  end if !present

end do !pft

!------------------------------------------
!calculate daily actual evapotranspiration and soil water balance
      
d = 1 !day of year

do m = 1,12

  mw1(m) = 0.

  do dm = 1,ndaymonth(m)

    do pft = 1,npft

      if (present(pft)) then

        !Use yesterday's potential water scalar to determine today's drought phenology

        if (d == 1) dwscal(pft) = dwscal365(pft)

        !Drought phenology and net phenology for today. Drought deciduous PFTs shed their leaves
        !when their water scalar falls below the PFT specific minimum value (minwscal).
        !Leaves are replaced immediately (i.e., daily) once the minimum water scalar is exceeded.

        !if (pft == 9) then
        !  write(0,'(a,2i5,3f10.3)')'gppmod',pft,d,dwscal(pft),minwscal(pft),dgp(d,9)
        !end if

        if (dwscal(pft) > minwscal(pft) .and. leafon(pft)) then

          dphen_w(d,pft)  = 1.
          dphen(d,pft)    = dphen_t(d,pft)
          leafondays(pft) = leafondays(pft) + 1

        else

          dphen_w(d,pft) = 0.
          dphen(d,pft)   = 0.

        end if

        !stop deciduous vegetation behaving like evergreen when climate permits

        if (raingreen(pft) .and. tree(pft)) then

          if (real(leafondays(pft)) >= (365. * longevity(pft))) then

            leafon(pft)      = .false.
            leafoffdays(pft) = leafoffdays(pft)+1

            if (real(leafoffdays(pft)) >= (365. * longevity(pft))) then

              leafoffdays(pft) = 0
              leafondays(pft)  = 0
              leafon(pft)      = .true.

            end if

          end if
        end if

      endif
    end do  !pft

    call waterbalance(d,present,rootprop,w,dgp,dpet,dphen,dgc,dmelt,dprec,ksat,awc,  &
                      drunoff_drain,drunoff_surf,dwscal,daet,fpc_grid,mat20,idx)

    !Store today's water content in soil layer 1

    dw1(d) = w(1)

    !Increment monthly runoff totals

    mrunoff_surf(m)  = mrunoff_surf(m)  + drunoff_surf
    mrunoff_drain(m) = mrunoff_drain(m) + drunoff_drain

    aaet = aaet + daet

    !Increment monthly w(1) total

    mw1(m) = mw1(m) + w(1)

    do pft = 1,npft
      if (present(pft)) then

        !Accumulate count of days with some leaf cover and pft-specific annual water scalar used in allocation

        if (dphen(d,pft) > 0.) then

          aleafdays(pft) = aleafdays(pft) + 1
          awscal(pft)    = awscal(pft) + dwscal(pft)

        end if 

        !Accumulate mean monthly fpc, actual (gc) and minimum (gmin) canopy conductances,
        !incorporating leaf phenology

        meangc(m,pft)   = meangc(m,pft)   + dgc(d,pft) / real(ndaymonth(m))
        meangmin(m,pft) = meangmin(m,pft) + gminp(pft) * dphen(d,pft) / real(ndaymonth(m))
        meanfpc(m,pft)  = meanfpc(m,pft)  + fpc_grid(pft) * dphen(d,pft) / real(ndaymonth(m))

        wscal_v(d,pft) = dwscal(pft)

        !Save final daily water scalar for next year

        if (d == 365) dwscal365(pft) = dwscal(pft)

      end if

    end do !pft
     
    d = d + 1

  end do !day of month

  !Increment annual runoff totals

  mrunoff(m)    = mrunoff_surf(m) ! + mrunoff_drain(m)
  arunoff_surf  = arunoff_surf  + mrunoff_surf(m)
  arunoff_drain = arunoff_drain + mrunoff_drain(m)

  !calculate gpp for each pft
  !Find water-limited daily net photosynthesis (And) and ratio of intercellular to ambient partial pressure of CO2
  !(lambda) by solving simultaneously Eqns 2, 18 and 19 (Haxeltine & Prentice 1996). 

  !Using a tailored implementation of the bisection method with a fixed 10 bisections, assuming root (f(lambda)=0)
  !bracketed by f(0.02) < 0 and f(lambdam + 0.05) > 0

  do pft = 1,npft 

    if (present(pft)) then

      !Convert canopy conductance assoc with photosynthesis (actual minus minimum gc) (gpd) from mm/sec to mm/day

      gpd = tsecs(m) * (meangc(m,pft) - meangmin(m,pft))

      fpar = meanfpc(m,pft)  !cover including phenology  

      if (gpd > 1.e-5) then  !canopy conductance
            
        !Implement numerical solution

        x1 = 0.02                  !minimum bracket of the root
        x2 = lambdam(pft) + 0.05   !maximum bracket of the root
        rtbis = x1                 !root of the bisection
        dx = x2 - x1

        b = 0  !number of tries towards solution

        fmid = epsilon + 1.

        do  !bisection root finding

          b    = b + 1
          dx   = dx * 0.5
          xmid = rtbis + dx

          !Calculate total daytime photosynthesis implied by canopy conductance from water balance routine and
          !current guess for lambda (xmid).  Units are mm/m2/day (mm come from gpd value, mm/day)
          !Eqn 18, Haxeltine & Prentice 1996

          adt1 = gpd / 1.6 * ca * (1. - xmid)

          !Call photosynthesis to determine alternative total daytime photosynthesis estimate (adt2) implied by
          !Eqns 2 & 19, Haxeltine & Prentice 1996, and current guess for lambda (xmid)

          call photosynthesis(ca,mtemp(m),fpar,mpar_day(m),mdayl(m),c4(pft),sla(pft),nmax(pft),xmid, &
                              rd,agd,adt2,inhibx1(pft),inhibx2(pft),inhibx3(pft),inhibx4(pft),pft)

          !Evaluate fmid at the point lambda = xmid fmid will be an increasing function with xmid,
          !with a solution (fmid=0) between x1 and x2
        
          fmid = adt2 - adt1

          if (fmid < 0.) rtbis = xmid

          !exit if solution found or > 10 iterations needed without solution
 
          if (abs(fmid) < epsilon .or. b > 10) exit

        end do !bisection

      else  !infinitesimal canopy conductance

        rd   = 0.
        agd  = 0.
        xmid = 0.

      end if  !canopy conductance

      !Estimate monthly gross photosynthesis and monthly leaf respiration from mid-month daily values
      !Agd = And + Rd (Eqn 2 Haxeltine & Prentice 1996)

      mgpp(m,pft,1) = real(ndaymonth(m)) * agd
      agpp(pft,1)   = agpp(pft,1) + mgpp(m,pft,1)

      mlresp(m,pft,1) = real(ndaymonth(m)) * rd
      alresp(pft,1)   = alresp(pft,1) + mlresp(m,pft,1)
      cratio(m,pft)   = xmid

    else

      mgpp(m,pft,1) = 0.
      agpp(pft,1)   = 0.

    end if
          
  end do  !pft
        
  mw1(m) = mw1(m) / real(ndaymonth(m))

end do  !month
       
!------------------------------------------
!Convert water scalar to average daily basis using phenology

do pft = 1,npft
  if (present(pft)) then

    if (aleafdays(pft) /= 0) then
      wscal(pft)=awscal(pft) / real(aleafdays(pft))
    else
      wscal(pft)=1.
    end if

  end if
end do

!------------------------------------------
!Calculate the carbon isotope fractionation in plants following Lloyd & Farquhar 1994

do pft=1,npft
   if (present(pft)) then

      if (agpp(pft,1).gt.0.0) then

         call isotope(pft,cratio,co2,mtemp,mlresp,c4,mgpp,agpp)

         !assign DELTA 14C value from atmospheric time series to plant production
         mgpp(:,pft,3) = co2(3)
         agpp(pft,3)   = co2(3)

      else

         mgpp(:,pft,2:3) = 0.
         agpp(pft,2:3)   = 0.

      end if

   end if
end do

!increement annual total runoff

arunoff = arunoff_surf + arunoff_drain

!write(0,*)sum(meanfpc,dim=1)
!write(0,*)sum(meangc,dim=1)
!write(0,*)sum(meangmin,dim=1)

!write(0,*)'  gpp: dgc    ',sum(dgc,dim=1)
!write(0,*)'  gpp: wscal  ',wscal   
!write(0,*)'  gpp: present',present!sum(dphen_t,dim=1)
!write(0,*)'  gpp: awscal ',awscal
!write(0,*)'  gpp: alfdays',aleafdays

!write(0,'(24f12.2)')mw1,mrunoff
    
end subroutine calcgpp

!------------------------------

end module gppmod
