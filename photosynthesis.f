c     SUBROUTINE PHOTOSYNTHESIS
c     Adapted from Farquhar (1982) photosynthesis model, as simplified by
c     Collatz et al 1991, Collatz et al 1992 and Haxeltine & Prentice 1996

      subroutine photosynthesis(ca,temp,fpar,par,dayl,c4,sla,nmax,
     *  lambda,rd,agd,adtmm,x1,x2,x3,x4,pfti)

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar,ncvar,pft
        parameter (npft=9,npftpar=51,nsoilpar=7,ncvar=3)
	  real pftpar(1:npft,1:npftpar)
	  real po2,p,bc3,bc4,theta,q10ko,q10kc,q10tau,ko25,kc25,tau25
      real alphaa,alphac3,alphac4,cmass,lambdamc3,lambdamc4,n0,m
      real e0,t0c3,t0c4,cq,tk25,tmc4,tmc3
	  real lambdam
	
      parameter (po2=20.9e3)!O2 partial pressure in Pa
      parameter (p=1.0e5)   !atmospheric pressure in Pa
      parameter (bc3=0.015) !leaf respiration as fraction of Vmax for C3 plants
      parameter (bc4=0.02)  !leaf respiration as fraction of vmax for C4 plants
      parameter (theta=0.7) !colimitation (shape) parameter
      parameter (q10ko=1.2) !q10 for temperature-sensitive parameter ko
      parameter (q10kc=2.1) !q10 for temperature-sensitive parameter kc
      parameter (q10tau=0.57) !q10 for temperature-sensitive parameter tau
      parameter (ko25=3.0e4)  !value of ko at 25 deg C
      parameter (kc25=30.0)   !value of kc at 25 deg C
      parameter (tau25=2600.0)!value of tau at 25 deg C
      parameter (alphaa=0.5)  !fraction of PAR assimilated at ecosystem level
                              !relative to leaf level
      parameter (alphac3=0.08)  !intrinsic quantum efficiency of CO2 uptake in
                                !C3 plants
      parameter (alphac4=0.053) !C4 intrinsic quantum efficiency
      parameter (lambdamc4=0.4) !optimal ratio of intercellular to ambient CO2
                               !concentration (lambda) in C4 plants
      parameter (lambdamc3=0.8)!optimal (maximum) lambda in C3 plants
      parameter (cmass=12.0)    !atomic mass of carbon
      parameter (cq=4.6e-6)   !conversion factor for solar radiation at 550 nm
                              !from J/m2 to E/m2 (E=mol quanta)
      parameter (n0=7.15)     !leaf N concentration (mg/g) not involved in
                              !photosynthesis
      parameter (m=25.0)      !corresponds to parameter p in Eqn 28, Haxeltine
                              !& Prentice 1996
      parameter (t0c3=250.0)  !base temperature (K) in Arrhenius temperature
                              !response function for C3 plants
      parameter (t0c4=260.0)  !base temperature in Arrhenius func for C4 plants
      parameter (e0=308.56)   !parameter in Arrhenius temp response function
      parameter (tk25=298.15) !25 deg C in Kelvin
      parameter (tmc3=45.0)   !maximum temperature for C3 photosynthesis
      parameter (tmc4=55.0)   !maximum temperature for C4 photosynthesis
     

c     ARGUMENTS
      real lambda,par,dayl,temp,rd,agd,adtmm,ca,fpar,sla,nmax
      logical c4
      real x1,x2,x3,x4
      integer pfti

c     LOCAL VARIABLES
      real apar,pi,ko,kc,tau,gammastar,c1,c2,b,phipi
      real s,sigma,vm,and,pa,je,jc,vmmax,cn,tk,t0,adt
      real tstress,k1,k2,k3,low,high
      integer i

       do pft=1,npft
        lambdam=pftpar(pft,26) 
      enddo

c     Return without performing calculations if daylength 0 hours

      if (dayl.lt.0.01) then
        agd=0.0
        adtmm=0.0
        rd=0.0
        return
      endif

c     APAR in J/m2/day
c     alphaa = scaling factor for absorbed PAR at ecosystem, versus leaf, scale
c     See Eqn 4, Haxeltine & Prentice 1996

      apar=par*fpar*alphaa

c     calculate temperate inhibition function

      if (temp.lt.x4) then
         k1=2.*ALOG((1./0.99)-1.)/(x1-x2)
         k2=(x1+x2)/2.
         low=1./(1.+ EXP(k1*(k2-temp)))
         k3=ALOG(0.99/0.01)/(x4-x3)
         high=1.-0.01*EXP(k3*(temp-x3))
         tstress=(low*high)
      else
         tstress=0.
      endif
      
c      write(0,*)'tstress',temp,x4,tstress
      if (tstress.lt.1e-2) tstress=0.  

c     First calculate catalytic capacity of rubisco, Vm, assuming optimal
c     (non-water-stressed) value for lambda, i.e. lambdamc3

      if (.not.c4) then  !C3 photosynthesis

c       Temperature-adjusted values of kinetic parameters, Eqn 22,
c       Haxeltine & Prentice 1996a

        ko=ko25*q10ko**((temp-25.0)/10.0) !Michaelis constant of rubisco for O2
        kc=kc25*q10kc**((temp-25.0)/10.0) !Michaelis constant for CO2
        tau=tau25*q10tau**((temp-25.0)/10.0) !CO2/O2 specificity ratio

c       CO2 compensation point (CO2 partial pressure, Pa)
c       Eqn 8, Haxeltine & Prentice 1996

        gammastar=po2/(2.0*tau)
       

c       Convert ambient CO2 level, ca, from mole fraction to partial pressure
c       in Pa

        !write(0,*) 'CO2 in photosynthesis: ', ca * 1e6
        !stop

        pa=ca*p

c       Non-water-stressed intercellular CO2 partial pressure in Pa
c       Eqn 7, Haxeltine & Prentice 1996

        pi=lambdam*pa

c       Calculation of C1C3, Eqn 4, Haxeltine & Prentice 1996
c       Notes: - there is an error in this equation in the above paper (missing
c                2.0* in denominator) which is fixed here (see Eqn A2, Collatz
c                et al 1991)
c              - There is no longer an explicit temperature inhibition function
c                (low-temperature inhibition is now done mechanistically
c                by imposing a temperature-dependent upper limit on Vm, see
c                below)
c              - There is no longer any reduction in maximum photosynthesis due
c                to leaf age (phic)
c              - alphaa, the upscaling parameter accounting for the reduction
c                in PAR utilisation in ecosystems compared with leaf level,
c                appears in the calculation of APAR instead of here
c              - Cmass, the atomic weight of carbon, used in unit conversion
c                from molC to g appears in the calculation of Vm instead of
c                here
 
        c1=tstress*alphac3*((pi-gammastar)/(pi+2.0*gammastar))

c       High temperature inhibition modelled primarily by suppression of LUE
c       by decreased relative affinity of rubisco for CO2 relative to O2 with
c       increasing temperature, but we also implement a step function to
c       prohibit any C3 photosynthesis above 45 degrees (Table 3.7, Larcher
c       1983)

        if (temp.gt.tmc3) c1=0.0

c       Calculation of C2C3, Eqn 6, Haxeltine & Prentice 1996

        c2=(pi-gammastar)/(pi+kc*(1.0+po2/ko))

        b=bc3   !Choose C3 value of b for Eqn 10, Haxeltine & Prentice 1996
        t0=t0c3 !base temperature for temperature response of rubisco

      else  !C4 photosynthesis

c       Specify C1C4, C2C4
c       Eqns 14,15, Haxeltine & Prentice 1996
c       Notes:
c              - alphaa, the upscaling parameter accounting for the reduction
c                in PAR utilisation in ecosystems compared with leaf level,
c                appears in the calculation of APAR instead of here
c              - Cmass, the atomic weight of carbon, used in unit conversion
c                from molC to g appears in the calculation of Vm instead of
c                here
c              - parameter phipi is not needed for calculation of optimal Vm
c                which assumes optimal intercellular CO2 concentration
c                (lambdamc4)

        c1=tstress*alphac4

c       High-temperature inhibition modelled conservatively as a step function
c       prohibiting photosynthesis above 55 deg C (Table 3.7, Larcher 1983)

        if (temp.gt.tmc4) c1=0.0

        c2=1.0
        b=bc4    !Choose C4 value of b for Eqn 10, Haxeltine & Prentice 1996
        t0=t0c4  !base temperature for temperature response of rubisco

      endif

c     Eqn 13, Haxeltine & Prentice 1996

      s=(24.0/dayl)*b

c     Eqn 12, Haxeltine & Prentice 1996

      sigma=sqrt(max(0.0,1.0-(c2-s)/(c2-theta*s)))

c     Calculation of optimal rubisco capacity, Vm, in gC/m2/day
c     Eqn 11, Haxeltine & Prentice 1996

      vm=(1.0/b)*(c1/c2)*((2.0*theta-1.0)*s-(2.0*theta*s-
     *  c2)*sigma)*apar*cmass*cq

c     Now use this Vm value to calculate actual photosynthesis

      if (.not.c4) then  !C3 photosynthesis

c       Intercellular CO2 partial pressure in Pa
c       Eqn 7, Haxeltine & Prentice 1996

        pi=lambda*pa

c       Recalculation of C1C3, C2C3 with actual pi

        c1=tstress*alphac3*((pi-gammastar)/(pi+2.0*gammastar))
        if (temp.gt.tmc3) c1=0.0  !high-temperature inhibition

        c2=(pi-gammastar)/(pi+kc*(1.0+po2/ko))

      else  !C4 photosynthesis

c       Parameter accounting for effect of reduced intercellular CO2
c       concentration on photosynthesis, Phipi.
c       Eqn 14,16, Haxeltine & Prentice 1996
c       Fig 1b, Collatz et al 1992
 
        phipi=min(lambda/lambdam,1.0)
        c1=tstress*phipi*alphac4
        if (temp.gt.tmc4) c1=0.0  !high-temperature inhibition

      endif

c     je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
c     Convert je from daytime to hourly basis

c     Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
c     Eqn 3, Haxeltine & Prentice 1996

      je=c1*apar*cmass*cq/dayl

c     Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h
c     Eqn 5, Haxeltine & Prentice 1996

      jc=c2*vm/24.0

      !write(0,'(i5,4f12.2,4f10.4)')pfti,temp,par,fpar,dayl,je,jc,tstress,c1

      if (je.lt.1e-10.or.jc.le.1e-10) then


        agd=0.0

      else
c      write(0,*)'daily',je,jc,(je+jc)**2.0,4.0*theta*je*jc

c     Calculation of daily gross photosynthesis, Agd, gC/m2/day
c     Eqn 2, Haxeltine & Prentice 1996
c     Note: - there is an error in this equation in the above paper (missing
c             theta in 4*theta*je*jc term) which is fixed here

      agd=(je+jc-sqrt((je+jc)**2.0-4.0*theta*je*jc))/(2.0*theta)*dayl
      
      end if

c     Daily leaf respiration, Rd, gC/m2/day
c     Eqn 10, Haxeltine & Prentice 1996

      rd=b*vm

c     Daily net photosynthesis (at leaf level), And, gC/m2/day

      and=agd-rd

c     Total daytime net photosynthesis, Adt, gC/m2/day
c     Eqn 19, Haxeltine & Prentice 1996

      adt=and+(1.0-dayl/24.0)*rd

c     Convert adt from gC/m2/day to mm/m2/day using
c     ideal gas equation

      adtmm=adt/cmass*8.314*(temp+273.3)/p*1000.0
      
      return
      end
