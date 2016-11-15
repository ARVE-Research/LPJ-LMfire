c     SUBROUTINE ISOTOPE
c     This subroutine is for calculating the total fractionation of 13C
c     as it goes from free air (as 13CO2) to fixed carbon in the leaf.
c     This program is based upon the model used by Lloyd and Farquhar (1994).

      subroutine isotope(pft,cratio,co2,mtemp,rd,c4,mgpp,agpp)


      implicit none

c     PARAMETERS
      integer npft,ncvar
      parameter (npft=9,ncvar=3)
      real a,es,a1,b,b3,f,e
      parameter (a=4.4,es=1.1,a1=0.7,b=27.5,b3=30.0,e=0.0,f=8.0)


   

c     ARGUMENTS
      integer pft
      logical c4(1:npft)
      real cratio(1:12,1:npft),co2(1:ncvar),mtemp(1:12)
      real rd(1:12,1:npft,1:ncvar),mgpp(1:12,1:npft,1:ncvar)
      real agpp(1:npft,1:ncvar)
      

c     LOCAL
      integer m     
      real deltaa,b4
      real leaftemp,gamma
      real catm,k,rdi
      real q,r,s,t
      real phi
      real c3gpp,c4gpp

      c3gpp=0.0
      c4gpp=0.0
      catm=0.0
c      phi=0.2                 !as suggested by Lloyd & Farquhar (1994)
   
      do m=1,12

        if (mgpp(m,pft,1).gt.0.0) then

           if (cratio(m,pft).lt.0.05) cratio(m,pft)=0.05
          
c          define fractionation parameters
     
           if (rd(m,pft,1).le.0.0) rd(m,pft,1)=0.01
           
           leaftemp = 1.05*(mtemp(m)+2.5)
           gamma = 1.54*leaftemp
           rdi = rd(m,pft,1)/(86400.0*12.0)
           catm = co2(1)*1.0e6
           k = rdi/11.0         
           b4=(26.19-(9483/(273.2+leaftemp))) !from Farquhar et al. 1982 p. 126
           
           if (c4(pft)) then
          
c             calculate PHI for C4 pathway 

              call calcphi(mgpp,pft,phi)

c             calculate fractionation -> C4 pathway           
              
              deltaa=a*(1-cratio(m,pft)+0.0125)+0.0375*(es+a1)+
     *             (b4+(b3-es-a1)*phi)*(cratio(m,pft)-0.05)
                            
              mgpp(m,pft,2)=deltaa-co2(2)
              agpp(pft,2)=agpp(pft,2)+mgpp(m,pft,2)*mgpp(m,pft,1)
              c4gpp=c4gpp+mgpp(m,pft,1)

           else

c             calculate the fractionation -> C3 pathway
              
              q = a*(1-cratio(m,pft)+0.025)
              r = 0.075*(es+a1)
              s = b*(cratio(m,pft)-0.1)
              t = (e*rdi/k+f*gamma)/catm

              deltaa = q+r+s-t

              mgpp(m,pft,2)=deltaa-co2(2)             
              agpp(pft,2)=agpp(pft,2)+mgpp(m,pft,2)*mgpp(m,pft,1)
              c3gpp=c3gpp+mgpp(m,pft,1)

           endif
           
        else       
           mgpp(m,pft,2)=0.0
        endif
      enddo                     !month

      agpp(pft,2)=agpp(pft,2)/(c3gpp+c4gpp)
            
      return
      end									

c--------------------------------------------------------------------------------

c     SUBROUTINE CALCPHI
c     This subroutine is for calculating the phi variable used in 
c     C4 photosynethsis isotope fractionation calculations

      subroutine calcphi(mgpp,pft,phi)

      implicit none

c     ARGUMENTS
      real mgpp(1:12,1:9,1:3),phi
      integer pft


C     LOCAL
      real snormavg(1:4),svar(1:4)
      real avar,a,gppm(1:12)
      real totgpp,meangpp,normgpp(1:12)
      integer m,s

      do m=1,12
         gppm(m)=mgpp(m,pft,1)
      enddo
      totgpp=0.0       !initialize a few variables
      do s=1,4
        svar(s)=0.0
      enddo

c     This first part of the subroutine estimates annual variability of
c     GPP first by normalizing and then summing seasonal variability
c     which compensates for amplitude and seasonal variation in GPP.

      do m=1,12
        totgpp=totgpp+gppm(m)
      enddo

      meangpp=totgpp/12.0

      do m=1,12
        normgpp(m)=gppm(m)/meangpp
      enddo 

      snormavg(1)=(normgpp(12)+normgpp(1)+normgpp(2))/3.0
      snormavg(2)=(normgpp(3)+normgpp(4)+normgpp(5))/3.0
      snormavg(3)=(normgpp(6)+normgpp(7)+normgpp(8))/3.0
      snormavg(4)=(normgpp(9)+normgpp(10)+normgpp(11))/3.0

c     calculate the population variances by season

      do m=1,2
        a=((normgpp(m)-snormavg(1))**2)/3
        svar(1)=svar(1)+a
      enddo
      svar(1)=svar(1)+(((normgpp(12)-snormavg(1))**2)/3)

      do m=3,5
        a=((normgpp(m)-snormavg(2))**2)/3
        svar(2)=svar(2)+a
      enddo

      do m=6,8
        a=((normgpp(m)-snormavg(3))**2)/3
        svar(3)=svar(3)+a
      enddo

      do m=9,11
        a=((normgpp(m)-snormavg(4))**2)/3
        svar(4)=svar(4)+a
      enddo

      avar=svar(1)+svar(2)+svar(3)+svar(4)

c     ------------------------------------------------

c     This part sets the phi value based upon the annual variability.
c     The equation is a simple regresion based upon hypothetical extreme
c     scenarios of phi.

      phi=0.3518717*avar+0.2552359


      return

      end 

