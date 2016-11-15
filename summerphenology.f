c     SUBROUTINE SUMMERPHENOLOGY
c     Temperature-based phenology for summergreen PFTs

      subroutine summerphenology(pftpar,mtemp,dtemp,gdd,dphen_t,
     *  summergreen,tree)

      implicit none

c     PARAMETERS
      integer npft,npftpar
        parameter (npft=9,npftpar=51)
      real gddbase
        parameter (gddbase=5.0)     !base temperature for growing degree days

c     ARGUMENTS
      real pftpar(1:npft,1:npftpar)
      real mtemp(1:12)
      real dtemp(1:365)
      real gdd
      real dphen_t(1:365,1:npft)
      logical summergreen(1:npft),tree(1:npft)
      integer year

c     LOCAL VARIABLES
      integer warmest,month,midsummer,firstday,pft,day
      real ramp
      real leafon
      integer ndaymonth(1:12)
      data (ndaymonth(month),month=1,12)
     *  / 31,28,31,30,31,30,31,31,30,31,30,31 /  
      integer cm1,cm2,cm3,coldest
      integer midday(1:12) 
      data (midday(month),month=1,12)
     *  / 16,44,75,105,136,166,197,228,258,289,319,350 /
      integer d
      real aphen       


c     First find warmest month

      warmest=1
      do month=1,12
        if (mtemp(month).gt.mtemp(warmest)) warmest=month
      enddo 

c     find the coldest month
      coldest=1
      do month=1,12
        if (mtemp(month).lt.mtemp(coldest)) coldest=month
      enddo

      midsummer=-15
      do month=1,warmest
        midsummer=midsummer+ndaymonth(month)
      enddo

c     Now midsummer is middle day (roughly) of warmest month
c     Find day of leaf abscission at end of summer

      firstday=midsummer+1
      do while (dtemp(firstday).ge.gddbase.and.
     *  firstday.ne.midsummer)
        firstday=firstday+1
        if (firstday.gt.365) firstday=1
      enddo

      do pft=1,npft
        
        if (summergreen(pft)) then  !summergreen taxa
          ramp=pftpar(pft,17)   !number of GDDs to attain full leaf cover

          if (firstday.eq.midsummer) then  !no leaf abscission

            do day=1,365
              dphen_t(day,pft)=1.0
            enddo

          else

            gdd=0.0     !accumulated growing degree days
            leafon=0.0  !proportional leaf-on today

            day=firstday+1
            if (day.gt.365) day=1
            do while (day.ne.firstday)
              if (dtemp(day).gt.gddbase) then  !growing day
                gdd=gdd+dtemp(day)-gddbase
                if (ramp.gt.0.0) then
                  leafon=min(gdd/ramp,1.0)
                else
                  leafon=1.0
                endif
              endif
              dphen_t(day,pft)=leafon
              day=day+1
              if (day.gt.365) day=1
            enddo              

          endif   


c         constrain woody deciduous phenology to max= 9 months 
          if (tree(pft)) then
             aphen=0.
             do day=1,365
                aphen=aphen+dphen_t(day,pft)  ! calculate the total number
                                              ! of days with foliage
             enddo 
             if (aphen.gt.210) then   ! limit summergreen to 210 
                                      ! days (approx 5 months) 
               do d=midday(coldest),midday(coldest)+75   ! set 75 days  
                                             ! after coldest with no foliage   
                  if (d.le.365) then
                     day=d
                  else
                     day=d-365      
                  endif
                  dphen_t(day,pft)=0.0
               enddo
               do d=midday(coldest)-75,midday(coldest)
                  if(d.ge.1) then
                     day=d
                  else
                     day=365+d
                  endif
                  dphen_t(day,pft)=0.0
               enddo
             endif
           endif          

        else  !non-summergreen taxa

          do day=1,365
            dphen_t(day,pft)=1.0
          enddo

        endif

      enddo   !pft

      return
      end
