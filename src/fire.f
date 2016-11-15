c     SUBROUTINE FIRE
c     Biomass destruction through disturbance by fire

      subroutine fire(pftpar,dtemp,litter_ag_fast,litter_ag_slow,
     *  acflux_fire,afire_frac,lm_ind,rm_ind,sm_ind,hm_ind,nind,dw1,
     *  present,tree,year)   

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar,ncvar
        parameter (npft=9,npftpar=51,nsoilpar=7,ncvar=3)
      real pi
        parameter (pi=3.14159265)
      real minfuel
        parameter (minfuel=200.0)  !fuel threshold to carry a fire (gC/m2)

c     ARGUMENTS
      real pftpar(1:npft,1:npftpar)
      real dtemp(1:365)
      real litter_ag_fast(1:npft,1:ncvar)
      real litter_ag_slow(1:npft,1:ncvar)
      real acflux_fire(1:ncvar)
      real lm_ind(1:npft,1:ncvar),rm_ind(1:npft,1:ncvar)
      real sm_ind(1:npft,1:ncvar),hm_ind(1:npft,1:ncvar)
      real nind(1:npft)
      real dw1(1:365)
      real afire_frac
      logical present(1:npft),tree(1:npft)
      integer year
    

c     LOCAL VARIABLES
      integer pft,d,c
      real fire_length,fuel,fire_prob
      real fire_index,disturb
      real resist(1:npft)
      real moistfactor,flam,litter_ag_total
      real fire_term
      real temp,acflux_fire_int(1:npft,1:ncvar)
      real all_ind,litter_inc

c     Calculate total above-ground litter

      litter_ag_total=0.0
      do pft=1,npft
        litter_ag_total=litter_ag_total+litter_ag_fast(pft,1)
     *        +litter_ag_slow(pft,1)
      enddo

c     Calculate litter moisture weighting factor

      moistfactor=0.0

      do pft=1,npft
        flam=pftpar(pft,39)
        if (litter_ag_total.gt.0.0) then
          moistfactor=moistfactor+((litter_ag_fast(pft,1)
     *          +litter_ag_slow(pft,1))/litter_ag_total)*flam    
        else
          moistfactor=0.0
        endif   
      enddo
c     First calculate the annual fraction of the grid cell affected by fire

c     Initialise

      fire_length=0.0
      fuel=0.0

      do c=1,ncvar
         acflux_fire(c)=0.0
         do pft=1,npft
            acflux_fire_int(pft,c)=0.0
         enddo
      enddo

c     Assign a minimum fire fraction (for presentational purposes)

      afire_frac=0.001

c     Assign PFT resistance to fire

      do pft=1,npft
        resist(pft)=pftpar(pft,40)
		
		write(*,*) 'resistance pft',resist(pft)
		
      enddo
       
c     Calculate the length of the fire season (units=days)

      do d=1,365

c       Calculate today's fire probability, fire_prob
c       Assume fire is only possible when temperature is above zero

        if (dtemp(d).gt.0.0.and.moistfactor.gt.0.0) then
          fire_prob=EXP((-pi/4.0)*(dw1(d)/moistfactor)**2)
        else 
          fire_prob=0.0
        endif

        fire_length=fire_length+fire_prob

      enddo
      
      !write(0,*)'fire length',fire_length

c     Calculate annual fire index

      fire_index=fire_length/365.0 


c     Calculate the fraction of the grid cell affected by fire

c      afire_frac=1.0-(EXP(-0.2*fire_index**1.4))**1.1 
      
      fire_term=fire_index-1.0
      afire_frac=fire_index*EXP( fire_term/(-0.13*fire_term**3 
     *    + 0.6*fire_term**2 + 0.8*fire_term + 0.45))
      
c     Calculate the available fuel (above-ground litter) to carry the fire

      fuel=litter_ag_total
            
c     Reduce fraction of grid cell affected by fire when fuel
c     becomes limiting (reduced carrying capacity)
      
      if (fuel.lt.minfuel) then
c        afire_frac=afire_frac*(1.0/(minfuel**2))*fuel**2   
        afire_frac=0.0
      endif
      
      if (afire_frac.lt.0.001) afire_frac=0.001   

c     Implement the effect of the fire on vegetation structure and litter 
c     in the disturbed fraction.  

c     Each PFT is assigned a resistance to fire, representing the fraction of
c     the PFT which survives a fire. Grasses assumed already to have completed
c     their life cycle and thus are not affected by fire, giving them
c     a competitive advantage against woody PFTs.

      do pft=1,npft

        if (present(pft).and.tree(pft)) then
     
        

c         Calculate the fraction of individuals in grid cell which die

           disturb=(1.0-resist(pft))*afire_frac
           
c          Calculate carbon flux to atmosphere (gC/m2) due to burnt biomass
           
           all_ind=lm_ind(pft,1)+sm_ind(pft,1)+hm_ind(pft,1)
     *          +rm_ind(pft,1)
           acflux_fire_int(pft,1)=disturb*(nind(pft)*all_ind)

           if (all_ind .gt. 0.) then

           do c=2,ncvar
              acflux_fire_int(pft,c)=(lm_ind(pft,c)*lm_ind(pft,1)
     *          +sm_ind(pft,c)*sm_ind(pft,1)+hm_ind(pft,c)*hm_ind(pft,1)
     *          +rm_ind(pft,c)*rm_ind(pft,1))/all_ind
           enddo
           end if
           
           
c          Update the individual density

           nind(pft)=nind(pft)*(1.0-disturb)
           
        endif

c       Add combusted litter to carbon flux to atmosphere term

        temp=acflux_fire_int(pft,1)
        acflux_fire_int(pft,1)=acflux_fire_int(pft,1)+(afire_frac
     *       *(litter_ag_fast(pft,1)+litter_ag_slow(pft,1)))
        if (acflux_fire_int(pft,1).gt.0.) then
           do c=2,ncvar
              litter_inc=litter_ag_fast(pft,c)*litter_ag_fast(pft,1)
     *             +litter_ag_slow(pft,c)*litter_ag_slow(pft,1)              
              acflux_fire_int(pft,c)=(acflux_fire_int(pft,c)*temp+
     *             afire_frac*litter_inc)/acflux_fire_int(pft,1)
           enddo
        else
           do c=2,ncvar
              acflux_fire_int(pft,c)=0.0
           enddo
        endif
        

c         Update the above ground litter term
         
        litter_ag_fast(pft,1)=(1.0-afire_frac)*litter_ag_fast(pft,1)
        litter_ag_slow(pft,1)=(1.0-afire_frac)*litter_ag_slow(pft,1)
        
c         Calculate carbon flux summed over all pft's
        
        temp=acflux_fire(1)
        acflux_fire(1)=acflux_fire(1)+acflux_fire_int(pft,1)        
        if (acflux_fire(1).gt.0.0) then         
           do c=2,ncvar
              acflux_fire(c)=(acflux_fire(c)*temp+acflux_fire_int(pft,c)
     *             *acflux_fire_int(pft,1))/acflux_fire(1)
           enddo
        else
           do c=2,ncvar
              acflux_fire(c)=0.0
           enddo
        endif
      
      enddo
      
      !write(0,*)'fire:',acflux_fire(1),acflux_fire_int(:,1)
      
      return
      end
