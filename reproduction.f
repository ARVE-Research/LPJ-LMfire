c     SUBROUTINE REPRODUCTION
c     Deduction of reproduction costs from annual biomass increment

      subroutine reproduction(bm_inc,lm_sapl,sm_sapl,hm_sapl,rm_sapl,
     *     litter_ag_fast,litter_ag_slow,present,tree,co2)

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar,ncvar
      parameter (npft=9,npftpar=51,nsoilpar=7,ncvar=3)
      real reprod_cost              !proportion of NPP lost to reproduction
        parameter (reprod_cost=0.1) !Harper 1977

c     ARGUMENTS
      real bm_inc(1:npft,1:ncvar)
      real lm_sapl(1:npft,1:ncvar),sm_sapl(1:npft,1:ncvar)
      real hm_sapl(1:npft,1:ncvar),rm_sapl(1:npft,1:ncvar)
      real litter_ag_fast(1:npft,1:ncvar)
      real litter_ag_slow(1:npft,1:ncvar)
      logical present(1:npft)
      logical tree(1:npft)
      real co2(1:ncvar)

c     LOCAL VARIABLES
      integer pft,c
      real reprod        !allocation to reproduction (gC/m2)
      real temp

      do pft=1,npft

        if (present(pft)) then

c         Calculate allocation to reproduction
c         Reproduction costs taken simply as a constant fraction of annual NPP

          reprod=max(bm_inc(pft,1)*reprod_cost,0.0)

c         assume the costs go to reproductive structures which will
c         eventually enter the litter pool

          temp=litter_ag_fast(pft,1)
          litter_ag_fast(pft,1)=litter_ag_fast(pft,1)+reprod
          if (litter_ag_fast(pft,1).gt.1.) then
             do c=2,ncvar
                litter_ag_fast(pft,c)=(litter_ag_fast(pft,c)*temp+
     *               bm_inc(pft,c)*reprod)/(litter_ag_fast(pft,1))
             enddo
          endif
     

c         Reduce biomass increment by reproductive cost

          bm_inc(pft,1)=bm_inc(pft,1)-reprod

          if (tree(pft)) then
            if (bm_inc(pft,1).gt.0.0) then
               do c=2,ncvar
                  lm_sapl(pft,c)=bm_inc(pft,c)
                  sm_sapl(pft,c)=bm_inc(pft,c)
                  hm_sapl(pft,c)=bm_inc(pft,c)
                  rm_sapl(pft,c)=bm_inc(pft,c)
               enddo
            else              
               lm_sapl(pft,2)=17.8-co2(2) !C3 value         
               sm_sapl(pft,2)=17.8-co2(2) !from lloyd & farquhar,1994
               hm_sapl(pft,2)=17.8-co2(2)
               rm_sapl(pft,2)=17.8-co2(2)
               lm_sapl(pft,3)=co2(3)
               sm_sapl(pft,3)=co2(3)
               hm_sapl(pft,3)=co2(3)
               rm_sapl(pft,3)=co2(3)
            endif
          else
            if (bm_inc(pft,1).gt.0.0) then
               do c=2,ncvar
                  lm_sapl(pft,c)=bm_inc(pft,c)
                  rm_sapl(pft,c)=bm_inc(pft,c)
               enddo
            else
               lm_sapl(pft,2)=3.6-co2(2)  !C4 value
               rm_sapl(pft,2)=3.6-co2(2)  !from lloyd & farquhar,1994 
               lm_sapl(pft,3)=co2(3)
               rm_sapl(pft,3)=co2(3)
            endif
          endif
        endif

      enddo

      return
      end
