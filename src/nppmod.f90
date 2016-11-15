module nppmod

implicit none

public :: calcnpp

contains

!----------------------------------------------------------------------------------------------------------------------

subroutine calcnpp(tair,tsoil,dphen,present,nind,lm_ind,sm_ind,hm_ind,rm_ind,cstore,mgpp,mnpp,anpp)

use parametersmod,   only : sp,npft,pft,pftpar,ndaymonth,firstday

!calculate NPP based on temperature and plant mass

implicit none

!arguments

real(sp), dimension(:),   intent(in)    :: tair    !per day
real(sp), dimension(:),   intent(in)    :: tsoil
real(sp), dimension(:,:), intent(in)    :: dphen   !per day,pft
logical,  dimension(:),   intent(in)    :: present
real(sp), dimension(:),   intent(in)    :: nind
real(sp), dimension(:),   intent(in)    :: lm_ind
real(sp), dimension(:),   intent(in)    :: sm_ind
real(sp), dimension(:),   intent(in)    :: hm_ind
real(sp), dimension(:),   intent(in)    :: rm_ind
real(sp), dimension(:),   intent(inout) :: cstore  !per pft
real(sp), dimension(:,:), intent(inout) :: mgpp    !per month,pft
real(sp), dimension(:,:), intent(out)   :: mnpp    !per month,pft
real(sp), dimension(:),   intent(out)   :: anpp    !per pft

!parameters

real(sp), parameter :: k  = 0.0548
real(sp), parameter :: tc = 1. / 56.02
real(sp), parameter :: min_npp = 0.02

!local variables

integer :: i
integer :: m
integer :: a
integer :: b

real(sp) :: l_n2c      ! leaf tissue C:N mass ratio
real(sp) :: s_n2c      ! sapwood tissue C:N mass ratio
real(sp) :: r_n2c      ! fineroot tissue C:N mass ratio 
real(sp) :: respcoeff  ! coefficient in respiration equations

real(sp) :: agpp

real(sp) :: aresp_ind

real(sp) :: nind0
real(sp) :: nind_kill

real(sp), dimension(365) :: gtemp_air  ! value of temperature response function given air temperature
real(sp), dimension(365) :: gtemp_soil ! value of temperature response function given soil temperature
real(sp), dimension(365) :: lresp_ind
real(sp), dimension(365) :: rresp_ind
real(sp), dimension(365) :: sresp_ind

real(sp), dimension(12,npft) :: gresp    !per month,pft

!----------------------------
!calculate temperature coefficients

where (tair > -40.) 
  gtemp_air = exp(308.56 * (tc - 1. / (tair  + 46.02)))
elsewhere
  gtemp_air = 0.
end where

where (tsoil > -40.) 
  gtemp_soil = exp(308.56 * (tc - 1. / (tsoil + 46.02)))
elsewhere
  gtemp_soil = 0.
end where

!---
!calculate NPP

do i = 1,npft
  if (present(i)) then
    
    !C to N ratios of leaf, sapwood & roots

    l_n2c = 1. / pftpar(i,11)
    r_n2c = 1. / pftpar(i,13)

    if (pft(i)%tree) then
      s_n2c = 1. / pftpar(i,12)
    else
      s_n2c = 0.
    end if

    respcoeff = pftpar(i,5)

    !---
    !daily individual respiration amounts

    lresp_ind = lm_ind(i) * l_n2c * respcoeff * k * gtemp_air  * dphen(:,i)
    rresp_ind = rm_ind(i) * r_n2c * respcoeff * k * gtemp_soil * dphen(:,i)
    sresp_ind = sm_ind(i) * s_n2c * respcoeff * k * gtemp_air

    !---
    !first calculate annual npp and make sure it is positive

    agpp = sum(mgpp(:,i))
    
    aresp_ind = sum(lresp_ind + rresp_ind + sresp_ind)

    anpp(i) = agpp - nind(i) * aresp_ind

    if (anpp(i) < 0.)  then
      
      !make up the negative npp with what is stored in cstore, if possible

      !write(0,*)'negative npp',i,anpp(i),cstore(i)
      
      mgpp(:,i) = mgpp(:,i) + cstore(i) / 12.
      cstore(i) = 0.

      !anpp(i) = cstore(i)
      !mnpp(:,i) = 0.
      
    end if
    
    !disaggregate npp into monthly values

    do m = 1,12

      a = firstday(m)
      b = a + ndaymonth(m) - 1

      mnpp(m,i) = mgpp(m,i) - nind(i) * sum(lresp_ind(a:b) + rresp_ind(a:b) + sresp_ind(a:b))

    end do
    
    gresp(:,i) = 0.25 * mnpp(:,i)
    
    cstore(i) = sum(gresp(:,i)) + 0.5 * cstore(i)
    
    mnpp(:,i) = mnpp(:,i) - gresp(:,i)
    
    anpp(i) = sum(mnpp(:,i))
  
  else

    anpp(i)   = 0.
    mnpp(:,i) = 0.

  end if

end do

end subroutine calcnpp

!----------------------------------------------------------------------------------------------------------------------

end module nppmod
