module bioclimmod

implicit none

public :: climate20
public :: bioclim

contains

!----------------------------------------------------------

subroutine climate20(mtemp,dtemp,gdd,mtemp_min_buf,gdd_buf,mtemp_min20,gdd20,mtemp_max,mat20,mat_buf)

use parametersmod, only : sp

implicit none

!parameters
real(sp), parameter :: gddbase = 5.
integer,  parameter :: climbuf = 20

!arguments

real(sp), intent(in), dimension(:) :: mtemp
real(sp), intent(in), dimension(:) :: dtemp

!arguments used locally by other core routines

real(sp), intent(out) :: mat20
real(sp), intent(out) :: gdd
real(sp), intent(out) :: gdd20
real(sp), intent(out) :: mtemp_max
real(sp), intent(out) :: mtemp_min20

!state variables that must be saved from one year to the next

real(sp), intent(inout), dimension(:) :: mat_buf         !
real(sp), intent(inout), dimension(:) :: gdd_buf
real(sp), intent(inout), dimension(:) :: mtemp_min_buf

!local variables
real(sp) :: mtemp_min
real(sp) :: yrs
real(sp) :: mat

!---------
!mean annual temperature

mat = sum(mtemp) / 12.

!write(0,'(13f7.2)')mtemp,mat

mat_buf = eoshift(mat_buf,-1,mat)
  
mat20 = sum(mat_buf,mask=mat_buf /= -9999.) / count(mat_buf /= -9999.)

!mean temperature of the coldest month

mtemp_min = minval(mtemp)

mtemp_min_buf = eoshift(mtemp_min_buf,-1,mtemp_min)

mtemp_min20 = sum(mtemp_min_buf,mask=mtemp_min_buf /= -9999.) / count(mtemp_min_buf /= -9999.)

!---------
!MTWM

mtemp_max = maxval(mtemp)

!---------
!calculate GDD5 and add to buffer

gdd = sum(dtemp-gddbase,mask = dtemp>=gddbase)

gdd_buf = eoshift(gdd_buf,-1,gdd)

gdd20 = sum(gdd_buf,mask=gdd_buf /= -9999.) / count(gdd_buf /= -9999.)

!write(0,'(a,22f8.1)')'bioclim',gdd,gdd20,gdd_buf

end subroutine climate20

!----------------------------------------------------------

subroutine bioclim(mtemp_min20,gdd,mtemp_max,survive,estab)

!Limits based on 20-year running averages of coldest-month mean temperature and growing degree days (5 degree base).
!For SURVIVAL, coldest month temperature and GDD should be at least as high as PFT-specific limits.
!For REGENERATION, PFT must be able to survive AND coldest month temperature should be no higher than a PFT-specific limit.

use parametersmod, only : sp,pftpar

implicit none

!arguments

real(sp), intent(in) :: mtemp_min20
real(sp), intent(in) :: gdd
real(sp), intent(in) :: mtemp_max

logical, dimension(:), intent(out) :: survive
logical, dimension(:), intent(out) :: estab

integer, parameter :: pfts = size(pftpar,dim=1)

!integer :: i

!local variables

real(sp), dimension(pfts) :: tcmin   !PFT-specific minimum coldest-month temperature
real(sp), dimension(pfts) :: tcmax   !PFT-specific maximum coldest-month temperature
real(sp), dimension(pfts) :: gddmin  !PFT-specific minimum GDD
real(sp), dimension(pfts) :: twmax   !PFT-specific upper limit of warmest-month temperature

!------------------------------------------

tcmin  = pftpar(:,27)
tcmax  = pftpar(:,28)
gddmin = pftpar(:,29)
twmax  = pftpar(:,30)

where (mtemp_min20 >= tcmin)

  survive = .true.

  where (gdd >= gddmin .and. mtemp_min20 <= tcmax .and. mtemp_max <= twmax)

    estab = .true.

  elsewhere

    estab = .false.

  end where

elsewhere

  survive = .false.

end where

!do i = 1,pfts
!  write(0,*)i,tcmin(i),tcmax(i),gddmin(i),twmax(i),gdd,mtemp_min20,mtemp_max,survive(i),estab(i)
!end do

end subroutine bioclim

!----------------------------------------------------------

end module bioclimmod
