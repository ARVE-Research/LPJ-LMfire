module utilitiesmod

use parametersmod, only : i4,sp

implicit none

public :: matsol
public :: findloc

contains

! ------------------------------------------------------------------------------------------------------------------

subroutine matsol(mat,sol)

! Provides matrix solution to X matrix in a A * X = B system using LU decomposition
! Code adapted from Press et al. (1996) Numerical Recipes in Fortran 90 : The Art of Parallel Scientific Computing
! 2nd Edition (P.1016-1017)

implicit none

real(sp), dimension(:,:), intent(inout) :: mat
real(sp), dimension(:)  , intent(inout) :: sol

integer(i4), dimension(size(sol)) :: indx
real(sp)   , dimension(size(sol)) :: mv
real(sp)   , dimension(:,:), allocatable   :: prod
integer(i4), dimension(1)   :: maxl

real(sp)   , parameter      :: tiny_sp = 1.0e-38_sp

integer(i4)                 :: i, n, k, ll
integer(i4)                 :: max
real(sp)                    :: summ

! ----------------------------------

n = size(sol)

mv = 1. / maxval(abs(mat), dim=2)

!---

do i = 1, n

  maxl = maxloc(mv(i:n) * abs(mat(i:n,i)))

  max = (i - 1) + maxl(1)

  indx(i) = max

  !---

  if (mat(i,i) == 0.) mat(i,i) = tiny_sp

  mat(i+1:n, i) = mat(i+1:n, i) / mat(i,i)

  !---

  allocate(prod(i+1:n, i+1:n))

  prod = spread(mat(i+1:n, i), dim=2, ncopies=size(mat(i, i+1:n)))

  prod = prod * spread(mat(i, i+1:n), dim=1, ncopies=size(mat(i+1:n, i)))

  !---

  mat(i+1:n, i+1:n) = mat(i+1:n, i+1:n) - prod

  deallocate(prod)

end do

!---

k = 0

do i = 1, n

  ll = indx(i)
  summ = sol(ll)
  sol(ll) = sol(i)

  if (k /= 0) then

    summ = summ - dot_product(mat(i, k:i-1), sol(k:i-1))

  else if (summ /= 0.) then

    k = i

  end if

  sol(i) = summ

end do

!---

do i = n, 1, -1

  sol(i) = (sol(i) - dot_product(mat(i, i+1:n), sol(i+1:n))) / mat(i,i)

end do

end subroutine matsol

! -------------------------------------------------------------------------------

subroutine findloc(mat,x,loc)

! returns the index (loc) of a value (x) in a 1D array (mat)
! adapted from Press et al. (1996) Numerical Recipes in Fortran 90 : The Art of Parallel Scientific Computing
! 2nd Edition

implicit none

real(sp), dimension(:), intent(in)  :: mat
real(sp),               intent(in)  :: x
integer(i4),            intent(out) :: loc

real(sp), dimension(:), allocatable :: diff
integer(i4) :: len

!----

len = size(mat)

allocate(diff(len))

diff = abs(mat - x)

loc = minloc(diff,dim=1)

end subroutine findloc

! ------------------------------------------------------------------------------------------------------------------

end module utilitiesmod
