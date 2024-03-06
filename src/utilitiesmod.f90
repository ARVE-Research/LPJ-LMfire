module utilitiesmod

use parametersmod, only : sp,dp

implicit none

public :: tunit2year
public :: matsol
public :: pos
public :: area

interface pos
  module procedure pos_sp
  module procedure pos_dp
end interface

contains

! ------------------------------------------------------------------------------------------------------------------

integer(i4) function tunit2year(timeunit) result(yr)

! extract the year from a netCDF-compliant time unit string

use parametersmod, only : i4

implicit none

character(*), intent(in) :: timeunit

integer(i4) :: i

if (timeunit(1:11) /= 'days since ') then

  write(0,*)'malformed unit for time in climate input file, must be in format: "days since YYYY-MM-DD HH:MM:SS"'
  stop

end if

i = verify(timeunit,'days since ')

read(timeunit(i:i+3),*)yr

end function tunit2year

! ------------------------------------------------------------------------------------------------------------------

subroutine matsol(mat,sol)

! Provides matrix solution to X matrix in a A * X = B system using LU decomposition
! Code adapted from Press et al. (1996) Numerical Recipes in Fortran 90 : The Art of Parallel Scientific Computing
! 2nd Edition (P.1016-1017)

use parametersmod, only : i4,sp

implicit none

! arguments

real(sp), dimension(:,:), intent(inout) :: mat
real(sp), dimension(:)  , intent(inout) :: sol

! parameter

real(sp), parameter :: tiny_sp = 1.e-38_sp

! local variables

integer(i4), dimension(size(sol))        :: indx
real(sp),    dimension(size(sol))        :: mv
real(sp),    dimension(:,:), allocatable :: prod
integer(i4), dimension(1)                :: maxl

integer(i4) :: i
integer(i4) :: n
integer(i4) :: k
integer(i4) :: ll
integer(i4) :: max
real(sp)    :: summ

! ----------------------------------

n = size(sol)

mv = 1. / maxval(abs(mat),dim=2)

! ---

do i = 1,n

  maxl = maxloc(mv(i:n) * abs(mat(i:n,i)))

  max = (i - 1) + maxl(1)

  indx(i) = max

  ! ---

  if (mat(i,i) == 0.) mat(i,i) = tiny_sp

  mat(i+1:n, i) = mat(i+1:n,i) / mat(i,i)

  ! ---

  allocate(prod(i+1:n,i+1:n))

  prod = spread(mat(i+1:n,i),dim=2,ncopies=size(mat(i,i+1:n)))

  prod = prod * spread(mat(i,i+1:n), dim=1, ncopies=size(mat(i+1:n, i)))

  ! ---

  mat(i+1:n, i+1:n) = mat(i+1:n,i+1:n) - prod

  deallocate(prod)

end do

!---

k = 0

do i = 1, n

  ll = indx(i)
  summ = sol(ll)
  sol(ll) = sol(i)

  if (k /= 0) then

    summ = summ - dot_product(mat(i,k:i-1),sol(k:i-1))

  else if (summ /= 0.) then

    k = i

  end if

  sol(i) = summ

end do

!---

do i = n,1,-1

  sol(i) = (sol(i) - dot_product(mat(i,i+1:n),sol(i+1:n))) / mat(i,i)

end do

end subroutine matsol

! -------------------------------------------------------------------------------

integer function pos_sp(arr,x)

  use parametersmod, only : sp
  
  implicit none
  
  real(sp), dimension(:), intent(in)  :: arr
  real(sp),               intent(in)  :: x
    
  ! ----
  
  pos_sp = minloc(abs(arr - x),dim = 1)

end function pos_sp

! -------------------------------------------------------------------------------

integer function pos_dp(arr,x)

  use parametersmod, only : dp
  
  implicit none
  
  real(dp), dimension(:), intent(in)  :: arr
  real(dp),               intent(in)  :: x
    
  ! ----
  
  pos_dp = minloc(abs(arr - x),dim = 1)

end function pos_dp

! ------------------------------------------------------------------------------------------------------------------

real(sp) function area(lat,minutes)

  ! this function returns the size of a regular grid cell in square meters.

  implicit none

  real(dp), intent(in) :: lat
  real(sp), intent(in), dimension(2) :: minutes

  real(dp), parameter :: pi      =    3.14159265359_dp
  real(dp), parameter :: radius  = 6378.137_dp ! km, WGS-84 spherical approximation
  real(dp), parameter :: deg2rad = pi / 180._dp

  real(dp) :: cellarea
  real(dp) :: deltalat
  real(dp) :: deltalon
  real(dp) :: elevation
  real(dp), dimension(2) :: resolution

  resolution = minutes / 60._dp

  elevation = deg2rad * (lat + 0.5_dp * resolution(2))

  deltalat = deg2rad * resolution(1)
  deltalon = deg2rad * resolution(2)

  cellarea = 2._dp * radius**2 * deltalon * cos(elevation) * sin(0.5_dp * deltalat)

  area = real(cellarea * 1.e6)

end function area

! ------------------------------------------------------------------------------------------------------------------

end module utilitiesmod
