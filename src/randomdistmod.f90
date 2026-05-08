module randomdistmod

! Random number generators and distributions module
!
! This module contains subroutines and functions to compute the gamma and
! normal distribution, as well as a random number generator. Most of the
! subroutines have been translated from other languages to FORTRAN. The main
! references are
!
! - http://www.johndcook.com/SimpleRNG.cpp for the random number generators
! - http://people.sc.fsu.edu/~jburkardt/f_src/prob/prob.html by John
!   Burkardt for the gamma and normal distribution functions

use parametersmod, only : sp,dp,i4

implicit none

public  :: ran_seed
public  :: ranu
public  :: ranur
public  :: ran_normal
public  :: ran_gamma
public  :: ran_gp
public  :: ran_gamma_gp

private :: refill

!----------------
!module variables

integer(i4), parameter :: qsiz  = 10   !41265_i4
integer(i4), parameter :: cmul  = 69609_i4
integer(i4), parameter :: coffs =   123_i4

real(sp), parameter :: rng1 = 1. / (2. * real(huge(i4)))  !scales a random integer to -0.5,0.5
real(sp), parameter :: rng2 = 1. / real(huge(i4))         !scales a random integer to -1,1

real(sp), parameter :: one    = 1.
real(sp), parameter :: half   = 0.5
real(sp), parameter :: vsmall = tiny(1.)
real(sp), parameter :: zero   = 0.

type randomstate

  integer(i4), dimension(qsiz) :: q
  integer(i4) :: carry =       362_i4
  integer(i4) :: xcng  =   1236789_i4
  integer(i4) :: xs    = 521288629_i4  !default seed
  integer(i4) :: indx  = qsiz + 1
  logical     :: have  = .false.

end type randomstate  !5+qsiz elements = 15 elements

contains

!-----------------------------------------------------------------------

subroutine ran_seed(sval,state)

! Set the seed of the random state

implicit none

integer(i4),               intent(in)    :: sval  ! state of the uniform random number generator
type(randomstate), target, intent(inout) :: state ! the random state

integer     :: i

integer(i4), pointer :: xcng
integer(i4), pointer :: xs
integer(i4), pointer, dimension(:) :: q

!---------------------

q    => state%q
xcng => state%xcng
xs   => state%xs

xs = sval

!---

do i = 1,qsiz

  xcng = xcng * cmul + coffs
  xs = ieor(xs,ishft(xs, 13))
  xs = ieor(xs,ishft(xs,-17))
  xs = ieor(xs,ishft(xs, -5))
  q(i) = xcng + xs

end do

end subroutine ran_seed

!-----------------------------------------------------------------------

integer(i4) function ranu(state)

! Generates a uniformly distributed random 4 byte integer with the range (-huge(i4),+huge(i4))
! based on the 32-bit super KISS random number generator by George Marsaglia, published online
! and translated to Fortran 90 by user "mecej4" and Marsaglia, http://forums.silverfrost.com/viewtopic.php?t=1480
! Further modifications to pass the complete state of the generator as an argument by J.O. Kaplan, 2011

implicit none

type(randomstate), target, intent(inout) :: state ! state of the uniform random number generator

integer(i4) :: supr

integer(i4), pointer :: indx
integer(i4), pointer :: xcng
integer(i4), pointer :: xs
integer(i4), pointer, dimension(:) :: q

!---------------------

indx => state%indx
q    => state%q
xcng => state%xcng
xs   => state%xs

!---

if (indx <= qsiz) then
  supr = q(indx)
  indx = indx + 1
else                     !reset the generator
  supr = refill(state)
end if

!---

xcng = xcng * cmul + coffs
xs   = ieor(xs,ishft(xs, 13))
xs   = ieor(xs,ishft(xs,-17))
xs   = ieor(xs,ishft(xs, -5))

ranu = xcng + xs + supr

end function ranu

!-----------------------------------------------------------------------

real(sp) function ranur(state)

! generate a random number in the range (0,1)

implicit none

type(randomstate), target, intent(inout) :: state ! state of the uniform random number generator

!----

ranur = real(ranu(state)) * rng1 + half

end function ranur

!-----------------------------------------------------------------------

function refill(state) result(s)

! reset a random state

implicit none

type(randomstate), target, intent(inout) :: state ! state of the uniform random number generator

integer(i4) :: s
integer(i4) :: z
integer(i4) :: h

integer     :: i

integer(i4), pointer :: indx
integer(i4), pointer :: carry
integer(i4), pointer, dimension(:) :: q

!---------------------

indx  => state%indx
carry => state%carry
q     => state%q

!---

do i = 1,qsiz

  h = iand(carry,1_i4)
  z = ishft(ishft(q(i),9),-1) + ishft(ishft(q(i),7),-1) + ishft(carry,-1)

  carry = ishft(q(i),-23) + ishft(q(i),-25) + ishft(z,-31)

  q(i) = not(ishft(z,1)+h)

end do

indx = 2
s = q(1)

end function refill

!-----------------------------------------------------------------------

subroutine ran_normal(state,nval)

! Sampler for the normal distribution centered at 0 with std. dev. of unity,
! based on Marsaglia polar method

implicit none

type(randomstate), intent(inout) :: state ! state of the uniform random number generator
real(sp),          intent(out)   :: nval  ! output: The random number from the normal distribution

!---

real(sp), dimension(2) :: vals

integer(i4), dimension(2) :: u
real(sp),    dimension(2) :: v

real(sp) :: s
real(sp) :: a

!---------------------

do

  u(1) = ranu(state)
  u(2) = ranu(state)

  v = real(u) * rng2   !convert integer (-huge,+huge) to (-1,+1)

  s = sum(v**2)

  if (s < 1.) exit

end do

a = sqrt(-2. * log(s) / s)

vals = v * a

if (state%have) then

  state%have = .false.
  nval = vals(2)

else

  nval = vals(1)
  state%have = .true.

end if

end subroutine ran_normal

!----------------------------------------------------------------------------------------------------------

real(sp) recursive function ran_gamma(state,first,shape,scale) result(ret)

! Select a random number from a Gamma distribution
!
! adapted from the cpp adaptation of the Marsaglia & Tsang random gamma algorithm in:
! http://www.johndcook.com/SimpleRNG.cpp
!
! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000), *A simple method for generating
! gamma variables*, Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

implicit none

!arguments

type(randomstate), intent(inout) :: state  ! state of the uniform random number generator
logical,           intent(in)    :: first  ! flag if this is the first call to the distribution with this shape
real(sp),          intent(in)    :: shape  ! shape parameter of the Gamma distribution (k or alpha, unitless) > 0
real(sp),          intent(in)    :: scale  ! scale parameter of the Gamma distribution (theta = 1/beta)

!local variables

real(sp), save  :: c
real(sp), save  :: d
real(sp)        :: u
real(sp)        :: v
real(sp)        :: x

!---------------------

if (shape <= 0.) then

  write(0,*) 'shape parameter value must be positive'
  stop

else if (shape >= 1.) then

  if (first) then
    d = shape - 1. / 3.
    c = 1. / sqrt(9. * d)
  end if

  do
    do  !generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

      call ran_normal(state,x)
      v = (1. + c * x)**3
      if (v > 0.) exit

    end do

    u = ranur(state)  ! generate uniform variable u in the range (0,1)

    if (u < 1. - 0.0331 * x**4 .or. log(u) < half * x**2 + d*(1. - v + log(v))) then
      ret = scale * d * v
      exit
    end if
  end do

else

  ret = scale * ran_gamma(state,first,shape + 1.,1.) * ranur(state)**(1. / shape)

end if

end function ran_gamma

!----------------------------------------------------------------------------------------------------------

real(sp) function ran_gamma_gp(state,first,shape,scale,thresh,shape_gp,scale_gp) result(ret)

! Select a random number from a hybrid Gamma-GP distribution

implicit none

!arguments

type(randomstate), intent(inout) :: state     ! state of the uniform random number generator
logical,           intent(in)    :: first     ! flag if this is the first call to the distribution with this shape
real(sp),          intent(in)    :: shape     ! shape parameter of the Gamma distribution (k or alpha, unitless) > 0
real(sp),          intent(in)    :: scale     ! scale parameter of the Gamma distribution (theta = 1/beta)
real(sp),          intent(in)    :: thresh    ! the threshold above which to choose the GP distribution
real(sp),          intent(in)    :: shape_gp  ! shape parameter of the GP distribution
real(sp),          intent(in)    :: scale_gp  ! scale parameter of the GP distribution

!---------------------------------------

ret = ran_gamma(state,first,shape,scale)

if (ret > thresh) then
  ret = ran_gp(state,shape_gp,scale_gp,thresh)
endif

end function ran_gamma_gp

!-------------------------------

real(sp) function ran_gp(state,shape,scale,loc)

! Select a random number from a generalized pareto (GP) distribution

implicit none

type(randomstate), intent(inout) :: state  ! state of the uniform random number generator
real(sp),          intent(in)    :: shape  ! shape parameter of the GP distribution (k or alpha, unitless) > 0
real(sp),          intent(in)    :: scale  ! scale parameter of the GP distribution (theta = 1/beta)
real(sp),          intent(in)    :: loc    ! the location of the GP distribution

real(sp)        :: u

!---------------------------------------

u = ranur(state)  ! generate uniform variable u in the range (0,1)

if (shape == 0.0_sp) then
  ran_gp = loc - scale * log(u)
else
  ran_gp = loc + scale * (U**(-shape) - 1) / shape
endif

end function ran_gp

!----------------------------------------------------------------------------------------------------------

end module randomdistmod
