module randomdistmod

implicit none

public  :: ranur
public  :: ranu
public  :: ran_seed
public  :: ran_normal
public  :: ran_gamma

private :: refill
private :: random_gamma1
private :: random_gamma2
private :: random_exponential

!----------------
!module variables

integer, parameter :: sp = selected_real_kind(4)      !4 byte real
integer, parameter :: i4 = selected_int_kind(9)       !10^9 fits into 4 bytes

integer(i4), parameter :: qsiz  =    10 !41265_i4
integer(i4), parameter :: cmul  = 69609_i4
integer(i4), parameter :: coffs =   123_i4

real(sp), parameter :: rng1 = 1. / (2. * real(huge(i4)))  !scales the random integer to -0.5,0.5
real(sp), parameter :: rng2 = 1. / real(huge(i4))         !scales the random integer to -1,1

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

real(sp) function ranur(state)

!generate a single precision random real on the interval [0,1]

implicit none

type(randomstate), target, intent(inout) :: state

!-----

ranur = rng1 * real(ranu(state)) + half

end function ranur

!-----------------------------------------------------------------------

integer(i4) function ranu(state)
  
!Generates a uniformly distributed random 4 byte integer with the range (-huge(i4),+huge(i4))
!based on the 32-bit super KISS random number generator by George Marsaglia, published online
!and translated to Fortran 90 by user "mecej4" and Marsaglia, http://forums.silverfrost.com/viewtopic.php?t=1480
!Further modifications to pass the complete state of the generator as an argument by J.O. Kaplan, 2011

implicit none

type(randomstate), target, intent(inout) :: state

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

if (indx <= qsiz) then   !reset the generator
  supr = q(indx)
  indx = indx + 1
else 
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

function refill(state) result(s)

implicit none

type(randomstate), target, intent(inout) :: state

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

subroutine ran_seed(sval,state)

implicit none

integer(i4),               intent(in)    :: sval
type(randomstate), target, intent(inout) :: state

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

subroutine ran_normal(state,nval)

!Sampler for the normal distribution based on Marsaglia polar method

implicit none

type(randomstate), intent(inout) :: state
real(sp),          intent(out)   :: nval

!---

real(sp), dimension(2) :: vals

integer(i4), dimension(2) :: u
real(sp),    dimension(2) :: v

real(sp) :: s
real(sp) :: a

!---------------------

if (state%have) then

  state%have = .false.
  nval = vals(2)

else
  
  do

    u(1) = ranu(state)
    u(2) = ranu(state)

    v = real(u) * rng2   !convert integer (-huge,+huge) to (-1,+1)

    s = sum(v**2)

    if (s < 1.) exit

  end do

  a = sqrt(-2. * log(s) / s)

  vals = v * a
  
  nval = vals(1)

  state%have = .true.

end if

end subroutine ran_normal

!----------------------------------------------------------------------------------------------------------

!The functions below are for randomly sampling the gamma distribution
  
!adapted from code by Alan J. Miller, URL,
!http://users.bigpond.net.au/amiller/random.html 

!----------------------------------------------------------------------------------------------------------

subroutine ran_gamma(state,s,b,first,fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     function generates a random gamma variate.
!     calls either random_gamma1 (s > 1.0)
!     or random_exponential (s = 1.0)
!     or random_gamma2 (s < 1.0).

!     s = shape parameter of distribution (0 < real)
!     b = scale parameter

implicit none

!arguments

type(randomstate), intent(inout) :: state
real(sp),          intent(in)    :: s       !shape parameter of the Gamma distribution (alpha, unitless)
real(sp),          intent(in)    :: b       !scale parameter of the Gamma distribution (Beta)
logical,           intent(in)    :: first   !flag if this is the first call to the distribution 
real(sp),          intent(out)   :: fn_val

!--------

if (s <= zero) then
  write(0, *) 'shape parameter value must be positive'
  stop
end if

if (s > one) then
  call random_gamma1(state,s,first,fn_val)
else if (s < one) then
  call random_gamma2(state,s,first,fn_val)
else
  call random_exponential(state,fn_val)
end if

!scale the random variable with Beta

fn_val = b * fn_val

end subroutine ran_gamma

!--------------------------------------------------------------------

subroutine random_gamma1(state,s,first,fn_val)

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s > 1.
  
implicit none

!arguments

type(randomstate), intent(inout) :: state
real(sp),          intent(in)    :: s
logical,           intent(in)    :: first
real(sp),          intent(out)   :: fn_val

!local variables

real(sp)   :: c
real(sp)   :: d
real(sp)        :: u
real(sp)        :: v
real(sp)        :: x

!--------

if (first) then
  d = s - one / 3.
  c = one / sqrt(9. * d)
end if

! start of main loop
do

  ! generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.
  do
    call ran_normal(state,x)
    v = (one + c * x)**3
    if (v > zero) exit
  end do

  ! generate uniform variable u in the range (0,1)

  u = real(ranu(state)) * rng1 + half

  if (u < one - 0.0331 * x**4) then
    fn_val = d*v
    exit
  else if (log(u) < half * x**2 + d*(one - v + log(v))) then
    fn_val = d*v
    exit
  end if
end do

end subroutine random_gamma1

!--------------------------------------------------------------------

subroutine random_gamma2(state,s,first,fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function generates a random variate in [0,infinity) from
! a gamma distribution with density proportional to
! gamma2**(s-1) * exp(-gamma2),
! using a switching method.

!    s = shape parameter of distribution
!          (real < 1.0)

implicit none

!arguments

type(randomstate), intent(inout) :: state
real(sp),          intent(in)    :: s
logical,           intent(in)    :: first
real(sp),          intent(out)   :: fn_val

!local variables

real(sp)       :: r
real(sp)       :: x
real(sp)       :: w

real(sp)  :: a
real(sp)  :: p
real(sp)  :: c
real(sp)  :: uf
real(sp)  :: vr
real(sp)  :: d

!--------

if (s <= zero .or. s >= one) then
  write(0, *) 'shape parameter value outside permitted range'
  stop
end if

if (first) then                        ! initialization, if necessary
  a = one - s
  p = a / (a + s * exp(-a))
  if (s < vsmall) then
    write(0,*) 'shape parameter value too small'
    stop
  end if
  c = one / s
  uf = p * (vsmall / a)**s
  vr = one - vsmall
  d = a * log(a)
end if

do
  r = real(ranu(state)) * rng1 + half !ranu(state)  !value between (0,1)

  if (r >= vr) then
    cycle
  else if (r > p) then
    x = a - log((one - r)/(one - p))
    w = a * log(x) - d
  else if (r > uf) then
    x = a * (r / p)**c
    w = x
  else
    fn_val = zero
    return
  end if

  r = real(ranu(state)) * rng1 + half !ranu(state)

  if (one - r <= w .and. r > zero) then
    if (r * (w + one) >= one) cycle
    if (-log(r) <= w) cycle
  end if

  exit

end do

fn_val = x

end subroutine random_gamma2

!--------------------------------------------------------------------

subroutine random_exponential(state,fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function generates a random variate in [0,infinity) from
! a negative exponential distribution wlth density proportional
! to exp(-random_exponential), using inversion.

implicit none

!arguments

type(randomstate), intent(inout) :: state
real(sp),          intent(out)   :: fn_val

!local variable

real(sp) :: r

!--------

do
  r = real(ranu(state)) * rng1 + half !ranu(state)
  if (r > zero) exit
end do

fn_val = -log(r)

end subroutine random_exponential

!-------------------------------

end module randomdistmod
