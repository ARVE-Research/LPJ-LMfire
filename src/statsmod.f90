module statsmod

! The following functions and subroutines are taken from R, version 3.3.2 and
! converted from the original src/nmath/qgamma.c to fortran by Philipp Sommer,
! 2016, November 26th
! updates for legibility by Jed Kaplan, summer 2019

! Functions and subroutines to calculate the gamma distribution
! function (gamma_cdf) and the gamma density (gamma_pdf) are taken from
! John Burkardt, Department of Scientific Computing at Florida State University (FSU)
! http://people.sc.fsu.edu/~jburkardt/f_src/prob/prob.html
! on February 23rd, 2016. Last revised on 21 August 2013.

! May 2020: 
!   Removed the hardcoded functions gamma and gamma_log in favor of the
!   (hopefully faster) Fortran 2008 intrinsics gamma and log_gamma.
!   Testing change calculations to single precision since that is all we need really.

implicit none

public :: normal_cdf_inv        ! inverse of the normal CDF
public :: gamma_pdf             ! Gamma probability density function (PDF)
public :: gamma_cdf             ! Gamma cumulative distribution function (CDF) 
public :: gamma_cdf_inv         ! inverse of the Gamma CDF

private :: normal_01_cdf_inv    ! one side of the inverse normal CDF, used by normal_01_cdf_inv
private :: evalpoly             ! function to evaluate a polynomial function of any order, used by normal_01_cdf_inv
private :: qchisq_appr          ! chi-square approximation, used by gamma_cdf_inv
private :: gamma_inc            ! incomplete Gamma function, used by gamma_cdf
private :: normal_01_cdf        ! evaluate the Normal 01 CDF, used by gamma_inc

interface gamma_cdf
  module procedure gamma_cdf_sp
  module procedure gamma_cdf_dp
end interface gamma_cdf

contains

! ------------------------------------------------------------------------------
! public routines
! ------------------------------------------------------------------------------

subroutine normal_cdf_inv(cdf,a,b,x)

! Invert the Normal CDF.

! Licensing:
!    This code is distributed under the GNU LGPL license.
! Modified:
!    - 23 February 1999
!    - Extracted: November, 2016
! Author:
!    - John Burkardt
!    - Extracted by Philipp Sommer

use parametersmod, only : dp

implicit none

! arguments

real(dp), intent(in)  :: cdf ! the value of the CDF. 0. <= CDF <= 1.
real(dp), intent(in)  :: a   ! the mean of the pdf
real(dp), intent(in)  :: b   ! the standard deviation of the pdf
real(dp), intent(out) :: x   ! the corresponding argument

! local variables

real(dp) :: x2

! ----

if (cdf < 0._dp .or. 1._dp < cdf) then
  write (0,'(a)')''
  write (0,'(a)')'NORMAL_CDF_INV - Fatal error!'
  write (0,'(a)')'CDF < 0 or 1 < CDF.'
  stop
end if

call normal_01_cdf_inv(cdf,x2)

x = a + b * x2

end subroutine normal_cdf_inv

! ------------------------------------------------------------------------------

subroutine gamma_pdf(x,a,b,c,pdf)

! Evaluate the Gamma PDF.

! .. math::

!    PDF(a,b,c;x) = \exp({-(x-a)/b}) \cdot ((x-a)/b)^{c-1} / (b \cdot \Gamma(c))

! - GAMMA_PDF(A,B,C;X), where C is an integer, is the Erlang PDF.
! - GAMMA_PDF(A,B,1;X) is the Exponential PDF.
! - GAMMA_PDF(0,2,C/2;X) is the Chi Squared PDF with C degrees of freedom.

! Licensing:
!    This code is distributed under the GNU LGPL license.
! Modified:
!    - 02 January 2000
!    - Extracted: June, 2016
! Author:
!    - John Burkardt
!    - Extracted by Philipp Sommer

use parametersmod, only : sp,dp

implicit none

! arguments

real(sp), intent(in)  :: x    ! the argument of the PDF. A <= X
real(sp), intent(in)  :: a    ! the location of the peak;  A is often chosen to be 0.
real(sp), intent(in)  :: b    ! the "scale" parameter; 0. < B, and is often 1.
real(sp), intent(in)  :: c    ! the "shape" parameter; 0. < C, and is often 1.
real(dp), intent(out) :: pdf  ! the returned value of the PDF

! local variable

real(dp) :: y
real(dp) :: c1

! ---------------------------

c1 = c

if (x <= a) then

  pdf = 0._dp

else

  y = (x - a) / b

  pdf = y**(c - 1.) / (b * gamma(c1) * exp(y))

end if

end subroutine gamma_pdf

! ------------------------------------------------------------------------------

subroutine gamma_cdf_sp(x,a,b,c,cdf)

! Evaluate the Gamma CDF.

! Licensing:
!   This code is distributed under the GNU LGPL license.
! Modified:
!   02 January 2000
!   Extracted: June, 2016
! Author:
!   John Burkardt
!   Extracted by Philipp Sommer

use parametersmod, only : sp,dp

implicit none

! arguments

real(sp), intent(in)  :: x   ! the input value for which to compute the CDF
real(sp), intent(in)  :: a   ! the location (< `x`) of the gamma distribution (usually 0)
real(sp), intent(in)  :: b   ! the shape (> 0.) of the distribution
real(sp), intent(in)  :: c   ! the scale (>0.) of the distribution
real(dp), intent(out) :: cdf ! the returned value of the CDF

! local variables

real(dp) :: p2
real(dp) :: x2

! ---------------

x2 = (x - a) / b
p2 = c

cdf = gamma_inc(p2,x2)

end subroutine gamma_cdf_sp

! ------------------------------------------------------------------------------

subroutine gamma_cdf_dp(x,a,b,c,cdf)

! Evaluate the Gamma CDF.

! Licensing:
!   This code is distributed under the GNU LGPL license.
! Modified:
!   02 January 2000
!   Extracted: June, 2016
! Author:
!   John Burkardt
!   Extracted by Philipp Sommer

use parametersmod, only : dp

implicit none

! arguments

real(dp), intent(in)  :: x   ! the input value for which to compute the CDF
real(dp), intent(in)  :: a   ! the location (< `x`) of the gamma distribution (usually 0)
real(dp), intent(in)  :: b   ! the shape (> 0.) of the distribution
real(dp), intent(in)  :: c   ! the scale (>0.) of the distribution
real(dp), intent(out) :: cdf ! the returned value of the CDF

! local variables

real(dp) :: p2
real(dp) :: x2

! ---------------

x2 = (x - a) / b
p2 = c

cdf = gamma_inc(p2,x2)

end subroutine gamma_cdf_dp

! ------------------------------------------------------------------------------

real(dp) function gamma_cdf_inv(p,alpha,scale)

! Compute the quantile function of the gamma distribution.

! This function is based on the Applied Statistics Algorithm AS 91
! ("ppchi2") and via pgamma(.) AS 239.

! References
!     Best, D. J. and D. E. Roberts (1975).
!     Percentage Points of the Chi-Squared Distribution.
!     Applied Statistics 24, page 385.

! .. note::

!    Compared to the original R function, we do not use the final
!    newton step which might lead to values going to infinity for
!    quantiles close to 1

use parametersmod, only : dp

implicit none

! arguments

real(dp), intent(in) :: p     ! the quantile between 0 and 1
real(dp), intent(in) :: alpha ! the shape of the gamma distribution
real(dp), intent(in) :: scale ! the scale of the gamma distribution

! parameters

integer,  parameter :: maxit   = 1000
real(dp), parameter :: eps1    = 1.e-2 
real(dp), parameter :: eps2    = 5.e-7
real(dp), parameter :: pmin    = 1.e-25
real(dp), parameter :: pmax    = 1. - 1.e-14
real(dp), parameter :: pneginf = -10.e34

! local variables

integer  :: i
real(dp) :: a
real(dp) :: b
real(dp) :: c 
real(dp) :: g 
real(dp) :: ch 
real(dp) :: ch0
real(dp) :: p1
real(dp) :: p2
real(dp) :: q
real(dp) :: s1 
real(dp) :: s2 
real(dp) :: s3
real(dp) :: s4
real(dp) :: s5
real(dp) :: s6
real(dp) :: t

! ---------------

if (alpha == 0) then

  gamma_cdf_inv = 0
  return
  
end if

g = log_gamma(alpha)

! ----- Phase I : Starting Approximation

ch = qchisq_appr(p,2.*alpha,g,eps1)

if ((ch < eps2) .or. (p > pmax) .or. (p < pmin)) return

! ----- Phase II: Iteration
! Call pgamma() [AS 239]  and calculate seven term taylor series

c = alpha - 1.

s6 = (120. + c * (346. + 127. * c)) / 5040.

ch0 = ch  ! save initial approx.

do i = 1,maxit

  q  = ch
  p1 = 0.5 * ch

  call gamma_cdf(p1*scale,0._dp,scale,alpha,p2)

  p2 = p - p2

  if (ch <= 0.) then
    ch = ch0
    exit
  end if

  if ((p2 < pNEGinf) .or. (ch <= 0)) then
    ch = ch0
    exit
  end if

  t = p2*exp(alpha * log(2.) + g + p1 - c * log(ch))

  b = t / ch

  a = 0.5 * t - b * c

  s1 = (210. + a * (140. + a * (105. + a * (84. + a * (70. + 60. * a))))) / 420.

  s2 = (420. +  a * (735. + a * (966. + a * (1141. + 1278. * a)))) / 2520.
  
  s3 = (210. + a * (462. + a * (707. + 932. * a))) / 2520.
 
  s4 = (252. + a * (672. + 1182. * a) + c * (294. +a * (889. + 1740. * a))) / 5040.

  s5 = (84. + 2264. * a + c*(1175. + 606. * a)) / 2520.

  ch = ch +  t * (1. + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))))

  if (abs(q - ch) < EPS2 * ch) exit

  if (abs(q - ch) > 0.1 * ch) then

    if (ch < q) then
      ch = 0.9 * q
    else
      ch = 1.1 * q
    end if

  end if

end do

gamma_cdf_inv = 0.5  * scale * ch

end function gamma_cdf_inv

! ------------------------------------------------------------------------------
! private (internal) routines
! ------------------------------------------------------------------------------

subroutine normal_01_cdf_inv (p,x)

! Invert the standard normal CDF.

! Licensing:
!    This code is distributed under the GNU LGPL license.
! Modified:
!   - 05 June 2007
!   - Extracted: November, 2016
! Author:
!    - Original FORTRAN77 version by Michael Wichura.
!    - FORTRAN90 version by John Burkardt.
!    - Extracted by Philipp Sommer
! Reference:
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
! .. note::
!
!    The result is accurate to about 1 part in 10^16.

use parametersmod, only : dp

implicit none

! arguments

real(dp), intent(in)  :: p ! the value of the cumulative probability densitity function.  0 < P < 1.
                           ! If P is outside this range, an "infinite" value will be returned.
                           
real(dp), intent(out) :: x ! the normal deviate value with the property that the probability of a
                           ! standard normal deviate being less than or equal to the value is P.

! parameters

real(dp), parameter :: const1 = 0.180625
real(dp), parameter :: const2 = 1.6

real(dp), parameter, dimension(8) :: a = &
  [ 3.3871328727963666080,     &
    1.3314166789178437745e2,   &
    1.9715909503065514427e3,   &
    1.3731693765509461125e4,   &
    4.5921953931549871457e4,   &
    6.7265770927008700853e4,   &
    3.3430575583588128105e4,   &
    2.5090809287301226727e3    ]
            
real(dp), parameter, dimension(8) :: b = &
  [ 1.,                        &
    4.2313330701600911252e1,   &
    6.8718700749205790830e2,   &
    5.3941960214247511077e3,   &
    2.1213794301586595867e4,   &
    3.9307895800092710610e4,   &
    2.8729085735721942674e4,   &
    5.2264952788528545610e3    ]

real(dp), parameter, dimension(8) :: c = &
  [ 1.42343711074968357734,    &
    4.63033784615654529590,    &
    5.76949722146069140550,    &
    3.64784832476320460504,    &
    1.27045825245236838258,    &
    2.41780725177450611770e-1, &
    2.27238449892691845833e-2, &
    7.74545014278341407640e-4  ]

real(dp), parameter, dimension(8) :: d = &
  [ 1.,                        &
    2.05319162663775882187,    &
    1.67638483018380384940,    &
    6.89767334985100004550e-1, &
    1.48103976427480074590e-1, &
    1.51986665636164571966e-2, &
    5.47593808499534494600e-4, &
    1.05075007164441684324e-9  ]

real(dp), parameter, dimension(8) :: e = &
  [ 6.65790464350110377720,    &
    5.46378491116411436990,    &
    1.78482653991729133580,    &
    2.96560571828504891230e-1, &
    2.65321895265761230930e-2, &
    1.24266094738807843860e-3, &
    2.71155556874348757815e-5, &
    2.01033439929228813265e-7  ]

real(dp), parameter, dimension(8) :: f = &
  [ 1.,                        &
    5.99832206555887937690e-1, &
    1.36929880922735805310e-1, &
    1.48753612908506148525e-2, &
    7.86869131145613259100e-4, &
    1.84631831751005468180e-5, &
    1.42151175831644588870e-7, &
    2.04426310338993978564e-15 ]

real(dp), parameter :: split1 = 0.425
real(dp), parameter :: split2 = 5.

! local variables

real(dp) :: q
real(dp) :: r

! -----------------------------------------------

if (p <= 0.) then
  x = -huge(x)
  return
end if

if (1. <= p) then
  x = huge(x)
  return
end if

q = p - 0.5

if (abs(q) <= split1) then

  r = const1 - q**2
  x = q * evalpoly(a,r) / evalpoly(b,r)

else

  if (q < 0.) then
    r = p
  else
    r = 1. - p
  end if

  if (r <= 0.) then

    x = huge(x)

  else

    r = sqrt(-log(r))

    if (r <= split2) then

      r = r - const2
      x = evalpoly(c,r) / evalpoly(d,r)

    else

      r = r - split2
      x = evalpoly(e,r) / evalpoly(f,r)

    end if

  end if

  if (q < 0.) then
    x = -x
  end if

end if

end subroutine normal_01_cdf_inv

! ------------------------------------------------------------------------------

real(dp) function evalpoly(coeff,x)

use parametersmod, only : dp

implicit none

real(dp), dimension(:), intent(in) :: coeff
real(dp),               intent(in) :: x

integer, allocatable, dimension(:) :: exponent
integer :: i
integer :: nc

! ----------------

nc = size(coeff)

allocate(exponent(nc))

exponent = [(i,i=0,nc-1)]

evalpoly = sum(coeff * x**exponent)

end function evalpoly

! ------------------------------------------------------------------------------

real(dp) function qchisq_appr(p,nu,g,tol)

! chi-square approximation of the "gamma_cdf_inv" function

use parametersmod, only : dp

implicit none

! arguments

real(dp), intent(in) :: p    ! the quantile
real(dp), intent(in) :: nu   ! twice the gamma shape
real(dp), intent(in) :: g    ! the logarithm of the gamma function at the gamma shape
real(dp), intent(in) :: tol  ! the tolerance for the approximation

! local variables

real(dp) :: alpha
real(dp) :: a
real(dp) :: c
real(dp) :: ch
real(dp) :: p1
real(dp) :: p2
real(dp) :: q
real(dp) :: t
real(dp) :: x
real(dp) :: lgam1pa

! ---------------------------------------

alpha = 0.5 * nu
c     = alpha - 1.

p1 = log(p)

if (nu < (-1.24) * p1) then

  ! for small chi-squared */
  !   log(alpha) + g = 
  !   log(alpha) + log(gamma(alpha)) =
  !   log(alpha*gamma(alpha)) =
  !   lgamma(alpha+1)
  ! suffers from catastrophic cancellation when alpha << 1

  if (alpha < 0.5) then
    lgam1pa = log_gamma(alpha + 1._dp)
  else
    lgam1pa = log(alpha) + g
  end if

  ch = exp((lgam1pa + p1)/alpha + log(2._dp))

else if (nu > 0.32) then ! using Wilson and Hilferty estimate

  call normal_cdf_inv(p,0._dp,1._dp,x)

  p1 = 2. / (9. * nu)
  ch = nu * (x * sqrt(p1) + 1. - p1) ** 3

  ! approximation for p tending to 1:

  if (ch > 2.2 * nu + 6) ch = -2. * (log(1 - p) - c * log(0.5 * ch) + g)

else

  ch = 0.4
  a = log(1 - p) + g + c * log(2.)

  do while (abs(q - ch) > tol * abs(ch))

    q = ch
    p1 = 1. / (1 + ch * (4.67 + ch))
    p2 = ch * (6.73 + ch * (6.66 + ch))
    t = -0.5 + (4.67 + 2 * ch) * p1 - (6.73 + ch*(13.32 + 3 * ch)) / p2
    ch = ch - (1 - exp(a + 0.5 * ch) * p2 * p1) / t

  end do

end if

qchisq_appr = ch

end function qchisq_appr

! ------------------------------------------------------------------------------

real(dp) function gamma_inc(p,x)

! Compute the incomplete Gamma function.
!
! Formulas:
!
!    .. math::
!
!        \Gamma_{inc}(P, 0) = 0
!
!    .. math::
!
!        \Gamma_{inc}(P, \infty) = 1.
!
!    .. math::
!
!        \Gamma_{inc}(P,X) = \int_0^x{T^{P-1} \exp{(-T)} \mathrm{d}t} / \Gamma(P)
!
! Licensing:
!    This code is distributed under the GNU LGPL license.
! Modified:
!    - 01 May 2001
!    - Extracted: June, 2016
! Author:
!    - Original FORTRAN77 version by B L Shea.
!    - FORTRAN90 version by John Burkardt
!    - Extracted by Philipp Sommer
! Reference:
!   BL Shea,
!   Chi-squared and Incomplete Gamma Integral,
!   Algorithm AS239,
!   Applied Statistics,
!   Volume 37, Number 3, 1988, pages 466-473.

use parametersmod, only : dp

implicit none

real(dp), intent(in) :: p  ! the exponent parameter (0. < P)
real(dp), intent(in) :: x  ! the integral limit parameter. If X is less than or equal to 0, GAMMA_INC is returned as 0.

real(dp), parameter :: exp_arg_min =  -88._dp
real(dp), parameter :: overflow    =    1.e37
real(dp), parameter :: plimit      = 1000._dp
real(dp), parameter :: tol         =    1.e-7
real(dp), parameter :: xbig        =    1.e8

real(dp) :: a
real(dp) :: arg
real(dp) :: b
real(dp) :: c
real(dp) :: cdf
real(dp) :: pn1
real(dp) :: pn2
real(dp) :: pn3
real(dp) :: pn4
real(dp) :: pn5
real(dp) :: pn6
real(dp) :: rn

! ----------------------------------------

gamma_inc = 0._dp

if (p <= 0._dp) then

  write(*,'(a)') ' '
  write(*,'(a)') 'GAMMA_INC - Fatal error!'
  write(*,'(a)') '  Parameter P <= 0.'
  stop

end if

if (x <= 0._dp) then

  gamma_inc = 0._dp
  return

end if

! Use a normal approximation if PLIMIT < P.

if (plimit < p) then

  pn1 = 3._dp * sqrt(p) * ((x / p)**(1._dp/3._dp) + 1._dp / (9._dp * p ) - 1._dp)

  call normal_01_cdf(pn1,cdf)

  gamma_inc = cdf

  return

end if

! Is X extremely large compared to P?

if (xbig < x) then

  gamma_inc = 1._dp
  return

end if

! Use Pearson's series expansion.
! (P is not large enough to force overflow in the log of Gamma.

if (x <= 1._dp .or. x < p) then

  arg = p * log(x) - x - log_gamma(p + 1._dp)
  c = 1._dp
  gamma_inc = 1._dp
  a = p

  do

    a = a + 1._dp
    c = c * x / a
    gamma_inc = gamma_inc + c

    if (c <= tol) exit

  end do

  arg = arg + log(gamma_inc)

  if ( exp_arg_min <= arg ) then
    gamma_inc = exp(arg)
  else
    gamma_inc = 0._dp
  end if

else

  ! Use a continued fraction expansion.

  arg = p * log(x) - x - log_gamma(p)

  a = 1._dp - p
  b = a + x + 1._dp
  c = 0._dp
  
  pn1 = 1._dp
  pn2 = x
  pn3 = x + 1._dp
  pn4 = x * b
  gamma_inc = pn3 / pn4

  do

    a = a + 1._dp
    b = b + 2._dp
    c = c + 1._dp
    pn5 = b * pn3 - a * c * pn1
    pn6 = b * pn4 - a * c * pn2

    if (0._dp < abs (pn6)) then

      rn = pn5 / pn6

      if (abs(gamma_inc - rn) <= min (tol,tol * rn)) then

        arg = arg + log(gamma_inc)

        if ( exp_arg_min <= arg ) then
          gamma_inc = 1._dp - exp(arg)
        else
          gamma_inc = 1._dp
        end if

      return

      end if

      gamma_inc = rn

    end if

    pn1 = pn3
    pn2 = pn4
    pn3 = pn5
    pn4 = pn6

    ! Rescale terms in continued fraction if terms are large.

    if (overflow <= abs(pn5)) then

      pn1 = pn1 / overflow
      pn2 = pn2 / overflow
      pn3 = pn3 / overflow
      pn4 = pn4 / overflow

    end if

  end do
end if

end function gamma_inc

! ------------------------------------------------------------------------------

subroutine normal_01_cdf(x,cdf)

! evaluate the Normal 01 CDF.

! Licensing:
!    This code is distributed under the GNU LGPL license.
! Modified:
!   - 10 February 1999
!   - Extracted: June, 2016
! Author:
!    - John Burkardt
!    - Extracted by Philipp Sommer
! Reference:
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.

use parametersmod, only : dp

implicit none

! arguments

real(dp), intent(in)  :: x   ! the argument of the CDF.
real(dp), intent(out) :: cdf ! the value of the CDF.

! parameters

real(dp), parameter :: a1  =  0.398942280444
real(dp), parameter :: a2  =  0.399903438504
real(dp), parameter :: a3  =  5.75885480458
real(dp), parameter :: a4  = 29.8213557808
real(dp), parameter :: a5  =  2.62433121679
real(dp), parameter :: a6  = 48.6959930692
real(dp), parameter :: a7  =  5.92885724438
real(dp), parameter :: b0  =  0.398942280385
real(dp), parameter :: b1  =  3.8052e-8
real(dp), parameter :: b2  =  1.00000615302
real(dp), parameter :: b3  =  3.98064794e-4
real(dp), parameter :: b4  =  1.98615381364
real(dp), parameter :: b5  =  0.151679116635
real(dp), parameter :: b6  =  5.29330324926
real(dp), parameter :: b7  =  4.8385912808
real(dp), parameter :: b8  = 15.1508972451
real(dp), parameter :: b9  =  0.742380924027
real(dp), parameter :: b10 = 30.789933034
real(dp), parameter :: b11 =  3.99019417011

! local variables

real(dp) :: q
real(dp) :: y

! -------------------------------------------------

if (abs(x) <= 1.28) then

  y = 0.5 * x**2

  q = 0.5 - abs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))))


else if ( abs ( x ) <= 12.7 ) then

  y = 0.5 * x**2

  q = exp(-y) * b0 / (abs(x) - b1 + b2 / (abs(x) + b3 &
              + b4 / (abs(x) - b5 + b6 / (abs(x) + b7 &
              - b8 / (abs(x) + b9 + b10 / (abs(x) + b11))))))

else     ! 12.7 < |X|

  q = 0._dp

end if

! Take account of negative X.

if (x < 0._dp) then
  cdf = q
else
  cdf = 1._dp - q
end if

end subroutine normal_01_cdf

! ------------------------------------------------------------------------------

end module statsmod
