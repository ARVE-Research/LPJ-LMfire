module newsplinemod

! This is newspline, a fast, means-preserving spline for interval data
! Leo O Lai and Jed O. Kaplan
! 2021

use parametersmod, only : i4,sp
use utilitiesmod,  only : matsol

implicit none

public  :: newspline          ! Newspline subroutine with all optional bounded adjustment schemes (version 2.0)

private :: alim_adjust        ! Bounded interpolation adjustment scheme for absolute limit (i.e., bound tolerance)
private :: plim_adjust        ! Bounded interpolation adjustmnet scheme for percentage limit (relative to each interval individually)
private :: ulim_adjust        ! Bounded interpolation adjustmnet scheme for maximum limit
private :: llim_adjust        ! Bounded interpolation adjustmnet scheme for minimum limit

contains

! ------------------------------------------------------------------------------------------------------------------

subroutine newspline(monthdata,nk,bcond,daydata,alim,llim,ulim,plim)

implicit none

! arguments

real(sp),    dimension(:), intent(in)  :: monthdata          ! Array of (monthly) interval data
integer(i4), dimension(:), intent(in)  :: nk                 ! Array of number of small time-steps for each interval (can be variable)
real(sp),    dimension(2), intent(in)  :: bcond              ! Boundary condition array, bcond(1) = bcond of first interval and bcond(2) = bcond of last interval
real(sp),    dimension(:), intent(out) :: daydata            ! Array of interpolated values (dimension must be equal to sum of nk)
real(sp),    optional    , intent(in)  :: alim               ! OPTIONAL :: Absolute limit for bounded interpolation
real(sp),    optional    , intent(in)  :: llim               ! OPTIONAL :: Minimum limit for bounded interpolation
real(sp),    optional    , intent(in)  :: ulim               ! OPTIONAL :: Maximum limit for bounded interpolation
real(sp),    optional    , intent(in)  :: plim               ! OPTIONAL :: Percetnage limit for bounded interpolation

! Local variables for wall control points

real(sp), dimension(:), allocatable :: wcp

! Local variables for linear system of mid control adjustments

real(sp),    dimension(:,:), allocatable :: mat
real(sp),    dimension(:),   allocatable :: solution

! Final vector of all control points

real(sp), dimension(:), allocatable :: all_cont

! Local variables for generating daily values

real(sp), dimension(:), allocatable :: m_cont

! Hermite cubic and quartic spline basis functions

real(sp) :: H_00, H_10, H_01, H_11
real(sp) :: G_00, G_10, G_01, G_11
real(sp) :: u
real(sp) :: del

integer :: len
integer :: len_cont
integer :: i
integer :: j
integer :: k
integer :: n

! -------------------------------------------------------------------------------
! Start of the spline routine

len = size(monthdata)

! -------------------------------------------------------------------------------
! Approximate first wall control points using linear interpolation by connecting yi with yi+1
! and calculating the interception with interval wall
! 
! New possible simplified method upon revision by Leo O Lai (Sep 2021)
! Since we assume the interval width to be arbitary 1 unit, the wall-CPs are simply the average of yi and yi+1
! First and last interval wall-CP (i.e. wcp(1) and wcp(N+1)) are determined by user specified boundary conditions

allocate(wcp(len+1))

wcp(1)     = (monthdata(1) + bcond(1)) / 2.
wcp(len+1) = (monthdata(len) + bcond(2)) / 2.

do i = 2, len

  wcp(i) = (monthdata(i-1) + monthdata(i)) / 2

end do

! -------------------------------------------------------------------------------
! Generate matrix for solution to mid-control points
allocate(mat(len,len))
allocate(solution(len))

! ====================
! Version 3 amendements - fixed G_00 integral u^4 term sign from '-' to '+'
! ====================

u = 1.

G_00 = u - (u**3) + (u**4) / 2.
G_10 = (u**2) * (3.*(u**2) - 8.*u + 6.) / 12.
G_01 = (u**3) * (1. - u / 2.)
G_11 = (u**3) * (3.*u - 4.) / 12.

! ---

mat = 0.

! Construct matrix coefficients and solution array

! ====================
! Version 3 amendements - elimited '+ 1.' terms in tridiagonal coefficients and solution array after new equation derivation
! ====================

! Consider two "buffer mid-CPs" outside of first and last interval

mat(1,1) = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11
mat(1,2) = 0.5 * G_11

mat(len,len-1) = -0.5 * G_10
mat(len,len)   = 0.5 * G_10 + G_01 + G_00 - 0.5 * G_11


n = 1
do i = 2, (len-1)

  mat(i,n)   = -0.5 * G_10

  mat(i,n+1) = G_00 + G_01 + 0.5 * G_10 - 0.5 * G_11

  mat(i,n+2) = 0.5 * G_11

  n = n + 1

end do

! ---

solution(1)   = 2 * monthdata(1) - &
                (G_00 - 0.5 * G_10 - 0.5 * G_11) * wcp(1) - &
                (G_01 + 0.5 * G_10 + 0.5 * G_11) * wcp(2)

solution(1)   = solution(1) + 0.5 * G_10 * bcond(1)

solution(len) = 2 * monthdata(len) - &
                (G_00 - 0.5 * G_10 - 0.5 * G_11) * wcp(len) - &
                (G_01 + 0.5 * G_10 + 0.5 * G_11) * wcp(len+1)

solution(len) = solution(len) - 0.5 * G_11 * bcond(2)


do i = 2, (len-1)

  solution(i) = 2 * monthdata(i) - &
                (G_00 - 0.5 * G_10 - 0.5 * G_11) * wcp(i) - &
                (G_01 + 0.5 * G_10 + 0.5 * G_11) * wcp(i+1)

end do


! -------------------------------------------------------------------------------
! Solve linear system to get final mid control points
call matsol(mat,solution)


! -------------------------------------------------------------------------------
! Compile wall control with newly adjusted mid control points (all_cont)
allocate(all_cont(size(wcp)+size(solution)))

n = 1
do i = 1, len

  all_cont(n)   = wcp(i)
  all_cont(n+1) = solution(i)

  n = n + 2

end do

all_cont(2*len+1) = wcp(len+1)


! -------------------------------------------------------------------------------
! Construct the spline for daily values based on all_cont

len_cont = size(all_cont)

! ---
! Calculate the slope for each control point

allocate(m_cont(len_cont))

m_cont(1)        = all_cont(2) - bcond(1)
m_cont(len_cont) = bcond(2) - all_cont(len_cont-1)

do i = 2, len_cont-1

  m_cont(i) = all_cont(i+1) - all_cont(i-1)

end do


! -------------------------------------------------------------------------------
! Find discrete daily values on the continuous function (i.e., series of u values in Hermite functions)

! ====================
! Version 3 amendements - include delta = 1/2 here instead of m_cont above to be consistent with equation in manuscript
! ====================

del = 0.5

n = 1
k = 1

do i = 1, len ! Outer loop start, for all monthly intervals N

  if(mod(nk(i),2) == 0) then ! Seperate into even or odd months, starting with EVEN

    u = 1. / nk(i)

    do j = 1, (nk(i) / 2)

      H_00 = 1. + (u**2) * (2*u - 3)
      H_01 = (u**2) * (3 - 2*u)
      H_10 = u * (u - 1) * (u - 1)
      H_11 = (u**2) * (u - 1)

      daydata(n) = all_cont(k) * H_00 + del * m_cont(k) * H_10 + all_cont(k+1) * H_01 + del * m_cont(k+1) * H_11

      u = u + (2. / nk(i))

      n = n + 1

    end do

    ! ---

    u = 1. / nk(i)

    do j = 1, (nk(i) / 2)

      H_00 = 1. + (u**2) * (2*u - 3)
      H_01 = (u**2) * (3 - 2*u)
      H_10 = u * (u - 1) * (u - 1)
      H_11 = (u**2) * (u - 1)

      daydata(n) = all_cont(k+1) * H_00 + del * m_cont(k+1) * H_10 + all_cont(k+2) * H_01 + del * m_cont(k+2) * H_11

      u = u + (2. / nk(i))

      n = n + 1

    end do

  else ! if odd months (or odd number of smaller time step)

    u = 1. / nk(i)

    do j = 1, ((nk(i)+1) / 2)

      H_00 = 1. + (u**2) * (2*u - 3)
      H_01 = (u**2) * (3 - 2*u)
      H_10 = u * (u - 1) * (u - 1)
      H_11 = (u**2) * (u - 1)

      daydata(n) = all_cont(k) * H_00 + del * m_cont(k) * H_10 + all_cont(k+1) * H_01 + del * m_cont(k+1) * H_11

      u = u + (2. / nk(i))

      n = n + 1

    end do

    ! ---

    u = 2. / nk(i)

    do j = 1, ((nk(i)-1) / 2)

      H_00 = 1. + (u**2) * (2*u - 3)
      H_01 = (u**2) * (3 - 2*u)
      H_10 = u * (u - 1) * (u - 1)
      H_11 = (u**2) * (u - 1)

      daydata(n) = all_cont(k+1) * H_00 + del * m_cont(k+1) * H_10 + all_cont(k+2) * H_01 + del * m_cont(k+2) * H_11

      u = u + (2. / nk(i))

      n = n + 1

    end do

  end if

  k = k + 2

end do ! End of outer loop


! -------------------------------------------------------------------------------
! Call bounded interpolation adjustment scheme if optional arguments are present

if (present(llim)) call llim_adjust(llim,monthdata,nk,bcond,all_cont,daydata)
if (present(ulim)) call ulim_adjust(ulim,monthdata,nk,bcond,all_cont,daydata)
if (present(alim)) call alim_adjust(alim,monthdata,nk,bcond,all_cont,daydata)
if (present(plim)) call plim_adjust(plim,monthdata,nk,bcond,all_cont,daydata)

end subroutine newspline

! ------------------------------------------------------------------------------------------------------------------

subroutine alim_adjust(alim,monthdata,nk,bcond,all_cont,daydata)

real(sp),                  intent(in)    :: alim        ! Absolute limit (e.g., no interpolated value can exceed +/- 0.05 of original interval mean input)
real(sp),    dimension(:), intent(in)    :: monthdata   ! Array of monthly (interval) input data
integer(i4), dimension(:), intent(in)    :: nk          ! Array of number of small time-steps for each interval (can be variable)
real(sp),    dimension(:), intent(in)    :: bcond       ! Boundary condition array
real(sp),    dimension(:), intent(in)    :: all_cont    ! Array of all control points (wall-CPs and mid-CPs)
real(sp),    dimension(:), intent(inout) :: daydata     ! Array of daily intepolated values

real(sp), allocatable, dimension(:) :: d_orig           ! Slope direction of the current interval (1 = positive, -1 = negative, 0 = local maxima/minima)
logical,  allocatable, dimension(:) :: osc_check        ! TRUE if interval require adjustment (dim = input data length)
real(sp), allocatable, dimension(:) :: c2               ! Array to store the amount of adjustment required
real(sp), allocatable, dimension(:) :: c_mon            ! Array to store current month (or interval) of values for bounded adjustment (dim = day in month)

real(sp) :: perc

real(sp) :: int_n
real(sp) :: int_nm1
real(sp) :: int_np1
real(sp) :: sip12

integer :: len
integer :: i
integer :: j
integer :: k
integer :: i0
integer :: i1

! -----

len = size(monthdata)

allocate(d_orig(len))
allocate(osc_check(len))
allocate(c2(len))

d_orig    = -9999.
osc_check = .false.
c2        = -9999.

! ------
! Assign -1 for negative slope, 1 for postive and 0 for turning point

do i = 1, len

  ! Assign intervals values to local variables
  int_n = monthdata(i)        ! Current interval

  if (i == 1) then              ! If first interval, apply bcond(1)
    int_nm1 = bcond(1)
    int_np1 = monthdata(i+1)
  else if (i == len) then       ! If last interval, apply bcond(2)
    int_nm1 = monthdata(i-1)
    int_np1 = bcond(2)
  else
    int_nm1 = monthdata(i-1)
    int_np1 = monthdata(i+1)
  end if

  ! Determine the slope direction of the current interval relative to previous and next interval

  if (int_np1 - int_n < 0) then

    d_orig(i) = -1

  else if (int_np1 - int_n > 0) then

    d_orig(i) = 1

  end if

  ! ---

  if (int_n > int_nm1 .and. int_n > int_np1) then

    d_orig(i) = 0

  else if (int_n < int_nm1 .and. int_n < int_np1) then

    d_orig(i) = 0

  end if

end do

! ------
! Assign TRUE if oscillation of turning point exceeds the predetermined threshold

do i = 1, len

  ! Assign intervals values to variables
  int_n = monthdata(i)        ! Current interval
  sip12 = all_cont(2*i)       ! Position of interval mid-control point (i.e. si+1/2)

  if (d_orig(i) == 0) then

    if (sip12 > int_n + alim .OR. sip12 < int_n - alim) then

      osc_check(i) = .TRUE.

    end if

  else if (d_orig(i) == 0 .and. int_n < 0) then

    if (sip12 < (1.0+perc) * int_n .OR. sip12 > (1.0-perc) * int_n) then

      osc_check(i) = .TRUE.

    end if

  end if

end do


! ------
! Calculate the amount of adjustment required and insert into the c2 variable

c2 = 0.

do i = 1, len

  ! Assign intervals values to local variables
  int_n   = monthdata(i)        ! Current interval

  if (i == 1) then              ! If first interval, apply bcond(1)
    int_nm1 = bcond(1)
    int_np1 = monthdata(i+1)
  else if (i == len) then       ! If last interval, apply bcond(2)
    int_nm1 = monthdata(i-1)
    int_np1 = bcond(2)
  else
    int_nm1 = monthdata(i-1)
    int_np1 = monthdata(i+1)
  end if

  ! ---

  if (osc_check(i)  .and. int_n > 0) then

    ! ---

    if (int_n > int_nm1 .and. int_n > int_np1) then

      c2(i) = alim

    else if (int_n < int_nm1 .and. int_n < int_np1) then

      c2(i) = -alim

    end if

    ! ---

  else if (osc_check(i)  .and. int_n < 0) then

    ! ---

    if (int_n > int_nm1 .and. int_n > int_np1) then

      c2(i) = alim

    else if (int_n < int_nm1 .and. int_n < int_np1) then

      c2(i) = -alim

    end if

    ! ---

  end if

end do


! ------
! Construct the spline for daily values based on all_cont

do i = 1, len

  if (osc_check(i)) then

    ! Assign intervals values to local variables
    int_n = monthdata(i)        ! Current interval

    if (i == 1) then              ! If first interval, apply bcond(1)
      int_nm1 = bcond(1)
      int_np1 = monthdata(i+1)
    else if (i == len) then       ! If last interval, apply bcond(2)
      int_nm1 = monthdata(i-1)
      int_np1 = bcond(2)
    else
      int_nm1 = monthdata(i-1)
      int_np1 = monthdata(i+1)
    end if

    ! ---

    if (i == 1) then
      i0 = 1
      i1 = i0 + nk(1) - 1
    else
      i0 = sum(nk(1:(i-1))) + 1
      i1 = i0 + nk(i) - 1
    end if

    ! ---
    ! Create seperate array of all the daily values within the interval being adjusted
    allocate(c_mon(nk(i)+2))

    if (i == 1) then
      c_mon(1)         = all_cont(1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = daydata(i1+1)
    else if (i == len) then
      c_mon(1)         = daydata(i0-1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = all_cont(2*len+1)
    else
      c_mon(1)         = daydata(i0-1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = daydata(i1+1)
    end if

    ! ---

    do j = 1, 1000

      if ((int_n > int_nm1 .and. int_n > int_np1) .and. &
          sum(c_mon(2:nk(i)+1)) / nk(i) - int_n < 0.01) exit

      if ((int_n < int_nm1 .and. int_n < int_np1) .and. &
          sum(c_mon(2:nk(i)+1)) / nk(i) - int_n > 0.01) exit


      do k = 2, nk(i)+1

        c_mon(k) = (c_mon(k-1) + c_mon(k) + c_mon(k+1)) / 3.

      end do

    end do

    ! ---

    daydata(i0:i1) = c_mon(2:nk(i)+1)

    deallocate(c_mon)

  end if

end do

! n_adjust = count(osc_check)

end subroutine alim_adjust

! ------------------------------------------------------------------------------------------------------------------

subroutine plim_adjust(plim,monthdata,nk,bcond,all_cont,daydata)

real(sp),                  intent(in)    :: plim        ! Percentage limit (e.g. no interpolated value can exceed 5% of original interval mean input)
real(sp),    dimension(:), intent(in)    :: monthdata   ! Array of monthly (interval) input data
integer(i4), dimension(:), intent(in)    :: nk          ! Array of number of small time-steps for each interval (can be variable)
real(sp),    dimension(:), intent(in)    :: bcond       ! Boundary condition array
real(sp),    dimension(:), intent(in)    :: all_cont    ! Array of all control points (wall-CPs and mid-CPs)
real(sp),    dimension(:), intent(inout) :: daydata     ! Array of daily intepolated values

real(sp), allocatable, dimension(:) :: d_orig           ! Slope direction of the current interval (1 = positive, -1 = negative, 0 = local maxima/minima)
logical,  allocatable, dimension(:) :: osc_check        ! TRUE if interval require adjustment (dim = input data length)
real(sp), allocatable, dimension(:) :: c2               ! Array to store the amount of adjustment required
real(sp), allocatable, dimension(:) :: c_mon            ! Array to store current month (or interval) of values for bounded adjustment (dim = day in month)

real(sp) :: perc

real(sp) :: int_n
real(sp) :: int_nm1
real(sp) :: int_np1
real(sp) :: sip12

integer :: len
integer :: i
integer :: j
integer :: k
integer :: i0
integer :: i1

! -----

len = size(monthdata)

allocate(d_orig(len))
allocate(osc_check(len))
allocate(c2(len))

d_orig = -9999.
osc_check = .false.
c2 = -9999.


! ------
! Assign -1 for negative slope, 1 for postive and 0 for turning point
do i = 1, len

  ! Assign intervals values to local variables
  int_n   = monthdata(i)        ! Current interval

  if (i == 1) then              ! If first interval, apply bcond(1)
    int_nm1 = bcond(1)
    int_np1 = monthdata(i+1)
  else if (i == len) then       ! If last interval, apply bcond(2)
    int_nm1 = monthdata(i-1)
    int_np1 = bcond(2)
  else
    int_nm1 = monthdata(i-1)
    int_np1 = monthdata(i+1)
  end if

  ! Determine the slope direction of the current interval relative to previous and next interval

  if (int_np1 - int_n < 0) then

    d_orig(i) = -1

  else if (int_np1 - int_n > 0) then

    d_orig(i) = 1

  end if

  ! ---

  if (int_n > int_nm1 .and. int_n > int_np1) then

    d_orig(i) = 0

  else if (int_n < int_nm1 .and. int_n < int_np1) then

    d_orig(i) = 0

  end if

end do


! ------
! Assign TRUE if oscillation of turning point exceeds the predetermined threshold
perc = plim / 100.

do i = 1, len

  ! Assign intervals values to variables
  int_n = monthdata(i)        ! Current interval
  sip12 = all_cont(2*i)       ! Position of interval mid-control point (i.e. si+1/2)

  ! ---

  if (d_orig(i) == 0 .and. monthdata(i) > 0) then

    if (sip12 > (1.0+perc) * int_n .OR. sip12 < (1.0-perc) * int_n) then

      osc_check(i) = .TRUE.

    end if

  else if (d_orig(i) == 0 .and. int_n < 0) then

    if (sip12 < (1.0+perc) * int_n .OR. sip12 > (1.0-perc) * int_n) then

      osc_check(i) = .TRUE.

    end if

  end if

end do


! ------
! Calculate the amount of adjustment required and insert into the c2 variable
c2 = 0.

do i = 1, len

  ! Assign intervals values to local variables
  int_n   = monthdata(i)        ! Current interval

  if (i == 1) then              ! If first interval, apply bcond(1)
    int_nm1 = bcond(1)
    int_np1 = monthdata(i+1)
  else if (i == len) then       ! If last interval, apply bcond(2)
    int_nm1 = monthdata(i-1)
    int_np1 = bcond(2)
  else
    int_nm1 = monthdata(i-1)
    int_np1 = monthdata(i+1)
  end if

  ! ---

  if (osc_check(i)  .and. int_n > 0) then

    ! ---

    if (int_n > int_nm1 .and. int_n > int_np1) then

      c2(i) = perc * int_n

    else if (int_n < int_nm1 .and. int_n < int_np1) then

      c2(i) = -perc * int_n

    end if

    ! ---

  else if (osc_check(i)  .and. int_n < 0) then

    ! ---

    if (int_n > int_nm1 .and. int_n > int_np1) then

      c2(i) = -perc * int_n

    else if (int_n < int_nm1 .and. int_n < int_np1) then

      c2(i) = perc * int_n

    end if

    ! ---

  end if

end do


! ------
! Construct the spline for daily values based on all_cont
do i = 1, len

  if (osc_check(i)) then

    ! Assign intervals values to local variables
    int_n = monthdata(i)        ! Current interval

    if (i == 1) then              ! If first interval, apply bcond(1)
      int_nm1 = bcond(1)
      int_np1 = monthdata(i+1)
    else if (i == len) then       ! If last interval, apply bcond(2)
      int_nm1 = monthdata(i-1)
      int_np1 = bcond(2)
    else
      int_nm1 = monthdata(i-1)
      int_np1 = monthdata(i+1)
    end if

    ! ---

    if (i == 1) then
      i0 = 1
      i1 = i0 + nk(1) - 1
    else
      i0 = sum(nk(1:(i-1))) + 1
      i1 = i0 + nk(i) - 1
    end if

    ! ---
    ! Create seperate array of all the daily values within the interval being adjusted
    allocate(c_mon(nk(i)+2))

    if (i == 1) then
      c_mon(1)         = all_cont(1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = daydata(i1+1)
    else if (i == len) then
      c_mon(1)         = daydata(i0-1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = all_cont(2*len+1)
    else
      c_mon(1)         = daydata(i0-1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = daydata(i1+1)
    end if

    ! ---

    do j = 1, 1000

      if ((int_n > int_nm1 .and. int_n > int_np1) .and. &
          sum(c_mon(2:nk(i)+1)) / nk(i) - int_n < 0.01) exit

      if ((int_n < int_nm1 .and. int_n < int_np1) .and. &
          sum(c_mon(2:nk(i)+1)) / nk(i) - int_n > 0.01) exit


      do k = 2, nk(i)+1

        c_mon(k) = (c_mon(k-1) + c_mon(k) + c_mon(k+1)) / 3.

      end do

    end do

    ! ---

    daydata(i0:i1) = c_mon(2:nk(i)+1)

    deallocate(c_mon)

  end if

end do


end subroutine plim_adjust

! ------------------------------------------------------------------------------------------------------------------

subroutine ulim_adjust(ulim,monthdata,nk,bcond,all_cont,daydata)

real(sp),                  intent(in)    :: ulim        ! Maximum limit (e.g. no interpolated value can exceed 1.0 in the ENTIRE interpolated series)
real(sp),    dimension(:), intent(in)    :: monthdata   ! Array of monthly (interval) input data
integer(i4), dimension(:), intent(in)    :: nk          ! Array of number of small time-steps for each interval (can be variable)
real(sp),    dimension(:), intent(in)    :: bcond       ! Boundary condition array
real(sp),    dimension(:), intent(in)    :: all_cont    ! Array of all control points (wall-CPs and mid-CPs)
real(sp),    dimension(:), intent(inout) :: daydata     ! Array of daily intepolated values


real(sp), allocatable, dimension(:) :: d_orig           ! Slope direction of the current interval (1 = positive, -1 = negative, 0 = local maxima/minima)
logical,  allocatable, dimension(:) :: osc_check        ! TRUE if interval require adjustment (dim = input data length)
real(sp), allocatable, dimension(:) :: c2               ! Array to store the amount of adjustment required
real(sp), allocatable, dimension(:) :: c_mon            ! Array to store current month (or interval) of values for bounded adjustment (dim = day in month)

real(sp) :: int_n
real(sp) :: int_nm1
real(sp) :: int_np1
real(sp) :: sip12

integer :: len
integer :: i
integer :: j
integer :: k
integer :: i0
integer :: i1

! -----

len = size(monthdata)

allocate(d_orig(len))
allocate(osc_check(len))
allocate(c2(len))

d_orig = -9999.
osc_check = .false.
c2 = -9999.


! ------
! Assign -1 for negative slope, 1 for postive and 0 for turning point
do i = 1, len

  ! Assign intervals values to local variables
  int_n = monthdata(i)        ! Current interval

  if (i == 1) then              ! If first interval, apply bcond(1)
    int_nm1 = bcond(1)
    int_np1 = monthdata(i+1)
  else if (i == len) then       ! If last interval, apply bcond(2)
    int_nm1 = monthdata(i-1)
    int_np1 = bcond(2)
  else
    int_nm1 = monthdata(i-1)
    int_np1 = monthdata(i+1)
  end if

  ! Determine the slope direction of the current interval relative to previous and next interval

  if (int_np1 - int_n < 0) then

    d_orig(i) = -1

  else if (int_np1 - int_n > 0) then

    d_orig(i) = 1

  end if

  ! ---

  if (int_n > int_nm1 .and. int_n > int_np1) then

    d_orig(i) = 0

  else if (int_n < int_nm1 .and. int_n < int_np1) then

    d_orig(i) = 0

  end if

end do


! ------
! Assign TRUE if oscillation of turning point exceeds the predetermined threshold
do i = 1, len

  sip12 = all_cont(2*i)

  if (sip12 > ulim .and. d_orig(i) == 0) then

    osc_check(i) = .TRUE.

  end if

end do


! ------
! Calculate the amount of adjustment required and insert into the c2 variable
c2 = 0.

do i = 1, len

  if (osc_check(i)) then

    c2(i) = ulim - monthdata(i)

  end if

end do


! ------
! Construct the spline for daily values based on all_cont
do i = 1, len

  if (osc_check(i)) then

    ! Assign intervals values to local variables
    int_n = monthdata(i)        ! Current interval

    if (i == 1) then              ! If first interval, apply bcond(1)
      int_nm1 = bcond(1)
      int_np1 = monthdata(i+1)
    else if (i == len) then       ! If last interval, apply bcond(2)
      int_nm1 = monthdata(i-1)
      int_np1 = bcond(2)
    else
      int_nm1 = monthdata(i-1)
      int_np1 = monthdata(i+1)
    end if

    ! ---

    if (i == 1) then
      i0 = 1
      i1 = i0 + nk(1) - 1
    else
      i0 = sum(nk(1:(i-1))) + 1
      i1 = i0 + nk(i) - 1
    end if

    ! ---
    ! Create seperate array of all the daily values within the interval being adjusted
    allocate(c_mon(nk(i)+2))

    if (i == 1) then
      c_mon(1)         = all_cont(1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = daydata(i1+1)
    else if (i == len) then
      c_mon(1)         = daydata(i0-1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = all_cont(2*len+1)
    else
      c_mon(1)         = daydata(i0-1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = daydata(i1+1)
    end if

    ! ---

    do j = 1, 1000

      if ((int_n > int_nm1 .and. int_n > int_np1) .and. &
          sum(c_mon(2:nk(i)+1)) / nk(i) - int_n < 0.01) exit

      if ((int_n < int_nm1 .and. int_n < int_np1) .and. &
          sum(c_mon(2:nk(i)+1)) / nk(i) - int_n > 0.01) exit


      do k = 2, nk(i)+1

        c_mon(k) = (c_mon(k-1) + c_mon(k) + c_mon(k+1)) / 3.

      end do

    end do

    ! ---

    daydata(i0:i1) = c_mon(2:nk(i)+1)

    deallocate(c_mon)

  end if

end do


end subroutine ulim_adjust

! ------------------------------------------------------------------------------------------------------------------

subroutine llim_adjust(llim,monthdata,nk,bcond,all_cont,daydata)

real(sp),                  intent(in)    :: llim        ! Minimum limit (e.g. no interpolated value can exceed 1.0 in the ENTIRE interpolated series)
real(sp),    dimension(:), intent(in)    :: monthdata   ! Array of monthly (interval) input data
integer(i4), dimension(:), intent(in)    :: nk          ! Array of number of small time-steps for each interval (can be variable)
real(sp),    dimension(:), intent(in)    :: bcond       ! Boundary condition array
real(sp),    dimension(:), intent(in)    :: all_cont    ! Array of all control points (wall-CPs and mid-CPs)
real(sp),    dimension(:), intent(inout) :: daydata     ! Array of daily intepolated values


real(sp), allocatable, dimension(:) :: d_orig           ! Slope direction of the current interval (1 = positive, -1 = negative, 0 = local maxima/minima)
logical,  allocatable, dimension(:) :: osc_check        ! TRUE if interval require adjustment (dim = input data length)
real(sp), allocatable, dimension(:) :: c2               ! Array to store the amount of adjustment required
real(sp), allocatable, dimension(:) :: c_mon            ! Array to store current month (or interval) of values for bounded adjustment (dim = day in month)

real(sp) :: int_n
real(sp) :: int_nm1
real(sp) :: int_np1
real(sp) :: sip12

integer :: len
integer :: i
integer :: j
integer :: k
integer :: i0
integer :: i1

! -----

len = size(monthdata)

allocate(d_orig(len))
allocate(osc_check(len))
allocate(c2(len))

d_orig = -9999.
osc_check = .false.
c2 = -9999.


! ------
! Assign -1 for negative slope, 1 for postive and 0 for turning point
do i = 1, len

  ! Assign intervals values to local variables
  int_n = monthdata(i)        ! Current interval

  if (i == 1) then              ! If first interval, apply bcond(1)
    int_nm1 = bcond(1)
    int_np1 = monthdata(i+1)
  else if (i == len) then       ! If last interval, apply bcond(2)
    int_nm1 = monthdata(i-1)
    int_np1 = bcond(2)
  else
    int_nm1 = monthdata(i-1)
    int_np1 = monthdata(i+1)
  end if

  ! Determine the slope direction of the current interval relative to previous and next interval

  if (int_np1 - int_n < 0) then

    d_orig(i) = -1

  else if (int_np1 - int_n > 0) then

    d_orig(i) = 1

  end if

  ! ---

  if (int_n > int_nm1 .and. int_n > int_np1) then

    d_orig(i) = 0

  else if (int_n < int_nm1 .and. int_n < int_np1) then

    d_orig(i) = 0

  end if

end do



! ------
! Assign TRUE if oscillation of turning point exceeds the predetermined threshold
do i = 1, len

  sip12 = all_cont(2*i)

  if (sip12 < llim .and. d_orig(i) == 0) then

    osc_check(i) = .TRUE.

  end if

end do


! ------
! Calculate the amount of adjustment required and insert into the c2 variable
c2 = 0.

do i = 1, len

  if (osc_check(i)) then

    c2(i) = llim - monthdata(i)

  end if

end do


! ------
! Construct the spline for daily values based on all_cont
do i = 1, len

  if (osc_check(i)) then

    ! Assign intervals values to local variables
    int_n = monthdata(i)        ! Current interval

    if (i == 1) then              ! If first interval, apply bcond(1)
      int_nm1 = bcond(1)
      int_np1 = monthdata(i+1)
    else if (i == len) then       ! If last interval, apply bcond(2)
      int_nm1 = monthdata(i-1)
      int_np1 = bcond(2)
    else
      int_nm1 = monthdata(i-1)
      int_np1 = monthdata(i+1)
    end if

    ! ---

    if (i == 1) then
      i0 = 1
      i1 = i0 + nk(1) - 1
    else
      i0 = sum(nk(1:(i-1))) + 1
      i1 = i0 + nk(i) - 1
    end if

    ! ---
    ! Create seperate array of all the daily values within the interval being adjusted
    allocate(c_mon(nk(i)+2))

    if (i == 1) then
      c_mon(1)         = all_cont(1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = daydata(i1+1)
    else if (i == len) then
      c_mon(1)         = daydata(i0-1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = all_cont(2*len+1)
    else
      c_mon(1)         = daydata(i0-1)
      c_mon(2:nk(i)+1) = int_n + c2(i)
      c_mon(nk(i)+2)   = daydata(i1+1)
    end if

    ! ---

    do j = 1, 1000

      if ((int_n > int_nm1 .and. int_n > int_np1) .and. &
          sum(c_mon(2:nk(i)+1)) / nk(i) - int_n < 0.01) exit

      if ((int_n < int_nm1 .and. int_n < int_np1) .and. &
          sum(c_mon(2:nk(i)+1)) / nk(i) - int_n > 0.01) exit


      do k = 2, nk(i)+1

        c_mon(k) = (c_mon(k-1) + c_mon(k) + c_mon(k+1)) / 3.

      end do

    end do

    ! ---

    daydata(i0:i1) = c_mon(2:nk(i)+1)

    deallocate(c_mon)

  end if

end do


end subroutine llim_adjust

! ------------------------------------------------------------------------------------------------------------------

end module newsplinemod
