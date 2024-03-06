module calendarmod

use iso_fortran_env, only : int16,int32,real32,real64
use ieee_arithmetic

implicit none

integer, parameter :: i2 = int16
integer, parameter :: i4 = int32
integer, parameter :: sp = real32
integer, parameter :: dp = real64

public  :: timestruct
public  :: ymdt2jd
public  :: jd2ymdt
private :: checkyear

type timestruct
  integer(i4) :: y    = -1
  integer(i4) :: m    = -1
  integer(i4) :: d    = -1
  integer(i4) :: hr   = -1
  integer(i4) :: min  = -1
  real(sp)    :: sec  = -1._sp
  real(dp)    :: jd   =  0._dp ! Julian date
end type timestruct

contains

! -----------------------------------------------------------------

integer(i4) function checkyear(y) result(y2)

implicit none

integer(i4), intent(in) :: y

! ---

if ( y < 0 ) then
  y2 = y + 1
else if ( y == 0 ) then
  write(0,*)'error, calendar does not have a year 0'
  stop
else
  y2 = y
end if

end function checkyear

! -----------------------------------------------------------------

subroutine ymdt2jd(dt)

! Retrieve a Julian date based on a (proleptic) Gregorian year-month-day-[time] combination
! Based on routines in calpak.f90 by John Burkardt

implicit none

! argument

type(timestruct), intent(inout) :: dt

! local variables

integer(i4) :: d_prime
integer(i4) :: g
integer(i4) :: ierror
integer(i4) :: j1
integer(i4) :: j2
integer(i4) :: m_prime
integer(i4) :: y2
integer(i4) :: y_prime
real(dp)    :: timefrac

! -----------------
! Check if the time is provided, if not fall back to default of 00:00:00.0

timefrac = 0.

if (dt%hr >= 0) then
  if (dt%min >= 0 .and. dt%sec >= 0.) then
    timefrac = (real(dt%hr) + real(dt%min) / 60. + dt%sec / 3600.) / 24.
  else
    write(0,*)'must include full time string or nothing at all, falling back to midnight default'
  end if
end if

! ---

! Convert the calendar date to a computational date

y2 = checkyear(dt%y)

y_prime = y2 + 4716 - (14 - dt%m) / 12
m_prime = mod (dt%m + 9,12)
d_prime = dt%d - 1

! Convert the computational date to a Julian date

j1 = (1461 * y_prime) / 4

j2 = (153 * m_prime + 2) / 5

g = (3 * ((y_prime + 184) / 100) / 4) - 38

dt%jd = real(j1 + j2 + d_prime - 1401 - g,kind=dp) - 0.5_dp + timefrac

end subroutine ymdt2jd

! -----------------------------------------------------------------

subroutine jd2ymdt(dt)

! Retrieve a (proleptic) Gregorian year-month-day-[time] combination based on a Julian date
! Based on routines in calpak.f90 by John Burkardt

implicit none

! argument

type(timestruct), intent(inout) :: dt

! local variables

integer(i4) :: g
integer(i4) :: j
integer(i4) :: j_prime
integer(i4) :: y_prime
integer(i4) :: m_prime
integer(i4) :: d_prime
integer(i4) :: t_prime
real(sp)    :: f
real(dp)    :: fpart
real(dp)    :: hr
real(dp)    :: mr

! Determine the computational date (Y'/M'/D').

j = int(dt%jd + 0.5_dp)
f = (dt%jd + 0.5_dp) - real(j,kind=dp)

g = (3 * ((4 * j + 274277) / 146097) / 4) - 38

j_prime = j + 1401 + g

y_prime =    (4 * j_prime + 3) / 1461
t_prime = mod(4 * j_prime + 3, 1461) / 4
m_prime =    (5 * t_prime + 2) / 153
d_prime = mod(5 * t_prime + 2, 153) / 5

!  Convert the computational date to a calendar date.

dt%d = d_prime + 1
dt%m = mod(m_prime + 2,12) + 1
dt%y = y_prime - 4716 + (14 - dt%m) / 12

! If there is a remainder calculate the time - 
! NB this algorithm resolves to noon on the selected computational date

fpart = dt%jd - floor(dt%jd)

if (fpart > 0.) then

  hr = 24. * fpart
  mr = 60. * (hr - floor(hr))
  dt%sec = 60. * (mr - floor(mr))

  dt%hr  = int(hr)
  dt%min = int(mr)

else

  dt%hr  = 0
  dt%min = 0
  dt%sec = 0.

end if

! Any year before 1 AD must be moved one year further back, since
! this calendar does not include a year 0.

if (dt%y <= 0) dt%y = dt%y - 1

end subroutine jd2ymdt

! -----------------------------------------------------------------

end module calendarmod
