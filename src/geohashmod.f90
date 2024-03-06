module geohashmod

use parametersmod, only : sp,dp,i4

implicit none

integer(i4), parameter :: mv = huge(i4)  !add offset to fit into values into signed number

contains

! ------

integer(i4) function geohash(lon,lat)

! encode a lon and lat pair into a signed 4-byte integer
! original use: to seed a random number generator in a predictable but geographically non-uniform way
! this should be precise and provide a unique number down to 30 arc-second resolution (1km)
! Jed O. Kaplan, 2011

implicit none

real(dp), intent(in) :: lon
real(dp), intent(in) :: lat

real(dp),    parameter :: scale  = 120._dp      !scale factor for geohash, the larger the number the more unique values
real(dp),    parameter :: offset =   0.5_dp     !offset to calculate pixel number assuming gridcell center coordinates are given
integer(i4), parameter :: rowlen = nint(scale * 360._dp)

integer(i4) :: i,j

! ---

i = nint(offset + scale * (lon + 180._dp))
j = nint(offset + scale * (lat + 90._dp))

geohash = i + rowlen * (j-1) - mv

end function geohash

! ------

end module geohashmod
