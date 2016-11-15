module geohashmod

implicit none

integer, parameter :: dp = selected_real_kind(9)      !8 byte real
integer, parameter :: sp = selected_real_kind(4)      !4 byte real
integer, parameter :: i4 = selected_int_kind(9)       !10^9 fits into 4 bytes

integer(i4), parameter :: mv = huge(i4)  !add offset to fit into values into signed number

contains

!------

integer(i4) function geohash(lon,lat)

!encode a lon and lat pair into a signed 4-byte integer
!original use: to seed a random number generator in a predictable but geographically non-uniform way
!this should be precise and provide a unique number down to 30 arc-second resolution (1km)
!Jed O. Kaplan, 2011

implicit none

real(dp), intent(in) :: lon
real(dp), intent(in) :: lat

real(dp),    parameter :: scale  = 120.d0      !scale factor for geohash, the larger the number the more unique values
real(dp),    parameter :: offset =   0.5d0     !offset to calculate pixel number assuming gridcell center coordinates are given
integer(i4), parameter :: rowlen = nint(scale * 360.d0)

integer(i4) :: i,j

!---

i = nint(offset + scale * (lon + 180.d0))
j = nint(offset + scale * (lat + 90.d0))

geohash = i + rowlen * (j-1) - mv

end function geohash

!------

end module geohashmod
