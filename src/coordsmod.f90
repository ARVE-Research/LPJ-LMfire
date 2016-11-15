module coordsmod

implicit none

public :: parsecoords
public :: calcpixels

contains

!-----------------------------------------------------------------------------------------------

subroutine parsecoords(coordstring,val)

!subroutine to parse a coordinate string

implicit none

character(45),      intent(in)  :: coordstring
real, dimension(4), intent(out) :: val

character(10), dimension(4) :: cval = '0'

integer :: i
integer :: lasti = 1
integer :: part  = 1

do i=1,len_trim(coordstring)
  if (coordstring(i:i) == '/') then
    cval(part) = coordstring(lasti:i-1)
    lasti=i+1
    part = part + 1
  end if
end do

cval(part) = coordstring(lasti:i-1)

read(cval,*)val

if (part < 4) then
  val(3)=val(2)
  val(4)=val(3)
  val(2)=val(1)
end if

end subroutine parsecoords

!-----------------------------------------------------------------------------------------------

subroutine calcpixels(bounds,srtx,srty,cntx,cnty,gridres,lboundlon,uboundlat)

implicit none

real, intent(in), dimension(4) :: bounds
real, intent(in), dimension(2) :: gridres
real, intent(in) :: lboundlon
real, intent(in) :: uboundlat

integer, intent(out) :: srtx
integer, intent(out) :: srty
integer, intent(out) :: cntx
integer, intent(out) :: cnty

!local variables

!these values are hardwired at the moment, but shouldn't be
!real :: grid_resx  =   30.   !grid resolution in minutes
!real :: grid_resy  =   30.   !grid resolution in minutes

real, dimension(2) :: pixfact

real :: minlon
real :: maxlon
real :: minlat
real :: maxlat

!-------------------------
!calculate the number of pixels to be calculated

pixfact = 60. / gridres

srtx = 1 + nint(pixfact(1) * (bounds(1) - lboundlon))
srty = 1 + nint(pixfact(2) * (bounds(3) - bounds(4)))
!srty = 1 + nint(pixfact(2) * (uboundlat - bounds(4)))

write(*,*)gridres
write(*,*)bounds
write(*,*)srtx,srty
read(*,*)

minlon = bounds(1)
if (bounds(2) == bounds(1)) then
  maxlon = bounds(1)  + (1./pixfact(1))
else
  maxlon = bounds(2)
end if

minlat = bounds(3)
if (bounds(4) == bounds(3)) then
  maxlat = bounds(3)  + (1./pixfact(2))
else
  maxlat = bounds(4)
end if

cntx = nint(pixfact(1) * (maxlon - minlon))
cnty = nint(pixfact(2) * (maxlat - minlat))

if (cntx <= 0 .or. cnty <= 0) then
  write(0,*)cntx,cnty
  stop 'invalid coordinates'
end if

end subroutine calcpixels

!-----------------------------------------------------------------------------------------------

end module coordsmod
