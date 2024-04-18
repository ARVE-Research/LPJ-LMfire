module landscape_geometrymod

use parametersmod, only : sp, pi

implicit none

public :: fragmentation

contains

! --------------------------------------------------
! NOTE this should also include total land fraction - 04.2024 have attempted something but need to check conceptually

subroutine fragmentation(landf,coverfrac,cellarea,avg_cont_area,totnat_area,avg_patch_number,median_distance,nbl)

implicit none

! Given a fractional land cover (a) present on a landscape, estimate the fraction of (a) that is contiguous with all other (a).
! Empirical relationship based on pixel counting of randomized distributions of (a) at concentrations.
! Discussed see Pfeiffer et al. (2013), see particularly Fig. S1.

! Landscape fragmentation, based on assumption that subgrid-kernels will be randomely distributed; 
! calculations based on a testgrid of 10000 sub-kernels. Calculating average contiguous area size, 
! depending on the tilefraction of tiles 1 + 3 the smallest size a natural patch can get is the kernel-size,
! no fractionation below kernel size allowed any more the kernel size defines the graininess of the landscape.


! arguments

real(sp),               intent(in)  :: landf            ! fraction of the gridcell that is dry land
real(sp), dimension(3), intent(in)  :: coverfrac        ! cover fractions of the three landuse tiles (1 = natural, 2 = cropland, 3 = recovering)
real(sp),               intent(in)  :: cellarea         ! total area of the gridcell, in m2
real(sp),               intent(out) :: avg_cont_area    ! average size of a natural patch at a given landuse fraction of the gridcell, in m2
real(sp),               intent(out) :: totnat_area      ! total natural area of a gridcell at a given landuse fraction, in m2   
real(sp),               intent(out) :: avg_patch_number ! average number of natural patches per gridcell at a given landuse fraction
real(sp),               intent(out) :: median_distance  ! patch radius minus the radius of circle whose area is exactly half the area of the natural patch; auxiliary for median distance to the edge of the natural patch, in (m)
real(sp),               intent(out) :: nbl              ! normalized boundary length; for boundary between natural and used part; normalized to the max. possible boundary length when having a chessboard-type distribution of kernels

! parameters

real(sp), parameter :: kernel_size   = 100000.    ! the smallest allowed size of a contiguous patch of homogenous land cover, in m2 (equals a size of 10 ha)
real(sp), parameter :: one_min_root  = 1. - sqrt(0.5)

! local variables

real(sp) :: nolanduse_frac
real(sp) :: patch_radius    ! radius of an average natural patch, assuming that the shape of the patch is a circle

real(sp) :: cont_area_frac

! ----

nolanduse_frac = landf * (1. - coverfrac(2))

cont_area_frac = (1.0025 + exp(16.607 - 41.503 * nolanduse_frac))**(-2.1694)

avg_cont_area = cont_area_frac * nolanduse_frac * cellarea

avg_cont_area = max(kernel_size,avg_cont_area)  ! average contiguous area size of non-used land (m2)

totnat_area = nolanduse_frac * cellarea

avg_patch_number = totnat_area / avg_cont_area

patch_radius = sqrt(avg_cont_area / pi)

median_distance = patch_radius * one_min_root

nbl = 1.9978 * nolanduse_frac - 1.9979 * nolanduse_frac**2

end subroutine fragmentation

! --------------------------------------------------

end module landscape_geometrymod
