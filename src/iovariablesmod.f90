module iovariablesmod

use parametersmod, only : sp,dp,i2,i4
use calendarmod,   only : timestruct

implicit none

integer :: ofid

real(sp) :: ocean_uptake

! required input files
character(220) :: cfile_spinup  ! a climate spinup file can be used alone for an equilibrium run
character(220) :: soilfile      ! file containing sand, clay, and lithology
character(220) :: topofile      ! file containing elevation, slope, and land fraction
character(220) :: pftparsfile   ! file containing pftparameters

! optional input files
character(220) :: cfile_transient = ''
character(220) :: co2file  = ''  ! or specify fixed CO2
character(220) :: poppfile = ''  ! file with potential population density of foragers
character(220) :: popdfile = ''  ! file with population density of different people

integer :: startyr_foragers   ! first year of spinup at which to introduce foragers (to allow vegetation to develop first)

! flag for projected grid

logical :: projgrid

integer, parameter :: maxoutvars = 100   ! maximum number of variables for output

character(40), dimension(maxoutvars) :: outputvar = 'null'

character(80) :: timeunit_climate
type(timestruct) :: timeunit_basedate

integer  :: cfid
integer  :: cfid1
integer  :: cfid2

integer  :: soilfid
integer  :: topofid

logical  :: calcforagers
logical  :: lucc
integer  :: popfid

integer  :: condfid

integer, allocatable, dimension(:) :: topotime

integer :: elvid
integer :: slopeid
integer :: landfid

real(sp) :: temp_sf
real(sp) :: prec_sf
real(sp) :: sunp_sf

real(sp) :: temp_ao
real(sp) :: prec_ao
real(sp) :: sunp_ao

real(sp) :: lu_turn_yrs = 0. ! absolute turnover time for anthopogenic land use (years)

logical  :: dospinup
logical  :: dotransient
logical  :: nolanduse

integer  :: climateyears
integer  :: nspinyrsout = 0   ! number of years of the spinup to write out
integer  :: cal_year    = 0   ! calendar year BP (1950) for the spinup and first year of the transient run

integer  :: inputlatlen
integer  :: inputlonlen
integer  :: inputclimlen
integer  :: climatemonths

integer  :: dtimelen    ! the length of the time dimension in the detrended climatology
integer  :: ttimelen    ! the length of the time dimension in the transient climatology

integer  :: timebuflen  ! the number of months of climate data to store in the input buffer

real(sp) :: maxmem      ! maximum amount of memory to use for input buffer (Mb),
                        ! specify in joboptions file or will be set to 512 Mb

!  real(sp), allocatable, save, dimension(:,:,:) :: cbuf ! a buffer for multi-year climate data

! maxibufsize sets the maximum number of gridcells the input and output buffers can hold.
! A maximum sized box of 259200 pixels (global 0.5 deg resolution), has some 67420 valid 
! (non-water non-ice pixels). Global soils data and 12 months of climate data results in 
! a memory requirement of ca. 13 MB

integer, parameter :: maxibufsize = 300000

real(dp), dimension(4) :: bounds

real(dp), dimension(2) :: gridres

real(dp) :: minlon
real(dp) :: maxlon
real(dp) :: minlat
real(dp) :: maxlat

integer :: xlen
integer :: ylen

real(sp), allocatable, dimension(:,:) :: landfrac  ! optional values that can be used when fractional land data is available
real(sp), allocatable, dimension(:,:) :: icefrac

real(sp) :: fixedco2

logical, allocatable, dimension(:,:)  :: cellmask

logical :: dosoilco2

! ----------------

type jobinfovars
  integer(i4) :: bufsize
  integer(i4) :: ntiles
  integer(i4) :: nlayers
  character(220) :: pftparsfile
end type jobinfovars

type(jobinfovars) :: jobinfo

! ----------------

type inputbuffer

  real(dp) :: lon
  real(dp) :: lat
  real(sp), dimension(3) :: popd                 ! human population density (person km-2) - three categories
  real(sp) :: lu_turnover                        ! fraction of the land use part of the gridcell that is turned over every year
  real(sp), allocatable, dimension(:) :: temp    ! mean monthly temperature (degC)
  real(sp), allocatable, dimension(:) :: prec    ! total monthly precipitation (mm)
  real(sp), allocatable, dimension(:) :: cldp    ! mean monthly cloud cover (%)
  real(sp), allocatable, dimension(:) :: wetd    ! total monthly days with precipitation (days)
  real(sp), allocatable, dimension(:) :: trng    ! diurnal temperature range (degC)
  real(sp), allocatable, dimension(:) :: temp0   ! mean monthly temperature of the previous year (degC)
  real(sp), allocatable, dimension(:) :: wind    ! mean monthly wind speed (m s-1)
  real(sp), allocatable, dimension(:) :: lght    ! total lightning flashes (flashes km-2 day-1)
  real(sp), allocatable, dimension(:) :: landuse ! fraction of gridcell used for the 5 types of land use: forage, crop, pasture, rangeland, urban
  real(sp), allocatable, dimension(:) :: co2     ! CO2 concentration (ppm) (total CO2, 13C, and 14C)

end type inputbuffer

type soiltype

  logical    :: water
  integer(2) :: elv
  real(sp)   :: slopeangle   ! mean slope angle (radians)
  real(sp)   :: landf       ! fraction of the grid cell that is land
  real(sp), allocatable, dimension(:) :: zpos
  real(sp), allocatable, dimension(:) :: whc
  real(sp), allocatable, dimension(:) :: cond
  real(sp), allocatable, dimension(:) :: sand
  real(sp), allocatable, dimension(:) :: clay
  real(sp), allocatable, dimension(:) :: orgm

end type soiltype

real,              allocatable, target, dimension(:)   :: co2vect
type(inputbuffer), allocatable, target, dimension(:,:) :: ibuf
type(soiltype),    allocatable, target, dimension(:,:) :: soil

real(sp),    allocatable, dimension(:,:)     :: foragerPD ! hunter-gatherer potential population density
integer(i2), allocatable, dimension(:,:,:)   :: input_i2
real(sp),    allocatable, dimension(:,:,:,:) :: input_sp

! ----------------

type inputvarinfo
  integer  :: varid
  real(sp) :: add_offset
  real(sp) :: scale_factor
end type inputvarinfo

integer, parameter :: nclimv = 7  ! number of climate variables read by the program

type(inputvarinfo), dimension(nclimv) :: varinfo

! ---

real(sp), allocatable, dimension(:) :: dco2
real(sp), allocatable, dimension(:) :: tco2

integer, allocatable, dimension(:,:) :: indexarray
integer, allocatable, dimension(:,:) :: indexmat

integer :: index         ! the index of the current gridcell in the input buffer

integer, parameter :: outyears    =   10
integer, parameter :: outmons     =   12 * outyears

integer :: loopyear  = 1
integer :: transyear = 1

real(dp), allocatable, dimension(:) :: lonvect
real(dp), allocatable, dimension(:) :: latvect

real(sp), allocatable, dimension(:,:) :: geolon  ! geodetic longitude and latitude for projected grids
real(sp), allocatable, dimension(:,:) :: geolat


real(sp), allocatable, dimension(:) :: timevect
real(sp), allocatable, dimension(:) :: outtimevect

real(sp), allocatable, dimension(:,:,:,:) :: obuf_pft   ! x,y,variable,pft (or pft)
real(sp), allocatable, dimension(:,:,:,:) :: obuf_time  ! x,y,variable,time

! integer, allocatable, dimension(:) :: lucctime
real(dp), allocatable, dimension(:) :: popdtime

integer :: srtx
integer :: srty
integer :: cntx
integer :: cnty
integer :: endx
integer :: endy

character(120) :: outputfile

integer, allocatable, save, dimension(:,:) :: cellindex

end module iovariablesmod
