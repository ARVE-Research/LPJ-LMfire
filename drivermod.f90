module drivermod

implicit none

contains

!------------------------------------------------------------------------------------------------------------

subroutine driver()

use mpistatevarsmod, only : in_master,sv_master
use mpimod,          only : initmpi,master
use initjobmod,      only : initjob
use getyrdatamod,    only : getdata
use iovariablesmod,  only : cal_year,dospinup,dotransient,cellmask,cfile_spinup,cfile_transient,nspinyrsout
use netcdfoutputmod, only : netcdf_output
use netcdfoutputmod, only : netcdf_close
use initclimatemod,  only : initclimate,closeclimate

implicit none

!local variables

integer :: ncells
integer :: ntiles
integer :: spinupyears
integer :: transientyears
integer :: firstyear

logical :: lastyear

integer :: year
integer :: i
integer :: tpos
integer :: time0

real :: time_begin
real :: time_end
real :: dt
real :: tsecs
integer :: thour
integer :: tmins

integer :: firstyrout

character(40) :: status_msg

!-------------------
!open input and output files, initialize gs

call initjob(ncells,ntiles,spinupyears,transientyears)

!-------------------
!open initial climate data files and allocate the climate input vector; allocate ibuf%

call initclimate(cfile_spinup,ncells)

!-------------------
!establish the size of the MPI buffer, etc.

call initmpi(ncells,ntiles)

lastyear  = .false.
firstyear = cal_year

time0 = 0

write(0,*)'calculating valid pixels'
call getdata(ncells,1,cal_year,firstyear,time0,in_master)        !returns gs filled with model input for first year as initial condition

write(0,*)'valid pixels at model start: ',count(cellmask)

tpos = 1

call cpu_time(time_begin)

!------------------
!spinup run

if (dospinup) then

  write(0,'(a,i6,a,f8.2)')'running spinup at calendar year: ',firstyear,' BP. CO2: ',in_master(1)%co2
  
  firstyrout = spinupyears - nspinyrsout
  
  do year = 1,spinupyears

    write(0,'(a18,2i8,f8.2)')' working on year: ',year,cal_year,in_master(1)%co2 !,year,lyear,co2(1)  !a,3i8,f8.2
    !call overprint(status_msg)

    in_master%spinup = .true.
    in_master%year = year

    do i = 1,12
      in_master%climate%temp0(i) = in_master%climate%temp(i)  !copy last year's temperature
    end do

    if (.not. dotransient .and. year == spinupyears) lastyear = .true.

    call getdata(ncells,year,cal_year,firstyear,time0,in_master)                   !returns gs filled with model input for this year

    call master(lastyear,ncells,in_master,sv_master)      !sends out and returns filled with model output
    
    if(year > firstyrout) then

       call netcdf_output(ncells,tpos,year-spinupyears,sv_master,in_master)
       
    end if   

  end do

  write(0,*)

  call closeclimate()

end if

!------------------
!transient run

if (dotransient) then

  write(0,*)'starting transient run'
  
  call initclimate(cfile_transient,ncells)
  time0 = 0

  do year = 1,transientyears

    write(0,'(a18,2i7,f8.2)')' working on year: ',year,cal_year,in_master(1)%co2 !,year,lyear,co2(1)  !a,3i8,f8.2
    !call overprint(status_msg)

    in_master%spinup = .false.
    in_master%year = year

    do i = 1,12
      in_master%climate%temp0(i) = in_master%climate%temp(i)  !copy last year's temperature
    end do

    if (year == transientyears) lastyear = .true.

    call getdata(ncells,year,cal_year,firstyear,time0,in_master)  !returns gs filled with model input for this year
    
    call master(lastyear,ncells,in_master,sv_master)        !sends out and returns filled with model output

    call netcdf_output(ncells,tpos,year,sv_master,in_master)
    
    cal_year = cal_year - 1

  end do

  write(0,*)

  call closeclimate()

end if

!------------------

write(0,*)'run finshed successfully'

call cpu_time(time_end)

dt = time_end - time_begin
thour = dt / 3600
tmins = mod(dt,3600.) / 60
tsecs = mod(dt,60.)

write(0,'(a,i5,a,i3,a,f5.1,a)') ' runtime:',thour,'h',tmins,'m',tsecs,'s'

call netcdf_close()

end subroutine driver

!------------------------------------------------------------------------------------------------------------

end module drivermod
