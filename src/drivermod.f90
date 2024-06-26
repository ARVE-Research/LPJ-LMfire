module drivermod

use parametersmod, only : stdout,stderr

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
integer :: nlayers
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

integer :: yrBP

integer :: firstyrout

character(100) :: status_msg

!-------------------
!open input and output files, initialize gs

call initjob(ncells,ntiles,nlayers,spinupyears,transientyears)

!-------------------
!open initial climate data files and allocate the climate input vector; allocate ibuf%

call initclimate(cfile_spinup,ncells)

!-------------------
!establish the size of the MPI buffer, etc.

call initmpi(ncells,ntiles,nlayers)

lastyear  = .false.
firstyear = cal_year

time0 = 0

write(stdout,*)'calculating valid pixels'
call getdata(ncells,1,cal_year,firstyear,time0,in_master)        !returns gs filled with model input for first year as initial condition

write(stdout,'(a,i20)')'valid pixels at model start: ',count(cellmask)

tpos = 1

call cpu_time(time_begin)

!------------------
!spinup run

if (dospinup) then

  write(stdout,'(a,i6,a,f8.2)')'starting spinup at calendar year: ',firstyear,' BP. CO2: ',in_master(1)%co2
  
  firstyrout = spinupyears - nspinyrsout
  
  do year = 1,spinupyears
  
    yrBP = cal_year ! + spinupyears - (year - 1) 

!    write(stdout,'(a18,3i8,f8.2)')' working on year: ',year,cal_year,yrBP,in_master(1)%co2  !,lyear,co2(1)  !a,3i8,f8.2
    write(stdout,'(a18,2i8,f8.2)')' working on year: ',year,cal_year,in_master(1)%co2 !,year,lyear,co2(1)  !a,3i8,f8.2
    write(status_msg,'(a,i6,a,i6)')' working on year',year,' out of',spinupyears
    ! call overprint(status_msg)

    in_master%spinup = .true.
    in_master%year = year

    do i = 1,12
      in_master%climate%temp0(i) = in_master%climate%temp(i)  !copy last year's temperature
    end do

    if (.not. dotransient .and. year == spinupyears) lastyear = .true.
    
    call getdata(ncells,year,cal_year,firstyear,time0,in_master)   !returns gs filled with model input for this year
    
     ! write(0,*)'driver',in_master(1:3)%cellarea
    
    call master(lastyear,ncells,in_master,sv_master)               !sends out and returns filled with model output
    
    if(year > firstyrout) then

       call netcdf_output(ncells,tpos,yrBP,sv_master,in_master)
       
       tpos = tpos + 1
 
    end if   

    cal_year = cal_year - 1

  end do

  write(stdout,*)

  call closeclimate()

end if

!------------------
!transient run

if (dotransient) then

  write(stdout,*)'starting transient run'
  
  call initclimate(cfile_transient,ncells)
  time0 = 0

  do year = 1,transientyears
  
    yrBP = cal_year
  
    in_master%spinup = .false.
    in_master%year = year

    do i = 1,12
      in_master%climate%temp0(i) = in_master%climate%temp(i)  !copy last year's temperature
    end do

    if (year == transientyears) lastyear = .true.

    call getdata(ncells,year,cal_year,firstyear,time0,in_master)  !returns gs filled with model input for this year

    write(stdout,'(a18,2i8,f8.2)')' working on year: ',year,cal_year,in_master(1)%co2 !,year,lyear,co2(1)  !a,3i8,f8.2
    write(status_msg,'(a,i6,a,i6)')' working on year',year,' out of',transientyears
    ! call overprint(status_msg)

    call master(lastyear,ncells,in_master,sv_master)        !sends out and returns filled with model output

    call netcdf_output(ncells,tpos,yrBP,sv_master,in_master)
    
    tpos = tpos + 1
    
    cal_year = cal_year - 1

  end do

  write(stdout,*)

  call closeclimate()

end if

!------------------

write(stderr,*)
write(stdout,*)'run finshed successfully'

call cpu_time(time_end)

dt = time_end - time_begin
thour = dt / 3600
tmins = mod(dt,3600.) / 60
tsecs = mod(dt,60.)

write(stdout,'(a,i5,a,i3,a,f5.1,a)') ' runtime:',thour,'h',tmins,'m',tsecs,'s'

call netcdf_close()

end subroutine driver

!------------------------------------------------------------------------------------------------------------

end module drivermod
