module drivermod_serial

implicit none

contains

!------------------------------------------------------------------------------------------------------------

subroutine driver()

use mpistatevarsmod, only : in_master,sv_master
use initjobmod,      only : initjob
use getyrdatamod,    only : getdata
use iovariablesmod,  only : cal_year,dospinup,dotransient,cellmask
use netcdfoutputmod, only : netcdf_output
use netcdfoutputmod, only : netcdf_close
use lpjmod,          only : lpjcore
use radiationmod,    only : initairmass
use parametersmod,   only : sp,lm_sapl,sm_sapl,rm_sapl,hm_sapl,pftpar,pft, &
                            allom1,allom2,allom3,allom4,sla,latosa,wooddens,reinickerp
!uses pftparameters.f

implicit none

real(sp), parameter, dimension(3) :: co2 = [ 273., -8., 0. ]

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
integer :: x,y
integer :: time0

real :: time_begin
real :: time_end

character(40) :: status_msg

!---------------------

call initjob(ncells,ntiles,spinupyears,transientyears)     !open input and output files, initialize gs

!-------------------------------
!initialize the airmass parameters

call initairmass()

!--------------------------
!initialize the pft parameters

call pftparameters(pftpar,sla,                                                 &
                   pft%tree,pft%evergreen,pft%summergreen,                     &
                   pft%raingreen,pft%needle,pft%boreal,                        &
                   lm_sapl,sm_sapl,hm_sapl,rm_sapl,                            &
                   latosa,allom1,allom2,allom3,allom4,wooddens,reinickerp,co2)

!year loop (equilibrium model means only one year)

lastyear  = .false.
firstyear = cal_year

time0 = 0

write(0,*)'reading initial input data'
call getdata(ncells,1,cal_year,firstyear,time0,in_master)        !returns gs filled with model input for first year as initial condition

write(0,*)'valid pixels at model start: ',count(cellmask)

!write(0,'(12f8.2)')in_master(1)%climate%temp

tpos = 1

call cpu_time(time_begin)

!------------------
!spinup run

if (dospinup) then

  write(0,'(a,i8,a,i6,a,f8.2)')'running ',spinupyears,' year spinup with calendar year: ',firstyear,' BP. CO2: ',in_master(1)%co2
  
  do year = 1,spinupyears

    write(0,'(a18,2i8)')' working on year: ',year,cal_year !,year,lyear,co2(1)  !a,3i8,f8.2
    !call overprint(status_msg)

    in_master%spinup = .true.
    in_master%year = year

    do i = 1,12
      in_master%climate%temp0(i) = in_master%climate%temp(i)  !copy last year's temperature
    end do

    if (.not. dotransient .and. year == spinupyears) lastyear = .true.

    call getdata(ncells,1,cal_year,firstyear,time0,in_master)        !returns gs filled with model input for first year as initial condition

    do i = 1,ncells  !loop over grid cells in the job
      
      x = in_master(i)%xpos
      y = in_master(i)%ypos
      !  write(0,*)i,in_master(i)%idx,'check lpjcore spinup'

      if (.not.cellmask(x,y)) cycle

      call lpjcore(in_master(i),sv_master(i))

    end do

    call netcdf_output(ncells,tpos,year-spinupyears,sv_master)

  end do

  write(0,*)

end if

!------------------
!transient run

if (dotransient) then

  write(0,*)'transient'

  do year = 1,transientyears

    write(0,'(a18,2i7,f8.2)')' working on year: ',year,cal_year,in_master(1)%co2 !,year,lyear,co2(1)  !a,3i8,f8.2
!    call overprint(status_msg)

    in_master%spinup = .false.
    in_master%year = year

    do i = 1,12
      in_master%climate%temp0(i) = in_master%climate%temp(i)  !copy last year's temperature
    end do

    if (year == transientyears) lastyear = .true.

    call getdata(ncells,1,cal_year,firstyear,time0,in_master)        !returns gs filled with model input for first year as initial condition

    do i = 1,ncells  !loop over grid cells in the job

      x = in_master(i)%xpos
      y = in_master(i)%ypos     
        
      if (.not.cellmask(x,y)) cycle	 

      call lpjcore(in_master(i),sv_master(i))

    end do

    call netcdf_output(ncells,tpos,year,sv_master)
    
    cal_year = cal_year - 1

  end do

  write(0,*)

end if

!------------------

call cpu_time(time_end)

write(0,*) 'Model was running for ', time_end - time_begin, 'seconds'

call netcdf_close()

write(0,*)'done'

!NB for the transient run, decrement cal_year in the year loop here

end subroutine driver

!------------------------------------------------------------------------------------------------------------

end module drivermod_serial
