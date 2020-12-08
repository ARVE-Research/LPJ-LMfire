module netcdfsetupmod

use parametersmod, only : maxoutvars
use parametersmod, only : stdout,stderr

implicit none

public  :: netcdf_create
private :: getvarinfo
private :: declvar

integer, parameter :: max_ndims  = 10

type variableinfo
  character(20) :: name
  character(80) :: longname
  character(20) :: units
  integer       :: ndims
  character(20), dimension(max_ndims) :: dimnames
end type variableinfo

type(variableinfo), target, dimension(maxoutvars) :: varinfo

character(50), dimension(maxoutvars) :: reqvarsout

contains

!------------------------------------------------------------------------------------------------------------------------------

subroutine netcdf_create(ncells)

use typesizes
use netcdf
use errormod, only : ncstat,netcdf_err

use parametersmod,  only : npft,lutype
use iovariablesmod, only : ofid,lonvect,latvect,srtx,cntx,endx,srty,cnty,endy,outputfile,outputvar,cellindex,cellmask,calchumans,cal_year

implicit none

integer, intent(in) :: ncells

integer :: dimid
integer :: varid
integer :: i
integer :: x,y

integer :: reqoutvars

character(8)  :: today
character(10) :: now

real, dimension(2) :: xrange
real, dimension(2) :: yrange

real :: xres
real :: yres

integer, allocatable, dimension(:) :: pftnum

character(40) :: yearstring

character(40), dimension(2), parameter :: varlabel = ['landf' ,'foragerPD']  !names of variables that will be automatically output

!----------

call getvarinfo()

!if (size(lonvect) > 1 .and. size(latvect) > 1) then
!  xres = (lonvect(endx) - lonvect(srtx)) / (size(lonvect) - 1)
!  yres = (latvect(endy) - latvect(srty)) / (size(latvect) - 1)
!else
  xres = 0.5
  yres = 0.5
!end if

write(stdout,*)'input grid resolution',xres,yres,' degrees'

xrange(1) = minval(lonvect(srtx:endx)) - 0.5 * xres
xrange(2) = maxval(lonvect(srtx:endx)) + 0.5 * xres

yrange(1) = minval(latvect(srty:endy)) - 0.5 * yres
yrange(2) = maxval(latvect(srty:endy)) + 0.5 * yres

allocate(cellmask(cntx,cnty))

cellmask = .false.

!----------------------------------
!dimensions and dimension variables

write(stdout,'(a,a)')'output file: ',trim(outputfile)

ncstat = nf90_create(outputfile,nf90_hdf5,ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'title','LPJ netCDF output file')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

call date_and_time(today,now)

ncstat = nf90_put_att(ofid,nf90_global,'timestamp',today//' '//now(1:4))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'Conventions','COARDS')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,nf90_global,'node_offset',1)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

write(stdout,*)'added global atts'

!----

ncstat = nf90_def_dim(ofid,'lon',cntx,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'lon',nf90_float,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','longitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_east')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',xrange)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'lat',cnty,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'lat',nf90_float,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','latitude')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','degrees_north')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',yrange)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'layer',2,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'layer',nf90_short,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','soil layer (0-30cm, 30cm-bottom)')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','layer')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'tile',size(lutype),dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'tile',nf90_short,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','land use tile')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','layer')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'pft',npft,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'pft',nf90_short,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','Plant Functional Type')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','PFT')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'month',12,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'month',nf90_short,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','month')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'units','month')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

ncstat = nf90_def_dim(ofid,'time',nf90_unlimited,dimid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_def_var(ofid,'time',nf90_int,dimid,varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'long_name','time')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

write(yearstring,'(a,i4,a)')'years since ',1950-cal_year,'-01-01'

ncstat = nf90_put_att(ofid,varid,'units',trim(yearstring))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'calendar','standard')
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!---------------------------------------------------------------------------
!regular variables

write(stdout,*)'done declaring dimension variables'

outputvar = eoshift(outputvar,-1,varlabel(1))
  
if (calchumans) then
  outputvar = eoshift(outputvar,-1,varlabel(2))
end if

reqoutvars = count(outputvar /= 'null')

if (reqoutvars < 1) then
  write(stdout,*)'WARNING, no variables were specified for output!'
    
else
  write(stdout,*)'the following variables will be written to netCDF output'
  
  do i = 1,reqoutvars

    write(stdout,*)outputvar(i)
    call declvar(ofid,outputvar(i))
  
  end do
end if

!---------------------------------------------------------------------------

ncstat = nf90_enddef(ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----

!write the dimension variables (except for time)

allocate(pftnum(npft))

forall (i=1:npft)
  pftnum(i) = i
end forall

ncstat = nf90_inq_varid(ofid,'lon',varid)
ncstat = nf90_put_var(ofid,varid,lonvect(srtx:endx))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'lat',varid)
ncstat = nf90_put_var(ofid,varid,latvect(srty:endy))
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'layer',varid)
ncstat = nf90_put_var(ofid,varid,[1,2])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'pft',varid)
ncstat = nf90_put_var(ofid,varid,pftnum)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'tile',varid)
ncstat = nf90_put_var(ofid,varid,[1,2,3])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(ofid,'month',varid)
ncstat = nf90_put_var(ofid,varid,[1,2,3,4,5,6,7,8,9,10,11,12])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!----------------

deallocate(pftnum)

end subroutine netcdf_create

!------------------------------------------------------------------------------------------------------------------------------

subroutine getvarinfo()

implicit none

namelist / varinfolist / varinfo

open(87,file='outvarsinfo.namelist',status='old')

read(87,varinfolist)

close(87)

end subroutine getvarinfo

!--------------------------------------------

subroutine declvar(fid,varname)

use typesizes
use netcdf
use errormod, only : ncstat,netcdf_err

implicit none

integer,      intent(in) :: fid
character(*), intent(in) :: varname

character(20), pointer :: name
character(80), pointer :: longname
character(20), pointer :: units

character(20), pointer, dimension(:) :: dimnames

real :: missing = -9999.

integer, pointer :: ndims

integer :: i,j

integer :: dimlen
integer :: varid

integer, allocatable, dimension(:) :: dimids
integer, allocatable, dimension(:) :: chunks

!------
!scan varinfo for the variable to be processed

do i = 1,maxoutvars
  if (trim(varinfo(i)%name) == trim(varname)) then
    j = i
    exit
  else
    j = 0
  end if
end do

if (j == 0) then
  write(stdout,*)'ERROR, a variable was requested for output that was not present'
  write(stdout,*)'in the attribute table: outvarsinfo.namelist'
  write(stdout,*)'Please check the variable name or add metadata to the table.'
  stop
end if  

name     => varinfo(i)%name
longname => varinfo(i)%longname
units    => varinfo(i)%units
ndims    => varinfo(i)%ndims
dimnames => varinfo(i)%dimnames

allocate(dimids(ndims))
allocate(chunks(ndims))

do i = 1,ndims
  
  ncstat = nf90_inq_dimid(fid,dimnames(i),dimids(i))
  if (ncstat/=nf90_noerr) call netcdf_err(ncstat)
    
  ncstat = nf90_inquire_dimension(fid,dimids(i),len=dimlen)
  if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

  if (i <= 2) then
    chunks(i) = dimlen
  else
    chunks(i) = 1
  end if
  
end do

ncstat = nf90_def_var(fid,varname,nf90_float,dimids,varid,chunksizes=chunks,deflate_level=1,shuffle=.false.)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(fid,varid,'long_name',longname)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(fid,varid,'units',units)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(fid,varid,'_FillValue',missing)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_att(fid,varid,'missing_value',missing)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

deallocate(dimids)
deallocate(chunks)

end subroutine declvar

!------------------------------------------------------------------------------------------------------------------------------

end module netcdfsetupmod
