module netcdfoutputmod

use parametersmod, only : sp,i2

implicit none

public  :: netcdf_output
public  :: netcdf_close
private :: putvar2d
private :: putvar3d

real(sp),    parameter :: rmissing = -9999.
integer(i2), parameter :: smissing = -32768

contains

!----------------------------------------------------------------------------------------------------------------------------

subroutine netcdf_output(ncells,tpos,year,sv,in_master)

use typesizes
use netcdf
use errormod, only : ncstat,netcdf_err

use iovariablesmod,  only : cntx,cnty,ofid,cellmask,calcforagers
use parametersmod,   only : npft
use mpistatevarsmod, only : statevars,inputdata

implicit none

integer, intent(in)    :: ncells
integer, intent(inout) :: tpos
integer, intent(in)    :: year
type(statevars), dimension(:), intent(in) :: sv
type(inputdata), dimension(:), intent(in) :: in_master

integer, dimension(1) :: tval
!integer, dimension(1) :: lyear

integer :: varid

integer :: i,j,x,y
integer :: ntiles

real(sp), allocatable, dimension(:) :: rvar1d

real(sp), allocatable, dimension(:,:)     :: rvar2d
real(sp), allocatable, dimension(:,:,:,:) :: rvar4d

real(sp), dimension(3) :: NBP

!----------------------------
!write the time variable

tval = year

ncstat = nf90_inq_varid(ofid,'time',varid)
ncstat = nf90_put_var(ofid,varid,tval,start=[tpos],count=[1])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

!-----------------------------------------------------
!write other variables
!-----------------------------------------------------

allocate(rvar1d(ncells))

!number of tiles used in that run

ntiles = count(in_master(1)%human%landuse >=0.)

!---
!land fraction

do i = 1,ncells
  rvar1d(i) = in_master(i)%landf
end do

call putvar2d(ofid,tpos,'landf',rvar1d)

!---

if (calcforagers) then

  !forager population density

  do i = 1,ncells
    rvar1d(i) = sum(sv(i)%tile%forager_pd * sv(i)%tile%coverfrac) 
  end do

  call putvar2d(ofid,tpos,'foragerPD',rvar1d)

end if

!---
!livebiomass

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%livebiomass * sv(i)%tile%coverfrac) 
end do

call putvar2d(ofid,tpos,'livebiomass',rvar1d)

!---
!albiomass Chaste 03.03.2016

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%albiomass * sv(i)%tile%coverfrac) 
end do

call putvar2d(ofid,tpos,'albiomass',rvar1d)

!---
!litterC_fast

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%litterC_fast * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'litterC_fast',rvar1d)

!---
!litterC_slow

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%litterC_slow * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'litterC_slow',rvar1d)

!---
!litterC_bg

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%litterC_bg * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'litterC_bg',rvar1d)

!---
!surface SOM

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%cpool_surf(1) * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'soilC_surf',rvar1d)

!---
!soilC_fast

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%cpool_fast(1) * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'soilC_fast',rvar1d)

!---
!soilC_slow

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%cpool_slow(1) * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'soilC_slow',rvar1d)

!--- ! Ajout par E.C. avril 2016
!GPP

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%grid_gpp(1) * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'GPP',rvar1d)

!---
!NPP

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%grid_npp(1) * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'NPP',rvar1d)


!---
!burned fraction

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%afire_frac * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'burnedf',rvar1d)

!---
!fire carbon flux

do i = 1,ncells
  !rvar1d(i) = sv(i)%tile(1)%acflux_fire(1)
  rvar1d(i) = sum(sv(i)%tile%acflux_fire(1) * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'acflux_fire',rvar1d)

!---
!NBP

do i = 1,ncells
  NBP = sv(i)%tile%arh(1) + sv(i)%tile%acflux_fire(1) + sv(i)%tile%acflux_conv(1)   &
                          - sv(i)%tile%grid_npp(1)    - sv(i)%tile%acflux_estab(1)
  !rvar1d(i) = sum(NBP * sv(i)%tile%coverfrac) + sum(sv(i)%carbon%crop_harvest) * sv(i)%tile(2)%coverfrac + sum(sv(i)%carbon%prod_flux)
  !rvar1d(i) = NBP(2) + sv(i)%carbon%crop_harvest ! + sum(sv(i)%carbon%prod_flux)
  rvar1d(i) = sv(i)%tile(1)%grasscover
end do

call putvar2d(ofid,tpos,'NBP',rvar1d)

!---
!fire CO2

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%aMx(1) * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'fireCO2',rvar1d)

!---
!fire CO

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%aMx(2) * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'fireCO',rvar1d)

!---
!fire CH4

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%aMx(3) * sv(i)%tile%coverfrac)
end do

call putvar2d(ofid,tpos,'fireCH4',rvar1d)

!---
!fire VOC

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%aMx(4) * sv(i)%tile%coverfrac)
  !rvar1d(i) = sv(i)%tile(1)%aMx(4)				!FLAG: temporary use aMX(4) to output annual number of fires per km2
end do

call putvar2d(ofid,tpos,'fireVOC',rvar1d)

!---
!fire TPM

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%aMx(5) * sv(i)%tile%coverfrac)
  !rvar1d(i) = sv(i)%tile(1)%aMx(5)				!FLAG: temporary use aMX(5) to output annual average FDI on burndays
  !rvar1d(i) = sv(i)%tile(1)%forager_pd				!FLAG: temporary use forager PD
end do

call putvar2d(ofid,tpos,'fireTPM',rvar1d)

!---
!fire NOx

do i = 1,ncells
  rvar1d(i) = sum(sv(i)%tile%aMx(6) * sv(i)%tile%coverfrac)
  !rvar1d(i) = sum(sv(i)%tile%aMx(6) * sv(i)%tile%coverfrac)				!FLAG: temporary use aMX(6) to output the total annual area burned in the gridcell
end do

call putvar2d(ofid,tpos,'fireNOx',rvar1d)

!---

deallocate(rvar1d)

!------------------------
!array values (3d output)

!-----
!fpc_grid

y = size(sv(1)%tile(1)%fpc_grid)

allocate(rvar2d(ncells,y))

do i = 1,ncells
  do j = 1,npft
    !rvar2d(i,j) = sv(i)%tile(1)%fpc_grid(j)
    rvar2d(i,j) = sum(sv(i)%tile%fpc_grid(j) * sv(i)%tile%coverfrac)  !sum across tiles
  end do
end do

call putvar3d(ofid,tpos,'cover',rvar2d)

deallocate(rvar2d)

!-----
!tree height

y = size(sv(1)%tile(1)%height)

allocate(rvar2d(ncells,y))

do i = 1,ncells
  do j = 1,npft
    rvar2d(i,j) = sv(i)%tile(1)%height(j)
    !rvar2d(i,j) = sum(sv(i)%tile%fpc_grid(j) * sv(i)%tile%coverfrac)  !sum across tiles
  end do
end do

call putvar3d(ofid,tpos,'height',rvar2d)

deallocate(rvar2d)

!-----
!nind

y = size(sv(1)%tile(1)%nind)

!write(0,*) 'size nind', y

allocate(rvar2d(ncells,y))

do i = 1,ncells
  do j = 1,npft
    rvar2d(i,j) = sv(i)%tile(1)%nind(j)
    !rvar2d(i,j) = sum(sv(i)%tile%fpc_grid(j) * sv(i)%tile%coverfrac)  !sum across tiles
  end do
end do

call putvar3d(ofid,tpos,'nind',rvar2d)

deallocate(rvar2d)

!!-----
! Je remplace le crown area par la biomasse aboveground par PFTs en output 
!!crownarea

!y = size(sv(1)%tile(1)%crownarea)

!allocate(rvar2d(ncells,y))

!do i = 1,ncells
!  do j = 1,npft
!    rvar2d(i,j) = sv(i)%tile(1)%height(j) !crownarea(j)
!    !rvar2d(i,j) = sum(sv(i)%tile%fpc_grid(j) * sv(i)%tile%coverfrac)  !sum across tiles
!  end do
!end do

!call putvar3d(ofid,tpos,'crownarea',rvar2d)

!deallocate(rvar2d)

!-----
!Aboveground Biomass by PFTs

y = size(sv(1)%tile(1)%pftalbiomass)

allocate(rvar2d(ncells,y))

do i = 1,ncells
  do j = 1,npft
    rvar2d(i,j) = sv(i)%tile(1)%above_carbon(j) ! 
  end do
end do

call putvar3d(ofid,tpos,'pftalbiomass',rvar2d)

deallocate(rvar2d)

!-----
!tile fraction

y = size(sv(1)%tile%coverfrac)

allocate(rvar2d(ncells,y))

do i = 1,ncells
  rvar2d(i,:) = sv(i)%tile%coverfrac
end do

call putvar3d(ofid,tpos,'coverfrac',rvar2d)

deallocate(rvar2d)

!-----
!tile carbon

y = size(sv(1)%tile%tilecarbon)

allocate(rvar2d(ncells,y))

do i = 1,ncells
  rvar2d(i,:) = sv(i)%tile%tilecarbon
end do

call putvar3d(ofid,tpos,'tilecarbon',rvar2d)

deallocate(rvar2d)

!-----
!monthly burned area fraction of grid cell

y = size(sv(1)%tile(1)%mburnedf)	!ATTENTION: tile integration has already been done at end of lpjmod and put in tile 1, hence we only need to output tile1 for the integrated

allocate(rvar2d(ncells,y))

do i = 1,ncells
  do j = 1,12
    rvar2d(i,j) = sv(i)%tile(1)%mburnedf(j)
  end do
end do

call putvar3d(ofid,tpos,'mburnedf',rvar2d)

deallocate(rvar2d)

!-----

!monthly mean LAI, per PFT

allocate(rvar4d(cntx,cnty,npft,12))

rvar4d = rmissing

do i = 1,ncells
  
  x = in_master(i)%xpos
  y = in_master(i)%ypos

  if (.not.cellmask(x,y)) cycle

  rvar4d(x,y,:,:) = sv(i)%tile(1)%mLAI  !unpack into a 4D array

end do

ncstat = nf90_inq_varid(ofid,'mLAI',varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,rvar4d,start=[1,1,1,1,tpos],count=[cntx,cnty,npft,12,1])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

deallocate(rvar4d)

!-----
!monthly biomass burned, per PFT

allocate(rvar4d(cntx,cnty,npft,12))

rvar4d = 0.

do i = 1,ncells
  
  x = in_master(i)%xpos
  y = in_master(i)%ypos

  if (.not.cellmask(x,y)) cycle

  do j = 1, ntiles

    rvar4d(x,y,:,:) = rvar4d(x,y,:,:) + sv(i)%tile(j)%mBBpft * sv(i)%tile(j)%coverfrac
   
  end do  

!  rvar4d(x,y,:,:) = sv(i)%tile(1)%mBBpft  !unpack into a 4D array

end do

ncstat = nf90_inq_varid(ofid,'mBBpft',varid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_put_var(ofid,varid,rvar4d,start=[1,1,1,1,tpos],count=[cntx,cnty,npft,12,1])
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

deallocate(rvar4d)

!--------------------------------------------
!flush the write buffer every 30 years of run

if (mod(tpos,30) == 0) ncstat = nf90_sync(ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

tpos = tpos + 1

end subroutine netcdf_output

!----------------------------------------------------------------------------------------------------------------------------

subroutine putvar2d(ofid,tpos,varname,var1d)

use netcdf
use errormod, only : ncstat,netcdf_err
use iovariablesmod, only : cntx,cnty,cellmask

implicit none

integer,                intent(in) :: ofid
character(*),           intent(in) :: varname
!class(*), dimension(:), intent(in) :: var1d
real, dimension(:), intent(in) :: var1d
integer,                intent(in) :: tpos

real(sp), allocatable, dimension(:,:) :: rvar2d

integer :: varid

!-----------------

!select type(var1d)
!type is (real)
    
  allocate(rvar2d(cntx,cnty))
  
  rvar2d = reshape(var1d,[cntx,cnty])

  where (.not. cellmask) rvar2d = rmissing
  
  ncstat = nf90_inq_varid(ofid,varname,varid)
  if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_put_var(ofid,varid,rvar2d,start=[1,1,tpos],count=[cntx,cnty,1])
  if (ncstat/=nf90_noerr) call netcdf_err(ncstat)
    
  deallocate(rvar2d)

!class default

!  write(0,*)'error an output data type was requested that is not supported'
!  stop

!end select
  
end subroutine putvar2d

!--------------------------------------------

subroutine putvar3d(ofid,tpos,varname,var2d)

use netcdf
use errormod, only : ncstat,netcdf_err
use iovariablesmod, only : cntx,cnty,cellmask

implicit none

integer,                  intent(in) :: ofid
character(*),             intent(in) :: varname
!class(*), dimension(:,:), intent(in) :: var2d
real(sp), dimension(:,:), intent(in) :: var2d
integer,                  intent(in) :: tpos

real(sp), allocatable, dimension(:,:,:) :: rvar3d

integer :: varid
integer :: i,l
integer :: x,y

!-----------------

l = size(var2d,dim=2)

!select type(var2d)
!type is (real)
    
  allocate(rvar3d(cntx,cnty,l))
  
  y = 1
  do i = 1,size(var2d,dim=1)

    x = mod(i,cntx)
    y = 1+ i / cntx

    if (x==0) then
      x = cntx
      y = i/cntx
    end if
    
    rvar3d(x,y,:) = var2d(i,:)
  end do

  do i = 1,l
    where (.not. cellmask) rvar3d(:,:,i) = rmissing
  end do

  ncstat = nf90_inq_varid(ofid,varname,varid)
  if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

  ncstat = nf90_put_var(ofid,varid,rvar3d,start=[1,1,1,tpos],count=[cntx,cnty,l,1])
  if (ncstat/=nf90_noerr) call netcdf_err(ncstat)
    
  deallocate(rvar3d)

!class default

!  write(0,*)'error an output data type was requested that is not supported'
!  stop

!end select

end subroutine putvar3d

!----------------------------------------------------------------------------------------------------------------------------

subroutine netcdf_close()

use typesizes
use netcdf
use errormod, only : ncstat,netcdf_err

use iovariablesmod, only : ofid

implicit none

ncstat = nf90_close(ofid)
if (ncstat/=nf90_noerr) call netcdf_err(ncstat)

end subroutine netcdf_close

!----------------------------------------------------------------------------------------------------------------------------

end module netcdfoutputmod
