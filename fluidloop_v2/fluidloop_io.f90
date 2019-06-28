module fluidloop_io
! Program fluidloop
!
! solves the 1D loop problem with point sources/sinks of heat and salt
! Parameters are set in namelist
! Circular and folded loop configurations
! Nonlinear equation of state
! Temperature relaxation and salinity fixed fluxes, applied at same height
! Possibility to initialize model with predefined temperature and salinity fields
!
! version 0.1: F. Pollmann and F. Roquet, 20/11/2013 (fixed-flux forcing, nonlinear EOS, leap-frog, namelist, netCDF)
! version 0.2: R. Lindqvist and F. Roquet,  13/04/2015 (mixed boundary condition, non-zero initial conditions)
! version 0.3: F. Roquet, 31/07/2015 (readability improvements, code reorganisation, non-inertialess case)
! version 0.4: F. Roquet, 20/08/2015 (elliptic shape)
! version 0.5: F. Roquet, 07/01/2016 (ascii outputs, netCDF capability removed)
! version 0.6: F. Roquet, 22/01/2016 (forcing inputs)
! version 1.0: F. Roquet, 22/01/2016 (simplification of the code: Remove elliptic shape, simplify forcing parameters, and remove convergence test)
!
! loop toolbox, distributed on GitHub: http://github.com/fabien-roquet/loop
! F. Roquet 2016
! GNU General Public License


USE netcdf
USE fluidloop_init
implicit none

PUBLIC open_output_files
PUBLIC close_output_files
PUBLIC write_results

PRIVATE netcdf_setup
PRIVATE netcdf_write
PRIVATE netcdf_close


integer :: ncid, status
integer :: rec_varid, height_varid
integer :: theta_varid, salt_varid, sigma_varid, N2_varid, w_varid
integer :: theta_adv_varid, theta_diff_varid, salt_adv_varid, salt_diff_varid, varid
integer :: theta_for_varid, salt_for_varid, time_varid
logical :: wr_adv


contains

!----------------------------------------------------------------------------------------!
subroutine open_output_files

  character (len=50) :: name_out_1, name_out_2

  ! open output files
  if (.NOT. use_netcdf) then
     write(name_out_1,'(a,"_out_0D.txt")') trim(name_exp)
     open(unit=11,file=name_out_1,action='write',status='new')
     write(11,*) '% niter time velocity mass'
     write(name_out_2,'(a,"_out_1D.txt")') trim(name_exp)
     open(unit=12,file=name_out_2,action='write',status='new')
     write(12,*) '% niter j theta salt sigma N2'
  else
     call netcdf_setup
  end if
  
end subroutine open_output_files


!----------------------------------------------------------------------------------------!
subroutine close_output_files

  ! finalize program
  print*, 'last velocity', w_n, 'at time=', time
  
  if (.NOT. use_netcdf) then
     close(unit=11)
     close(unit=12)
  else
     call netcdf_close
  end if
  
end subroutine close_output_files


!----------------------------------------------------------------------------------------!
subroutine write_results
  integer            :: jk
  
	! Write data 1D
	if (.NOT. use_netcdf) then
     write (unit=11,fmt="(i10,3(f10.4))") niter,time,w_n,mass
     do jk = 1,nl
	     write (unit=12,fmt="(2(i10),4(f10.3))") niter,jk,theta_n(jk),salt_n(jk),sigma(jk),N2(jk)
	   end do
  else
	   CALL netcdf_write
  end if
  
end subroutine write_results


!----------------------------------------------------------------------------------------!
subroutine netcdf_setup

	integer :: rec_dimid,pos_dimid,pos_varid,dimids(2)

  character (len=50) :: file_name
	character (len = *), parameter :: REC_NAME = "record"
	character (len = *), parameter :: POS_NAME = "index"
	character (len = *), parameter :: UNITS = "units"
	character (len = *), parameter :: POS_UNITS = "loop index (starting from top)"
	character (len = *), parameter :: TIME_UNITS = "record index"

  wr_adv = .FALSE.
  
  
	! Create output file
  write(file_name,'(a,"_output.nc")') trim(name_exp)
	status = nf90_create(file_name,nf90_clobber,ncid)

	! Define dimensions
	status = nf90_def_dim(ncid, POS_NAME, nl, pos_dimid)
	status = nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid)

	! Define coordinate variable  
	status = nf90_def_var(ncid, POS_NAME, NF90_INT, pos_dimid, pos_varid)
	status = nf90_def_var(ncid, REC_NAME, NF90_INT, rec_dimid, rec_varid)

	! Assign unit attributes
	status = nf90_put_att(ncid, pos_varid, UNITS, POS_UNITS)
	status = nf90_put_att(ncid, rec_varid, UNITS, TIME_UNITS)

	dimids = (/ pos_dimid, rec_dimid/)
	
	! Define netCDF variables
	status = nf90_def_var(ncid, 'w', NF90_REAL8, rec_dimid, w_varid)
	status = nf90_def_var(ncid, 'z', NF90_REAL8, pos_dimid, height_varid)
	status = nf90_def_var(ncid, 'time', NF90_INT, rec_dimid, time_varid)
	status = nf90_def_var(ncid, 'theta', NF90_REAL8, dimids, theta_varid)
	status = nf90_def_var(ncid, 'salt', NF90_REAL8, dimids, salt_varid)
	status = nf90_def_var(ncid, 'sigma', NF90_REAL8, dimids, sigma_varid)
	status = nf90_def_var(ncid, 'N2', NF90_REAL8, dimids, N2_varid)

	if (wr_adv)then
		status = nf90_def_var(ncid, 'theta_adv', NF90_REAL8, dimids, theta_adv_varid)
		status = nf90_def_var(ncid, 'theta_diff', NF90_REAL8, dimids, theta_diff_varid)
		status = nf90_def_var(ncid, 'theta_for', NF90_REAL8, dimids, theta_for_varid)
		status = nf90_def_var(ncid, 'salt_adv', NF90_REAL8, dimids, salt_adv_varid)
		status = nf90_def_var(ncid, 'salt_diff', NF90_REAL8, dimids, salt_diff_varid)
		status = nf90_def_var(ncid, 'salt_for', NF90_REAL8, dimids, salt_for_varid)
	endif

	! End define mode
	status = nf90_enddef(ncid)

	! Write coordinate variable data
	status = nf90_put_var(ncid, height_varid, z)
	status = nf90_put_var(ncid, pos_varid, location)

end subroutine netcdf_setup


!----------------------------------------------------------------------------------------!
subroutine netcdf_write

	integer :: count1(2), start1(2), start2(1)

	count1 = (/ nl, 1/)
	start1 = (/ 1, nwrite /)
	start2 = (/ nwrite /)

	! Write data
	status = nf90_put_var(ncid, rec_varid, nwrite, start = start2)
	status = nf90_put_var(ncid, w_varid, w_a, start = start2 )
	status = nf90_put_var(ncid, time_varid, time, start = start2)
	status = nf90_put_var(ncid, theta_varid, theta_a, start = start1, count = count1)
	status = nf90_put_var(ncid, salt_varid, salt_a, start = start1, count = count1)
	status = nf90_put_var(ncid, sigma_varid, sigma, start = start1, count = count1)
  status = nf90_put_var(ncid, N2_varid, N2, start = start1, count = count1)

	if(wr_adv)then
		status = nf90_put_var(ncid, theta_adv_varid, theta_adv, start = start1, count = count1)
		status = nf90_put_var(ncid, theta_diff_varid, theta_diff, start = start1, count = count1)
		status = nf90_put_var(ncid, theta_for_varid, theta_for, start = start1, count = count1)
		status = nf90_put_var(ncid, salt_adv_varid, salt_adv, start = start1, count = count1)
		status = nf90_put_var(ncid, salt_diff_varid, salt_diff, start = start1, count = count1)
		status = nf90_put_var(ncid, salt_for_varid, salt_for, start = start1, count = count1)
	endif

end subroutine netcdf_write


!----------------------------------------------------------------------------------------!
subroutine netcdf_close

   status = nf90_close(ncid)
   
end subroutine netcdf_close



!----------------------------------------------------------------------------------------!
end module fluidloop_io

