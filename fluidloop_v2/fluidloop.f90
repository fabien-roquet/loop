program fluidloop
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


!USE netcdf
USE fluidloop_init
USE fluidloop_io

implicit none

real(kind=8)       :: del_t


CALL initialization
CALL open_output_files

! run iterations
do while (niter .lt. niter_max)

	del_t = 2.*dt
	if(niter == 1) del_t = dt

	call advection
	call forcing
	if (.NOT. implicit_diffusion) then
    call diffusion_explicit
  end if
  call time_step
    
  call compute_N2
  if (.NOT. constant_diffusion) then
    call Rz_mixing
  end if


	! Writing results every niter_write timesteps
	if((modulo(niter,niter_write)==0.) .OR. (niter==niter_max))then
	  nwrite = nwrite + 1
		call write_results
	endif
	
enddo 

CALL close_output_files
CALL finalization



contains


!----------------------------------------------------------------------------------------!
subroutine time_step

	! Temperature and salinity field
	if (implicit_diffusion) then
	   call diffusion_implicit
  else
	   theta_a(:) = theta_b(:)+del_t*(theta_adv(:) + theta_diff(:) + theta_for(:))
	   salt_a(:)  = salt_b(:) +del_t*(salt_adv(:)  + salt_diff(:)  + salt_for(:) )
  endif
  
	! Density and velocity
	sigma(:) = -(1.+0.5*lambda*theta_a(:)-mu*z(:))*theta_a(:) + salt_a(:)
	torque  = sum(sigma(:)*curv(:))/nl
	mass = sum(sigma(:))/nl
  w_a = tau + torque
  
  ! RAW filter
	call filter
  
	w_b = w_n
	theta_b(:) = theta_n(:)
	salt_b(:) = salt_n(:)
  
	w_n = w_a
	theta_n(:) = theta_a(:)
	salt_n(:) = salt_a(:)
	
	niter = niter + 1
  time  = niter*dt

end subroutine time_step


!----------------------------------------------------------------------------------------!
subroutine filter
	! Robert-Asselin-Williams filter

	theta_n(:) = theta_n(:) + 0.5*alpRAW*gamRAW*(theta_b(:) - 2.0*theta_n(:) + theta_a(:))
	salt_n(:) = salt_n(:) + 0.5*alpRAW*gamRAW*(salt_b(:) - 2.0*salt_n(:) + salt_a(:))
	theta_a(:) = theta_a(:) - 0.5*gamRAW*(1.0-alpRAW)*(theta_b(:) - 2.0*theta_n(:) + theta_a(:)) 
	salt_a(:) = salt_a(:) - 0.5*gamRAW*(1.0-alpRAW)*(salt_b(:) - 2.0*salt_n(:) + salt_a(:)) 

end subroutine filter


!----------------------------------------------------------------------------------------!
subroutine Rz_mixing

  integer      :: jk
  real(kind=8) :: I
    
  ! integrate calculation
  I = 0.
  do jk = 1,nl-1
	  if (N2(jk)>0) then
	    I = I + N2(jk)*dl*curv_w(jk)**2
	  end if
  end do

  do jk = 0,nl-1
    Rz(jk) = Emix/I*curv_w(jk)**2
  end do

end subroutine Rz_mixing


!----------------------------------------------------------------------------------------!
subroutine advection

  integer      :: jk
	real(kind=8) :: aux
	real(kind=8), dimension(0:nl+1) :: dummy
	
	! centered scheme
	aux = w_n / (2*dl)
	
	! temperature
	dummy(1:nl) = aux*theta_n(1:nl)
	dummy(0)    = aux*theta_n(nl)
	dummy(nl+1) = aux*theta_n(1)
	do jk = 1,nl
	   theta_adv(jk) = dummy(jk-1) - dummy(jk+1)
	enddo
	
	! salinity
	dummy(1:nl) = aux*salt_n(1:nl)
	dummy(0)    = aux*salt_n(nl)
	dummy(nl+1) = aux*salt_n(1)
	do jk = 1,nl
	   salt_adv(jk)  = dummy(jk-1) - dummy(jk+1)
	enddo

end subroutine advection


!----------------------------------------------------------------------------------------!
subroutine forcing

  if(implicit_diffusion)then
	  theta_for(jsource) =  nl*xi_t*ref_t
	  theta_for(jsink)   = -nl*xi_t*ref_t
	else
	  theta_for(jsource) =  nl*xi_t*(ref_t-theta_b(jsource))
	  theta_for(jsink)   = -nl*xi_t*(ref_t+theta_b(jsink)  )
  end if
  
	salt_for(jsource)  =  nl*F_s
	salt_for(jsink)    = -nl*F_s

end subroutine forcing


!----------------------------------------------------------------------------------------!
subroutine diffusion_explicit

  integer            :: jk
	real(kind=8), dimension(0:nl-1) :: aux
	real(kind=8), dimension(0:nl) :: dummy

  aux(0:nl-1) = (Rz(0:nl-1)+Rh(0:nl-1)) / (dl*dl)

	! theta
	dummy(0) = aux(0) * (theta_b(1)-theta_b(nl))
	do jk = 1,nl-1
	  dummy(jk) = aux(jk) * (theta_b(jk+1)-theta_b(jk))
	end do
	dummy(nl) = dummy(0)
	do jk = 1,nl
	  theta_diff(jk) = dummy(jk) - dummy(jk-1)
	end do	
	
	! salt
	dummy(0) = aux(0) * (salt_b(1)-salt_b(nl))
	do jk = 1,nl-1
	  dummy(jk) = aux(jk) * (salt_b(jk+1)-salt_b(jk))
	end do
	dummy(nl) = dummy(0)
	do jk = 1,nl
	  salt_diff(jk) = dummy(jk) - dummy(jk-1)
	end do	
  
end subroutine diffusion_explicit


!----------------------------------------------------------------------------------------!
subroutine diffusion_implicit

  integer            :: jk
  real(kind=8),dimension(1:nl)  :: zwi, zwd, zws, zwy
  real(kind=8),dimension(1:nl)  :: zwx, zwz
  real(kind=8),dimension(1:nl)  :: y, q
  real(kind=8),dimension(0:nl-1):: aux
  real(kind=8)                  :: u1, uN, vN

  aux(:) = del_t*(Rz(:)+Rh(:))/(dl*dl)
  zwi(1) = -aux(0)
  zws(1) = -aux(1)
  do jk = 2,nl-1
	  zwi(jk) = -aux(jk-1)
	  zws(jk) = -aux(jk)
	end do
  zwi(nl) = -aux(nl-1)
  zws(nl) = -aux(0)
  zwd(:) = 1. - zwi(:) - zws(:)
	
	! for temperature only (forcing)
  zwd(jsource) = zwd(jsource) + del_t*nl*xi_t
  zwd(jsink  ) = zwd(jsink  ) + del_t*nl*xi_t
  
  ! modified matrix (Sherman-Morrison formula)
  u1 = -zwd(1)
  uN = zws(nl)
  vN = zwi(1) / u1
	zwd(1) = 2.*zwd(1)
	zwd(nl)= zwd(nl) - uN * vN
	zwi(1) = 0.
	zws(nl)= 0.
	
  zwy(:) = theta_b(:)+del_t*(theta_adv(:) + theta_diff(:) + theta_for(:))
	y = solve_tridiag(zwi, zwd, zws, zwy)
  zwy(:) = 0.
  zwy(1) = u1
  zwy(nl)= uN
	q = solve_tridiag(zwi, zwd, zws, zwy)
	theta_a(:) = y(:) - q(:) * (y(1)+vN*y(nl))/(1.+q(1)+vN*q(nl))
	
  aux(:) = del_t*(Rz(:)+Rh(:))/(dl*dl)
  zwi(1) = -aux(0)
  zws(1) = -aux(1)
  do jk = 2,nl-1
	  zwi(jk) = -aux(jk-1)
	  zws(jk) = -aux(jk)
	end do
  zwi(nl) = -aux(nl-1)
  zws(nl) = -aux(0)
  zwd(:) = 1. - zwi(:) - zws(:)
	
  ! modified matrix (Sherman-Morrison formula)
  u1 = -zwd(1)
  uN = zws(nl)
  vN = zwi(1) / u1
	zwd(1) = 2.*zwd(1)
	zwd(nl)= zwd(nl)- uN * vN
	zwi(1) = 0.
	zws(nl)= 0.
	
  zwy(:) = salt_b(:)+del_t*(salt_adv(:) + salt_diff(:) + salt_for(:))
	y = solve_tridiag(zwi, zwd, zws, zwy)
  zwy(:) = 0.
  zwy(1) = u1
  zwy(nl)= uN
	q = solve_tridiag(zwi, zwd, zws, zwy)
	salt_a(:) = y(:) - q(:) * (y(1)+vN*y(nl))/(1.+q(1)+vN*q(nl))
	
end subroutine diffusion_implicit


function solve_tridiag(zwi, zwd, zws, zwy) result(zwx)
  !! Matrix inversion from the first level
  !!----------------------------------------------------------------------
  !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
  !
  !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
  !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
  !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
  !        (        ...               )( ...  ) ( ...  )
  !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
  !
  !   m is decomposed in the product of an upper and lower triangular matrix.
  !   The 3 diagonal terms are in 1d arrays: zwd, zws, zwi.
  !   Suffices i,s and d indicate "inferior" (below diagonal), diagonal
  !   and "superior" (above diagonal) components of the tridiagonal system.
  integer                                  :: jk
  real(kind=8),dimension(1:nl),intent(in)  :: zwi, zwd, zws, zwy
  real(kind=8),dimension(1:nl)             :: zwx
  real(kind=8),dimension(1:nl)             :: zwz,zwt

  ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
  zwt(1) = zwd(1)
  DO jk = 2, nl
    zwt(jk) = zwd(jk) - zwi(jk) * zws(jk-1) / zwt(jk-1)
  END DO
  ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1 (using zwd memory to store Z)
  zwz(1) = zwy(1)
  DO jk = 2, nl
    zwz(jk) = zwy(jk) - zwi(jk) / zwt(jk-1) * zwz(jk-1)
  END DO
  ! third recurrence:    Xk = (Zk - Sk Xk+1 ) / Tk
  zwx(nl) = zwz(nl) / zwt(nl)
  DO jk = nl-1, 1, -1
    zwx(jk) = ( zwz(jk) - zws(jk) * zwx(jk+1) ) / zwt(jk)
  END DO

 end function solve_tridiag
 
 
!----------------------------------------------------------------------------------------!
end program

