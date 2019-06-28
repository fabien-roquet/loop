program fluidloop
implicit none
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

! Declaration parameters namelist namrun
character (len=50) :: name_exp
integer            :: nl, niter_t1, niter_max, niter_write_0D, niter_write_1D
logical            :: init_from_file
character (len=50) :: init_state

! Declaration parameters namelist namconfig
real(kind=8)       :: Zf, R, lambda, mu
logical            :: foldtrue

! Declaration parameters namelist namforcing
real(kind=8)       :: tau, ref_t, xi_t, F_s

! Running parameters
integer            :: niter, nwrite
real(kind=8)       :: pi, dt, del_t

real(kind=8)       :: w_b, w_a, w_n
real(kind=8)       :: dl, torque, mass, time
real(kind=8)       :: gamRAW, alpRAW
integer            :: jsink, jsource, jk, ji
character (len=50) :: name_out_0, name_out_1, name_out_2
  
real(kind=8), dimension(:), allocatable :: z, location, sigma, curv, mfold
real(kind=8), dimension(:), allocatable :: theta_b,theta_n,theta_a,salt_b,salt_n,salt_a
real(kind=8), dimension(:), allocatable :: theta_diff, theta_adv, theta_for
real(kind=8), dimension(:), allocatable :: salt_diff, salt_adv, salt_for

! Reading parameters
namelist / namrun / name_exp, nl, niter_t1, niter_max, niter_write_0D, niter_write_1D,  &
  & init_from_file, init_state
namelist / namconfig / Zf,R,lambda,mu,foldtrue
namelist / namforcing / tau,ref_t,xi_t,F_s

! 1) read namelist
open (unit=33,file='namelist',action='read',status='old')
read (33,namrun)
read (33,namconfig)
read (33,namforcing)
close(unit=33)

! open output files

write(name_out_0,'(a,"_config.txt")') trim(name_exp)
open(unit=10,file=name_out_0,action='write',status='new')

write(name_out_1,'(a,"_out_0D.txt")') trim(name_exp)
open(unit=11,file=name_out_1,action='write',status='new')
write(11,*) '% niter time velocity mass'

write(name_out_2,'(a,"_out_1D.txt")') trim(name_exp)
open(unit=12,file=name_out_2,action='write',status='new')
write(12,*) '% niter j theta salt sigma'

! Initialization
pi              = 4.0*atan(1.0)
dt              = 1. / niter_t1
dt = 0.0001
call allocation  
call initialization

! run iterations
do while (niter .le. niter_max)

	del_t = 2.*dt
	if(niter == 1) del_t = dt

	call advection
	call diffusion
	call forcing
	call time_step

	! Writing results every niter_write timesteps
	if((modulo(niter,niter_write_0D)==0.) .OR. (niter==niter_max))then
		call write_results_0D
	  nwrite = nwrite + 1
	endif
	
	if((modulo(niter,niter_write_1D)==0.) .OR. (niter==niter_max))then 
		call write_results_1D
	endif
	
enddo 

! finalize program
print*, 'last velocity', w_n, 'at time=', time
write(10,*) 'last velocity', w_n, 'at time=', time

close(unit=10)
close(unit=11)
close(unit=12)

call deallocation






contains


!----------------------------------------------------------------------------------------!
subroutine initialization

  write (unit=10,fmt=*) 'fluidloop.f90, version 1.0'
  write (unit=10,fmt=*) '    available on GitHub: http://github.com/fabien-roquet/loop'
  write (unit=10,fmt=*) '    F. Roquet 2016'
  write (unit=10,fmt=*) '    GNU General Public License'
  write (unit=10,fmt=*) ''
  write (unit=10,fmt=*) ''
  write (unit=10,fmt=*) 'Namelist parameters:'
  write (10,namrun)
  write (10,namconfig)
  write (10,namforcing)
  write (unit=10,fmt=*) ''

	dl   = 2*pi/nl
	location = (/ (jk, jk = 1, nl) /)
  
  if (MOD(nl,2) .EQ. 1) then
    write(10,*) 'NL must be an even number: nl=', nl
    stop
	endif
  
  if(Zf .LT. -1 .OR. Zf .GT. 1)then
     write(10,*) 'Source/sink height Zf must be a number between -1 and 1: Zf=',Zf
     stop
  endif

	! Forcing position source and sink
	z(:) = cos(dl*location(:))
	jsink  = 0
  do while ((z(jsink+1) - Zf > 1.e-10) .AND. (jsink < nl))
    jsink = jsink + 1
  end do
	jsource = nl - jsink
  write(10,*) 'Loop index of sink  : jsink  =', jsink
  write(10,*) 'Loop index of source: jsource=', jsource
  write(10,*) 'Actual height of source/sink: z(jsink)=', z(jsink)

	! Fold mask
	mfold(:)=1.0
	if(foldtrue)then
		mfold(1:jsink)=0.
		mfold(jsource:nl)=0.
	endif

	! Height and curvature
	z(:) = z(:)*mfold(:)+(1.0-mfold(:))*z(jsink)
  curv(:)  = sin(dl*location(:))*mfold(:)
  
  ! initialize tables
	theta_b(:) = 0.
	salt_b(:) = 0.
	w_b = 0.
	theta_n(:) = 0.
	salt_n(:) = 0.
	w_n = 0.
	theta_a(:) = 0.
	salt_a(:) = 0.
	w_a = 0.

	sigma(:) = 0.
	mass = 0.
	torque = 0.

	theta_diff(:) = 0.
	theta_adv(:) = 0.
	theta_for(:) = 0.
	salt_diff(:) = 0.
	salt_adv(:) = 0.
	salt_for(:) = 0.

  ! Start the model from input file
  if (init_from_file) then
    open (unit=34,file=init_state,action='read',status='old')
  	read (unit=34,fmt=*) w_n
    do jk = 1,nl
	    read (unit=34,fmt=*) ji,theta_n(jk),salt_n(jk)
	  end do
	  close (unit=34)
    w_b = w_n
    theta_b = theta_n
    salt_b = salt_n
	endif
	write(10,*) 'w_n init =',w_n
	
  ! Parameters of RAW filter for the leapfrog
  gamRAW = 0.1                      ! first Robert-Asselin-Williams filter constant
  alpRAW = 0.53                     ! second Robert-Asselin-Williams filter constant

  ! Stability check
  call stability

  ! count number of iteration
  niter = 1
  nwrite= 1
  time  = niter*dt
  
  ! output
  write(10,*) 'Initialization complete'
  write (unit=10,fmt=*) ''



end subroutine initialization


!----------------------------------------------------------------------------------------!
subroutine stability

	real(kind=8) :: sc

	! Stability criterion from pure diffusion equation: 4 R dt/(dl**2) .le. 1
	sc = 4.*R*dt/(dl**2.)
	if (sc >= 1.)then
		write(10,*) 'Diffusion scheme unstable!'
		write(10,*) 'Stability criterion: 4 R dt/(dl**2) .le.1'
		write(10,*) 'Factor is', sc,'; adjust R, dt or dl'
		stop
	else
		write(10,*) 'Diffusion scheme stable'
		write(10,*) 'Stability factor is',sc,'(must be less than 1, should be less than 0.5)'
	endif

	! Weak stability criterion for temperature relaxation xi*nl*dt .le. 1
  sc = xi_t*nl*dt
	if (xi_t .GT. 0) then
    if (sc >= 0.90)then
      write(10,*) 'Temperature relaxation unstable!'
      write(10,*) 'Stability criterion: xi*nl*dt .le. 0.90'
      write(10,*) 'Factor is', sc,'; adjust xi, dt or nl'
      stop
    else
      write(10,*) 'Temperature relaxation stable'
      write(10,*) 'Stability factor is',sc,'(must be less than 0.90)'
    endif
	else if (xi_t .LT. 0) then
		write(10,*) 'Temperature Restoring constant must be positive'
		stop
  endif

end subroutine stability


!----------------------------------------------------------------------------------------!
subroutine allocation

	allocate(theta_a(1:nl),theta_n(1:nl),theta_b(1:nl))
	allocate(salt_a(1:nl),salt_n(1:nl),salt_b(1:nl))
	allocate(salt_diff(1:nl),salt_adv(1:nl),salt_for(1:nl))
	allocate(theta_diff(1:nl),theta_adv(1:nl),theta_for(1:nl))
	allocate(z(1:nl),location(1:nl),mfold(1:nl))
	allocate(sigma(1:nl),curv(1:nl))
	
end subroutine allocation


!----------------------------------------------------------------------------------------!
subroutine deallocation

	deallocate(theta_a,theta_b,theta_n, salt_a, salt_b, salt_n)
	deallocate(theta_diff, theta_adv, salt_diff, salt_adv, theta_for, salt_for)
	deallocate(sigma,curv,z,location,mfold)

end subroutine deallocation


!----------------------------------------------------------------------------------------!
subroutine advection

	real(kind=8) :: aux
	real(kind=8), dimension(0:nl+1) :: dummy_s, dummy_t

	aux = 1./(2.*dl)
	
	dummy_t(1:nl) = theta_n(:)
	dummy_t(0)    = theta_n(nl)
	dummy_t(nl+1) = theta_n(1)
	
	dummy_s(1:nl) = salt_n(:)
	dummy_s(0)    = salt_n(nl)
	dummy_s(nl+1) = salt_n(1)
	
	do jk = 1,nl
		theta_adv(jk) = -w_n*aux*(dummy_t(jk+1) - dummy_t(jk-1))
		salt_adv(jk)  = -w_n*aux*(dummy_s(jk+1) - dummy_s(jk-1))
	enddo

end subroutine advection


!----------------------------------------------------------------------------------------!
subroutine diffusion

	real(kind=8), dimension(0:nl) :: dummy

	dummy(0) = R*(theta_b(1)-theta_b(nl))/(dl*dl)
  do jk = 1,nl-1
	  dummy(jk) = R*(theta_b(jk+1)-theta_b(jk))/(dl*dl)
	end do
	dummy(nl) = dummy(0)
	do jk = 1,nl
	  theta_diff(jk) = dummy(jk)-dummy(jk-1)
  end do	
	
	dummy(0) = R*(salt_b(1)-salt_b(nl))/(dl*dl)
  do jk = 1,nl-1
	  dummy(jk) = R*(salt_b(jk+1)-salt_b(jk))/(dl*dl)
	end do
	dummy(nl) = dummy(0)
	do jk = 1,nl
	  salt_diff(jk) = dummy(jk)-dummy(jk-1)
  end do	

end subroutine diffusion


!----------------------------------------------------------------------------------------!
subroutine forcing

	theta_for(jsource) =  nl*xi_t*(ref_t-theta_b(jsource))
	theta_for(jsink)   = -nl*xi_t*(ref_t+theta_b(jsink)  )
 
	salt_for(jsource)  =  nl*F_s
	salt_for(jsink)    = -nl*F_s

end subroutine forcing


!----------------------------------------------------------------------------------------!
subroutine time_step

	! Temperature and salinity field
	theta_a(:) = theta_b(:)+del_t*(theta_adv(:) + theta_diff(:) + theta_for(:))
	salt_a(:)  = salt_b(:) +del_t*(salt_adv(:)  + salt_diff(:)  + salt_for(:) )
  
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
subroutine write_results_0D

	! Write data 0D
	write (unit=11,fmt="(i10,3(f10.4))") niter,time,w_n,mass

end subroutine write_results_0D


!----------------------------------------------------------------------------------------!
subroutine write_results_1D

	! Write data 1D
  do jk = 1,nl
	  write (unit=12,fmt="(2(i10),3(f10.3))") niter,jk,theta_n(jk),salt_n(jk),sigma(jk)
	end do

end subroutine write_results_1D


!----------------------------------------------------------------------------------------!
end program

