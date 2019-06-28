module fluidloop_init
implicit none
! module fluidloop_init
!
! loop toolbox, distributed on GitHub: http://github.com/fabien-roquet/loop
! F. Roquet 2016
! GNU General Public License


PRIVATE allocation
PRIVATE mixing_energy
PRIVATE control_stability
PRIVATE deallocation

PUBLIC  initialization
PUBLIC  compute_N2
PUBLIC  finalization



! Declaration parameters namelist namrun
character (len=50) :: name_exp
integer            :: nl, niter_t1, niter_max, niter_write
logical            :: use_netcdf,init_from_file
character (len=50) :: init_state

! Declaration parameters namelist namconfig
real(kind=8)       :: Zf, R, lambda, mu
logical            :: foldtrue, constant_diffusion, implicit_diffusion

! Declaration parameters namelist namforcing
real(kind=8)       :: tau, ref_t, xi_t, F_s

! Running parameters
integer            :: niter, nwrite
real(kind=8)       :: pi, dt, time

real(kind=8)       :: w_b, w_a, w_n, Emix
real(kind=8)       :: dl, torque, mass
real(kind=8)       :: gamRAW, alpRAW
integer            :: jsink, jsource
  
! Reading parameters
namelist / namrun / name_exp, nl, niter_t1, niter_max, niter_write,  &
  &    use_netcdf,init_from_file, init_state
namelist / namconfig / Zf, R,lambda,mu,foldtrue,constant_diffusion,  &
  &  implicit_diffusion
namelist / namforcing / tau,ref_t,xi_t,F_s

real(kind=8), dimension(:), allocatable :: z, location, sigma, curv, curv_w, mfold
real(kind=8), dimension(:), allocatable :: theta_b,theta_n,theta_a,salt_b,salt_n,salt_a
real(kind=8), dimension(:), allocatable :: theta_diff, theta_adv, theta_for
real(kind=8), dimension(:), allocatable :: salt_diff, salt_adv, salt_for
real(kind=8), dimension(:), allocatable :: Rz, Rh, N2




contains


!----------------------------------------------------------------------------------------!
subroutine initialization
implicit none

  integer            :: jk, ji
  character (len=50) :: name_out_0

  ! 1) read namelist
  open (unit=33,file='namelist',action='read',status='old')
  read (33,namrun)
  read (33,namconfig)
  read (33,namforcing)
  close(unit=33)

  ! allocate memory
  CALL allocation  
  
  ! open output file
  write(name_out_0,'(a,"_config.txt")') trim(name_exp)
  open(unit=10,file=name_out_0,action='write')  
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

  ! Initialization
  pi              = 4.0*atan(1.0)
  dt              = 1. / niter_t1
	dl   = 2*pi/nl
	w_b = 0.
	w_n = 0.
	w_a = 0.
 
  if (MOD(nl,2) .EQ. 1) then
    write(10,*) 'NL must be an even number: nl=', nl
    stop
	endif
  
  if(Zf .LT. -1 .OR. Zf .GT. 1)then
     write(10,*) 'Source/sink height Zf must be a number between -1 and 1: Zf=',Zf
     stop
  endif

	! Forcing position source and sink
	location = (/ (jk, jk = 1, nl) /)
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
	mass = 0.
	torque = 0.
  
  ! diffusivities
  do jk = 0,nl-1
     Rh(jk) = R*cos(dl*(jk+0.5))**2.
     curv_w(jk) = -sin(dl*(jk+0.5))
     Rz(jk) = R*curv_w(jk)**2.
  end do
  
	if(foldtrue)then
	   Rh(0:jsink-1) = R
	   Rh(jsource:nl-1) = R
     curv_w(0:jsink-1) = 0.
	   curv_w(jsource:nl-1) = 0.
	   Rz(0:jsink-1) = 0.
	   Rz(jsource:nl-1) = 0.
  end if
  
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
	
  ! compute total energy of mixing Emix in initial state
  call mixing_energy
  
  ! Parameters of RAW filter for the leapfrog
  gamRAW = 0.1                      ! first Robert-Asselin-Williams filter constant
  alpRAW = 0.53                     ! second Robert-Asselin-Williams filter constant

  ! Stability check
  if (.NOT. implicit_diffusion) then
     call control_stability
  end if
  
  ! count number of iteration
  niter     = 0
  nwrite    = 0
  time  = niter*dt
  
  ! output
  write(10,*) 'Initialization complete'
	write(10,*) 'w_n init =',w_n
	write(10,*) 'Emix  init =',Emix
	
  write (unit=10,fmt=*) ''

end subroutine initialization

!----------------------------------------------------------------------------------------!
subroutine mixing_energy

  integer      :: jk
	
  call compute_N2
  Emix = 0.
  do jk = 1,nl-1
	  if (N2(jk)>0) then
	    Emix = Emix + Rz(jk)*N2(jk)*dl
	  end if
  end do

end subroutine mixing_energy


!----------------------------------------------------------------------------------------!
subroutine compute_N2

  integer      :: jk
	   
  sigma(:) = -(1.+0.5*lambda*theta_n(:)-mu*z(:))*theta_n(:) + salt_n(:)
  N2(0) = - curv_w(0)*(sigma(1)-sigma(nl))/dl
  do jk = 1,nl-1
    N2(jk) = - curv_w(jk)*(sigma(jk+1)-sigma(jk))/dl
  end do
    
end subroutine compute_N2


!----------------------------------------------------------------------------------------!
subroutine finalization

  character (len=50) :: name_out
  integer            :: jk

  write(10,*) 'last velocity', w_n, 'at time=', time
  close(unit=10)
  
  ! write restart
  write(name_out,'(a,"_restart.txt")') trim(name_exp)
  open (unit=34,file=name_out,action='write',status='new')
  write (unit=34,fmt="(f15.7)") w_n
  do jk = 1,nl
    write (unit=34,fmt="(i10, 2(f15.7))") jk,theta_n(jk),salt_n(jk)
  end do
  close (unit=34)

  CALL deallocation

end subroutine finalization


!----------------------------------------------------------------------------------------!
subroutine control_stability

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

end subroutine control_stability


!----------------------------------------------------------------------------------------!
subroutine allocation

	allocate(theta_a(1:nl),theta_n(1:nl),theta_b(1:nl))
	allocate(salt_a(1:nl),salt_n(1:nl),salt_b(1:nl))
	allocate(theta_diff(1:nl),theta_adv(1:nl),theta_for(1:nl))
	allocate(salt_diff(1:nl),salt_adv(1:nl),salt_for(1:nl))
	allocate(z(1:nl),location(1:nl),mfold(1:nl))
	allocate(sigma(1:nl),curv(1:nl),curv_w(0:nl-1))
	allocate(N2(0:nl-1),Rz(0:nl-1),Rh(0:nl-1))

  ! initialize tables
	theta_a(:) = 0.
	theta_n(:) = 0.
	theta_b(:) = 0.
	salt_b(:) = 0.
	salt_n(:) = 0.
	salt_a(:) = 0.
	
	theta_diff(:) = 0.
	theta_adv(:) = 0.
	theta_for(:) = 0.
	salt_diff(:) = 0.
	salt_adv(:) = 0.
	salt_for(:) = 0.
	
	z(:) = 0.
	location(:) = 0.
	mfold(:) =0.
	sigma(:) = 0.
	curv(:) = 0.
	curv_w(:) = 0.
	N2(:) = 0.
	Rz(:) = 0.
	Rh(:) = 0.
		
end subroutine allocation


!----------------------------------------------------------------------------------------!
subroutine deallocation

	deallocate(theta_a,theta_b,theta_n, salt_a, salt_b, salt_n)
	deallocate(theta_diff, theta_adv, salt_diff, salt_adv, theta_for, salt_for)
	deallocate(z,location,mfold,sigma,curv,curv_w,N2,Rz,Rh)

end subroutine deallocation


!----------------------------------------------------------------------------------------!
end module fluidloop_init

