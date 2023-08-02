module lorenz
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp) :: sigma, beta, rho
  real(dp), dimension(4) :: current
  contains

  function dxdt(x)
    real(dp), dimension(3) :: x
    real(dp), dimension(3) :: dxdt
    dxdt = 0d0

    dxdt(1) = sigma*(x(2)-x(1))
    dxdt(2) = x(1)*(rho-x(3)) - x(2)
    dxdt(3) = x(1)*x(2) - beta * x(3)

  end function dxdt

  subroutine timestepLorenz(in, out, dt)
    real(dp), dimension(3) :: in, out
    real(dp) :: dt

    out = in + dt*dxdt(in)
  end subroutine timestepLorenz
  
end module lorenz


program lorenzRPO
  use lorenz
  use RPOSolve, only : RPONewtonSetTimestep, RPONewtonSetDotprod, RPONewtonSetSaveGuess, stepOrbit, RPONewtonVecSetup, RPONewtonParams, RPONewtonSetup, RPOSolveNewton, dotprod
  integer :: info
  type(RPONewtonParams) :: rpoParams
  real(dp), dimension(:,:), allocatable :: guesses
  real(dp), dimension(:,:), allocatable :: states
  integer :: spf


  ! Lorenz model parameters
  sigma = 10d0
  beta  = 8d0/3d0
  rho   = 28d0

  ! Initial guess
  current = (/3.084d0, -10d0,-20.0d0,31d0 /)

  ! RPO Solver parameters encapsulated in RPONEWTONPARAMS derived type (RPOSolve Module)
  rpoParams%mgmres = 4          ! Iterations of GMRES per Newton iteration
  rpoParams%nits   = 100         ! Newton Iteration max
  rpoParams%rel_err   = 1d-7     ! Target relative error for residual
  rpoParams%del    = -1d0        ! Initial hookstep size (-1 for automatic calculation)
  rpoParams%mndl   = 1d-14       ! Minimum hookstep size. Program currently terminates if hookstep gets too small (TODO: Deal with this better)
  rpoParams%mxdl   = 1d+14       ! Maximum hookstep size
  rpoParams%gtol   = 1d-3        ! GMRES tolerance
  rpoParams%epsJ   = 1d-5        ! epsilon used in Jacobian calculation
  rpoParams%info   = 0           ! Used for output - currentky unused on input
  rpoParams%ndts   = 10000000    ! Number of timesteps. Specifying this keeps this (roughly) fixed and adjusts timestep size accordingly.
                                 ! Could also specify rpoParams%dt in which the inverse is true

  spf = 100000 ! Number of timesteps per output frame

  allocate(guesses(4,rpoParams%nits+1))       ! Guesses gives the initial guess at every newton iteration (including 0th and last)
  allocate(states(4,rpoParams%ndts/spf + 1)) ! States stores the best solution


  call RPONewtonSetTimestep( timestepLorenzRPO ) ! Tells RPO solver what the timestepping function is. Must be of the form: function
                                                 ! function name(n, in, ndts, dt)
                                                 ! where name      returns an array of size n containing final state
                                                 !       n         is the dimension of the system (including parameters)
                                                 !       in        is the initial guess
                                                 !       ndts      number of timesteps
                                                 !       dt        size of timestep
  call RPONewtonSetDotprod( lorenzDotprod )      ! Tells RPO solver what to use as dotprod function
                                                 ! used for magnitude and projection calculations. Takes an argument n which, if -1 only includes the real space components
                                                 ! (in this case if n==-1 the period, in(0), is not included
  call RPONewtonSetSaveGuess( lorenzSaveGuess )  ! Tells RPO solver what to call every iteration to save the current guess (might rework how this works)

  call RPONewtonVecSetup((/ 0 /),(/ .false./), info)  ! Tells RPO solver what properties the parameters have
                                                      ! First argument is array of types
                                                      ! 0    Period Type
                                                      ! 1    Translation length like
                                                      ! 2    Phase Like
                                                      ! 3    Control like (e.g if you wanted a model parameter to be solved for)
  call RPONewtonSetup(rpoParams)                      ! Setup function for RPOSolve - pass RPOParams
  call RPOSolveNewton(4,current, info)                ! Tells RPOSolve to run
  current = timestepLorenzRPOSave(4, current, rpoParams%ndts, current(1)/(rpoParams%ndts), states_=states, spf=spf)
  call writeFile("./guesses.dat", guesses)
  call writeFile("./soln.dat", states)

contains

 
  function timestepLorenzRPO(n, in, ndts, dt)
    implicit none
    real(dp), dimension(:), intent(in) :: in
    integer, intent(in) :: ndts, n
    real(dp), intent(in) :: dt
    real(dp), dimension(:), allocatable :: timestepLorenzRPO

    integer :: j
    real(dp), dimension(3) :: current
    allocate(timestepLorenzRPO(n))
    timestepLorenzRPO = 0d0
    
    current(:) = in(2:)
    do j = 1, ndts
       call timestepLorenz(current, current, dt)
    end do

    timestepLorenzRPO(2:) = current(:)
  end function timestepLorenzRPO

  function timestepLorenzRPOSave(n, in, ndts, dt, states_, spf)
    implicit none
    real(dp), dimension(:), intent(in) :: in
    integer, intent(in) :: ndts, n
    real(dp), intent(in) :: dt
    real(dp) :: time
    real(dp), dimension(:), allocatable :: timestepLorenzRPOSave
    integer :: spf
    real(dp), dimension(4,ndts/spf + 1) :: states_


    integer :: j
    real(dp), dimension(3) :: current
    allocate(timestepLorenzRPOSave(n))
    timestepLorenzRPOSave = 0d0
    current(:) = in(2:)
   
    states_(1 ,1) = 0d0
    states_(2:,1) = current(:)
    do j = 1, ndts
       time = time + dt
       call timestepLorenz(current, current, dt)
       if (mod(j,spf)==0) then

          states_(1, j/spf + 1) = time
          states_(2:,j/spf + 1) = current(:)
       end if
    end do
!    write(*,*) states_
    timestepLorenzRPOSave(2:) = current(:)
  end function timestepLorenzRPOSave
  
  double precision function lorenzDotprod(n_,a,b)
    implicit none
    integer,          intent(in) :: n_
    double precision,       dimension(:), intent(in) ::  a
    double precision,       dimension(:), intent(in) ::  b
    !real(dp),       dimension(:) :: a1(4), b1(4)
    !   double precision :: d,d_
    integer :: n1

    n1 = 1
    !a1 = a
    !b1 = b
    if (n_==-1) n1 = 2
    lorenzDotprod = dble(dot_product(a(n1:),b(n1:)))

  end function lorenzDotprod

  subroutine lorenzSaveGuess()
    use newton, only : new_x, new_nits ! Ideally want to do this without depending on newton module
    implicit none
    
    guesses(:,new_nits) = new_x
    
  end subroutine lorenzSaveGuess

  subroutine writeFile(fileName, mat)
    character(*), intent(in) :: fileName
    character(len = 21) :: outfmt
    real(dp), intent(in), dimension(1:, 1:) :: mat
    integer, dimension(2) :: shp
    integer :: cols, rows
    integer :: j

    shp = shape(mat)
    cols = shp(1)
    rows = shp(2)
    open(12, file = fileName)
    write(outfmt, '(A1, I9, A11)') '(', cols, '(E30.15E3))'

    do j = 1, rows
       write(12, fmt=outfmt) mat(:, j)
    end do

    close(12)

  end subroutine writeFile
  
end program lorenzRPO
