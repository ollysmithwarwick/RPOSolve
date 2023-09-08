module RPOSolve
  use parameters
  use Matrix
  implicit none
  type RPONewtonParams
     ! Derived type for concisely storing and passing RPOSolve parameters
     integer          :: mgmres, nits 	 !m for gmres(m), max newton its
     real(dp)         :: rel_err	 !relative error
     real(dp)         :: del, mndl, mxdl !delta for hookstep-trust-region
     real(dp)         :: gtol, epsJ	 !gmres tolerance, eps in eqn above
     integer          :: info            !Allows for error retuns
     integer          :: ndts = -1       ! If period-like parameter exists,
                                         ! these will be recalculated each iteration
     real(dp)         :: dt = -1d0
     real(dp)         :: shortDt = 1d-5  ! Used to calculate du/dt
  end type RPONewtonParams

  interface
     ! Types for timestep and dotprod functions to be passed by user.
     function timestepFunctionType(n_, in_, ndts_, dt_)
       use parameters
       implicit none
       real(dp), dimension(:), intent(in) :: in_
       integer, intent(in) :: ndts_, n_
       real(dp), intent(in) :: dt_
       real(dp), dimension(:), allocatable :: timestepFunctionType
     end function timestepFunctionType

     double precision function metricType(n_, a, b)
       implicit none
       integer, intent(in) :: n_
       double precision, dimension(:), intent(in) :: a
       double precision, dimension(:), intent(in) :: b
     end function metricType

     subroutine saveGuessType()
       implicit none
     end subroutine saveGuessType
  end interface
  
  ! Pointers to functions passed by user (set to null until user specifies functions)
  procedure(timestepFunctionType), pointer :: timestepFunction => null ()
  procedure(saveGuessType), pointer :: saveGuess => saveGuessDefault
  procedure(metricType), pointer :: dotprod => null ()

  ! Store current operation params
  type(RPONewtonParams) :: rpoParams

  ! Labels for types of params (make enum?)
  integer :: RPO_PARAM_PERIOD_LIKE = 0
  integer :: RPO_PARAM_TRANS_LIKE  = 1
  integer :: RPO_PARAM_PHASE_LIKE  = 2
  integer :: RPO_PARAM_CONTROL_LIKE = 3 ! TODO: Make some kind of C-bind enum here?

  ! Stores info about the parameters given by user
  integer, dimension(:), allocatable, private  :: vecTypes
  logical, dimension(:), allocatable, private  :: fixedTypes
 
  ! In fixedNDTS mode (T) or fixed dt mode (F)
  logical :: fixedNDTS = .true.

contains

   subroutine assignRPORun(RPORunIn, RPORunOut)
     type(RPONewtonParams), intent(in) :: RPORunIn
     type(RPONewtonParams), intent(out) :: RPORunOut

     RPORunOut%mgmres = RPORunIn%mgmres
     RPORunOut%nits = RPORunIn%nits
     RPORunOut%rel_err = RPORunIn%rel_err
     RPORunOut%del = RPORunIn%del
     RPORunOut%mndl = RPORunIn%mndl
     RPORunOut%mxdl = RPORunIn%mxdl
     RPORunOut%gtol = RPORunIn%gtol
     RPORunOut%epsJ = RPORunIn%epsJ
     RPORunOut%info = RPORunIn%info
     RPORunOut%ndts = RPORunIn%ndts
     RPORunOut%dt = RPORunIn%dt
     RPORunOut%shortDt = RPORunIn%shortDt

   end subroutine assignRPORun

   subroutine file2RPORun(fileName, RPORunOut)
     character(len = 100), intent(in)  :: fileName
     type(RPONewtonParams) , intent(out) :: RPORunOut
     
     integer          :: mgmres, nits 	 !m for gmres(m), max newton its
     real(dp)         :: rel_err	 !relative error
     real(dp)         :: del, mndl, mxdl !delta for hookstep-trust-region
     real(dp)         :: gtol, epsJ	 !gmres tolerance, eps in eqn above
     integer          :: info            !Allows for error retuns
     integer          :: ndts = -1       ! If period-like parameter exists,
                                         ! these will be recalculated each iteration
     real(dp)         :: dt = -1d0
     real(dp)         :: shortDt = 1d-5  ! Used to calculate du/dt
          
     namelist /RPONewtonParams/ mgmres, nits, rel_err, &
          del, mndl, mxdl, gtol, epsJ, info, ndts, dt, shortDt
     
     open(12, file = fileName)
     read(12, nml = RPONewtonParams)
     close(12)

     RPORunOut%mgmres = mgmres
     RPORunOut%nits = nits
     RPORunOut%rel_err = rel_err
     RPORunOut%del = del
     RPORunOut%mndl = mndl
     RPORunOut%mxdl = mxdl
     RPORunOut%gtol = gtol
     RPORunOut%epsJ = epsJ
     RPORunOut%info = info
     RPORunOut%ndts = ndts
     RPORunOut%dt = dt
     RPORunOut%shortDt = shortDt

   end subroutine file2RPORun

   subroutine RPORun2file(RPORunIn, fileName)
     character(len = 100), intent(in)  :: fileName
     type(RPONewtonParams) , intent(out) :: RPORunIn
     
     integer          :: mgmres, nits 	 !m for gmres(m), max newton its
     real(dp)         :: rel_err	 !relative error
     real(dp)         :: del, mndl, mxdl !delta for hookstep-trust-region
     real(dp)         :: gtol, epsJ	 !gmres tolerance, eps in eqn above
     integer          :: info            !Allows for error retuns
     integer          :: ndts = -1       ! If period-like parameter exists,
                                         ! these will be recalculated each iteration
     real(dp)         :: dt = -1d0
     real(dp)         :: shortDt = 1d-5  ! Used to calculate du/dt

     namelist /RPONewtonParams/ mgmres, nits, rel_err, &
          del, mndl, mxdl, gtol, epsJ, info, ndts, dt, shortDt
     

     mgmres = RPORunIn%mgmres
     nits = RPORunIn%nits
     rel_err = RPORunIn%rel_err
     del = RPORunIn%del
     mndl = RPORunIn%mndl
     mxdl = RPORunIn%mxdl
     gtol = RPORunIn%gtol
     epsJ = RPORunIn%epsJ
     info = RPORunIn%info
     ndts = RPORunIn%ndts
     dt = RPORunIn%dt
     shortDt = RPORunIn%shortDt

     open(12, file = fileName)
     write(12, nml = RPONewtonParams)
     close(12)

   end subroutine RPORun2file


  subroutine RPONewtonSetDotProd( func )
    ! Set dotprod function
    implicit none
    procedure(metricType) :: func

    dotprod => func

    return
  end subroutine RPONewtonSetDotProd

  subroutine RPONewtonSetTimestep( func )
    ! Set timestep function
    implicit none
    procedure(timestepFunctionType) :: func
    timestepFunction => func
    return
  end subroutine RPONewtonSetTimestep

  subroutine RPONewtonSetSaveGuess( func )
    ! Set saveguess function (called once per iteration)
    implicit none
    procedure(saveGuessType) :: func

    saveGuess => func

    return
  end subroutine RPONewtonSetSaveGuess

  subroutine RPONewtonSetup(solveParams)
    ! Setup parameters for RPOSolve via RPONewtonParams derived type
    type(RPONewtonParams) :: solveParams
    rpoParams = solveParams
  end subroutine RPONewtonSetup

  subroutine RPONewtonVecSetup(vecTypes_, fixedTypes_, info)
    ! Tells RPOSolve what form parameters are
    integer, dimension(:) :: vecTypes_
    logical, dimension(:) :: fixedTypes_
    integer, intent(out) :: info

    if (allocated(vecTypes)) then
       deallocate(vecTypes)
    end if
    if (allocated(fixedTypes)) then
       deallocate(fixedTypes)
    end if

    if (size(vecTypes_) /= size(fixedTypes_) ) then
       write(*,*) 'RPONewtonInitVec: Error establishing vec types due to incompatible sizes'
       info = 1
       return
    end if

    allocate(vecTypes(size(vecTypes_)))
    allocate(fixedTypes(size(fixedTypes_)))

    vecTypes(:) = vecTypes_(:)
    fixedTypes(:) = fixedTypes_(:)

    info = 0
  end subroutine RPONewtonVecSetup

  subroutine RPONewtonClearVec()
    deallocate(vecTypes)
    deallocate(fixedTypes)
  end subroutine RPONewtonClearVec

  subroutine RPOSolveNewton(n, initGuess, info)
    use newton
    implicit none
    integer :: n
    real(dp), dimension(n) :: initGuess
   ! real(dp), dimension(:,:), allocatable :: guesses
    integer :: info
    real(dp), dimension(n) :: tmp
    real(dp) :: d, del, mndl, mxdl, tol, gtol

    
    ! Sort out ndts/dt
    if (rpoParams%ndts == -1) then
       fixedNDTS = .false.
    else
       fixedNDTS = .true.
    end if
    
    ! Allocate memory
    !if (allocated(guesses)) deallocate(guesses)
    !allocate(guesses(n,rpoParams%nits))

    if (allocated(new_x)) deallocate(new_x)
    allocate(new_x(n))

    if (allocated(new_fx)) deallocate(new_fx)
    allocate(new_fx(n))

    new_x = initGuess
    call getrhs(n,new_x, new_fx)

    ! Assign tolerances
    d = dotprod(-1, new_x, new_x) !Dotprod of phase space only

    tol  = RPOParams%rel_err * dsqrt(d)
    del  = RPOParams%del     * dsqrt(d)
    mndl = RPOParams%mndl    * dsqrt(d)
    mxdl = RPOParams%mxdl    * dsqrt(d)

    info = 1

    ! Get best parameters (LATER)

    ! Store trajectory (LATER)
    
    ! Call newtonhook
    if (RPOParams%nits > 0) then
       write(*,*) 'RPOSolve: Calling NewtonHook'
       call newtonhook(getrhs, multJ, multJp, saveGuess, dotprod, &
       RPOParams%mgmres, n, RPOParams%gtol, tol, del, mndl, mxdl, RPOParams%nits, info)
    end if
    ! deallocate
    
    initGuess = new_x
  end subroutine RPOSolveNewton

  subroutine stepOrbit(n_, x, y, ndts_, dt_)
    integer, intent(in) :: n_                         ! Size of vec
    integer, optional, intent(in) :: ndts_            ! Number of timesteps
    real(dp), optional, intent(in) :: dt_             ! Timestep size
    real(dp), dimension(:), intent(in) :: x           ! Input vec
    real(dp), dimension(:), intent(out) :: y          ! Output vec

    real(dp), dimension(:), allocatable :: currentY   ! Temporary y
    integer :: ndts
    real(dp) :: dt

    integer :: j
    real(dp) :: duration = 1.0 ! Default duration = 1.0
    
    do j = 1, size(vecTypes)
       if (vecTypes(j) == RPO_PARAM_PERIOD_LIKE) then
          duration = x(j)   !If period given in vec then use this as period 
          ! - will mean last period-like param is chosen
       end if
    end do

    !! Determine dt, ndts etc. If both given, use these, otherwise determine using period
    if (.not. present(dt_)) then
       if (present(ndts_)) then
          ndts = ndts_
       else
          ndts = 1
       end if
       dt = duration/ndts !Default to a period of 1.0
    else
       dt = dt_
       if (present(ndts_)) then
          ndts = ndts_
          duration = dt * ndts
       else
          ndts = ceiling(duration/dt)
          dt = duration/ndts
       end if
    end if

    allocate(currentY(size(x)))

    currentY = x
    currentY = timestepFunction(n_, currentY, ndts, dt) !ndts and dt can be overriden in timestep function to account for different periodicities (e.g PI model)
    y = currentY
    do j=1, size(vecTypes)
       y(j) = 0d0
    end do

  end subroutine stepOrbit
  
  subroutine getrhs(n_,x, y)
!    use orbit
    implicit none
    integer,          intent(in)  :: n_
    double precision, intent(in)  :: x(n_)
    double precision, intent(out) :: y(n_)
    double precision :: x_(n_), y_(n_)
    integer :: i

    x_ = x
    if (fixedNDTS) then
       call steporbit(n_, x_, y_, ndts_ = rpoParams%ndts)
    else
       call steporbit(n_, x_, y_, dt_ = rpoParams%dt)
    end if
    y = y_ - x
    
    do i = 1, size(vecTypes)
       y(i) = 0d0              ! One constraint per parameter
    end do
    
  end subroutine getrhs


  subroutine multJ(n_, x, y)
    use newton
    integer, intent(in)   :: n_
    real(dp), intent(in)  :: x(n_)
    real(dp), intent(out) :: y(n_)

    integer  :: i
    double precision :: ss(n_), ts(n_)
    real(dp) :: eps
    real(dp) :: shortDt 
    shortDt = rpoParams%shortDt
    
    eps = dsqrt(dotprod(1,x,x))
    if (eps==0d0) stop 'multJ: eps=0'
    eps = rpoParams%epsJ * dsqrt(dotprod(1, new_x, new_x)) / eps !eps is a proportion of the size of new_x
    if(eps==0d0) then
       stop 'multJ: eps=0'
    end if
    
    y = new_x + eps*x
    call getrhs(n_, y, ss)
    ss = (ss - new_fx)/eps
    y = new_x - eps*x
    call getrhs(n_, y, ts)
    ts = (ts - new_fx)/eps
    y = (ss - ts)/2d0

    ! No update in symmetry directions, and tracjectory direction
    do i=1, size(vecTypes)
       if (vecTypes(i) == RPO_PARAM_PERIOD_LIKE .or. vecTypes(i) == RPO_PARAM_CONTROL_LIKE) then
          call steporbit(n_, new_x, ss, 1, shortDt)
          ss = (ss - new_x)/shortDt
          y(i) = dotprod(-1, ss, x)
       else if (vecTypes(i) == RPO_PARAM_PHASE_LIKE .or. vecTypes(i) == RPO_PARAM_TRANS_LIKE) then
          ts = new_x
          ts(i) = new_x(i) + 1d-4
          call steporbit(n_, ts, ss, 0, 0d0)
          ss = (ss - new_x)/1d-4
          y(i) = dotprod(-1, ss, x)
       end if
    end do
  end subroutine multJ

 !-------------------------------------------------------------------------
 !  preconditioner for multJ.  Empty - no preconditioner required
 !-------------------------------------------------------------------------
  subroutine multJp(n, x)
    implicit none
    integer,          intent(in)    :: n
    double precision, intent(inout) :: x(n)
  end subroutine multJp

  ! double precision function dotprod(n_, a, b)
  !   implicit none
  !   integer, intent(in) :: n_
  !   double precision, dimension(:), intent(in) :: a
  !   double precision, dimension(:), intent(in) :: b
  !   dotprod = dotprodIn(n_, a, b)
  ! end function dotprod

  subroutine saveGuessDefault()
    implicit none  
  end subroutine saveGuessDefault

  subroutine RPOGetMatrix(n, ignore, m)
    integer :: n, ignore
    real(dp), dimension(n,n) :: m
    call GetMatrix(n, ignore, multJ, m)
  end subroutine RPOGetMatrix

  subroutine RPOGetEigenvalues(n, A, eigVals, eigvecs)
    integer :: n
    real(dp), dimension(n,n) :: A
    complex(dp), dimension(n) :: eigVals
    complex(dp), dimension(n,n) :: eigVecs

    write(*,*) 'RPOSolve: Getting eigenvalues, n', n
    call GetEigenvalues(n, A, eigvals, eigvecs)
!    call GetEigenvalues(n, A, eigvals)
    write(*,*) 'RPOSolve: Found eigenvalues'
  end subroutine RPOGetEigenvalues

end module RPOSolve
