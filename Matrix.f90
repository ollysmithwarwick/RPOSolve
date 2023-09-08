module matrix
  implicit none
  ! ------------------------------------------
  ! Finds the matrix A representing a function matvec
  ! ------------------------------------------
contains
  subroutine GetMatrix(n, ignore, matvec, A)
    use parameters
    implicit none

    integer,          intent(in)  :: n, ignore
    double precision              :: x(n)
    double precision              :: y1(n), y2(n)
    external                      :: matvec
    double precision, intent(out) :: A(n - ignore, n - ignore)
    integer :: i1

    do i1 = (ignore + 1), n
       x(:) = 0.0_dp
       x(i1) = 1.0_dp
       call matvec(n, x, y1)
       call matvec(n,-x, y2)
       A(1:n-ignore, i1-ignore) = (y1(ignore + 1: n) - y2(ignore+1: n))/2.0_dp
    end do
    
    do i1 = 1, n-ignore
       A(i1, i1) = A(i1,i1)  + 1.0_dp
    end do
  end subroutine GetMatrix

  ! subroutine GetSym(sym, run0)
  !   use runType
  !   double precision, intent(out) :: sym(1:)
  !   double precision, dimension(:), allocatable :: symL, x
  !   type(POPIRun), intent(in) :: RPORun0
  !   type(run), intent(in) :: run0

  !   external steporbit
    
  !   allocate(symL(size(sym)+ 3))
  !   allocate(   x(size(sym)+ 3))
  !   write(*,*) 'GetSym'
  !   x(1)  = RPORun0%SGuess
  !   x(2)  = RPORun0%ShiftGuess
  !   x(3)  = RPORun0%PhaseGuess
  !   x(4:) = RPORun0%initGuess
    
  !   call steporbit(n, 0, 1d-2, 0d0, run0, x, symL)
  !   sym = (symL(4:) - new_x(4:))/(1d-2)
  ! end subroutine GetSym
  
  subroutine GetEigenvalues(n, A, eigVals, eigvecs)
    use parameters
    
    integer, intent(in) :: n
    integer :: lda, ldvl, ldvr
    integer, parameter :: lwmax = 10000
    integer :: info, lwork
!    logical, optional :: ignoreSym
    logical :: ignore
    logical :: returnVecs

    complex(dp) :: eigVals(:)
    complex(dp), optional :: eigvecs(:,:)
    complex(dp), dimension(:,:), allocatable :: eigVecs_
    real(dp) :: A(1:, 1:)
    double precision, allocatable :: mat(:,:)
!    double precision, intent(in), optional :: sym(1:,1:)
!    real(dp), allocatable :: sym_(:,:)
    double precision, allocatable :: vl(:,:), vr(:,:), wr(:), wi(:), work(:), coss(:,:)
    double precision :: rpart(n), ipart(n)

    integer :: i1, i2, j
    integer :: shp(2)

    external dgeev

    write(*,*) 'GetEigenvalues: Start'
       

    ! if(present(ignoreSym).and.present(sym)) then
    !    shp = shape(sym)
    !    ignore = ignoreSym
    !    allocate(sym_(shp(1),shp(2)))
    !    sym_ = sym
    ! else
    !    ignore = .false.
    ! end if

    if (present(eigvecs)) then
       returnVecs=.true.
    end if
    allocate(eigvecs_(n,n))
       
    lda=N
    ldvl=N
    ldvr=N
    lwork = -1

    shp=shape(A)
    allocate(vl(ldvl,n))
    allocate(vr(ldvr,n))
    allocate(wr(n))
    allocate(wi(n))
    allocate(work(lwmax))
    allocate(coss(4, n))
    allocate(mat(shp(1),shp(2)))

    mat(:,:)=A(:,:)

    write(*,*) 'calling dgeev (1)'
    call dgeev('N','V',n, mat ,n,rpart,ipart,vl,ldvl,vr,ldvr, work, lwork, info)
    lwork = min(lwmax, int(work(1)))
!    write(*,*) 'lwork = ', lwork
    write(*,*) 'calling dgeev (2)'

    call dgeev('N','V',n, mat,n,rpart,ipart,vl,ldvl,vr,ldvr, work, lwork, info)

    eigVals(:) = rpart(:) + im*(ipart(:))
    write(*,*) 'completed dgeev (2)'

    if (info > 0) then
       write(*,*) 'Failed to compute eigenvalues'
       stop
    end if

    do i1 = 1, n
       J = 1
       do while (J .LE. n)
          if(ipart(j) .EQ. 0d0) then
             eigVecs_(i1, j) = vr(i1,j)
             j=j+1
          else
             eigVecs_(i1,j  ) = vr(i1, j) + im* vr(i1,j+1)
             eigVecs_(i1,j+1) = vr(i1, j) - im* vr(i1,j+1)
             j=j+2
          end if
       end do
    end do
    
    
    ! if (ignore == .true.) then
    !    write(*,*) 'ignore=True'
    !    do i1 = 1, n
    !       do i2 =1,4
    !          coss(i2,i1) = real(dot_product(sym_(i2,4:), eigVecs_(:,i1)))
    !          coss(i2,i1) = coss(i2,i1)!/sqrt((dot_product(sym_(i2,4:),sym_(i2,4:))))
    !       end do
    !    end do
    !    do i1 = 1,n
    !       if ((coss(1,i1)**2+ coss(2, i1)**2+ coss(3,i1)**2 + coss(4,i1)**2)>=0.99) then
    !          eigVals(i1)=eigVals(i1)-1d0
    !       end if
    !    end do
    ! end if

    if (returnVecs) then
       eigvecs=eigvecs_
    end if

  end subroutine GetEigenValues

end module matrix
! subroutine CheckMatrix(n, ignore, matvec, A)
!     implicit none
!   integer,          intent(in)  :: n, ignore
!   integer, parameter            :: dp = selected_real_kind(15,307)
!   double precision              :: x(n)
!   double precision              :: y1(n), y2(n)
!   external                      :: matvec
!   double precision, intent(in) :: A(n - ignore, n - ignore)
!   integer :: i

!   x(:) = 0.0_dp
!   x(5:6) = 1.0_dp
!   call multj(n, x, y1)
!   y2(4:n) = A(:, 2) + A(:, 3)
!   write(*,*) 'CheckMatrix', (sum(abs(y1(ignore+1:)-y2(ignore+1:)))/(n-ignore))

! end subroutine CheckMatrix
