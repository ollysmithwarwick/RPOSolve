!----------------------------------------------------------------------

! Find root f(x) = 0 by Newton method with gmres-hookstep:
!    invert  df(x_n)/dx . s = f(x_n)  
! for s subject to constraint |s| < del  (size of 'trust region'), 
!    then     x_{n+1} = x_n - s
! To define subroutines f,df below, write as
!    df . s = f
! where f is a function to be minimised and operator df provides the
! action of the Jacobian, but may also contain any constraints on the 
! update s (e.g. no shifts in a homogenous direction), arranged such 
! that the rhs of the constraint (in f) is 0 when converged.
!----------------------------------------------------------------------
!
! put initial x(1:n) in newtondata module, updated each iteration
!
! f       evaluates y=f(x):  call f(n,x, y)
! df      evaluates y=Jacobian.x, see GMRESm.f90:  call df(n,x, y)
! dfp     preconditioner for df:  call df(n, x)
! sub     parameterless subroutine called at end of each iteration,
!         so that data could be saved etc.
! dp      dot product:  d = dp(n,a,b)
! m	  gmres dimension (also max num gmres its)
! n 	  dimension of x
! gtol    tolerence for gmres, typically 0.001
! tol     request |f(x)|<tol
! del     initial size of trust region for usefulness of df
!         set = 0d0 for no hookstep, <0d0 for del = |f(x0)|/1d2
! mndl    min size of trust region
! mxdl    max size of trust region, if del<0 set no upper limit
! nits	  max num Newton its 
! info	  on input if =1 print* details
!         on exit if =0 sucessful
!                    =1 gmres failed
!                    =2 reached max iterations, nits
!                    =3 trust region got too small
!                    =4 unable to allocate memory
!							A.P.Willis 2008
!----------------------------------------------------------------------
 module newton
   implicit none
   save
   double precision, allocatable :: new_x(:), new_fx(:)
   double precision              :: new_tol,  new_del
   double precision, allocatable :: J(:, :)
   integer                       :: new_nits, new_gits
!   integer, parameter :: RL = 1, IM = 2

 contains
!    subroutine makepos(n_)
!      implicit none
!      integer, intent(in) :: n_
!      integer :: i1
!      do i1 = 1, n_                                 !     To help visualise - each row follows on from the previous. The tilde quantities take twice as much space as they have real and imag comps.
!         pos(i1,    1, RL) =             i1   + 3   !           S        |    shift      |     phase    | NBAR(1)| NBAR(2) | NBAR(3) | NBAR(4) | ...
!         pos(i1,    2, RL)     =   n_ + i1   + 3   !         NBAR(N_x-2)| NBAR(N_x - 1) |   NBAR(N_x)  | E(1)   |   E(2)  |  E(3)   |  E(4)   | ...
!         pos(i1,        3, RL) = 2*n_ + i1*2 + 2   !            E(N_x-2)|   E(N_x-1)    |     E(N_x)   |   PHITIL(1)      |     PHITIL(2)     | ...
!         pos(i1,        3, IM) = 2*n_ + i1*2 + 3   !    PHITIL(N_x/2-1) |           PHITIL(N_x/2)      | PHITIL(N_x/2+1)  |  PHITIL(N_x/2+2)  | ...
!         pos(i1,        4, RL) = 4*n_ + i1*2 + 2   !    PHITIL(N_x - 1) |            PHITIL(N_x)       |
!         pos(i1,        4, IM) = 4*n_ + i1*2 + 3   !

! !      inv_pos(          i1 + 2, :) = (i1,     NBAR, RL) ! This ended up not being used but might be useful somehow
!    end do
   

!   end subroutine makepos
   
 end module newton


 subroutine newtonhook(f, df, dfp, sub, dotprod, m, n, gtol, &
                       tol, del, mndl, mxdl, nits, info)
   use newton

   implicit none
   external                        :: f, df, dfp, sub
   interface
      double precision function metricType(n_, a, b)
        implicit none
        integer, intent(in) :: n_
        double precision, dimension(:), intent(in) :: a
        double precision, dimension(:), intent(in) :: b
      end function metricType
   end interface
   procedure(metricType)     :: dotprod
   integer,          intent(in)    :: m, n
   double precision, intent(in)    :: gtol, tol, del, mndl, mxdl
   integer,          intent(inout)    :: nits
   integer,          intent(inout) :: info
   integer :: i2
   double precision, allocatable :: v(:,:)
   double precision :: tol_,x_(n),fx_(n), tol__,x__(n),fx__(n),del__
   double precision :: mxdl_, ared, pred, snrm, s(n)
   double precision :: gres, gdel, h((m+1)*m)
   double precision, dimension(:, :), allocatable :: hfull
   integer :: ginfo
!   logical :: JacOut
!   character(16) :: AFmt
!   character(16) :: hessFmt
!   character(50) :: AFile = 'POPI/Output/matrix.dat'
!   character(50) :: hessFile = 'POPI/Output/hessenberg.dat'   
   
   write(*,*) 'NewtonHook: Called. Info =', info
   write(*,*) 'NewtonHook: Called. tol =', tol
 
   allocate(v(n,m+1), stat=ginfo)
   allocate(hfull(n, n), stat = ginfo)

 !  write(hessFmt, "(A1, I0, A7)") "(",n,"E30.18)"
!   write(AFmt, "(A1, I0, A7)") "(",n-3,"E30.18)"
   if(ginfo/=0) then 
      if(info==1) print*, 'newton: unable to allocate memory'
      info = 4
      nits = 0
      return
   end if

   write(*,*) 'NewtonHook: Allocated Memory'

   new_nits = 0
   new_gits = 0
   new_del  = del
   mxdl_    = mxdl
   ginfo    = info
   call f(n,new_x, new_fx)

   new_tol  = dsqrt(dotprod(n,new_fx,new_fx))
   if(del<0d0)  new_del = new_tol / 10d0
   if(del<0d0)  mxdl_   = 1d99
   if(info==1)  print*,'newton: nits=',new_nits,' res=',real(new_tol), 'tol=',real(tol)


   call sub()
   x_   = new_x
   fx_  = new_fx
   tol_ = new_tol
   tol__ = 1d99

   if(new_tol<tol) then
      if(info==1) print*, 'newton: input already converged'
      info = 0
      nits = new_nits
      deallocate(v)
      return
   end if 
   write(*,*) 'NewtonHook: Starting Main Loop'
     				! - - - - Start main loop - - - - -  -
   do while(.true.)

      if(new_del<mndl) then
         if(info==1) print*, 'newton: trust region too small. new_del = ', new_del, 'mndl = ',mndl
         info = 3
         nits = new_nits
         deallocate(v)
         return      
      end if
         			! find hookstep s and update x
      s        = 0d0
      gres     = gtol * new_tol
      gdel     = new_del
      if(ginfo/=2) new_gits = m
      print*, 'max_gits = ', new_gits, 'gdel= ', gdel, 'gres= ',gres
      call gmresm(m,n,s,fx_,df,dfp,dotprod,h,v,gres,gdel,new_gits,ginfo)

      ginfo = info

      new_x = x_ - s
         			! calc new norm, compare with prediction
      call f(n,new_x, new_fx)
      new_tol = dsqrt(dotprod(n,new_fx,new_fx))
      snrm = dsqrt(dotprod(n,s,s))
      ared = tol_ - new_tol
      pred = tol_ - gdel
         
      if(info==1) then 
         print*,'newton: nits=',new_nits,' res=',real(new_tol), ' tol=',real(tol)
         print*,'newton: gits=',new_gits,' del=',real(new_del)
         print*,'newton: |s|=',real(snrm),' pred=',real(pred)
         print*, 'newton: ared=', real(ared)
         print*, 'newton: pred=', real(pred)
         print*,'newton: ared/pred=',real(ared/pred)
      end if

      if(new_tol>tol__) then
         if(info==1) print*, 'newton: accepting previous step'
         new_x   = x__
         new_fx  = fx__
         new_tol = tol__
         new_del = del__
      else if(ared<0d0) then
         if(info==1) print*, 'newton: norm increased, try smaller step'
      !   new_del = snrm * 0.5d0 ORIGINAL
         new_del = snrm * 0.9d0
         ginfo   = 2
      !   x__     = new_x   !TEST
      !   fx__    = new_fx  !TEST
      !   tol__   = new_tol !TEST
      !   del__   = snrm    !TEST
      else if(ared/pred<0.75d0) then
         if(info==1)  print*, 'newton: step ok, trying smaller step'
         x__     = new_x
         fx__    = new_fx
         tol__   = new_tol
         if(ared/pred> 0.1d0)  del__ = snrm
         if(ared/pred<=0.1d0)  del__ = snrm*0.5d0
         new_del = snrm * 0.7d0
         ginfo   = 2
      else if(snrm<new_del*0.9d0) then
         if(info==1) print*, 'newton: step good, took full newton step'
         new_del = min(mxdl_,snrm*2d0)
      else if(new_del<mxdl_*0.9d0) then
         if(info==1) print*, 'newton: step good, trying larger step'
         x__     = new_x
         fx__    = new_fx
         tol__   = new_tol
         del__   = new_del
         new_del = min(mxdl_,snrm*2d0)
         ginfo   = 2
      end if
         				! check if need to try another s
      if(ginfo==2) cycle
         				! end of iteration
      new_nits = new_nits + 1
      call sub()
      x_   = new_x
      fx_  = new_fx
      tol_ = new_tol
      tol__ = 1d99

      if(new_tol<tol) then
         if(info==1) print*, 'newton: converged'
         info = 0
         nits = new_nits
         new_gits = n-1 ! Want a full Krylov subspace generated below
         gres = 0d0
         gdel = 0d0
!         open(12, file = AFile)
         ! Print matrix in reduced (without 3 constraints) form and Upper Hessenberg form
!         if (JacOut) then
!            write(*,*) 'Getting Jacobian'
!            call GetMatrix(n, 3, df, J)
!            write(*,*) 'Jacobian complete'
!         end if
         deallocate(v)
         return
      else if(new_nits==nits) then
         if(info==1) print*, 'newton: reached max its'
         info = 2
         nits = new_nits
         new_gits = n-1 ! Want a full Krylov subspace generated below
         gres = 0d0
         gdel = 0d0
  !       open(12, file = AFile)
         ! Print matrix in reduced (without 3 constraints) form and Upper Hessenberg form
!         if (JacOut) then
!            write(*,*) 'Getting Jacobian'
!            call GetMatrix(n, 3, df, J)
!            write(*,*) 'Jacobian complete'
!         end if
 !        call CheckMatrix(n, 3, df, A)
  !       do i2 = 1, n-3
  !          write(12, fmt=AFmt) A(i2,:)
  !       end do
  !       close(12)
  !       open(12, file = hessFile)
  !       do i2 = 1, n !Only n rows needed for full square hessenberg matrix
  !          write(12, fmt=hessFmt) hfull(i2,:)
  !       end do
  !       close(12)
         deallocate(v)
         return
      end if
   
   end do

 end subroutine newtonhook

