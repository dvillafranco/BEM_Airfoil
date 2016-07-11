module BEM2d_mod

   use TypeKinds
   use BEM2d_typ
   use Globals
   use IO
   use Greens_mod
   use Vortex_mod
   !use f77_lapack
   implicit none
   private

   ! Module procedures
   public :: getCollPts, WakeTE, EvolveWake
   public :: Solve!, CalcCP
   public :: sourceVector, wakeSource
   public :: ludcmp, lubksb

   private :: Solve0, SolveM

   ! TEcond=1 parameters
   real(kind=dp), private :: d1,d2
   integer      , private :: m1,m2

contains

subroutine getCollPts(W)
   type(Wing)                   , intent(inout) :: W
   integer                                      :: i,N1,N2, idx

   idx = 0
   do i = 1, Npieces
      idx = idx + size(W%B(i,1)%X,1) - 1
      W%idx(i) = idx
   end do

   if(associated(W%Xc)) then
      deallocate(W%Xc)
   end if
   nullify(W%Xc)
   allocate(W%Xc(idx,2))

   N1 = 1
   do i = 1, Npieces
      N2 = W%idx(i)
      W%Xc(N1:N2,:) = 0.5_dp * (W%B(i,1)%X(N1:N2,:) + W%B(i,1)%X(N1+1:N2+1,:))
      N1 = N2 + 1
   end do
end subroutine getCollPts

subroutine WakeTE(B,W)
   type(Surf), intent(inout) :: B, W

   W%X(1,:) = B%X(1,:)
   if(Steady) then 
      W%X(2,1) = W%X(1,1) + D_Start_Vortex
      W%X(2,2) = W%X(1,2)
   else
      if(Fixed_Wake .or. (.not. BisectorWake)) then
         W%X(2,1) = W%X(1,1) + 0.5_dp*Time(1)
         W%X(2,2) = W%X(1,2)
      !if(TEcond==1) then  ! Mangler-Smith parallel to lower edge
      !   W%X(2,:) = B%X(Nx,:) - B%X(Nx-1,:)
      !   W%X(2,:) = W%X(1,:) + 0.5_dp * Time(1) * &
      !                     W%X(2,:) / sqrt(W%X(2,1)**2 + W%X(2,2)**2)
      !else if(TEcond==2) then  ! Mangler-Smith parallel to upper edge
      !   W%X(2,:) = B%X(1,:) - B%X(2,:)
      !   W%X(2,:) = W%X(1,:) + 0.5_dp * Time(1) * &
      !                     W%X(2,:) / sqrt(W%X(2,1)**2 + W%X(2,2)**2)
      else! if(BisectorWake) then
         W%X(2,:) = 0.5_dp * (B%X(2,:) + B%X(size(B%X,1)-1,:))
         W%X(2,:) = W%X(1,:) - W%X(2,:)
         W%X(2,:) = W%X(1,:) + 0.5_dp * Time(1) * &
                           W%X(2,:) / sqrt(W%X(2,1)**2 + W%X(2,2)**2)
      end if
   end if
end subroutine WakeTE

subroutine EvolveWake(W,dt)
   ! NB: Evolution assumes incompressible !
   type(Wing)                , intent(inout)  :: W
   real(kind=dp)             , intent(in)     :: dt

   type(Greens)                               :: Gb,Gw
   real(kind=dp), dimension(:,:), allocatable :: V
   integer :: N, i, j, k

   call nullGreens(Gb)
   call nullGreens(Gw)
   do i = 1, Npieces
      N = size(W%W(i,Timeiter)%X,1)-1
      !N = Timeiter
      allocate(V(N,2))

      nullify(W%W(i,Timeiter+1)%X,W%W(i,Timeiter+1)%phi,W%B(i,Timeiter+1)%phi)
      allocate(W%W(i,Timeiter+1)%X(N+2,2),W%W(i,Timeiter+1)%phi(N+1))
      allocate(W%B(i,Timeiter+1)%phi(size(W%B(i,Timeiter)%phi)))

      ! Freestream velocity * twopi
      V(:,1) = twopi
      V(:,2) = 0.0_dp

      if(.not.Fixed_Wake)then
         do j = 1, Npieces
            call CalcGreens(W%B(j,    1   )%X, W%W(i,Timeiter)%X(2:,:), Gb)
            call CalcGreens(W%W(j,Timeiter)%X, W%W(i,Timeiter)%X(2:,:), Gw)
            do k = 1, 2
               V(:,k) = V(:,k) + matmul(Gb%vd(:,:,k),W%B(i,Timeiter)%phi ) + &
                                 matmul(Gb%vs(:,:,k),W%B(i,Timeiter)%phin) + &
                                 matmul(Gw%vd(:,:,k),W%W(i,Timeiter)%phi )
            end do
            if(DoVort)then
               call CalcVortVel(W%W(i,Timeiter)%X(2:,:), V)
            end if
            call killGreens(Gb)
            call killGreens(Gw)
         end do
      end if

      call WakeTE(W%B(i,1),W%W(i,Timeiter+1))
      W%W(i,Timeiter+1)%X(3:,:) = W%W(i,Timeiter)%X(2:,:) + dt*V/twopi
      W%W(i,Timeiter+1)%phi(1)  = W%W(i,Timeiter)%phi(1)
      W%W(i,Timeiter+1)%phi(2:) = W%W(i,Timeiter)%phi(1:N)

      deallocate(V)
   end do

end subroutine EvolveWake

subroutine Solve(W,recalc)
   type(Wing)                   , intent(inout) :: W
   logical                      , intent(in)    :: recalc

   if(Mach > 0.0_dp)then
      call SolveM(W,recalc)
   else
      call Solve0(W,recalc)
   end if   
end subroutine Solve

subroutine Solve0(W, recalc)
   ! Incompressible routine
   type(Wing)                   , intent(inout) :: W
   logical                      , intent(in)    :: recalc

   type(Greens) , dimension(:)  , allocatable, save :: Gb,Gw

   real(kind=dp), dimension(:,:), allocatable, save :: array, lufac, BC
   real(kind=dp), dimension(:)  , allocatable, save :: rhs
   integer      , dimension(:)  , allocatable, save :: lupiv

   integer       :: i,j,k,M,N,N1,N2,ii,jj
   real(kind=dp) :: tmp, dt

   N = size(W%Xc,1)  ! no coll pts
   M = N + Npieces   ! no equations

   if(.not.allocated(Gb))then
      allocate(Gb(Npieces),Gw(Npieces))
   end if
      
   if(Timeiter==1)then
      dt = Time(1)
   else
      dt = Time(Timeiter) - Time(Timeiter-1)
   end if

   if(Timeiter==1 .or. recalc)then

      if(Timeiter==1)then
         allocate(array(M,M),lufac(M,M),rhs(M),lupiv(M))
         array = 0.0_dp
         lufac = 0.0_dp
         lupiv = 0.0_dp
         rhs   = 0.0_dp
      end if

      N1 = 1
      do j = 1, Npieces
         call nullGreens(Gb(j))
         call nullGreens(Gw(j))
         allocate(BC(size(W%B(j,1)%X,1)-1,size(Time,1)))
         call CalcGreens(W%B(j,1)%X, W%Xc, Gb(j), BC)

         ! Store BC data in W structure
         do k = 1, size(Time,1)
            if(associated(W%B(j,k)%phin))then
               deallocate(W%B(j,k)%phin)
            end if
            nullify(W%B(j,k)%phin)
            allocate(W%B(j,k)%phin(size(W%B(j,1)%X,1)-1))
            W%B(j,k)%phin = BC(:,k)
         end do
         deallocate(BC)

         N2 = W%idx(j)

         !  Morino formulation... Take collocation pts just outside the body
         ! A = 2pi H I - Gd  = 2pi - Gd    where diag(Gd) = pi
         !  Katz formulation... Take collocation pts just inside the body
         ! A = 2pi H I - Gd  = 0 - Gd    where diag(Gd) = -pi

         array(1:N,N1:N2) = -Gb(j)%d
         N1 = N2 + 1
      end do
      do i = 1, N
         array(i,i)  = pi
      end do

   end if

   rhs   = 0.0_dp
   lufac = array

   if(DoVort)then
      call VortInfl(W)
   end if

   N1 = 1
   do i = 1, Npieces
      N2 = W%idx(i)
      call CalcGreens(W%W(i,Timeiter)%X, W%Xc, Gw(i))

      ! assemble forcing vector, b = -gsb*bc + gdw(:,2>)*phiw(2>)
      ! Morino formulation
      ! rhs(1:N) = -matmul(Gb%s,BC(:,Timeiter))
      ! Katz formulation
      rhs(1:N)  = rhs(1:N) + matmul(Gb(i)%s,W%B(i,Timeiter)%phin)  &
                + matmul(Gw(i)%d(:,2:),W%W(i,Timeiter)%phi(2:))
                !+ matmul(Gw(i)%d(:,2:Timeiter),W%W(i,Timeiter)%phi(2:Timeiter))
      rhs(N+1:) = 0.0_dp

      ! assemble implicit matrix, A = (2pi)H(f)I - Gdb - Gdw*(s+ - s-)
      lufac(1:N,N+i) = -Gw(i)%d(:,1)

      select case (TEcond)
         case(1)           ! dphiTE(t) = phi+(t-x+/U) - phi-(t-x+/U)
		if (Timeiter ==1 .and. ot) then 
            d1 = abs( 0.5_dp*(W%W(i,Timeiter)%X(1,1)+W%W(i,Timeiter)%X(2,1)) &
                      - W%Xc(N2,1) ) / oldtime 
            d2 = abs( 0.5_dp*(W%W(i,Timeiter)%X(1,1)+W%W(i,Timeiter)%X(2,1)) &
                      - W%Xc(N1,1) ) /oldtime 
		else
            d1 = abs( 0.5_dp*(W%W(i,Timeiter)%X(1,1)+W%W(i,Timeiter)%X(2,1)) &
                      - W%Xc(N2,1) ) / dt
            d2 = abs( 0.5_dp*(W%W(i,Timeiter)%X(1,1)+W%W(i,Timeiter)%X(2,1)) &
                      - W%Xc(N1,1) ) / dt
		endif
            !print *, "TE wake delays: d1=",d1," , d2=",d2
            ! take average delay
            d1 = 0.5_dp * (d1 + d2)
            d2 = d1

            m1 = int(d1)
            d1 = d1 - real(m1,kind=dp)
            !m2 = int(d2)
            !d2 = d2 - real(m2,kind=dp)

            lufac(N+i,N+i) =  1.0_dp
            if(m1==0)then
               lufac(N+i,N2) = d1 - 1.0_dp
               lufac(N+i,N1) = 1.0_dp - d1
            else if(m1 < Timeiter)then
               rhs(N+i) = rhs(N+i) + (1.0_dp-d1)*W%B(i,Timeiter-m1)%phi(N2-N1+1)
               rhs(N+i) = rhs(N+i) - (1.0_dp-d1)*W%B(i,Timeiter-m1)%phi(1)
            end if
            if(m1 < Timeiter-1)then
               rhs(N+i) = rhs(N+i) + d1*W%B(i,Timeiter-m1-1)%phi(N2-N1+1)
               rhs(N+i) = rhs(N+i) - d1*W%B(i,Timeiter-m2-1)%phi(1)
            end if
            if(Timeiter ==1 .and. ot)then
               rhs(N+i) = rhs(N+i) + d1*W%B(i,Timeiter-m1)%phi(N2-N1+1)
               rhs(N+i) = rhs(N+i) - d1*W%B(i,Timeiter-m2)%phi(1)
            end if
            !if(m2==0)then
            !   lufac(N+i,N1) = 1.0_dp - d2
            !else if(m2 < Timeiter)then
            !   rhs(N+i) = rhs(N+i) - (1.0_dp-d2)*W%B(i,Timeiter-m2)%phi(1)
            !end if
            !if(m2 < Timeiter-1)then
            !   rhs(N+i) = rhs(N+i) - d2*W%B(i,Timeiter-m2-1)%phi(1)
            !end if

         case default      ! default Kutta condition
            lufac(N+i, N1) =  1.0_dp
            lufac(N+i, N2) = -1.0_dp
            lufac(N+i,N+i) =  1.0_dp
            rhs(N+i)       =  0.0_dp
      end select

      N1 = N2 + 1

   end do    ! Npieces

   ! compute lu factorization
   call ludcmp(lufac,M,M,lupiv,tmp)          ! NR routine
   !call la_getrf(M,M,lufac,M,lupiv, ier)    ! LAPACK routine

   ! solve for phi
   call lubksb(lufac,M,M,lupiv,rhs)                    ! NR routine
   !call la_getrs("n",M,1,lufac,M,lupiv,rhs,M, ier)    ! LAPACK routine
   !   if(ier/=0) then
   !      print *, "error: la_getrs error in Solve"
   !      stop
   !   end if

   N1 = 1
   do i = 1, Npieces
      N2 = W%idx(i)
      W%B(i,Timeiter)%phi    = rhs(N1:N2)
      W%W(i,Timeiter)%phi(1) = rhs(N+i)
      !call CalcCP(Xb,Pb,Cp,CL)
      N1 = N2 + 1
   end do

end subroutine Solve0


subroutine SolveM(W, recalc)
   ! Compressible routine
   ! NB: TE condition still incompressible !
   type(Wing)                   , intent(inout) :: W
   logical                      , intent(in)    :: recalc

   type(Greens) , dimension(:)  , allocatable, save :: Gb,Gw

   real(kind=dp), dimension(:,:), allocatable, save :: array, lufac, BC
   real(kind=dp), dimension(:)  , allocatable, save :: rhs
   integer      , dimension(:)  , allocatable, save :: lupiv

   real(kind=dp), dimension(:,:), allocatable       :: alf, src,dbl
   integer      , dimension(:,:), allocatable       :: ith

   integer       :: i,j,k,M,N,N1,N2, Nel
   real(kind=dp) :: tmp, dt

	write(*,*) 'getting to here for mach', Mach
   N = size(W%Xc,1)  ! no coll pts
   M = N + Npieces   ! no equations

   if(.not.allocated(Gb))then
      allocate(Gb(Npieces),Gw(Npieces))
   end if
      
   if(Timeiter==1)then
      dt = Time(1)
   else
      dt = Time(Timeiter) - Time(Timeiter-1)
   end if

   if(Timeiter==1 .or. recalc)then

      if(Timeiter==1)then
         allocate(array(M,M),lufac(M,M),rhs(M),lupiv(M))
         array = 0.0_dp
         lufac = 0.0_dp
         lupiv = 0.0_dp
         rhs   = 0.0_dp
      end if

	write(*,*) 'getting to here for before nullgreens' 
      N1 = 1
      do j = 1, Npieces
         Nel = size(W%B(j,1)%X,1) - 1
         call nullGreens(Gb(j))
         call nullGreens(Gw(j))
         allocate(BC(Nel,size(Time,1)))
         call CalcGreens(W%B(j,1)%X, W%Xc, Gb(j), BC)

         ! Store BC data in W structure
         do k = 1, size(Time,1)
            if(associated(W%B(j,k)%phin))then
               deallocate(W%B(j,k)%phin)
            end if
            nullify(W%B(j,k)%phin)
            allocate(W%B(j,k)%phin(Nel))
            W%B(j,k)%phin = BC(:,k)
         end do
         deallocate(BC)
	write(*,*) 'getting to here for after BC deallo'

         N2 = W%idx(j)

         !  Morino formulation... Take collocation pts just outside the body
         ! A = 2pi H I - Gd  = 2pi - Gd    where diag(Gd) = pi
         !  Katz formulation... Take collocation pts just inside the body
         ! A = 2pi H I - Gd  = 0 - Gd    where diag(Gd) = -pi

         allocate(alf(N,Nel), ith(N,Nel))
         alf = Gb(j)%th/dt
         ith = int(alf)
         alf = alf - real(ith,kind=dp)

	write(*,*) 'getting to here after allo alf'
         do i = 1, Nel
            Gb(j)%d(i,i)  = -pi
         end do
         where(ith==0)
            array(1:N,N1:N2) = -(1.0_dp - alf) * (Gb(j)%d + Gb(j)%r)
         end where

         ! Add Euler forward time differential
         alf = Gb(j)%thp/dt
         ith = int(alf)
         alf = alf - real(ith,kind=dp)
         where(ith==0)
            array(1:N,N1:N2) =  array(1:N,N1:N2) + (1.0_dp - alf) * Gb(j)%r
         end where

         ! Add Figueiredo discretization
         if(Figueiredo)then
            alf = Gb(j)%thpF/dt
            ith = int(alf)
            alf = alf - real(ith,kind=dp)
            where(ith==0)
               array(1:N,N1:N2) = array(1:N,N1:N2) + (1.0_dp - alf) * Gb(j)%rF
            end where
         end if

         deallocate(alf, ith)

	write(*,*) 'getting to here for after deallo alf'
         N1 = N2 + 1
      end do

   end if

   rhs   = 0.0_dp
   lufac = array

	write(*,*) 'getting to here before calcgreens'
   N1 = 1
   do i = 1, Npieces
      N2 = W%idx(i)
      call CalcGreens(W%W(i,Timeiter)%X, W%Xc, Gw(i))

      ! assemble forcing vector, b = -gsb*bc + gdw(:,2>)*phiw(2>)
      ! Morino formulation
      ! rhs(1:N) = -matmul(Gb%s,BC(:,Timeiter))
      ! Katz formulation
      ! rhs(1:N) =  matmul(Gb%s,BC(:,Timeiter))

	write(*,*) 'after calcgreens'
      Nel = size(W%B(i,1)%X,1) - 1
      allocate(src(N,Nel), dbl(N,Nel))

      allocate(alf(N,Nel), ith(N,Nel))

	write(*,*) 'after allocates src etc'
      alf = Gb(i)%thp/dt
      ith = int(alf)
      alf = alf - real(ith,kind=dp)
      call sourceVector(W%B,i,ith,alf, src,dbl, .false.)
      rhs(1:N)  = rhs(1:N) - sum(Gb(i)%r*dbl, dim=2)

	write(*,*) 'after sourcevector1'
      alf = Gb(i)%th/dt
      ith = int(alf)
      alf = alf - real(ith,kind=dp)
      call sourceVector(W%B,i,ith,alf, src,dbl, .false.)
      rhs(1:N)  = rhs(1:N) + sum(Gb(i)%s*src + (Gb(i)%d + Gb(i)%r)*dbl, dim=2)

      if(Figueiredo)then
         rhs(1:N) = rhs(1:N) + sum(Gb(i)%rF*dbl, dim=2)
         
         alf = Gb(i)%thpF/dt
         ith = int(alf)
         alf = alf - real(ith,kind=dp)
         call sourceVector(W%B,i,ith,alf, src,dbl, .false.)
         rhs(1:N) = rhs(1:N) - sum(Gb(i)%rF*dbl, dim=2)
      end if

      deallocate(alf,ith, src,dbl)
	write(*,*) 'after deallocate src etc'

      !Nel = size(W%W(i,1)%X,1) - 1
      Nel = Timeiter
      if(Nel > 1) then
	write(*,*) 'in if allocate'
         allocate(alf(N,Nel), ith(N,Nel), dbl(N,Nel))
	write(*,*) 'allocate workds'
	write(*,*) dt, Gw(i)%thp 
         alf = Gw(i)%thp/dt
         ith = int(alf)
         alf = alf - real(ith,kind=dp)
	write(*,*) 'before wakesource', i, ith, alf, dbl
         call wakeSource(W%W,i,ith,alf, dbl)
	write(*,*) 'after wkae course'
         rhs(1:N)  = rhs(1:N) - sum(Gw(i)%r(:,2:Nel)*dbl(:,2:Nel), dim=2)

         alf = Gw(i)%th/dt
         ith = int(alf)
         alf = alf - real(ith,kind=dp)
         call wakeSource(W%W,i,ith,alf, dbl)
	write(*,*) 'after 2nd wakesource'
         rhs(1:N)  = rhs(1:N) + sum( &
            (Gw(i)%d(:,2:Nel)+Gw(i)%r(:,2:Nel))*dbl(:,2:Nel), dim=2)

         if(Figueiredo)then
	write(*,*) 'in fig if'
            rhs(1:N) = rhs(1:N) + sum(Gb(i)%rF(:,2:Nel)*dbl(:,2:Nel),dim=2)

            alf = Gb(i)%thpF/dt
            ith = int(alf)
            alf = alf - real(ith,kind=dp)
            call wakeSource(W%W,i,ith,alf, dbl)
	write(*,*) 'after wake source call in figuer'
            rhs(1:N) = rhs(1:N) - sum(Gb(i)%rF(:,2:Nel)*dbl(:,2:Nel), dim=2)
         end if

         deallocate(alf,ith, dbl)
	write(*,*) 'wake source part deallo alf'
      end if

      rhs(N+1:) = 0.0_dp

      ! assemble implicit matrix, A = (2pi)H(f)I - Gdb - Gdw*(s+ - s-)
      lufac(1:N,N+i) = -Gw(i)%d(:,1)

	write(*,*) 'getting to here before select case'
      select case (TEcond)
         case(1)           ! dphiTE(t) = phi+(t-x+/U) - phi-(t-x+/U)
            d1 = abs( 0.5_dp*(W%W(i,Timeiter)%X(1,1)+W%W(i,Timeiter)%X(2,1)) &
                      - W%Xc(N2,1) ) / (dt*Beta)
            d2 = abs( 0.5_dp*(W%W(i,Timeiter)%X(1,1)+W%W(i,Timeiter)%X(2,1)) &
                      - W%Xc(N1,1) ) / (dt*Beta)
            !print *, "TE wake delays: d1=",d1," , d2=",d2
            m1 = int(d1)
            d1 = d1 - real(m1,kind=dp)
            m2 = int(d2)
            d2 = d2 - real(m2,kind=dp)

            lufac(N+i,N+i) =  1.0_dp
            if(m1==0)then
               lufac(N+i,N2) = d1 - 1.0_dp
            else if(m1 < Timeiter)then
               rhs(N+i) = rhs(N+i) + (1.0_dp-d1)*W%B(i,Timeiter-m1)%phi(N2-N1+1)
            end if
            if(m1 < Timeiter-1)then
               rhs(N+i) = rhs(N+i) + d1*W%B(i,Timeiter-m1-1)%phi(N2-N1+1)
            end if
            if(m2==0)then
               lufac(N+i,N1) = 1.0_dp - d2
            else if(m2 < Timeiter)then
               rhs(N+i) = rhs(N+i) - (1.0_dp-d2)*W%B(i,Timeiter-m2)%phi(1)
            end if
            if(m2 < Timeiter-1)then
               rhs(N+i) = rhs(N+i) - d2*W%B(i,Timeiter-m2)%phi(1)
            end if

         case default      ! default Kutta condition
            lufac(N+i, N1) =  1.0_dp
            lufac(N+i, N2) = -1.0_dp
            lufac(N+i,N+i) =  1.0_dp
            rhs(N+i)       =  0.0_dp
      end select

      N1 = N2 + 1

   end do    ! Npieces

   ! compute lu factorization
   call ludcmp(lufac,M,M,lupiv,tmp)          ! NR routine
   !call la_getrf(M,M,lufac,M,lupiv, ier)    ! LAPACK routine

   ! solve for phi
   call lubksb(lufac,M,M,lupiv,rhs)                    ! NR routine
   !call la_getrs("n",M,1,lufac,M,lupiv,rhs,M, ier)    ! LAPACK routine
   !   if(ier/=0) then
   !      print *, "error: la_getrs error in Solve"
   !      stop
   !   end if

   N1 = 1
   do i = 1, Npieces
      N2 = W%idx(i)
      W%B(i,Timeiter)%phi    = rhs(N1:N2)
      W%W(i,Timeiter)%phi(1) = rhs(N+i)
      !call CalcCP(Xb,Pb,Cp,CL)
      N1 = N2 + 1
   end do

end subroutine SolveM

subroutine sourceVector(B,ipiece,ith,alf, src,dbl, known)
   ! Assumes zero source strength for t <= 0
   type(Surf)   , dimension(:,:), pointer     :: B
   integer                      , intent(in)  :: ipiece
   integer      , dimension(:,:), intent(in)  :: ith
   real(kind=dp), dimension(:,:), intent(in)  :: alf
   real(kind=dp), dimension(:,:), intent(out) :: src, dbl
   logical                      , intent(in)  :: known
   integer :: i,j, t

   src = 0.0_dp
   dbl = 0.0_dp
   do j = 1, size(ith,2)
   do i = 1, size(ith,1)
      t = Timeiter - ith(i,j)
      if(t > 0)then
         src(i,j) = src(i,j) + (1.0_dp - alf(i,j)) * B(ipiece,t)%phin(j)
         if(ith(i,j) > 0 .or. known) then
            dbl(i,j) = dbl(i,j) + (1.0_dp - alf(i,j)) * B(ipiece,t)%phi(j)
         end if
      end if
      if(t > 1)then
         src(i,j) = src(i,j) + alf(i,j) * B(ipiece,t-1)%phin(j)
         dbl(i,j) = dbl(i,j) + alf(i,j) * B(ipiece,t-1)%phi(j)
      end if
   end do
   end do
end subroutine sourceVector

subroutine wakeSource(W,ipiece,ith,alf, dbl)
   ! Assumes zero source strength for t <= 0 
   ! Assumes fixed-wake model !!!
   type(Surf)   , dimension(:,:), pointer     :: W
   integer                      , intent(in)  :: ipiece
   integer      , dimension(:,:), intent(in)  :: ith
   real(kind=dp), dimension(:,:), intent(in)  :: alf
   real(kind=dp), dimension(:,:), intent(out) :: dbl
   integer :: M, i,j, idx

   M = Timeiter
   dbl = 0.0_dp
   do j = 1, M
   do i = 1, size(ith,1)
      idx = j + ith(i,j)
      if(idx <= M)then
         dbl(i,j) = dbl(i,j) + (1.0_dp - alf(i,j)) * W(ipiece,Timeiter)%phi(idx)
      end if
      if(idx < M)then
         dbl(i,j) = dbl(i,j) + alf(i,j) * W(ipiece,Timeiter)%phi(idx+1)
      end if
   end do
   end do
end subroutine wakeSource

! subroutine CalcCP(Xb,Phi,Cp,CL)
!    real(kind=dp), dimension(:,:), intent(in)  :: Xb,Phi
!    real(kind=dp), dimension(:,:), intent(out) :: Cp
!    real(kind=dp), dimension(:)  , intent(out) :: CL
!    real(kind=dp), dimension(size(Phi,1)+1)    :: phin,lpn
!    real(kind=dp), dimension(size(Phi,1))      :: dphi, lp, vt, t1
!    real(kind=dp), save                        :: c
!    integer :: N, i
!    
!    N  = size(Xb,1)-1
!    lp = sqrt( (Xb(2:N+1,1)-Xb(1:N,1))**2 + (Xb(2:N+1,2)-Xb(1:N,2))**2 )
! 
!    ! Time derivative of potential
!    if(Timeiter==1)then
!       dphi = Phi(:,1) / Time(1)
!       !dphi = 0.0_dp
!       c = maxval(sqrt( (Xb(2:N,1)-Xb(1,1))**2 + (Xb(2:N,2)-Xb(1,2))**2 ))
!       print *, "Chord length = ", c
!    else
!       dphi = ( Phi(:,Timeiter) - Phi(:,Timeiter-1) ) / &
!              ( Time(Timeiter)  - Time(Timeiter-1)  )
!    end if
! 
!    ! determine nodal potential via weighted averaging by panel size
!    lpn         = 0.0_dp
!    lpn(1: N )  = lp
!    lpn(2:N+1)  = lp + lpn(2:N+1)
!    phin        = 0.0_dp
!    phin(1: N ) = Phi(:,Timeiter)*lp
!    phin(2:N+1) = Phi(:,Timeiter)*lp + phin(2:N+1)
!    phin        = phin / lpn
!    
!    ! tangential velocity (dphi/dx1 + U.t = dphi/dx1 + t1)
!    t1 = (Xb(2:N+1,1) - Xb(1:N,1)) / lp
!    vt = (phin(2:N+1) - phin(1:N)) / lp + t1
!    Cp(:,Timeiter) = 1.0_dp - vt**2 - 2.0_dp*dphi
! 
!    ! CL = -int Cp dx n3 / c,  n3 = t1
!    CL(Timeiter) = sum(-Cp(:,Timeiter)*lp * t1) / c
! 
! end subroutine CalcCP


! |----------------------------|
! |  Matrix Solution Routines  |
! |----------------------------|

subroutine ludcmp(a,n,np,indx,d)
! Given a matrix a(1:n,1:n), with physical dimension np by np, this routine
! replaces it by the LU decomposition of a rowwise permutation of itself.  
! a and n are input.  a is output, arranged as in equation (2.3.14 above); 
! indx(1:n) is an output vector that records the row permutation effected
! by the partial pivoting; d is output as +/- 1 depending on whether the 
! number of row interchanges was even or odd, respectively.  This routine
! is used in combination with lubksb to solve linear equation or invert a 
! matrix.
   real(kind=dp), dimension(:,:), intent(inout) :: a
   integer                      , intent(in)    :: n,np
   integer      , dimension(:)  , intent(out)   :: indx
   real(kind=dp)                , intent(out)   :: d

   real(kind=dp), dimension(:), allocatable :: vv  !! implicit row scaling
   real(kind=dp)                            :: aamax,dum,ssum
   integer                                  :: i,imax,j,k, err

   err = 0
   if(size(a,1) /= np .or. size(a,2) /= np .or. size(indx) /= n) then
      print *, "ludcmp error: wrong input sizes"
      stop
   end if
   allocate(vv(n), stat=err)
   if(err /= 0) then
      print *, "ludcmp error: cannot allocate vv"
      stop
   end if

   d = 1.0_dp    ! No row interchanges yet
   do i = 1, n
      aamax = 0.0_dp
      do j = 1, n
         if(abs(a(i,j)) > aamax) then
            aamax = abs(a(i,j))
         end if
      end do
      if(aamax == 0.0_dp) then
         print *, "ludcmp error: singular matrix"
      end if
      vv(i) = 1.0_dp / aamax
   end do
   do j = 1, n
      do i = 1, j-1
         ssum = a(i,j)
         do k = 1, i-1
            ssum = ssum - a(i,k)*a(k,j)
         end do
         a(i,j) = ssum
      end do
      aamax = 0.0_dp
      do i = j, n
         ssum = a(i,j)
         do k = 1, j-1
            ssum = ssum - a(i,k)*a(k,j)
         end do
         a(i,j) = ssum
         dum = vv(i)*abs(ssum)
         if(dum >= aamax) then
            imax  = i
            aamax = dum
         end if
      end do
      if(j /= imax) then
         do k = 1, n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
         end do
         d = -d
         vv(imax) = vv(j)
      end if
      indx(j) = imax
      if(a(j,j) == 0.0_dp) then
      ! if the pivot element is zero, the matrix is singular to the precision
      ! of the algorithm.  For some applications, it is desirable to
      ! substitute TINY for zero.  
         a(j,j) = tiny(1.0_dp)
      end if
      if(j /= n) then
         dum = 1.0_dp/a(j,j)
         do i = j+1, n
            a(i,j) = a(i,j)*dum
         end do
      end if
   end do
   deallocate(vv)
   return

end subroutine ludcmp

subroutine lubksb(a,n,np,indx,b)
! Solves the set of n linear equations AX=B.  Here a is input, not as the
! matrix A, but rather as its LU decomposition, determined by the routine
! ludcmp.  indx is input as the permutation vector returned by ludcmp; b(1:n)
! is input as the right-hand side vector B, and returns with the solution
! vector X.  a, n, np, and indx are not modified by this routine and can be
! left in place for successive calls with difference right-hand sides b.  This
! routine takes into account the possibility that b will begin with many zero
! elements, so it is efficient for use in matrix inversion.
   real(kind=dp), dimension(:,:), intent(in)    :: a
   integer                      , intent(in)    :: n,np
   integer      , dimension(:)  , intent(in)    :: indx
   real(kind=dp), dimension(:)  , intent(inout) :: b

   integer       :: i,ii,j,ll
   real(kind=dp) :: ssum

   if(size(a,1)/=np .or. size(a,2)/=np .or. size(indx)/=n .or. size(b)/=n) then
      print *, "ludcmp error: wrong input sizes"
      stop
   end if

   ii = 0
   do i = 1, n
      ll = indx(i)
      ssum = b(ll)
      b(ll) = b(i)
      if(ii /= 0) then
         do j = ii, i-1
            ssum = ssum - a(i,j)*b(j)
         end do
      else if(ssum /= 0) then
         ii = i
      end if
      b(i) = ssum
   end do
   do i = n, 1, -1
      ssum = b(i)
      do j = i+1, n
         ssum = ssum - a(i,j) * b(j)
      end do
      b(i) = ssum / a(i,i)
   end do
   return

end subroutine lubksb

end module BEM2d_mod


! **********************************
! **** OLD TEcond investigation ****
! **********************************
!
! ! Initialize TEcond parameters for Kutta Cond
! if(Timeiter==1)then
!    d1  = sqrt( (Xb(1,1)-Xc(N,1))**2 + (Xb(1,2)-Xc(N,2))**2 )
!    d2  = sqrt( (Xb(1,1)-Xc(1,1))**2 + (Xb(1,2)-Xc(1,2))**2 )
!    t11 = (Xb(1,1) - Xc(N,1)) / d2
!    t21 = (Xb(1,1) - Xc(1,1)) / d1
!    a   = 2.0_dp / dt
! 
!    p1p = 0.0_dp
!    p2p = 0.0_dp
!    dpp = 0.0_dp
! else
!    p1p = Pb(N,Timeiter-1)
!    p2p = Pb(1,Timeiter-1)
!    dpp = Pw(1)
! end if
! select case (TEcond)
!    case(1)   ! Mangler-Smith parallel to lower edge
!       if(Timeiter==1)then
!          db = sqrt( (Xc(1,1)-Xc(2,1))**2 + (Xc(1,2)-Xc(2,2))**2 )
!          dw = sqrt( ( 0.5_dp*(Xw(1,1)+Xw(2,1)) - Xc(1,1))**2 + &
!                     ( 0.5_dp*(Xw(1,2)+Xw(2,2)) - Xc(1,2))**2 )
!       end if
!       lufac(N+1, 1 ) =  1.0_dp * (1.0_dp - dw/db)
!       lufac(N+1, 2 ) =  1.0_dp *           dw/db
!       lufac(N+1, N ) = -1.0_dp
!       lufac(N+1,N+1) =  1.0_dp
!       rhs(N+1)       =  0.0_dp
! 
!    case(2)   ! Mangler-Smith parallel to upper edge
!       if(Timeiter==1)then
!          db = sqrt( (Xc(N,1)-Xc(N-1,1))**2 + (Xc(N,2)-Xc(N-1,2))**2 )
!          dw = sqrt( ( 0.5_dp*(Xw(1,1)+Xw(2,1)) - Xc(N,1))**2 + &
!                     ( 0.5_dp*(Xw(1,2)+Xw(2,2)) - Xc(N,2))**2 )
!       end if
!       lufac(N+1, 1 ) =  1.0_dp
!       lufac(N+1,N-1) =  1.0_dp *           dw/db
!       lufac(N+1, N ) = -1.0_dp * (1.0_dp + dw/db)
!       lufac(N+1,N+1) =  1.0_dp
!       rhs(N+1)       =  0.0_dp
! 
!    case(3)   ! Bose-style upper TE stagnation
!       if(Timeiter==1)then
!          pwp = p1p - dpp
!          b   = ( pwp - p2p ) / d2**2 + t21/d2
!          c   = t11**2 - t21**2 + 2.0_dp*(pwp-p1p)/dt
!          d   = b + a
!       end if
!       lufac(N+1, 1 ) =  1.0_dp *           b/d
!       lufac(N+1, N ) = -1.0_dp * (1.0_dp - a/d)
!       lufac(N+1,N+1) =  1.0_dp
!       rhs(N+1)       = -c/d
! 
!    case(4)   ! Bose-style lower TE stagnation
!       if(Timeiter==1)then
!          pwp = dpp + p2p
!          b   = ( pwp - p1p ) / d1**2 + t11/d1
!          c   = t21**2 - t11**2 + 2.0_dp*(pwp-p2p)/dt
!          d   = b + a
!       end if
!       lufac(N+1, 1 ) =  1.0_dp * (1.0_dp - a/d)
!       lufac(N+1, N ) = -1.0_dp *           b/d
!       lufac(N+1,N+1) =  1.0_dp
!       rhs(N+1)       =  c/d
! 
!    case(5)   ! Specify vn = 0 at TE explicitly
!       !Xtmp(:,1) = Xb(1,1)
!       !Xtmp(1,2) = Xb(1,2) + 1.0e-3
!       !Xtmp(1,2) = Xb(1,2) - 1.0e-3
!       !!! DID I use Xtmp in GbTE and Xb(1,:) in GwTE !!!
!       call CalcGreens(Xb,Xb((/1/),:), GbTE)
!       call CalcGreens(Xw,Xb((/1/),:), GwTE)
!       dwn = sqrt( (Xw(1,1)-Xw(2,1))**2 + (Xw(1,2)-Xw(2,2))**2 )
!       n1  = -(Xw(1,2)-Xw(2,2)) / dwn
!       n2  =  (Xw(1,1)-Xw(2,1)) / dwn
!       ! Fix GbTE and GwTE
!       GbTE%vs(1,1,:) = 0.0_dp
!       GbTE%vs(1,N,:) = 0.0_dp
!       GwTE%vs(1,1,:) = 0.0_dp
! 
!       lufac(N+1,1:N) = -GbTE%vd(1,:,1)*n1 -GbTE%vd(1,:,2)*n2
!       lufac(N+1,N+1) = -GwTE%vd(1,1,1)*n1 -GwTE%vd(1,1,2)*n2
!       rhs(N+1) = +sum((GbTE%vs(1,:,1)*n1+GbTE%vs(1,:,2)*n2)*BC(:,Timeiter)) &
!                  +sum((GbTE%vd(1,2:Nw,1)*n1+GwTE%vd(1,2:Nw,2)*n2)*Pw(2:Nw))
!       !lufac(N+1,1:N) = -GbTE%vd(1,:,2)
!       !lufac(N+1,N+1) = -GwTE%vd(1,1,2)
!       !rhs(N+1) = -sum(GbTE%vs(1,:,2)*BC(:,Timeiter)) &
!       !           +sum(GbTE%vd(1,2:Nw,2)*Pw(2:Nw))
! 
!    case(6)           ! dphiTE(t) = phi+(t-x+/U) - phi-(t-x+/U)
!       if(Timeiter==1)then
!          d1 = abs( 0.5_dp*(Xw(1,1)+Xw(2,1))-Xc(N,1) ) / dt
!          d2 = abs( 0.5_dp*(Xw(1,1)+Xw(2,1))-Xc(1,1) ) / dt
!          print *, "TE wake delays: d1=",d1," , d2=",d2
!          m1 = int(d1)
!          d1 = d1 - real(m1,kind=dp)
!          m2 = int(d2)
!          d2 = d2 - real(m2,kind=dp)
!       end if
!       lufac(N+1,N+1) =  1.0_dp
!       if(m1==0)then
!          lufac(N+1,N) = d1 - 1.0_dp
!       else if(m1 < Timeiter)then
!          rhs(N+1) = rhs(N+1) + (1.0_dp - d1)*Pb(N,Timeiter-m1)
!       end if
!       if(m1 < Timeiter-1)then
!          rhs(N+1) = rhs(N+1) + d1*Pb(N,Timeiter-m1-1)
!       end if
!       if(m2==0)then
!          lufac(N+1,1) = 1.0_dp - d2
!       else if(m2 < Timeiter)then
!          rhs(N+1) = rhs(N+1) - (1.0_dp - d2)*Pb(1,Timeiter-m2)
!       end if
!       if(m2 < Timeiter-1)then
!          rhs(N+1) = rhs(N+1) - d2*Pb(1,Timeiter-m2-1)
!       end if
! 
!    case(7)           ! phiTE(t) determined by quadratic fit
!       if(Timeiter==1)then
!          s1 = 0.5_dp*sqrt((Xb(1,1)-Xc(N:N-2:-1,1))**2+(Xb(1,2)-Xc(N:N-2:-1,2))**2)
!          s2 = 0.5_dp*sqrt((Xb(1,1)-Xc(1:3     ,1))**2+(Xb(1,2)-Xc(1:3     ,2))**2)
!          s1(3) = s1(3) + 2.0_dp*s1(2) + 2.0_dp*s1(1)
!          s1(2) = s1(2) + 2.0_dp*s1(1)
!          s2(3) = s2(3) + 2.0_dp*s2(2) + 2.0_dp*s2(1)
!          s2(2) = s2(2) + 2.0_dp*s2(1)
!          a1(1) = s1(3)*s1(2)*(s1(3)-s1(2))
!          a1(2) = s1(1)*s1(3)*(s1(1)-s1(3))
!          a1(3) = s1(2)*s1(1)*(s1(2)-s1(1))
!          a2(1) = s2(3)*s2(2)*(s2(3)-s2(2))
!          a2(2) = s2(1)*s2(3)*(s2(1)-s2(3))
!          a2(3) = s2(2)*s2(1)*(s2(2)-s2(1))
!          a1 = a1/sum(a1)
!          a2 = a2/sum(a2)
!       end if
!       lufac(N+1,1:3)      =  a2
!       lufac(N+1,N:N-2:-1) = -a1
!       lufac(N+1,N+1)      =  1.0_dp
!       rhs(N+1)            =  0.0_dp
! 
!    case default      ! default Kutta condition
!       lufac(N+1, 1 ) =  1.0_dp
!       lufac(N+1, N ) = -1.0_dp
!       lufac(N+1,N+1) =  1.0_dp
!       rhs(N+1)       =  0.0_dp
! end select
! if(MorinoCorr)then
!    dcorr = sqrt( (Xc(N,1)-Xc(1,1))**2 + (Xc(N,2)-Xc(1,2))**2 )
!    rhs(N+1) = rhs(N+1) + dcorr
! end if
