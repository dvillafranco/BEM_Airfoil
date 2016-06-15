module Vortex_mod
   use TypeKinds
   use BEM2d_typ
   use Globals
   use IO
   use Greens_mod

   implicit none
   private

   ! Module procedures
   public  :: InitVortex, EvolVort, VortInfl, CalcVortVel, OutputVort, Velfield
   !public  :: OutputPhi

   ! Module data
   real(kind=dp), private, dimension(:,:), pointer :: XVort  ! location matrix
   real(kind=dp), private                          :: StVort ! strength of vortex

contains

subroutine InitVortex()
   ! Reads vortex position from file "Vort.in"

   allocate(XVort(2,size(Time)))
   open(unit=2,file="Vort."//trim(filetag)//"in",action="read",status="old")
      read(unit=2,fmt=*) XVort(1,1), XVort(2,1), StVort
   close(unit=2)
   !xvort(:,1) = (/ -2.0_dp, 0.0_dp /)
   !stvort     = 0.5_dp
end subroutine Initvortex


subroutine EvolVort(W,dt)
   ! Evolves position of vortex Xvort at current timestep
   ! NB: Evolution assumes incompressible !
   type(Wing)                , intent(inout)  :: W
   real(kind=dp)             , intent(in)     :: dt


   type(Greens)                               :: Gb,Gw
   real(kind=dp), dimension(1,2)              :: V,Xo
   real(kind=dp)	            :: dtt, xvx, xvy 
   integer :: N, i, j, k, nummer, ii

	if (DoVort) then 
   Xo(1,:) = Xvort(:,Timeiter)

   call nullGreens(Gb)
   call nullGreens(Gw)
	
   do i = 1, Npieces
      !N = size(W%W(i,Timeiter)%X,1)-1

      ! Freestream velocity * twopi
      V(:,1) = twopi
      V(:,2) = 0.0_dp

      do j = 1, Npieces
         call CalcGreens(W%B(j,    1   )%X, Xo, Gb)
         call CalcGreens(W%W(j,Timeiter)%X, Xo, Gw)
         do k = 1, 2
            V(:,k) = V(:,k) + matmul(Gb%vd(:,:,k),W%B(i,Timeiter)%phi ) + &
                              matmul(Gb%vs(:,:,k),W%B(i,Timeiter)%phin) + &
                              matmul(Gw%vd(:,:,k),W%W(i,Timeiter)%phi )
         end do
         call killGreens(Gb)
         call killGreens(Gw)
      end do

      Xvort(:,Timeiter+1) = Xvort(:,Timeiter) + dt*V(1,:)/twopi
   end do


	else

	nummer = 3200
  	xvx = -4.5
	xvy = -.165
	dtt = 10.0/nummer
	write(19,*) xvx,xvy,dtt

	do ii = 1,nummer
	Xo(1,1) = xvx
	Xo(1,2) = xvy
   call nullGreens(Gb)
   call nullGreens(Gw)
	
   do i = 1, Npieces
      ! Freestream velocity * twopi
      V(:,1) = twopi
      V(:,2) = 0.0_dp

      do j = 1, Npieces
         call CalcGreens(W%B(j,    1   )%X, Xo, Gb)
         call CalcGreens(W%W(j,size(Time))%X, Xo, Gw)
         do k = 1, 2
            V(:,k) = V(:,k) + matmul(Gb%vd(:,:,k),W%B(i,size(Time))%phi ) + &
                              matmul(Gb%vs(:,:,k),W%B(i,size(Time))%phin) + &
                              matmul(Gw%vd(:,:,k),W%W(i,size(Time))%phi )
         end do
         call killGreens(Gb)
         call killGreens(Gw)
      end do
	xvx = xvx+ dtt*V(1,1)/twopi
	xvy = xvy + dtt*V(1,2)/twopi
	write(19,*) xvx, xvy, dtt
   end do
   end do
	endif

end subroutine EvolVort


subroutine VortInfl(W)
   ! Modifies source strength phin on the airfoil surface to account for
   ! presence of vortex at Xvort
   type(Wing)                   , intent(inout) :: W

   integer                     :: i,j,m,n
   real(kind=dp)               :: mr2, n1,n2,mn, dt
   real(kind=dp), dimension(2) :: r,Xc

   do i = 1, Npieces

      m = size(W%B(i,1)%X,1) - 1      ! no of panels

      !
      ! get the collocation points
      !
      do j = 1, m    ! Source Loop
         Xc(:) = (W%B(i,1)%X(j+1,:) + W%B(i,1)%X(j,:)) * 0.5_dp

         mn   = sqrt((W%B(i,1)%X(j+1,1) - W%B(i,1)%X(j,1))**2 + &
                     (W%B(i,1)%X(j+1,2) - W%B(i,1)%X(j,2))**2)
         n1   = -(W%B(i,1)%X(j+1,2) - W%B(i,1)%X(j,2)) / mn
         n2   =  (W%B(i,1)%X(j+1,1) - W%B(i,1)%X(j,1)) / mn

         r    = Xc - Xvort(:,Timeiter)
         r(1) = r(1) / Beta        ! PGtrans
         mr2  = r(1)**2 + r(2)**2

         W%B(i,Timeiter)%phin(j) = W%B(i,Timeiter)%phin(j) + &
            stvort/twopi/mr2 * (r(2)*n1 - r(1)*n2)

      end do !panels
   end do !Npieces

end subroutine VortInfl 

subroutine CalcVortVel(Xo, V)
   ! Modifies velocity V at observers Xo due to presence of vortex at XVort
   real(kind=dp), dimension(:,:), intent(in)    :: Xo
   real(kind=dp), dimension(:,:), intent(inout) :: V

   integer                     :: i
   real(kind=dp), dimension(2) :: r
   real(kind=dp)               :: mr2

   if(size(Xo,1)/=size(V,1))then
      print *, "CalcVortVel Error: wrong size for first dimension"
      stop
   end if

   do i = 1, size(V,1)
      r = Xo(i,:) - Xvort(:,Timeiter)
      r(1) = r(1) / Beta     !PGtrans
      mr2  = r(1)**2 + r(2)**2 

!print *, 'before', v(i,1), v(i,2)
      ! N.B.: twopi taken out in EvolveWake
      V(i,1) = V(i,1) - StVort / mr2 * r(2)

      V(i,2) = V(i,2) + StVort / mr2 * r(1)

!print *, 'vort infl', v(i,1), v(i,2)
!print *, ' on point ', Xo(i,:)
!print *,  '  ' 

   end do

end subroutine CalcVortVel

subroutine OutputVort()
   ! Writes vortex position in time in file "Vort.m"

   integer :: fid
   call openmatfile("Vort_"//trim(Filetag),"w",fid)
   call savematfile(fid,"Xvort",transpose(Xvort))
   close(unit=fid)
end subroutine OutputVort
 
!subroutine OutputPhi(W,Nx)
!   type(Wing)                , intent(in)  :: W
!    integer :: Nx,i
!   ! Writes phi at ending time
!
!   integer :: fid
!   call openmatfile("Phiend_"//trim(Filetag),"w",fid)
!      do i = 1,Nx
!	      write(fid,*) W%B(1,Timeiter)%phi(i)
!      enddo
!      write(fid,*) W%W(1,Timeiter)%phi(1)
!   close(unit=fid)
!end subroutine OutputPhi

subroutine Velfield(W,dt)
   ! gets velocity at field points for contour plot
   type(Wing)                , intent(inout)  :: W
   real(kind=dp)             , intent(in)     :: dt

   type(Greens)                                 :: Gb,Gw
   real(kind=dp), dimension(:,:)  , allocatable :: Xo, V
   real(kind=dp), dimension(:,:,:), allocatable :: Vel
   real(kind=dp), dimension(:,:,:), allocatable :: Xpoint
   real(kind=dp)   :: dx,dy
   integer :: N, i, j, k, mm,nn,mi,mj, fid

   mm = 50
   nn = 50
   dx = 3.00_dp / real(mm-1,kind=dp)
   dy = 0.55_dp / real(nn-1,kind=dp)
   allocate(Xpoint(mm,nn,2),Vel(mm,nn,2), Xo(mm*nn,2),V(mm*nn,2) )

   do mj = 1, nn
   do mi = 1, mm
      Xpoint(mi,mj,1) = -1.50_dp + dx*(mi-1)
      Xpoint(mi,mj,2) = -0.15_dp + dy*(mj-1)
   end do
   end do
   Xo  = reshape(Xpoint,(/mm*nn,2/))
   V   = 0.0_dp
   Vel = 0.0_dp
   
   N = Timeiter
   do i  = 1, Npieces
      call nullGreens(Gb)
      call nullGreens(Gw)

      ! Freestream velocity * twopi
      V(:,1) = twopi
      V(:,2) = 0.0_dp

      do j = 1, Npieces
         call CalcGreens(W%B(j,    1   )%X, Xo, Gb)
         call CalcGreens(W%W(j,Timeiter)%X, Xo, Gw)
         do k = 1, 2
            V(:,k) = V(:,k) +  &
                     matmul(Gb%vd(:,:,k),W%B(i,Timeiter)%phi ) + &
                     matmul(Gb%vs(:,:,k),W%B(i,Timeiter)%phin) + &
                     matmul(Gw%vd(:,:,k),W%W(i,Timeiter)%phi )
         end do
         call killGreens(Gb)
         call killGreens(Gw)
      end do   ! Npieces - j

   end do   !Npieces - i

   Vel = reshape(V / twopi, (/ mm,nn,2 /) )

   call openmatfile("velfield.","w",fid)
   call savematfile(fid,"X",Xpoint)
   call savematfile(fid,"V",Vel)
   close(unit=fid)

   deallocate(Xpoint,Vel)

end subroutine Velfield



end module Vortex_mod
