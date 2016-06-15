module Greens_mod

   use TypeKinds
   use Globals
   implicit none
   private

   ! Module procedures
   public :: nullGreens,killGreens,CalcGreens

   type, public :: Greens
   ! idx 1 => observer; idx 2 => source; idx 3 => spatial component
      real(kind=dp), dimension(:,:)  , pointer :: s,d,r,rF,th,thp,thpF
      real(kind=dp), dimension(:,:,:), pointer :: vs,vd,vr,vrF
   end type Greens

contains

subroutine nullGreens(G)
   type(Greens), intent(inout) :: G
   nullify(G%s,G%d,G%r,G%rF,G%vs,G%vd,G%vr,G%vrF, G%th,G%thp,G%thpF)
end subroutine nullGreens

subroutine killGreens(G)
   type(Greens), intent(inout) :: G
   if(associated(G%s))then
      deallocate(G%s,G%d,G%r,G%rF,G%vs,G%vd,G%vr,G%vrF, G%th,G%thp,G%thpF)
   end if
   nullify(G%s,G%d,G%r,G%rF,G%vs,G%vd,G%vr,G%vrF, G%th,G%thp,G%thpF)
end subroutine killGreens

subroutine CalcGreens(Xs,Xo, G, BC)
!  ASSUMES: Xs encloses the surface in CLOCKWISE fashion such that the normal
!           vector points into the fluid.  The wake sources must also begin at 
!           the trailing edge and progress downtream (in x1 direction).
   real(kind=dp), dimension(:,:), intent(in)    :: Xs,Xo
   type(Greens)                 , intent(inout) :: G
   real(kind=dp), dimension(:,:), intent(out), optional :: BC

   integer                     :: i,j,m,n
   real(kind=dp)               :: x,z,x2,r1,r12,r2,r22,mrc,mrcp,th1,th2, xc
   real(kind=dp)               :: n1,n2,mn, Mn1,camp,dt
   real(kind=dp), dimension(2) :: ex,ez,r,rc,rcp,gvs!,gvd

   m = size(Xs,1) - 1
   n = size(Xo,1)
   if(size(Xs,2)/=2 .or. size(Xo,2)/=2)then
      print *, "CalcGreens Error: wrong size for second dimension"
      stop
   end if
   if(Timeiter > 1)then
      dt = Time(Timeiter) - Time(Timeiter-1)
   else
      dt = Time(1)
   end if
   
   if(associated(G%s))then
      deallocate(G%s,G%d,G%r,G%rF,G%vs,G%vd,G%vr,G%vrF, G%th,G%thp,G%thpF)
   end if
   nullify(G%s,G%d,G%r,G%rF,G%vs,G%vd,G%vr,G%vrF, G%th,G%thp,G%thpF)
   allocate(G%s(n,m)   ,G%d(n,m)   ,G%r(n,m)   , G%th(n,m),  & 
            G%vs(n,m,2),G%vd(n,m,2),G%vr(n,m,2), G%thp(n,m), &
            G%rF(n,m)  ,G%thpF(n,m),G%vrF(n,m,2) )

   do j = 1, m    ! Source Loop
      ex = Xs(j+1,:) - Xs(j,:)
      ex(1) = ex(1) / Beta      ! PGtrans
      x2 = sqrt(ex(1)**2 + ex(2)**2)
      ex = ex / x2
      ! in P-G coordinates
      ez = (/ -ex(2), ex(1) /)

      ! in real coordinates
      mn   = sqrt((Xs(j+1,1) - Xs(j,1))**2 + (Xs(j+1,2) - Xs(j,2))**2)
      n1   = -(Xs(j+1,2) - Xs(j,2)) / mn
      n2   =  (Xs(j+1,1) - Xs(j,1)) / mn
      Mn1  = Mach * n1
      camp = sqrt(1.0_dp - Mn1*Mn1)

      xc = 0.5_dp * (Xs(j,1)+Xs(j+1,1))
      if(present(BC))then
         select case (BCmode)
            case (1) ! Steady freestream in \hat i direction
               BC(j,:) = -n1
            case (2) ! Sinusoidal vertical (\hat j) gust travelling in \hat i
               BC(j,:) = -n2 * cos(BCfreq*(Time - xc))
            case (3) ! Sinusoidal vertical (\hat j) gust travelling in \hat i
               BC(j,:) = -n2 * sin(BCfreq*(Time - xc))
            case (4) ! Sinusoidal pitch oscillation about x1=0
               BC(j,:) = -n2 * cos(BCfreq*Time) * xc
            case (5) ! Sinusoidal pitch oscillation about x1=0
               BC(j,:) = -n2 * sin(BCfreq*Time) * xc
            case (6) ! Sinusoidal heave
               BC(j,:) = -n2 * cos(BCfreq*Time)
            case (7) ! Sinusoidal heave
               BC(j,:) = -n2 * sin(BCfreq*Time)
            case (8) ! Bose calc, th=0.4=alf/(f h), alf axis=1, alf leads h pi/2
               BC(j,:) = -n2 * (                                    &
                  0.4_dp*BCfreq**2*(xc-1.0_dp) * cos(BCfreq*Time) - &
                  0.6_dp*BCfreq                * sin(BCfreq*Time) )
            case default
               print *, "Error: invalid BCmode"
               stop
         end select
      end if

      do i = 1, n    ! Observer Loop
         r   = Xo(i,:) - Xs(j,:)
         r(1) = r(1) / Beta        ! PGtrans
         x   =  r(1)*ex(1) + r(2)*ex(2)
         z   =  r(1)*ez(1) + r(2)*ez(2)
         r12 =  x**2     + z**2
         r22 = (x-x2)**2 + z**2
         r1  = sqrt(r12)
         r2  = sqrt(r22)
         rc(1) = (Xo(i,1) - 0.5_dp*(Xs(j,1)+Xs(j+1,1))) / Beta
         rc(2) = (Xo(i,2) - 0.5_dp*(Xs(j,2)+Xs(j+1,2)))
         mrc = sqrt(rc(1)**2 + rc(2)**2)
         th1 = atan2(z,x)
         th2 = atan2(z,x-x2)

         G%d(i,j)     = 0.0_dp
         G%s(i,j)     = 0.0_dp
         G%r(i,j)     = 0.0_dp
         G%rF(i,j)    = 0.0_dp
         G%th(i,j)    = 0.0_dp
         G%thp(i,j)   = 0.0_dp
         G%thpF(i,j)  = 0.0_dp
         G%vd(i,j,:)  = 0.0_dp
         G%vs(i,j,:)  = 0.0_dp
         G%vr(i,j,:)  = 0.0_dp
         G%vrF(i,j,:) = 0.0_dp
         if(r1 > tol)then
            G%s(i,j)    = G%s(i,j)    +  x    *log(r1)
            G%vs(i,j,1) = G%vs(i,j,1) + log(r1)
            !G%vd(i,j,1) = G%vd(i,j,1) -  z    /(r12+delta**2)
            !G%vd(i,j,2) = G%vd(i,j,2) +  x    /(r12+delta**2)
         else
            th1 = 0.0_dp
            th2 = pi
         end if
         if(r2 > tol)then
            G%s(i,j)    = G%s(i,j)    - (x-x2)*log(r2)
            G%vs(i,j,1) = G%vs(i,j,1) - log(r2)
            !G%vd(i,j,1) = G%vd(i,j,1) +  z    /(r22+delta**2)
            !G%vd(i,j,2) = G%vd(i,j,2) - (x-x2)/(r22+delta**2)
         else
            th1 = 0.0_dp
            th2 = pi
         end if
         G%d(i,j)    = th2 - th1
         G%s(i,j)    = G%s(i,j)    + z*G%d(i,j)

         G%th(i,j)   = Mach/Beta*(mrc - Mach*rc(1))
         if(Mach > 0.0_dp)then
            G%thp(i,j)  = G%th(i,j) + dn*dt
            G%r(i,j)    = -Mach*Mn1/camp * G%s(i,j) / (dn*dt)

            if(Figueiredo)then
               rcp  = rc - ez*dn*x2
               mrcp = sqrt(rcp(1)**2 + rcp(2)**2)
               G%thpF(i,j)  = Mach/Beta*(mrcp - Mach*rc(1))
               G%rF(i,j)    = G%s(i,j)/(dn*x2)
            else
               G%r(i,j) = G%r(i,j) + Mach*mrc/Beta * G%d(i,j) / (dn*dt)
            end if
         end if

         G%vs(i,j,2) = G%d(i,j)
         G%vd(i,j,1) = G%vd(i,j,1)                      - &
            (Xo(i,2) - Xs(j+1,2))/(Beta*(r22+delta**2)) + &
            (Xo(i,2) - Xs( j ,2))/(Beta*(r12+delta**2))
         G%vd(i,j,2) = G%vd(i,j,2)                      + &
            (Xo(i,1) - Xs(j+1,1))/(Beta*(r22+delta**2)) - &
            (Xo(i,1) - Xs( j ,1))/(Beta*(r12+delta**2))

         G%s(i,j)    = G%s(i,j)    * camp
         G%vs(i,j,:) = G%vs(i,j,:) * camp


         ! Transform velocities back into physical plane
         gvs(1) =  n2*G%vs(i,j,1) + n1*G%vs(i,j,2)
         gvs(2) = -n1*G%vs(i,j,1) + n2*G%vs(i,j,2)
         !gvd(1) =  n2*G%vd(i,j,1) + n1*G%vd(i,j,2)
         !gvd(2) = -n1*G%vd(i,j,1) + n2*G%vd(i,j,2)
         G%vs(i,j,1) = gvs(1)
         G%vs(i,j,2) = gvs(2)
         !G%vd(i,j,1) = gvd(1)
         !G%vd(i,j,2) = gvd(2)

         if(Mach > 0.0_dp)then
            G%vr(i,j,:) = -Mach*Mn1/(camp**2) * G%vs(i,j,:) / (dn*dt)
            if(Figueiredo)then
               G%vrF(i,j,:) = G%vs(i,j,:) / (dn*x2)
            else
               G%vr(i,j,:) = G%vr(i,j,:) + Mach*mrc/Beta * G%vd(i,j,:) / (dn*dt)
            end if
         end if
      end do
   end do

end subroutine CalcGreens


end module Greens_mod
