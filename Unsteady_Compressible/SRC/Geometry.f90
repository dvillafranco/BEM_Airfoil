module Geometry

   use TypeKinds
   use Globals
   use IO
   implicit none
   private
   
   private :: getafgeom,naca,nacathick
   public  :: afoil, afgeom, aforient
   !public  :: winggeom, orient

   ! required parameters
   !
   !    ysymmetry: declared in global.f90; initialized in bem.f90
   !    pi       : declared/initialized in global.f90
   ! 

contains

   subroutine afoil(tok,n, xu,xl,xc,AFspacing)

      character(len=*), intent(inout)    :: tok
      integer, intent(in)                :: n
      real(kind=dp), dimension(:,:) , pointer :: xu,xl,xc
      character(len=*), intent(in)       :: AFspacing

      character(len=max_str_len)         :: tag, tok1, tok2, tok3, tok4
      integer, dimension(1)              :: idx
      integer                            :: ii
      integer, dimension(2)              :: err
      real(kind=dp), dimension(n)        :: x
      real(kind=dp)                      :: tau, lo,hi
      logical                            :: docut

      err = 0
      docut = .false.
      if(associated(xu)) then
         deallocate( xu, xc, xl, stat=err(1) )
      end if
      nullify(xu,xc,xl)
      allocate( xu(n,2), xl(n,2), xc(n,2), stat=err(2) )
      if(any(err/=0)) then
         print *, "error: cannot de/allocate in afoil", err
         stop
      end if

      if(AFspacing(1:6)=="linear" .or. AFspacing(1:6)=="LINEAR")then
         x = (/ ( ( 1.0_dp - real(ii,kind=dp)/real(n-1,kind=dp) ), ii=0,n-1) /)
      else
         x = (/ ( ( cos(ii*pi/(n-1)) + 1 ) / 2 , ii=0,n-1) /)
      end if

      call strtok(tok ," ",tag,tok1)
      tok2 = tag
      call strtok(tok2," ",tag,tok2)
      tok3 = tag
      call strtok(tok3," ",tag,tok3)
      docut = .false.
      if (index(tok3,"cut") /= 0) then     ! cutout for te flap
         docut = .true.
         tok3  = tag
         call strtok(tok3," ",tag,tok3)
         if(tok3=="") then
            lo = 0.6_dp
            hi = 0.8_dp
         else
            lo    = str2num(tok3)
            tok4  = tag
            call strtok(tok4," ",tag,tok4)
            if(tok4=="") then
               hi = 0.8_dp
            else
               hi = str2num(tok4)
            end if
         end if

         idx         = minloc(lo - x, mask = x < lo) - 1
         x(1:idx(1)) = hi - (hi-lo) * (/ (real(ii)/real(idx(1)),ii=0,idx(1)-1) /)
      end if

      if(tok1(1:4) == "PARA") then
            tau = str2num(tok2)/200

            xu(:,1) = x
            xu(:,2) = sqrt(((1+tau**2)/(2*tau))**2-(2*x-1)**2) - (1-tau*tau)/(2*tau)
            xl(:,1) = x
            xl(:,2) = -xu(:,2)
            xc(:,1) = x
            xc(:,2) = 0.0_dp

      elseif (tok1(1:4) == "FLAT") then
            ! tau = str2num(tok2)/200
            ! x = (/ (real(ii)/real(n-1),ii=n-1,0,-1) /)

            xu(:,1) = x
            xu(:,2) = 0.0_dp
            xl      = xu
            xc      = xu

      elseif (tok1(1:4) == "NACA") then
            call naca(x,n, tok2, xu,xl,xc)
      else
            call naca(x,n, tok1, xu,xl,xc)
      end if

      if (docut) then
         xl(1:idx(1),2) = xl(1:idx(1),2) + (xu(1:idx(1),2)-xl(1:idx(1),2)) * &
              0.5_dp*(1.0_dp - cos(pi*(x(1:idx(1))-lo)/(hi-lo)))
      end if

   end subroutine afoil
   
   subroutine afgeom(tok, n, x, idxwake, AFspacing)

      character(len=*), intent(inout)        :: tok
      integer, intent(in)                    :: n
      real(kind=dp), dimension(:,:), pointer :: x
      integer, dimension(:,:), pointer       :: idxwake
      character(len=*), intent(in)           :: AFspacing
      character(len=max_str_len)             :: tag,tok1
      real(kind=dp), dimension(:,:), pointer :: xu, xl, xc
   
      nullify(xu,xl,xc,x,idxwake)
      call strtok(tok ," ",tag,tok1)
      if(tok1(1:4) == "PARA" .or. tok1(1:4) == "NACA") then
         call afoil(tok,n,xu,xl,xc,AFspacing)
         allocate(x(2*n-1,3),idxwake(1,3))
         x = 0.0_dp
         x(1:  n  ,(/ 1,3 /)) = xu
         x(n:2*n-1,(/ 1,3 /)) = xl(n:1:-1,:)
         idxwake(1,:) = (/ 1,2*n-2,1 /)      ! Kutta shedding
         !idxwake(1,:) = (/ 1,2*n-2,2 /)      ! bruggeman shedding
         deallocate(xu,xl,xc)
         nullify(xu,xl,xc)
      elseif(tok1(1:4) == "FLAT") then
         call afoil(tok,n,xu,xl,xc,AFspacing)
         allocate(x(n,3),idxwake(1,2))
         x = xu
         idxwake(1,:) = (/ 1,1 /)
         deallocate(xu,xl,xc)
         nullify(xu,xl,xc)
      else
         call getafgeom(tok1,x,idxwake)
      end if

   end subroutine afgeom

   subroutine getafgeom(tok,x,idxwake)
   ! reads airfoil geomtry data from file.  the geometry is split up into 
   ! vector segments.  the number of segments corresponds to the number of 
   ! vortex sheet wakes that will be generated.  the location of the wake
   ! generation is determined by the break between segments.  
   ! idxwake(i,j) indexes the j=1|2 panels adjacent to the ith wake.  
   !    idxwake(:,3) determines the mode of wake shedding
   !                 = 1 -> kutta condition
   !                 = 2 -> bruggeman method

      character(len=*) , intent(in)          :: tok
      real(kind=dp), dimension(:,:), pointer :: x
      integer, dimension(:,:), pointer       :: idxwake
      real(kind=dp), dimension(:,:), allocatable :: xtmp
      integer                        :: nseg, np, err, i,j, lastidx
   
      lastidx = 0
      open(unit=2, file=trim(GeomPath)//"/"//trim(tok), &
           action="read", status="old", iostat=err)
      if(err/=0) then
         print *, "getafgeom error: cannot open ",trim(GeomPath)//"/"//trim(tok)
         stop
      end if

      read(unit=2,fmt=*) nseg
      allocate(idxwake(nseg,3))
      do i = 1, nseg
         read(unit=2,fmt=*,iostat=err) np, idxwake(i,3)
         if(err/=0) then
            print *, "getafgeom i/o error: read np, wake mode"
            stop
         end if
         if(associated(x)) then
            deallocate(x)
         end if
         nullify(x)
         allocate(x(lastidx+np,3))
         x = 0.0_dp
         if(allocated(xtmp)) then
            x(1:size(xtmp,1),:) = xtmp
         end if
         do j = 1, np
            read(unit=2,fmt=*,iostat=err) x(j+lastidx,1), x(j+lastidx,3) 
            if(err/=0) then
               print *, "getafgeom i/o error: read x, z"
               stop
            end if
         end do
         lastidx      = lastidx + np
         idxwake(i,1) = lastidx
         idxwake(i,2) = lastidx - 1
         if(allocated(xtmp)) then
            deallocate(xtmp)
         end if
         allocate(xtmp(lastidx,3))
         xtmp    = x
      end do
      idxwake(nseg,1) = 1
      deallocate(xtmp)
      close(unit=2)

   end subroutine getafgeom

   subroutine naca(x,n,tag, xu,xl,xc)
      real(kind=dp), dimension(:)  , intent(in)  :: x
      integer                      , intent(in)  :: n
      character(len=*)             , intent(in)  :: tag
      real(kind=dp), dimension(:,:), intent(out) :: xu,xl,xc
      real(kind=dp), dimension(n)      :: yc, th, yt
      real(kind=dp)                    :: m, p, t, rt, k1

      if(size(x)/=n .or. size(xu,1)/=n .or. size(xl,1)/=n .or. size(xc,1)/=n)then
         print *, "naca ERROR: wrong sizes for x"
         stop
      end if

      select case (len_trim(tag))

         case (4)
            m = str2num(tag(1:1))/100.0_dp    ! max ordinate of mean line (scaled on chord)
            p = str2num(tag(2:2))/10.0_dp     ! chordwise position of max ordinate
            t = str2num(tag(3:4))/100.0_dp    ! max thickness

            call nacathick(x,t, yt,rt)

            if (p > 0.0_dp) then
               where (x < p)
                  yc = m/p**2*(2*p-x)*x
                  th = 2*m/p**2*(p-x)
               elsewhere
                  yc  = m/(1-p)**2*(1-2*p+x*(2*p-x))
                  th  = 2*m/(1-p)**2*(p-x)
               end where
            else
               yc = 0*x
               th = 0*x
            end if

         case (5)
            if     (tag(1:3) == "210") then
               p=0.05_dp
               m=0.0580_dp
               k1=361.4_dp
            else if(tag(1:3) == "220") then
               p=0.10_dp
               m=0.1260_dp
               k1=51.64_dp
            else if(tag(1:3) == "230") then
               p=0.15_dp
               m=0.2025_dp
               k1=15.957_dp
            else if(tag(1:3) == "240") then
               p=0.20_dp
               m=0.2900_dp
               k1=6.643_dp
            else if(tag(1:3) == "250") then
               p=0.25_dp
               m=0.3910_dp
               k1=3.230_dp
            else
               print *, "incorrect mean-line designation for naca 5-series"
               stop
            end if

            t = str2num(tag(4:5))/100   ! max thickness

            call nacathick(x,t, yt,rt)

            if (p > 0) then
               where (x <  p)
                  yc = k1/6*(((x-3*m)*x+m*m*(3-m))*x)
                  th =atan( k1/6*((3*x-6*m)*x+m*m*(3-m)))
               elsewhere
                  yc  = k1/6*m**3*(1-x)
                  th  =atan(-k1/6*m**3 + 0*x)
               end where
            else
               yc = 0*x
               th = yc
            end if

         case default
            print *, "error: naca series tag"
            stop

      end select

      xu(:,1) = (x-yt*sin(th))
      xu(:,2) = (yc+yt*cos(th))
      xl(:,1) = (x+yt*sin(th))
      xl(:,2) = (yc-yt*cos(th))
      xc(:,1) = x
      xc(:,2) = yc

   end subroutine naca


   subroutine nacathick(x,t, yt,rt)

      real(kind=dp), dimension(:), intent(in)  :: x
      real(kind=dp)              , intent(in)  :: t
      real(kind=dp), dimension(:), intent(out) :: yt
      real(kind=dp)              , intent(out) :: rt
      real(kind=dp), dimension(5)              :: c

      if(size(yt)/=size(x)) then
         print *, "nacathick error: yt wrong size"
         stop
      end if

      !c = [0.29690,-0.12600,-0.35160,0.28430,-0.10150];
      c = (/  0.2953721665951587_dp , &
             -0.1192678795153606_dp , &
             -0.3742385808026736_dp , &
              0.3229485273880801_dp , &
             -0.1248142336652090_dp /)
      yt = t/0.2*(c(1)*sqrt(x)+x*(c(2)+x*(c(3)+x*(c(4)+c(5)*x))))
      rt = 1.1019_dp*t*t
      where (x==1) 
         yt = 0.0_dp
      end where

   end subroutine nacathick


!    subroutine winggeom(wg,idxwake,aftype,ls,nc,ns,taper,sweep,dihedral,twist,mode)
! 
!    ! generates wing geometry for a wing element with linear sweep, taper, 
!    ! dihedral, and twist variations in the span.  
!    !
!    ! aftype   -- airfoil geometry of the root and hub sections
!    ! ls       -- span / chord length
!    ! nc       -- no. of elements in the chord direction
!    ! ns       -- no. of elements in the span direction
!    ! taper    -- taper ratio.  tip airfoil is scaled as (taper) * (root geometry)
!    !             note that the sweep is consequently affected.
!    ! sweep    -- downstream shift of tip airfoil geometry
!    ! dihedral -- upwards shift of tip airfoil geometry
!    ! twist    -- angle of attack (in degrees) of the tip section 
!    !             (relative to the root)
!    ! mode     -- "cosine' or 'linear" spacing of the spanwise airfoil sections
!    !
!    !
!    ! the output is a cell array of the wing geometry and a matrix of indices of 
!    ! panels adjacent to wakes.
!    ! "wg" is a rank-3 array of the geometry of patch i.
!    ! i.e.,
!    !     x = wg(:,:,1)
!    !     y = wg(:,:,2)
!    !     z = wg(:,:,3)
!    !
!    ! note: wing outputs the wing geometry relative to a unit length root chord
!    !       section at 0 degrees aoa.  to arbitrarily scale and orient the wing, 
!    !       use the function "orient(wing,...)"
!    !
! 
! 
!       real(kind=dp), dimension(:,:,:), pointer :: wg
!       integer      , dimension(:,:)  , pointer :: idxwake
!       character(len=*), optional, intent(in)   :: aftype
!       real(kind=dp)   , optional, intent(in)   :: ls
!       integer         , optional, intent(in)   :: nc, ns
!       real(kind=dp)   , optional, intent(in)   :: taper, sweep, dihedral, twist
!       character(len=*), optional, intent(in)   :: mode
! 
!       character(len=max_str_len)  :: aftype2, mode2
!       real(kind=dp)               :: ls2, taper2, sweep2, dihedral2, twist2
!       integer                     :: nc2, ns2, ns3
! 
!       real(kind=dp), dimension(:,:)  , pointer :: xr, xt, rrt
!       real(kind=dp), dimension(:)    , pointer :: t
!       real(kind=dp)          :: ang
!       integer                :: n, i,j, err
! 
!       nullify(xr,xt,rrt,t)
!       err = 0
!       if(associated(wg)) then
!          deallocate(wg, stat=err)
!       end if
!       nullify(wg)
!       if(err/=0) then
!          print *, "error: cannot deallocate in wing"
!          stop
!       end if
! 
!       nullify(xr)
! 
!         aftype2 = "NACA 0012"
!             ls2 = 6
!             nc2 = 20
!             ns2 = 10
!          taper2 = 1.0_dp
!          sweep2 = 0.0_dp
!       dihedral2 = 0.0_dp
!          twist2 = 0.0_dp
!           mode2 = "cosine"
! 
!       if (present(  aftype)) then
!            aftype2 =   aftype
!       endif
!       if (present(      ls)) then
!                ls2 =       ls
!       end if
!       if (present(      nc)) then
!                nc2 =       nc
!       end if
!       if (present(      ns)) then
!                ns2 =       ns
!       end if
!       if (present(   taper)) then
!             taper2 =    taper
!       end if
!       if (present(   sweep)) then
!             sweep2 =    sweep
!       end if
!       if (present(dihedral)) then
!          dihedral2 = dihedral
!       end if
!       if (present(   twist)) then
!             twist2 =    twist
!       end if
!       if (present(    mode)) then
!              mode2 =     mode
!       end if
! 
!       ang = -twist2 * pi/180
! 
!       call afgeom(aftype2, nc2, xr, idxwake)
! 
!       n = size(xr,1)
!       if(ysymmetry) then
!          ns3 = 0
!       else
!          ns3 = -ns2
!       end if
!       allocate( xt(n,3), rrt(n,3), t(0:ns2) , wg(n, ns3:ns2, 3), stat=err )
!       if(err/=0) then
!          print *, "error: cannot allocate xr in wing"
!          stop
!       end if
! 
!       xt = (xr * taper2) +                     &
!            ( real(ones(n),kind=dp)     .o.     &
!               (/ sweep2+0.5_dp*(1.0_dp-taper2), 1.0_dp, dihedral2 /) )
! 
!       xt = matmul(xt , reshape(                   &
!             (/ cos(ang) , 0.0_dp , -sin(ang)    , &
!                0.0_dp   , 1.0_dp , 0.0_dp       , &
!                sin(ang) , 0.0_dp ,  cos(ang) /) , (/ 3, 3 /) ))
!       
!       rrt  = xt - xr
!       t    = real((/ (i,i=0,ns2) /), kind=dp) / ns2
! 
!       if (mode(1:3)=="cos") then
!          t = sin(pi*t/2)
!       end if
! 
!       do j = 0, ns2
!          wg(:,j,:) = xr * (1.0_dp - t(j)) + xt * t(j)
!       end do
!       if(.not.ysymmetry) then
!          wg(:,-ns2:-1,:) =  wg(:, ns2:1:-1,:)
!          wg(:,-ns2:-1,2) = -wg(:,-ns2:-1,2)
!       end if
!       
!       wg(:,:,2) = ls / 2.0_dp * wg(:,:,2)
! 
!       deallocate( xr,xt,rrt,t, stat=err )
!       nullify(xr,xt,rrt,t)
!       if(err/=0) then
!          print *, "error: cannot deallocate xr in wing"
!          stop
!       end if
! 
!    end subroutine winggeom
!    
!    
!    subroutine orient(wg,trans,scales,euler)
! 
!    ! translates, scales and orients the wing geometry using the 3-vectors, 
!    ! trans = [dx,dy,dz], scale = [sx,sy,sz], and euler = [ax,ay,az].
!    
! 
!       real(kind=dp), dimension(:,:,:), pointer :: wg
!       real(kind=dp), dimension(:), intent(in)  :: trans, scales, euler
!       real(kind=dp), dimension(:,:)  , pointer :: w
!       real(kind=dp), dimension(3)   :: ang
!       real(kind=dp), dimension(3,3) :: rx, ry, rz, r
!       integer      :: i, dims, err
! 
!       err = 0
!       if(size(trans)/=3 .or. size(scales)/=3 .or. size(euler)/=3) then
!          print *, "orient error: wrong sizes of inputs"
!          stop
!       end if
! 
!       ang = euler * pi/180
! 
!       rx  = reshape( (/ 1.0_dp      ,  0.0_dp       , 0.0_dp        , &
!                         0.0_dp      ,  cos(ang(1))  ,-sin(ang(1))   , &
!                         0.0_dp      ,  sin(ang(1))  , cos(ang(1)) /), (/3,3/) )
! 
!       ry  = reshape( (/ cos(ang(2)) ,  0.0_dp       ,-sin(ang(2))   , &
!                         0.0_dp      ,  1.0_dp       , 0.0_dp        , &
!                         sin(ang(2)) ,  0.0_dp       , cos(ang(2)) /), (/3,3/) )
! 
!       rz  = reshape( (/ cos(ang(3)) , -sin(ang(3))  , 0.0_dp        , &
!                         sin(ang(3)) ,  cos(ang(3))  , 0.0_dp        , &
!                         0.0_dp      ,  0.0_dp       , 1.0_dp      /), (/3,3/) )
! 
!       r   = matmul(rx, matmul(ry,rz))
! 
!       dims = size(wg,1)*size(wg,2)
!       allocate( w(dims,3), stat=err )
!       if(err/=0) then
!          print *, "error: cannot allocate w in orient"
!          stop
!       end if
!       w    = reshape(wg, (/ dims,3 /))
! 
!       do i = 1,dims
!          w(i,:) = w(i,:) * scales
!       end do
!       w    = matmul(w , r)
!       do i = 1,dims
!          w(i,:) = w(i,:) + trans
!       end do
! 
!       wg = reshape( w , (/ size(wg,1),size(wg,2),3 /) )
! 
!       deallocate( w, stat=err )
!       nullify(w)
!       if(err/=0) then
!          print *, "error: cannot deallocate w in orient"
!          stop
!       end if
! 
!    end subroutine orient

   subroutine aforient(X,t,s,a)
      real(kind=dp), dimension(:,:), intent(inout) :: X
      real(kind=dp), dimension(:)  , intent(in)    :: t,s
      real(kind=dp)                , intent(in)    :: a
      real(kind=dp), dimension(size(X,1),2)        :: Xt

      Xt = X
      Xt(:,1) = Xt(:,1) * s(1) + t(1)
      Xt(:,2) = Xt(:,2) * s(2) + t(2)
      X(:,1) = cos(a)*Xt(:,1) - sin(a)*Xt(:,2)
      X(:,2) = sin(a)*Xt(:,1) + cos(a)*Xt(:,2)
   end subroutine aforient

end module Geometry
