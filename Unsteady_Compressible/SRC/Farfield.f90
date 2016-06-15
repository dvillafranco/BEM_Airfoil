! |----------------------|
! |  Far-Field Routines  |
! |----------------------|

module Farfield
   use TypeKinds
   use BEM2d_typ
   use Globals
   use IO
   use BEM2d_mod
   use Greens_mod

   implicit none
   private

   ! Module procedures
   public  :: InitFField, CalcFField, WriteFField
   private :: CalcFField0, CalcFFieldM, CalcDeadTime

   ! Module data
   type(Greens) , private, dimension(:)    , pointer :: Gbff,Gwff
   real(kind=dp), private, dimension(:,:)  , pointer :: Xff, Phiff
   real(kind=dp), private, dimension(:,:,:), pointer :: Uff
   real(kind=dp), private                            :: Rff, DeadTime
   integer      , private                            :: Nff
   logical      , private                            :: WriteFFdat

contains

subroutine InitFField()
   integer :: err, i, N, fid
   
   open(unit=11,file="field.in",status="old",action="read",iostat=err)
   if(err/=0)then
      print *, "InitFField Error: cannot open field.in"
      stop
   end if
   read(unit=11,fmt=*) Rff, Nff, WriteFFdat
   close(unit=11)

   N = size(Time,1)

   if(associated(Xff))then
      deallocate(Xff)
   end if
   nullify(Xff)
   allocate(Xff(Nff,2),Phiff(Nff,N),Uff(Nff,2,N), stat=err)
   if(err/=0)then
      print *, "InitFField Error: cannot allocate ff data"
      stop
   end if

   ! Set far-field observer locations on a circle of radius Rff indexed in
   ! counter-clockwise orientation starting with the flyover direction
   Xff(:,1)=(/ (cos(2.0_dp*pi*(real(i,kind=dp)/real(Nff,kind=dp) - 0.25_dp)),i=0,Nff-1) /)
   Xff(:,2)=(/ (sin(2.0_dp*pi*(real(i,kind=dp)/real(Nff,kind=dp) - 0.25_dp)),i=0,Nff-1) /)
   Xff     = Rff * Xff

   Phiff = 0.0_dp
   Uff   = 0.0_dp

   call openmatfile("field_data."//trim(Filetag),"w", fid, WriteFFdat)
      call savematfile(fid, "time"   , Time, WriteFFdat)
      call savematfile(fid, "Xff"    , Xff , WriteFFdat)
   call closematfile(fid)

end subroutine InitFField

subroutine CalcDeadTime(W)
   type(Wing), intent(inout) :: W
   integer :: i

   if(associated(Gbff))then
      deallocate(Gbff,Gwff)
   end if
   nullify(Gbff,Gwff)
   allocate(Gbff(Npieces),Gwff(Npieces))

   DeadTime = huge(1.0_dp)
   do i = 1, Npieces
      call nullGreens(Gbff(i))
      call nullGreens(Gwff(i))
      call CalcGreens(W%B(i,1)%X, Xff, Gbff(i))
      DeadTime = min(DeadTime, minval(Gbff(i)%th))
   end do
   print *, "   Far-Field DeadTime = ", DeadTime
end subroutine CalcDeadTime

subroutine CalcFField(W, recalc)
   type(Wing), intent(inout) :: W
   logical   , intent(in)    :: recalc

   if(Timeiter==1)then
      call CalcDeadTime(W)
   end if
   if(Mach > 0.0_dp)then
      call CalcFFieldM(W, recalc)
   else
      call CalcFField0(W, recalc)
   end if
end subroutine CalcFField

subroutine CalcFField0(W, recalc)
   type(Wing), intent(inout) :: W
   logical   , intent(in)    :: recalc
   integer :: i, k

   do i = 1, Npieces
      if(Timeiter>1 .and. recalc)then
         call CalcGreens(W%B(i,Timeiter)%X, Xff, Gbff(i))
      end if
      call CalcGreens(W%W(i,Timeiter)%X, Xff, Gwff(i))

      Phiff(:,Timeiter) = 0.5_dp * ( &
         matmul(Gbff(i)%d,W%B(i,Timeiter)%phi ) + &
         matmul(Gbff(i)%s,W%B(i,Timeiter)%phin) + &
         matmul(Gwff(i)%d,W%W(i,Timeiter)%phi ) ) / twopi
      do k = 1, 2
         Uff(:,k,Timeiter) = 0.5_dp * ( &
            matmul(Gbff(i)%vd(:,:,k),W%B(i,Timeiter)%phi ) + &
            matmul(Gbff(i)%vd(:,:,k),W%B(i,Timeiter)%phin) + &
            matmul(Gwff(i)%vd(:,:,k),W%W(i,Timeiter)%phi ) ) / twopi
      end do
   end do

end subroutine CalcFField0

subroutine CalcFFieldM(W, recalc)
   type(Wing), intent(inout) :: W
   logical   , intent(in)    :: recalc
   real(kind=dp), dimension(:,:), allocatable    :: alf, src,dbl
   integer      , dimension(:,:), allocatable    :: ith
   integer :: i, k, Nel

   Phiff(:,Timeiter) = 0.0_dp
   Uff(:,:,Timeiter) = 0.0_dp
   do i = 1, Npieces
      if(Timeiter>1 .and. recalc)then
         call CalcGreens(W%B(i,Timeiter)%X, Xff, Gbff(i))
      end if
      call CalcGreens(W%W(i,Timeiter)%X, Xff, Gwff(i))

      Nel = size(W%B(i,1)%X,1) - 1
      allocate(alf(Nff,Nel), ith(Nff,Nel), src(Nff,Nel), dbl(Nff,Nel))
      alf = (Gbff(i)%th-DeadTime) / Time(1)
      ith = int(alf)
      alf = alf - real(ith,kind=dp)
      call sourceVector(W%B,i, ith,alf, src,dbl, .true.)
      Phiff(:,Timeiter) =  Phiff(:,Timeiter) + sum( &
            (Gbff(i)%d+Gbff(i)%r)*dbl - Gbff(i)%s*src, dim=2)
      do k = 1, 2
         Uff(:,k,Timeiter) = Uff(:,k,Timeiter) + sum( &
            (Gbff(i)%vd(:,:,k)+Gbff(i)%vr(:,:,k))*dbl - Gbff(i)%vs(:,:,k)*src, dim=2)
      end do
      if(Figueiredo)then
         Phiff(:,Timeiter) =  Phiff(:,Timeiter) + sum( &
               Gbff(i)%rF*dbl, dim=2)
         do k = 1, 2
            Uff(:,k,Timeiter) = Uff(:,k,Timeiter) + sum( &
               Gbff(i)%vrF(:,:,k)*dbl, dim=2)
         end do
      end if

      ! Add Euler forward differential
      alf = (Gbff(i)%thp-DeadTime) / Time(1)
      ith = int(alf)
      alf = alf - real(ith,kind=dp)
      call sourceVector(W%B,i, ith,alf, src,dbl, .true.)
      Phiff(:,Timeiter) = Phiff(:,Timeiter) - sum(Gbff(i)%r*dbl, dim=2)
      do k = 1, 2
         Uff(:,k,Timeiter) = Uff(:,k,Timeiter) - &
            sum(Gbff(i)%vr(:,:,k)*dbl, dim=2)
      end do

      ! Add Figueiredo influences
      if(Figueiredo)then
         alf = (Gbff(i)%thpF-DeadTime) / Time(1)
         ith = int(alf)
         alf = alf - real(ith,kind=dp)
         call sourceVector(W%B,i, ith,alf, src,dbl, .true.)
         Phiff(:,Timeiter) = Phiff(:,Timeiter) - sum(Gbff(i)%rF*dbl, dim=2)
         do k = 1, 2
            Uff(:,k,Timeiter) = Uff(:,k,Timeiter) - &
               sum(Gbff(i)%vrF(:,:,k)*dbl, dim=2)
         end do
      end if

      deallocate(alf,ith, src,dbl)

      ! Add wake contribution
      Nel = Timeiter
      allocate(alf(Nff,Nel), ith(Nff,Nel), dbl(Nff,Nel))
      alf = (Gwff(i)%th-DeadTime) / Time(1)
      ith = int(alf)
      alf = alf - real(ith,kind=dp)
      call wakeSource(W%W,i, ith,alf, dbl)
      Phiff(:,Timeiter) = Phiff(:,Timeiter) + sum((Gwff(i)%d+Gwff(i)%r)*dbl, dim=2)
      do k = 1, 2
         Uff(:,k,Timeiter) = Uff(:,k,Timeiter) + &
            sum((Gwff(i)%vd(:,:,k)+Gwff(i)%vr(:,:,k))*dbl, dim=2)
      end do

      alf = (Gwff(i)%thp-DeadTime) / Time(1)
      ith = int(alf)
      alf = alf - real(ith,kind=dp)
      call wakeSource(W%W,i, ith,alf, dbl)
      Phiff(:,Timeiter) = Phiff(:,Timeiter) - sum(Gwff(i)%r*dbl, dim=2)
      do k = 1, 2
         Uff(:,k,Timeiter) = Uff(:,k,Timeiter) - &
            sum(Gwff(i)%vr(:,:,k)*dbl, dim=2)
      end do

      if(Figueiredo)then
         alf = (Gwff(i)%thpF-DeadTime) / Time(1)
         ith = int(alf)
         alf = alf - real(ith,kind=dp)
         call wakeSource(W%W,i, ith,alf, dbl)
         Phiff(:,Timeiter) = Phiff(:,Timeiter) - sum(Gwff(i)%rF*dbl, dim=2)
         do k = 1, 2
            Uff(:,k,Timeiter) = Uff(:,k,Timeiter) - &
               sum(Gwff(i)%vrF(:,:,k)*dbl, dim=2)
         end do
      end if

      deallocate(alf,ith, dbl)

   end do

   Phiff(:,Timeiter) = 0.5_dp * Phiff(:,Timeiter) / twopi
   Uff(:,:,Timeiter) = 0.5_dp * Uff(:,:,Timeiter) / twopi

end subroutine CalcFFieldM


subroutine WriteFField()
   integer :: fid

   call openmatfile("field_data."//trim(Filetag),"u", fid, WriteFFdat)
   call savematfile(fid,"Phiff",Phiff   , WriteFFdat)
   call savematfile(fid,"Uff",Uff(:,1,:), WriteFFdat)
   call savematfile(fid,"Vff",Uff(:,2,:), WriteFFdat)
   call closematfile(fid)
end subroutine WriteFField

end module Farfield
