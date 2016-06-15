module IO

   use TypeKinds
   implicit none
   private
   
   private :: savematfile1, savematfile2, savematfile3
   private :: readmatfile1, readmatfile2, readmatfile3, numcols, getMat
   public  :: openmatfile, closematfile, savematfile, readmatfile
   public  :: elapsed, strcmp, str2num, num2str, strtok, sizecell

   integer             , parameter, private :: lform = 14, lfield = 13
   character(len=lform), parameter, private :: form = "(1p     e13.5)"
   character(len=1)    , parameter, private :: ldel="{", rdel="}", fs=","
   character(len=1)    , parameter, private :: vi="[", vo="]"

   interface savematfile
      module procedure savematfile1, savematfile2, savematfile3
   end interface 
   interface readmatfile
      module procedure readmatfile1, readmatfile2, readmatfile3
   end interface 

contains

   subroutine openmatfile(name,perm,fid,isasc)
      character(len=*), intent(in)  :: name
      character(len=*), intent(in)  :: perm
      integer, intent(out)          :: fid
      logical, intent(in), optional :: isasc

      integer                       :: stat
      character(len=max_str_len)    :: pos
      logical                       :: doBin
      
      doBin = .false.
      if(present(isasc))then
         if(.not.isasc)then
            doBin = .true.
         end if
      end if
      
      if(perm=="u") then
         pos = "append"
      else
         pos = "rewind"
      end if         

      fid = 50
      do
         if(doBin)then           ! BINARY mode
            open(unit=fid,file=trim(name)//"bin",action="readwrite", &
                 position=trim(pos),access="sequential",iostat=stat, &
                 status="unknown",form="unformatted")
         else                    ! ASCII mode
            open(unit=fid,file=trim(name)//"m",action="readwrite",   &
                 position=trim(pos),access="sequential",iostat=stat, &
                 status="unknown",form="formatted")
         end if
         if(stat==0) then
            exit
         end if
         fid = fid+1
         if(fid > 100)then
            print *, "OPENMATFILE WARNING: Cannot get unit between 50 and ",fid
         end if
      end do
   end subroutine openmatfile

   subroutine closematfile(fid) 
      integer, intent(in)      :: fid
      close(unit=fid)
   end subroutine closematfile

   subroutine savematfile1(fid,name,t, isasc)
      integer, intent(in)                     :: fid
      character(len=*), intent(in)            :: name
      real(kind=dp), dimension(:), intent(in) :: t
      logical, intent(in), optional           :: isasc

      logical                                 :: doBin
      
      doBin = .false.
      if(present(isasc))then
         if(.not.isasc)then
            doBin = .true.
         end if
      end if
      
      if(doBin)then
         write(unit=fid) size(t,1),t
      else
         write(unit=fid,fmt=*   ) trim(name)," = "//vi
         write(unit=fid,fmt=form) t
         write(unit=fid,fmt=*   ) " "//vo//";"
      end if
   end subroutine savematfile1

   subroutine savematfile2(fid,name,t, isasc)
      integer, intent(in)                       :: fid
      character(len=*), intent(in)              :: name
      real(kind=dp), dimension(:,:), intent(in) :: t
      logical, intent(in), optional             :: isasc

      character(len=lform)                      :: tform
      logical                                   :: doBin
      
      doBin = .false.
      if(present(isasc))then
         if(.not.isasc)then
            doBin = .true.
         end if
      end if
      
      if(doBin)then
         write(unit=fid) size(t,1),size(t,2),t
      else
         tform = form
         write(unit=tform(4:8),fmt="(i5)") size(t,2)
         write(unit=fid       ,fmt=*     ) trim(name)," = "//vi
         write(unit=fid       ,fmt=tform ) transpose(t)
         write(unit=fid       ,fmt=*     ) " "//vo//";"
      end if
   end subroutine savematfile2

   subroutine savematfile3(fid,name,t, isasc)
      integer, intent(in)                         :: fid
      character(len=*), intent(in)                :: name
      real(kind=dp), dimension(:,:,:), intent(in) :: t
      logical, intent(in), optional               :: isasc

      integer                    :: i,n
      character(len=lform)       :: tform
      character(len=max_str_len) :: name2
      logical                    :: doBin
      
      doBin = .false.
      if(present(isasc))then
         if(.not.isasc)then
            doBin = .true.
         end if
      end if
      
      if(doBin)then
         write(unit=fid) size(t,1),size(t,2),size(t,3),t
      else
         tform = form
         n = len_trim(name)
         name2            = ""
         name2(1   :n   ) = name
         name2(n+1 :n+5 ) = "(:,:,"
         name2(n+11:n+11) = ")"

         do i = 1, size(t,3)
            write(unit=name2(n+6:n+10),fmt="(i5)") i
            write(unit=tform(4:8),fmt="(i5)") size(t,2)
            write(unit=fid       ,fmt=*     ) trim(name2)," = "//vi
            write(unit=fid       ,fmt=tform ) transpose(t(:,:,i))
            write(unit=fid       ,fmt=*     ) " "//vo//";"
         end do
      end if
   end subroutine savematfile3

   subroutine readmatfile1(fid,name,t, isasc,cell)
      integer                      , intent(in)           :: fid
      character(len=*)             , intent(in)           :: name
      real(kind=dp), dimension(:)  , pointer              :: t
      logical                      , intent(in), optional :: isasc
      integer      , dimension(:)  , intent(in), optional :: cell

      real(kind=dp), dimension(:,:), pointer              :: tmp
      integer                                             :: maxrecl, N1
      logical                                             :: doBin
      
      doBin = .false.
      if(present(isasc))then
         if(.not.isasc)then
            doBin = .true.
         end if
      end if
      
      if(associated(t))then
         deallocate(t)
      end if
      nullify(t)
      if(doBin)then
         read(unit=fid) N1
         backspace(unit=fid)
         allocate(t(N1))
         read(unit=fid) N1,t
      else
         inquire(unit=fid,recl=maxrecl)
         if(present(cell))then
            call getMat(fid, maxrecl, name, tmp, cell)
         else
            call getMat(fid, maxrecl, name, tmp)
         end if
         if(size(tmp,2)/=1)then
            print *, "Error: quantity read in readmatfile1 is not a vector"
            stop
         end if
         allocate(t(size(tmp,1)))
         t = tmp(:,1)
         deallocate(tmp)
      end if
   end subroutine readmatfile1

   subroutine readmatfile2(fid,name,t, isasc,cell)
      integer                      , intent(in)           :: fid
      character(len=*)             , intent(in)           :: name
      real(kind=dp), dimension(:,:), pointer              :: t
      logical                      , intent(in), optional :: isasc
      integer      , dimension(:)  , intent(in), optional :: cell

      integer                                             :: maxrecl,N1,N2
      logical                                             :: doBin
      
      doBin = .false.
      if(present(isasc))then
         if(.not.isasc)then
            doBin = .true.
         end if
      end if
      if(associated(t))then
         deallocate(t)
      end if
      nullify(t)
      if(doBin)then
         read(unit=fid) N1,N2
         backspace(unit=fid)
         allocate(t(N1,N2))
         read(unit=fid) N1,N2,t
      else
         inquire(unit=fid,recl=maxrecl)
         if(present(cell))then
            call getMat(fid, maxrecl, name, t, cell)
         else
            call getMat(fid, maxrecl, name, t)
         end if
      end if
   end subroutine readmatfile2

   subroutine readmatfile3(fid,name,t,N, isasc,cell)
      integer                        , intent(in)           :: fid
      character(len=*)               , intent(in)           :: name
      real(kind=dp), dimension(:,:,:), pointer              :: t
      integer                        , intent(in)           :: N
      logical                        , intent(in), optional :: isasc
      integer      , dimension(:)    , intent(in), optional :: cell

      integer                                               :: i
      real(kind=dp), dimension(:,:)  , pointer              :: tmp
      integer                                               :: maxrecl,N1,N2,N3
      logical                                               :: doBin
      
      doBin = .false.
      if(present(isasc))then
         if(.not.isasc)then
            doBin = .true.
         end if
      end if
      if(associated(t))then
         deallocate(t)
      end if
      nullify(t)
      if(doBin)then
         read(unit=fid) N1,N2,N3
         backspace(unit=fid)
         allocate(t(N1,N2,N3))
         read(unit=fid) N1,N2,N3,t
      else
         inquire(unit=fid,recl=maxrecl)
         if(present(cell))then
            call getMat(fid, maxrecl, name, tmp, cell)
         else
            call getMat(fid, maxrecl, name, tmp)
         end if
         allocate(t(size(tmp,1),size(tmp,2),N))
         t(:,:,1) = tmp
         deallocate(tmp)
         do i = 2, N
            if(present(cell))then
               call getMat(fid, maxrecl, name, tmp, cell)
            else
               call getMat(fid, maxrecl, name, tmp)
            end if
            t(:,:,i) = tmp
            deallocate(tmp)
         end do
      end if
   end subroutine readmatfile3

   subroutine getMat(fid,maxl, name, mat, incell)
      integer                      , intent(in)           :: fid, maxl
      character(len=*)             , intent(in)           :: name
      real(kind=dp), dimension(:,:), pointer              :: mat
      integer      , dimension(:)  , intent(in), optional :: incell

      integer      , dimension(:)  , pointer              :: vec
      integer                                             :: i, M, N, idx, ios
      character(len=maxl)                                 :: line, cell,tok,rem
      character(len=lform)                                :: eform

      eform = form
      nullify(mat)
      do   ! find desired matrix
         read(unit=fid,fmt="(a)",iostat=ios) line
         if(ios/=0)then
            print *, "getMat Error: could not find ",name,"{",incell,"} in file #",fid
            stop
         end if
         !if(index(line,name) > 0)then
         call strtok(adjustl(line)," ",rem,tok)
         idx = scan(tok,"{}(),:[]=")
         if(idx > 0)then
            tok = tok(1:idx-1)
         endif
         if(strcmp(tok,name))then
            if(present(incell))then
               cell = line(index(line,ldel)+1:index(line,rdel)-1)
               N    = 1
               do   ! count number of elements in cell vector
                  idx = index(cell,fs)
                  if(idx > 0)then
                     N    = N + 1
                     cell = cell(:idx-1)//cell(idx+1:)   ! remove field separator
                  else
                     exit
                  end if
               end do
               allocate(vec(N))
               read(unit=cell,fmt=*) vec
               if(N==size(incell))then
                  if(all(vec==incell))then
                     deallocate(vec)
                     exit   ! matrix found
                  end if
               end if
               deallocate(vec)
            else
               exit   ! matrix found
            end if
         end if
      end do

      read(unit=fid,fmt="(a)") line
      if(index(line,"];") > 0)then
         print *, "Error: matrix ",name," is empty in fid ",fid
         stop
      end if
      N = numcols(line," ")     ! find no. of columns
      
      M = 1
      do   ! find no. of rows
         M = M + 1
         read(unit=fid,fmt=*) line
         if(index(line,"];") > 0)then
            M = M - 1
            exit
         end if
      end do
      do i = 1, M+1
         backspace(unit=fid)
      end do

      allocate(mat(M,N))
      if(N>1)then
         write(unit=eform(4:8),fmt="(i5)") N
      end if
      do i = 1, M
         read(unit=fid,fmt=eform) mat(i,:)
      end do

   end subroutine getMat

   function numcols(str,fs) result(N)
      character(len=*), intent(in)  :: str, fs
      integer                       :: N
      character(len=len_trim(str))  :: s
      integer                       :: idx
      
      s = trim(adjustl(str))
      N = 1
      do 
         idx = index(s,fs)
         if(idx <= 1 .or. idx >= len_trim(s))then
            exit
         end if
         N = N + 1
         s = adjustl(s(idx+1:))
      end do
   end function numcols

   subroutine sizecell(fid,ls,str,N)
      integer                , intent(in)    :: fid, ls
      character(len=*)       , intent(in)    :: str
      integer, dimension(:,:), intent(inout) :: N
      character(len=ls)                      :: line, rem,tok
      integer                                :: ios, ip
      integer, dimension(size(N,2))          :: tmp
      
      N = 0
      do
         do
            read(unit=fid,fmt="(a)",iostat=ios) line
            if(ios /= 0)then
               rewind(unit=fid)
               return
            end if
            call strtok(adjustl(line)," ",rem,tok)
            tok = tok(1:scan(tok,"{}(),:[]=")-1)
            if(strcmp(tok,str))then
               exit
            end if
         end do
         line = line(index(line,ldel)+1:index(line,rdel)-1)
         do ! remove field separator marks
            ip   = index(line,fs)
            if(ip > 1)then
               line = line(:ip-1)//line(ip+1:)
            else
               exit
            end if
         end do
         ip   = 0
         tmp  = 0
         if(size(N,2)==1)then    ! limits do not depend on patch index
            read(unit=line,fmt=*) ip
            where(ip > N(1,:))
               N(1,:) = ip
            end where
         else                    ! limits depend on patch index
            read(unit=line,fmt=*) ip,tmp
            where(tmp > N(ip,:))
               N(ip,:) = tmp
            end where
         end if
      end do
   end subroutine sizecell

   subroutine elapsed(t)
      real(kind=dp), intent(out) :: t
      real(kind=dp), save        :: tmp
      !integer, dimension(8)      :: vals
      integer                    :: ic,ir
      
      call system_clock(count = ic, count_rate = ir)
      t = real(ic,kind=dp) / real(ir,kind=dp) - tmp
      !call data_and_time(values=vals)
      !t = real(vals(8),kind=dp)/1000.0_dp + &
      !    real(vals(7),kind=dp)           + &
      !    real(vals(6),kind=dp)*60        + &
      !    real(vals(5),kind=dp)*3600      + &
      !    real(vals(4),kind=dp)*86400     - tmp
      tmp = t
             
   end subroutine elapsed

   function strcmp(s1,s2) result(same)
      character(len=*), intent(in) :: s1, s2
      logical                      :: same
      if(len_trim(adjustl(s1))==len_trim(adjustl(s2)))then
         same = trim(adjustl(s1))==trim(adjustl(s2))
      else
         same = .false.
      end if
   end function strcmp

   function str2num(str) result(t)
      character(len=*), intent(in)     :: str
      real(kind=dp)                    :: t
      character(len=12)                :: form
      integer                          :: l, d
      
      form = "(g0000.0000)"
      l = len_trim(str)
      d = l - scan(str,".")
      if (d==l) then
         d=0
      end if
      
      write(unit=form,fmt="('(g',i4.4,'.',i4.4,')')") l, d
      read(unit=str,fmt=form) t
      
      return
      
   end function str2num

   function num2str(num, form) result(t)
      real(kind=dp), intent(in)        :: num
      character(len=*), intent(in)     :: form
      character(len=max_str_len)       :: t

      write(unit=t,fmt=form) num
      
   end function num2str

   subroutine strtok(instr,tok,rem,outstr)
      character(len=*), intent(in)    :: instr
      character(len=*), intent(in)    :: tok
      character(len=*), intent(out)   :: rem
      character(len=*), intent(out)   :: outstr
      integer                         :: n
      
      !outstr = ""
      !rem    = ""
      n      = index(instr,tok)
      outstr = trim(adjustl(instr(1:n-1)))
      rem    = trim(adjustl(instr(n:)))

   end subroutine strtok
   
end module io

! test programs
! =============
!
! str2num
! -------
! program t
!    use vecutils
!    character(max_str_len) :: line
! 
!    do
!      read(*,*) line
!      print *, str2num(line)
!    end do
! 
! end program t
!
!
! strtok
! ------
! program t
!    use vecutils
!    character(max_str_len) :: line, tok, str,rem
! 
!    do
!       read(*,*) line, tok
!       call strtok(trim(line),trim(tok),str,rem)
!       print *, 'strtok(',trim(line),',',trim(tok),'):  ',&
!       &         '"',trim(str),'"    "', trim(rem),'"'
!    end do
! 
! end program t
