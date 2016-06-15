module TypeKinds
   implicit none
   
! standard f90 kind types
   intrinsic selected_real_kind
   integer, public, parameter :: i8b  = 8
   integer, public, parameter :: i4b  = 4
   integer, public, parameter :: sp   = selected_real_kind(4)
   integer, public, parameter :: dp   = selected_real_kind(8)
   integer, public, parameter :: csp  = 8
   integer, public, parameter :: cdp  = 16

   integer, public, parameter :: max_str_len = 80     ! maximum string length

! fplus (f95) kind types
!  integer, parameter :: i8b = 3
!  integer, parameter :: i4b = 3
!  integer, parameter :: sp  = 1
!  integer, parameter :: dp  = 2
!  integer, parameter :: csp = 1
!  integer, parameter :: cdp = 2

end module TypeKinds
