module Globals

   use TypeKinds
   implicit none
   private

   ! Global parameters
   real(kind=dp), public, parameter :: pi = 3.14159265358979_dp, twopi = 6.28318530717959_dp
   real(kind=dp), public, parameter :: tol = 1.0e-7_dp        ! epsilon
   real(kind=dp), public, parameter :: D_Start_Vortex = 50.0_dp ! length of infinite vortex
                                                              ! segments for steady runs
   integer, public       :: Timeiter           ! current time iteration
   integer, public       :: Npieces            ! no solid body surfaces
   integer, public       :: BCmode             ! type of boundary condition
   integer, public       :: TEcond             ! alt Kutta Condition mode

   real(kind=dp), public, pointer, dimension(:) :: Time   ! time vector for entire calculation
   real(kind=dp), public :: Mach      ! freestream mach number
   real(kind=dp), public :: Beta      ! compressibility factor = sqrt(1-m^2)
   real(kind=dp), public :: dn        ! for d\phi/dn in timepole influence
   real(kind=dp), public :: BCfreq    ! frequency of oscillating boundary condition
   real(kind=dp), public :: delta     ! regularization parameter
   real(kind=dp), public :: oldtime    ! dt ofsteady state run - for restart 

   logical, public, save :: ot      = .false.   ! flag for restart 
   logical, public, save :: DoField       = .true.   ! perform acoustic field evaluation
   logical, public, save :: DoVort        = .true.   ! perform vortex source calc
   logical, public, save :: FreeVort        = .true.  ! perform real vortex evolv 
   logical, public, save :: Fixed_Wake    = .false.  ! sets v_wake = \hat x_1
   logical, public, save :: BisectorWake  = .true.   ! wake sheds with bisector (t) or freestream (f)
   logical, public, save :: Figueiredo    = .true.   ! Alt Ratelet discretization
   logical, public, save :: Steady        = .false.  ! pushes starting vortex "infinitely" far away
   logical, public, save :: WriteDataFile = .false.  ! Write BEM_data.m or not

   character(len=max_str_len), public :: Filetag  ! for differentiating i/o files of diff runs
   character(len=max_str_len), public :: GeomPath ! directory of non-analytical airfoil geometries

end module Globals
