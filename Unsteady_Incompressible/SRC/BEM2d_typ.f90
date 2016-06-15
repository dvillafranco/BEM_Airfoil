module BEM2d_typ
   use TypeKinds

   type, public :: Wing
   ! A Wing comprises a (B)ody and (W)ake surface and coll pts (Xc)
   ! idx is a running total of the no. of coll pts
   ! idx 1 => wing element ; idx 2 => Time
      type(Surf)   , dimension(:,:)  , pointer :: B, W
      real(kind=dp), dimension(:,:)  , pointer :: Xc
      integer      , dimension(:)    , pointer :: idx
   end type Wing

   type, public :: Surf
   ! A Surf is a clockwise oriented contour with nodes X and perturbation 
   ! potential phi and normal flux velocity phin at central collocation pts.
      real(kind=dp), dimension(:,:)  , pointer :: X
      real(kind=dp), dimension(:)    , pointer :: phi, phin
   end type Surf
end module BEM2d_typ

