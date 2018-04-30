!> Module containing the definition of the sub-domain type.
module subdomain_mod
  use region_mod
  implicit none

  type :: subdomain_type
     !> Position of bottom-left and top-right corners of this sub-domain
     !! in the global domain.
     integer :: xstart, ystart, xstop, ystop
     !> Full width and height of this subdomain - includes all halo and
     !! boundary points. nx != xstop - xstart +1 because xstart and xstop
     !! are the coordinates of the internal region in the global domain.
     integer :: nx, ny
     !> The internal region of this subdomain (excluding halo and
     !! boundary points)
     type(region_type) :: internal
  end type subdomain_type

end module subdomain_mod
