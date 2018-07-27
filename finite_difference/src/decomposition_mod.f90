module decomposition_mod
  use region_mod, only: region_type
  implicit none

  !> Type encapsulating the information required to define a single
  !! sub-domain
  type :: subdomain_type
     !> The definition of this subdomain in terms of the global domain
     type(region_type) :: global
     !> The internal region of this subdomain (excluding halo and
     !! boundary points)
     type(region_type) :: internal
  end type subdomain_type

  !> Type encapsulating all information regarding a regular, 2D
  !! domain decomposition
  type :: decomposition_type
     !> Dimensions of the global domain that is decomposed
     integer :: global_nx, global_ny
     !> Dimensions of the grid of sub-domains
     integer :: nx, ny
     !> Number of sub-domains (=nx*ny)
     integer :: ndomains
     !> Max dimensions of any sub-domain
     integer :: max_width, max_height
     !> Array of the sub-domain definitions
     type(subdomain_type), allocatable :: subdomains(:)
     !> An MPI process may have more than one sub-domain
     integer, allocatable :: proc_subdomains(:,:)
  end type decomposition_type

end module decomposition_mod
