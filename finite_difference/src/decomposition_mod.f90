module decomposition_mod
  use region_mod, only: region_type
  implicit none

  ! The subdomain type captures everything we need to know about
  ! the definition of a subdomain within a larger, global domain.
  ! e.g. if '.' indicates a grid point in the global domain that
  ! does not belong to the current subdomain, 'h' is a halo grid
  ! point in the current subdomain and 'i' is an 'internal' grid
  ! point in the current subdomain, then:
  !
  ! 14 .  .  .  .  .  .  .  .  .  .  .  .
  ! 13 .  .  .  .  .  .  .  .  .  .  .  .
  ! 12 .  .  .  .  .  .  .  .  .  .  .  .
  ! 11 .  .  .  h  h  h  h  h  h  .  .  .
  ! 10 .  .  .  h  i  i  i  i  h  .  .  .
  !  9 .  .  .  h  i  i  i  i  h  .  .  .
  !  8 .  .  .  h  i  i  i  i  h  .  .  .
  !  7 .  .  .  h  i  i  i  i  h  .  .  .
  !  6 .  .  .  h  i  i  i  i  h  .  .  .
  !  5 .  .  .  h  i  i  i  i  h  .  .  .
  !  4 .  .  .  h  i  i  i  i  h  .  .  .
  !  3 .  .  .  h  h  h  h  h  h  .  .  .
  !  2 .  .  .  .  .  .  .  .  .  .  .  .
  !  1 .  .  .  .  .  .  .  .  .  .  .  .
  !    1  2  3  4  5  6  7  8  9  10 11 12
  !             Global x index -->
  !
  ! Position of the internal part of the subdomain in the global domain
  ! global%(xstart, ystart) = (5, 4), global%(xstop, ystop) = (8, 10)
  ! Width and height of the *whole* subdomain
  ! global%nx = 9 - 4 + 1 = 6
  ! global%ny = 11 - 3 + 1 = 9
  !
  ! Position of the 'internal' region (i.e. excluding halos, boundary
  ! points) within the sub-domain.
  ! internal%(xstart, ystart) = (2, 2), internal%(xstop, ystop) = (5, 8)
  ! Width and height of the internal region
  ! internal%nx = 5 - 2 + 1 = 4
  ! internal%ny = 8 - 1 + 1 = 7

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
