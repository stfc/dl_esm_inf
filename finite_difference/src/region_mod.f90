module region_mod
  implicit none

  !> Specify a region on the simulation grid
  type :: region_type
     !> Extent of this region in x, y and z
     integer :: nx, ny, nz
     integer :: xstart, xstop
     integer :: ystart, ystop
     integer :: zstart, zstop
  end type region_type

  interface region_type
     module procedure region_constructor
  end interface region_type

contains

  function region_constructor() result(self)
    implicit none
    type(region_type) :: self

    self%nx = 0
    self%ny = 0
    self%nz = 0

    self%xstart = 0
    self%xstop  = 0
    self%ystart = 0
    self%ystop  = 0
    self%zstart = 0
    self%zstop  = 0

  end function region_constructor

end module region_mod
