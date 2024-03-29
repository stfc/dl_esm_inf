!------------------------------------------------------------------------------
! BSD 2-Clause License
! 
! Copyright (c) 2017-2019, Science and Technology Facilities Council
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
! 
! * Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
! 
! * Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!------------------------------------------------------------------------------
! Author: A. R. Porter, STFC Daresbury Laboratory
! Modified: J. Henrichs, Australian Bureau of Meteorology

!> Module for describing all aspects of a field (which exists on some
!! grid).
module field_mod
  use kind_params_mod
  use region_mod
  use halo_mod
  use grid_mod
  use gocean_mod, only: gocean_stop
  use tile_mod
  use parallel_mod, only: decomposition_type
  implicit none

  private

  ! Enumeration of grid-point types on the Arakawa C grid. A
  ! field lives on one of these types.
  integer, public, parameter :: GO_U_POINTS   = 0
  integer, public, parameter :: GO_V_POINTS   = 1
  integer, public, parameter :: GO_T_POINTS   = 2
  integer, public, parameter :: GO_F_POINTS   = 3
  !> A field that lives on all grid-points of the grid
  integer, public, parameter :: GO_ALL_POINTS = 4

  !> Abstract interfaces for C and Fortran subroutines to be implemented by
  ! the infrastructure user to retrieve data from, or push data to, devices
  ! using the desired programming model.
  !
  ! All interfaces include:
  !  - A host and a device pointer specifying where the data is copied from/to.
  !  - startx, starty, nx, and ny integers to specify if copying just a sub-
  !  region of the field is enough (the implementation may copy a larger part).
  !  - A blocking boolean to specify if the data copy operation should have
  !  finished on the routine return or just been dispatched.
  !
  abstract interface
    subroutine read_from_device_c_interface(from, to, startx, starty, nx, ny, blocking)
        use iso_c_binding, only: c_ptr, c_int, c_bool
        type(c_ptr), intent(in), value :: from
        type(c_ptr), intent(in), value :: to
        integer(c_int), intent(in), value :: startx, starty, nx, ny
        logical(c_bool), intent(in), value :: blocking
    end subroutine read_from_device_c_interface
  end interface

  abstract interface
    subroutine read_from_device_f_interface(from, to, startx, starty, nx, ny, blocking)
        use iso_c_binding, only: c_ptr
        use kind_params_mod, only: go_wp
        type(c_ptr), intent(in) :: from
        real(go_wp), dimension(:,:), target, intent(inout) :: to
        integer, intent(in) ::  startx, starty, nx, ny
        logical, intent(in) :: blocking
    end subroutine read_from_device_f_interface
  end interface

  abstract interface
    subroutine write_to_device_c_interface(from, to, startx, starty, nx, ny, blocking)
        use iso_c_binding, only: c_ptr, c_int, c_bool
        type(c_ptr), intent(in), value :: from
        type(c_ptr), intent(in), value :: to
        integer(c_int), intent(in), value :: startx, starty, nx, ny
        logical(c_bool), intent(in), value :: blocking
    end subroutine write_to_device_c_interface
  end interface

  abstract interface
    subroutine write_to_device_f_interface(from, to, startx, starty, nx, ny, blocking)
        use iso_c_binding, only: c_ptr
        use kind_params_mod, only: go_wp
        real(go_wp), dimension(:,:), target, intent(in) :: from
        type(c_ptr), intent(in) :: to
        integer, intent(in) ::  startx, starty, nx, ny
        logical, intent(in) :: blocking
    end subroutine write_to_device_f_interface
  end interface


  !> The base field type. Intended to represent a global field
  !! such as the x-component of velocity.
  type, public :: field_type
     !> Which mesh points the field is defined upon
     integer :: defined_on
     !> The grid on which this field is defined
     type(grid_type), pointer :: grid
     !> The internal region of this field
     type(region_type) :: internal
     !> The whole region covered by this field - includes 
     !! boundary points
     type(region_type) :: whole
     !> The number of halo regions that this field has.
     !! Halo region values are not computed but copied
     !! from elsewhere.
     integer :: num_halos
     !> Array of objects describing the halos belonging to this field.
     type(halo_type), dimension(:), allocatable :: halo
     !> Whether the data for this field lives in a remote memory space
     !! (e.g. on a GPU)
     logical :: data_on_device
     !> Function pointers to the functions that read/write data from/to
     ! devices.
     procedure(read_from_device_c_interface), POINTER, nopass :: read_from_device_c
     procedure(read_from_device_f_interface), POINTER, nopass :: read_from_device_f
     procedure(write_to_device_c_interface), POINTER, nopass :: write_to_device_c
     procedure(write_to_device_f_interface), POINTER, nopass :: write_to_device_f

  end type field_type

  !> A real, 2D field.
  type, public, extends(field_type) :: r2d_field
     integer :: ntiles
     !> The dimensions of the tiles into which the field
     !! is sub-divided.
     type(tile_type), dimension(:), allocatable :: tile
     !> Array holding the actual field values. Eventually this should
     !! be private so that users are forced to access it via get_data().
     !! This will allow us to guarantee that the data is valid (e.g. in
     !! the case where values are computed on a separate accelerator
     !! device).
     real(go_wp), dimension(:,:), allocatable :: data
     !> Pointer to corresponding buffer on remote device (if any).
     !! This requires variables that are declared to be of this type
     !! have the 'target' attribute.
     type(c_ptr) :: device_ptr
   contains
     !> Setter for the data associated with this field
     procedure, pass :: set_data
     !> Getter for the data associated with this field. Fetches data
     !! from remote accelerator if necessary.
     procedure, pass :: get_data
     procedure, pass :: read_halo_from_device
     procedure, pass :: write_halo_to_device
     procedure, pass :: read_from_device
     procedure, pass :: write_to_device
     procedure, public :: halo_exchange
     procedure, public :: gather_inner_data
  end type r2d_field

  !> Interface for the copy_field operation. Overloaded to take
  !! an array or an r2d_field type.
  !! \todo Remove support for raw arrays from this interface.
  interface copy_field
     module procedure copy_2dfield_array, copy_2dfield_array_patch, &
                      copy_2dfield, copy_2dfield_patch
  end interface copy_field

  ! User-defined constructor for r2d_field type objects
  interface r2d_field
     module procedure r2d_field_constructor
  end interface r2d_field

  !> Interface for the field checksum operation. Overloaded to take either
  !! a field object or a 2D, real(go_wp) array.
  interface field_checksum
     module procedure fld_checksum, array_checksum
  end interface field_checksum

  interface free_field
     module procedure r2d_free_field
  end interface

  public copy_field
  public set_field
  public field_checksum
  public free_field

! Grid points on an Arakawa C grid with NE offset (i.e. the U,V and F pts
! immediately to the North and East of a T point share its grid indices) 
! are arranged like so:
!
!v(1,ny)----f(1,ny)---v(i-1,ny)--f(i-1,ny)--v(i,ny)----f(i,ny)--v(nx,ny)---f(nx,ny)  
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!T[1,ny]----u(1,ny)---T(i-1,ny)--u(i-1,ny)--T(i,ny)----u(i,ny)--T(nx,ny)---u(nx,ny)  
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!v(1,j)-----f(1,j)----v(i-1,j)---f(i-1,j)---v(i,j)-----f(i,j)---v(nx,j)----f(nx,j)   
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!T[1,j]-----u(1,j)----T(i-1,j)---u(i-1,j)---T(i,j)-----u(i,j)---T(nx,j)----u(nx,j)   
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!v(1,j-1)---f(1,j-1)--v(i-1,j-1)-f(i-1,j-1)-v(i,j-1)---f(i,j-1)-v(nx,j-1)--f(nx,j-1) 
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!T[1,j-1]---u(1,j-1)--T(i-1,j-1)-u(i-1,j-1)-T(i,j-1)---u(i,j-1)-T(nx,j-1)--u(nx,j-1) 
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!v(1,1)-----f(1,1)----v(i-1,1)---f(i-1,1)---v(i,1)-----f(i,1)---v(nx,1)----f(nx,1)   
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!T[1,1]     u(1,1)    T(i-1,1)---u(i-1,1)---T(i,1)-----u(i,1)---T(nx,1)----u(nx,1)   

  !> The no. of cols/rows used to define boundary data in the absence
  !! of periodic BCs.
  !! \todo remove this parameter and determine it from
  !! the supplied T mask.
  integer, public, parameter :: NBOUNDARY = 1

  !> Whether or not to 'tile' (sub-divide) field arrays
  logical, public, parameter :: TILED_FIELDS = .FALSE.

  !> The tiling decomposition to use for each field
  type(decomposition_type), save :: decomp

  !> Whether or not the tiling decomposition has been created yet
  logical, save :: tiling_initialised = .false.

contains

  !===================================================

  function r2d_field_constructor(grid,        &
                                 grid_points, &
                                 do_tile,     &
                                 init_global_data) result(self)
    use parallel_mod, only: go_decompose, on_master
!$    use omp_lib, only : omp_get_max_threads
    implicit none
    ! Arguments
    !> Pointer to the grid on which this field lives
    type(grid_type), intent(in), target  :: grid
    !> Which grid-point type the field is defined on
    integer,         intent(in)          :: grid_points
    !> If the field should be tiled among all threads, or if only
    !> a single field should be allocated (which is not currently
    !> supported by PSyclone)
    logical, intent(in), optional :: do_tile
    !> An optional global array with which the field in each domain
    !> will be  initialsed
    real(go_wp), dimension(:,:), intent(in), optional :: init_global_data
    ! Local declarations
    type(r2d_field), target :: self
    logical :: local_do_tiling = .false.
    integer :: ierr
    character(len=8) :: fld_type
    integer :: ji, jj
    !> The upper bounds actually used to allocate arrays (as opposed
    !! to the limits carried around with the field)
    integer :: upper_x_bound, upper_y_bound
    integer :: itile, nthreads, ntilex, ntiley, dx, dy

    if (present(do_tile)) then
       local_do_tiling = do_tile
    endif

    ! Set this field's grid pointer to point to the grid pointed to
    ! by the supplied grid_ptr argument
    self%grid => grid

    !> The data associated with this device is currently local
    !! to where we're executing
    self%data_on_device = .FALSE.

    !> The function pointers are initialized with NULL
    nullify(self%read_from_device_c)
    nullify(self%read_from_device_f)
    nullify(self%write_to_device_c)
    nullify(self%write_to_device_f)

    ! Set-up the limits of the 'internal' region of this field
    !
    call set_field_bounds(self,fld_type,grid_points)

    ! Dimensions of the grid of tiles. 
    if(.not. get_grid_dims(ntilex, ntiley) )then
       ntilex = 1
       ntiley = 1
    end if

    if (local_do_tiling) then
       nthreads = 1
       !> \TODO #38 We don't actually support building dl_esm_inf
       !! with OpenMP enabled!
!$     nthreads = omp_get_max_threads()

       if (.not. tiling_initialised) then
!$        write (*,"(/'Have ',I3,' OpenMP threads available.')") nthreads
          decomp = go_decompose(self%internal%nx, self%internal%ny, &
               nthreads, ntilex, ntiley)
          tiling_initialised = .true.
       end if

       self%ntiles = nthreads
       allocate(self%tile(self%ntiles), Stat=ierr)
       if(ierr /= 0)then
          call gocean_stop('r2d constructor failed to allocate tiling structures')
       end if
       do itile = 1, self%ntiles
          self%tile(itile)%whole = decomp%subdomains(itile)%global
          self%tile(itile)%internal = decomp%subdomains(itile)%internal
       end do
    else
       ! We're not tiling this field so set ntiles to zero
       self%ntiles = 0
    end if

    ! We allocate *all* fields to have the same extent as that
    ! of the grid. This enables the (Cray) compiler
    ! to safely evaluate code within if blocks that are
    ! checking for conditions at the boundary of the domain.
    ! Hence we use self%whole%{x,y}stop + 1...
    upper_x_bound = self%grid%nx
    upper_y_bound = self%grid%ny

    if (on_master()) then
        write(*, "('Allocating ',(A),' field with bounds: (',I1,':',I4, "// &
            "',',I1,':',I4,'), internal region is (',I1,':',I4, ',',I1,':',I4,')' )") &
            TRIM(ADJUSTL(fld_type)), 1, upper_x_bound, 1, upper_y_bound, &
             self%internal%xstart, self%internal%xstop, &
             self%internal%ystart, self%internal%ystop
    endif

    ! Allocating with a lower bound != 1 causes problems whenever
    ! array passed as assumed-shape dummy argument because lower
    ! bounds default to 1 in called unit (unless array declared as
    ! allocatable). However, all loops will be in the generated, middle
    ! layer and the generator knows the array bounds. This may give us the
    ! ability to solve this problem (by passing array bounds to the
    ! kernels).
    allocate(self%data(1:upper_x_bound, 1:upper_y_bound), &
                       Stat=ierr)
    if(ierr /= 0)then
       call gocean_stop('r2d_field_constructor: ERROR: failed to '// &
                        'allocate field')
    end if

    ! Since we're allocating the arrays to be larger than strictly
    ! required we explicitly set all elements to 0 in case the code
    ! does access 'out-of-bounds' elements during speculative
    ! execution. If we're running with OpenMP this also gives
    ! us the opportunity to do a 'first touch' policy to aid with
    ! memory<->thread locality...
    if (self%ntiles > 0) then
!$OMP PARALLEL DO schedule(runtime), default(none), &
!$OMP private(itile,ji,jj), shared(self)
       do itile = 1, self%ntiles
          do jj = self%tile(itile)%whole%ystart, self%tile(itile)%whole%ystop
             do ji = self%tile(itile)%whole%xstart, self%tile(itile)%whole%xstop
                self%data(ji,jj) = 0.0_go_wp
             end do
          end do
       end do
!$OMP END PARALLEL DO
    else
       self%data(:,:) = 0.0_go_wp
    end if

    if (present(init_global_data)) then
        dx = grid%subdomain%global%xstart - self%internal%xstart
        dy = grid%subdomain%global%ystart - self%internal%ystart

!$OMP PARALLEL DO schedule(runtime), default(none), &
!$OMP private(ji,jj), shared(self, grid, init_global_data, dx, dy)
        do jj = grid%subdomain%internal%ystart, grid%subdomain%internal%ystop
            do ji = grid%subdomain%internal%xstart, grid%subdomain%internal%xstop
                self%data(ji, jj) = init_global_data(ji+dx, jj+dy)
            end do
        end do
    end if
  end function r2d_field_constructor
   
  !===================================================
  ! Frees the data allocated for this array. Does nothing
  ! if the data was never allocated.
  subroutine r2d_free_field(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld
    ! Arguments
    !> fld: Pointer to the fld which data is to be freed
    if (allocated(fld%data)) then
       deallocate(fld%data)
    end if
  end subroutine r2d_free_field

  !===================================================

  subroutine read_from_device(self, startx, starty, nx, ny, blocking)
    !> Update the host data with the device data.
    ! The optional startx, starty, nx and ny parameters can be provided to
    ! read just from a sub-region of the field.
    ! The blocking optional parameter specifies if the data transfer should
    ! have finished on the return from this routine.

    class(r2d_field), target :: self
    integer, optional, intent(in) :: startx, starty, nx, ny
    logical, optional, intent(in) :: blocking
    integer :: local_startx, local_starty, local_nx, local_ny
    logical :: local_blocking

    if(self%data_on_device)then

       if (present(startx)) then
           local_startx = startx
       else
           local_startx = 1
       endif

       if (present(starty)) then
           local_starty = starty
       else
           local_starty = 1
       endif

       if (present(nx)) then
           local_nx = nx
       else
           local_nx = size(self%data, 1)
       endif

       if (present(ny)) then
           local_ny = ny
       else
           local_ny = size(self%data, 2)
       endif

       if (present(blocking)) then
           local_blocking = blocking
       else
           local_blocking = .true.
       endif

       if(associated(self%read_from_device_c))then
          call self%read_from_device_c(self%device_ptr, C_LOC(self%data), &
                int(local_startx, c_int), int(local_starty, c_int), &
                int(local_nx, c_int), int(local_ny, c_int), &
                logical(local_blocking, c_bool))
        else if(associated(self%read_from_device_f))then
          call self%read_from_device_f(self%device_ptr, self%data, &
                local_startx, local_starty, local_nx, local_ny, local_blocking)
        else
          call gocean_stop("ERROR: Data is on a device but no instructions " // &
              "about how to retrieve the data have been provided.")
       endif
    end if
  end subroutine read_from_device

  subroutine write_to_device(self, startx, starty, nx, ny, blocking)
    !> Update the device data with host data.
    ! The optional startx, starty, nx and ny parameters can be provided to
    ! write just a sub-region of the field to the device.
    ! The blocking optional parameter specifies if the data transfer should
    ! have finished on the return from this routine.

    class(r2d_field), target :: self
    integer, optional, intent(in) :: startx, starty, nx, ny
    logical, optional, intent(in) :: blocking
    integer :: local_startx, local_starty, local_nx, local_ny
    logical :: local_blocking

    if(self%data_on_device)then

       if (present(startx)) then
           local_startx = startx
       else
           local_startx = 1
       endif

       if (present(starty)) then
           local_starty = starty
       else
           local_starty = 1
       endif

       if (present(nx)) then
           local_nx = nx
       else
           local_nx = size(self%data, 1)
       endif

       if (present(ny)) then
           local_ny = ny
       else
           local_ny = size(self%data, 2)
       endif

       if (present(blocking)) then
           local_blocking = blocking
       else
           local_blocking = .true.
       endif

       if(associated(self%write_to_device_c))then
          call self%write_to_device_c(C_LOC(self%data), self%device_ptr, &
                int(local_startx, c_int), int(local_starty, c_int), &
                int(local_nx, c_int), int(local_ny, c_int), &
                logical(local_blocking, c_bool))
       else if(associated(self%write_to_device_f))then
          call self%write_to_device_f(self%data, self%device_ptr, &
                local_startx, local_starty, local_nx, local_ny, local_blocking)
        else
          call gocean_stop("ERROR: Data is on a device but no instructions " // &
              "about how to write new data have been provided.")
       endif
    endif
  end subroutine write_to_device


  !===================================================

  function get_data(self) result(dptr)
    !> Getter for the data associated with a field. If the data is on a
    ! device it ensures that the host copy is up-to-date with that on
    ! any accelerator device.
    class(r2d_field), target :: self
    real(go_wp), dimension(:,:), pointer :: dptr

    ! Synchronise the whole array
    call self%read_from_device()

    ! Return a reference to the data
    dptr => self%data
  end function get_data

  !===================================================

  function set_data(self, array) result(flag)
    !> Setter for the data associated with a field. If data is on a
    !! remote accelerator device then the device copy is updated too.
    implicit none
    class(r2d_field) :: self
    integer :: flag
    real(go_wp), dimension(:,:) :: array
    self%data = array

    ! Synchronise the whole array
    call self%write_to_device()

    flag = 0
  end function set_data

  !===================================================

  subroutine set_field_bounds(fld, fld_type, grid_points)
    implicit none
    !> The field we're working on. We use class as we want this
    ! to be polymorphic as this routine doesn't care whether the
    ! field is tiled or not.
    class(field_type), intent(inout) :: fld
    !> Which grid-point type the field is defined on
    integer,          intent(in)    :: grid_points
    !> Character string describing the grid-point type
    character(len=8), intent(out)   :: fld_type

    select case(grid_points)

    case(GO_U_POINTS)
       write(fld_type, "('C-U')")
       call cu_field_init(fld)
    case(GO_V_POINTS)
       write(fld_type, "('C-V')")
       call cv_field_init(fld)
    case(GO_T_POINTS)
       write(fld_type, "('C-T')")
       call ct_field_init(fld)
    case(GO_F_POINTS)
       write(fld_type, "('C-F')")
       call cf_field_init(fld)
    case(GO_ALL_POINTS)
       write(fld_type, "('C-All')")
       call field_init(fld)
    case default
       call gocean_stop('r2d_field_constructor: ERROR: invalid '//&
                        'specifier for type of mesh points')
    end select

    ! Compute and store dimensions of internal region of field
    fld%internal%nx = fld%internal%xstop - fld%internal%xstart + 1
    fld%internal%ny = fld%internal%ystop - fld%internal%ystart + 1

    ! In addition to the 'internal region' of the field, we may have
    ! external points that define B.C.'s or that act as halos. Here
    ! we store the full extent of the field, inclusive of such
    ! points.
    !> \todo Replace the use of NBOUNDARY here with info. computed
    !! from the T-point mask.
    if(fld%grid%boundary_conditions(1) /= GO_BC_PERIODIC)then
       fld%whole%xstart = fld%internal%xstart - NBOUNDARY
       fld%whole%xstop  = fld%internal%xstop  + NBOUNDARY
    else
       fld%whole%xstart = fld%internal%xstart - NBOUNDARY
       fld%whole%xstop  = fld%internal%xstop  + NBOUNDARY
    end if
    if(fld%grid%boundary_conditions(2) /= go_BC_PERIODIC)then
       fld%whole%ystart = fld%internal%ystart - NBOUNDARY
       fld%whole%ystop  = fld%internal%ystop  + NBOUNDARY
    else
       fld%whole%ystart = fld%internal%ystart - NBOUNDARY
       fld%whole%ystop  = fld%internal%ystop  + NBOUNDARY
    end if

    fld%whole%nx = fld%whole%xstop - fld%whole%xstart + 1
    fld%whole%ny = fld%whole%ystop - fld%whole%ystart + 1

  end subroutine set_field_bounds

  !===================================================

  subroutine field_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    fld%defined_on = GO_ALL_POINTS

    M = fld%grid%nx
    N = fld%grid%ny

    ! An 'all points' field is defined upon every point in the grid
    fld%internal%xstart = 1
    fld%internal%xstop  = M
    fld%internal%ystart = 1
    fld%internal%ystop  = N

    ! We have no halo regions
    fld%num_halos = 0

  end subroutine field_init
  
  !===================================================

  subroutine cu_field_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    fld%defined_on = GO_U_POINTS

    select case(fld%grid%offset)

    case(GO_OFFSET_SW)
       call cu_sw_init(fld)

    case(GO_OFFSET_NE)
       call cu_ne_init(fld)

    case default
       call gocean_stop('cu_field_init: ERROR - unsupported grid offset!')

    end select

  end subroutine cu_field_init

  !===================================================

  subroutine cu_sw_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    ! Set up a field defined on U points when the grid has
    ! a South-West offset:

    !   vi-1j+1--fij+1---vij+1---fi+1j+1
    !   |        |       |       |
    !   |        |       |       |
    !   Ti-1j----uij-----Tij-----ui+1j
    !   |        |       |       |
    !   |        |       |       |
    !   vi-1j----fij-----vij-----fi+1j
    !   |        |       |       |
    !   |        |       |       |
    !   Ti-1j-1--uij-1---Tij-1---ui+1j-1

    !
    if(fld%grid%boundary_conditions(1) == GO_BC_PERIODIC)then
       ! When implementing periodic boundary conditions, all mesh
       ! point types have the same extents as the grid of T points. We
       ! then have a halo of width grid_mod::HALO_WIDTH_X on either
       ! side of the domain.
       ! When updating a quantity on U points we write to:
       ! (using 'x' to indicate a location that is written and 'a' an
       ! additional point that shallow dispenses with):
       !
       ! i=istart i=istop
       !  o  o  o  o  a
       !  o  x  x  x  a  j=ystop
       !  o  x  x  x  a
       !  o  x  x  x  a  j=ystart
       !  a  a  a  a  a

       fld%internal%xstart = fld%grid%subdomain%internal%xstart
       fld%internal%xstop  = fld%grid%subdomain%internal%xstop
    else
       ! When updating a quantity on U points with this offset convention
       ! we write to (using 'x' to indicate a location that is written,
       !                    'b' a boundary point and
       !                    'o' a point that is external to the domain):
       !
       ! i=1         i=M
       !  o  b  b  b  b   j=N 
       !  o  b  x  x  b
       !  o  b  x  x  b
       !  o  b  x  x  b
       !  o  b  b  b  b   j=1
       fld%internal%xstart = fld%grid%subdomain%internal%xstart + 1
       fld%internal%xstop  = fld%grid%subdomain%internal%xstop
    end if

    fld%internal%ystart = fld%grid%subdomain%internal%ystart
    fld%internal%ystop  = fld%grid%subdomain%internal%ystop

    ! When applying periodic (wrap-around) boundary conditions (PBCs)
    ! we must fill the regions marked with 'b' above.
    ! This looks like (using x to indicate a location that is written
    ! first and y a location that is written second):
    !
    !  i=2      i=M
    ! _ y  y  y  y   j=N  
    ! /|x  o  o  o
    ! | x  o  o  o
    ! \ x  o  o  o
    !  \x_ o  o  o   j=1
    !   |\______/ 

    ! In array notation this looks like:
    !
    ! (2    , 1:N-1) = (M  , 1:N-1)
    ! (2:M, N) = (2:M, 1)

    call init_periodic_bc_halos(fld)

  end subroutine cu_sw_init

  !===================================================

  subroutine cu_ne_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    ! Set up a field defined on U points when the grid types have 
    ! a North-East offset relative to the T point.

    if(fld%grid%boundary_conditions(1) /= GO_BC_PERIODIC)then
       ! If we do not have periodic boundary conditions then we do
       ! not need to allow for boundary points here - they are
       ! already contained within the region defined by T mask.
       ! The T mask has been used to determine the grid%subdomain
       ! which describes the area on the grid that is actually being
       ! modelled (as opposed to having values supplied from B.C.'s etc.)
       fld%internal%xstart = fld%grid%subdomain%internal%xstart
       ! Now that all fields are allocated to be the same size, we always have
       ! as many U-points as T-points
       fld%internal%xstop  = fld%grid%subdomain%internal%xstop
    else
      call gocean_stop('ERROR: cu_ne_init: implement periodic boundary conditions!')
    end if
    if(fld%grid%boundary_conditions(2) /= GO_BC_PERIODIC)then
      fld%internal%ystart = fld%grid%subdomain%internal%ystart
      fld%internal%ystop  = fld%grid%subdomain%internal%ystop
    else
      call gocean_stop('ERROR: cu_ne_init: implement periodic BCs!')
    end if

!> \todo Is this concept of halo definitions useful?
    fld%num_halos = 0

  end subroutine cu_ne_init

  !===================================================

  subroutine cv_field_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    fld%defined_on = GO_V_POINTS

    select case(fld%grid%offset)

    case(GO_OFFSET_SW)
       call cv_sw_init(fld)

    case(GO_OFFSET_NE)
       call cv_ne_init(fld)

    case default
       call gocean_stop('cv_field_init: ERROR - unsupported grid offset!')

    end select

  end subroutine cv_field_init

  !===================================================

  subroutine cv_sw_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    if(fld%grid%boundary_conditions(2) == GO_BC_PERIODIC)then
       ! When implementing periodic boundary conditions, all
       ! mesh point types have the same extents as the grid of
       ! T points. We then have a halo of width 1 on either side
       ! of the domain.
       fld%internal%xstart = fld%grid%subdomain%internal%xstart
       fld%internal%xstop  = fld%grid%subdomain%internal%xstop

       fld%internal%ystart = fld%grid%subdomain%internal%ystart
       fld%internal%ystop  = fld%grid%subdomain%internal%ystop
    else
       ! When updating a quantity on V points we write to:
       ! (using x to indicate a location that is written):
       !
       ! i=1      i=M
       !  b  b  b  b   j=N 
       !  b  x  x  b
       !  b  x  x  b
       !  b  b  b  b
       !  o  o  o  o   j=1
       ! We're not offset from the T points in the x dimension so
       ! we have the same x bounds.
       fld%internal%xstart = fld%grid%subdomain%internal%xstart
       fld%internal%xstop  = fld%grid%subdomain%internal%xstop

       fld%internal%ystart = fld%grid%subdomain%internal%ystart + 1
       fld%internal%ystop  = fld%grid%subdomain%internal%ystop
       call gocean_stop('cv_sw_init: IMPLEMENT non-periodic BCs!')
    endif

    ! When applying periodic (wrap-around) boundary conditions (PBCs)
    ! we must fill the regions marked with 'b' above.
    ! This looks like (using x to indicate a location that is written
    ! first and y a location that is written second):
    !
    !  i=1     i=M
    !  -o  o  o  y   j=N  
    ! / o  o  o  y
    ! | o  o  o  y
    ! \ o  o  o  y
    !  \x  x  x _y   j=2
    !    \______/|

    ! In array notation this looks like:
    ! First row = last internal row
    ! field(1:M    ,1:1  ) = field(1:M,N-1:N-1)
    ! Last col = first internal col
    ! field(M:M,1:N) = field(2:2,  1:N)

    call init_periodic_bc_halos(fld)

  end subroutine cv_sw_init

  !===================================================

  subroutine cv_ne_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    if(fld%grid%boundary_conditions(1) /= GO_BC_PERIODIC)then
      ! If we do not have periodic boundary conditions then we do
      ! not need to allow for boundary points here - they are
      ! already contained within the region.
      fld%internal%xstart = fld%grid%subdomain%internal%xstart
      fld%internal%xstop  = fld%grid%subdomain%internal%xstop
    else
      call gocean_stop('ERROR: cv_ne_init: implement periodic BCs!')
    end if

    if(fld%grid%boundary_conditions(2) /= GO_BC_PERIODIC)then
       fld%internal%ystart = fld%grid%subdomain%internal%ystart
       ! Since we allocate all fields with the same extents, we always have
       ! as many V points as T points
       fld%internal%ystop  = fld%grid%subdomain%internal%ystop
    else
      call gocean_stop('ERROR: cv_ne_init: implement periodic BCs!')
    end if

  end subroutine cv_ne_init

  !===================================================

  subroutine ct_field_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    fld%defined_on = GO_T_POINTS

    select case(fld%grid%offset)

    case(GO_OFFSET_SW)
       call ct_sw_init(fld)

    case(GO_OFFSET_NE)
       call ct_ne_init(fld)

    case default
       call gocean_stop('ct_field_init: ERROR - unsupported grid offset!')

    end select

  end subroutine ct_field_init

  !===================================================

  subroutine ct_sw_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    ! When updating a quantity on T points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1      i=M
    !  b  b  b  b   j=N 
    !  b  x  x  b
    !  b  x  x  b
    !  b  b  b  b   j=1

    fld%internal%xstart = fld%grid%subdomain%internal%xstart
    fld%internal%xstop  = fld%grid%subdomain%internal%xstop
    fld%internal%ystart = fld%grid%subdomain%internal%ystart
    fld%internal%ystop  = fld%grid%subdomain%internal%ystop

    ! When applying periodic (wrap-around) boundary conditions
    ! (PBCs) we must fill the regions marked with 'b' above.
    ! This looks like (using x to indicate a location that is 
    ! written first and y a location that is written second):
    !
    !  i=1      i=M
    ! _ y  y  y  y   j=N  
    ! /|o  o  o  x
    ! | o  o  o  x
    ! \ o  o  o  x
    !  \o  o  o _x   j=1
    !    \______/|

    ! In array notation this looks like:
    ! Last col = first col
    ! field(M:M,  1:N-1  ) = field(1:1  ,1:N-1)
    ! Last row = first row
    ! field(1:M,N:N) = field(1:M,1:1)

    call init_periodic_bc_halos(fld)

  end subroutine ct_sw_init

  !===================================================

  subroutine ct_ne_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    ! When updating a quantity on T points with a NE offset
    ! we write to (using x to indicate a location that is written):
    !
    ! i=1          Nx
    !  b  b  b  b  b Ny
    !  b  x  x  x  b
    !  b  x  x  x  b 
    !  b  x  x  x  b
    !  b  b  b  b  b j=1

    if(fld%grid%boundary_conditions(1) /= GO_BC_PERIODIC)then
      ! If we do not have periodic boundary conditions then we do
      ! not need to allow for boundary points here - they are
      ! already contained within the region.
      ! Start and stop are just the same as those calculated from the T mask
      ! earlier because this is a field on T points.
      fld%internal%xstart = fld%grid%subdomain%internal%xstart
      fld%internal%xstop  = fld%grid%subdomain%internal%xstop
    else
      call gocean_stop('ERROR: ct_ne_init: implement periodic BCs!')
    end if

    if(fld%grid%boundary_conditions(2) /= GO_BC_PERIODIC)then
      ! Start and stop are just the same as those calculated from the T mask
      ! earlier because this is a field on T points.
      fld%internal%ystart = fld%grid%subdomain%internal%ystart
      fld%internal%ystop  = fld%grid%subdomain%internal%ystop
    else
      call gocean_stop('ERROR: ct_ne_init: implement periodic BCs!')
    end if

  end subroutine ct_ne_init

  !===================================================

  subroutine cf_field_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    fld%defined_on = GO_F_POINTS

    select case(fld%grid%offset)

    case(GO_OFFSET_SW)
       call cf_sw_init(fld)

    case(GO_OFFSET_NE)
       call cf_ne_init(fld)

    case default
       call gocean_stop('cf_field_init: ERROR - unsupported grid offset!')

    end select

  end subroutine cf_field_init

  !===================================================

  subroutine cf_sw_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    ! When updating a quantity on F points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1         i=M
    !  o  b  b  b  b   j=N 
    !  o  b  x  x  b
    !  o  b  x  x  b
    !  o  b  b  b  b
    !  o  o  o  o  o   j=1
    if(fld%grid%boundary_conditions(1) == GO_BC_PERIODIC)then
       fld%internal%xstart = fld%grid%subdomain%internal%xstart
       fld%internal%xstop  = fld%grid%subdomain%internal%xstop !internal%xstart + fld%grid%simulation_domain%nx - 1
    else
       fld%internal%xstart = fld%grid%subdomain%internal%xstart + 1
       fld%internal%xstop  = fld%grid%subdomain%internal%xstop
       ! I think these are correct but we stop because I've not properly
       ! gone through the coding.
       call gocean_stop('cf_sw_init: CHECK non-periodic BCs!')
    end if

    if(fld%grid%boundary_conditions(2) == GO_BC_PERIODIC)then
       fld%internal%ystart = fld%grid%subdomain%internal%ystart
       fld%internal%ystop  = fld%grid%subdomain%internal%ystop !fld%internal%ystart + fld%grid%simulation_domain%ny - 1
    else
       fld%internal%ystart = fld%grid%subdomain%internal%ystart + 1
       fld%internal%ystop  = fld%grid%subdomain%internal%ystop
       ! I think these are correct but we stop because I've not properly
       ! gone through the coding.
       call gocean_stop('cf_sw_init: CHECK non-periodic BCs!')
    end if


    ! When applying periodic (wrap-around) boundary conditions
    ! (PBCs) we must fill the regions marked with 'b' above.
    ! This looks like (using x to indicate a location that is 
    ! written first and y a location that is written second):
    !
    !  i=2      i=M
    ! .-x  o  o  o   j=N  
    ! | x  o  o  o
    ! | x  o  o  o
    ! | x  o  o  o
    ! ->y_ y  y  y   j=2
    !   |\______/ 

    ! In array notation this looks like:
    ! First col = last col
    ! field(2:2, 2:N) = field(M:M, 2:N)
    ! First row = last row
    ! field(2:M, 2:2) = field(2:M, N:N)

    call init_periodic_bc_halos(fld)

  end subroutine cf_sw_init

  !===================================================

  subroutine cf_ne_init(fld)
    implicit none
    class(field_type), intent(inout) :: fld

    ! When updating a quantity on F points we write to:
    ! (using x to indicate a location that is written
    !        b a boundary point - defined by ext. b.c.
    !        o a point that is external to the domain):
    !
    ! i=1       Nx
    !  o  o  o  o   Ny
    !  b  b  b  o   
    !  b  x  b  o   
    !  b  x  b  o
    !  b  b  b  o   j=1

    if(fld%grid%boundary_conditions(1) /= GO_BC_PERIODIC)then
      ! If we do not have periodic boundary conditions then we do
      ! not need to allow for boundary points here - they are
      ! already contained within the region.
      fld%internal%xstart = fld%grid%subdomain%internal%xstart
      fld%internal%xstop  = fld%grid%subdomain%internal%xstop
    else
      call gocean_stop('ERROR: cf_ne_init: implement periodic BCs!')
      stop
    end if

    if(fld%grid%boundary_conditions(2) /= GO_BC_PERIODIC)then
      fld%internal%ystart = fld%grid%subdomain%internal%ystart
      fld%internal%ystop  = fld%grid%subdomain%internal%ystop
    else
      call gocean_stop('ERROR: cf_ne_init: implement periodic BCs!')
    end if

  end subroutine cf_ne_init

  !===================================================

  SUBROUTINE copy_2dfield_array(field_in, field_out)
    IMPLICIT none
    REAL(go_wp), INTENT(in),  DIMENSION(:,:) :: field_in
    REAL(go_wp), INTENT(out), DIMENSION(:,:) :: field_out
        
    field_out(:,:) = field_in(:,:)
        
  end subroutine copy_2dfield_array

  !===================================================

  !> Copy from one patch in an array to another patch
  !! Will be superseded by interface using field object
  !! instead of array data.
  subroutine copy_2dfield_array_patch(field, src, dest)
    implicit none
    real(go_wp),       intent(inout), dimension(:,:) :: field
    type(region_type), intent(in)                    :: src, dest

    field(dest%xstart:dest%xstop,dest%ystart:dest%ystop) = &
     field(src%xstart:src%xstop ,src%ystart:src%ystop)
        
  end subroutine copy_2dfield_array_patch

  !===================================================

  SUBROUTINE copy_2dfield(field_in, field_out)
    IMPLICIT none
    type(r2d_field), intent(in)    :: field_in
    type(r2d_field), intent(inout) :: field_out
    integer :: it, ji, jj

    if (field_out%ntiles == 0) then
       field_out%data(:,:) = field_in%data(:,:)
    else
       !> \TODO #38 We don't actually support building dl_esm_inf
       !! with OpenMP enabled!
       !$OMP PARALLEL DO SCHEDULE(RUNTIME), DEFAULT(none), &
       !$OMP SHARED(field_out, field_in), PRIVATE(ji,jj)
       do it = 1, field_out%ntiles, 1
          do jj= field_out%tile(it)%whole%ystart, field_out%tile(it)%whole%ystop
             do ji = field_out%tile(it)%whole%xstart, field_out%tile(it)%whole%xstop
                field_out%data(ji,jj) = field_in%data(ji,jj)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine copy_2dfield

  !===================================================

  !> Copy from one patch to another in a field
  subroutine copy_2dfield_patch(field, src, dest)
    implicit none
    type(r2d_field), intent(inout) :: field
    type(region_type),    intent(in)    :: src, dest

    field%data(dest%xstart:dest%xstop,dest%ystart:dest%ystop) = &
     field%data(src%xstart:src%xstop ,src%ystart:src%ystop)
        
  end subroutine copy_2dfield_patch

  !===================================================

  SUBROUTINE set_field(fld, val)
    implicit none
    class(field_type), INTENT(out) :: fld
    real(go_wp), INTENT(in) :: val

    select type(fld)
    type is (r2d_field)
       fld%data = val
    class default
    end select

  END SUBROUTINE set_field

  !===================================================

  !> Compute the checksum of the internal values of the supplied field
  !> object. If the data is in a remote device, the field%get_data()
  !> call will first bring it to the host.
  function fld_checksum(field) result(val)
    implicit none
    type(r2d_field), intent(in) :: field
    real(go_wp) :: val

    val = array_checksum(field%get_data(),           &
                         field%internal%xstart, field%internal%xstop, &
                         field%internal%ystart, field%internal%ystop)
    return

  end function fld_checksum
  
  !===================================================

  !> Provide access to exchange_generic for halo swaps for this field.
  !! If the data is in a remote device, it will be synchronised with the host
  !! before and after the exchange.
  !! The depth value is currently ignored and hardwired to 1: the device
  !! synchonization and the halo_exchange communicators need to be updated
  !! with the provided depth.
  !!
  !! TODO #58: This method could benefit from an asynchronous execution model.
  subroutine halo_exchange(self, depth)
    use parallel_comms_mod, only: Iplus, Iminus, Jplus, Jminus, NONE, &
         exchange_generic
    implicit none
    class(r2d_field), target, intent(inout) :: self
    integer, intent(in) :: depth
    ! Locals
    integer :: exch  !> Handle for exchange

    ! Make sure the data from each halo is in the host before the communication
    call self%read_halo_from_device(comm=JPlus, blocking=.false.)
    call self%read_halo_from_device(comm=Jminus, blocking=.false.)
    call self%read_halo_from_device(comm=IPlus, blocking=.false.)
    call self%read_halo_from_device(comm=IMinus, blocking=.true.)

    ! Execute the halo exchange
    call exchange_generic(b2=self%data, handle=exch, &
                          comm1=JPlus, comm2=Jminus, comm3=IPlus, comm4=IMinus)

    ! Synchronize the data back to the accelerator device if necessary
    call self%write_halo_to_device(comm=JPlus, blocking=.false.)
    call self%write_halo_to_device(comm=Jminus, blocking=.false.)
    call self%write_halo_to_device(comm=IPlus, blocking=.false.)
    call self%write_halo_to_device(comm=IMinus, blocking=.true.)

  end subroutine halo_exchange
  
  subroutine read_halo_from_device(self, comm, blocking)
    use parallel_comms_mod, only: isrcsend, jsrcsend, nxsend, nysend
    implicit none
    class(r2d_field), target, intent(inout) :: self
    integer, intent(in) :: comm
    logical, intent(in) :: blocking

    ! Get location and dimension to read data for the given communicator
    if (isrcsend(comm) > 0) then
        call self%read_from_device(isrcsend(comm), jsrcsend(comm), &
            nxsend(comm), nysend(comm), blocking)
    endif
  end subroutine read_halo_from_device

  subroutine write_halo_to_device(self, comm, blocking)
    use parallel_comms_mod, only: idesrecv, jdesrecv, nxrecv, nyrecv
    implicit none
    class(r2d_field), target, intent(inout) :: self
    integer, intent(in) :: comm
    logical, intent(in) :: blocking

    ! Get location and dimension to write data for the given communicator
    if (idesrecv(comm) > 0) then
        call self%write_to_device(idesrecv(comm), jdesrecv(comm), &
            nxrecv(comm), nyrecv(comm), blocking)
    endif
  end subroutine write_halo_to_device
  !===================================================

  !> Compute the checksum of ALL of the elements of supplied array. Performs a
  !! global sum if built with distributed-memory support.
  function array_checksum(field, &
                          xstart, xstop, &
                          ystart, ystop) result(val)
    use parallel_comms_mod, only: global_sum
    implicit none
    real(go_wp), dimension(:,:), intent(in) :: field
    integer, optional, intent(in) :: xstart, xstop, ystart, ystop
    real(go_wp) :: val

    if( present(xstart) )then
       val = SUM( ABS(field(xstart:xstop,ystart:ystop) ) )
    else
       val = SUM( ABS(field(:,:)) )
    end if

    ! Sum over all processors if running with distributed memory
    call global_sum(val)

  end function array_checksum

  !===================================================

  !> Collect a (distributed) field into a global array
  !! on the master node.
  subroutine gather_inner_data(self, global_data)
    use parallel_utils_mod, only: get_num_ranks, gather
    use parallel_mod, only: on_master

    class(r2d_field), intent(in) :: self
    real(go_wp), dimension(:,:),           &
        allocatable, intent(out) :: global_data

    real(go_wp), dimension(:), allocatable :: send_buffer, recv_buffer

    integer :: dx, dy, ji, jj, i, n, rank, halo_x, halo_y
    integer :: x_start, x_stop, y_start, y_stop, ierr

    allocate(global_data(self%grid%global_nx, self%grid%global_ny), &
             stat=ierr)
    if(ierr /= 0)then
       call gocean_stop('gather_inner_data failed to allocate global result array')
    end if

    ! No MPI (or single process), just copy the data out
    if (get_num_ranks() == 1) then
        ! Compute size of inner area only
        dx = self%internal%xstart - 1
        dy = self%internal%ystart - 1
        do jj= self%internal%ystart, self%internal%ystop
           do ji = self%internal%xstart, self%internal%xstop
              global_data(ji-dx,jj-dy) = self%data(ji,jj)
           end do
        end do
        return
    endif

    ! Determine maximum size of data to be sent. We don't need
    ! to sent the halo, so reduce max_width and max_height by
    ! 2*halo size.
    halo_x = self%internal%xstart-1
    halo_y = self%internal%ystart-1
    n = (self%grid%decomp%max_width - 2*halo_x) *  &
        (self%grid%decomp%max_height - 2*halo_y)
    allocate(send_buffer(n), stat=ierr)
    if(ierr /= 0)then
       call gocean_stop('gather_inner_data failed to allocate send buffer')
    end if
    allocate(recv_buffer(n*get_num_ranks()), stat=ierr)
    if(ierr /= 0)then
       call gocean_stop('gather_inner_data failed to allocate receive buffer')
    end if

    ! Copy data into 1D send buffer.
    i = 0
    do jj= self%internal%ystart, self%internal%ystop
        do ji = self%internal%xstart, self%internal%xstop
            i = i + 1
            send_buffer(i) = self%data(ji,jj)
        end do
    end do

    ! Collect all send_buffers on the master:
    call gather(send_buffer, recv_buffer)

    if (on_master()) then
        ! Copy the data from each process into the global array
        do rank=1, get_num_ranks()
            x_start = self%grid%decomp%subdomains(rank)%global%xstart
            x_stop = self%grid%decomp%subdomains(rank)%global%xstop
            y_start = self%grid%decomp%subdomains(rank)%global%ystart
            y_stop = self%grid%decomp%subdomains(rank)%global%ystop
            i = (rank-1) * n
            do jj= y_start, y_stop
                do ji = x_start, x_stop
                    i = i + 1
                    global_data(ji, jj) = recv_buffer(i)
                end do
            end do
        enddo
    endif

  end subroutine gather_inner_data

  !===================================================

  subroutine init_periodic_bc_halos(fld)
    implicit none
    class(field_type), intent(inout) :: fld
    ! Locals
    integer :: ihalo

    ! Check whether we have PBCs in x AND y dimensions
    fld%num_halos = 0
    if( fld%grid%boundary_conditions(1) == GO_BC_PERIODIC )then
       fld%num_halos = fld%num_halos + 2
    end if
    if( fld%grid%boundary_conditions(2) == GO_BC_PERIODIC )then
       fld%num_halos = fld%num_halos + 2
    end if

    allocate( fld%halo(fld%num_halos) )

    ihalo = 0
    if( fld%grid%boundary_conditions(1) == GO_BC_PERIODIC )then
       ! E-most column set to W-most internal column
       ihalo = ihalo + 1
       fld%halo(ihalo)%dest%xstart = fld%internal%xstop + 1
       fld%halo(ihalo)%dest%xstop  = fld%internal%xstop + 1
       fld%halo(ihalo)%dest%ystart = fld%internal%ystart   
       fld%halo(ihalo)%dest%ystop  = fld%internal%ystop

       fld%halo(ihalo)%source%xstart = fld%internal%xstart   
       fld%halo(ihalo)%source%xstop  = fld%internal%xstart
       fld%halo(ihalo)%source%ystart = fld%internal%ystart  
       fld%halo(ihalo)%source%ystop  = fld%internal%ystop

       ! W-most column set to E-most internal column
       ihalo = ihalo + 1
       fld%halo(ihalo)%dest%xstart = fld%internal%xstart-1   
       fld%halo(ihalo)%dest%xstop  = fld%internal%xstart-1
       fld%halo(ihalo)%dest%ystart = fld%internal%ystart   
       fld%halo(ihalo)%dest%ystop  = fld%internal%ystop

       fld%halo(ihalo)%source%xstart = fld%internal%xstop 
       fld%halo(ihalo)%source%xstop  = fld%internal%xstop
       fld%halo(ihalo)%source%ystart = fld%internal%ystart
       fld%halo(ihalo)%source%ystop  = fld%internal%ystop
    end if

    if( fld%grid%boundary_conditions(2) == GO_BC_PERIODIC )then
       ! N-most row set to S-most internal row
       ihalo = ihalo + 1
       fld%halo(ihalo)%dest%xstart = fld%internal%xstart - 1   
       fld%halo(ihalo)%dest%xstop  = fld%internal%xstop  + 1
       fld%halo(ihalo)%dest%ystart = fld%internal%ystop+1   
       fld%halo(ihalo)%dest%ystop  = fld%internal%ystop+1

       fld%halo(ihalo)%source%xstart = fld%internal%xstart - 1
       fld%halo(ihalo)%source%xstop  = fld%internal%xstop  + 1
       fld%halo(ihalo)%source%ystart = fld%internal%ystart 
       fld%halo(ihalo)%source%ystop  = fld%internal%ystart

       ! S-most row set to N-most internal row
       ihalo = ihalo + 1
       fld%halo(ihalo)%dest%xstart = fld%internal%xstart - 1   
       fld%halo(ihalo)%dest%xstop  = fld%internal%xstop  + 1
       fld%halo(ihalo)%dest%ystart = fld%internal%ystart - 1   
       fld%halo(ihalo)%dest%ystop  = fld%internal%ystart - 1

       fld%halo(ihalo)%source%xstart = fld%internal%xstart - 1 
       fld%halo(ihalo)%source%xstop  = fld%internal%xstop  + 1
       fld%halo(ihalo)%source%ystart = fld%internal%ystop 
       fld%halo(ihalo)%source%ystop  = fld%internal%ystop
    end if

  end subroutine init_periodic_bc_halos

  !==============================================

  !> Routine to get the dimensions of the OpenMP tiling grid.
  !! Reads the GOCEAN_OMP_GRID environment variable which
  !! should have the format "NxM" where N is nx and M is ny.
  !! Returns false if the environment variable is not set
  !! or does not conform to this format.
  function get_grid_dims(nx, ny) result(success)
    implicit none
    integer, intent(inout) :: nx, ny
    logical :: success
    character(len=20) :: lstr
    integer :: idx, ierr

    success = .FALSE.

    call get_environment_variable(NAME='GOCEAN_OMP_GRID', VALUE=lstr, &
                                  STATUS=ierr)

    if(ierr /= 0)return

    ! We expect the string to have the format 'AxB' where A and B are
    ! integers.
    idx = index(lstr, 'x')
    if(idx == 0)then
       write (*,"(/'shallow_omp_mod::get_grid_dims: failed to parse ' &
                 &  'GOCEAN_OMP_GRID string: ',(A))") TRIM(lstr)
       write (*,"('   -  will use defaults for dimensions of tiling grid')")
       return
    endif

    read(lstr(1:idx-1),*,iostat=ierr) nx
    if(ierr /= 0)return

    read(lstr(idx+1:),*,iostat=ierr) ny
    if(ierr == 0)success = .TRUE.

  end function get_grid_dims

end module field_mod
