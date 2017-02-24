module grid_mod
  use kind_params_mod
  use region_mod
  use gocean_mod
  implicit none

  private

  ! Enumeration of possible grid types (we only actually
  ! support ARAKAWA_C at the moment)
  integer, public, parameter :: ARAKAWA_C = 0
  integer, public, parameter :: ARAKAWA_B = 1

  ! Enumeration of the four possible choices for 
  ! offsetting the grid-point types relative to the T point.
  !> Points to South and West of T point have same 
  !! i,j index (e.g. 'shallow' code)
  integer, public, parameter :: OFFSET_SW = 0
  integer, public, parameter :: OFFSET_SE = 1
  integer, public, parameter :: OFFSET_NW = 2
  !> Points to North and East of T point have same
  !! i,j index (e.g. NEMO code).
  integer, public, parameter :: OFFSET_NE = 3
  !> Value to signify no dependence on the relative offset of
  !! the different grid-point types.
  integer, public, parameter :: OFFSET_ANY = 4

  ! Enumeration of boundary-condition types
  !> Grid (model domain) has periodic boundary condition
  integer, public, parameter :: BC_PERIODIC = 0
  !> Grid (model domain) has external boundary conditions. This is a
  !! placeholder really as this is a complex area.
  integer, public, parameter :: BC_EXTERNAL = 1
  !> Grid (model domain) has no boundary conditions
  integer, public, parameter :: BC_NONE = 2

  !> The width of the halos we set up for implementing PBCs
  integer, parameter :: HALO_WIDTH_X = 1
  integer, parameter :: HALO_WIDTH_Y = 1

  type, public :: grid_type
     !> The type of grid this is (e.g. Arakawa C Grid)
     integer :: name
     !> Specifies the convention by which grid-point
     !! types are indexed relative to a T point.
     integer :: offset
     !> Total number of grid points
     integer :: npts
     !> Extent of T-point grid in x. Note that this is the whole grid,
     !! not just the region that is simulated.
     integer :: nx
     !> Extent of T-point grid in y. Note that this is the whole grid,
     !! not just the region that is simulated.
     integer :: ny
     !> Number of vertical levels that the grid has
     integer :: nlevels
     !> Uniform grid spacing in x (m)
     real(wp) :: dx
     !> Uniform grid spacing in y (m)
     real(wp) :: dy
     !> Uniform grid spacing in z (m)
     real(wp) :: dz

     !> Nature of each T point: 1.0 == wet inside simulated region
     !!                         0.0 == land
     !!                        -1.0 == wet outside simulated region
     !! This is the key quantity that determines the region that
     !! is actually simulated. However, we also support the
     !! specification of a model consisting entirely of wet points
     !! with periodic boundary conditions. Since this does not
     !! require a T-mask, we do not allocate this array for that
     !! case.
     real(wp), allocatable :: tmask(:,:,:)
     real(wp), allocatable :: umask(:,:,:)
     real(wp), allocatable :: vmask(:,:,:)

     !> The type of boundary conditions applied to the model domain
     !! in the x, y and z dimensions. Note that at this stage
     !! this is really only required for Periodic Boundary
     !! conditions.
     integer, dimension(3) :: boundary_conditions

     !> Where on the grid our simulated domain sits.
     !! \todo Decide whether this is useful.
     type(region_type) :: simulation_domain

     !> Horizontal scale factors at t point (m)
     real(wp), allocatable :: dx_t(:,:), dy_t(:,:)
     !> Horizontal scale factors at u point (m)
     real(wp), allocatable :: dx_u(:,:), dy_u(:,:)
     !> Horizontal scale factors at v point (m)
     real(wp), allocatable :: dx_v(:,:), dy_v(:,:)  
     !> Horizontal scale factors at f point (m)
     real(wp), allocatable :: dx_f(:,:), dy_f(:,:)
     !> Unknown \todo Name these fields!
     real(wp), allocatable :: area_t(:,:), area_u(:,:), area_v(:,:)
     !> Latitude of u points
     real(wp), allocatable :: gphiu(:,:)
     !> Latitude of v points
     real(wp), allocatable :: gphiv(:,:)
     !> Latitude of f points
     real(wp), allocatable :: gphif(:,:)

     !> Coordinates of grid (T) points in horizontal plane
     real(wp), allocatable :: xt(:,:), yt(:,:)
   contains
     procedure :: get_tmask

  end type grid_type

  interface grid_type
     module procedure grid_constructor
  end interface grid_type

  !> Single interface to the grid-initialisation routines
  !! for both 2D and 3D
  interface grid_init
     module procedure grid_init2d, grid_init3d
  end interface grid_init

  public grid_init

contains

  !============================================
  function get_tmask(self) result(tmask)
    implicit none
    class (grid_type), target, intent(in) :: self
    real(wp), pointer :: tmask(:,:,:)

    tmask => self%tmask

    return
  end function get_tmask

  !> Basic constructor for the grid type. Full details
  !! are fleshed-out by the grid_init() routine.
  function grid_constructor(grid_name, &
                            boundary_conditions, grid_offsets) result(self)
    implicit none
    integer, intent(in) :: grid_name
    !> The boundary conditions that will be applied to all fields that
    !! are defined on this grid.
    integer, dimension(3), intent(in) :: boundary_conditions
    !> The choice of the way the indices of the various grid-point types
    !! are offset from the T point. This is an optional argument as it
    !! will be supplied by PSyclone-modified call of this constructor.
    !! PSyclone will obtain the value of this offset from the kernel
    !! metadata. Despite being optional, it is actually required
    !! and we have a run-time check for this below.
    integer, optional, intent(in) :: grid_offsets
    type(grid_type), target :: self
    ! Locals
    integer :: grid_stagger

    ! We must be told the offset expected by the kernels. In a manual
    ! implementation it is up to the algorithm writer to inspect the
    ! kernels and provide us with the correct value. PSyclone will
    ! automate this process. This argument is only optional in order
    ! to ensure that the code to be processed by PSyclone will
    ! compile.
    if(present(grid_offsets))then
       grid_stagger = grid_offsets
    else
       call gocean_stop('ERROR: grid offset not specified in call to '//&
                        'grid_constructor.')
    end if

    ! This case statement is mainly to check that the caller
    ! has specified a valid value for grid_name.
    select case(grid_name)

    case(ARAKAWA_C)
       self%name = ARAKAWA_C
    case(ARAKAWA_B)
       self%name = ARAKAWA_B
    case default
       write(*,*) 'grid_constructor: ERROR: unsupported grid type: ', &
                  grid_name
       call gocean_stop('')
    end select

    ! Ditto for the choice of how the grid-point types are indexed
    ! relative to the T point
    select case(grid_stagger)

    case(OFFSET_NE)
       self%offset = OFFSET_NE
    case(OFFSET_NW)
       self%offset = OFFSET_NW
    case(OFFSET_SE)
       self%offset = OFFSET_SE
    case(OFFSET_SW)
       self%offset = OFFSET_SW
    case default
       write(*,*) 'grid_constructor: ERROR: unsupported relative offsets of grid types: ', &
                  grid_stagger
       call gocean_stop('')
    end select

    ! Store the boundary conditions that the model domain is
    ! subject to.
    self%boundary_conditions(1:3) = boundary_conditions(1:3)

  end function grid_constructor

  !============================================

  !> Initialise the supplied grid object for a 2D model
  !! consisting of m x n points. Ultimately, this routine should be
  !! general purpose but it is not there yet.
  !! N.B. the definition of m and n (the grid extents) depends on
  !! the type of boundary conditions that the model is subject to.
  !! For periodic boundary conditions they specify the extent of the
  !! simulated region (since we don't require the user to specify 
  !! the halos required to *implement* the PBCs). However, when a
  !! T-mask is used to define the model domain, m and n give the
  !! extents of that mask/grid. This, of necessity, includes boundary
  !! points. Therefore, the actual simulated region has an extent
  !! which is less than m x n.
  !! @param[inout] grid The object to initialise
  !! @param[in] m Extent in x of domain for which we have information
  !! @param[in] n Extent in y of domain for which we have information
  !! @param[in] dxarg Grid spacing in x dimension
  !! @param[in] dyarg Grid spacing in y dimension
  !! @param[in] tmask Array holding the T-point mask which defines
  !!                  the contents of the domain. Need not be
  !!                  supplied if domain is all wet and has PBCs.
  subroutine grid_init2d(grid, m, n, dxarg, dyarg, tmask)
    implicit none
    type(grid_type), intent(inout) :: grid
    integer,         intent(in)    :: m, n
    real(wp),        intent(in)    :: dxarg, dyarg
    integer, dimension(m,n), intent(in), optional :: tmask
    ! Locals
    integer :: mlocal
    integer :: ierr(5)
    integer :: ji, jj
    integer :: xstart, ystart ! Start of internal region of T-pts
    integer :: xstop, ystop ! End of internal region of T-pts

    ! Set-up the global x-y dimensions of the grid. If no T-mask has
    ! been supplied then we assume that we have periodic boundary conditions
    call set_2d_global_dims(grid, m, n, pbcs=(.not. present(tmask)))

    ! A 2D model has just a single level
    grid%nlevels = 1

    ! Copy-in the externally-supplied T-mask, if any. If using OpenMP
    ! then apply first-touch policy for data locality.
    if( present(tmask) )then
       allocate(grid%tmask(grid%nx,grid%ny,1), stat=ierr(1))
       if( ierr(1) /= 0 )then
          call gocean_stop('grid_init: failed to allocate array for T mask')
       end if
!$OMP PARALLEL DO schedule(runtime), default(none), private(ji,jj), &
!$OMP shared(m, n, grid, tmask)
       do jj = 1, n
          grid%tmask(1:m,jj,1) = REAL(tmask(1:m,jj), wp)
          ! Our saved mask is padded for alignment purposes so set
          ! any additional points to be outside the domain
          do ji = m+1, grid%nx
             grid%tmask(ji,jj,1) = REAL(tmask(m,jj), wp)
          end do
       end do
!$OMP END PARALLEL DO
       ! Additional rows
       do jj = n+1, grid%ny
          do ji = 1, m
             grid%tmask(ji,jj,1) = REAL(tmask(ji, n), wp)
          end do
       end do
       ! Additional corner points
       do jj = n+1, grid%ny
          do ji = m+1, grid%nx
             grid%tmask(ji,jj,1) = REAL(tmask(m,n), wp)
          end do
       end do

       ! Compute masks that are derived from the T-point mask
       !> \todo Rather than carry around extra arrays, it might
       !! be better to just compute these quantities as
       !! required.
       call compute_derived_masks(grid)

    else
       ! No T-mask supplied. Check that grid has PBCs in both
       ! x and y dimensions otherwise we won't know what to do.
       if( .not. ( (grid%boundary_conditions(1) == BC_PERIODIC) .and. &
                   (grid%boundary_conditions(2) == BC_PERIODIC) ) )then
          call gocean_stop('grid_init: ERROR: No T-mask supplied and '// &
                           'grid does not have periodic boundary conditions!')
       end if
    end if ! T-mask supplied

    ! Use the T mask to determine the dimensions of the
    ! internal, simulated region of the grid.
    ! This call sets grid%simulation_domain.
    call compute_internal_region(grid, m, n)

    ! For a regular, orthogonal mesh the spatial resolution is constant
    call setup_xy_scalefactors(grid, dxarg, dyarg)

    ! Co-ordinates of the T points
    call set_tpt_coords(grid)

  end subroutine grid_init2d

  !============================================

  !> Initialise the supplied grid object for a 3D model
  !! consisting of m x n points and nlevel vertical levels. Ultimately,
  !! this routine should be general purpose but it is not there yet.
  !! N.B. the definition of m and n (the grid extents) depends on
  !! the type of boundary conditions that the model is subject to.
  !! For periodic boundary conditions they specify the extent of the
  !! simulated region (since we don't require the user to specify 
  !! the halos required to *implement* the PBCs). However, when a
  !! T-mask is used to define the model domain, m and n give the
  !! extents of that mask/grid. This, of necessity, includes boundary
  !! points. Therefore, the actual simulated region has an extent
  !! which is less than m x n.
  !! @param[inout] grid The object to initialise
  !! @param[in] m Extent in x of domain for which we have information
  !! @param[in] n Extent in y of domain for which we have information
  !! @param[in] nlevel No. of vertical levels in the domain
  !! @param[in] dxarg Grid spacing in x dimension
  !! @param[in] dyarg Grid spacing in y dimension
  !! @param[in] dzarg Grid spacing in z dimension
  !! @param[in] tmask Array holding the T-point mask which defines
  !!                  the contents of the domain. Need not be
  !!                  supplied if domain is all wet and has PBCs.
  subroutine grid_init3d(grid, m, n, nlevel, dxarg, dyarg, dzarg, tmask)
    implicit none
    type(grid_type), intent(inout) :: grid
    integer,         intent(in)    :: m, n, nlevel
    real(wp),        intent(in)    :: dxarg, dyarg, dzarg
    real(wp), dimension(m,n,nlevel), intent(in), optional :: tmask
    ! Locals
    integer :: mlocal
    integer :: ierr(5)
    integer :: ji, jj, jk
    integer :: xstart, ystart ! Start of internal region of T-pts
    integer :: xstop, ystop ! End of internal region of T-pts

    ! Set-up the global x-y dimensions of the grid. If no T-mask has
    ! been supplied then we assume that we have periodic boundary conditions
    call set_2d_global_dims(grid, m, n, pbcs=(.not. present(tmask)))

    ! Store the global vertical dimension of the grid.
    grid%nlevels = nlevel

    ! Copy-in the externally-supplied T-mask, if any. If using OpenMP
    ! then apply first-touch policy for data locality.
    if( present(tmask) )then
       allocate(grid%tmask(grid%nx,grid%ny,grid%nlevels), stat=ierr(1))
       if( ierr(1) /= 0 )then
          call gocean_stop('grid_init: failed to allocate array for T mask')
       end if
!$OMP PARALLEL DO schedule(runtime), default(none), private(ji,jj,jk), &
!$OMP shared(m, n, nlevel, grid, tmask)
       do jk = 1, nlevel
          do jj = 1, n
             grid%tmask(1:m,jj,jk) = tmask(1:m,jj,jk)
             ! Our saved mask is padded for alignment purposes so set
             ! any additional points to be outside the domain
             do ji = m+1, grid%nx
                grid%tmask(ji,jj,jk) = tmask(m,jj,jk)
             end do
          end do
       end do
!$OMP END PARALLEL DO
       do jk = 1, nlevel
          ! Additional rows
          do jj = n+1, grid%ny
             do ji = 1, m
                grid%tmask(ji,jj,jk) = tmask(ji,n,jk)
             end do
          end do
          ! Additional corner points
          do jj = n+1, grid%ny
             do ji = m+1, grid%nx
                grid%tmask(ji,jj,jk) = tmask(m,n,jk)
             end do
          end do
       end do

       ! Compute masks that are derived from the T-point mask
       !> \todo Rather than carry around extra arrays, it might
       !! be better to just compute these quantities as
       !! required.
       call compute_derived_masks(grid)

    else
       ! No T-mask supplied. Check that grid has PBCs in both
       ! x and y dimensions otherwise we won't know what to do.
       if( .not. ( (grid%boundary_conditions(1) == BC_PERIODIC) .and. &
                   (grid%boundary_conditions(2) == BC_PERIODIC) ) )then
          call gocean_stop('grid_init: ERROR: No T-mask supplied and '// &
                           'grid does not have periodic boundary conditions!')
       end if
    end if ! T-mask supplied

    ! Use the T mask to determine the dimensions of the
    ! internal, simulated region of the grid.
    ! This call sets grid%simulation_domain.
    call compute_internal_region(grid, m, n)

    ! For a regular, orthogonal mesh the spatial resolution is constant
    call setup_xy_scalefactors(grid, dxarg, dyarg)
    grid%dz = dzarg

    ! Co-ordinates of the T points
    call set_tpt_coords(grid)

  end subroutine grid_init3d

  !================================================

  subroutine setup_xy_scalefactors(grid, dx, dy)
    !> Computes x/y-related scale factors for the grid
    implicit none
    type(grid_type), intent(inout) :: grid
    real(wp), intent(in) :: dx, dy
    ! Locals
    integer :: ji, jj, ierr(5)

    grid%dx = dx
    grid%dy = dy

    allocate(grid%dx_t(grid%nx,grid%ny), grid%dy_t(grid%nx,grid%ny), &
             grid%dx_u(grid%nx,grid%ny), grid%dy_u(grid%nx,grid%ny), &
             stat=ierr(1))
    allocate(grid%dx_f(grid%nx,grid%ny), grid%dy_f(grid%nx,grid%ny), &
             grid%dx_v(grid%nx,grid%ny), grid%dy_v(grid%nx,grid%ny), &
             stat=ierr(2)) 
    allocate(grid%area_t(grid%nx,grid%ny), grid%area_u(grid%nx,grid%ny), &
             grid%area_v(grid%nx,grid%ny), stat=ierr(3))
    allocate(grid%gphiu(grid%nx,grid%ny), grid%gphiv(grid%nx,grid%ny), &
             grid%gphif(grid%nx,grid%ny), stat=ierr(4))
    allocate(grid%xt(grid%nx,grid%ny), grid%yt(grid%nx,grid%ny), stat=ierr(5))

    if( any(ierr /= 0, 1) )then
       call gocean_stop('setup_xy_scalefactors: failed to allocate arrays')
    end if

    ! Initialise the horizontal scale factors for a regular,
    ! orthogonal mesh. (Constant spatial resolution.)
!$OMP PARALLEL DO schedule(runtime), default(none), private(ji,jj), &
!$OMP shared(grid)
    do jj = 1, grid%ny
       do ji = 1, grid%nx
          grid%dx_t(ji, jj)   = grid%dx
          grid%dy_t(ji, jj)   = grid%dy

          grid%dx_u(ji, jj)   = grid%dx
          grid%dy_u(ji, jj)   = grid%dy

          grid%dx_v(ji, jj)   = grid%dx
          grid%dy_v(ji, jj)   = grid%dy

          grid%dx_f(ji, jj)   = grid%dx
          grid%dy_f(ji, jj)   = grid%dy
       end do
    end do
!$OMP END PARALLEL DO

    ! calculate t,u,v cell area
!$OMP PARALLEL DO schedule(runtime), default(none), private(ji,jj), &
!$OMP shared(grid)
    do jj = 1, grid%ny
       do ji = 1, grid%nx
          grid%area_t(ji,jj) = grid%dx_t(ji,jj) * grid%dy_t(ji,jj)

          grid%area_u(ji,jj) = grid%dx_u(ji,jj) * grid%dy_u(ji,jj)

          grid%area_v(ji,jj) = grid%dx_v(ji,jj) * grid%dy_v(ji,jj)
       END DO
    END DO
!$OMP END PARALLEL DO

    ! -here is an f-plane testing case
    ! i.e. the Coriolis parameter is set to a constant value.
!$OMP PARALLEL DO schedule(runtime), default(none), private(ji,jj), &
!$OMP shared(grid)
    do jj = 1, grid%ny
       do ji = 1, grid%nx
          grid%gphiu(ji, jj) = 50._wp
          grid%gphiv(ji, jj) = 50._wp
          grid%gphif(ji, jj) = 50._wp
       end do
    end do
!$OMP END PARALLEL DO

  end subroutine setup_xy_scalefactors

  !================================================

  subroutine set_tpt_coords(grid)
    !> Set the locations of the T points in 2D
    implicit none
    type(grid_type), intent(inout) :: grid
    ! Locals
    integer :: ji, jj
    integer :: xstart, xstop, ystart, ystop

    ! Do first-touch initialisation before setting actual values
!$OMP PARALLEL DO schedule(runtime), default(none), private(ji,jj), &
!$OMP shared(grid)
    do jj = 1, grid%ny
       do ji = 1, grid%nx
          grid%xt(ji,jj) = 0.0
          grid%yt(ji,jj) = 0.0
       end do
    end do

    xstart = grid%simulation_domain%xstart
    xstop  = grid%simulation_domain%xstop
    ystart = grid%simulation_domain%ystart
    ystop  = grid%simulation_domain%ystop
    grid%xt(xstart, :) = 0.0_wp + 0.5_wp * grid%dx_t(xstart,:)
    grid%yt(:,ystart)  = 0.0_wp + 0.5_wp * grid%dy_t(:,ystart)

    DO ji = xstart+1, xstop
      grid%xt(ji,ystart:ystop) = grid%xt(ji-1, ystart:ystop) + grid%dx
    END DO
            
    DO jj = ystart+1, ystop
      grid%yt(xstart:xstop,jj) = grid%yt(xstart:xstop, jj-1) + grid%dy
    END DO
  end subroutine set_tpt_coords

  !================================================

  subroutine set_2d_global_dims(grid, m, n, pbcs)
    use global_parameters_mod, only: ALIGNMENT
    implicit none
    !> Store the global x-y dimensions of the grid (the nx and ny
    !! components of the grid object).
    type(grid_type), intent(inout) :: grid
    integer, intent(in) :: m, n
    !> Whether or not we apply periodic boundary conditions
    logical, intent(in) :: pbcs
    ! Locals
    integer :: mlocal, nlocal

    if( .not. pbcs )then
       ! A T-mask has been supplied and that tells us everything
       ! about the extent of this model.
       ! Extend the domain by unity in each dimension to allow
       ! for staggering of variables. All fields will be
       ! allocated with extent (nx,ny).
       mlocal = m + 1
       if( mod(mlocal, ALIGNMENT) > 0 )then
          ! Since this is the dimension of the array and not that of
          ! the internal region, we add two lots of 'ALIGNMENT'. This
          ! allows us to subsequently extend the loop over the internal
          ! region so that it too is aligned without array accesses of
          ! the form a(i+1,j) going out of bounds.
          grid%nx = (mlocal/ALIGNMENT + 2)*ALIGNMENT
       else
          grid%nx = mlocal
       end if
       grid%ny = n + 1
    else
       ! No T-mask has been supplied so we assume we're implementing
       ! periodic boundary conditions and allow for halos of width
       ! HALO_WIDTH_{X,Y} here.  Currently we put a halo on all four
       ! sides of our rectangular domain. This is actually unnecessary
       ! - depending on the variable staggering used only one of the
       ! E/W halos and one of the N/S halos are required. However,
       ! that is an optimisation and this framework must be developed
       ! in such a way that that optimisation is supported.
       mlocal = m + 2*HALO_WIDTH_X
       if( mod(mlocal, ALIGNMENT) > 0 )then
          ! Since this is the dimension of the array and not that of
          ! the internal region, we add two lots of 'ALIGNMENT'. This
          ! allows us to subsequently extend the loop over the internal
          ! region so that it too is aligned without array accesses of
          ! the form a(i+1,j) going out of bounds.
          grid%nx = (mlocal/ALIGNMENT + 2)*ALIGNMENT
       else
          grid%nx = mlocal
       end if

       grid%ny = n + 2*HALO_WIDTH_Y
    end if

  end subroutine set_2d_global_dims

  !================================================

  !> Use the T-mask to deduce the inner or simulated region
  !! of the supplied grid.
  subroutine compute_internal_region(grid, morig, norig)
    implicit none
    type(grid_type), intent(inout) :: grid
    !> The original dimensions of the supplied T-mask (before we
    !! padded for alignment). This is a bit of a hack to the 
    !! interface but it will do for now (in the absence of an 
    !! algorithm to determine the internal region).
    integer,         intent(in) :: morig, norig

    if( allocated(grid%tmask) )then

       ! Here we will loop over the grid points, looking for
       ! the first occurrence of wet points.
       ! However, for the moment we just hardwire the routine
       ! to return results appropriate for a T mask that has
       ! a shell of unit depth of boundary/external points:

       ! i= 1           nx 
       !    b   b   b   b   ny
       !    b   x   x   b  
       !    b   x   x   b  
       !    b   x   x   b   
       !    b   b   b   b   1
       !                    j

       ! The actual part of this domain that is simulated. The outer-most 
       ! rows and columns of T points are not in the domain but are needed
       ! to specify the type of boundary (whether hard or open).
       !> \todo Generate the bounds of the simulation domain by
       !! examining the T-point mask.
       ! This defines the internal region of any T-point field.
       grid%simulation_domain%xstart = 2
       grid%simulation_domain%xstop  = morig - 1
       grid%simulation_domain%ystart = 2
       grid%simulation_domain%ystop  = norig - 1

    else

       ! We don't have a T mask so we must have PBCs in both x and y
       ! dimensions. In this case, the grid dimensions stored in grid%{nx,ny}
       ! may have been padded for alignment and so we use morig,norig from
       ! the namelist file. These are taken to specify the dimension of the
       ! *simulated* domain, excluding halos.
       grid%simulation_domain%xstart = 2
       grid%simulation_domain%xstop  = grid%simulation_domain%xstart + &
                                         morig - 1
       grid%simulation_domain%ystart = 2
       grid%simulation_domain%ystop  = grid%simulation_domain%ystart + &
                                         norig - 1

    end if

    grid%simulation_domain%nx =  grid%simulation_domain%xstop -  &
                                 grid%simulation_domain%xstart + 1
    grid%simulation_domain%ny =  grid%simulation_domain%ystop -  &
                                 grid%simulation_domain%ystart + 1

    write(*,"('GRID: Simulation domain = [',I4,':',I4,',',I4,':',I4,']')") &
         grid%simulation_domain%xstart, &
         grid%simulation_domain%xstop,  &
         grid%simulation_domain%ystart, &
         grid%simulation_domain%ystop

  end subroutine compute_internal_region

  !================================================

  subroutine compute_derived_masks(grid)
    implicit none
    type(grid_type), intent(inout) :: grid
    ! Locals
    integer :: ji, jj, jk, ierr
    integer :: nx, ny, nz

    if(.not. allocated(grid%tmask))then
       call gocean_stop('compute_derived_masks: T-mask not yet allocated.')
    end if

    nx = grid%nx
    ny = grid%ny
    nz = grid%nlevels

    allocate(grid%umask(nx, ny, nz), grid%vmask(nx,ny,nz), Stat=ierr)
    if(ierr /= 0)then
       call gocean_stop('compute_derived_masks: failed to allocate masks')
    end if

    ! This code is taken from dommsk.F90 in NEMO
    DO jk = 1, nz
       DO jj = 1, ny-1
          DO ji = 1, nx-1
             grid%umask(ji,jj,jk) = grid%tmask(ji,jj,jk)*grid%tmask(ji+1,jj  ,jk)
             grid%vmask(ji,jj,jk) = grid%tmask(ji,jj,jk)*grid%tmask(ji  ,jj+1,jk)
          END DO
! We don't support an f-mask just yet
!          DO ji = 1, jpi-1
!             grid%fmask(ji,jj,jk) = grid%tmask(ji,jj,jk)*grid%tmask(ji+1,jj  ,jk)   &
!                  &            * grid%tmask(ji,jj+1,jk)*grid%tmask(ji+1,jj+1,jk)
!          END DO
       END DO
    END DO

  end subroutine compute_derived_masks

end module grid_mod
