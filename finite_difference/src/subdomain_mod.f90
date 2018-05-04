!> Module containing the definition of the sub-domain type.
module subdomain_mod
  use region_mod
  implicit none

  !> Type encapsulating the information required to define a single
  !! sub-domain
  type :: subdomain_type
     !> Position of bottom-left and top-right corners of this sub-domain
     !! in the global domain.
     integer :: xstart, ystart, xstop, ystop
     !> Full width and height of this subdomain - includes all halo and
     !! boundary points. n{x,y} != {x,y}stop - {x,y}start +1 because
     !! {x,y}start and {x,y}stop are the coordinates of the internal
     !! region in the global domain.
     integer :: nx, ny
     !> The internal region of this subdomain (excluding halo and
     !! boundary points)
     type(region_type) :: internal
  end type subdomain_type

  !> Type encapsulating all information regarding a regular, 2D
  !! domain decomposition
  type :: decomposition_type
     !> Dimensions of the grid of sub-domains
     integer :: nx, ny
     !> Number of sub-domains (=nx*ny)
     integer :: ndomains
     !> Max dimensions of any sub-domain
     integer :: max_width, max_height
     !> Array of the sub-domain definitions
     type(subdomain_type), allocatable :: subdomains(:)
  end type decomposition_type

contains

  function decompose(domainx, domainy,             &
                     ndomains, ndomainx, ndomainy, &
                     halo_width) result(decomp)
    !> Decompose a domain consisting of domainx x domainy points
    !! into a 2D grid.
    !! Returns an array of tiles describing each of the subdomains.    
    use parallel_mod, only: get_num_ranks, parallel_abort
    implicit none
    !> The decomposition that this function will return
    type(decomposition_type), target :: decomp
    !> Dimensions of the domain to decompose
    integer, intent(in) :: domainx, domainy
    !> No. of domains to decompose into. If not supplied will use the
    !! number of MPI ranks (or 1 if serial).
    integer, intent(in), optional :: ndomains
    !> Optional specification of the dimensions of the tiling grid
    integer, intent(in), optional :: ndomainx, ndomainy
    !> Width of halo to allow for when constructing sub-domains.
    !! Default value is 1. Must be > 0 for a multi-processor run.
    integer, intent(in), optional :: halo_width
    ! Locals
    integer :: ival, jval ! For tile extent calculation
    integer :: internal_width, internal_height
    integer :: ierr, nwidth
    integer :: ji,jj, ith
    integer :: jover, junder, height
    integer :: iover, iunder, width
    ! For doing stats on tile sizes
    integer :: nvects, nvects_sum, nvects_min, nvects_max 
    logical, parameter :: print_tiles = .TRUE.
    ! Whether to automatically compute the dimensions of the tiling grid
    logical :: auto_tile = .TRUE.
    integer :: xlen, ylen
    integer :: ntilex, ntiley
    ! Local var to hold number of sub-domains
    integer :: ndom = 1
    ! Local var to hold halo width
    integer :: hwidth = 1
    ! Pointer to the current subdomain being defined
    type(subdomain_type), pointer :: subdomain

    ! Deal with optional arguments to this routine
    if(.not. present(ndomains))then
       if(.not. present(ndomainx) .and. .not. present(ndomainy))then
          ! Automatically use the number of MPI processes
          ndom = get_num_ranks()
          auto_tile = .TRUE.
       else if(present(ndomainx) .and. present(ndomainy))then
          ndom = ndomainx * ndomainy
          auto_tile = .FALSE.
       else
          call parallel_abort('decompose: invalid arguments supplied')
       end if
    else
       ndom = ndomains
       auto_tile = .TRUE.
    end if

    if(present(halo_width))then
       if(halo_width < 1 .and. get_num_ranks() > 1)then
          call parallel_abort('decompose: halo width must be > 0 if '// &
                              'running on more than one process')
       end if
       hwidth = halo_width
    end if

    decomp%ndomains = ndom
    allocate(decomp%subdomains(ndom), Stat=ierr)
    if(ierr /= 0)then
       call parallel_abort('decompose: failed to allocate tiles array')
    end if

    xlen = domainx
    ylen = domainy
    
    if(auto_tile)then

       ntilex = INT( SQRT(REAL(ndom)) )
       do WHILE(MOD(ndom, ntilex) /= 0)
          ntilex = ntilex - 1
       end do
       ntiley = ndom/ntilex

       ! Match longest dimension of domain to longest dimension of 
       ! process grid
       if(xlen > ylen)then
          if( ntilex < ntiley )then
             ierr   = ntiley
             ntiley = ntilex
             ntilex = ierr
          end if
       else
          ! N >= M so want nprocy >= nprocx
          if( ntiley < ntilex )then
             ierr   = ntiley
             ntiley = ntilex
             ntilex = ierr
          end if
       end if
    else
       ntilex = ndomainx
       ntiley = ndomainy
    end if ! automatic determination of tiling grid

    WRITE (*,"('decompose: using grid of ',I3,'x',I3)") ntilex, ntiley
    decomp%nx = ntilex
    decomp%ny = ntiley

    ! If we think about the internal regions of the tiles,
    ! then they should share the domain between them:
    internal_width = NINT(REAL(xlen) / REAL(ntilex))
    internal_height = NINT(REAL(ylen) / REAL(ntiley))

    ! Integer arithmetic means that ntiley tiles of height idy might
    ! actually span a height greater or less than N. If so, we try and
    ! reduce the height of each row by just one until we've accounted
    ! for the <jover> extra rows.
    !nwidth = (ntiley-2)*(idy-2) + 2*(idy-1)
    nwidth = ntiley * internal_height
    if(nwidth > ylen)then
       jover  = nwidth - ylen
       junder = 0
    else if(nwidth < ylen)then
       jover  = 0
       junder = ylen - nwidth
    else
       jover  = 0
       junder = 0
    end if
    ! Ditto for x dimension
    !nwidth = (ntilex-2)*(idx-2) + 2*(idx-1)
    nwidth = ntilex * internal_width
    if(nwidth > xlen)then
       iover  = nwidth - xlen
       iunder = 0
    else if(nwidth < xlen)then
       iover  = 0
       iunder = xlen - nwidth
    else
       iover  = 0
       iunder = 0
    end if

    if(print_tiles)then
       WRITE(*,"('Tile width = ',I4,', tile height = ',I4)") &
            internal_width, internal_height
       WRITE(*,"('iover = ',I3,', iunder = ',I3)") iover, iunder
       WRITE(*,"('jover = ',I3,', junder = ',I3)") jover, junder
    end if

    ith = 1
    ! The starting point of the tiles in y
    jval = 1

    nvects_max = 0
    nvects_min = 1000000
    nvects_sum = 0
    decomp%max_width  = 0
    decomp%max_height = 0

    if(print_tiles)WRITE(*,"(/'Sub-domains:')")

    do jj = 1, ntiley, 1

       ! Point to subdomain we are about to define
       subdomain => decomp%subdomains(ith)

       ! If necessary, correct the height of this tile row
       if(jover > 0)then
          height = internal_height - 1
          jover = jover - 1
       else if(junder > 0)then
          height = internal_height + 1
          junder = junder - 1
       else
          height = internal_height
       end if

       !  . . . . . . . .
       !  o o o o o . . .
       !  o h h h h
       !  o h h h h
       !  o h h i i
       !  o h h i i
       !  o h h i i
       !  .
       !  .

       ! The starting point of the tiles in x
       ival = 1

       do ji = 1, ntilex, 1
         
          ! If necessary, correct the width of this tile column
          if(iover > 0)then
             width = internal_width - 1
             iover = iover - 1
          else if(iunder > 0)then
             width = internal_width + 1
             iunder = iunder - 1
          else
             width = internal_width
          end if

          subdomain => decomp%subdomains(ith)
          subdomain%internal%xstart = hwidth+1
          subdomain%internal%xstop = subdomain%internal%xstart+width-1
          subdomain%internal%nx = subdomain%internal%xstop - &
                                  subdomain%internal%xstart + 1
          ! Which part of the global domain this represents
          subdomain%xstart = ival
          subdomain%xstop = ival + subdomain%internal%nx - 1
          subdomain%nx  = 2*hwidth + subdomain%internal%nx

          subdomain%internal%ystart = hwidth+1
          subdomain%internal%ystop = subdomain%internal%ystart+height-1
          subdomain%internal%ny = subdomain%internal%ystop - &
                                  subdomain%internal%ystart + 1
          ! Which part of the global domain this represents
          subdomain%ystart = jval
          subdomain%ystop = jval + subdomain%internal%ny - 1
          subdomain%ny = 2*hwidth + subdomain%internal%ny

          if(print_tiles)then
             WRITE(*,"('subdomain[',I4,'](',I4,':',I4,')(',I4,':',I4,'), "// &
                  & "interior:(',I4,':',I4,')(',I4,':',I4,') ')")       &
                  ith,                                                  &
                  subdomain%xstart, subdomain%xstop,        &
                  subdomain%ystart, subdomain%ystop,        &
                  subdomain%internal%xstart, subdomain%internal%xstop, &
                  subdomain%internal%ystart, subdomain%internal%ystop
          end if

          ! Collect some data on the distribution of tile sizes for 
          ! loadbalance info
          nvects = subdomain%internal%nx * subdomain%internal%ny
          nvects_sum = nvects_sum + nvects
          nvects_min = MIN(nvects_min, nvects)
          nvects_max = MAX(nvects_max, nvects)

          ! For use when allocating tile-'private' work arrays
          decomp%max_width  = MAX(decomp%max_width, subdomain%nx)
          decomp%max_height = MAX(decomp%max_height, subdomain%ny)

          ival = subdomain%xstop + 1
          ith = ith + 1
       end do
       jval = subdomain%ystop + 1
    end do

    ! Print tile-size statistics
    if(print_tiles)then
       write(*,"(/'Mean sub-domain size = ',F8.1,' pts = ',F7.1,' KB')") &
                                     REAL(nvects_sum)/REAL(decomp%ndomains), &
                              REAL(8*nvects_sum)/REAL(decomp%ndomains*1024)
       write(*,"('Min,max sub-domain size (pts) = ',I6,',',I6)") &
                                                      nvects_min, nvects_max
       write(*,"('Domain load imbalance (%) =',F6.2)") &
                               100.0*(nvects_max-nvects_min)/REAL(nvects_min)
       write(*,"('Max sub-domain dims are ',I4,'x',I4/)") decomp%max_width, &
                                                          decomp%max_height
    end if

  end function decompose

end module subdomain_mod
