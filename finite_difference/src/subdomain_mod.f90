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

contains

  function decompose(domainx, domainy,   &
                     ndomains, ndomainx, &
                     ndomainy) result(subdomains)
    use parallel_mod, only: get_num_ranks, parallel_abort
    implicit none
    !! Decompose a domain consisting of domainx x domainy points
    !! into a 2D grid.
    !! Returns an array of tiles describing each of the subdomains.
    
    !> The array of tiles this function will return
    type(subdomain_type), allocatable :: subdomains(:)
    !> Dimensions of the domain to decompose
    integer, intent(in) :: domainx, domainy
    !> No. of domains to decompose into
    integer, intent(in), optional :: ndomains
    !> Optional specification of the dimensions of the tiling grid
    integer, intent(in), optional :: ndomainx, ndomainy
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
    ! TODO these need to be stored with the generated decomposition
    ! maybe we need a decomposition_type?
    integer :: max_tile_height
    integer :: max_tile_width
    integer :: halo_width = 2

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

    allocate(subdomains(ndom), Stat=ierr)
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
    max_tile_width  = 0
    max_tile_height = 0

    if(print_tiles)WRITE(*,"(/'Sub-domains:')")

    do jj = 1, ntiley, 1

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

          subdomains(ith)%internal%xstart = halo_width+1
          subdomains(ith)%internal%xstop = &
                                 subdomains(ith)%internal%xstart+width-1
          subdomains(ith)%internal%nx = subdomains(ith)%internal%xstop - &
                                        subdomains(ith)%internal%xstart + 1
          ! Which part of the global domain this represents
          subdomains(ith)%xstart = ival
          subdomains(ith)%xstop = ival + subdomains(ith)%internal%nx - 1
          subdomains(ith)%nx  = 2*halo_width + subdomains(ith)%internal%nx

          subdomains(ith)%internal%ystart = halo_width+1
          subdomains(ith)%internal%ystop = &
                    subdomains(ith)%internal%ystart+height-1
          subdomains(ith)%internal%ny = subdomains(ith)%internal%ystop - &
                                        subdomains(ith)%internal%ystart + 1
          ! Which part of the global domain this represents
          subdomains(ith)%ystart = jval
          subdomains(ith)%ystop = jval + subdomains(ith)%internal%ny - 1
          subdomains(ith)%ny = 2*halo_width + subdomains(ith)%internal%ny

          if(print_tiles)then
             WRITE(*,"('subdomain[',I4,'](',I4,':',I4,')(',I4,':',I4,'), "// &
                  & "interior:(',I4,':',I4,')(',I4,':',I4,') ')")       &
                  ith,                                                  &
                  subdomains(ith)%xstart, subdomains(ith)%xstop,        &
                  subdomains(ith)%ystart, subdomains(ith)%ystop,        &
                  subdomains(ith)%internal%xstart, subdomains(ith)%internal%xstop, &
                  subdomains(ith)%internal%ystart, subdomains(ith)%internal%ystop
          end if

          ! Collect some data on the distribution of tile sizes for 
          ! loadbalance info
          nvects = subdomains(ith)%internal%nx * subdomains(ith)%internal%ny
          nvects_sum = nvects_sum + nvects
          nvects_min = MIN(nvects_min, nvects)
          nvects_max = MAX(nvects_max, nvects)

          ! For use when allocating tile-'private' work arrays
          max_tile_width  = MAX(max_tile_width, subdomains(ith)%nx)
          max_tile_height = MAX(max_tile_height, subdomains(ith)%ny)

          ival = subdomains(ith)%xstop + 1
          ith = ith + 1
       end do
       jval = subdomains(ith-1)%ystop + 1
    end do

    ! Print tile-size statistics
    if(print_tiles)then
       WRITE(*,"(/'Mean sub-domain size = ',F8.1,' pts = ',F7.1,' KB')") &
            REAL(nvects_sum)/REAL(ndom), &
            REAL(8*nvects_sum)/REAL(ndom*1024)
       WRITE(*,"('Min,max sub-domain size (pts) = ',I6,',',I6)") &
            nvects_min, nvects_max
       WRITE(*,"('Domain load imbalance (%) =',F6.2)") &
            100.0*(nvects_max-nvects_min)/REAL(nvects_min)
       WRITE (*,"('Max sub-domain dims are ',I4,'x',I4/)") max_tile_width, &
            max_tile_height
    end if

  end function decompose

end module subdomain_mod
