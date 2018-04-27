!------------------------------------------------------------------------------
! BSD 2-Clause License
! 
! Copyright (c) 2017-2018, Science and Technology Facilities Council
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

module parallel_mod
  use mpi
  implicit none

  private

  integer :: mpierr

  !> Our MPI communicator
  integer :: comm
  !> MPI rank of current process
  integer :: rank
  !> Total no. of MPI processes
  integer :: nranks

contains

  !================================================

  subroutine parallel_init()

    call mpi_init(mpierr)
    call mpi_comm_dup(MPI_COMM_WORLD, comm, mpierr)
    call mpi_comm_size(comm, nranks, mpierr)
    call mpi_comm_rank(comm, rank, mpierr)

    if(rank == 0)then
       write (*,*) "Number of MPI ranks: ", nranks
    end if
  end subroutine parallel_init

  !================================================

  subroutine parallel_finalise()
    call mpi_finalize(mpierr)
  end subroutine parallel_finalise

  !================================================

  function decompose(domainx, domainy,   &
                     ndomains, ndomainx, &
                     ndomainy) result(tiles)
    use tile_mod, only: tile_type
    implicit none
    !! Decompose a domain consisting of domainx x domainy points
    !! into a 2D grid.
    !! Returns an array of tiles describing each of the subdomains.
    
    !> The array of tiles this function will return
    type(tile_type), allocatable :: tiles(:)
    !> Dimensions of the domain to decompose
    integer, intent(in) :: domainx, domainy
    !> No. of domains to decompose into
    integer, intent(in), optional :: ndomains
    !> Optional specification of the dimensions of the tiling grid
    integer, intent(in), optional :: ndomainx, ndomainy
    integer :: nx, ny
    INTEGER :: ival, jval ! For tile extent calculation
    integer :: internal_width, internal_height
    INTEGER :: ierr, nwidth
    INTEGER :: ji,jj, ith
    INTEGER :: nthreads       ! No. of OpenMP threads being used
    INTEGER :: jover, junder, idytmp
    INTEGER :: iover, iunder, idxtmp
    ! For doing stats on tile sizes
    INTEGER :: nvects, nvects_sum, nvects_min, nvects_max 
    LOGICAL, PARAMETER :: print_tiles = .FALSE.
    ! Whether to automatically compute the dimensions of the tiling grid
    logical :: auto_tile
    integer :: xlen, ylen
    integer :: ntilex, ntiley
    ! Local var to hold number of sub-domains
    integer :: ndom

    if(.not. present(ndomains))then
       if(.not. present(ndomainx) .or. .not. present(ndomainy))then
          call gocean_stop('decompose: invalid arguments supplied')
       end if
       ndom = ndomainx * ndomainy
       auto_tile = .FALSE.
    else
       ndom = ndomains
       auto_tile = .TRUE.
    end if

    allocate(tiles(ndom), Stat=ierr)
    if(ierr /= 0)then
       call gocean_stop('decompose: failed to allocate tiles array')
    end if

!!$    xlen = fld%internal%nx
!!$    ylen = fld%internal%ny
!!$    
!!$    IF(auto_tile)THEN
!!$
!!$       ntilex = INT( SQRT(REAL(fld%ntiles)) )
!!$       DO WHILE(MOD(fld%ntiles,ntilex) /= 0)
!!$          ntilex = ntilex - 1
!!$       END DO
!!$       ntiley = fld%ntiles/ntilex
!!$
!!$       ! Match longest dimension of MPI domain to longest dimension of 
!!$       ! thread grid
!!$       IF(xlen > ylen)THEN
!!$          IF( ntilex < ntiley )THEN
!!$             ierr   = ntiley
!!$             ntiley = ntilex
!!$             ntilex = ierr
!!$          END IF
!!$       ELSE
!!$          ! N >= M so want nthready >= nthreadx
!!$          IF( ntiley < ntilex )THEN
!!$             ierr   = ntiley
!!$             ntiley = ntilex
!!$             ntilex = ierr
!!$          END IF
!!$       END IF
!!$
!!$    END IF ! automatic determination of tiling grid
!!$
!!$    WRITE (*,"('OpenMP thread tiling using grid of ',I3,'x',I3)") ntilex,ntiley
!!$    write(*,*) 'ntiles for this field = ',fld%ntiles
!!$
!!$    ! Tiles at left and right of domain only have single
!!$    ! overlap. Every other tile has two overlaps. So: 
!!$    ! xlen = (ntilex-2)*(idx-2) + 2*(idx-1)
!!$    !      = ntilex.idx - 2.ntilex - 2.idx + 4 + 2.idx - 2
!!$    !      = ntilex.idx + 2 - 2.ntilex
!!$    !=> idx = (xlen - 2 + 2.ntilex)/ntilex
!!$    ! where idx is the whole width of a tile.
!!$    !idx = NINT(REAL(xlen - 2 + 2*ntilex)/REAL(ntilex))
!!$    !idy = NINT(REAL(ylen - 2 + 2*ntiley)/REAL(ntiley))
!!$    ! Alternatively, if we think about the internal regions of the tiles,
!!$    ! then they should share the domain between them:
!!$    internal_width = NINT(REAL(xlen) / REAL(ntilex))
!!$    internal_height = NINT(REAL(ylen) / REAL(ntiley))
!!$
!!$    ! Integer arithmetic means that ntiley tiles of height idy might
!!$    ! actually span a height greater or less than N. If so, we try and
!!$    ! reduce the height of each row by just one until we've accounted
!!$    ! for the <jover> extra rows.
!!$    !nwidth = (ntiley-2)*(idy-2) + 2*(idy-1)
!!$    nwidth = ntiley * internal_height
!!$    IF(nwidth > ylen)THEN
!!$       jover  = nwidth - ylen
!!$       junder = 0
!!$    ELSE IF(nwidth < ylen)THEN
!!$       jover  = 0
!!$       junder = ylen - nwidth
!!$    ELSE
!!$       jover  = 0
!!$       junder = 0
!!$    END IF
!!$    ! Ditto for x dimension
!!$    !nwidth = (ntilex-2)*(idx-2) + 2*(idx-1)
!!$    nwidth = ntilex * internal_width
!!$    IF(nwidth > xlen)THEN
!!$       iover  = nwidth - xlen
!!$       iunder = 0
!!$    ELSE IF(nwidth < xlen)THEN
!!$       iover  = 0
!!$       iunder = xlen - nwidth
!!$    ELSE
!!$       iover  = 0
!!$       iunder = 0
!!$    END IF
!!$
!!$    ! For AVX (256-bit vector) instructions, I think we want
!!$    ! MOD(idx,4) == 0 idx = idx + (4 - MOD(idx,4))
!!$
!!$    if(print_tiles)then
!!$       WRITE(*,"('Tile width = ',I4,', tile height = ',I4)") &
!!$            internal_width, internal_height
!!$       WRITE(*,"('iover = ',I3,', iunder = ',I3)") iover, iunder
!!$       WRITE(*,"('jover = ',I3,', junder = ',I3)") jover, junder
!!$    end if
!!$
!!$    ith = 1
!!$    ! The starting point of the tiles in y
!!$    jval = fld%internal%ystart
!!$
!!$    nvects_max = 0
!!$    nvects_min = 1000000
!!$    nvects_sum = 0
!!$    max_tile_width  = 0
!!$    max_tile_height = 0
!!$
!!$    IF(print_tiles)WRITE(*,"(/'Tile dimensions:')")
!!$
!!$    DO jj = 1, ntiley, 1
!!$
!!$       ! If necessary, correct the height of this tile row
!!$       IF(jover > 0)THEN
!!$          idytmp = internal_height - 1
!!$          jover = jover - 1
!!$       ELSE IF(junder > 0)THEN
!!$          idytmp = internal_height + 1
!!$          junder = junder - 1
!!$       ELSE
!!$          idytmp = internal_height
!!$       END IF
!!$
!!$       ! The starting point of the tiles in x
!!$       ival = fld%internal%xstart
!!$
!!$       DO ji = 1, ntilex, 1
!!$         
!!$          ! If necessary, correct the width of this tile column
!!$          IF(iover > 0)THEN
!!$             idxtmp = internal_width - 1
!!$             iover = iover - 1
!!$          ELSE IF(iunder > 0)THEN
!!$             idxtmp = internal_width + 1
!!$             iunder = iunder - 1
!!$          ELSE
!!$             idxtmp = internal_width
!!$          END IF
!!$
!!$          if(ji == 1)then
!!$             fld%tile(ith)%whole%xstart    = fld%whole%xstart
!!$             fld%tile(ith)%internal%xstart = ival
!!$          else
!!$             fld%tile(ith)%internal%xstart = ival
!!$             fld%tile(ith)%whole%xstart    = ival
!!$          end if
!!$          
!!$          IF(ji == ntilex)THEN
!!$             fld%tile(ith)%internal%xstop = fld%internal%xstop
!!$             fld%tile(ith)%whole%xstop = fld%whole%xstop
!!$          ELSE
!!$             fld%tile(ith)%internal%xstop =  MIN(fld%internal%xstop-1, &
!!$                                     fld%tile(ith)%internal%xstart + idxtmp - 1)
!!$             fld%tile(ith)%whole%xstop = fld%tile(ith)%internal%xstop
!!$          END IF
!!$          
!!$          if(jj == 1)then
!!$             fld%tile(ith)%whole%ystart    = fld%whole%ystart
!!$             fld%tile(ith)%internal%ystart = jval
!!$          else
!!$             fld%tile(ith)%whole%ystart    = jval
!!$             fld%tile(ith)%internal%ystart = jval
!!$          end if
!!$
!!$          IF(jj /= ntiley)THEN
!!$             fld%tile(ith)%internal%ystop =  MIN(fld%tile(ith)%internal%ystart+idytmp-1, &
!!$                                             fld%internal%ystop-1)
!!$             fld%tile(ith)%whole%ystop = fld%tile(ith)%internal%ystop
!!$          ELSE
!!$             fld%tile(ith)%internal%ystop = fld%internal%ystop
!!$             fld%tile(ith)%whole%ystop = fld%whole%ystop
!!$          END IF
!!$
!!$          IF(print_tiles)THEN
!!$             WRITE(*,"('tile[',I4,'](',I4,':',I4,')(',I4,':',I4,'), "// &
!!$                  & "interior:(',I4,':',I4,')(',I4,':',I4,') ')")       &
!!$                  ith,                                                  &
!!$                  fld%tile(ith)%whole%xstart, fld%tile(ith)%whole%xstop,       &
!!$                  fld%tile(ith)%whole%ystart, fld%tile(ith)%whole%ystop,       &
!!$                  fld%tile(ith)%internal%xstart, fld%tile(ith)%internal%xstop, &
!!$                  fld%tile(ith)%internal%ystart, fld%tile(ith)%internal%ystop
!!$          END IF
!!$
!!$          ! Collect some data on the distribution of tile sizes for 
!!$          ! loadbalance info
!!$          nvects = (fld%tile(ith)%internal%xstop - fld%tile(ith)%internal%xstart + 1) &
!!$                  * (fld%tile(ith)%internal%ystop - fld%tile(ith)%internal%ystart + 1)
!!$          nvects_sum = nvects_sum + nvects
!!$          nvects_min = MIN(nvects_min, nvects)
!!$          nvects_max = MAX(nvects_max, nvects)
!!$
!!$          ! For use when allocating tile-'private' work arrays
!!$          max_tile_width  = MAX(max_tile_width, &
!!$                  (fld%tile(ith)%whole%xstop - fld%tile(ith)%whole%xstart + 1) )
!!$          max_tile_height = MAX(max_tile_height, &
!!$                  (fld%tile(ith)%whole%ystop - fld%tile(ith)%whole%ystart + 1) )
!!$
!!$          ival = fld%tile(ith)%whole%xstop
!!$          ith = ith + 1
!!$       END DO
!!$       jval = fld%tile(ith-1)%whole%ystop
!!$    END DO
!!$
!!$    ! Print tile-size statistics
!!$    if(print_tiles)then
!!$       WRITE(*,"(/'Mean tile size = ',F8.1,' pts = ',F7.1,' KB')") &
!!$            REAL(nvects_sum)/REAL(fld%ntiles), &
!!$            REAL(8*nvects_sum)/REAL(fld%ntiles*1024)
!!$       WRITE(*,"('Min,max tile size (pts) = ',I6,',',I6)") nvects_min,nvects_max
!!$       WRITE(*,"('Tile load imbalance (%) =',F6.2)") &
!!$            100.0*(nvects_max-nvects_min)/REAL(nvects_min)
!!$       WRITE (*,"('Max tile dims are ',I4,'x',I4/)") max_tile_width, &
!!$            max_tile_height
!!$    end if

  end function decompose

end module parallel_mod
