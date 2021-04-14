!------------------------------------------------------------------------------
! BSD 2-Clause License
! 
! Copyright (c) 2019, Science and Technology Facilities Council
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
! Authors: M. Ashworth, S. Pickles and A. R. Porter, STFC Daresbury Laboratory

module parallel_comms_mod
  use kind_params_mod, only: go_wp
  use parallel_utils_mod, only: get_num_ranks, get_rank, parallel_abort,     &
       MSG_UNDEFINED, MSG_REQUEST_NULL, msg_wait, msg_wait_all, get_max_tag, &
       post_receive, post_send, global_sum
  use decomposition_mod, only: subdomain_type, decomposition_type
  implicit none

  private

  integer, parameter :: JPK = 1 !< Only 1 vertical level
  logical, parameter :: DEBUG = .False.
  !> En/Disable debug output on a per-message basis
  logical, parameter :: DEBUG_COMMS = .False.
  logical :: lwp !< Whether or not to write out from this rank

  ! Set by mapcomms
  integer :: halo_depthx, halo_depthy
  integer, parameter :: MAX_HALO_DEPTH = 1

  ! Communications lists (one for sending, one for receiving)

  !     nsend         Number of messages to be sent.
  !     dirsend       Direction.
  !     destination   Destination process id.
  !     isrcsend      X coordinate of source data.
  !     jsrcsend      Y coordinate of source data.
  !     idessend      X coordinate of destination.
  !     jdessend      Y coordinate of destination.
  !     nxsend        Size in X of data to be sent.
  !     nysend        Size in Y of data to be sent.

  !     nrecv         Number of messages to be received.
  !     dirrecv       Direction.
  !     source        Source process id.
  !     idesrecv      X coordinate of destination.
  !     jdesrecv      Y coordinate of destination.
  !     nxrecv        Size in X of data to be received.
  !     nyrecv        Size in Y of data to be received.

  integer, parameter :: MaxComm=16
  integer, save, dimension(MaxComm)         :: dirsend,destination, &
                                               dirrecv,source
  integer, save, dimension(MaxComm) :: isrcsend,jsrcsend,    &
                                       isrcrecv,jsrcrecv,    &
                                       idessend,jdessend,    &
                                       nxsend,nysend,nzsend, &
                                       idesrecv,jdesrecv,    &
                                       nxrecv,nyrecv,nzrecv
  integer, save :: nsend, nrecv

  ! Total number of points in each message
  integer, save, dimension(MaxComm,MAX_HALO_DEPTH) :: nsendp, nsendp2d, &
                                                   nrecvp, nrecvp2d

  ! Process dependent partitioning information.

  !     ielb       Lower (west) longitude bound index.
  !     ieub       Upper (east) longitude bound index.
  !     iesub      Number of longitude gridpoints.
  !     jelb       Lower (south) latitude bound index.
  !     jeub       Upper (north) latitude bound index.
  !     jesub      Number of latitude gridpoints.

  integer, save :: ielb, ieub, iesub
  integer, save :: jelb, jeub, jesub

  ! Global definitions for parallel execution.

  ! Direction flags for communications.
  ! Listed so that opposite directions are given values maximally spaced
  integer, parameter :: NONE=0         &
                       ,Iplus=1        & 
                       ,Iminus=2       &
                       ,Jplus=3        &
                       ,Jminus=4       &
                       ,IplusJplus=5   &
                       ,IminusJminus=6 &
                       ,IplusJminus=7  &
                       ,IminusJplus=8  &
                       ,MaxCommDir=8
  ! Array to hold direction flags for looking up the
  ! direction that is opposite to one we have
  integer, dimension(MaxCommDir) :: opp_dirn

  ! Set up indices indicating the north-south and east-west
  ! attributes of the eight basic communication directions:
  ! four edges: W, E, S, N;
  ! four corners: SW, SE, NW, NE.
  integer, parameter, dimension(MaxCommDir) ::       &
                west  = (/ 1, 0, 0, 0, 1, 0, 1, 0 /) &
               ,east  = (/ 0, 1, 0, 0, 0, 1, 0, 1 /) &
               ,south = (/ 0, 0, 1, 0, 1, 0, 0, 1 /) &
               ,north = (/ 0, 0, 0, 1, 0, 1, 1, 0 /)
                         ! 1  2  3  4  5  6  7  8
                         ! W  E  S  N SW NE NW  SE
  
  ! exch_flags    -  Array of flag arrays for exchanges
  ! exch_flags1d  -  Array of only the current MPI receive operations
  ! exch_tag      -  The tag value associated with this exchange
  ! exch_busy     -  Indicates whether a slot in the flag array is being used

  integer, parameter :: indexs=1,indexr=2
  integer, parameter :: max_flags=40
  integer, parameter :: min_tag=0
  integer, save :: current_tag,max_tag_used,max_tag,n_tag_cycles=0
  logical, save :: first_mod=.TRUE.

  integer, ALLOCATABLE, dimension(:,:,:), save :: exch_flags
  integer, ALLOCATABLE, dimension(:),     SAVE :: exch_tag, exch_flags1d
  logical, ALLOCATABLE, dimension(:),     SAVE :: exch_busy

  integer, save :: nextFreeExchItem, maxExchItems
  
  ! Buffer for doing halo-exchange.
  ! For a 3D array, halos are 2D slabs but copied into these buffers 
  ! as 1D vectors. 2nd dimension refers to the direction of the 
  ! communication.
  ! For a 2D array, halos are 1D vectors anyway.
  real(go_wp), dimension(:,:), allocatable, save :: sendBuff, recvBuff
  integer , dimension(:,:), allocatable, save :: sendIBuff, recvIBuff

  ! Public routines
  public :: map_comms, iprocmap, exchmod_alloc, exchange_generic, global_sum

  ! Public variables
  public :: MaxComm,nsend,nrecv,nxsend,nysend,destination,dirrecv, &
            dirsend,isrcsend,jsrcsend,idesrecv, jdesrecv,          &
            nxrecv, nyrecv, source, idessend, jdessend

  public :: nsendp,nsendp2d,nrecvp,nrecvp2d
  public :: ielb, ieub, jeub, jelb

  public :: NONE         &
           ,Iplus        & 
           ,Iminus       &
           ,Jplus        &
           ,Jminus       &
           ,IplusJplus   &
           ,IminusJminus &
           ,IplusJminus  &
           ,IminusJplus  &
           ,MaxCommDir

  public :: opp_dirn

contains

  subroutine map_comms (decomp, tmask, pbc, halo_depths, ierr)
    !!------------------------------------------------------------------
    ! Maps out the communications requirements for the partitioned
    ! domain, adding communications descriptions to the list.
    !
    !     Mike Ashworth, CLRC Daresbury Laboratory, July 1999
    !     Andy Porter, STFC Daresbury Laboratory, March 2019
    !!------------------------------------------------------------------
    use parallel_utils_mod, only: DIST_MEM_ENABLED
    implicit none

    ! Subroutine arguments.
    type(decomposition_type), target, intent(in) :: decomp
    integer, intent(in), allocatable :: tmask(:,:)  ! decomposed mask:
                                                    ! 0 for land, 1 for ocean
    logical, intent(in) :: pbc ! Whether mesh is periodic in x dimension
    integer, intent(in) :: halo_depths(2)
    integer, intent(out):: ierr

    ! Local variables.
    integer :: idirn, i1, i2, ihalo, iproc, iprocx, &
         iprocy, j, j1, j2, nadd, naddmaxr, naddmaxs, &
         iinside, ioutside
    integer :: ldiff0, ldiff1 ! Local vars for coping with wrapping of coords
    integer :: imax, imin ! Max/min value of i that a halo strip can run 
                          ! to/from (used to avoid including E/W halo cells
                          ! in N/S transfers)
    integer :: ielb_iproc, ieub_iproc ! Lower and upper bounds for proc
                                      ! iproc corrected for E & W halos
    integer :: jelb_iproc, jeub_iproc
    integer, dimension(MAX_HALO_DEPTH) :: idesr, jdesr, idess, jdess, &
                                isrcr, jsrcr, isrcs, jsrcs, &
                                nxr, nyr, nxs, nys
    logical :: addcorner
    integer :: nprocp, irank
    type(subdomain_type), pointer :: subdomain

    ! Do nothing if distributed-memory support is not enabled
    if(.not. DIST_MEM_ENABLED) return

    halo_depthx = halo_depths(1)
    halo_depthy = halo_depths(2)
    if(halo_depthx > MAX_HALO_DEPTH .OR. halo_depthy > MAX_HALO_DEPTH)then
       call parallel_abort("map_comms: specified halo depth exceeds "// &
                           "MAX_HALO_DEPTH limit in parallel_comms_mod")
    end if
    if( pbc )then
       call parallel_abort("map_comms: periodic-boundary conditions are "// &
                           "not yet supported")
    end if
    
    nprocp = get_num_ranks()
    irank = get_rank()
    lwp = (irank == 1)
    
    ! Clear the error code.
    ierr = 0

    ! Populate the look-up table of opposite directions
    opp_dirn(Iplus       ) = Iminus
    opp_dirn(Iminus      ) = Iplus 
    opp_dirn(Jplus       ) = Jminus
    opp_dirn(Jminus      ) = Jplus
    opp_dirn(IplusJplus  ) = IminusJminus
    opp_dirn(IminusJminus) = IplusJplus
    opp_dirn(IplusJminus ) = IminusJplus
    opp_dirn(IminusJplus ) = IplusJminus

    dirsend = -999
    destination = -999
    isrcsend = -999
    jsrcsend = -999
    idessend = -999
    jdessend = -999
    nxsend = -999
    nysend = -999
    nzsend = -999
    dirrecv = -999
    source = -999
    idesrecv = -999
    jdesrecv = -999
    nxrecv = -999
    nyrecv = -999
    nzrecv = -999

    ! For each of the eight communication directions on a 2d grid of
    ! processes, add a communication to the list.

    ! Note: the parameters Iplus, Iminus, etc. refer to array subscript
    ! references in the code like (I+1), (I-1), etc. which represent
    ! a requirement for communications in the OPPOSITE direction. So
    ! Iplus associates with sending data in the minus I direction,
    ! Iminus associates with sending data in the plus I direction, etc.

    nsend = 0
    nrecv = 0

    ! For convenience take copies of the upper and lower bounds of the
    ! subdomain owned by this process.
    subdomain => decomp%subdomains(irank)
    jelb = subdomain%global%ystart
    jeub = subdomain%global%ystop
    ielb = subdomain%global%xstart
    ieub = subdomain%global%xstop
    iesub = subdomain%internal%nx
    jesub = subdomain%internal%ny

    ! =================================================================
    ! Looking at the border where we will
    ! send data in the minus I direction (Iplus) and
    ! receive data that has been sent in the plus I direction (Iminus).
    ! =================================================================

    ! Start from the lower bound of the sub-domain, and carry on looking
    ! for communications with neighbours until we have reached
    ! the upper bound.

    ! \TODO support multiple subdomains per rank
    j1 = jelb
    do WHILE (j1.LE.jeub)

       ! Look for the process which owns the neighbouring point in the
       ! minus I direction.

       iproc = iprocmap(decomp, ielb-1, j1)
       if ( iproc.GT.0 ) then

          ! Find where in the j direction the common border between these
          ! sub-domains ends.
          j2 = MIN(jeub, decomp%subdomains(iproc)%global%ystop)

          ! ( {i,I}nternal cells, {h,H}alo cells )
          !
          !                    PE=iproc                PE=myself
          !
          !                 nleit(iproc)      data        nldi
          !                                    s          
          !                      |         <--------       |
          !                                    r
          !                                -------->
          !          -----------------------     ---------------------
          !          | In-1 | In | H1 | H2 |     | h2 | h1 | i1 | i2 |
          !          -----------------------     ---------------------
          !                                In -> h1
          !                              In-1 -> h2
          !                                H1 <- i1
          !                                H2 <- i2

          ! Construct the rest of the data describing the zone,
          ! convert to local indices and extend to multiple halo widths.
          isrcs(:) = subdomain%internal%xstart

          do ihalo=1,halo_depthx
            ! Source for the receive must be within internal domain of the
            ! (sending) PE, iproc
            isrcr(ihalo) = decomp%subdomains(iproc)%internal%xstop - ihalo + 1
            idesr(ihalo) = ihalo ! Halo goes from 1..halo_depthx
            nxr(ihalo) = ihalo
            nxs(ihalo) = ihalo
          enddo

          ! Destination for a send must be a halo region
          idess(:) = decomp%subdomains(iproc)%internal%xstop+1
          ! Source for a send is always within internal region
          jsrcs(:) = j1 - subdomain%global%ystart + subdomain%internal%ystart
          jdess(:) = j1 - decomp%subdomains(iproc)%global%ystart + &
                          decomp%subdomains(iproc)%internal%ystart
          jdesr(:) = jsrcs(:)
          jsrcr(:) = jdess(:)
          nyr(:) = j2-j1+1
          nys(:) = nyr(:)

          ! Examine whether corner points should be added to the start.

          naddmaxr = 0
          naddmaxs = 0
          ! ARPDBG - why loop below when naddmaxs and naddmaxr are scalars?
          do ihalo=1,halo_depthx,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.

             if ( j1-ihalo.GE.jelb .AND.  &
                  iprocmap(decomp,ielb-ihalo,j1).GT.0 ) then
                naddmaxs = ihalo
             endif

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            if ( j1.EQ.jelb .AND.  &
                 iprocmap(decomp,ielb-ihalo,j1-ihalo).EQ.iproc ) then
              naddmaxr = ihalo
            endif
          enddo

          ! Add the extra points.

          do ihalo=1,halo_depthx
            nadd = MIN(ihalo,naddmaxs)
            jdess(ihalo) = jdess(ihalo) - nadd
            jsrcs(ihalo) = jsrcs(ihalo) - nadd
            nys(ihalo) = nys(ihalo)+nadd
            if(DEBUG)then
               if ( nadd > 0 ) then
                  write (*,"(I3,': Adding starting corner to send for halo '" &
                         //",I2,' with ',I3,' points')") irank, ihalo, nadd
               end if
            end if
            nadd = MIN(ihalo,naddmaxr)
            jdesr(ihalo) = jdesr(ihalo) - nadd
            jsrcr(ihalo) = jsrcr(ihalo) - nadd
            nyr(ihalo) = nyr(ihalo)+nadd
            if(DEBUG)then
               if ( nadd > 0 ) then
                  write (*,"(I3,': Adding starting corner to receive for " &
                         //"halo ',I2,' with ',I3,' points')") irank,ihalo, nadd
               endif
            end if
          enddo ! Loop over 'overlap' points in i direction

          ! Examine whether corner points should be added to the end.

          naddmaxr = 0
          naddmaxs = 0
          do ihalo=1,halo_depthx,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.

            if ( j2+ihalo.LE.jeub .AND. &
                 iprocmap(decomp,ielb-ihalo,j2).GT.0 ) then
              naddmaxs = ihalo
            endif

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            if ( j2.EQ.jeub .AND. &
                 iprocmap(decomp,ielb-ihalo,j2+ihalo).EQ.iproc ) then
              naddmaxr = ihalo
            endif
          enddo 

          ! Add the extra points.

          do ihalo=1,halo_depthx,1
            nadd = MIN(ihalo,naddmaxs)
            nys(ihalo) = nys(ihalo)+nadd
            if(DEBUG)then
               if ( nadd.GT.0 ) then
                  write (*,"(I3,': Adding starting corner to send for halo '," &
                       //"I2,' with ',I3,' points')") irank,ihalo, nadd
               endif
            end if
            nadd = MIN(ihalo,naddmaxr)
            nyr(ihalo) = nyr(ihalo)+nadd
            if(DEBUG)then
               if ( nadd.GT.0  ) then
                  write (*,"(I3,': Adding starting corner to receive for halo " &
                       //"',I2,' with ',I3,' points')") irank,ihalo, nadd
               endif
            end if
          enddo

          ! Add a send and a receive to the lists for this section
          ! of border.
          CALL addsend (nsend,Iplus,iproc-1, &
                        isrcs,jsrcs,idess,jdess,   &
                        nxs,nys,ierr)
          if ( ierr.NE.0 ) RETURN

          ! This is a receive for data _from_ the neighbouring domain
          ! - it is NOT the corresponding receive for the above send.
          CALL addrecv (nrecv,Iminus,iproc-1, &
                        isrcr,jsrcr,idesr,jdesr,    &
                        nxr,nyr,ierr)
          if ( ierr.NE.0 ) RETURN
          if(DEBUG)then
             if ( lwp ) then
                write (*,'(a,7i6)') 'Adding receive ',iproc-1 &
                   ,isrcr(1),jsrcr(1),idesr(1),jdesr(1),nxr(1),nyr(1)
             endif
          end if

          ! Move the start point to one beyond this strip.
          j1 = j2+1

        else

           ! No process found, continue searching up the boundary.
          j1 = j1+1
        endif
      enddo

      ! =================================================================
      ! Looking at the border where we will
      ! send data in the plus I direction (Iminus) and
      ! receive data that has been sent in the minus I direction (Iplus).
      ! =================================================================

      ! Start from the lower bound of the sub-domain, and carry on looking
      ! for communications with neighbours until we have reached
      ! the upper bound.

      j1 = jelb
      do WHILE (j1.LE.jeub)

         ! Look for the process which owns the neighbouring point in the
         ! plus I direction.

        iproc = iprocmap(decomp,ieub+1,j1)
        if ( iproc.GT.0 ) then

           ! Find where in the j direction the common border between these
           ! sub-domains ends.
          j2 = MIN(jeub,decomp%subdomains(iproc)%global%ystop)

          if(DEBUG)then
             write (*,FMT="(I3,': ARPDBG strip for plus I is ',I3,',',I3)") &
               irank,j1,j2
          end if

          ! ( {i,I}nternal cells, {h,H}alo cells )
          !
          !                    PE=myself                PE=iproc
          !
          !                     nlei          data      nldit(iproc)
          !                                    s          
          !                      |         -------->       |
          !                                    r
          !                                <--------
          !          -----------------------     ---------------------
          !          | In-1 | In | H1 | H2 |     | h2 | h1 | i1 | i2 |
          !          -----------------------     ---------------------
          !                                In -> h1
          !                              In-1 -> h2
          !                                H1 <- i1
          !                                H2 <- i2

          ! Construct the rest of the data describing the zone,
          ! convert to local indexes and extend to multiple halo widths.
          do ihalo=1,halo_depthx
             isrcr(ihalo) = decomp%subdomains(iproc)%internal%xstart + ihalo
             ! Outermost halo -> innermost col.
             isrcs(ihalo) = subdomain%internal%xstop - ihalo + 1
             idess(ihalo) = ihalo ! Halo runs from 1..halo_depthx
             nxr(ihalo) = ihalo
             nxs(ihalo) = ihalo
             ! Destination for a receive is in local halo region
             idesr(ihalo) = subdomain%internal%xstop + ihalo
          enddo
          ! Destination for a send is within halo on remote domain
          jdess(:) = j1 - decomp%subdomains(iproc)%global%ystart  + &
                          decomp%subdomains(iproc)%internal%ystart
          ! Source for a send is within local internal domain
          jsrcs(:) = j1-subdomain%global%ystart + subdomain%internal%ystart

          jdesr(:) = jsrcs(:)
          jsrcr(:) = jdess(:)
          nyr(:) = j2-j1+1
          nys(:) = nyr(:)
          ! Examine whether corner points should be added to the START.

          naddmaxr = 0
          naddmaxs = 0
          do ihalo=1,halo_depthx,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.
             if ( j1-ihalo.GE.jelb .AND.  &
                 iprocmap(decomp,ieub+ihalo,j1).GT.0 ) then
                naddmaxs = ihalo
             endif

             ! Receive corner data only when we are at the corner
             ! and while the sending sub-domain is the same as for the edge.

             if ( j1.EQ.jelb .AND.  &
                 iprocmap(decomp,ieub+ihalo,j1-ihalo).EQ.iproc ) then
                naddmaxr = ihalo
             endif
          enddo

          ! Add the extra points.
          do ihalo=1,halo_depthx,1
            nadd = MIN(ihalo,naddmaxs)
            jdess(ihalo) = jdess(ihalo) - nadd
            jsrcs(ihalo) = jsrcs(ihalo) - nadd
            nys(ihalo) = nys(ihalo)+nadd
            nadd = MIN(ihalo,naddmaxr)
            jdesr(ihalo) = jdesr(ihalo) - nadd
            jsrcr(ihalo) = jsrcr(ihalo) - nadd
            nyr(ihalo) = nyr(ihalo)+nadd
          enddo

          ! Examine whether corner points should be added to the end.
          naddmaxr = 0
          naddmaxs = 0
          do ihalo=1,halo_depthx,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.

            if ( j2+ihalo.LE.jeub .AND. &
                 iprocmap(decomp,ieub+ihalo,j2).GT.0 ) then
              naddmaxs = ihalo
            endif

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            if ( j2.EQ.jeub .AND. &
                 iprocmap(decomp,ieub+ihalo,j2+ihalo).EQ.iproc ) then
              naddmaxr = ihalo
            endif
          enddo 

          ! Add the extra points.

          do ihalo=1,halo_depthx,1
            nadd = MIN(ihalo,naddmaxs)
            nys(ihalo) = nys(ihalo)+nadd
            if(DEBUG)then
               if ( nadd.GT.0 .AND. lwp ) then
                  write (*,*) 'Adding starting corner to send' &
                          ,' for halo ',ihalo,' with ',nadd,' points'
               endif
            end if
            nadd = MIN(ihalo,naddmaxr)
            nyr(ihalo) = nyr(ihalo)+nadd
            if(DEBUG)then
               if ( nadd.GT.0 .AND. lwp ) then
                  write (*,*) 'Adding starting corner to receive' &
                  ,' for halo ',ihalo,' with ',nadd,' points'
               endif
            end if
          enddo

          ! Add a send and a receive to the lists for this section
          ! of border.
          CALL addsend (nsend,Iminus,iproc-1,            &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        ierr)
          if ( ierr.NE.0 ) RETURN

          CALL addrecv (nrecv,Iplus,iproc-1,       &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        ierr)
          if ( ierr.NE.0 ) RETURN

          ! Move the start point to one beyond this strip.
          j1 = j2+1
        else

          ! No process found, continue searching up the boundary.
          j1 = j1+1
        endif
      enddo

      ! =================================================================
      ! Looking at the border where we will send data in the minus J
      ! direction (Jplus) and receive data that has been sent in the plus
      ! J direction (Jminus).
      ! =================================================================

      ! Start from the lower bound of the sub-domain (in global coords), 
      ! and carry on looking
      ! for communications with neighbours until we have reached
      ! the upper bound.
      imin = ielb
      imax = ieub

      i1 = imin
      do WHILE (i1.LE.imax)

         ! Look for the process which owns the neighbouring point in the
         ! minus J direction.

         iproc = iprocmap(decomp,i1,jelb-1)
         if ( iproc.GT.0 ) then

            ! Ensure we don't include halos from the global domain borders if
            ! we have cyclic bc's.
            ieub_iproc = decomp%subdomains(iproc)%global%xstop

            ! Find where in the i direction the common border between these
            ! sub-domains ends (in the global domain).
            i2 = MIN(imax,ieub_iproc)

            if(DEBUG)then
               write (*,FMT="(I3,': ARPDBG strip for minus J is ',I3,',',I3)") &
                 irank-1,i1,i2
            end if

            !  |        |        |        |        ||
            !  |        |        |        |        ||
            !  |        |        |        |        ||
            !  ------------------------------------------------
            !          ||        |        |        |
            !          ||        |        |        |
            !          ||        |        |        |

            ! Construct the rest of the data describing the zone,
            ! convert to local indexes and extend to multiple halo widths.
            ! Convert to local coords:
            ! Convert from global i1 to local i in current domain
            isrcs(:) = i1 - subdomain%global%xstart + subdomain%internal%xstart

            ! Convert from global i1 to local i in the *destination* domain
            idess(:) = i1- decomp%subdomains(iproc)%global%xstart + &
                 decomp%subdomains(iproc)%internal%xstart
            idesr(:) = isrcs(:)
            isrcr(:) = idess(:)
            nxr(:) = i2-i1+1
            nxs(:) = nxr(:)

            jsrcs(:) = subdomain%internal%ystart ! First row of 'internal'
                                                 ! region of domain
            do ihalo=1,halo_depthy,1
               ! Source for a receive must be within internal region
               jsrcr(ihalo) = decomp%subdomains(iproc)%internal%ystop-ihalo+1
               jdesr(ihalo) = ihalo ! Halo runs from 1..halo_depthy
               nyr(ihalo) = ihalo
               nys(ihalo) = ihalo
            enddo
            ! Destination for a send must be a halo region
            jdess(:) = decomp%subdomains(iproc)%internal%ystop + 1

            ! Examine whether corner points should be added to the START.
            naddmaxr = 0
            naddmaxs = 0
            do ihalo=1,halo_depthy,1

               ! Send corner data while we have data to send
               ! and while there is a point that depends on it.
               if ( i1-ihalo.GE.imin .AND.  &
                    iprocmap(decomp,i1-ihalo,jelb-ihalo).GT.0 ) then
                  naddmaxs = ihalo
               endif

               ! Receive corner data only when we are at the corner and
               ! while the sending sub-domain is the same as for the edge.
               if ( i1.EQ.imin .AND. &
                    iprocmap(decomp,i1-ihalo,jelb-ihalo).EQ.iproc ) then
                  naddmaxr = ihalo
               endif
            enddo

            ! Add the extra points.

            do ihalo=1,halo_depthy,1
               nadd = MIN(ihalo,naddmaxs)
               idess(ihalo) = idess(ihalo)-nadd
               isrcs(ihalo) = isrcs(ihalo)-nadd
               nxs(ihalo) = nxs(ihalo)+nadd
               if(DEBUG)then
                  if ( nadd.GT.0 ) then
                     write (*,"(I3,': Adding starting corner to send for &
                       & halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
                  endif
               end if
               nadd = MIN(ihalo,naddmaxr)
               idesr(ihalo) = idesr(ihalo)-nadd
               isrcr(ihalo) = isrcr(ihalo)-nadd
               nxr(ihalo) = nxr(ihalo)+nadd

               if(DEBUG)then
                  if ( nadd.GT.0 ) then
                     write (*,"(I3,': Adding starting corner to receive for &
                       & halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
                  endif
               end if
            enddo

            ! Examine whether corner points should be added to the end.

          naddmaxr = 0
          naddmaxs = 0
          do ihalo=1,halo_depthy,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.

            if ( i2+ihalo.LE.imax .AND. &
                 iprocmap(decomp,i2,jelb-ihalo).GT.0 ) then
              naddmaxs = ihalo
            endif

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            if ( i2.EQ.imax .AND. & 
                 iprocmap(decomp,i2+ihalo,jelb-ihalo).EQ.iproc ) then
              naddmaxr = ihalo
            endif
          enddo 

          ! Add the extra points.

          do ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            nxs(ihalo) = nxs(ihalo)+nadd
            if(DEBUG)then
               if ( nadd.GT.0 ) then
                  write (*,"(I3,': Adding starting corner to send for &
                       & halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
               endif
            end if
            nadd = MIN(ihalo,naddmaxr)
            nxr(ihalo) = nxr(ihalo)+nadd
            if(DEBUG)then
               if ( nadd.GT.0 ) then
                  write (*,"(I3,': Adding starting corner to receive for & 
                       & halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
               endif
            end if
          enddo

          ! Add a send and a receive to the lists for this section
          ! of border.

          CALL addsend (nsend,Jplus,iproc-1,       &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        ierr)
          if ( ierr.NE.0 ) RETURN

          CALL addrecv (nrecv,Jminus,iproc-1,      &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        ierr)
          if ( ierr.NE.0 ) RETURN

          ! Move the start point to one beyond this strip.
          i1 = i2+1

        else

           ! No process found, continue searching up the boundary.
          i1 = i1+1
        endif
      enddo

      ! =================================================================
      ! Looking at the border where we will
      ! send data in the plus J direction (Jminus) and
      ! receive data that has been sent in the minus J direction (Jplus).
      ! =================================================================

      ! Start from the lower bound of the sub-domain, and carry on looking
      ! for communications with neighbours until we have reached
      ! the upper bound.
      imin = ielb
      imax = ieub
      i1 = imin

      do WHILE (i1.LE.imax)

         ! Look for the process which owns the neighbouring point in the
         ! plus J direction.

         iproc = iprocmap(decomp,i1,jeub+1)
         if ( iproc.GT.0 ) then
           ielb_iproc = decomp%subdomains(iproc)%global%xstart
           ieub_iproc = decomp%subdomains(iproc)%global%xstop

           ! Find where in the i direction the common border between these
           ! sub-domains ends (in the global domain).
           i2 = MIN(imax, ieub_iproc)

           if(DEBUG)then
              write (*,FMT="(I3,': ARPDBG strip for plus J is ',I3,',',I3)") &
                 irank-1,i1,i2
           end if

           ! Construct the rest of the data describing the zone,
           ! convert to local indexes and extend to multiple halo widths.
           isrcs(:) = i1 - subdomain%global%xstart + subdomain%internal%xstart
           idess(:) = i1 - decomp%subdomains(iproc)%global%xstart + &
                decomp%subdomains(iproc)%internal%xstart
           idesr(:) = isrcs(:)
           isrcr(:) = idess(:)
           nxr(:) = i2-i1+1
           nxs(:) = nxr(:)

            ! Source for a receive must be within an internal region
            jsrcr(:) = decomp%subdomains(iproc)%internal%ystart

            do ihalo=1,halo_depthy,1
               ! innermost row -> outermost halo
               jsrcs(ihalo) = subdomain%internal%ystop-ihalo+1
               jdess(ihalo) = ihalo        ! Halo runs from 1..halo_depthy
               nyr(ihalo)   = ihalo
               nys(ihalo)   = ihalo
            enddo
            ! Destination for receive must be in halo
            jdesr(:) = subdomain%internal%ystop + 1

          ! Examine whether corner points should be added to the START.

          naddmaxr = 0
          naddmaxs = 0
          do ihalo=1,halo_depthy,1

            ! Send corner data while we have data to send
            ! and while there is a point that depends on it.
            if(i1-ihalo >= imin .and. &
               iprocmap(decomp, i1, jeub+ihalo) > 0)then
              naddmaxs = ihalo
            endif

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.
            if ( i1.EQ.imin .AND. &
                 iprocmap(decomp,i1-ihalo,jeub+ihalo).EQ.iproc ) then
              naddmaxr = ihalo
            endif
          enddo

          ! Add the extra points.

          do ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            idess(ihalo) = idess(ihalo) -nadd
            isrcs(ihalo) = isrcs(ihalo) -nadd
            nxs(ihalo) = nxs(ihalo)+nadd
            if(DEBUG)then
               if ( nadd.GT.0 ) then
                  write (*,"(I3,': Adding starting corner to send for halo '," &
                  //"I2,' with ',I3,' points')") irank-1,ihalo, nadd
               endif
            end if
            nadd = MIN(ihalo,naddmaxr)
            idesr(ihalo) = idesr(ihalo) - nadd
            isrcr(ihalo) = isrcr(ihalo) - nadd
            nxr(ihalo) = nxr(ihalo)+nadd

            if(DEBUG)then
               if ( nadd.GT.0 ) then
                  write (*,"(I3,': Adding starting corner to receive for " &
                       //"halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
               endif
            end if
          enddo

          ! Examine whether corner points should be added to the end.
          naddmaxr = 0
          naddmaxs = 0
          do ihalo=1,halo_depthy,1

            ! Send corner data while we have data to send
            ! and while there is a point that depends on it.
            if ( i2+ihalo.LE.imax .AND. &
                 iprocmap(decomp,i2,jeub+ihalo).GT.0 ) then
              naddmaxs = ihalo
            endif

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            if ( i2.EQ.imax .AND.       & ! Are we at the corner?
                 iprocmap(decomp,i2+ihalo,jeub+ihalo).EQ.iproc ) then
              naddmaxr = ihalo
            endif
          enddo 

          ! Add the extra points.

          do ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            nxs(ihalo) = nxs(ihalo)+nadd

            if(DEBUG)then
               if ( nadd.GT.0 .AND. lwp ) then
                  write (*,*) irank-1,': Adding starting corner to send' &
                       ,' for halo ',ihalo,' with ',nadd,' points'
               endif
            end if
            
            nadd = MIN(ihalo,naddmaxr)
            nxr(ihalo) = nxr(ihalo)+nadd
            if(DEBUG)then
               if ( nadd.GT.0 .and. lwp ) then
                  write (*,*) irank-1,': Adding starting corner to receive' &
                  ,' for halo ',ihalo,' with ',nadd,' points'
               endif
            end if
          enddo

          ! Add a send and a receive to the lists for this section
          ! of border.
          CALL addsend (nsend,Jminus,iproc-1,      &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        ierr)
          if ( ierr.NE.0 ) RETURN

          CALL addrecv (nrecv,Jplus,iproc-1,       &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        ierr)
          if ( ierr.NE.0 ) RETURN

          ! Move the start point to one beyond this strip.
          i1 = i2+1

       else ! iproc < 0

           ! No process found, continue searching up the boundary.
          i1 = i1+1
        endif
      enddo

      ! =================================================================
      ! Corner points are sent with the edge data where possible.
      ! Check for cases where corner data resides on a processor which
      ! is not participating in edge communications and set up extra
      ! diagonal communications for these cases.
      ! =================================================================

      ! Loop over the four corner directions
      ! i = 1  2  3  4  5  6  7 8
      !     W  E  S  N SW NE NW SE
      do idirn = 5, 8

        ! At first assume there is to be no corner communication
        addcorner = .FALSE.

        ! ioutside is to be x-coord just OUTSIDE our domain
        ! iinside is to be x-coord just WITHIN our domain
        ioutside = west(idirn)*(subdomain%global%xstart-1) + &
             east(idirn)*(subdomain%global%xstop+1)
        iinside = west(idirn)*subdomain%global%xstart + &
             east(idirn)*subdomain%global%xstop

        ! For e.g. a NW corner:
        !               | 
        !       iproc   |  iprocy
        !       ________|______
        !               |
        !       iprocx  |  Me
        !               |

        ! x coord just OUTSIDE our domain but y INSIDE
        iprocx = iprocmap(decomp, ioutside, &
             south(idirn)*jelb+north(idirn)*jeub)
        ! x coord just INSIDE our domain but y OUTSIDE
        iprocy = iprocmap(decomp,iinside, &
             south(idirn)*(jelb-1)+north(idirn)*(jeub+1))

        ! Loop over all required halo widths

        do ihalo=halo_depthx,1,-1

          ! Look at the processor in the corner at width ihalo from the
          ! corner of the sub-domain. We want both x and y to be just
          ! outside our domain. ioutside is already just outside our domain
          ! so we subtract one from ihalo below:
          iproc = iprocmap(decomp, &
                    ioutside - west(idirn)*(ihalo-1)+ east(idirn)*(ihalo-1), &
                    south(idirn)*(jelb-ihalo)+north(idirn)*(jeub+ihalo))

          ! If the corner processor is different from those to X and Y
          ! we will need a corner communication.
          if ( iproc.GT.0 .AND. iprocx.GT.0 .AND. iprocy.GT.0 .AND. &
               iproc.NE.iprocx .AND. iproc.NE.iprocy ) then

             ielb_iproc = decomp%subdomains(iproc)%global%xstart
             ieub_iproc = decomp%subdomains(iproc)%global%xstop
             jelb_iproc = decomp%subdomains(iproc)%global%ystart
             jeub_iproc = decomp%subdomains(iproc)%global%ystop

             if(DEBUG)then
                write(*, &
                     FMT="(I3,': adding corner as ',I3,' differs from ',2I3)")&
                    irank-1, iproc,iprocx,iprocy
             end if

             ! If the furthest corner point needs a communication,
             ! we will need them all.
             if (ihalo .EQ. halo_depthx) then
               ielb_iproc = decomp%subdomains(iproc)%global%xstart
               ieub_iproc = decomp%subdomains(iproc)%global%xstop

               ! Set the flag to add everything to the communications list.
               addcorner = .TRUE.
             endif

             ! Set the parameters for the communication.
             !     __|__
             !  __|__|
             !    |
             !    |
             ldiff0 = ielb_iproc - ieub
             ldiff1 = ielb - ieub_iproc

             nxs(ihalo) = ihalo -  east(idirn)*(ldiff0-1) &
                                -  west(idirn)*(ldiff1-1)
             ldiff0 = jelb_iproc - jeub
             ldiff1 = jelb - jeub_iproc
             nys(ihalo) = ihalo - north(idirn)*(ldiff0-1) &
                                - south(idirn)*(ldiff1-1)

             ! Source for a send must be within the *internal* region of
             ! the LOCAL domain
             isrcs(ihalo) = east(idirn)*(subdomain%internal%xstop - ihalo + 1) + &
                  west(idirn)*(subdomain%internal%xstart + ihalo - 1)
             jsrcs(ihalo) = north(idirn)*(subdomain%internal%ystop - ihalo + 1) + &
                  south(idirn)*(subdomain%internal%ystart + ihalo - 1)

             ! Destination for a send must be within a halo region on the
             ! REMOTE domain
             idess(ihalo) =  west(idirn)*( &
                  decomp%subdomains(iproc)%internal%xstop + ihalo) + &
                  east(idirn)*(decomp%subdomains(iproc)%internal%xstart-ihalo)

             jdess(ihalo) = south(idirn)*( &
                  decomp%subdomains(iproc)%internal%ystop + ihalo) + &
               north(idirn)*(decomp%subdomains(iproc)%internal%ystart - ihalo)

             ! Source for a receive must be in an internal region of the
             ! REMOTE domain
             isrcr(ihalo) = west(idirn)* &
                 (decomp%subdomains(iproc)%internal%xstop - ihalo + 1) + &
                 east(idirn)*(decomp%subdomains(iproc)%internal%xstart + ihalo -1)
             jsrcr(ihalo) = south(idirn)* &
                  (decomp%subdomains(iproc)%internal%ystop - ihalo + 1) + &
                  north(idirn)*(decomp%subdomains(iproc)%internal%ystart + ihalo - 1)

             ! Destination for a receive must be in a halo region (on LOCAL
             ! domain) and therefore:
             !  1 <= {i,j}desr <= jprec{i,j} or >= nle{i,j} + 1
             idesr(ihalo) =  east(idirn)*(subdomain%internal%xstop + ihalo) &
                         + west(idirn)*(subdomain%internal%xstart-ihalo)

             jdesr(ihalo) = north(idirn)*(subdomain%internal%ystop + ihalo) &
                         + south(idirn)*(subdomain%internal%ystart - ihalo)

          else

             if(DEBUG)then
                write (*,FMT="(I3,': skipping corner for halo ',I3,' PE ',I3)")&
                       irank-1,ihalo,iproc
             end if

            ! No communication for this halo width.
            ! Clear all the parameters
            ! in case there are comms for other halo widths.
            isrcs(ihalo) = 0
            jsrcs(ihalo) = 0
            idess(ihalo) = 0
            jdess(ihalo) = 0
            isrcr(ihalo) = 0
            jsrcr(ihalo) = 0
            idesr(ihalo) = 0
            jdesr(ihalo) = 0
            nxs  (ihalo) = 0
            nys  (ihalo) = 0
          endif

       enddo ! Loop over halo depths

       ! The size of the received data is always the same as 
       ! that of the sent data.
       nxr(:) = nxs(:)
       nyr(:) = nys(:)

        ! Add the data to the communications lists.
        if ( addcorner ) then
          if(DEBUG)then
             write (*,FMT="(I3,': ARPDBG adding corner send to ',I4, &
                  & ', dir = ',I1)") irank-1, iproc-1,idirn
          end if
          call addsend (nsend,idirn,iproc-1,            &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        ierr)
          if ( ierr.NE.0 ) then
            call parallel_abort("Call to addsend failed.")
          end if

          ! Manually reverse the direction indicator for the receive.
          j = opp_dirn(idirn)

          if(DEBUG)then
             write (*,FMT="(I3,': ARPDBG adding corner recv. from ',I4,', &
                & old dir = ',I1,' new dir = ',I1)") irank-1, iproc-1,idirn, j
          end if
          call addrecv (nrecv,j,iproc-1,          &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        ierr)
          if ( ierr.NE.0 ) then
             call parallel_abort("Call to addrecv failed")
          end if
        endif

      enddo

    end subroutine map_comms

  subroutine addsend (icomm, dir, proc, isrc, jsrc, &
                      ides, jdes, nx, ny, ierr )
    !!------------------------------------------------------------------
    !   Adds a send communication specified by the parameters dir through 
    !   to ny to the send communication list at the next position. 
    !   icomm points to the last entry and is incremented and returned 
    !   if successful.
    !
    !   icomm           int   in/out    Location in comms list.
    !   dir             int   input     Direction.
    !   proc            int   input     Process id.
    !   isrc            int   input     X coordinate of source data.
    !   jsrc            int   input     Y coordinate of source data.
    !   ides            int   input     X coordinate of destination data.
    !   jdes            int   input     Y coordinate of destination data.
    !   nx              int   input     Size in X of data to be sent.
    !   ny              int   input     Size in Y of data to be sent.
    !   ierr            int   output    Error flag.
    !
    !             Mike Ashworth, CLRC Daresbury Laboratory, March 1999
    !             Stephen Pickles, STFC Daresbury Laboratory
    !                - Sep 2009: added run-length encoding
    !!------------------------------------------------------------------
    implicit none

    ! Subroutine arguments.
    integer,                    intent(inout) :: icomm
    integer,                    intent( in  ) :: dir, proc
    integer,                    intent( out ) :: ierr
    integer, dimension(:), intent( in  ) :: isrc, jsrc, ides, jdes, nx, ny
    integer :: irank
         
    irank = get_rank()
    ! Clear the error flag.
    ierr = 0

    ! Return if the process id is not set.
    if ( proc.LT.0 ) then
       RETURN
    endif

    icomm = icomm+1

    ! Check that the comms list has space for another entry.
    if ( icomm.GT.MaxComm ) then
       if ( lwp ) then
          write (*,*) 'ERROR: Number of separate ', &
               'communications exceeds maximum of ',MaxComm
       endif
       ierr = -12
       RETURN
    endif
            
    ! Add the data into the comms list at the new location.
    dirsend(icomm)     = dir
    destination(icomm) = proc
    isrcsend(icomm)    = isrc(1)
    jsrcsend(icomm)    = jsrc(1)
    idessend(icomm)    = ides(1)
    jdessend(icomm)    = jdes(1)
    nxsend(icomm)      = nx(1)
    nysend(icomm)      = ny(1)
    nsendp2d(icomm, 1) = nxsend(icomm)*nysend(icomm)
    nzsend(icomm)      = 1
    nsendp(icomm, 1) = nsendp2d(icomm,1)*nzsend(icomm)
    
    if(DEBUG)then
       write(*,FMT="(I4,': ARPDBG adding Send:')") irank-1  
       write(*,FMT="(I4,': ARPDBG: icomm = ',I2)") irank-1,icomm
       write(*,FMT="(I4,': ARPDBG:   dir = ',I2)") irank-1,dirsend(icomm)
       write(*,FMT="(I4,': ARPDBG:  proc = ',I4)") irank-1,destination(icomm)
       write(*,FMT="(I4,': ARPDBG:  isrc = ',I4)") irank-1,isrcsend(icomm)
       write(*,FMT="(I4,': ARPDBG:  jsrc = ',I4)") irank-1,jsrcsend(icomm)
       write(*,FMT="(I4,': ARPDBG:  ides = ',I4)") irank-1,idessend(icomm)
       write(*,FMT="(I4,': ARPDBG:  jdes = ',I4)") irank-1,jdessend(icomm)
       write(*,FMT="(I4,': ARPDBG:    nx = ',I4)") irank-1,nxsend(icomm)
       write(*,FMT="(I4,': ARPDBG:    ny = ',I4)") irank-1,nysend(icomm)
       write(*,FMT="(I4,': ARPDBG:    nz = ',I4)") irank-1,nzsend(icomm) 
       write (*,FMT="(I4,': ARPDBG:nsendp = ',I4)") irank-1,nsendp(icomm,1)
       write (*,FMT="(I4,': ARPDBG Send ends')")    irank-1
    end if

    if(nxsend(icomm) < 1)then
       write(*,"(I4,': ERROR in addsend: nxsend = ',I4)") irank-1, &
            nxsend(icomm)
       call parallel_abort("")
    end if
    if(nysend(icomm) < 1)then
       write(*,"(I4,': ERROR in addsend: nysend = ',I4)") irank-1, &
            nysend(icomm)
       call parallel_abort("")
    end if

  end subroutine addsend

  subroutine addrecv ( icomm, dir, proc, isrc, jsrc, &
                       ides, jdes, nx, ny, ierr )
    !!------------------------------------------------------------------
    !   Adds a recv communication specified by the parameters dir through 
    !   to ny to the recv communication list at the next position. 
    !   icomm points to the last entry and is incremented and returned 
    !   if successful.
    !
    !   icomm                   int   in/out    Location in comms list.
    !   dir                     int   input     Direction.
    !   proc                    int   input     Process id.
    !   isrc                    int   input     X coordinate of source data.
    !   jsrc                    int   input     Y coordinate of source data.
    !   ides                    int   input     X coordinate of dest. data.
    !   jdes                    int   input     Y coordinate of dest. data.
    !   nx                      int   input     Size in X of data to be sent.
    !   ny                      int   input     Size in Y of data to be sent.
    !   ierr                    int   output    Error flag.
    !
    !          Mike Ashworth, CLRC Daresbury Laboratory, March 1999
    !!------------------------------------------------------------------
    implicit none

    ! Subroutine arguments.
    integer,                 intent(inout) :: icomm
    integer,                 intent( in  ) :: dir, proc
    integer,                 intent(out)   :: ierr
    integer, dimension(halo_depthx)        :: isrc, jsrc, ides, jdes, nx, ny
    ! Local variables.
    integer :: irank

    irank = get_rank()
    ! Clear the error flag.
    ierr = 0

    ! Return if the process id is not set.
    if ( proc.LT.0 ) then
       RETURN
    endif

    icomm = icomm + 1

    ! Check that the comms list has space for another entry.

    if ( icomm.GT.MaxComm ) then
       if ( lwp ) then
          write (*,*) 'ERROR: Number of separate ', &
               'communications exceeds maximum of ',MaxComm
       endif
       ierr = -12
       RETURN
    endif

    ! Add the data into the comms list at the new location.
    dirrecv(icomm)  = dir
    source(icomm)   = proc
    isrcrecv(icomm) = isrc(1)
    jsrcrecv(icomm) = jsrc(1)
    idesrecv(icomm) = ides(1)
    jdesrecv(icomm) = jdes(1)
    nxrecv(icomm)   = nx(1)
    nyrecv(icomm)   = ny(1)
    !> TOdo handle multiple halo depths
    nrecvp2d(icomm, 1) = nxrecv(icomm)*nyrecv(icomm)
    nzrecv(icomm)   = 1
    nrecvp(icomm, 1) = nzrecv(icomm) * nrecvp2d(icomm, 1)
    
    if(DEBUG)then
       write (*,FMT="(I4,': ARPDBG adding RECV:')") irank-1
       write (*,FMT="(I4,': ARPDBG: icomm = ',I2)") irank-1,icomm
       write (*,FMT="(I4,': ARPDBG:   dir = ',I2)") irank-1,dir
       write (*,FMT="(I4,': ARPDBG:  proc = ',I4)") irank-1,proc
       write (*,FMT="(I4,': ARPDBG:  isrc = ',I4)") irank-1,isrcrecv(icomm)
       write (*,FMT="(I4,': ARPDBG:  jsrc = ',I4)") irank-1,jsrcrecv(icomm)
       write (*,FMT="(I4,': ARPDBG:  ides = ',I4)") irank-1,idesrecv(icomm)
       write (*,FMT="(I4,': ARPDBG:  jdes = ',I4)") irank-1,jdesrecv(icomm)
       write (*,FMT="(I4,': ARPDBG:    nx = ',I4)") irank-1,nxrecv(icomm)
       write (*,FMT="(I4,': ARPDBG:    ny = ',I4)") irank-1,nyrecv(icomm)
       write (*,FMT="(I4,': ARPDBG:nrecvp = ',I4)") irank-1,nrecvp(icomm,1)
       write (*,FMT="(I4,': ARPDBG:nrecvp2d = ',I4)") irank-1,nrecvp2d(icomm,1)
       write (*,FMT="(I4,': ARPDBG RECV ends')")    irank-1
    end if

    if(nxrecv(icomm) < 1)then
       write(*,"(I4,': ERROR in addrecv: nxrecv = ',I4)") irank-1, &
            nxrecv(icomm)
       call parallel_abort("")
    end if
    if(nyrecv(icomm) < 1)then
       write(*,"(I4,': ERROR in addrecv: nyrecv = ',I4)") irank-1, &
            nyrecv(icomm)
       call parallel_abort("")
    end if

  end subroutine addrecv

  function iprocmap (decomp, ia, ja)
    !!------------------------------------------------------------------
    !> Returns the process number (1...N) of the process whose sub-domain
    !! contains the point with global coordinates (ia, ja).
    !! If no process owns the point, returns zero.

    !         Mike Ashworth, CLRC Daresbury Laboratory, July 1999
    !         Andrew Porter, STFC Daresbury Laboratory, May  2008
    !!------------------------------------------------------------------
    implicit none
    !> The process owning the supplied point or zero if none.
    integer                              :: iprocmap
    !> The domain decomposition
    type(decomposition_type), intent(in) :: decomp
    !> Global x, y grid indices of the point to search for
    integer,                  intent(in) :: ia, ja
    ! Local variables.
    integer :: iproc, nprocp

    nprocp = get_num_ranks()
    iprocmap = 0

    ! Search the processes for the one owning point (i,j).
    do iproc=1,nprocp
       if (decomp%subdomains(iproc)%global%xstart <= ia .and. &
           ia <= decomp%subdomains(iproc)%global%xstop  .and. &
           decomp%subdomains(iproc)%global%ystart <= ja .and. &
           ja <= decomp%subdomains(iproc)%global%ystop) then
          iprocmap = iproc
          exit
       end if
    enddo

  end function iprocmap

  ! ---------------------------------------------------------------
  !> Allocate the various arrays required to hold the state of
  !! halo exchanges.
  integer function exchmod_alloc()
    implicit none
    ! Locals
    integer :: ierr
    maxExchItems = 20
    allocate(exch_flags(max_flags,MaxComm,2),   &
             exch_flags1d(MaxComm),             &
             exch_busy(max_flags),              &
             exch_tag(max_flags),               &
             STAT=ierr)

    if(ierr .eq. 0)then
       exch_busy   = .FALSE.
    else
       maxExchItems = 0
    end if

    nextFreeExchItem = 1

    ! Pass back the allocation status flag
    exchmod_alloc = ierr

  end function exchmod_alloc

  ! ---------------------------------------------------------------
  !> Gets a new exchange handle
  integer function get_exch_handle ( )
    implicit none
    ! Local variables.
    integer :: h

    if ( first_mod ) then

       ! First time in the module (i.e. exch or glob) set up the tags.
       max_tag = get_max_tag()
       if(DEBUG)then
          if ( lwp ) write (*,*) 'MAX_TAG: set to ', max_tag
       end if

       ! Set the current tag to the minimum.
       current_tag = min_tag
       max_tag_used = current_tag

       first_mod = .FALSE.
    endif

    ! Look for a free location in the flags array.

    flag_search : do h=1,max_flags
       if ( .NOT.exch_busy(h) ) EXIT flag_search
    enddo flag_search

    if ( h.GT.max_flags ) then

       ! If no free flags array was found, flag an error.
       call parallel_abort('ERROR: get_exch_handle: no free flag array')
    else

       ! Assign a new tag.
       exch_busy(h) = .TRUE.

       if ( current_tag.GE.(max_tag-MaxCommDir) ) then

          ! Wrap around.
          current_tag = min_tag
          n_tag_cycles = n_tag_cycles+1
       else
          current_tag = current_tag + MaxCommDir
       endif
       max_tag_used = MAX(max_tag_used,current_tag)
       exch_tag(h) = current_tag

    endif

    get_exch_handle = h

  end function get_exch_handle

  !================================================================

  subroutine free_exch_handle ( h )
    !> Frees exchange handle, h.
    implicit none

    ! Subroutine arguments.
    integer :: h ! Handle to be free'd

    ! Free the flags array.
    if ( h.GT.0 .AND. h.LE.max_flags ) then
       exch_busy(h) = .FALSE.
    else
       write (*,*) 'free_exch_handle: invalid handle ',h
    endif

  end subroutine free_exch_handle

  !================================================================

  subroutine exchange_generic ( b2, ib2, b3, ib3, &
                                handle, comm1, comm2, comm3, comm4)

    ! ******************************************************************
    ! Send and receive boundary data elements to adjacent sub-domains.

    ! b2(:,:)                real   inout       2D real*8 local array.
    ! ib2(:,:)               int    inout       2D integer local array.
    ! b3(:,:,:)              real   inout       3D real*8 local array.
    ! ib3(:,:,:)             int    inout       3D integer local array.
    ! rows/cols to exchange.
    ! handle                 int    output      Exchange handle.
    ! comm1                  int    input       Send in direction comm1.
    ! comm2                  int    input       Send in direction comm2.
    ! comm3                  int    input       Send in direction comm3.
    ! comm4                  int    input       Send in direction comm4.
    !
    ! Mike Ashworth, CCLRC, March 2005.
    ! Andrew Porter, STFC,  January 2008
    ! ******************************************************************
    use parallel_utils_mod, only: DIST_MEM_ENABLED
    IMPLICIT none

    ! Subroutine arguments.
    integer,               intent(out)                     :: handle
    real(go_wp), optional, intent(inout), dimension(:,:)   :: b2
    integer,     optional, intent(inout), dimension(:,:)   :: ib2
    real(go_wp), optional, intent(inout), dimension(:,:,:) :: b3
    integer,     optional, intent(inout), dimension(:,:,:) :: ib3
    integer,               intent(in) :: comm1, comm2, comm3, comm4

    ! Local variables.

    logical :: enabled(0:MaxCommDir)
    integer :: ierr, irecv, isend, &
               isrc, jsrc, nxs, nys, tag, tag_orig
    integer :: maxrecvpts, maxsendpts ! Max no. of grid points involved in 
                                      ! any one halo exchange
    integer :: i, j, k, ic ! Loop counters
    integer :: istart, iend, jstart, jend
    logical, save :: first_time = .TRUE.
    integer, parameter :: index_z = 3
    !!--------------------------------------------------------------------

    ! Do nothing if distributed-memory support is not enabled
    if(.not. DIST_MEM_ENABLED) return
    
    ierr = 0

    ! Allocate a communications tag/handle and a flags array.

    handle   = get_exch_handle()
    tag_orig = exch_tag(handle)

    ! Set enabled flags according to the subroutine arguments.

    enabled(Iplus ) = .FALSE.
    enabled(Jplus ) = .FALSE.
    enabled(Iminus) = .FALSE.
    enabled(Jminus) = .FALSE.
    enabled(comm1) = comm1 > 0
    enabled(comm2) = comm2 > 0
    enabled(comm3) = comm3 > 0
    enabled(comm4) = comm4 > 0

    ! Set diagonal communications according to the non-diagonal flags.

    enabled(IplusJplus ) = enabled(Iplus ).AND.enabled(Jplus )
    enabled(IminusJminus)= enabled(Iminus).AND.enabled(Jminus)
    enabled(IplusJminus) = enabled(Iplus ).AND.enabled(Jminus)
    enabled(IminusJplus )= enabled(Iminus).AND.enabled(Jplus )

    ! Main communications loop.
    maxrecvpts = maxval(nrecvp(1:nrecv,1))
    maxsendpts = maxval(nsendp(1:nsend,1))

    if(present(b2) .OR. present(b3))then
       if(.NOT. allocated(sendBuff))then
          ! Only allocate the sendBuff once
          allocate(recvBuff(maxrecvpts,nrecv), &
                   sendBuff(maxsendpts,nsend),stat=ierr)
       else
          allocate(recvBuff(maxrecvpts,nrecv),stat=ierr)
       end if
    else if(present(ib2) .OR. present(ib3))then
       if(.NOT. allocated(sendIBuff))then
          allocate(recvIBuff(maxrecvpts,nrecv), &
                   sendIBuff(maxsendpts,nsend),stat=ierr)
       else
          allocate(recvIBuff(maxrecvpts,nrecv),stat=ierr)
       end if
    end if

    if (ierr /= 0) then
       call parallel_abort('exchange_generic: unable to allocate send/recvBuffs')
    end if

    ! Initiate receives in case posting them first improves 
    ! performance.

    do irecv=1,nrecv

       if ( enabled(dirrecv(irecv)) .AND. &
            source(irecv) >= 0 .AND. nxrecv(irecv) > 0 ) then

          tag = tag_orig + dirrecv(irecv)

          ! ARPDBG - nrecvp second rank is for multiple halo widths but
          !          that isn't used
          call post_receive(nrecvp2d(irecv,1), source(irecv), tag, &
               exch_flags(handle,irecv,indexr), rbuff=recvBuff(:,irecv))

          if(DEBUG_COMMS)then
             write (*,FMT="(I4,': exchs post recv : hand = ',I2,' dirn = '," &
             //"I1,' src = ',I3,' tag = ',I4,' npoints = ',I6)") &
                  get_rank(),handle,dirrecv(irecv), &
                  source(irecv), tag, nrecvp(irecv,1)
          endif

       else
          exch_flags(handle,irecv,indexr) = MSG_REQUEST_NULL
       end if

    end do

    if (.not. first_time) then        

       ! Check that all sends from previous call have completed before 
       ! we continue to modify the send buffers
       call msg_wait_all(nsend, exch_flags1d)
    else
        first_time = .FALSE.
    end if ! .not. first_time

    ! Send all messages in the communications list.

    do isend=1,nsend

       if ( enabled(dirsend(isend)) .AND. &
            destination(isend) >= 0 .AND. nxsend(isend) > 0 ) then

          isrc = isrcsend(isend)
          jsrc = jsrcsend(isend)
          nxs  =   nxsend(isend)
          nys  =   nysend(isend)

          tag = tag_orig + dirsend(isend)

          if( DEBUG_COMMS)then
             if(present(b3))then
                write (*,FMT="(I4,': handle ',I4,' tag ',I4,' sending to '," &
                     //"I4,' data ',I4,' direction ',I3)")                   &  
                     get_rank(), handle, tag, destination(isend),            &
                     nsendp(isend,1),dirsend(isend)
             else if(present(b2))then
                write (*,FMT="(I4,': handle ',I4,' tag ',I4,' sending to '," &
                     //"I4,' data ',I4,' direction ',I3)")                   &  
                     get_rank(), handle, tag, destination(isend),            &
                     nsendp2d(isend,1),dirsend(isend)
             end if
          endif

          ! Copy the data into the send buffer and send it...
          if ( present(b2) )then

             ic = 0
             ! Stored patch coordinates are for T points
             istart = isrcsend(isend)
             iend = istart+nxsend(isend)-1
             jstart = jsrcsend(isend)
             jend = jstart+nysend(isend)-1

             if(DEBUG_COMMS)then
                write(*,"(I3,': packing from:',I3,':',I3,',',I3,':',I3)") &
                  get_rank(),istart,iend,jstart,jend
             end if
             
             do j=jstart, jend, 1
                do i=istart, iend, 1
                   ic = ic + 1
                   sendBuff(ic,isend) = b2(i,j)
                end do
             end do

             if(DEBUG_COMMS)then
                write(*,"(I3,': packed: ',6E15.4)") get_rank(), &
                     sendBuff(1:min(ic,6),isend)
             end if
             
             call post_send(sendBuff(1,isend),ic,destination(isend),tag, &
                            exch_flags(handle,isend,indexs))

          else if ( present(ib2) ) then

             call parallel_abort('exchange_generic: halo-swaps for 2D integer ' &
                             & //'fields not implemented.')
!!$             ic = 0
!!$             pack_patches2i: do ipatch=1, npatchsend(isend,1), 1
!!$                jstart = jsrcsendp(ipatch,isend,1)
!!$                istart = isrcsendp(ipatch,isend,1)
!!$                jend   = jstart+nysendp(ipatch,isend,1)-1
!!$                iend   = istart+nxsendp(ipatch,isend,1)-1
!!$
!!$                do j=jstart, jend, 1
!!$                   do i=istart, iend, 1
!!$                      ic = ic + 1
!!$                      sendIBuff(ic,isend) = ib2(i,j)
!!$                   end do
!!$                end do
!!$             end do pack_patches2i
!!$
!!$             CALL MPI_Isend(sendIBuff(1,isend),ic, MPI_INTEGER,   &
!!$                            destination(isend),tag,mpi_comm_world,&
!!$                            exch_flags(handle,isend,indexs),ierr)

          else if ( present(b3) )then

             call parallel_abort('exchange_generic: halo-swaps for 3D real ' &
                             & //'fields not implemented.')
!!$             ic = 0
!!$             pack_patches3r: do ipatch=1,npatchsend(isend,1)
!!$
!!$                jstart = jsrcsendp(ipatch,isend,1)
!!$                istart = isrcsendp(ipatch,isend,1)
!!$                jend   = jstart+nysendp(ipatch,isend,1)-1
!!$                iend   = istart+nxsendp(ipatch,isend,1)-1
!!$                do k=1,nzsendp(ipatch,isend,1),1
!!$                   do j=jstart, jend, 1
!!$                      do i=istart, iend, 1
!!$                         ic = ic + 1
!!$                         sendBuff(ic, isend) = b3(i,j,k)
!!$                      end do
!!$                   end do
!!$                end do
!!$             end do pack_patches3r


             !CALL MPI_Isend(sendBuff(1,isend),ic,                    &
             !               MPI_DOUBLE_PRECISION,                    &
             !               destination(isend), tag, MPI_COMM_WORLD, &
             !               exch_flags(handle,isend,indexs),ierr)

          end if

       else

          exch_flags(handle,isend,indexs) = MSG_REQUEST_NULL
       end if ! direction is enabled and have something to send

    enddo ! Loop over sends

    if(DEBUG_COMMS)then
       write (*,FMT="(I3,': exch tag ',I4,' finished all ',I4,' sends')") &
            get_rank(), tag, nsend
    end if

    ! Wait on the receives that were posted earlier

    ! Copy just the set of flags we're interested in for passing 
    ! to MPI_waitany
    exch_flags1d(1:nrecv) = exch_flags(handle, 1:nrecv, indexr)

    if(DEBUG_COMMS)then
       write(*,"(I3,': Doing waitany: nrecv =',I3,' handle = ',I3)") &
          get_rank(), nrecv,handle
    endif

    ! Get the first available message that we've received
    call msg_wait(nrecv, exch_flags1d, irecv)

    do while(irecv /= MSG_UNDEFINED)

       if ( present(b2) ) then

          ! Copy received data back into array
          ic = 0
          jstart = jdesrecv(irecv)
          jend   = jstart+nyrecv(irecv)-1

          istart = idesrecv(irecv)
          iend   = istart+nxrecv(irecv)-1

          if(DEBUG_COMMS)then
             write(*,"(I3,': unpacking to:',I3,':',I3,',',I3,':',I3)") &
                  get_rank(), istart, iend, jstart, jend
          end if
             
          do j=jstart, jend, 1
             do i=istart, iend, 1
                ic = ic + 1
                b2(i,j) = recvBuff(ic,irecv)
             end do
          end do
          
          if(DEBUG_COMMS)then
             write(*,"(I3,': unpacked: ',6E15.4)") get_rank(), &
                  recvBuff(1:min(6,ic),irecv)
          end if
             
       else if ( present(ib2) ) then

          ! Copy received data back into array
          ic = 0
          jstart = jdesrecv(irecv)
          jend   = jstart+nyrecv(irecv)-1
          istart = idesrecv(irecv)
          iend   = istart+nxrecv(irecv)-1
          do j=jstart, jend, 1
             do i=istart, iend, 1
                ic = ic + 1
                ib2(i,j) = recvIBuff(ic,irecv)
             end do
          end do

       else if (present(b3) ) then

          ic = 0
          
          jstart = jdesrecv(irecv)
          jend   = jstart+nyrecv(irecv)-1
          istart = idesrecv(irecv)
          iend   = istart+nxrecv(irecv)-1

          do k=1,nzrecv(irecv),1
             do j=jstart, jend, 1
                do i=istart, iend, 1
                   ic = ic + 1
                   b3(i,j,k) = recvBuff(ic,irecv)
                end do
             end do
          end do
             
       end if

       ! Wait for the next message
       call msg_wait(nrecv, exch_flags1d, irecv)
    end do ! while irecv != MSG_UNDEFINED

    if(DEBUG_COMMS)then
       write(*,"(I3,': Finished all ',I3,' receives for handle ',I3)") &
             get_rank(), nrecv, handle
    end if

    ! Copy just the set of flags we're interested in for passing to  
    ! MPI_waitall next time around  
    exch_flags1d(1:nsend) = exch_flags(handle, 1:nsend, indexs)

    ! Free the exchange communications handle.
    call free_exch_handle(handle)

    ! All receives done so we can safely free the MPI receive buffers
    if( allocated(recvBuff) ) deallocate(recvBuff)
    if( allocated(recvIBuff) )deallocate(recvIBuff)

  end subroutine exchange_generic

end module parallel_comms_mod
