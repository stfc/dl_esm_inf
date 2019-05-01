module parallel_comms_mod
  use kind_params_mod, only: go_wp
  use parallel_utils_mod, only: get_num_ranks, get_rank, parallel_abort, &
       MPI_UNDEFINED, MPI_REQUEST_NULL, msg_wait, get_max_tag, post_receive, &
       post_send
  use decomposition_mod, only: subdomain_type, decomposition_type
  implicit none

  private

  integer, parameter :: jpk = 1 !< Only 1 vertical level
  logical, parameter :: DEBUG = .true.
  !> En/Disable debug output on a per-message basis
  logical, parameter :: DEBUG_COMMS = .False.
  logical :: lwp !< Whether or not to write out from this rank

  ! Set by mapcomms
  integer :: halo_depthx, halo_depthy
  integer, parameter :: MAX_HALO_DEPTH = 1

  ! Information held by the ensemble leader about all processes.
  !     pnactive   Number of active points in each sub-domain.
  !     pielb      Lower (west) longitude bound index.
  !     pieub      Upper (east) longitude bound index.
  !     piesub     Number of longitude gridpoints.
  !     pilbext    True if the lower longitude boundary is external.
  !     piubext    True if the upper longitude boundary is external.
  !     pjelb      Lower (south) latitude bound index.
  !     pjeub      Upper (north) latitude bound index.
  !     pjesub     Number of latitude gridpoints.
  !     pjlbext    True if the lower latitude boundary is external.
  !     pjubext    True if the upper latitude boundary is external.

  !INTEGER, ALLOCATABLE, DIMENSION(:) :: pnactive,              &
  !     piesub, pjelb, pjeub, pjesub
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: pilbext, piubext, pjlbext, pjubext

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

  INTEGER, PARAMETER :: MaxComm=16
  INTEGER, SAVE, DIMENSION(MaxComm)         :: dirsend,destination, &
                                               dirrecv,source
!  INTEGER, SAVE, DIMENSION(MaxComm,halo_depthx ) :: isrcsend,jsrcsend &
  INTEGER, SAVE, DIMENSION(MaxComm) :: isrcsend,jsrcsend, &
                                       isrcrecv,jsrcrecv, &
                                       idessend,jdessend, &
                                       nxsend,nysend,nzsend, &
                                       idesrecv,jdesrecv, &
                                       nxrecv,nyrecv,nzrecv
  INTEGER, SAVE :: nsend,nrecv

  ! Total number of points in each message
  INTEGER, SAVE, DIMENSION(MaxComm,MAX_HALO_DEPTH) :: nsendp, nsendp2d, &
                                                   nrecvp, nrecvp2d

  ! Process dependent partitioning information.

  !     ielb       Lower (west) longitude bound index.
  !     ieub       Upper (east) longitude bound index.
  !     iesub      Number of longitude gridpoints.
  !     ilbext     True if the lower longitude boundary is external.
  !     iubext     True if the upper longitude boundary is external.
  !     jelb       Lower (south) latitude bound index.
  !     jeub       Upper (north) latitude bound index.
  !     jesub      Number of latitude gridpoints.
  !     jlbext     True if the lower latitude boundary is external.
  !     jubext     True if the upper latitude boundary is external.

  INTEGER, SAVE :: ielb, ieub, iesub
  INTEGER, SAVE :: jelb, jeub, jesub
  LOGICAL, SAVE :: ilbext, iubext, jlbext, jubext

  ! Global definitions for parallel execution.

  ! Direction flags for communications.
  ! Listed so that opposite directions are given values maximally spaced
  INTEGER, PARAMETER :: NONE=0         &
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
  INTEGER, DIMENSION(MaxCommDir) :: opp_dirn

  ! Set up indices indicating the north-south and east-west
  ! attributes of the eight basic communication directions:
  ! four edges: W, E, S, N;
  ! four corners: SW, SE, NW, NE.
  INTEGER, PARAMETER, DIMENSION(MaxCommDir) ::       &
                west  = (/ 1, 0, 0, 0, 1, 0, 1, 0 /) &
               ,east  = (/ 0, 1, 0, 0, 0, 1, 0, 1 /) &
               ,south = (/ 0, 0, 1, 0, 1, 0, 0, 1 /) &
               ,north = (/ 0, 0, 0, 1, 0, 1, 1, 0 /)
                         ! 1  2  3  4  5  6  7  8
                         ! W  E  S  N SW NE NW  SE

  ! Value representing land in the mask used for partitioning and message
  ! trimming
  INTEGER, PARAMETER :: LAND = 0

  ! Rather than trim to a point immediately next to a wet point, we 
  ! back off nextra points. If we don't do this then the sea-ice
  ! computation goes wrong because it does use values over the land
  ! that immediately border the ocean.
  INTEGER, SAVE :: nextra
  
  ! exch_flags    -  Array of flag arrays for exchanges
  ! exch_flags1d  -  Array of only the current MPI receive operations
  ! exch_tag      -  The tag value associated with this exchange
  ! exch_busy     -  Indicates whether a slot in the flag array is being used

  INTEGER, PARAMETER :: indexs=1,indexr=2
  INTEGER, PARAMETER :: max_flags=40
  INTEGER, PARAMETER :: min_tag=0
  INTEGER, SAVE :: current_tag,max_tag_used,max_tag,n_tag_cycles=0
  LOGICAL, SAVE :: first_mod=.TRUE.

  INTEGER, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: exch_flags
  INTEGER, ALLOCATABLE, DIMENSION(:),     SAVE :: exch_tag, exch_flags1d
  LOGICAL, ALLOCATABLE, DIMENSION(:),     SAVE :: exch_busy

  TYPE exch_item
     INTEGER               :: halo_width
     INTEGER, DIMENSION(4) :: dirn
     INTEGER               :: isgn
     CHARACTER(LEN=1)      :: grid
     LOGICAL               :: lfill
     integer,  dimension(:,:),   pointer :: i2dptr
     integer,  dimension(:,:,:), pointer :: i3dptr
     real(go_wp), dimension(:,:),   pointer :: r2dptr
     real(go_wp), dimension(:,:,:), pointer :: r3dptr
  END TYPE exch_item

  TYPE (exch_item), ALLOCATABLE, SAVE :: exch_list(:)
  INTEGER, SAVE :: nextFreeExchItem, maxExchItems
  
  ! Buffer for doing halo-exchange.
  ! For a 3D array, halos are 2D slabs but copied into these buffers 
  ! as 1D vectors. 2nd dimension refers to the direction of the 
  ! communication.
  ! For a 2D array, halos are 1D vectors anyway.
  real(go_wp), dimension(:,:), allocatable, save :: sendBuff, recvBuff
  integer , dimension(:,:), allocatable, save :: sendIBuff, recvIBuff

  ! Public routines
  PUBLIC :: map_comms, iprocmap, exchmod_alloc, exchs_generic

  ! Public variables
  PUBLIC :: MaxComm,nsend,nrecv,nxsend,nysend,destination,dirrecv, &
            dirsend,isrcsend,jsrcsend,idesrecv, jdesrecv,          &
            nxrecv, nyrecv, source, idessend, jdessend

  PUBLIC :: nsendp,nsendp2d,nrecvp,nrecvp2d

  PUBLIC :: ielb,  ieub,                      &
            jeub, ilbext, iubext, jubext, jlbext, &
            jelb, pilbext, pjlbext, pjubext, piubext

  PUBLIC :: NONE         &
           ,Iplus        & 
           ,Iminus       &
           ,Jplus        &
           ,Jminus       &
           ,IplusJplus   &
           ,IminusJminus &
           ,IplusJminus  &
           ,IminusJplus  &
           ,MaxCommDir

  PUBLIC :: opp_dirn, land

  PUBLIC :: nextra

  ! Switch for trimming dry points from halo swaps
  LOGICAL, PARAMETER :: msgtrim   = .TRUE.

  ! Switch for trimming points below ocean floor from halo swaps
  ! Defaults to true unless set via NEMO_MSGTRIM_Z environment var.
  LOGICAL, PUBLIC, SAVE :: msgtrim_z

contains

  subroutine map_comms (decomp, tmask, pbc, halo_depths, ierr)
    !!------------------------------------------------------------------
    ! Maps out the communications requirements for the partitioned
    ! domain, adding communications descriptions to the list.
    !
    !     Mike Ashworth, CLRC Daresbury Laboratory, July 1999
    !     Andy Porter, STFC Daresbury Laboratory, March 2019
    !!------------------------------------------------------------------
    IMPLICIT NONE

    ! Subroutine arguments.
    type(decomposition_type), target, intent(in) :: decomp
    integer, intent(in), allocatable :: tmask(:,:)  ! decomposed mask: 0 for land, 1 for ocean
    logical, intent(in) :: pbc ! Whether mesh is periodic in x dimension
    integer, intent(in) :: halo_depths(2)
    integer, intent(out):: ierr

    ! Local variables.
    INTEGER :: idirn, i1, i2, ihalo, iproc, iprocx, &
         iprocy, j, j1, j2, nadd, naddmaxr, naddmaxs, &
         iinside, ioutside
    INTEGER :: ldiff0, ldiff1 ! Local vars for coping with wrapping of coords
    INTEGER :: imax, imin ! Max/min value of i that a halo strip can run 
                          ! to/from (used to avoid including E/W halo cells
                          ! in N/S transfers)
    INTEGER :: ielb_iproc, ieub_iproc ! Lower and upper bounds for proc
                                      ! iproc corrected for E & W halos
    INTEGER :: jelb_iproc, jeub_iproc
    INTEGER, DIMENSION(MAX_HALO_DEPTH) :: idesr, jdesr, idess, jdess, &
                                isrcr, jsrcr, isrcs, jsrcs, &
                                nxr, nyr, nxs, nys
    logical :: addcorner
    integer :: nprocp, irank
    type(subdomain_type), pointer :: subdomain

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
    DO WHILE (j1.LE.jeub)

       ! Look for the process which owns the neighbouring point in the
       ! minus I direction.

       iproc = iprocmap(decomp, ielb-1, j1)
       IF ( iproc.GT.0 ) THEN

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

          DO ihalo=1,halo_depthx
            ! Source for the receive must be within internal domain of the
            ! (sending) PE, iproc
            isrcr(ihalo) = decomp%subdomains(iproc)%internal%xstop - ihalo + 1
            idesr(ihalo) = ihalo ! Halo goes from 1..halo_depthx
            nxr(ihalo) = ihalo
            nxs(ihalo) = ihalo
          ENDDO

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
          DO ihalo=1,halo_depthx,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.

             IF ( j1-ihalo.GE.jelb .AND.  &
                  iprocmap(decomp,ielb-ihalo,j1).GT.0 ) THEN
                naddmaxs = ihalo
             ENDIF

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            IF ( j1.EQ.jelb .AND.  &
                 iprocmap(decomp,ielb-ihalo,j1-ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO

          ! Add the extra points.

          DO ihalo=1,halo_depthx
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
          DO ihalo=1,halo_depthx,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.

            IF ( j2+ihalo.LE.jeub .AND. &
                 iprocmap(decomp,ielb-ihalo,j2).GT.0 ) THEN
              naddmaxs = ihalo
            ENDIF

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            IF ( j2.EQ.jeub .AND. &
                 iprocmap(decomp,ielb-ihalo,j2+ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO 

          ! Add the extra points.

          DO ihalo=1,halo_depthx,1
            nadd = MIN(ihalo,naddmaxs)
            nys(ihalo) = nys(ihalo)+nadd
            if(DEBUG)then
               IF ( nadd.GT.0 ) THEN
                  WRITE (*,"(I3,': Adding starting corner to send for halo '," &
                       //"I2,' with ',I3,' points')") irank,ihalo, nadd
               ENDIF
            end if
            nadd = MIN(ihalo,naddmaxr)
            nyr(ihalo) = nyr(ihalo)+nadd
            if(DEBUG)then
               IF ( nadd.GT.0  ) THEN
                  WRITE (*,"(I3,': Adding starting corner to receive for halo " &
                       //"',I2,' with ',I3,' points')") irank,ihalo, nadd
               ENDIF
            end if
          ENDDO

          ! Add a send and a receive to the lists for this section
          ! of border.
          CALL addsend (nsend,Iplus,iproc-1, &
                        isrcs,jsrcs,idess,jdess,   &
                        nxs,nys,tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          ! This is a receive for data _from_ the neighbouring domain
          ! - it is NOT the corresponding receive for the above send.
          CALL addrecv (nrecv,Iminus,iproc-1, &
                        isrcr,jsrcr,idesr,jdesr,    &
                        nxr,nyr,tmask,ierr)
          IF ( ierr.NE.0 ) RETURN
          if(DEBUG)then
             IF ( lwp ) THEN
                WRITE (*,'(a,7i6)') 'Adding receive ',iproc-1 &
                   ,isrcr(1),jsrcr(1),idesr(1),jdesr(1),nxr(1),nyr(1)
             ENDIF
          end if

          ! Move the start point to one beyond this strip.
          j1 = j2+1

        ELSE

           ! No process found, continue searching up the boundary.
          j1 = j1+1
        ENDIF
      ENDDO

      ! =================================================================
      ! Looking at the border where we will
      ! send data in the plus I direction (Iminus) and
      ! receive data that has been sent in the minus I direction (Iplus).
      ! =================================================================

      ! Start from the lower bound of the sub-domain, and carry on looking
      ! for communications with neighbours until we have reached
      ! the upper bound.

      j1 = jelb
      DO WHILE (j1.LE.jeub)

         ! Look for the process which owns the neighbouring point in the
         ! plus I direction.

        iproc = iprocmap(decomp,ieub+1,j1)
        IF ( iproc.GT.0 ) THEN

           ! Find where in the j direction the common border between these
           ! sub-domains ends.
          j2 = MIN(jeub,decomp%subdomains(iproc)%global%ystop)

          if(DEBUG)then
             WRITE (*,FMT="(I3,': ARPDBG strip for plus I is ',I3,',',I3)") &
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
          DO ihalo=1,halo_depthx
             isrcr(ihalo) = decomp%subdomains(iproc)%internal%xstart + ihalo
             ! Outermost halo -> innermost col.
             isrcs(ihalo) = subdomain%internal%xstop - ihalo + 1
             idess(ihalo) = ihalo ! Halo runs from 1..halo_depthx
             nxr(ihalo) = ihalo
             nxs(ihalo) = ihalo
             ! Destination for a receive is in local halo region
             idesr(ihalo) = subdomain%internal%xstop + ihalo
          ENDDO
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
          DO ihalo=1,halo_depthx,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.
             IF ( j1-ihalo.GE.jelb .AND.  &
                 iprocmap(decomp,ieub+ihalo,j1).GT.0 ) THEN
                naddmaxs = ihalo
             ENDIF

             ! Receive corner data only when we are at the corner
             ! and while the sending sub-domain is the same as for the edge.

             IF ( j1.EQ.jelb .AND.  &
                 iprocmap(decomp,ieub+ihalo,j1-ihalo).EQ.iproc ) THEN
                naddmaxr = ihalo
             ENDIF
          ENDDO

          ! Add the extra points.
          DO ihalo=1,halo_depthx,1
            nadd = MIN(ihalo,naddmaxs)
            jdess(ihalo) = jdess(ihalo) - nadd
            jsrcs(ihalo) = jsrcs(ihalo) - nadd
            nys(ihalo) = nys(ihalo)+nadd
            nadd = MIN(ihalo,naddmaxr)
            jdesr(ihalo) = jdesr(ihalo) - nadd
            jsrcr(ihalo) = jsrcr(ihalo) - nadd
            nyr(ihalo) = nyr(ihalo)+nadd
          ENDDO

          ! Examine whether corner points should be added to the end.
          naddmaxr = 0
          naddmaxs = 0
          DO ihalo=1,halo_depthx,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.

            IF ( j2+ihalo.LE.jeub .AND. &
                 iprocmap(decomp,ieub+ihalo,j2).GT.0 ) THEN
              naddmaxs = ihalo
            ENDIF

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            IF ( j2.EQ.jeub .AND. &
                 iprocmap(decomp,ieub+ihalo,j2+ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO 

          ! Add the extra points.

          DO ihalo=1,halo_depthx,1
            nadd = MIN(ihalo,naddmaxs)
            nys(ihalo) = nys(ihalo)+nadd
            if(DEBUG)then
               IF ( nadd.GT.0 .AND. lwp ) THEN
                  WRITE (*,*) 'Adding starting corner to send' &
                          ,' for halo ',ihalo,' with ',nadd,' points'
               ENDIF
            end if
            nadd = MIN(ihalo,naddmaxr)
            nyr(ihalo) = nyr(ihalo)+nadd
            if(DEBUG)then
               IF ( nadd.GT.0 .AND. lwp ) THEN
                  WRITE (*,*) 'Adding starting corner to receive' &
                  ,' for halo ',ihalo,' with ',nadd,' points'
               ENDIF
            end if
          ENDDO

          ! Add a send and a receive to the lists for this section
          ! of border.
          CALL addsend (nsend,Iminus,iproc-1,            &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          CALL addrecv (nrecv,Iplus,iproc-1,       &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          ! Move the start point to one beyond this strip.
          j1 = j2+1
        ELSE

          ! No process found, continue searching up the boundary.
          j1 = j1+1
        ENDIF
      ENDDO

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
      DO WHILE (i1.LE.imax)

         ! Look for the process which owns the neighbouring point in the
         ! minus J direction.

         iproc = iprocmap(decomp,i1,jelb-1)
         IF ( iproc.GT.0 ) THEN

            ! Ensure we don't include halos from the global domain borders if
            ! we have cyclic bc's.
            ieub_iproc = decomp%subdomains(iproc)%global%xstop

            ! Find where in the i direction the common border between these
            ! sub-domains ends (in the global domain).
            i2 = MIN(imax,ieub_iproc)

            if(DEBUG)then
               WRITE (*,FMT="(I3,': ARPDBG strip for minus J is ',I3,',',I3)") &
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
            DO ihalo=1,halo_depthy,1
               ! Source for a receive must be within internal region
               jsrcr(ihalo) = decomp%subdomains(iproc)%internal%ystop-ihalo+1
               jdesr(ihalo) = ihalo ! Halo runs from 1..halo_depthy
               nyr(ihalo) = ihalo
               nys(ihalo) = ihalo
            ENDDO
            ! Destination for a send must be a halo region
            jdess(:) = decomp%subdomains(iproc)%internal%ystop + 1

            ! Examine whether corner points should be added to the START.
            naddmaxr = 0
            naddmaxs = 0
            DO ihalo=1,halo_depthy,1

               ! Send corner data while we have data to send
               ! and while there is a point that depends on it.
               IF ( i1-ihalo.GE.imin .AND.  &
                    iprocmap(decomp,i1-ihalo,jelb-ihalo).GT.0 ) THEN
                  naddmaxs = ihalo
               ENDIF

               ! Receive corner data only when we are at the corner and
               ! while the sending sub-domain is the same as for the edge.
               IF ( i1.EQ.imin .AND. &
                    iprocmap(decomp,i1-ihalo,jelb-ihalo).EQ.iproc ) THEN
                  naddmaxr = ihalo
               ENDIF
            ENDDO

            ! Add the extra points.

            DO ihalo=1,halo_depthy,1
               nadd = MIN(ihalo,naddmaxs)
               idess(ihalo) = idess(ihalo)-nadd
               isrcs(ihalo) = isrcs(ihalo)-nadd
               nxs(ihalo) = nxs(ihalo)+nadd
               if(DEBUG)then
                  IF ( nadd.GT.0 ) THEN
                     WRITE (*,"(I3,': Adding starting corner to send for &
                       & halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
                  ENDIF
               end if
               nadd = MIN(ihalo,naddmaxr)
               idesr(ihalo) = idesr(ihalo)-nadd
               isrcr(ihalo) = isrcr(ihalo)-nadd
               nxr(ihalo) = nxr(ihalo)+nadd

               if(DEBUG)then
                  IF ( nadd.GT.0 ) THEN
                     WRITE (*,"(I3,': Adding starting corner to receive for &
                       & halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
                  ENDIF
               end if
            ENDDO

            ! Examine whether corner points should be added to the END.

          naddmaxr = 0
          naddmaxs = 0
          DO ihalo=1,halo_depthy,1

             ! Send corner data while we have data to send
             ! and while there is a point that depends on it.

            IF ( i2+ihalo.LE.imax .AND. &
                 iprocmap(decomp,i2,jelb-ihalo).GT.0 ) THEN
              naddmaxs = ihalo
            ENDIF

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            IF ( i2.EQ.imax .AND. & 
                 iprocmap(decomp,i2+ihalo,jelb-ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO 

          ! Add the extra points.

          DO ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            nxs(ihalo) = nxs(ihalo)+nadd
            if(DEBUG)then
               IF ( nadd.GT.0 ) THEN
                  WRITE (*,"(I3,': Adding starting corner to send for &
                       & halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
               ENDIF
            end if
            nadd = MIN(ihalo,naddmaxr)
            nxr(ihalo) = nxr(ihalo)+nadd
            if(DEBUG)then
               IF ( nadd.GT.0 ) THEN
                  WRITE (*,"(I3,': Adding starting corner to receive for & 
                       & halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
               ENDIF
            end if
          ENDDO

          ! Add a send and a receive to the lists for this section
          ! of border.

          CALL addsend (nsend,Jplus,iproc-1,       &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          CALL addrecv (nrecv,Jminus,iproc-1,      &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          ! Move the start point to one beyond this strip.
          i1 = i2+1

        ELSE

           ! No process found, continue searching up the boundary.
          i1 = i1+1
        ENDIF
      ENDDO

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

      DO WHILE (i1.LE.imax)

         ! Look for the process which owns the neighbouring point in the
         ! plus J direction.

         iproc = iprocmap(decomp,i1,jeub+1)
         IF ( iproc.GT.0 ) THEN
           ielb_iproc = decomp%subdomains(iproc)%global%xstart
           ieub_iproc = decomp%subdomains(iproc)%global%xstop

           ! Find where in the i direction the common border between these
           ! sub-domains ends (in the global domain).
           i2 = MIN(imax, ieub_iproc)

           if(DEBUG)then
              WRITE (*,FMT="(I3,': ARPDBG strip for plus J is ',I3,',',I3)") &
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

            DO ihalo=1,halo_depthy,1
               ! innermost row -> outermost halo
               jsrcs(ihalo) = subdomain%internal%ystop-ihalo+1
               jdess(ihalo) = ihalo        ! Halo runs from 1..halo_depthy
               nyr(ihalo)   = ihalo
               nys(ihalo)   = ihalo
            ENDDO
            ! Destination for receive must be in halo
            jdesr(:) = subdomain%internal%ystop + 1

          ! Examine whether corner points should be added to the START.

          naddmaxr = 0
          naddmaxs = 0
          DO ihalo=1,halo_depthy,1

            ! Send corner data while we have data to send
            ! and while there is a point that depends on it.
            if(i1-ihalo >= imin .and. &
               iprocmap(decomp, i1, jeub+ihalo) > 0)then
              naddmaxs = ihalo
            endif

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.
            IF ( i1.EQ.imin .AND. &
                 iprocmap(decomp,i1-ihalo,jeub+ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO

          ! Add the extra points.

          DO ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            idess(ihalo) = idess(ihalo) -nadd
            isrcs(ihalo) = isrcs(ihalo) -nadd
            nxs(ihalo) = nxs(ihalo)+nadd
            if(DEBUG)then
               IF ( nadd.GT.0 ) THEN
                  WRITE (*,"(I3,': Adding starting corner to send for halo '," &
                  //"I2,' with ',I3,' points')") irank-1,ihalo, nadd
               ENDIF
            end if
            nadd = MIN(ihalo,naddmaxr)
            idesr(ihalo) = idesr(ihalo) - nadd
            isrcr(ihalo) = isrcr(ihalo) - nadd
            nxr(ihalo) = nxr(ihalo)+nadd

            if(DEBUG)then
               IF ( nadd.GT.0 ) THEN
                  WRITE (*,"(I3,': Adding starting corner to receive for " &
                       //"halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
               ENDIF
            end if
          ENDDO

          ! Examine whether corner points should be added to the END.
          naddmaxr = 0
          naddmaxs = 0
          DO ihalo=1,halo_depthy,1

            ! Send corner data while we have data to send
            ! and while there is a point that depends on it.
            IF ( i2+ihalo.LE.imax .AND. &
                 iprocmap(decomp,i2,jeub+ihalo).GT.0 ) THEN
              naddmaxs = ihalo
            ENDIF

            ! Receive corner data only when we are at the corner
            ! and while the sending sub-domain is the same as for the edge.

            IF ( i2.EQ.imax .AND.       & ! Are we at the corner?
                 iprocmap(decomp,i2+ihalo,jeub+ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO 

          ! Add the extra points.

          DO ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            nxs(ihalo) = nxs(ihalo)+nadd

            if(DEBUG)then
               IF ( nadd.GT.0 .AND. lwp ) THEN
                  WRITE (*,*) irank-1,': Adding starting corner to send' &
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
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          CALL addrecv (nrecv,Jplus,iproc-1,       &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          ! Move the start point to one beyond this strip.
          i1 = i2+1

       ELSE ! iproc < 0

           ! No process found, continue searching up the boundary.
          i1 = i1+1
        ENDIF
      ENDDO

      ! =================================================================
      ! Corner points are sent with the edge data where possible.
      ! Check for cases where corner data resides on a processor which
      ! is not participating in edge communications and set up extra
      ! diagonal communications for these cases.
      ! =================================================================

      ! Loop over the four corner directions
      ! i = 1  2  3  4  5  6  7 8
      !     W  E  S  N SW NE NW SE
      DO idirn = 5, 8

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

        DO ihalo=halo_depthx,1,-1

          ! Look at the processor in the corner at width ihalo from the
          ! corner of the sub-domain. We want both x and y to be just
          ! outside our domain. ioutside is already just outside our domain
          ! so we subtract one from ihalo below:
          iproc = iprocmap(decomp, &
                    ioutside - west(idirn)*(ihalo-1)+ east(idirn)*(ihalo-1), &
                    south(idirn)*(jelb-ihalo)+north(idirn)*(jeub+ihalo))

          ! If the corner processor is different from those to X and Y
          ! we will need a corner communication.
          IF ( iproc.GT.0 .AND. iprocx.GT.0 .AND. iprocy.GT.0 .AND. &
               iproc.NE.iprocx .AND. iproc.NE.iprocy ) THEN

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
             ENDIF

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

          ELSE

             if(DEBUG)then
                WRITE (*,FMT="(I3,': skipping corner for halo ',I3,' PE ',I3)")&
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
          ENDIF

       ENDDO ! Loop over halo depths

       ! The size of the received data is always the same as 
       ! that of the sent data.
       nxr(:) = nxs(:)
       nyr(:) = nys(:)

        ! Add the data to the communications lists.
        IF ( addcorner ) THEN
          if(DEBUG)then
             write (*,FMT="(I3,': ARPDBG adding corner send to ',I4, &
                  & ', dir = ',I1)") irank-1, iproc-1,idirn
          end if
          call addsend (nsend,idirn,iproc-1,            &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        tmask, ierr)
          if ( ierr.NE.0 ) then
            call parallel_abort("Call to addsend failed.")
          end if

          ! Manually reverse the direction indicator for the receive.
          j = opp_dirn(idirn)

          if(DEBUG)then
             WRITE (*,FMT="(I3,': ARPDBG adding corner recv. from ',I4,', &
                & old dir = ',I1,' new dir = ',I1)") irank-1, iproc-1,idirn, j
          end if
          call addrecv (nrecv,j,iproc-1,          &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        tmask, ierr)
          if ( ierr.NE.0 ) then
             call parallel_abort("Call to addrecv failed")
          end if
        ENDIF

      ENDDO

    END SUBROUTINE map_comms

  subroutine addsend (icomm, dir, proc, isrc, jsrc, &
                      ides, jdes, nx, ny, tmask, ierr )
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
    INTEGER,                    INTENT(inout) :: icomm
    INTEGER,                    INTENT( in  ) :: dir, proc
    INTEGER, DIMENSION(:,:),    INTENT( in  ) :: tmask
    INTEGER,                    INTENT( out ) :: ierr
    INTEGER, DIMENSION(:), INTENT( in  ) :: isrc, jsrc, ides, jdes, nx, ny
    INTEGER :: irank
         
    irank = get_rank()
    ! Clear the error flag.
    ierr = 0

    ! Return if the process id is not set.
    IF ( proc.LT.0 ) THEN
       RETURN
    ENDIF

    icomm = icomm+1

    ! Check that the comms list has space for another entry.
    IF ( icomm.GT.MaxComm ) THEN
       IF ( lwp ) THEN
          WRITE (*,*) 'ERROR: Number of separate ', &
               'communications exceeds maximum of ',MaxComm
       ENDIF
       ierr = -12
       RETURN
    ENDIF
            
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
       WRITE(*,FMT="(I4,': ARPDBG adding SEND:')") irank-1  
       WRITE(*,FMT="(I4,': ARPDBG: icomm = ',I2)") irank-1,icomm
       WRITE(*,FMT="(I4,': ARPDBG:   dir = ',I2)") irank-1,dirsend(icomm)
       WRITE(*,FMT="(I4,': ARPDBG:  proc = ',I4)") irank-1,destination(icomm)
       WRITE(*,FMT="(I4,': ARPDBG:  isrc = ',I4)") irank-1,isrcsend(icomm)
       WRITE(*,FMT="(I4,': ARPDBG:  jsrc = ',I4)") irank-1,jsrcsend(icomm)
       WRITE(*,FMT="(I4,': ARPDBG:  ides = ',I4)") irank-1,idessend(icomm)
       WRITE(*,FMT="(I4,': ARPDBG:  jdes = ',I4)") irank-1,jdessend(icomm)
       WRITE(*,FMT="(I4,': ARPDBG:    nx = ',I4)") irank-1,nxsend(icomm)
       WRITE(*,FMT="(I4,': ARPDBG:    ny = ',I4)") irank-1,nysend(icomm)
       WRITE(*,FMT="(I4,': ARPDBG:    nz = ',I4)") irank-1,nzsend(icomm) 
       WRITE (*,FMT="(I4,': ARPDBG:nsendp = ',I4)") irank-1,nsendp(icomm,1)
       WRITE (*,FMT="(I4,': ARPDBG SEND ends')")    irank-1
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
                       ides, jdes, nx, ny, tmask, ierr )
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
    INTEGER,                 INTENT(inout) :: icomm
    INTEGER,                 INTENT( in  ) :: dir, proc
    INTEGER,                 INTENT(out)   :: ierr
    INTEGER, DIMENSION(:,:), INTENT( in  ) :: tmask
    INTEGER, DIMENSION(halo_depthx)        :: isrc, jsrc, ides, jdes, nx, ny
    ! Local variables.
    INTEGER :: irank

    irank = get_rank()
    ! Clear the error flag.
    ierr = 0

    ! Return if the process id is not set.
    IF ( proc.LT.0 ) THEN
       RETURN
    ENDIF

    icomm = icomm + 1

    ! Check that the comms list has space for another entry.

    IF ( icomm.GT.MaxComm ) THEN
       IF ( lwp ) THEN
          WRITE (*,*) 'ERROR: Number of separate ', &
               'communications exceeds maximum of ',MaxComm
       ENDIF
       ierr = -12
       RETURN
    ENDIF

    ! Add the data into the comms list at the new location.
    dirrecv(icomm)  = dir
    source(icomm)   = proc
    isrcrecv(icomm) = isrc(1)
    jsrcrecv(icomm) = jsrc(1)
    idesrecv(icomm) = ides(1)
    jdesrecv(icomm) = jdes(1)
    nxrecv(icomm)   = nx(1)
    nyrecv(icomm)   = ny(1)
    !> TODO handle multiple halo depths
    nrecvp2d(icomm, 1) = nxrecv(icomm)*nyrecv(icomm)
    nzrecv(icomm)   = 1
    nrecvp(icomm, 1) = nzrecv(icomm) * nrecvp2d(icomm, 1)
    
    if(DEBUG)then
       WRITE (*,FMT="(I4,': ARPDBG adding RECV:')") irank-1
       WRITE (*,FMT="(I4,': ARPDBG: icomm = ',I2)") irank-1,icomm
       WRITE (*,FMT="(I4,': ARPDBG:   dir = ',I2)") irank-1,dir
       WRITE (*,FMT="(I4,': ARPDBG:  proc = ',I4)") irank-1,proc
       WRITE (*,FMT="(I4,': ARPDBG:  isrc = ',I4)") irank-1,isrcrecv(icomm)
       WRITE (*,FMT="(I4,': ARPDBG:  jsrc = ',I4)") irank-1,jsrcrecv(icomm)
       WRITE (*,FMT="(I4,': ARPDBG:  ides = ',I4)") irank-1,idesrecv(icomm)
       WRITE (*,FMT="(I4,': ARPDBG:  jdes = ',I4)") irank-1,jdesrecv(icomm)
       WRITE (*,FMT="(I4,': ARPDBG:    nx = ',I4)") irank-1,nxrecv(icomm)
       WRITE (*,FMT="(I4,': ARPDBG:    ny = ',I4)") irank-1,nyrecv(icomm)
       WRITE (*,FMT="(I4,': ARPDBG:nrecvp = ',I4)") irank-1,nrecvp(icomm,1)
       WRITE (*,FMT="(I4,': ARPDBG:nrecvp2d = ',I4)") irank-1,nrecvp2d(icomm,1)
       WRITE (*,FMT="(I4,': ARPDBG RECV ends')")    irank-1
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
  
  INTEGER FUNCTION exchmod_alloc()
    IMPLICIT none
    ! Locals
    INTEGER :: ierr, ii
    ! Since halos are broken up into wet-point-only patches we
    ! allocate the send and receive buffers  on a per-PE basis once we
    ! know the sizes of the patches (in exchs_generic).
    maxExchItems = 20
    ALLOCATE(exch_list(maxExchItems),           &
             exch_flags(max_flags,MaxComm,2),   &
             exch_flags1d(MaxComm),             &
             exch_busy(max_flags),              &
             exch_tag(max_flags),               &
             STAT=ierr)

    IF(ierr .eq. 0)THEN

       DO ii=1,maxExchItems,1
          NULLIFY(exch_list(ii)%r2dptr, exch_list(ii)%r3dptr, &
                  exch_list(ii)%i2dptr, exch_list(ii)%i3dptr)
       END DO

       exch_busy   = .FALSE.
    ELSE
       maxExchItems = 0
    END IF

    nextFreeExchItem = 1

    ! Pass back the allocation status flag
    exchmod_alloc = ierr

  END FUNCTION exchmod_alloc

  INTEGER FUNCTION get_exch_handle ( )
    ! ---------------------------------------------------------------
    ! Gets a new exchange handle
    ! ---------------------------------------------------------------
    IMPLICIT NONE
    ! Local variables.
    INTEGER :: h

    IF ( first_mod ) THEN

       ! First time in the module (i.e. exch or glob) set up the tags.
       max_tag = get_max_tag()
       if(DEBUG)then
          IF ( lwp ) WRITE (*,*) 'MAX_TAG: set to ',max_tag
       end if

       ! Set the current tag to the minimum.
       current_tag = min_tag
       max_tag_used = current_tag

       first_mod = .FALSE.
    ENDIF

    ! Look for a free location in the flags array.

    flag_search : DO h=1,max_flags
       IF ( .NOT.exch_busy(h) ) EXIT flag_search
    ENDDO flag_search

    IF ( h.GT.max_flags ) THEN

       ! If no free flags array was found, flag an error.
       call parallel_abort('ERROR: get_exch_handle: no free flag array')
    ELSE

       ! Assign a new tag.
       exch_busy(h) = .TRUE.

       IF ( current_tag.GE.(max_tag-MaxCommDir) ) THEN

          ! Wrap around.
          current_tag = min_tag
          n_tag_cycles = n_tag_cycles+1
       ELSE
          current_tag = current_tag + MaxCommDir
       ENDIF
       max_tag_used = MAX(max_tag_used,current_tag)
       exch_tag(h) = current_tag

    ENDIF

    get_exch_handle = h

  END FUNCTION get_exch_handle

  !================================================================

  SUBROUTINE free_exch_handle ( h )
    ! Frees exchange handle, h.
    IMPLICIT NONE

    ! Subroutine arguments.
    INTEGER :: h ! Handle to be free'd

    ! Free the flags array.
    IF ( h.GT.0 .AND. h.LE.max_flags ) THEN
       exch_busy(h) = .FALSE.
    ELSE
       WRITE (*,*) 'free_exch_handle: invalid handle ',h
    ENDIF

  END SUBROUTINE free_exch_handle

  !================================================================

  SUBROUTINE exchs_generic ( b2, ib2, b3, ib3, &
                             handle, comm1, comm2, comm3, comm4)

    ! *******************************************************************
    ! Send boundary data elements to adjacent sub-domains.

    ! b2(:,:)                real   input       2D real*8 local array.
    ! ib2(:,:)               int    input       2D integer local array.
    ! b3(:,:,:)              real   input       3D real*8 local array.
    ! ib3(:,:,:)             int    input       3D integer local array.
    ! rows/cols to exchange.
    ! handle                 int    output      Exchange handle.
    ! comm1                  int    input       Send in direction comm1.
    ! comm2                  int    input       Send in direction comm2.
    ! comm3                  int    input       Send in direction comm3.
    ! comm4                  int    input       Send in direction comm4.
    !
    ! Mike Ashworth, CCLRC, March 2005.
    ! Andrew Porter, STFC,  January 2008
    ! *******************************************************************
    use parallel_utils_mod, only: DIST_MEM_ENABLED
    IMPLICIT none

    ! Subroutine arguments.
    INTEGER, INTENT(out) :: handle
    REAL(go_wp),OPTIONAL, INTENT(inout), DIMENSION(:,:)   :: b2
    INTEGER, OPTIONAL, INTENT(inout), DIMENSION(:,:)   :: ib2
    REAL(go_wp),OPTIONAL, INTENT(inout), DIMENSION(:,:,:) :: b3
    INTEGER, OPTIONAL, INTENT(inout), DIMENSION(:,:,:) :: ib3

    INTEGER,           INTENT(in) :: comm1, comm2, comm3, comm4

    ! Local variables.

    LOGICAL :: enabled(0:MaxCommDir)
    INTEGER :: ierr, irecv, isend, &
               isrc, jsrc, nxs, nys, tag, tag_orig
    INTEGER :: maxrecvpts, maxsendpts ! Max no. of grid points involved in 
                                      ! any one halo exchange
    INTEGER :: i, j, k, ic ! Loop counters
    INTEGER :: istart, iend, jstart, jend
    LOGICAL, SAVE :: first_time = .TRUE.
    INTEGER, PARAMETER :: index_z = 3
    !!--------------------------------------------------------------------

    ! Do nothing if distributed-memory support is not enabled
    if(.not. DIST_MEM_ENABLED)return
    
    ierr = 0

    ! Allocate a communications tag/handle and a flags array.

    handle   = get_exch_handle()
    tag_orig = exch_tag(handle)

    ! Set enabled flags according to the subroutine arguments.

    enabled(Iplus ) = .FALSE.
    enabled(Jplus ) = .FALSE.
    enabled(Iminus) = .FALSE.
    enabled(Jminus) = .FALSE.
    enabled(comm1) = comm1.GT.0
    enabled(comm2) = comm2.GT.0
    enabled(comm3) = comm3.GT.0
    enabled(comm4) = comm4.GT.0

    ! Set diagonal communications according to the non-diagonal flags.

    enabled(IplusJplus ) = enabled(Iplus ).AND.enabled(Jplus )
    enabled(IminusJminus)= enabled(Iminus).AND.enabled(Jminus)
    enabled(IplusJminus) = enabled(Iplus ).AND.enabled(Jminus)
    enabled(IminusJplus )= enabled(Iminus).AND.enabled(Jplus )

    ! Main communications loop.
    maxrecvpts = MAXVAL(nrecvp(1:nrecv,1))
    maxsendpts = MAXVAL(nsendp(1:nsend,1))

    IF(PRESENT(b2) .OR. PRESENT(b3))THEN
       IF(.NOT. ALLOCATED(sendBuff))THEN
          ! Only allocate the sendBuff once
          ALLOCATE(recvBuff(maxrecvpts,nrecv), &
                   sendBuff(maxsendpts,nsend),stat=ierr)
       ELSE
          ALLOCATE(recvBuff(maxrecvpts,nrecv),stat=ierr)
       END IF
    ELSE IF(PRESENT(ib2) .OR. PRESENT(ib3))THEN
       IF(.NOT. ALLOCATED(sendIBuff))THEN
          ALLOCATE(recvIBuff(maxrecvpts,nrecv), &
                   sendIBuff(maxsendpts,nsend),stat=ierr)
       ELSE
          ALLOCATE(recvIBuff(maxrecvpts,nrecv),stat=ierr)
       END IF
    END IF

    IF (ierr .ne. 0) THEN
       CALL parallel_abort('exchs_generic: unable to allocate send/recvBuffs')
    END IF

    ! Initiate receives in case posting them first improves 
    ! performance.

    DO irecv=1,nrecv

       IF ( enabled(dirrecv(irecv)) .AND. &
            source(irecv).GE.0 .AND. nxrecv(irecv).GT.0 ) THEN

          tag = tag_orig + dirrecv(irecv)

          ! ARPDBG - nrecvp second rank is for multiple halo widths but
          !          that isn't used
          call post_receive(nrecvp2d(irecv,1), source(irecv), tag, &
               exch_flags(handle,irecv,indexr), rbuff=recvBuff(:,irecv))

          if(DEBUG_COMMS)then
             WRITE (*,FMT="(I4,': exchs post recv : hand = ',I2,' dirn = '," &
             //"I1,' src = ',I3,' tag = ',I4,' npoints = ',I6)") &
                  get_rank(),handle,dirrecv(irecv), &
                  source(irecv), tag, nrecvp(irecv,1)
          endif

       ELSE
          exch_flags(handle,irecv,indexr) = MPI_REQUEST_NULL
       END IF

    ENDDO

    if (.not. first_time) then        

       ! Check that all sends from previous call have completed before 
       ! we continue to modify the send buffers
       call msg_wait(nsend, exch_flags1d, irecv, all=.True.)
    else
        first_time = .FALSE.
    end if ! .not. first_time

    ! Send all messages in the communications list.

    DO isend=1,nsend

       IF ( enabled(dirsend(isend)) .AND. &
            destination(isend) >= 0 .AND. nxsend(isend) > 0 ) THEN

          isrc = isrcsend(isend)
          jsrc = jsrcsend(isend)
          nxs  =   nxsend(isend)
          nys  =   nysend(isend)

          tag = tag_orig + dirsend(isend)

          if( DEBUG_COMMS)then
             IF(PRESENT(b3))THEN
                WRITE (*,FMT="(I4,': handle ',I4,' tag ',I4,' sending to '," &
                     //"I4,' data ',I4,' direction ',I3)")                   &  
                     get_rank(), handle, tag, destination(isend),            &
                     nsendp(isend,1),dirsend(isend)
             ELSE IF(PRESENT(b2))THEN
                WRITE (*,FMT="(I4,': handle ',I4,' tag ',I4,' sending to '," &
                     //"I4,' data ',I4,' direction ',I3)")                   &  
                     get_rank(), handle, tag, destination(isend),            &
                     nsendp2d(isend,1),dirsend(isend)
             END IF
          endif

          ! Copy the data into the send buffer and send it...

          IF ( PRESENT(b2) )THEN

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
             
             DO j=jstart, jend, 1
                DO i=istart, iend, 1
                   ic = ic + 1
                   sendBuff(ic,isend) = b2(i,j)
                END DO
             END DO

             if(DEBUG_COMMS)then
                write(*,"(I3,': packed: ',6E15.4)") get_rank(), &
                     sendBuff(1:min(ic,6),isend)
             end if
             
             call post_send(sendBuff(1,isend),ic,destination(isend),tag, &
                            exch_flags(handle,isend,indexs))

          ELSEIF ( PRESENT(ib2) ) THEN

             ic = 0
!!$             pack_patches2i: DO ipatch=1, npatchsend(isend,1), 1
!!$                jstart = jsrcsendp(ipatch,isend,1)
!!$                istart = isrcsendp(ipatch,isend,1)
!!$                jend   = jstart+nysendp(ipatch,isend,1)-1
!!$                iend   = istart+nxsendp(ipatch,isend,1)-1
!!$
!!$                DO j=jstart, jend, 1
!!$                   DO i=istart, iend, 1
!!$                      ic = ic + 1
!!$                      sendIBuff(ic,isend) = ib2(i,j)
!!$                   END DO
!!$                END DO
!!$             END DO pack_patches2i
!!$
!!$             CALL MPI_Isend(sendIBuff(1,isend),ic, MPI_INTEGER, &
!!$                            destination(isend),tag,mpi_comm_world,&
!!$                            exch_flags(handle,isend,indexs),ierr)

          ELSEIF ( PRESENT(b3) )THEN

             ! CALL timing_start('3dr_pack')
!!$             ic = 0
!!$             pack_patches3r: DO ipatch=1,npatchsend(isend,1)
!!$
!!$                jstart = jsrcsendp(ipatch,isend,1)
!!$                istart = isrcsendp(ipatch,isend,1)
!!$                jend   = jstart+nysendp(ipatch,isend,1)-1
!!$                iend   = istart+nxsendp(ipatch,isend,1)-1
!!$if defined key_z_first
!!$                DO j=jstart, jend, 1
!!$                   DO i=istart, iend, 1
!!$                      DO k=1,mbkmax(i,j),1
!!$                      !DO k=1, nzsendp(ipatch,isend,1),1
!!$else
!!$                DO k=1,nzsendp(ipatch,isend,1),1
!!$                   DO j=jstart, jend, 1
!!$                      DO i=istart, iend, 1
!!$endif
!!$                         ic = ic + 1
!!$                         sendBuff(ic, isend) = b3(i,j,k)
!!$                      END DO
!!$                   END DO
!!$                END DO
!!$             END DO pack_patches3r

             ! CALL timing_stop('3dr_pack')

             !CALL MPI_Isend(sendBuff(1,isend),ic,                  &
             !               MPI_DOUBLE_PRECISION,                  &
             !               destination(isend), tag, mpi_comm_world, &
             !               exch_flags(handle,isend,indexs),ierr)

          ENDIF

       ELSE

          exch_flags(handle,isend,indexs) = MPI_REQUEST_NULL

       ENDIF ! direction is enabled and have something to send

    ENDDO ! Loop over sends

    if(DEBUG_COMMS)then
       write (*,FMT="(I3,': exch tag ',I4,' finished all ',I4,' sends')") &
            get_rank(), tag, nsend
    end if

    ! Wait on the receives that were posted earlier

    ! Copy just the set of flags we're interested in for passing 
    ! to MPI_waitany
    exch_flags1d(1:nrecv) = exch_flags(handle, 1:nrecv, indexr)

    if(DEBUG_COMMS)then
       WRITE(*,"(I3,': Doing waitany: nrecv =',I3,' handle = ',I3)") &
          get_rank(), nrecv,handle
    endif

    ! Get the first available message that we've received
    call msg_wait(nrecv, exch_flags1d, irecv)

    DO WHILE(irecv .ne. MPI_UNDEFINED)

       IF ( PRESENT(b2) ) THEN

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
             
          DO j=jstart, jend, 1
             DO i=istart, iend, 1
                ic = ic + 1
                b2(i,j) = recvBuff(ic,irecv)
             END DO
          END DO
          
          if(DEBUG_COMMS)then
             write(*,"(I3,': unpacked: ',6E15.4)") get_rank(), &
                  recvBuff(1:min(6,ic),irecv)
          end if
             
       ELSE IF ( PRESENT(ib2) ) THEN

          ! Copy received data back into array
          ic = 0
          jstart = jdesrecv(irecv)
          jend   = jstart+nyrecv(irecv)-1
          istart = idesrecv(irecv)
          iend   = istart+nxrecv(irecv)-1
          DO j=jstart, jend, 1
             DO i=istart, iend, 1
                ic = ic + 1
                ib2(i,j) = recvIBuff(ic,irecv)
             END DO
          END DO

       ELSE IF (PRESENT(b3) ) THEN

          ic = 0
          
          jstart = jdesrecv(irecv)
          jend   = jstart+nyrecv(irecv)-1
          istart = idesrecv(irecv)
          iend   = istart+nxrecv(irecv)-1

          DO k=1,nzrecv(irecv),1
             DO j=jstart, jend, 1
                DO i=istart, iend, 1
                   ic = ic + 1
                   b3(i,j,k) = recvBuff(ic,irecv)
                END DO
             END DO
          END DO
             
       END IF

       ! Wait for the next message
       call msg_wait(nrecv, exch_flags1d, irecv)
    END DO ! while irecv != MPI_UNDEFINED

    if(DEBUG_COMMS)then
       WRITE(*,"(I3,': Finished all ',I3,' receives for handle ',I3)") &
             get_rank(), nrecv, handle
    end if

    ! Copy just the set of flags we're interested in for passing to  
    ! MPI_waitall next time around  
    exch_flags1d(1:nsend) = exch_flags(handle, 1:nsend, indexs)

    ! Free the exchange communications handle.
    CALL free_exch_handle(handle)

    ! All receives done so we can safely free the MPI receive buffers
    IF( ALLOCATED(recvBuff) ) DEALLOCATE(recvBuff)
    IF( ALLOCATED(recvIBuff) )DEALLOCATE(recvIBuff)

  end subroutine exchs_generic

end module parallel_comms_mod