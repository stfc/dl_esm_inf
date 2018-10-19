module parallel_comms_mod
  use kind_params_mod, only: go_wp
  use parallel_utils_mod, only: get_num_ranks, get_rank, parallel_abort
  use decomposition_mod, only: subdomain_type, decomposition_type
  use mpi !> TODO move MPI-specific stuff out of here
  implicit none

  private

  integer, parameter :: jpk = 1 !< Only 1 vertical level
  logical, parameter :: DEBUG = .true.
  logical, parameter :: DEBUG_COMMS = .true.
  logical :: lwp !< Whether or not to write out from this rank

  integer, parameter :: halo_depthx = 2, halo_depthy = 2

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

  ! SMP 22 Sep 2009
  ! Alternative, run-length encoded communications lists
  ! omitting permanently dry points.
  ! Of these, isrcrecp, jsrcrecvp
  ! are set up but not currently used,
  ! and could be eliminated.

  ! Maximum number of patches a single halo communication can be broken 
  ! into when trimming dry points in msg_trim()
  INTEGER, PARAMETER :: MaxPatch=8
  INTEGER, SAVE, DIMENSION(MaxPatch,MaxComm,halo_depthx) :: &
                                    isrcsendp, jsrcsendp,&
                                    nxsendp, nysendp, nzsendp, &
                                    isrcrecvp, jsrcrecvp,&
                                    idesrecvp, jdesrecvp,&
                                    nxrecvp, nyrecvp, nzrecvp
  INTEGER, SAVE, DIMENSION(MaxComm,halo_depthx) :: npatchsend, npatchrecv
  ! Total number of points in each message
  INTEGER, SAVE, DIMENSION(MaxComm,halo_depthx) :: nsendp, nsendp2d, nrecvp, nrecvp2d

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

  ! cyclic_bc     True if a cyclic boundary condition is to be applied
  !               Set using the value of jperio.
  LOGICAL :: cyclic_bc

  ! Stores whether a domain's boundaries have been trimmed as
  ! trimmed(boundary, PE) where boundary is one of {n,e,s,w}idx
  ! for the Northern, Eastern, Southern or Western boundary, respectively.
  ! Allocated in finish_partition(), deallocated in...ARPDBG
  LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: trimmed
  INTEGER, PARAMETER :: nidx = 1, eidx= 2, sidx = 3, widx = 4

  ! Value representing land in the mask used for partitioning and message
  ! trimming
  INTEGER, PARAMETER :: LAND = 0

  ! Rather than trim to a point immediately next to a wet point, we 
  ! back off nextra points. If we don't do this then the sea-ice
  ! computation goes wrong because it does use values over the land
  ! that immediately border the ocean.
  INTEGER, SAVE :: nextra

  ! Public routines
  PUBLIC :: map_comms, iprocmap, exchmod_alloc, exchs_generic

  ! Public variables
  PUBLIC :: MaxComm,nsend,nrecv,nxsend,nysend,destination,dirrecv, &
            dirsend,isrcsend,jsrcsend,idesrecv, jdesrecv,          &
            nxrecv, nyrecv, source, cyclic_bc, idessend, jdessend

  PUBLIC :: nsendp,nsendp2d,nrecvp,nrecvp2d,npatchsend,npatchrecv, &
            nxsendp,nysendp,nzsendp,nxrecvp,nyrecvp,nzrecvp,       &
            idesrecvp,jdesrecvp,isrcsendp,jsrcsendp

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

  PUBLIC :: trimmed, nidx, eidx, sidx, widx, nextra

  ! Switch for trimming dry points from halo swaps
  LOGICAL, PARAMETER :: msgtrim   = .TRUE.

  ! Switch for trimming points below ocean floor from halo swaps
  ! Defaults to true unless set via NEMO_MSGTRIM_Z environment var.
  LOGICAL, PUBLIC, SAVE :: msgtrim_z

  ! exch_flags    -  Array of flag arrays for exchanges
  ! exch_flags1d  -  Array of only the current MPI receive operations
  ! exch_tag      -  The tag value associated with this exchange
  ! exch_busy     -  Indicates whether a slot in the flag array is being used

  INTEGER, PARAMETER :: indexs=1,indexr=2
  INTEGER, PARAMETER :: max_flags=40
  INTEGER, PARAMETER :: min_tag=0
  INTEGER :: current_tag,max_tag_used,max_tag,n_tag_cycles=0
  LOGICAL :: first_mod=.TRUE.

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

contains

  subroutine map_comms (decomp, tmask, pbc, ierr )
    !!------------------------------------------------------------------
    ! Maps out the communications requirements for the partitioned
    ! domain, adding communications descriptions to the list.
    !
    !     Mike Ashworth, CLRC Daresbury Laboratory, July 1999
    !!------------------------------------------------------------------
    IMPLICIT NONE

    ! Subroutine arguments.
    type(decomposition_type), target, intent(in) :: decomp
    INTEGER, INTENT(in), allocatable :: tmask(:,:)  ! decomposed mask: 0 for land, 1 for ocean
    logical, intent(in) :: pbc ! Whether mesh is periodic in x dimension
    INTEGER, INTENT(out):: ierr

    ! Local variables.
    INTEGER :: i, i1, i2, ihalo, iproc, iprocc, iprocx, &
               iprocy, j, j1, j2, nadd, naddmaxr, naddmaxs
    INTEGER :: ldiff0, ldiff1 ! Local vars for coping with wrapping of coords
    INTEGER :: imax, imin ! Max/min value of i that a halo strip can run 
                          ! to/from (used to avoid including E/W halo cells
                          ! in N/S transfers)
    INTEGER :: ielb_iproc, ieub_iproc ! Lower and upper bounds for proc
                                      ! iproc corrected for E & W halos
    INTEGER :: ielb_no_halo, ieub_no_halo ! Lower and upper bounds for local
                                          ! domain corrected for E & W halos
    INTEGER, DIMENSION(halo_depthx) :: idesr, jdesr, idess, jdess &
                                , isrcr, jsrcr, isrcs, jsrcs &
                                , nxr, nyr, nxs, nys
    logical :: addcorner
    integer :: nprocp, irank
    type(subdomain_type), pointer :: subdomain
    integer :: nimpp

    cyclic_bc = pbc ! \todo remove need to set this module variable
    
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

    ! =================================================================
    ! Looking at the border where we will
    ! send data in the minus I direction (Iplus) and
    ! receive data that has been sent in the plus I direction (Iminus).
    ! =================================================================

    ! Start from the lower bound of the sub-domain, and carry on looking
    ! for communications with neighbours until we have reached
    ! the upper bound.

    subdomain => decomp%subdomains(irank)
    jelb = subdomain%global%ystart
    jeub = subdomain%global%ystop
    ielb = subdomain%global%xstart
    ieub = subdomain%global%xstop

    ! \TODO support multiple subdomains per rank
    j1 = jelb
    DO WHILE (j1.LE.jeub)

       ! Look for the process which owns the neighbouring point in the
       ! minus I direction.

       iproc = iprocmap(decomp,ielb-1,j1)
       IF ( iproc.GT.0 ) THEN

          ! Find where in the j direction the common border between these
          ! sub-domains ends.

          j2 = MIN(jeub,decomp%subdomains(iproc)%internal%ystop)

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
            isrcr(ihalo) = decomp%subdomains(iproc)%internal%xstop - halo_depthx - ihalo + 1 ! nleit(iproc)-ihalo+1 
            idesr(ihalo) = ihalo ! Halo goes from 1..halo_depthx
            nxr(ihalo) = ihalo
            nxs(ihalo) = ihalo
          ENDDO

          ! MIN below allows for fact that NEMO sets nlci==nlei at the E 
          ! boundary of the global domain when using cyclic bc's
          idess(:) = MIN(decomp%subdomains(iproc)%internal%xstop+1,decomp%subdomains(iproc)%internal%xstop)  ! Send _to_ E halo column of iproc
          ! Source for a send is always within internal region
          jsrcs(:) = j1-jelb+subdomain%internal%ystart !Add nldj 'cos jsrcs is local address in domain
          jdess(:) = j1-decomp%subdomains(iproc)%global%ystart+decomp%subdomains(iproc)%internal%ystart ! ditto
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
            IF ( nadd.GT.0 ) THEN
               WRITE (*,"(I3,': Adding starting corner to send for halo ',I2,' with ',I3,' points')") irank, ihalo, nadd
            END IF
            end if
            nadd = MIN(ihalo,naddmaxr)
            jdesr(ihalo) = jdesr(ihalo) - nadd
            jsrcr(ihalo) = jsrcr(ihalo) - nadd
            nyr(ihalo) = nyr(ihalo)+nadd
            if(DEBUG)then
            IF ( nadd.GT.0 ) THEN
              WRITE (*,"(I3,': Adding starting corner to receive for halo ',I2,' with ',I3,' points')") irank,ihalo, nadd
            ENDIF
            end if
          ENDDO ! Loop over 'overlap' points in i direction

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
               WRITE (*,"(I3,': Adding starting corner to send for halo ',I2,' with ',I3,' points')") irank,ihalo, nadd
            ENDIF
end if
            nadd = MIN(ihalo,naddmaxr)
            nyr(ihalo) = nyr(ihalo)+nadd
if(DEBUG)then
            IF ( nadd.GT.0  ) THEN
              WRITE (*,"(I3,': Adding starting corner to receive for halo ',I2,' with ',I3,' points')") irank,ihalo, nadd
            ENDIF
end if
          ENDDO

!         Add a send and a receive to the lists for this section
!         of border.

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
!            WRITE (*,'(21x,6i6)') &
!                   isrcr(2),jsrcr(2),idesr(2),jdesr(2),nxr(2),nyr(2)
          ENDIF
end if

!         Move the start point to one beyond this strip.

          j1 = j2+1

        ELSE

!         No process found, continue searching up the boundary.

          j1 = j1+1
        ENDIF
      ENDDO

!     =================================================================
!     Looking at the border where we will
!     send data in the plus I direction (Iminus) and
!     receive data that has been sent in the minus I direction (Iplus).
!     =================================================================

!     Start from the lower bound of the sub-domain, and carry on looking
!     for communications with neighbours until we have reached
!     the upper bound.

      j1 = jelb
      DO WHILE (j1.LE.jeub)

!       Look for the process which owns the neighbouring point in the
!       plus I direction.

        iproc = iprocmap(decomp,ieub+1,j1)
        IF ( iproc.GT.0 ) THEN

!         Find where in the j direction the common border between these
!         sub-domains ends.

          j2 = MIN(jeub,decomp%subdomains(iproc)%internal%ystop)!pjeub(iproc))

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

!         Construct the rest of the data describing the zone,
!         convert to local indexes and extend to multiple halo widths.

          isrcr(:) = 1 + halo_depthx ! nldit(iproc) ARPDBG because NEMO sets nldi
                                ! to unity if nbondi == -1 (W boundary of 
                                ! global domain)
          DO ihalo=1,halo_depthx
             ! NEMO sets nlei = nlci if nbondi == 1 (E boundary of 
             ! global domain). Normally, nlci = jpi = halo_depthx + iesub + halo_depthx
             ! \TODO fix the next line as I don't think blah%nx includes halos
             ! at the minute
             isrcs(ihalo) = subdomain%global%nx - halo_depthx - ihalo + 1 ! Outermost halo -> innermost col.
             idess(ihalo) = ihalo ! Halo runs from 1..halo_depthx
             nxr(ihalo) = ihalo
             nxs(ihalo) = ihalo
          ENDDO
          idesr(:) = subdomain%internal%xstop !MIN(nlei + 1, subdomain%internal%xstop) ! Allow for case when on boundary
                                         ! of global domain and thus have 
                                         ! no (explicit) halo
          ! Source for a send is within local internal domain
          jsrcs(:) = j1-jelb+subdomain%internal%ystart
          ! Destination for a send is within halo on remote domain
          jdess(:) = j1-decomp%subdomains(iproc)%global%ystart+decomp%subdomains(iproc)%internal%ystart

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

!         Examine whether corner points should be added to the end.

          naddmaxr = 0
          naddmaxs = 0
          DO ihalo=1,halo_depthx,1

!           Send corner data while we have data to send
!           and while there is a point that depends on it.

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

!         Add a send and a receive to the lists for this section
!         of border.

          CALL addsend (nsend,Iminus,iproc-1,     &
                        isrcs,jsrcs,idess,jdess,nxs,nys,&
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN
if(DEBUG)then
          IF ( lwp ) THEN
            WRITE (*,'(a,7i6)') 'Adding send -1111   ',iproc-1, &
                  isrcs(1),jsrcs(1),idess(1),jdess(1),nxs(1),nys(1)
          ENDIF
end if

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
      ! Looking at the border where we will
      ! send data in the minus J direction (Jplus) and
      ! receive data that has been sent in the plus J direction (Jminus).
      ! =================================================================

      ! Ensure we don't include any values from the model domain W and E
      ! halos if we have cyclic b.c.'s
      ielb_no_halo = ielb
      ieub_no_halo = ieub
      IF(cyclic_bc)THEN
         ! West
         IF((.NOT. trimmed(widx,irank)) .AND. &
                        pilbext(irank)) ielb_no_halo = ielb_no_halo + halo_depthx
         ! East
         IF((.NOT. trimmed(eidx,irank)) .AND. &
                        piubext(irank)) ieub_no_halo = ieub_no_halo - halo_depthx
      END IF

      ! Start from the lower bound of the sub-domain (in global coords), 
      ! and carry on looking
      ! for communications with neighbours until we have reached
      ! the upper bound.

      imin = ielb_no_halo
      imax = ieub_no_halo

      i1 = imin
      DO WHILE (i1.LE.imax)

         ! Look for the process which owns the neighbouring point in the
         ! minus J direction.

         iproc = iprocmap(decomp,i1,jelb-1)
         IF ( iproc.GT.0 ) THEN

            ! Ensure we don't include halos from the global domain borders if
            ! we have cyclic bc's.
!          ielb_iproc = pielb(iproc)
            ieub_iproc = decomp%subdomains(iproc)%global%xstop
            IF(cyclic_bc)THEN
               !             IF(pilbext(iproc))ielb_iproc = pielb(iproc)+halo_depthx
               IF( (.NOT. trimmed(eidx,iproc)) .AND. &
                    piubext(iproc)) ieub_iproc = decomp%subdomains(iproc)%global%xstop-halo_depthx
            END IF

!           Find where in the i direction the common border between these
!           sub-domains ends.

            i2 = MIN(imax,ieub_iproc)

if(DEBUG)then
            WRITE (*,FMT="(I3,': ARPDBG strip for minus J is ',I3,',',I3)") &
                 irank-1,i1,i2
end if

!           |        |        |        |        ||
!           |        |        |        |        ||
!           |        |        |        |        ||
!           ------------------------------------------------
!                   ||        |        |        |
!                   ||        |        |        |
!                   ||        |        |        |


!           Construct the rest of the data describing the zone,
!           convert to local indexes and extend to multiple halo widths.
            ! Convert to local coords:
            ! Dist into zone = ipos - start + 1
            ! Pos in zone in local = (start of internal region) + (dist into zone) - 1
            ! Convert from global i1 to local i in current domain
            ! if i1 == nimpp then we must start at i=1 (because nimpp is absolute position
! of starting edge of domain including overlap regions)
            nimpp = subdomain%global%xstart
            isrcs(:) = i1 - nimpp + 1
            ! Convert from global i1 to local i in the destination domain
            idess(:) = i1- decomp%subdomains(iproc)%global%xstart + 1
            idesr(:) = isrcs(:)
            isrcr(:) = idess(:)
            nxr(:) = i2-i1+1
            nxs(:) = nxr(:)

            jsrcs(:) = subdomain%internal%ystart ! First row of 'internal' region of domain
            DO ihalo=1,halo_depthy,1
               ! Source for a receive must be within internal region
               jsrcr(ihalo) = decomp%subdomains(iproc)%internal%ystop-ihalo+1
               jdesr(ihalo) = ihalo ! Halo runs from 1..halo_depthy
               nyr(ihalo) = ihalo
               nys(ihalo) = ihalo
            ENDDO
            ! Destination for a send must be a halo region
            ! nlcjt(iproc) is always in a halo region. Not sure what 
            ! happens when halo wider than 1.
            jdess(:) = decomp%subdomains(iproc)%internal%ystop

!           Examine whether corner points should be added to the START.

            naddmaxr = 0
            naddmaxs = 0
            DO ihalo=1,halo_depthy,1

!             Send corner data while we have data to send
!             and while there is a point that depends on it.

            IF ( i1-ihalo.GE.imin .AND.  &
                 iprocmap(decomp,i1-ihalo,jelb-ihalo).GT.0 ) THEN
              naddmaxs = ihalo
            ENDIF

!           Receive corner data only when we are at the corner
!           and while the sending sub-domain is the same as for the edge.

            IF ( i1.EQ.imin .AND. &
                 iprocmap(decomp,i1-ihalo,jelb-ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO

!         Add the extra points.

          DO ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            idess(ihalo) = idess(ihalo)-nadd
            isrcs(ihalo) = isrcs(ihalo)-nadd
            nxs(ihalo) = nxs(ihalo)+nadd
if(DEBUG)then
            IF ( nadd.GT.0 ) THEN
               WRITE (*,"(I3,': Adding starting corner to send for halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
            ENDIF
end if
            nadd = MIN(ihalo,naddmaxr)
            idesr(ihalo) = idesr(ihalo)-nadd
            isrcr(ihalo) = isrcr(ihalo)-nadd
            nxr(ihalo) = nxr(ihalo)+nadd

if(DEBUG)then
            IF ( nadd.GT.0 ) THEN
              WRITE (*,"(I3,': Adding starting corner to receive for halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
            ENDIF
end if
          ENDDO

!         Examine whether corner points should be added to the END.

          naddmaxr = 0
          naddmaxs = 0
          DO ihalo=1,halo_depthy,1

!           Send corner data while we have data to send
!           and while there is a point that depends on it.

            IF ( i2+ihalo.LE.imax .AND. &
                 iprocmap(decomp,i2,jelb-ihalo).GT.0 ) THEN
              naddmaxs = ihalo
            ENDIF

!           Receive corner data only when we are at the corner
!           and while the sending sub-domain is the same as for the edge.

            IF ( i2.EQ.imax .AND. & 
                 iprocmap(decomp,i2+ihalo,jelb-ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO 

!         Add the extra points.

          DO ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            nxs(ihalo) = nxs(ihalo)+nadd
if(DEBUG)then
            IF ( nadd.GT.0 ) THEN
               WRITE (*,"(I3,': Adding starting corner to send for halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
            ENDIF
end if
            nadd = MIN(ihalo,naddmaxr)
            nxr(ihalo) = nxr(ihalo)+nadd
if(DEBUG)then
            IF ( nadd.GT.0 ) THEN
              WRITE (*,"(I3,': Adding starting corner to receive for halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
            ENDIF
end if
          ENDDO

!         Add a send and a receive to the lists for this section
!         of border.

          CALL addsend (nsend,Jplus,iproc-1,       &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          CALL addrecv (nrecv,Jminus,iproc-1,      &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

!         Move the start point to one beyond this strip.

          i1 = i2+1

        ELSE

!         No process found, continue searching up the boundary.

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

      imin = ielb_no_halo
      imax = ieub_no_halo
      i1 = imin

      DO WHILE (i1.LE.imax)

         ! Look for the process which owns the neighbouring point in the
         ! plus J direction.

         iproc = iprocmap(decomp,i1,jeub+1)
         IF ( iproc.GT.0 ) THEN
            ! Ensure we don't include halos from the global borders if we
            ! have cyclic b.c.'s
!           ielb_iproc = pielb(iproc)
            ieub_iproc = decomp%subdomains(iproc)%global%xstop! pieub(iproc)
            IF(cyclic_bc)THEN
!              IF(pilbext(iproc))ielb_iproc = pielb(iproc)+halo_depthx
               IF((.NOT. trimmed(eidx,iproc)) .AND. &
                            piubext(iproc))ieub_iproc = decomp%subdomains(iproc)%global%xstop-halo_depthx
            END IF

!           Find where in the i direction the common border between these
!           sub-domains ends.

            i2 = MIN(imax, ieub_iproc)

if(DEBUG)then
            WRITE (*,FMT="(I3,': ARPDBG strip for plus J is ',I3,',',I3)") &
                 irank-1,i1,i2
end if

            ! Construct the rest of the data describing the zone,
            ! convert to local indexes and extend to multiple halo widths.

            isrcs(:) = i1 - nimpp + 1
            idess(:) = i1 - decomp%subdomains(iproc)%global%xstart + 1
            idesr(:) = isrcs(:)
            isrcr(:) = idess(:)
            nxr(:) = i2-i1+1
            nxs(:) = nxr(:)

            ! Source for a receive must be within an internal region
            ! nldj incorporates whether or not lower halo exists
            jsrcr(:) = decomp%subdomains(iproc)%internal%ystart

            DO ihalo=1,halo_depthy,1
               jsrcs(ihalo) = subdomain%internal%ystop-ihalo+1 ! innermost row -> outermost halo
               jdess(ihalo) = ihalo        ! Halo runs from 1..halo_depthy
               nyr(ihalo)   = ihalo
               nys(ihalo)   = ihalo
            ENDDO
            jdesr(:) = subdomain%internal%ystop

!         Examine whether corner points should be added to the START.

          naddmaxr = 0
          naddmaxs = 0
          DO ihalo=1,halo_depthy,1

!           Send corner data while we have data to send
!           and while there is a point that depends on it.

            IF ( i1-ihalo.GE.imin .AND. iprocmap(decomp,i1,jeub+ihalo).GT.0 ) THEN
              naddmaxs = ihalo
            ENDIF

!           Receive corner data only when we are at the corner
!           and while the sending sub-domain is the same as for the edge.

            IF ( i1.EQ.imin .AND. &
                 iprocmap(decomp,i1-ihalo,jeub+ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO

!         Add the extra points.

          DO ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            idess(ihalo) = idess(ihalo) -nadd
            isrcs(ihalo) = isrcs(ihalo) -nadd
            nxs(ihalo) = nxs(ihalo)+nadd
if(DEBUG)then
            IF ( nadd.GT.0 ) THEN
               WRITE (*,"(I3,': Adding starting corner to send for halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
            ENDIF
end if
            nadd = MIN(ihalo,naddmaxr)
            idesr(ihalo) = idesr(ihalo) - nadd
            isrcr(ihalo) = isrcr(ihalo) - nadd
            nxr(ihalo) = nxr(ihalo)+nadd

if(DEBUG)then
            IF ( nadd.GT.0 ) THEN
               WRITE (*,"(I3,': Adding starting corner to receive for halo ',I2,' with ',I3,' points')") irank-1,ihalo, nadd
            ENDIF
end if
          ENDDO

!         Examine whether corner points should be added to the END.

          naddmaxr = 0
          naddmaxs = 0
          DO ihalo=1,halo_depthy,1

!           Send corner data while we have data to send
!           and while there is a point that depends on it.

            IF ( i2+ihalo.LE.imax .AND. &
                 iprocmap(decomp,i2,jeub+ihalo).GT.0 ) THEN
              naddmaxs = ihalo
            ENDIF

!           Receive corner data only when we are at the corner
!           and while the sending sub-domain is the same as for the edge.

            IF ( i2.EQ.imax .AND.       & ! Are we at the corner?
                 iprocmap(decomp,i2+ihalo,jeub+ihalo).EQ.iproc ) THEN
              naddmaxr = ihalo
            ENDIF
          ENDDO 

!         Add the extra points.

          DO ihalo=1,halo_depthy,1
            nadd = MIN(ihalo,naddmaxs)
            nxs(ihalo) = nxs(ihalo)+nadd
if(DEBUG)then
            IF ( nadd.GT.0 .AND. lwp ) THEN
              WRITE (*,*) irank-1,': Adding starting corner to send' &
                  ,' for halo ',ihalo,' with ',nadd,' points'
            ENDIF
end if
            nadd = MIN(ihalo,naddmaxr)
            nxr(ihalo) = nxr(ihalo)+nadd
if(DEBUG)then
            IF ( nadd.GT.0 .AND. lwp ) THEN
              WRITE (*,*) irank-1,': Adding starting corner to receive' &
                  ,' for halo ',ihalo,' with ',nadd,' points'
            ENDIF
end if
          ENDDO

!         Add a send and a receive to the lists for this section
!         of border.

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


      DO i=5,8

        ! At first assume there is to be no corner communication

        addcorner = .FALSE.

        ! i1 is to be x-coord just OUTSIDE our domain
        ! i2 is to be x-coord just WITHIN our domain
        ! All the following complexity is to allow for fact that the first 
        ! and last columns of the global domain are actually halos when
        ! we have cyclic E/W boundary conditions.
        IF( (iubext .OR. ilbext) .AND. cyclic_bc) THEN
           i1 = ielb
           i2 = ielb
           IF(ilbext .AND. (.NOT. trimmed(widx,irank)) )THEN
              i2 = i2+west(i) ! If on W boundary with cyclic bc's, ielb _is_ the halo column
                              ! so add 1 to move inside domain
           ELSE
              i1 = i1-west(i)
           END IF
           IF(iubext .AND. (.NOT. trimmed(eidx,irank)) )THEN
              ! If upper bound is on domain boundary then iesub already
              ! includes the halo column
              i1 = i1+east(i)*(iesub-1)
              i2 = i2+east(i)*(iesub-2)
           ELSE
              i1 = i1+east(i)*iesub
              i2 = i2+east(i)*(iesub-1)
           END IF

        ELSE
           i1 = ielb-west(i)+east(i)*iesub
           i2 = ielb+east(i)*(iesub-1)
        END IF

        ! For a NW corner:
        !               | 
        !       iproc   |  iprocy
        !       ________|______
        !               |
        !       iprocx  |  Me
        !               |

        ! x coord just OUTSIDE our domain but y INSIDE
        iprocx = iprocmap(decomp,i1, jelb+north(i)*(jesub-1))
        ! x coord just INSIDE our domain but y OUTSIDE
        iprocy = iprocmap(decomp,i2, jelb-south(i)+north(i)*jesub)

        iprocc = 0

if(DEBUG)then
        WRITE(*,FMT="(I3,': ARPDBG: i, i1 (outside), i2 (inside), iprocx, iprocy = ',5I4)") &
             irank-1, i,i1,i2,iprocx,iprocy
end if

        ! Loop over all required halo widths

        DO ihalo=halo_depthx,1,-1

          ! Look at the processor in the corner at width ihalo from the
          ! corner of the sub-domain. We want both x and y to be just
          ! outside our domain. i1 is already just outside our domain
          ! so we subtract one from ihalo below:
          iproc = iprocmap(decomp, i1 - west(i)*(ihalo-1)+ east(i)*(ihalo-1) &
                          ,south(i)*(jelb-ihalo)+north(i)*(jeub+ihalo))
!          iproc = iprocmap(decomp, west(i)*(ielb-ihalo)+ east(i)*(ieub+ihalo) &
!                          ,south(i)*(jelb-ihalo)+north(i)*(jeub+ihalo))

          ! If the corner processor is different from those to X and Y
          ! we will need a corner communication.

          IF ( iproc.GT.0 .AND. iprocx.GT.0 .AND. iprocy.GT.0 .AND. &
               iproc.NE.iprocx .AND. iproc.NE.iprocy ) THEN

             ! Ensure we don't include halos from the global borders if we
             ! have cyclic E/W boundaries.
             ielb_iproc = decomp%subdomains(iproc)%global%xstart
             ieub_iproc = decomp%subdomains(iproc)%global%xstop
             IF( cyclic_bc )THEN
                IF(pilbext(iproc))ielb_iproc = decomp%subdomains(iproc)%global%xstart+ihalo
                IF(piubext(iproc))ieub_iproc = decomp%subdomains(iproc)%global%xstop-ihalo
             END IF

if(DEBUG)then
             WRITE (*,FMT="(I3,': adding corner as ',I3,' differs from ',2I3)")&
                    irank-1, iproc,iprocx,iprocy
end if
            ! If the furthest corner point needs a communication,
            ! we will need them all.

            IF ( ihalo.EQ.halo_depthx ) THEN
              iprocc = iproc

              ! Ensure we don't include halos from the global borders if we
              ! have cyclic E/W boundaries.
              ielb_iproc = decomp%subdomains(iprocc)%global%xstart
              ieub_iproc = decomp%subdomains(iprocc)%global%xstop
              IF( cyclic_bc )THEN
                 IF(pilbext(iprocc))ielb_iproc = decomp%subdomains(iprocc)%global%xstart+halo_depthx
                 IF(piubext(iprocc))ieub_iproc = decomp%subdomains(iprocc)%global%xstop-halo_depthx
              END IF
              ! Set the flag to add everything to the communications list.

              addcorner = .TRUE.

            ENDIF

            ! Set the parameters for the communication.
            ldiff0 = ielb_iproc - ieub_no_halo
            ldiff1 = ielb_no_halo - ieub_iproc
            ! Allow for wrap-around if necessary
            IF(cyclic_bc)THEN
               IF(ldiff0 < 1)THEN
                  !ARPDBG -2 for consistency with procmap
                  ldiff0 = ldiff0 + (decomp%global_nx - 2)
               END IF
               IF(ldiff1 < 1)THEN
                  !ARPDBG -2 for consistency with procmap
                  ldiff1 = ldiff1 + (decomp%global_nx - 2)
               END IF
            END IF
            nxs  (ihalo) = ihalo -  east(i)*(ldiff0-1) &
                                 -  west(i)*(ldiff1-1)
            ! Have no cyclic b.c.'s in N/S direction so probably don't need
            ! the following checks on ldiff{0,1}
            ldiff0 = decomp%subdomains(iprocc)%internal%xstart - jeub
            IF(ldiff0 < 1) ldiff0 = ldiff0 + decomp%global_ny
            ldiff1 = jelb - decomp%subdomains(iprocc)%internal%ystop !pjeub(iprocc)
            IF(ldiff1 < 1) ldiff1 = ldiff1 + decomp%global_ny
            nys  (ihalo) = ihalo - north(i)*(ldiff0-1) &
                                 - south(i)*(ldiff1-1)

            ! Source for a send must be within the internal region of
            ! the LOCAL domain
            isrcs(ihalo) = east(i) *(iesub-nxs(ihalo)) + subdomain%internal%xstart
            jsrcs(ihalo) = north(i)*(jesub-nys(ihalo)) + subdomain%internal%ystart
            IF( cyclic_bc )THEN
               IF( ilbext )THEN
                  ! nldi is still within halo for domains on W edge of
                  ! global domain
                  isrcs(ihalo) = isrcs(ihalo) + west(i)
               ELSE IF( iubext )THEN
                  ! Final column is actually halo for domains on E edge of
                  ! global domain
                  isrcs(ihalo) = isrcs(ihalo) - east(i)
               END IF
            END IF

            ! Destination for a send must be within a halo region on the
            ! REMOTE domain
            ! MAX and MIN below to allow for edge domains that do not have
            ! explicit halo
            idess(ihalo) =  west(i)*(decomp%subdomains(iprocc)%internal%xstop & !(MIN(nleit(iprocc)+ihalo, decomp%subdomains(iprocc)%internal%xstop) &
                         - nxs(ihalo)+1) &
                         + east(i)*MAX(decomp%subdomains(iprocc)%internal%xstart-ihalo, 1)

            ! MAX and MIN below to allow for edge domains that do not have
            ! explicit halo
            jdess(ihalo) = south(i)*(decomp%subdomains(iprocc)%internal%ystop-nys(ihalo)+1) &
                         + north(i)*MAX(decomp%subdomains(iprocc)%internal%ystart - ihalo,1)

            ! Source for a receive must be in an internal region of the REMOTE domain
            isrcr(ihalo) = west(i)*(decomp%subdomains(iprocc)%internal%nx-nxs(ihalo)) + decomp%subdomains(iprocc)%internal%xstart
            IF( cyclic_bc )THEN

               ! This _could_ be a corner exchange wrapped around by the cyclic
               ! boundary conditions:
               !
               !  ||------|                      ||
               !  ||      |           |          ||
               !  ||a_____|__ _ _     |          ||
               !  ||                  -----------||
               !  ||                    |       a||
               !  ||                    |________||

               IF(pilbext(iprocc))THEN
                  ! nldi is still within halo for domains on E edge of
                  ! global domain
                  isrcr(ihalo) = isrcr(ihalo) + east(i)
               ELSE IF(piubext(iprocc))THEN
                  ! Final column is actually halo for domains on W edge of
                  ! global domain
                  isrcr(ihalo) = isrcr(ihalo) - west(i)
               END IF
            END IF
            jsrcr(ihalo) = south(i)*(decomp%subdomains(iprocc)%internal%ny-&
                        nys(ihalo)) + decomp%subdomains(iprocc)%internal%ystart

            ! Destination for a receive must be in a halo region (on LOCAL
            ! domain) and therefore:
            !  1 <= {i,j}desr <= jprec{i,j} or >= nle{i,j} + 1
            idesr(ihalo) =  east(i)*(subdomain%internal%xstop-nxs(ihalo)+1) & !ARPDBG s/iesub/nlei/
                         + west(i)*MAX(1,subdomain%internal%xstart-ihalo) !ARPDBG incl. nldi-

            jdesr(ihalo) = north(i)*(subdomain%internal%ystop-nys(ihalo) + 1) &
                         + south(i)*MAX(1,subdomain%internal%ystart - ihalo)

          ELSE

if(DEBUG)then
            IF ( iprocc.GT.0 ) THEN
                WRITE (*,FMT="(I3,': skipping corner for halo ',I3,' PE ',I3)")&
                       irank-1,ihalo,iprocc-1
            ENDIF
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

        ENDDO

        ! The size of the received data is always the same as 
        ! that of the sent data.

        nxr(:) = nxs(:)
        nyr(:) = nys(:)

        ! Add the data to the communications lists.

        IF ( addcorner ) THEN
if(DEBUG)then
          WRITE (*,FMT="(I3,': ARPDBG adding corner send to ',I4,', dir = ',I1)") &
                 irank-1, iprocc-1,i
end if
          CALL addsend (nsend,i,iprocc-1,          &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        tmask, ierr)
          IF ( ierr.NE.0 ) RETURN

          ! Manually reverse the direction indicator for the receive.
          j = opp_dirn(i)

if(DEBUG)then
          WRITE (*,FMT="(I3,': ARPDBG adding corner recv. from ',I4,', old dir = ',I1,' new dir = ',I1)") &
                 irank-1, iprocc-1,i, j
end if
          CALL addrecv (nrecv,j,iprocc-1,          &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        tmask, ierr)
          IF ( ierr.NE.0 ) RETURN

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
    ! Run-length encoded versions corresponding to above
    INTEGER, DIMENSION(MaxPatch,halo_depthx) :: risrc,rjsrc,rides,rjdes,rnx,rny,rnz
    ! Number of patches in run-length encoded message
    INTEGER, DIMENSION(halo_depthx) :: npatches
    INTEGER :: ihalo, ipatch, irank
    INTEGER :: nsendp_untrimmedz ! How many pts we'd be sending without
                                 ! trimming in z direction
    ! Whether there is still a message after clipping
    LOGICAL :: something_left
         
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
    nsendp(icomm, 1)   = nxsend(icomm)*nysend(icomm)
    nzsend(icomm)      = jpk
    
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
       WRITE(*,FMT="(I4,': ARPDBG:npatch = ',I3)") irank-1,npatches(1)
 
       WRITE (*,FMT="(I4,': ARPDBG:nsendp = ',I4)") irank-1,nsendp(icomm,1)
       WRITE (*,FMT="(I4,': ARPDBG SEND ends')")    irank-1
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

    !     Subroutine arguments.
    INTEGER,                 INTENT(inout) :: icomm
    INTEGER,                 INTENT( in  ) :: dir, proc
    INTEGER,                 INTENT(out)   :: ierr
    INTEGER, DIMENSION(:,:), INTENT( in  ) :: tmask
    INTEGER, DIMENSION(halo_depthx)        :: isrc, jsrc, ides, jdes, nx, ny

    ! Local variables.

    ! Run-length encoded versions corresponding to above
    INTEGER, dimension(MaxPatch,halo_depthx) :: risrc,rjsrc,rides,rjdes,rnx,rny,rnz
    ! Number of patches in run-length encoded message
    INTEGER, DIMENSION(halo_depthx) :: npatches
    INTEGER :: ihalo, ipatch, irank
    ! Whether there is still a message after clipping
    LOGICAL :: something_left

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

    dirrecv(icomm)  = dir
    source(icomm)   = proc
    isrcrecv(icomm) = isrc(1)
    jsrcrecv(icomm) = jsrc(1)
    idesrecv(icomm) = ides(1)
    jdesrecv(icomm) = jdes(1)
    nxrecv(icomm)   = nx(1)
    nyrecv(icomm)   = ny(1)
    !> TODO handle multiple halo depths
    nrecvp(icomm, 1) = nxrecv(icomm)*nyrecv(icomm)
    nzrecv(icomm)   = jpk
    
    if(DEBUG)then
       WRITE (*,FMT="(I3,': ARPDBG adding RECV:')") irank-1
       WRITE (*,FMT="(I3,': ARPDBG: icomm = ',I2)") irank-1,icomm
       WRITE (*,FMT="(I3,': ARPDBG:   dir = ',I2)") irank-1,dir
       WRITE (*,FMT="(I3,': ARPDBG:  proc = ',I4)") irank-1,proc
       WRITE (*,FMT="(I3,': ARPDBG:  isrc = ',I4)") irank-1,isrcrecv(icomm)
       WRITE (*,FMT="(I3,': ARPDBG:  jsrc = ',I4)") irank-1,jsrcrecv(icomm)
       WRITE (*,FMT="(I3,': ARPDBG:  ides = ',I4)") irank-1,idesrecv(icomm)
       WRITE (*,FMT="(I3,': ARPDBG:  jdes = ',I4)") irank-1,jdesrecv(icomm)
       WRITE (*,FMT="(I3,': ARPDBG:    nx = ',I4)") irank-1,nxrecv(icomm)
       WRITE (*,FMT="(I3,': ARPDBG:    ny = ',I4)") irank-1,nyrecv(icomm)
       WRITE (*,FMT="(I3,': ARPDBG:    nz = ',I4)") irank-1,nzrecv(icomm)
       WRITE (*,FMT="(I3,': ARPDBG RECV ends')")    irank-1
    end if

  END SUBROUTINE addrecv

    function iprocmap (decomp, ia, ja )
!!------------------------------------------------------------------
! Returns the process number (1...N) of the process whose sub-domain
! contains the point with global coordinates (i,j).
! If no process owns the point, returns zero.

!     i                       int   input     global x-coordinate
!     j                       int   input     global y-coordinate

!         Mike Ashworth, CLRC Daresbury Laboratory, July 1999
!         Andrew Porter, STFC Daresbury Laboratory, May  2008
      !!------------------------------------------------------------------
      implicit none

      !  Function arguments.
      integer                              :: iprocmap
      type(decomposition_type), intent(in) :: decomp
      integer,                  intent(in) :: ia, ja
      ! Local variables.
      integer :: iproc, i, j, iwidth, nprocp

      nprocp = get_num_ranks()
      iprocmap = 0

      ! Make sure we don't change variable values in calling
      ! routine...
      i = ia
      j = ja
      IF(cyclic_bc)THEN
         ! Allow for fact that first and last columns in global domain
         ! are actually halos
         iwidth = decomp%global_nx - 2*halo_depthx
         IF(i >= decomp%global_nx) i = ia - iwidth
         IF(i <= 1     ) i = ia + iwidth
         ! No cyclic BCs in North/South direction
         !IF(j > jpjglo) j = ja - jpjglo
         !IF(j < 1     ) j = ja + jpjglo
      END IF

        ! Search the processes for the one owning point (i,j).

      DO iproc=1,nprocp
         IF (decomp%subdomains(iproc)%global%xstart.LE.i .AND. &
             i.LE.decomp%subdomains(iproc)%global%xstop .AND.  &
             decomp%subdomains(iproc)%global%ystart.LE.j .AND. &
             j.LE.decomp%subdomains(iproc)%global%ystop) THEN
           iprocmap = iproc
           EXIT
         END IF
      ENDDO

! ARP - for debugging only
!!$        IF(iprocmap == 0)THEN
!!$           WRITE(*,"('iprocmap: failed to find owner PE for (',I3,1x,I3,')')") ia, ja
!!$           WRITE(*,*) 'PE domains are [xmin:xmax][ymin:ymax]:'
!!$           DO iproc=1,nprocp,1
!!$              WRITE(*,"(I3,': [',I3,':',I3,'][',I3,':',I3,']')") &
!!$                    iproc, pielb(iproc), pieub(iproc), pjelb(iproc), pjeub(iproc)
!!$           END DO
!!$        END IF

      END FUNCTION iprocmap

      
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

    INTEGER :: h,ierr
    LOGICAL :: got

    IF ( first_mod ) THEN

       ! First time in the module (i.e. exch or glob) set up the tags.

       ! Set the maximum tag value.

       got = .FALSE.
       CALL MPI_attr_get(MPI_comm_world,MPI_tag_ub,max_tag,got,ierr)
       IF ( ierr.NE.0 ) CALL abort ()

       IF ( .NOT.got ) THEN

          ! If no value was returned use the minimum possible tag max.
          ! (p. 28 of Version 2.1 of the MPI standard or p. 19 of V.1 of the standard.)
          max_tag = 32767
       ENDIF
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

       STOP 'ERROR: get_exch_handle: no free flag array'
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

       if (  DEBUG .or. DEBUG_COMMS)then
          IF ( lwp ) THEN
             WRITE (*,'(1x,a,i6,a,i8,a,i3,a,i3,a)')  &
               'Process ',get_rank(),' exch tag ',exch_tag(h) &
               ,' assigned flags ',h,' (',COUNT(exch_busy),' busy)'
          ENDIF
       endif
    ENDIF

    get_exch_handle = h

    RETURN

  END FUNCTION get_exch_handle

  ! ---------------------------------------------------------------

  SUBROUTINE free_exch_handle ( h )
    ! Frees exchange handle, h.
    IMPLICIT NONE

    ! Subroutine arguments.
    INTEGER :: h ! Handle to be free'd

    ! Free the flags array.
    
    IF ( h.GT.0 .AND. h.LE.max_flags ) THEN
       exch_busy(h) = .FALSE.
       if( DEBUG .or. DEBUG_COMMS)then
          IF ( lwp ) THEN
             WRITE (*,'(1x,a,i6,a,i8,a,i3)') 'Process ',get_rank(), &
                  ' exch tag ',exch_tag(h), ' freed    flags ',h
          endif
       endif
    ELSE
       WRITE (*,*) 'free_exch_handle: invalid handle ',h
    ENDIF

  END SUBROUTINE free_exch_handle

  !================================================================

  SUBROUTINE exchs_generic ( b2, ib2, b3, ib3, nhalo, nhexch, &
                             handle, comm1, comm2, comm3, comm4)

    ! *******************************************************************
    ! Send boundary data elements to adjacent sub-domains.

    ! b2(:,:)                real   input       2D real*8 local array.
    ! ib2(:,:)               int    input       2D integer local array.
    ! b3(:,:,:)              real   input       3D real*8 local array.
    ! ib3(:,:,:)             int    input       3D integer local array.
    ! nhalo                  int    input       Width of halo.
    ! nhexch                 int    input       Number of halo
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
    IMPLICIT none

    ! Subroutine arguments.
    INTEGER, INTENT(in)  :: nhalo,nhexch
    INTEGER, INTENT(out) :: handle

    REAL(go_wp),OPTIONAL, INTENT(inout), DIMENSION(:,:)   :: b2
    INTEGER, OPTIONAL, INTENT(inout), DIMENSION(:,:)   :: ib2
    REAL(go_wp),OPTIONAL, INTENT(inout), DIMENSION(:,:,:) :: b3
    INTEGER, OPTIONAL, INTENT(inout), DIMENSION(:,:,:) :: ib3

    INTEGER,           INTENT(in) :: comm1, comm2, comm3, comm4

    ! Local variables.

    LOGICAL :: enabled(0:MaxCommDir)
    INTEGER :: ierr, irecv, ircvdt, isend, isnddt, &
               isrc, jsrc, kdim1, &  ! ides, jdes, nxr, nyr,        &
               nxs, nys, tag, tag_orig
    INTEGER :: maxrecvpts, maxsendpts ! Max no. of grid points involved in 
                                      ! any one halo exchange
    INTEGER :: i, j, k, ic, ipatch ! Loop counters
    INTEGER :: istart, iend, jstart, jend
    INTEGER :: index  ! To hold index returned from MPI_waitany
    INTEGER, DIMENSION(3) :: isubsizes, istarts ! isizes
    INTEGER :: status(MPI_status_size)
    INTEGER :: astatus(MPI_status_size,MaxComm)
    LOGICAL, SAVE :: first_time = .TRUE.
    INTEGER, PARAMETER :: index_z = 3
    !!--------------------------------------------------------------------

    !CALL prof_region_begin(ARPEXCHS_GENERIC, "Exchs_indiv", iprofStat)
    !CALL timing_start('exchs_generic')

    ierr = 0

    ! Check nhexch is in range.

    if(nhexch > halo_depthx .or. nhexch > halo_depthy)then
       CALL parallel_abort('exchs: halo width greater than maximum')
    ENDIF

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

          if (DEBUG .or. DEBUG_COMMS)then
             WRITE (*,FMT="(I4,': tag ',I4,' ireceiving from ',I4,' data ',I4)") get_rank(),tag ,source(irecv), nrecvp(irecv,1)
          end if
          ! ARPDBG - nrecvp second rank is for multiple halo widths but
          !          that isn't used
          IF ( PRESENT(b2) ) THEN
             CALL MPI_irecv (recvBuff(1,irecv),nrecvp2d(irecv,1), &
                             MPI_DOUBLE_PRECISION, source(irecv), &
                             tag, mpi_comm_world,                   &
                             exch_flags(handle,irecv,indexr), ierr)
          ELSEIF ( PRESENT(ib2) ) THEN
             CALL MPI_irecv (recvIBuff(1,irecv),nrecvp2d(irecv,1), &
                             MPI_INTEGER, source(irecv),         &
                             tag, mpi_comm_world,                  &
                             exch_flags(handle,irecv,indexr),ierr)
          ELSEIF ( PRESENT(b3) ) THEN
             CALL MPI_irecv (recvBuff(1,irecv),nrecvp(irecv,1),   &
                             MPI_DOUBLE_PRECISION, source(irecv), &
                             tag, mpi_comm_world,                   &
                             exch_flags(handle,irecv,indexr),ierr)
          ELSEIF ( PRESENT(ib3) ) THEN
             CALL MPI_irecv (recvIBuff(1,irecv),nrecvp(irecv,1), &
                             MPI_INTEGER, source(irecv),         &
                             tag, mpi_comm_world,                  &
                             exch_flags(handle,irecv,indexr),ierr)
          ENDIF

          if(DEBUG_COMMS)then
             WRITE (*,FMT="(I4,': exchs post recv : hand = ',I2,' dirn = ',I1,' src = ',I3,' tag = ',I4,' npoints = ',I6)") &
                  get_rank(),handle,dirrecv(irecv), &
                  source(irecv), tag, nrecvp(irecv,1)
          endif

       ELSE
          exch_flags(handle,irecv,indexr) = MPI_REQUEST_NULL
       END IF

    ENDDO

    IF (.not. first_time) THEN        

       ! Check that all sends from previous call have completed before 
       ! we continue to modify the send buffers
       CALL MPI_waitall(nsend, exch_flags1d, astatus, ierr)
       IF ( ierr.NE.0 ) CALL MPI_abort(mpi_comm_world,1,ierr)

     ELSE
        first_time = .FALSE.
    END IF ! .not. first_time


    ! Send all messages in the communications list.

!    CALL timing_start('mpi_sends')

    DO isend=1,nsend

       IF ( enabled(dirsend(isend)) .AND. &
            destination(isend) >= 0 .AND. nxsend(isend) > 0 ) THEN

          isrc = isrcsend(isend)
          jsrc = jsrcsend(isend)
          nxs  =   nxsend(isend)
          nys  =   nysend(isend)

          tag = tag_orig + dirsend(isend)

          if( DEBUG .or. DEBUG_COMMS)then
             IF(PRESENT(b3))THEN
                WRITE (*,FMT="(I4,': handle ',I4,' tag ',I4,' sending to ',I4,' data ',I4,' direction ',I3)") &  
               get_rank(), handle, tag, destination(isend),nsendp(isend,1),dirsend(isend)
             ELSE IF(PRESENT(b2))THEN
                WRITE (*,FMT="(I4,': handle ',I4,' tag ',I4,' sending to ',I4,' data ',I4,' direction ',I3)") &  
               get_rank(), handle, tag, destination(isend),nsendp2d(isend,1),dirsend(isend)
             END IF
          endif

          ! Copy the data into the send buffer and send it...

          IF ( PRESENT(b2) )THEN

!             CALL timing_start('2dr_pack')
             ic = 0
             istart = isrcsend(isend)
             iend   = istart+nxsend(isend)-1
             jstart = jsrcsend(isend)
             jend   = jstart+nysend(isend)-1

             DO j=jstart, jend, 1
                DO i=istart, iend, 1
                   ic = ic + 1
                   sendBuff(ic,isend) = b2(i,j)
                END DO
             END DO

!             CALL timing_stop('2dr_pack')

             CALL MPI_Isend(sendBuff(1,isend),ic,MPI_DOUBLE_PRECISION, &
                            destination(isend),tag,mpi_comm_world, &
                            exch_flags(handle,isend,indexs),ierr)

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
!!$#if defined key_z_first
!!$                DO j=jstart, jend, 1
!!$                   DO i=istart, iend, 1
!!$                      DO k=1,mbkmax(i,j),1
!!$                      !DO k=1, nzsendp(ipatch,isend,1),1
!!$#else
!!$                DO k=1,nzsendp(ipatch,isend,1),1
!!$                   DO j=jstart, jend, 1
!!$                      DO i=istart, iend, 1
!!$#endif
!!$                         ic = ic + 1
!!$                         sendBuff(ic, isend) = b3(i,j,k)
!!$                      END DO
!!$                   END DO
!!$                END DO
!!$             END DO pack_patches3r

             ! CALL timing_stop('3dr_pack')

             CALL MPI_Isend(sendBuff(1,isend),ic,                  &
                            MPI_DOUBLE_PRECISION,                  &
                            destination(isend), tag, mpi_comm_world, &
                            exch_flags(handle,isend,indexs),ierr)

!!$#if defined DEBUG_COMMS
!!$             WRITE (*,FMT="(I4,': Isend of ',I3,' patches, ',I6,' points, to ',I3)") &
!!$                     narea-1, npatchsend(isend,1),ic, &
!!$                     destination(isend)
!!$#endif

           ELSEIF ( PRESENT(ib3) ) THEN

              ic = 0
!!$              jstart = jsrcsend(isend,1) !+nhalo
!!$              istart = isrcsend(isend,1) !+nhalo
!!$              jend   = jstart+nysend(isend,1)-1
!!$              iend   = istart+nxsend(isend,1)-1
!!$              DO k=1,nzsend(isend,1),1
!!$                 DO j=jstart, jend, 1
!!$                    DO i=istart, iend, 1
!!$                       ic = ic + 1
!!$                       sendIBuff(ic, isend) = ib3(i,j,k)
!!$                    END DO
!!$                 END DO
!!$              END DO
!!$
!!$             CALL MPI_Isend(sendIBuff(1,isend),ic,               &
!!$                            MPI_INTEGER,                         &
!!$                            destination(isend),tag,mpi_comm_world, &
!!$                            exch_flags(handle,isend,indexs),ierr)
          ENDIF

          !IF ( ierr.NE.0 ) CALL MPI_abort(mpi_comm_world,1,ierr)

       ELSE

          exch_flags(handle,isend,indexs) = MPI_REQUEST_NULL

       ENDIF ! direction is enabled and have something to send

    ENDDO ! Loop over sends

    ! CALL timing_stop('mpi_sends')

#if ( defined DEBUG && defined DEBUG_EXCHANGE ) || defined DEBUG_COMMS
    WRITE (*,FMT="(I3,': exch tag ',I4,' finished all sends')") get_rank(),tag
#endif

    ! Wait on the receives that were posted earlier

    ! CALL timing_start('mpi_recvs')

    ! Copy just the set of flags we're interested in for passing 
    ! to MPI_waitany
    exch_flags1d(1:nrecv) = exch_flags(handle, 1:nrecv, indexr)

    if( DEBUG_COMMS)then
       WRITE(*,"(I3,': Doing waitany: nrecv =',I3,' handle = ',I3)") &
          get_rank(), nrecv,handle
    endif

    CALL MPI_waitany (nrecv, exch_flags1d, irecv, status, ierr)
    IF ( ierr .NE. MPI_SUCCESS ) THEN

       IF(ierr .EQ. MPI_ERR_REQUEST)THEN
          WRITE (*,"(I3,': ERROR: exchs_generic: MPI_waitany returned MPI_ERR_REQUEST')") get_rank()
       ELSE IF(ierr .EQ. MPI_ERR_ARG)THEN
          WRITE (*,"(I3,': ERROR: exchs_generic: MPI_waitany returned MPI_ERR_ARG')") get_rank()
       ELSE
          WRITE (*,"(I3,': ERROR: exchs_generic: MPI_waitany returned unrecognised error')") get_rank()
       END IF
       CALL parallel_abort('exchs_generic: MPI_waitany returned error')
    END IF

    DO WHILE(irecv .ne. MPI_UNDEFINED)

          IF ( PRESENT(b2) ) THEN

             ! CALL timing_start('2dr_unpack')

             ! Copy received data back into array
             ic = 0
             unpack_patches2r: DO ipatch=1,npatchrecv(irecv,nhexch)

                jstart = jdesrecvp(ipatch,irecv,1)!+nhalo
                jend   = jstart+nyrecvp(ipatch,irecv,1)-1
                istart = idesrecvp(ipatch,irecv,1)!+nhalo
                iend   = istart+nxrecvp(ipatch,irecv,1)-1
                DO j=jstart, jend, 1
                   DO i=istart, iend, 1
                      ic = ic + 1
                      b2(i,j) = recvBuff(ic,irecv)
                   END DO
                END DO
             END DO unpack_patches2r

             ! CALL timing_stop('2dr_unpack')

          ELSE IF ( PRESENT(ib2) ) THEN

             ! Copy received data back into array
             ic = 0
             unpack_patches2i: DO ipatch=1,npatchrecv(irecv,nhexch),1

                jstart = jdesrecvp(ipatch,irecv,1)
                jend   = jstart+nyrecvp(ipatch,irecv,1)-1
                istart = idesrecvp(ipatch,irecv,1)
                iend   = istart+nxrecvp(ipatch,irecv,1)-1
                DO j=jstart, jend, 1
                   DO i=istart, iend, 1
                      ic = ic + 1
                      ib2(i,j) = recvIBuff(ic,irecv)
                   END DO
                END DO
             END DO unpack_patches2i

           ELSE IF (PRESENT(b3) ) THEN

              ! CALL timing_start('3dr_unpack')
             ic = 0
             unpack_patches3r: DO ipatch=1,npatchrecv(irecv,nhexch)

                jstart = jdesrecvp(ipatch,irecv,1)!+nhalo
                jend   = jstart+nyrecvp(ipatch,irecv,1)-1
                istart = idesrecvp(ipatch,irecv,1)!+nhalo
                iend   = istart+nxrecvp(ipatch,irecv,1)-1

                DO k=1,nzrecvp(ipatch,irecv,1),1
                   DO j=jstart, jend, 1
                      DO i=istart, iend, 1
                         ic = ic + 1
                         b3(i,j,k) = recvBuff(ic,irecv)
                      END DO
                   END DO
                END DO

                ! ARPDBG - wipe anything below the ocean bottom
                DO k=nzrecvp(ipatch,irecv,1)+1,jpk,1
                   DO j=jstart, jend, 1
                      DO i=istart, iend, 1
                         b3(i,j,k) = 0.0_go_wp
                      END DO
                   END DO
                END DO
             END DO unpack_patches3r

!             CALL timing_stop('3dr_unpack')

          ELSEIF ( PRESENT(ib3) ) THEN

             ic = 0
             unpack_patches3i: DO ipatch=1,npatchrecv(irecv,nhexch),1

                jstart = jdesrecvp(ipatch,irecv,1)!+nhalo
                jend   = jstart+nyrecvp(ipatch,irecv,1)-1
                istart = idesrecvp(ipatch,irecv,1)!+nhalo
                iend   = istart+nxrecvp(ipatch,irecv,1)-1
                DO k=1,nzrecvp(ipatch,irecv,1),1
                   DO j=jstart, jend, 1
                      DO i=istart, iend, 1
                         ic = ic + 1
                         ib3(i,j,k) = recvIBuff(ic,irecv)
                      END DO
                   END DO
                END DO
             END DO unpack_patches3i

          END IF

       CALL MPI_waitany (nrecv, exch_flags1d, irecv, status, ierr)
       !IF ( ierr.NE.0 ) CALL MPI_abort(mpi_comm_world,1,ierr)

    END DO ! while irecv != MPI_UNDEFINED

    ! CALL timing_stop('mpi_recvs')

    ! All receives done and unpacked so can deallocate the associated
    ! buffers
    !IF(ALLOCATED(recvBuff ))DEALLOCATE(recvBuff)
    !IF(ALLOCATED(recvIBuff))DEALLOCATE(recvIBuff)

#if defined DEBUG_COMMS
    WRITE(*,"(I3,': Finished all ',I3,' receives for handle ',I3)") &
             get_rank(), nrecv, handle
#endif

    ! Copy just the set of flags we're interested in for passing to  
    ! MPI_waitall next time around  
    exch_flags1d(1:nsend) = exch_flags(handle, 1:nsend, indexs)

    ! Free the exchange communications handle.
    CALL free_exch_handle(handle)

    ! All receives done so we can safely free the MPI receive buffers
    IF( ALLOCATED(recvBuff) ) DEALLOCATE(recvBuff)
    IF( ALLOCATED(recvIBuff) )DEALLOCATE(recvIBuff)

    !CALL timing_stop('exchs_generic','section')
    !CALL prof_region_end(ARPEXCHS_GENERIC, iprofStat)

  END SUBROUTINE exchs_generic

end module parallel_comms_mod
