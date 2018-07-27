module parallel_comms_mod
  use parallel_utils_mod, only: get_num_ranks, get_rank
  use decomposition_mod, only: subdomain_type, decomposition_type
  implicit none

  private

  integer, parameter :: jpk = 1 !< Only 1 vertical level
  logical, parameter :: DEBUG = .true.
  logical :: lwp !< Whether or not to write out from this rank

  integer, parameter :: halo_depthx = 2, halo_depthy = 2

  ! Process ids of ensemble member processes in a linear list.
  INTEGER, ALLOCATABLE, DIMENSION(:) :: procid

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

  INTEGER, ALLOCATABLE, DIMENSION(:) :: pnactive,              &
       piesub, pjelb, pjeub, pjesub
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
  PUBLIC :: map_comms, iprocmap

  ! Public variables
  PUBLIC :: MaxComm,nsend,nrecv,nxsend,nysend,destination,dirrecv, &
            dirsend,isrcsend,jsrcsend,idesrecv, jdesrecv,          &
            nxrecv, nyrecv, source, cyclic_bc, idessend, jdessend

  PUBLIC :: nsendp,nsendp2d,nrecvp,nrecvp2d,npatchsend,npatchrecv, &
            nxsendp,nysendp,nzsendp,nxrecvp,nyrecvp,nzrecvp,       &
            idesrecvp,jdesrecvp,isrcsendp,jsrcsendp

  PUBLIC :: ielb,  ieub,  pjelb, pjeub,                    &
            iesub, jesub, jeub, ilbext, iubext, jubext, jlbext, pnactive,&
            piesub, pjesub, jelb, pilbext, pjlbext, pjubext, piubext

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

!    ALLOCATE(procid(nprocp), Stat=ierr)
!    IF (ierr > 0) THEN
!       WRITE (*,*) 'ERROR: mapcomms: Allocate failed for iproc'
!       RETURN
!    END IF
!    ! Create ordered list of process ids
!    DO i=1,nprocp,1
!       procid(i) = i-1
!    END DO

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

          j2 = MIN(jeub,pjeub(iproc))

! nldit(iproc) = decomp%subdomains(iproc)%internal%xstart
! nleit(iproc) = decomp%subdomains(iproc)%internal%xstop
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

          ! NEMO sets nldi=1 for nbondi==-1 but then calculates the source
          ! of data to send as follows...
          ! Since this is not == nldi, we don't use nldi here as one
          ! might have expected
          isrcs(:) = halo_depthx + 1

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

          CALL addsend (nsend,Iplus,procid(iproc), &
                        isrcs,jsrcs,idess,jdess,   &
                        nxs,nys,tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          ! This is a receive for data _from_ the neighbouring domain
          ! - it is NOT the corresponding receive for the above send.
          CALL addrecv (nrecv,Iminus,procid(iproc), &
                        isrcr,jsrcr,idesr,jdesr,    &
                        nxr,nyr,tmask,ierr)
          IF ( ierr.NE.0 ) RETURN
if(DEBUG)then
          IF ( lwp ) THEN
            WRITE (*,'(a,7i6)') 'Adding receive ',procid(iproc) &
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

          j2 = MIN(jeub,pjeub(iproc))

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

          CALL addsend (nsend,Iminus,procid(iproc),     &
                        isrcs,jsrcs,idess,jdess,nxs,nys,&
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN
if(DEBUG)then
          IF ( lwp ) THEN
            WRITE (*,'(a,7i6)') 'Adding send -1111   ',procid(iproc) &
                  ,isrcs(1),jsrcs(1),idess(1),jdess(1),nxs(1),nys(1)
          ENDIF
end if

          CALL addrecv (nrecv,Iplus,procid(iproc),       &
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

          CALL addsend (nsend,Jplus,procid(iproc),       &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          CALL addrecv (nrecv,Jminus,procid(iproc),      &
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

          CALL addsend (nsend,Jminus,procid(iproc),      &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        tmask,ierr)
          IF ( ierr.NE.0 ) RETURN

          CALL addrecv (nrecv,Jplus,procid(iproc),       &
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
            ldiff0 = pjelb(iprocc) - jeub
            IF(ldiff0 < 1) ldiff0 = ldiff0 + decomp%global_ny
            ldiff1 = jelb - pjeub(iprocc)
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
            isrcr(ihalo) = west(i)*(piesub(iprocc)-nxs(ihalo)) + decomp%subdomains(iprocc)%internal%xstart
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
            jsrcr(ihalo) = south(i)*(pjesub(iprocc)-nys(ihalo)) + decomp%subdomains(iprocc)%internal%ystart

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
                 irank-1, procid(iprocc),i
end if
          CALL addsend (nsend,i,procid(iprocc),          &
                        isrcs,jsrcs,idess,jdess,nxs,nys, &
                        tmask, ierr)
          IF ( ierr.NE.0 ) RETURN

          ! Manually reverse the direction indicator for the receive.
          j = opp_dirn(i)

if(DEBUG)then
          WRITE (*,FMT="(I3,': ARPDBG adding corner recv. from ',I4,', old dir = ',I1,' new dir = ',I1)") &
                 irank-1, procid(iprocc),i, j
end if
          CALL addrecv (nrecv,j,procid(iprocc),          &
                        isrcr,jsrcr,idesr,jdesr,nxr,nyr, &
                        tmask, ierr)
          IF ( ierr.NE.0 ) RETURN

        ENDIF

      ENDDO

      DEALLOCATE(procid)
      
    END SUBROUTINE map_comms

  subroutine addsend (icomm, dir, proc, isrc, jsrc, &
                      ides, jdes, nx, ny, tmask, ierr )
!!------------------------------------------------------------------
!     Adds a send communication specified by the parameters dir through 
!     to ny to the send communication list at the next position. 
!     icomm points to the last entry and is incremented and returned 
!     if successful.
!
!     icomm                   int   in/out    Location in comms list.
!     dir                     int   input     Direction.
!     proc                    int   input     Process id.
!     isrc                    int   input     X coordinate of source data.
!     jsrc                    int   input     Y coordinate of source data.
!     ides                    int   input     X coordinate of destination data.
!     jdes                    int   input     Y coordinate of destination data.
!     nx                      int   input     Size in X of data to be sent.
!     ny                      int   input     Size in Y of data to be sent.
!     ierr                    int   output    Error flag.
!
!               Mike Ashworth, CLRC Daresbury Laboratory, March 1999
!               Stephen Pickles, STFC Daresbury Laboratory
!                  - Sep 2009: added run-length encoding
!!------------------------------------------------------------------
         IMPLICIT NONE

         ! Subroutine arguments.
         INTEGER,                    INTENT(inout) :: icomm
         INTEGER,                    INTENT( in  ) :: dir, proc
         INTEGER, DIMENSION(:,:),    INTENT( in  ) :: tmask
         INTEGER,                    INTENT( out ) :: ierr
         INTEGER, DIMENSION(halo_depthx), INTENT( in  ) :: isrc, jsrc, &
                                                      ides, jdes, nx, ny
         ! Values of corresponding input arguments after clipping
         INTEGER, DIMENSION(halo_depthx) :: cisrc,cjsrc,cides,cjdes,cnx,cny,cnz
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
            isrcsend(icomm)    = cisrc(1)
            jsrcsend(icomm)    = cjsrc(1)
            idessend(icomm)    = cides(1)
            jdessend(icomm)    = cjdes(1)
            nxsend(icomm)      = cnx(1)
            nysend(icomm)      = cny(1)
            IF(msgtrim_z)THEN
               nzsend(icomm)   = cnz(1)
            ELSE
               nzsend(icomm)   = jpk
            END IF

            ! Zero count of untrimmed pts to send
            nsendp_untrimmedz = 0

            ! Also set up the comms lists encoded as the start points and
            ! lengths of the contiguous runs of wet points.
            DO ihalo=1,halo_depthx

               nsendp2d(icomm,ihalo)   = 0
               nsendp(icomm,ihalo)     = 0
               npatchsend(icomm,ihalo) = npatches(ihalo)

               DO ipatch=1,npatches(ihalo)

                  isrcsendp(ipatch,icomm,ihalo) = risrc(ipatch,ihalo)
                  jsrcsendp(ipatch,icomm,ihalo) = rjsrc(ipatch,ihalo)

                  nxsendp(ipatch,icomm,ihalo)   = rnx(ipatch,ihalo)
                  nysendp(ipatch,icomm,ihalo)   = rny(ipatch,ihalo)
                  IF(msgtrim_z)THEN
                     nzsendp(ipatch,icomm,ihalo)= rnz(ipatch,ihalo)
                  ELSE
                     nzsendp(ipatch,icomm,ihalo) = jpk
                  END IF

                  ! Sum the no. of points to be sent over all
                  ! patches for both 2D-array halos and 3D-array halos
                  nsendp2d(icomm,ihalo) = nsendp2d(icomm,ihalo) +      &
                                          nxsendp(ipatch,icomm,ihalo)* &
                                          nysendp(ipatch,icomm,ihalo)
                  nsendp(icomm,ihalo) = nsendp(icomm,ihalo) +          &
                                          nxsendp(ipatch,icomm,ihalo)* &
                                          nysendp(ipatch,icomm,ihalo)* &
                                          nzsendp(ipatch,icomm,ihalo)
                  IF(msgtrim_z)THEN
                     nsendp_untrimmedz = nsendp_untrimmedz +           &
                                          nxsendp(ipatch,icomm,ihalo)* &
                                          nysendp(ipatch,icomm,ihalo)* &
                                          jpk
                  END IF
               END DO
            END DO

if(DEBUG)then
            WRITE (*,FMT="(I4,': ARPDBG adding SEND:')") irank-1  
            WRITE (*,FMT="(I4,': ARPDBG: icomm = ',I2)") irank-1,icomm
            WRITE (*,FMT="(I4,': ARPDBG:   dir = ',I2)") irank-1,dirsend(icomm)
            WRITE (*,FMT="(I4,': ARPDBG:  proc = ',I4)") irank-1,destination(icomm)
            WRITE (*,FMT="(I4,': ARPDBG:  isrc = ',I4)") irank-1,isrcsend(icomm)
            WRITE (*,FMT="(I4,': ARPDBG:  jsrc = ',I4)") irank-1,jsrcsend(icomm)
            WRITE (*,FMT="(I4,': ARPDBG:  ides = ',I4)") irank-1,idessend(icomm)
            WRITE (*,FMT="(I4,': ARPDBG:  jdes = ',I4)") irank-1,jdessend(icomm)
            WRITE (*,FMT="(I4,': ARPDBG:    nx = ',I4)") irank-1,nxsend(icomm)
            WRITE (*,FMT="(I4,': ARPDBG:    ny = ',I4)") irank-1,nysend(icomm)
            WRITE (*,FMT="(I4,': ARPDBG:    nz = ',I4)") irank-1,nzsend(icomm)
            WRITE (*,FMT="(I4,': ARPDBG:npatch = ',I3)") irank-1,npatches(1)
 
            DO ipatch=1,npatches(1)
               WRITE (*,FMT="(I4,': ARPDBG:  patch ',I2,': isrc = ',I4)") &
                                  irank-1,ipatch,isrcsendp(ipatch,icomm,1)
               WRITE (*,FMT="(I4,': ARPDBG:  patch ',I2,': jsrc = ',I4)") &
                                  irank-1,ipatch,jsrcsendp(ipatch,icomm,1)
               WRITE (*,FMT="(I4,': ARPDBG:  patch ',I2,':   nx = ',I4)") &
                                  irank-1,ipatch,nxsendp(ipatch,icomm,1)  
               WRITE (*,FMT="(I4,': ARPDBG:  patch ',I2,':   ny = ',I4)") &
                                  irank-1,ipatch,nysendp(ipatch,icomm,1)  
               WRITE (*,FMT="(I4,': ARPDBG:  patch ',I2,':   nz = ',I4)") &
                                  irank-1,ipatch,nzsendp(ipatch,icomm,1)  
            END DO

            WRITE (*,FMT="(I4,': ARPDBG:nsendp = ',I4)") irank-1,nsendp(icomm,1)
            IF(msgtrim_z)THEN
               WRITE (*,FMT="(I4,': ARPDBG:nsendp WITHOUT z trim = ',I4)") &
                                                       irank-1,nsendp_untrimmedz
            END IF
            WRITE (*,FMT="(I4,': ARPDBG SEND ends')")    irank-1
end if

    END SUBROUTINE addsend

    SUBROUTINE addrecv ( icomm, dir, proc, isrc, jsrc, &
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
      IMPLICIT NONE

!     Subroutine arguments.
      INTEGER,                 INTENT(inout) :: icomm
      INTEGER,                 INTENT( in  ) :: dir, proc
      INTEGER,                 INTENT(out)   :: ierr
      INTEGER, DIMENSION(:,:), INTENT( in  ) :: tmask
      INTEGER, DIMENSION(halo_depthx)             :: isrc, jsrc, ides, jdes, nx, ny

      ! Local variables.

      ! Values of corresponding input arguments after clipping
      INTEGER, DIMENSION(halo_depthx) :: cisrc,cjsrc,cides,cjdes,cnx,cny,cnz
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
         isrcrecv(icomm) = cisrc(1)
         jsrcrecv(icomm) = cjsrc(1)
         idesrecv(icomm) = cides(1)
         jdesrecv(icomm) = cjdes(1)

         nxrecv(icomm)   = cnx(1)
         nyrecv(icomm)   = cny(1)
         IF(msgtrim_z)THEN
            nzrecv(icomm)   = cnz(1)
         ELSE
            nzrecv(icomm)   = jpk
         END IF

         DO ihalo=1,halo_depthx

            nrecvp2d(icomm,ihalo) = 0
            nrecvp(icomm,ihalo)   = 0
            npatchrecv(icomm,ihalo) = npatches(ihalo)

            DO ipatch=1,npatches(ihalo)
               isrcrecvp(ipatch,icomm,ihalo) = risrc(ipatch,ihalo)
               jsrcrecvp(ipatch,icomm,ihalo) = rjsrc(ipatch,ihalo)
               idesrecvp(ipatch,icomm,ihalo) = rides(ipatch,ihalo)
               jdesrecvp(ipatch,icomm,ihalo) = rjdes(ipatch,ihalo)
               nxrecvp(ipatch,icomm,ihalo)   = rnx(ipatch,ihalo)
               nyrecvp(ipatch,icomm,ihalo)   = rny(ipatch,ihalo)
               IF(msgtrim_z)THEN
                  nzrecvp(ipatch,icomm,ihalo) = rnz(ipatch,ihalo)
               ELSE
                  nzrecvp(ipatch,icomm,ihalo) = jpk
               END IF

               ! Sum the no. of points to be received over all
               ! patches
               nrecvp2d(icomm,ihalo) = nrecvp2d(icomm,ihalo) +           &
                                            nxrecvp(ipatch,icomm,ihalo)* &
                                            nyrecvp(ipatch,icomm,ihalo)
                    
               nrecvp(icomm,ihalo) = nrecvp(icomm,ihalo) +               &
                                            nxrecvp(ipatch,icomm,ihalo)* &
                                            nyrecvp(ipatch,icomm,ihalo)* &
                                            nzrecvp(ipatch,icomm,ihalo)
            END DO
         END DO

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
         WRITE (*,FMT="(I3,': ARPDBG:npatch = ',I3)") irank-1,npatches(1)
         DO ipatch=1,npatches(1)
            WRITE (*,FMT="(I3,': ARPDBG:  patch ',I2,': isrc = ',I4)") &
                                  irank-1,ipatch,isrcrecvp(ipatch,icomm,1)
            WRITE (*,FMT="(I3,': ARPDBG:  patch ',I2,': jsrc = ',I4)") &
                                  irank-1,ipatch,jsrcrecvp(ipatch,icomm,1)
            WRITE (*,FMT="(I3,': ARPDBG:  patch ',I2,': ides = ',I4)") &
                                  irank-1,ipatch,idesrecvp(ipatch,icomm,1)
            WRITE (*,FMT="(I3,': ARPDBG:  patch ',I2,': jdes = ',I4)") &
                                  irank-1,ipatch,jdesrecvp(ipatch,icomm,1)
            WRITE (*,FMT="(I3,': ARPDBG:  patch ',I2,':   nx = ',I4)") &
                                  irank-1,ipatch,nxrecvp(ipatch,icomm,1)  
            WRITE (*,FMT="(I3,': ARPDBG:  patch ',I2,':   ny = ',I4)") &
                                  irank-1,ipatch,nyrecvp(ipatch,icomm,1)  
            WRITE (*,FMT="(I3,': ARPDBG:  patch ',I2,':   nz = ',I4)") &
                                  irank-1,ipatch,nzrecvp(ipatch,icomm,1)  
         END DO
         WRITE (*,FMT="(I3,': ARPDBG:nrecvp = ',I4)") irank-1,nrecvp(icomm,1)
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
           IF ( decomp%subdomains(iproc)%global%xstart.LE.i .AND. i.LE.decomp%subdomains(iproc)%global%xstop .AND. &
                pjelb(iproc).LE.j .AND. j.LE.pjeub(iproc) ) THEN
              iprocmap = iproc
              EXIT
           ENDIF
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

end module parallel_comms_mod
