!------------------------------------------------------------------------------
! BSD 2-Clause License
! 
! Copyright (c) 2018, Science and Technology Facilities Council.
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

!> Implementation of parallel (distributed-memory) support using MPI.
!! Requires that the compiler be able to find the 'mpi' module.
module parallel_utils_mod
  use mpi
  use kind_params_mod, only: go_wp
  implicit none

  private

  !> Error flag for MPI calls. Note that this module makes no attempt to
  !! check these since we assume that MPI will abort if an error occurs
  !! (which is the default behaviour on the majority of systems).
  integer :: mpierr

  !> Our MPI communicator
  integer :: comm
  !> MPI rank + 1 of current process
  integer :: rank
  !> Total no. of MPI processes
  integer :: nranks

  !> TODO get rid of these
  logical, parameter :: DEBUG=.True., lwp=.True.

  !> If this module has been compiled in to the library then we have
  !! been built with distributed-memory support (MPI).
  logical, parameter :: DIST_MEM_ENABLED = .True.

  public parallel_init, parallel_finalise, parallel_abort, get_max_tag
  public get_rank, get_num_ranks, post_receive, post_send, msg_wait
  ! Re-export some MPI constants
  public MPI_UNDEFINED, MPI_REQUEST_NULL, DIST_MEM_ENABLED

contains

  !================================================

  !> Initialise the parallel program. Calls mpi_init() and then
  !! queries and stores the size of the communicator and the
  !! rank of the current process.
  subroutine parallel_init()

    call mpi_init(mpierr)
    call mpi_comm_dup(MPI_COMM_WORLD, comm, mpierr)
    call mpi_comm_size(comm, nranks, mpierr)
    call mpi_comm_rank(comm, rank, mpierr)

    rank = rank + 1
    if(rank == 1)then
       write (*,*) "Number of MPI ranks: ", nranks
    end if

  end subroutine parallel_init

  !================================================

  !> Finish the parallel part of the program. Calls mpi_finalize().
  subroutine parallel_finalise()
    call mpi_finalize(mpierr)
  end subroutine parallel_finalise

  !================================================

  !> Stop a parallel model run. Currently simply does
  !! an MPI abort.
  !! @param[in] msg Message to print - reason we're stopping
  subroutine parallel_abort(msg)
    use iso_fortran_env, only : error_unit ! access computing environment
    implicit none
    character(len=*), intent(in) :: msg

    write(error_unit, *) msg
    call mpi_abort(comm, 1, mpierr)

  end subroutine parallel_abort
  
  !================================================

  function get_rank()
    integer :: get_rank
    get_rank = rank
  end function get_rank

  !================================================

  function get_num_ranks() result(num)
    integer :: num
    num = nranks
  end function get_num_ranks
  
  !================================================

  integer function get_max_tag()
    !> Returns maximum available MPI tag value.
    integer :: ierr
    logical :: got

    got = .FALSE.
    call MPI_attr_get(MPI_comm_world,MPI_tag_ub,get_max_tag,got,ierr)
    if ( ierr /= 0 ) call parallel_abort("MPI_attr_get() failed.")

    if ( .not. got ) then
       ! If no value was returned use the minimum possible tag max.
       ! (p. 28 of Version 2.1 of the MPI standard or p. 19 of V.1 of
       ! the standard.)
       get_max_tag = 32767
    endif
  end function get_max_tag
  
  !================================================

  subroutine post_receive(nrecv, source, tag, exch_flag, rbuff, ibuff)
    real(kind=go_wp), dimension(:), optional, intent(inout) :: rbuff
    integer, dimension(:), optional, intent(inout) :: ibuff
    integer, intent(in) :: nrecv
    integer, intent(in) :: tag
    integer, intent(in) :: source
    integer :: exch_flag
    ! Locals
    integer :: ierr

    if(present(rbuff))then
       CALL MPI_irecv (rbuff, nrecv, MPI_DOUBLE_PRECISION, source, &
            tag, MPI_COMM_WORLD, exch_flag, ierr)
    else if(present(ibuff))then
       CALL MPI_irecv (ibuff, nrecv, MPI_INTEGER, source, &
            tag, MPI_COMM_WORLD, exch_flag, ierr)
    else
       call parallel_abort("post_recv: one of ibuf or rbuf must be supplied")
    end if
  end subroutine post_receive
  
  !================================================
  
  subroutine post_send(sendBuff,nsend,destination,tag, &
                       exch_flag)
    integer, intent(in) :: nsend, destination
    real(kind=go_wp), dimension(nsend), intent(in) :: sendBuff
    integer :: tag ! intent is?
    integer :: exch_flag ! intent is??
    ! Locals
    integer :: ierr
    
    call MPI_Isend(sendBuff,nsend,MPI_DOUBLE_PRECISION,     &
                            destination,tag,MPI_COMM_WORLD, &
                            exch_flag,ierr)
  end subroutine post_send
  
  !================================================

  subroutine msg_wait(nmsg, flags, irecv, all)
    integer :: nmsg
    integer, dimension(:) :: flags
    integer :: irecv
    logical, optional,  intent(in) :: all
    ! Locals
    integer :: ierr
    logical :: l_all
    integer :: status(MPI_status_size)
    integer :: astatus(MPI_status_size, nmsg)

    if(present(all))then
       l_all = all
    else
       l_all = .False.
    end if

    if(l_all)then
       CALL MPI_waitall(nmsg, flags, astatus, ierr)
       IF ( ierr.NE.0 ) CALL MPI_abort(MPI_COMM_WORLD,1,ierr)
    else
       call MPI_waitany(nmsg, flags, irecv, status, ierr)
       if ( ierr /= MPI_SUCCESS ) then

          IF(ierr .EQ. MPI_ERR_REQUEST)THEN
             WRITE (*,"(I3,': ERROR: msg_wait: MPI_waitany returned MPI_ERR_REQUEST')") get_rank()
          ELSE IF(ierr .EQ. MPI_ERR_ARG)THEN
             WRITE (*,"(I3,': ERROR: msg_wait: MPI_waitany returned MPI_ERR_ARG')") get_rank()
          ELSE
             WRITE (*,"(I3,': ERROR: msg_wait: MPI_waitany returned unrecognised error')") get_rank()
          END IF
          CALL parallel_abort('exchs_generic: MPI_waitany returned error')
       END IF

    end if
  end subroutine msg_wait

end module parallel_utils_mod
