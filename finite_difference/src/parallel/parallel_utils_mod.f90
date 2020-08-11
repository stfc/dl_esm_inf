!------------------------------------------------------------------------------
! BSD 2-Clause License
! 
! Copyright (c) 2018-2020, Science and Technology Facilities Council.
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
  !> \todo #39. Ideally we would name the various symbols with an 'only' clause
  !! but this causes compilation errors with MPICH.
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

  !> Copy some MPI parameters into our namespace (because they're used
  !! when looping over messages in the layer above this one)
  integer, parameter :: MSG_UNDEFINED = MPI_UNDEFINED
  integer, parameter :: MSG_REQUEST_NULL = MPI_REQUEST_NULL

  public parallel_init, parallel_finalise, parallel_abort, get_max_tag
  public get_rank, get_num_ranks, post_receive, post_send, global_sum, on_master
  public msg_wait, msg_wait_all
  public MSG_UNDEFINED, MSG_REQUEST_NULL, DIST_MEM_ENABLED

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
       write (*,"('Number of MPI ranks:', I4)") nranks
       write (*,*)
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

  function on_master()
    logical :: on_master
    on_master = get_rank() == 1
  end function on_master

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
    integer, intent(in)  :: nrecv
    integer, intent(in)  :: tag
    integer, intent(in)  :: source
    integer, intent(out) :: exch_flag
    ! Locals
    integer :: ierr

    if(present(rbuff))then
       call MPI_irecv (rbuff, nrecv, MPI_DOUBLE_PRECISION, source, &
            tag, MPI_COMM_WORLD, exch_flag, ierr)
    else if(present(ibuff))then
       call MPI_irecv (ibuff, nrecv, MPI_INTEGER, source, &
            tag, MPI_COMM_WORLD, exch_flag, ierr)
    else
       call parallel_abort("post_recv: one of ibuf or rbuf must be supplied")
    end if
  end subroutine post_receive
  
  !================================================
  
  subroutine post_send(sendBuff, nsend, destination, tag, &
                       exch_flag)
    integer, intent(in) :: nsend, destination
    real(kind=go_wp), dimension(nsend), intent(in) :: sendBuff
    integer, intent(in)  :: tag
    integer, intent(out) :: exch_flag
    ! Locals
    integer :: ierr
    
    call MPI_Isend(sendBuff, nsend, MPI_DOUBLE_PRECISION, &
                   destination, tag, MPI_COMM_WORLD,      &
                   exch_flag, ierr)
  end subroutine post_send
  
  !================================================

  subroutine msg_wait(nmsg, flags, irecv)
    !> Number of requests to wait on
    integer, intent(in) :: nmsg
    !> List of requests
    integer, dimension(:), intent(inout) :: flags
    !> Which message we've received
    integer, intent(out) :: irecv
    ! Locals
    integer :: ierr
    integer :: status(MPI_status_size)
    integer :: astatus(MPI_status_size, nmsg)

    call MPI_waitany(nmsg, flags, irecv, status, ierr)

    if ( ierr /= MPI_SUCCESS ) then
       if(ierr == MPI_ERR_REQUEST)then
          write (*,"(I3,': ERROR: msg_wait: MPI_waitany returned MPI_ERR_REQUEST')") get_rank()
       else if(ierr == MPI_ERR_ARG)then
          write (*,"(I3,': ERROR: msg_wait: MPI_waitany returned MPI_ERR_ARG')") get_rank()
       else
          write (*,"(I3,': ERROR: msg_wait: MPI_waitany returned unrecognised error')") get_rank()
       end if
       call parallel_abort('exchs_generic: MPI_waitany returned error')
    end if

  end subroutine msg_wait
  
  !================================================

  subroutine msg_wait_all(nmsg, flags)
    !> Number of messages to wait on
    integer, intent(in)                  :: nmsg
    integer, dimension(:), intent(inout) :: flags
    ! Locals
    integer :: ierr
    integer :: astatus(MPI_status_size, nmsg)

    call MPI_waitall(nmsg, flags, astatus, ierr)
    if ( ierr /= MPI_SUCCESS ) call MPI_abort(MPI_COMM_WORLD, 1, ierr)

  end subroutine msg_wait_all

  !================================================

  subroutine global_sum(var)
    !> Performs a global sum on a single, double-precision scalar.
    real(go_wp), intent(inout) :: var
    ! Locals
    integer :: ierr

    call MPI_allreduce(MPI_IN_PLACE, var, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)

  end subroutine global_sum

end module parallel_utils_mod
