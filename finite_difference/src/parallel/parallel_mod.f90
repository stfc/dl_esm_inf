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
  !> MPI rank + 1 of current process
  integer :: rank
  !> Total no. of MPI processes
  integer :: nranks

  public parallel_init, parallel_finalise, parallel_abort
  public get_rank

contains

  !================================================

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

  subroutine parallel_finalise()
    call mpi_finalize(mpierr)
  end subroutine parallel_finalise

  !================================================

  function get_rank()
    integer :: get_rank
    get_rank = rank
  end function get_rank

  !================================================

  function get_num_ranks()
    integer :: get_num_ranks
    get_num_ranks = nranks
  end function get_rank

  !================================================

  !> Stop a parallel model run. Currently simply does
  !! an MPI abort.
  !! @param[in] msg Message to print - reason we're stopping
  subroutine parallel_abort(msg)
    use iso_fortran_env, only : error_unit ! access computing environment
    implicit none
    character(len=*), intent(in) :: msg
    integer :: ierr

    write(error_unit, *) msg
    call mpi_abort(comm, 1, ierr)

  end subroutine parallel_abort

end module parallel_mod
