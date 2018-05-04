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
module parallel_mod
  use mpi
  use parallel_common_mod
  implicit none

  private

  !> Error flag for MPI calls. Note that this module makes no attempt to
  !! check these since we assume that MPI will abort if an error occurs
  !! (which is the default behaviour on the majority of systems).
  integer :: mpierr

  !> Our MPI communicator
  integer :: comm

  public parallel_init, parallel_finalise, parallel_abort
  !> Make certain routines from parallel_common available to
  !! USE'rs of this module.
  public get_rank, get_num_ranks
  public set_proc_grid, get_proc_grid

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

    ! Initialise to meaningless values since we don't
    ! yet know what our processor grid will look like.
    nprocx = 0
    nprocy = 0

  end subroutine parallel_init

  !================================================

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

end module parallel_mod
