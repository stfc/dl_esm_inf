!------------------------------------------------------------------------------
! BSD 2-Clause License
! 
! Copyright (c) 2021, Science and Technology Facilities Council
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
! Author: S. Siso, STFC Daresbury Laboratory

!> Tests for the field read and write to device functionality.
program test_device_io
    use kind_params_mod
    use parallel_mod
    use grid_mod
    use field_mod
    use gocean_mod
    implicit none
    ! Total size of the model domain
    integer :: jpiglo = 5, jpjglo = 5
    ! (Uniform) grid resolution
    real(go_wp) :: dx = 1.0
    real(go_wp) :: dy = 1.0
    !> The grid on which our fields are defined
    type(grid_type), target :: model_grid
    !> An example field
    type(r2d_field) :: test_field
    ! Local definition of the T-point mask which defines whether T points are
    ! wet (1), dry (0) or outside (-1) the simulation domain.
    integer, allocatable :: tmask(:,:)
    integer :: ierr
    integer :: i

    call gocean_initialise()

    !> Create our grid object
    model_grid = grid_type(GO_ARAKAWA_C, &
                           (/GO_BC_EXTERNAL,GO_BC_EXTERNAL,GO_BC_NONE/), &
                           GO_OFFSET_NE)

    !> Generate a domain decomposition. This automatically uses the number
    !! of MPI ranks available.
    call model_grid%decompose(jpiglo, jpjglo)

    !> Create a T-point mask describing the (local) domain
    allocate(tmask(model_grid%subdomain%global%nx,  &
             model_grid%subdomain%global%ny), Stat=ierr)
    if(ierr /= 0)then
        call gocean_stop('Failed to allocate T-mask')
    end if
    ! To keep things simple we set all points to be wet and within the domain
    tmask = 1

    !> Complete the initialisation of the grid using the T-mask and
    !! grid resolution
    call grid_init(model_grid, dx, dy, tmask)

    !> Create fields on U,V,T,F-points of the grid
    test_field = r2d_field(model_grid, GO_U_POINTS)
    call init_device_memory(test_field)

    test_field%data = 1
    call test_field%write_to_device()

    do i=1, 8
    write(*, "(10f5.1)") test_field%data(i, :)
    enddo

    call gocean_finalise()

contains

    ! Mock device implementation using separate host memory
    subroutine init_device_memory(field)
        type(r2d_field), intent(inout) :: field
        real, dimension(:), allocatable, target :: device_memory

        allocate(device_memory(size(field%data,1) * size(field%data,2)))
        field%device_ptr = C_LOC(device_memory)
        field%read_from_device_f => read_from_device_impl
        field%write_to_device_f => write_to_device_impl

    end subroutine

    subroutine read_from_device_impl(from, to, offset, nx, ny, stride_gap)
        use iso_c_binding, only: c_intptr_t, c_int
        use kind_params_mod, only: go_wp
        integer(c_intptr_t), intent(in) :: from
        real(go_wp), dimension(:,:), intent(inout) :: to
        integer, intent(in) :: offset, nx, ny, stride_gap

    end subroutine read_from_device_impl

    subroutine write_to_device_impl(from, to, offset, nx, ny, stride_gap)
        use iso_c_binding, only: c_intptr_t, c_int
        use kind_params_mod, only: go_wp
        real(go_wp), dimension(:,:), intent(in) :: from
        integer(c_intptr_t), intent(in) :: to
        integer, intent(in) :: offset, nx, ny, stride_gap
        real, dimension(:), pointer :: device_memory
        integer :: i, startx, starty

        call C_F_POINTER(to, device_memory)
        
        startx = 1
        starty = 1

        do i = starty, ny
            device_memory(offset+startx:offset+startx+nx) = from(i,startx:startx+nx) 
        enddo

    end subroutine write_to_device_impl

end program test_device_io
