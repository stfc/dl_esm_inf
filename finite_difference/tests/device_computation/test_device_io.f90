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

! Mock device implementation using separate host memory
module virtual_device
    use field_mod
    use kind_params_mod, only: go_wp
    use iso_c_binding
    implicit none

contains

    subroutine init_device_memory(field)
        type(r2d_field), intent(inout) :: field
        real(go_wp), dimension(:), allocatable, target :: device_memory

        allocate(device_memory(size(field%data,1) * size(field%data,2)))
        field%device_ptr = C_LOC(device_memory)
        field%data_on_device = .true.
        field%read_from_device_f => read_from_device_impl
        field%write_to_device_f => write_to_device_impl
    end subroutine

    subroutine read_from_device_impl(from, to, startx, starty, nx, ny)
        type(c_ptr), intent(in) :: from
        real(go_wp), dimension(:,:), target, intent(inout) :: to
        integer, intent(in) :: startx, starty, nx, ny
        real(go_wp), dimension(:), pointer :: device_memory
        integer :: next_offset, gap, i

        write(*, *) "Read operation", startx, starty, nx, ny
        call C_F_POINTER(from, device_memory, [size(to,1) * size(to,2)])

        next_offset = ((starty - 1) * size(to,1) + (startx-1)) + 1

        do i = starty, starty + ny - 1
            ! Copy next contiguous chunk
            to(startx:startx+nx, i) = device_memory(next_offset:next_offset+nx)
            next_offset = next_offset + size(to,1)
        enddo

    end subroutine read_from_device_impl

    subroutine write_to_device_impl(from, to, startx, starty, nx, ny)
        real(go_wp), dimension(:,:), target, intent(in) :: from
        type(c_ptr), intent(in) :: to
        integer, intent(in) ::  startx, starty, nx, ny
        real(go_wp), dimension(:), pointer :: device_memory
        integer :: i, next_offset, gap

        write(*, *) "Write operation", startx, starty, nx, ny
        call C_F_POINTER(to, device_memory, [size(from,1) * size(from,2)])

        next_offset = ((starty - 1) * size(from,1) + (startx-1)) + 1

        do i = starty, starty + ny - 1
            ! Copy next contiguous chunk
            device_memory(next_offset:next_offset+nx) = from(startx:startx+nx, i)
            next_offset = next_offset + size(from,1)
        enddo

    end subroutine write_to_device_impl

    subroutine simulate_device_computation(buffer, total_size)
        type(c_ptr), intent(in) :: buffer
        integer, intent(in) :: total_size
        real(go_wp), dimension(:), pointer :: device_memory

        write(*, *) "Device computation"
        call C_F_POINTER(buffer, device_memory, [total_size])
        device_memory = device_memory * 2

    end subroutine simulate_device_computation

end module virtual_device

!> Tests for the field read and write to device functionality.
program test_device_io
    use kind_params_mod
    use parallel_mod
    use grid_mod
    use field_mod
    use gocean_mod
    use virtual_device
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

    test_field%data = 0
    call test_field%write_to_device()  ! All device data is 0
    test_field%data = 1
    call test_field%write_to_device(2, 2, 5, 5)  ! A 5x5 block starting at 2,2 is 1
    call simulate_device_computation(test_field%device_ptr, 8*8) ! Double device values
    call test_field%read_from_device(5, 5, 4, 4) ! Read the bottom-right quadrant

    write(*,*) "Resulting array:"
    do i=1, 8
        write(*, "(10f5.1)") test_field%data(:, i)
    enddo

    if ( &
        ! Data not read back must be 1.0
        test_field%data(1, 1) /= 1.0 .or. &
        ! Data written, doubled and read is 2.0
        test_field%data(5, 5) /= 2.0 .or. &
        ! Data not written but read back is still 0.0
        test_field%data(7, 7) /= 0.0 ) then
        stop "Error - Results are incorrect"
    endif

    write(*,*) "Test passed successfuly! Note that this test should still ", &
               "succeed when changing the DL_ESM_ALIGNMENT environment variable."
    call gocean_finalise()

end program test_device_io
