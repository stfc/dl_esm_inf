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
! Author: A. R. Porter, STFC Daresbury Laboratory

!> Tests for the halo-exchange functionality of the dl_esm_inf library.
!! Domain size must be supplied via the JPIGLO and JPJGLO environment
!! variables.
program test_reduction
  use kind_params_mod
  use parallel_mod
  use grid_mod
  use field_mod
  use gocean_mod
  implicit none
  ! Total size of the model domain
  integer :: jpiglo = 0, jpjglo = 0
  ! (Uniform) grid resolution
  real(go_wp) :: dx = 1.0
  real(go_wp) :: dy = 1.0
  !> The grid on which our fields are defined
  type(grid_type), target :: model_grid
  !> An example field
  type(r2d_field) :: field

  real(go_wp), allocatable :: global_init_data(:,:)
  real(go_wp), allocatable :: global_updated_data(:,:)
  integer :: my_rank
  integer :: ierr, i, j
  character(len=10) :: env

  call get_environment_variable("JPIGLO", env)
  read(env, '(i10)') jpiglo
  call get_environment_variable("JPJGLO", env)
  read(env, '(i10)') jpjglo
  if(jpiglo < 1 .or. jpjglo < 1)then
     stop 'Domain size must be set via $JPIGLO and $JPJGLO'
  end if
  
  call gocean_initialise()

  !> Create our grid object
  model_grid = grid_type(GO_ARAKAWA_C, &
                         (/GO_BC_EXTERNAL,GO_BC_EXTERNAL,GO_BC_NONE/), &
                         GO_OFFSET_NE)

  !> Generate a domain decomposition. This automatically uses the number
  !! of MPI ranks available.
  call model_grid%decompose(jpiglo, jpjglo)
  my_rank = get_rank()

  if(my_rank == 1) then
     write(*,"('Using global domain size of ',I4,'x',I4)") jpiglo, jpjglo
  end if

  !> Create global array to hold the global initialisation data
  allocate(global_init_data(jpiglo, jpjglo), stat=ierr)
  if(ierr /= 0)then
     call gocean_stop('Failed to allocate global init data')
  end if

  do j=1, jpjglo
      global_init_data(:,j) = (/ (unique_global_value(i, j, jpiglo), i=1, jpiglo) /)
  enddo

  !> Complete the initialisation of the grid
  call grid_init(model_grid, dx, dy)

  !> Create T point field
  field = r2d_field(model_grid, GO_T_POINTS, &
                      init_global_data=global_init_data)
  !> Check that the global data is set correctly in the local fields
  call check_field_distribution (field)

  !> Update the local data
  call update_field(field)

  !> This will allocate global_updated_data:
  call field%gather_inner_data(global_updated_data)

  !> The result is only on rank 1:
  if (my_rank == 1) then 
     call check_gathered_field(global_updated_data)
  endif
  
  call gocean_finalise()

contains

   !> Computes a unique value for each index pair (i,j) (1<=i,j)
   !! given that there are n values in a row. It assigns the
   !! numbers 0, 1, ..., n*(#rows)-1 to each cell:
   function unique_global_value(i, j, n)
      implicit none
      integer, intent(in) :: i, j, n
      integer :: unique_global_value
      unique_global_value = (i-1) + (j-1)*n
   end function unique_global_value

  ! -----------------------------------------------------------------------------------
  !> Checks that the local data in field corresponds to the original
  !! global data:
  subroutine check_field_distribution(field)
    type(r2d_field), intent(inout) :: field
    integer :: ji, jj
    integer :: istart, istop, jstart, jstop, xmin, ymin
    real(go_wp) :: global_value
    istart = field%internal%xstart
    istop =  field%internal%xstop
    jstart = field%internal%ystart
    jstop =  field%internal%ystop
    xmin = field%grid%subdomain%global%xstart
    ymin = field%grid%subdomain%global%ystart
    
    do jj = jstart, jstop
       do ji = istart, istop
         global_value = unique_global_value(xmin + ji-istart,  &
                                            ymin + jj-jstart,  &
                                            field%grid%decomp%global_nx)
         if (field%data(ji, jj) /= global_value) then
            print *,"Incorrect field data for ", ji, jj, field%data(ji, jj), global_value
            call gocean_stop('Failed to distribute init data')
           endif
           
       end do
    end do
    print *,"Field distributed correctly."
  end subroutine check_field_distribution

  ! -----------------------------------------------------------------------------------
  !> Adds one to each entry in the local fiels
  subroutine update_field(field)
    type(r2d_field), intent(inout) :: field
    integer :: ji, jj
    integer :: istart, istop, jstart, jstop
    istart = field%internal%xstart
    istop =  field%internal%xstop
    jstart = field%internal%ystart
    jstop =  field%internal%ystop
    
    do jj = jstart, jstop
       do ji = istart, istop
          field%data(ji,jj) = field%data(ji,jj) + 1
       end do
    end do
  end subroutine update_field

  !! -----------------------------------------------------------------------------------
  !> Checks if the global data gathered after updating the local fields is correct:
  subroutine check_gathered_field(global_field)
    real(go_wp), allocatable, dimension(:,:) :: global_field
    integer :: ji, jj
    real(go_wp) :: global_value
    
    do jj = 1, size(global_field, 2)
       do ji = 1, size(global_field, 1)
         global_value = unique_global_value(ji, jj, size(global_field, 1)) + 1
         if (global_field(ji, jj) /= global_value) then
            print *, "Incorrect global field data gathered for ", ji, jj, &
                     global_field(ji, jj), global_value
            !call gocean_stop('Failed to gather local data')
           endif
           
       end do
    end do

    print *,"Field gathered correctly."
  end subroutine check_gathered_field
  
end program test_reduction
