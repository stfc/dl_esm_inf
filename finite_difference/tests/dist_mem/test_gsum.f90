!------------------------------------------------------------------------------
! BSD 2-Clause License
! 
! Copyright (c) 2020, Science and Technology Facilities Council.
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

!> Tests for the global-sum functionality of the dl_esm_inf library.
!! Domain size must be supplied via the JPIGLO and JPJGLO environment
!! variables.
program test_gsum
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
  type(r2d_field) :: u_field, v_field, t_field, f_field
  ! Local definition of the T-point mask which defines whether T points are
  ! wet (1), dry (0) or outside (-1) the simulation domain.
  integer, allocatable :: tmask(:,:)
  integer :: my_rank
  integer :: ierr
  character(len=10) :: env
  real(go_wp) :: expected_sum, actual_sum
  real(go_wp) :: TOL_ZERO = 1.0D-8

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
  u_field = r2d_field(model_grid, GO_U_POINTS)
  v_field = r2d_field(model_grid, GO_V_POINTS)
  t_field = r2d_field(model_grid, GO_T_POINTS)
  f_field = r2d_field(model_grid, GO_F_POINTS)

  !> Set internal field points to 1.0, external to -100.0
  call init_field(u_field)
  call init_field(v_field)
  call init_field(t_field)
  call init_field(f_field)

  ! Global sum should be equal to the number of points in the domain
  expected_sum = jpiglo*jpjglo
  actual_sum = field_checksum(u_field)
  if(ABS(expected_sum - actual_sum) > TOL_ZERO)then
     call gocean_stop('Checksum of U field incorrect!')
  else
     write(*,"(I3,' Global sum for u: ', E15.8)") my_rank, actual_sum
  end if
  actual_sum = field_checksum(v_field)
  if(ABS(expected_sum - actual_sum) > TOL_ZERO)then
     call gocean_stop('Checksum of V field incorrect!')
  else
     write(*,"(I3,' Global sum for v: ', E15.8)") my_rank, actual_sum
  end if
  actual_sum = field_checksum(t_field)
  if(ABS(expected_sum - actual_sum) > TOL_ZERO)then
     call gocean_stop('Checksum of T field incorrect!')
  else
     write(*,"(I3,' Global sum for t: ', E15.8)") my_rank, actual_sum
  end if
  actual_sum = field_checksum(f_field)
  if(ABS(expected_sum - actual_sum) > TOL_ZERO)then
     call gocean_stop('Checksum of F field incorrect!')
  else
     write(*,"(I3,' Global sum for t: ', E15.8)") my_rank, actual_sum
  end if

  call gocean_finalise()

contains

  subroutine init_field(field)
    type(r2d_field), intent(inout) :: field
    integer :: ji, jj
    integer :: istart, istop, jstart, jstop
    istart = field%internal%xstart
    istop =  field%internal%xstop
    jstart = field%internal%ystart
    jstop =  field%internal%ystop
    
    field%data(:,:) = -100.0
    do jj = jstart, jstop
       do ji = istart, istop
          field%data(ji,jj) = 1.0
       end do
    end do

  end subroutine init_field

end program test_gsum
