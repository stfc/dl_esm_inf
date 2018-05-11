!------------------------------------------------------------------------------
! BSD 2-Clause License
! 
! Copyright (c) 2018, Science and Technology Facilities Council
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

!> A simple example of the expected usage of the dl_esm_inf library in
!! constructing a finite-difference model.
program model
  use kind_params_mod
  use parallel_mod
  use subdomain_mod
  use grid_mod
  use field_mod
  use gocean_mod
  implicit none
  ! Total size of the model domain
  integer :: jpiglo = 100, jpjglo = 100
  ! (Uniform) grid resolution
  real(wp) :: dx = 1.0
  real(wp) :: dy = 1.0
  !> Our domain decomposition
  type(decomposition_type) :: decomp
  !> The grid on which our fields are defined
  type(grid_type), target :: model_grid
  !> An example field
  type(r2d_field) :: a_field
  ! Local definition of the T-point mask which defines whether T points are
  ! wet (1), dry (0) or outside (-1) the simulation domain.
  integer, allocatable :: tmask(:,:)
  integer :: my_rank
  integer :: ierr

  call gocean_initialise()

  !> Create our grid object
  model_grid = grid_type(ARAKAWA_C, &
                         (/BC_EXTERNAL,BC_EXTERNAL,BC_NONE/), &
                         OFFSET_NE)

  !> Generate a domain decomposition
  decomp = decompose(jpiglo, jpjglo)
  my_rank = get_rank()

  !> Create a T-point mask describing the (local) domain
  allocate(tmask(decomp%subdomains(my_rank)%global%nx,  &
                 decomp%subdomains(my_rank)%global%ny), Stat=ierr)
  if(ierr /= 0)then
     call gocean_stop('Failed to allocate T-mask')
  end if
  ! To keep things simple We set all points to be wet and within the domain
  tmask = 1

  !> Complete the initialisation of the grid using the T-mask and
  !! grid resolution
  call grid_init(model_grid, decomp, dx, dy, tmask)

  !> Create a field on the grid
  a_field = r2d_field(model_grid, U_POINTS)

  ! All done!
  if (my_rank == 1) write(*,'(/"Example model set-up complete."/)')

  call gocean_finalise()

end program model
