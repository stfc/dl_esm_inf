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
  use decomposition_mod, only: decomposition_type
  use subdomain_mod
  use grid_mod
  use field_mod
  use gocean_mod
  implicit none
  ! Total size of the model domain
  integer :: jpiglo = 100, jpjglo = 50
  ! (Uniform) grid resolution
  real(go_wp) :: dx = 1.0
  real(go_wp) :: dy = 1.0
  !> Our domain decomposition
  type(decomposition_type) :: decomp
  !> The grid on which our fields are defined
  type(grid_type), target :: model_grid
  !> An example field
  type(r2d_field) :: u_field, v_field, t_field, f_field
  ! Local definition of the T-point mask which defines whether T points are
  ! wet (1), dry (0) or outside (-1) the simulation domain.
  integer, allocatable :: tmask(:,:)
  integer :: my_rank
  integer :: ierr

  call gocean_initialise()

  !> Create our grid object
  model_grid = grid_type(GO_ARAKAWA_C, &
                         (/GO_BC_EXTERNAL,GO_BC_EXTERNAL,GO_BC_NONE/), &
                         GO_OFFSET_NE)

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

  !> \TODO put these inside library initialisation
  call map_comms(decomp, tmask, .false., ierr)
  
  !> Create a field on U-points of the grid
  u_field = r2d_field(model_grid, GO_U_POINTS)
  v_field = r2d_field(model_grid, GO_V_POINTS)
  t_field = r2d_field(model_grid, GO_T_POINTS)
  f_field = r2d_field(model_grid, GO_F_POINTS)

  call init_field_hill(u_field)
  call init_field_hill(v_field)
  call init_field_hill(t_field)
  call init_field_hill(f_field)

  write(*,*) "Halo exchange for u:"
  call u_field%halo_exch(1)
  write(*,*) "Halo exchange for v:"
  call v_field%halo_exch(1)
  !call t_field%halo_exch(1)
  !call f_field%halo_exch(1)

  call dump_field(u_field, "u_fld", halo_depth=1)
  call dump_field(v_field, "v_fld")

  ! Check the halo regions
  call check_hill_halos(u_field)
  
  ! All done!
  if (my_rank == 1) write(*,'(/"Example model set-up complete."/)')

  call gocean_finalise()

contains

  subroutine init_field_by_rank(field)
    !> Initialise a field with the MPI rank of this process
    use parallel_mod, only: get_rank
    type(r2d_field), intent(inout) :: field
    ! Locals
    integer :: my_rank
    my_rank = get_rank()
    field%data(:,:) = real(my_rank)
    
  end subroutine init_field_by_rank

  subroutine init_field_hill(field)
    type(r2d_field), intent(inout) :: field
    integer :: ji, jj
    integer :: istart, istop, jstart, jstop
    istart = field%internal%xstart
    istop =  field%internal%xstop
    jstart = field%internal%ystart
    jstop =  field%internal%ystop
    
    field%data(:,:) = 0.0
    do jj = jstart, jstop
       do ji = istart, istop
          field%data(ji,jj) = hill(field%grid%subdomain%global%xstart +ji-1 , &
                                   field%grid%subdomain%global%ystart +jj-1)
       end do
    end do
    ! Set external points to plausible but wrong values
    do jj = 1, field%grid%ny
       field%data(1:istart-1,jj) = field%data(istart,jj)
       field%data(istop+1:,jj) = field%data(istop,jj)
    end do
    do ji = 1, field%grid%nx
       field%data(ji,1:jstart-1) = field%data(ji,jstart)
       field%data(ji,jstop+1:) = field%data(ji,jstop)
    end do
  end subroutine init_field_hill

  function hill(ji, jj)
    integer, intent(in) :: ji, jj
    real(go_wp) :: hill
    hill = real(ji + jj)
  end function hill
  
  subroutine check_hill_halos(field)
    type(r2d_field), intent(in) :: field
    ! Locals
    integer :: ji, jj, iglobal, jglobal
    real(go_wp), parameter :: TOL_ZERO = 1.0E-8

    if(field%grid%subdomain%global%xstart > 1)then
       ! This subdomain has a halo region in the -x direction
       ji = field%internal%xstart - 1
       do jj = field%internal%ystart, field%internal%ystop
          iglobal = field%grid%subdomain%global%xstart + ji - 1
          jglobal = field%grid%subdomain%global%ystart + jj - 1 
          if( abs(field%data(ji, jj) - hill(iglobal, jglobal)) > TOL_ZERO)then
             write(*,*) "ERROR in halo value: ", ji, jj, iglobal, jglobal, &
                  field%data(ji,jj), hill(iglobal,jglobal)
          end if
       end do
    end if
    
  end subroutine check_hill_halos
  
  subroutine dump_field(field, name, halo_depth)
    use parallel_mod, only: get_rank
    type(r2d_field), intent(in) :: field
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: halo_depth
    ! Locals
    integer, parameter :: unit_no = 21
    integer :: ji, jj
    character(len=5) :: rank_str
    integer :: lhalo

    ! Each rank creates its own file
    my_rank = get_rank()
    write(rank_str, "(I5.5)") my_rank
    open(unit=unit_no, file=name//"_"//TRIM(rank_str)//".dat", &
         status="unknown", action="write")

    if(present(halo_depth))then
       lhalo = halo_depth
       if(lhalo >= field%internal%xstart .or. &
            lhalo >= field%internal%ystart)then
          write (*,*) "field_dump: error: requested halo depth is &
               & greater than the halo depth of supplied field - not &
               & writing data"
          return
       end if
    else
       lhalo = 0
    end if
    
    do jj = field%internal%ystart-lhalo, field%internal%ystop+lhalo, 1
       do ji = field%internal%xstart-lhalo, field%internal%xstop+lhalo, 1

          write(unit_no,'(3e16.7e3)') field%grid%xt(ji,jj), &
                                      field%grid%yt(ji,jj), &
                                      field%data(ji,jj)
       end do
       ! Blank line at end of scan for gnuplot
       write(unit_no,*)
    end do
    close(unit_no)
       
  end subroutine dump_field
end program model
