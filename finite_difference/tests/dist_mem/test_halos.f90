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

!> A simple example of the expected usage of the dl_esm_inf library in
!! constructing a finite-difference model.
program test_halos
  use kind_params_mod
  use parallel_mod
  use decomposition_mod, only: decomposition_type
  use subdomain_mod
  use grid_mod
  use field_mod
  use gocean_mod
  implicit none
  ! Total size of the model domain
  integer :: jpiglo = 0, jpjglo = 0
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

  !> Generate a domain decomposition
  decomp = decompose(jpiglo, jpjglo)
  my_rank = get_rank()

  if(my_rank == 1) then
     write(*,"('Using global domain size of ',I4,'x',I4)") jpiglo, jpjglo
  end if

  !> Create a T-point mask describing the (local) domain
  allocate(tmask(decomp%subdomains(my_rank)%global%nx,  &
                 decomp%subdomains(my_rank)%global%ny), Stat=ierr)
  if(ierr /= 0)then
     call gocean_stop('Failed to allocate T-mask')
  end if
  ! To keep things simple we set all points to be wet and within the domain
  tmask = 1

  !> Complete the initialisation of the grid using the T-mask and
  !! grid resolution
  call grid_init(model_grid, decomp, dx, dy, tmask)

  !> \TODO put these inside library initialisation
  call map_comms(decomp, tmask, .false., ierr)
  
  !> Create fields on U,V,T,F-points of the grid
  u_field = r2d_field(model_grid, GO_U_POINTS)
  v_field = r2d_field(model_grid, GO_V_POINTS)
  t_field = r2d_field(model_grid, GO_T_POINTS)
  f_field = r2d_field(model_grid, GO_F_POINTS)

  call init_field_hill(u_field)
  call init_field_hill(v_field)
  call init_field_hill(t_field)
  call init_field_hill(f_field)

  write(*,"(I3,' Halo exchange for u:')") my_rank
  call u_field%halo_exch(1)
  write(*,"(I3,' Halo exchange for v:')") my_rank
  call v_field%halo_exch(1)
  write(*,"(I3,' Halo exchange for t:')") my_rank
  call t_field%halo_exch(1)
  write(*,"(I3,' Halo exchange for f:')") my_rank
  call f_field%halo_exch(1)

  !call dump_field(u_field, "u_fld", halo_depth=1)
  !call dump_field(v_field, "v_fld", halo_depth=1)
  !call dump_field(f_field, "f_fld", halo_depth=1)

  ! Check the halo regions
  call check_hill_halos(t_field, "t_fld")
  call check_hill_halos(u_field, "u_fld")
  call check_hill_halos(v_field, "v_fld")
  call check_hill_halos(f_field, "f_fld")
  
  call gocean_finalise()

contains

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
          field%data(ji,jj) = hill(field, ji, jj)
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

  function hill(field, ji, jj)
    !> Utility function that returns a unique value constructed
    !! from the global location corresponding to the local ji,jj
    !! point.
    use field_mod, only: r2d_field
    type(r2d_field), intent(in) :: field
    integer, intent(in) :: ji, jj
    real(go_wp) :: hill, xpos, ypos
    ! Global position of the T point
    xpos = field%grid%xt(ji,jj)
    ypos = field%grid%yt(ji,jj)
    select case(field%grid%offset)
    case(GO_OFFSET_SW)
       select case(field%defined_on)
       case(GO_U_POINTS)
          xpos = xpos - 0.5*field%grid%dx_u(ji, jj)
       case(GO_F_POINTS)
          xpos = xpos - 0.5*field%grid%dx_f(ji, jj)
          ypos = ypos - 0.5*field%grid%dy_f(ji, jj)
       case(GO_V_POINTS)
          ypos = ypos - 0.5*field%grid%dy_v(ji, jj)
       end select
    case(GO_OFFSET_NE)
       select case(field%defined_on)
       case(GO_U_POINTS)
          xpos = xpos + 0.5*field%grid%dx_u(ji, jj)
       case(GO_F_POINTS)
          xpos = xpos + 0.5*field%grid%dx_f(ji, jj)
          ypos = ypos + 0.5*field%grid%dy_f(ji, jj)
       case(GO_V_POINTS)
          ypos = ypos + 0.5*field%grid%dy_v(ji, jj)
       end select
    case default
       call gocean_stop("hill: unsupported grid offset")
    end select
    hill = real(10000.0*xpos + ypos)
  end function hill
  
  subroutine check_hill_halos(field, name)
    use parallel_mod, only: get_rank
    type(r2d_field), intent(in) :: field
    character(len=*), intent(in) :: name
    ! Locals
    integer :: ji, jj, iglobal, jglobal
    real(go_wp), parameter :: TOL_ZERO = 1.0E-8
    integer :: my_rank

    my_rank = get_rank()
    
    if(field%grid%subdomain%global%xstart > 1)then
       write(*,"(I4,' Subdomain has a halo region in the -x direction')") &
            my_rank
       ! Local and global x location of the L1 halo
       ji = field%internal%xstart - 1
       iglobal = field%grid%subdomain%global%xstart - 1
       do jj = field%internal%ystart, field%internal%ystop
          jglobal = field%grid%subdomain%global%ystart + jj - &
               field%internal%ystart 
          if( abs(field%data(ji, jj) - hill(field, ji, jj)) > TOL_ZERO)then
             write(*,'(I4," ERROR in halo value for ",(A),": ",4I5,2F9.1)') &
                  my_rank, TRIM(name), ji, jj, iglobal, jglobal, &
                  field%data(ji,jj), hill(field, ji, jj)
          end if
       end do
    end if
    if(field%grid%subdomain%global%xstop < field%grid%global_nx)then
       write(*,"(I4,' Subdomain has a halo region in the +x direction')") &
            my_rank
       ! Local and global x location of the L1 halo
       ji = field%internal%xstop + 1
       iglobal = field%grid%subdomain%global%xstop + 1
       do jj = field%internal%ystart, field%internal%ystop
          jglobal = field%grid%subdomain%global%ystart + jj - &
               field%internal%ystart 
          if( abs(field%data(ji, jj) - hill(field, ji, jj)) > TOL_ZERO)then
             write(*,'(I4," ERROR in halo value for ",(A),": ",4I5,2F9.1)') &
                  my_rank, TRIM(name), ji, jj, iglobal, jglobal, &
                  field%data(ji,jj), hill(field, ji, jj)
          end if
       end do
    end if
    if(field%grid%subdomain%global%ystart > 1)then
       write(*,"(I4,' Subdomain has a halo region in the -y direction')") &
            my_rank
       jj = field%internal%ystart - 1
       jglobal = field%grid%subdomain%global%ystart - 1
       do ji = field%internal%xstart, field%internal%xstop
          iglobal = field%grid%subdomain%global%xstart + ji - &
               field%internal%xstart
          if( abs(field%data(ji, jj) - hill(field, ji, jj)) > TOL_ZERO)then
             write(*,'(I4," ERROR in halo value for ",(A),": ",4I5,2F9.1)') &
                  my_rank, TRIM(name), ji, jj, iglobal, jglobal, &
                  field%data(ji,jj), hill(field, ji, jj)
          end if
       end do
    end if
    if(field%grid%subdomain%global%ystop < field%grid%global_ny)then
       write(*,"(I4,' Subdomain has a halo region in the +y direction')") &
            my_rank
       jj = field%internal%ystop + 1
       jglobal = field%grid%subdomain%global%ystop + 1
       do ji = field%internal%xstart, field%internal%xstop
          iglobal = field%grid%subdomain%global%xstart + ji - &
               field%internal%xstart
          if( abs(field%data(ji, jj) - hill(field, ji, jj)) > TOL_ZERO)then
             write(*,'(I4," ERROR in halo value for ",(A),": ",4I5,2F9.1)') &
                  my_rank, TRIM(name), ji, jj, iglobal, jglobal, &
                  field%data(ji,jj), hill(field, ji, jj)
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
    real :: xpos, ypos
    
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
    ! Global position of the T point
    xpos = field%grid%xt(ji,jj)
    ypos = field%grid%yt(ji,jj)
    select case(field%grid%offset)
    case(GO_OFFSET_SW)
       select case(field%defined_on)
       case(GO_U_POINTS)
          xpos = xpos - 0.5*field%grid%dx_u(ji, jj)
       case(GO_F_POINTS)
          xpos = xpos - 0.5*field%grid%dx_f(ji, jj)
          ypos = ypos - 0.5*field%grid%dy_f(ji, jj)
       case(GO_V_POINTS)
          ypos = ypos - 0.5*field%grid%dy_v(ji, jj)
       end select
    case(GO_OFFSET_NE)
       select case(field%defined_on)
       case(GO_U_POINTS)
          xpos = xpos + 0.5*field%grid%dx_u(ji, jj)
       case(GO_F_POINTS)
          xpos = xpos + 0.5*field%grid%dx_f(ji, jj)
          ypos = ypos + 0.5*field%grid%dy_f(ji, jj)
       case(GO_V_POINTS)
          ypos = ypos + 0.5*field%grid%dy_v(ji, jj)
       end select
    case default
       call gocean_stop("dump_field: unsupported grid offset")
    end select

          write(unit_no,'(3e16.7e3)') xpos, &
                                      ypos, &
                                      field%data(ji,jj)
       end do
       ! Blank line at end of scan for gnuplot
       write(unit_no,*)
    end do
    close(unit_no)
       
  end subroutine dump_field
end program test_halos
