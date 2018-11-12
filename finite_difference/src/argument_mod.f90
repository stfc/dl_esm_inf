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


!> Module containing definitions that enable kernel meta-data to be
!! written as valid Fortran.
module argument_mod
use iso_c_binding
use global_parameters_mod
implicit none
private

enum, bind(c) 
   ! The following value is valid for any arg.
   enumerator :: GO_READ
   ! The following values are only valid for fields.
   enumerator :: GO_WRITE, GO_READWRITE, GO_INC
   ! The following values are only valid for globals.
   enumerator :: GO_MIN, GO_MAX, GO_SUM
end enum

type :: go_stencil
   integer :: first_row
   integer :: second_row
   integer :: third_row
end type go_stencil

public go_stencil

  !args(fs,stencil,arg_intent) ! this need defining
type :: go_arg
  integer(kind(GO_READ)) :: arg_intent
  integer :: element
  type(go_stencil) :: stencil_type = go_stencil(0,0,0)
end type go_arg

!-------------------------------------------------------------------------------
! Expose public types
!-------------------------------------------------------------------------------

! Types to enable declarations of elements.
integer, public, parameter :: GO_R_SCALAR=0, GO_I_SCALAR=1
integer, public, parameter :: GO_EVERY=1
! The four types of grid-point on an Arakawa C-grid
integer, public, parameter :: GO_CU=1, GO_CV=2, GO_CT=3, GO_CF=4
! Arguments that a kernel can request that are supported/provided by
! the infrastructure
!> Kernel requires the model time-step
integer, public, parameter :: GO_TIME_STEP   = 1
!> Kernel requires the cell areas of the T-point grid
integer, public, parameter :: GO_GRID_AREA_T = 2
!> Kernel requires the cell areas of the U-point grid
integer, public, parameter :: GO_GRID_AREA_U = 3
!> Kernel requires the cell areas of the V-point grid
integer, public, parameter :: GO_GRID_AREA_V = 4
!> Kernel requires the land/sea mask at T points
integer, public, parameter :: GO_GRID_MASK_T = 5
!> Kernel requires the horizontal grid spacings of the T-point grid
integer, public, parameter :: GO_GRID_DX_T   = 6
!> Kernel requires the horizontal grid spacings of the U-point grid
integer, public, parameter :: GO_GRID_DX_U   = 7
!> Kernel requires the horizontal grid spacings of the V-point grid
integer, public, parameter :: GO_GRID_DX_V   = 8
!> Kernel requires the vertical grid spacings of the T-point grid
integer, public, parameter :: GO_GRID_DY_T   = 9
!> Kernel requires the vertical grid spacings of the U-point grid
integer, public, parameter :: GO_GRID_DY_U   = 10
!> Kernel requires the vertical grid spacings of the V-point grid
integer, public, parameter :: GO_GRID_DY_V   = 11
!> Kernel requires the geographical latitude of U points
integer, public, parameter :: GO_GRID_LAT_U  = 12
!> Kernel requires the geographical latitude of V points
integer, public, parameter :: GO_GRID_LAT_V  = 13
!> Kernel requires the horizontal grid spacing of the grid.
!! Requires/assumes that this quantity is constant.
integer, public, parameter :: GO_GRID_DX_CONST = 14
!> Kernel requires the vertical grid spacing of the grid.
!! Requires/assumes that this quantity is constant.
integer, public, parameter :: GO_GRID_DY_CONST = 15

public :: go_arg
public :: GO_READ, GO_WRITE, GO_READWRITE, GO_INC
public :: GO_SUM, GO_MIN, GO_MAX

!-------------------------------------------------------------------------------
! Member subroutines
!-------------------------------------------------------------------------------
!contains

end module argument_mod
