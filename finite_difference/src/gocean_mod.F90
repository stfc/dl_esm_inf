!> Module containing basic utilities
module gocean_mod
  use iso_c_binding, only: c_intptr_t
  use ocl_utils_mod, only: CL_UTIL_STR_LEN, init_device
  use kind_params_mod
  implicit none

  !> Interface to logging routines
  interface model_write_log
     module procedure write_log_a, write_log_ir, &
                      write_log_i, write_log_r
  end interface

  !> Whether or not this model will run using OpenCL
  logical :: use_opencl

  !> Pointer to the OpenCL device being used (if any)
  integer(c_intptr_t) :: cl_device
  !> The OpenCL context
  integer(c_intptr_t) :: cl_context
  !> Version of OpenCL supported by the device
  character(len=CL_UTIL_STR_LEN) :: version_str

  !> Number of OpenCL command queues
  integer, save :: num_queues
  !> Array of command queues - used to achieve concurrent execution
  integer(c_intptr_t), allocatable, target :: cl_cmd_queues(:)

  !> The OpenCL kernels used by the model
  integer, save :: num_kernels
  character(len=CL_UTIL_STR_LEN), allocatable :: cl_kernel_names(:)
  integer(c_intptr_t), target, allocatable :: cl_kernels(:)

contains

  !===================================================

  !> Initialise the GOcean environment
  subroutine gocean_init(opencl)
    use ocl_utils_mod, only: init_device
#if _OPENACC
    use openacc
#endif
    implicit none
    !> Whether or not to use OpenCL
    logical, optional :: opencl
    integer :: ierr

#if _OPENACC
    call acc_init(acc_device_nvidia)
#endif

    if(present(opencl))then
       use_opencl = opencl
    else
       use_opencl = .False.
    end if

    if(use_opencl)then
       ! Initialise the OpenCL device
       call init_device(cl_device, version_str, cl_context)
       ! Create command queue(s)
       num_queues = 4
       allocate(cl_cmd_queues(num_queues), Stat=ierr)
       if(ierr /= 0)then
          call gocean_stop("Failed to allocate list for OpenCL command queues")
       end if
       call init_cmd_queues(num_queues, cl_cmd_queues, cl_context, cl_device)
    end if

  end subroutine gocean_init

  !===================================================

  subroutine gocean_add_kernels(nkernels, kernel_names, filename)
    use iso_c_binding, only: c_intptr_t
    use ocl_utils_mod, only: get_program, get_kernel, release_program
    integer, intent(in) :: nkernels
    character(len=*), intent(in) :: kernel_names(nkernels)
    character(len=*), intent(in) :: filename
    ! Locals
    integer :: ik, ierr
    integer(c_intptr_t), target :: prog

    ! Get a program object containing all of our kernels
    prog = get_program(cl_context, cl_device, version_str, filename)

    num_kernels = nkernels

    allocate(cl_kernels(num_kernels), cl_kernel_names(num_kernels), &
             Stat=ierr)
    if(ierr /= 0)then
       call gocean_stop("Failed to allocate memory for kernel table")
    end if

    do ik = 1, num_kernels
       cl_kernels(ik) = get_kernel(prog, kernel_names(ik))
    end do

    ! Release the program now that we've created the kernels
    call release_program(prog)

  end subroutine gocean_add_kernels

  !===================================================

  function get_kernel_by_name(name) result(kern)
    integer(c_intptr_t), target :: kern
    character(len=*), intent(in) :: name
    ! Locals
    integer :: ik, match
    character(len=256) :: msg

    !> \TODO is there a better way to do this that reduces the need for
    !! string comparisons?
    match = 0
    do ik = 1, num_kernels
       if(name == cl_kernel_names(ik))then
          ! We can't just return out of this loop because this is a
          ! function
          match = ik
          exit
       end if
    end do

    if(match == 0)then
       !> \TODO add check that we don't go out of bounds when writing to msg
       write(msg, "('get_kernel_by_name: no kernel with name ',(A),' found')")&
            name
       call gocean_stop(msg)
    end if

    kern = cl_kernels(match)

  end function get_kernel_by_name

  !===================================================

  !> Stop the model run. Currently simply does
  !! a Fortran STOP.
  !! @param[in] msg Message to print - reason we're stopping
  subroutine gocean_stop(msg)
    use iso_fortran_env, only : error_unit ! access computing environment
    implicit none
    character(len=*), intent(in) :: msg

    write(error_unit, *) msg
    stop 1

  end subroutine gocean_stop

  !===================================================

  !> Write log entry with one integer and one real arg
  SUBROUTINE write_log_ir(fmtstr, istep, fvar)
    use iso_fortran_env, only : output_unit ! access computing environment
    IMPLICIT none
    CHARACTER(LEN=*), INTENT(in) :: fmtstr
    INTEGER,          INTENT(in) :: istep
    REAL(wp),         INTENT(in) :: fvar

    WRITE(output_unit,FMT=fmtstr) istep, fvar

  END SUBROUTINE write_log_ir

  !===================================================

  !> Write log entry with one integer arg
  SUBROUTINE write_log_i(fmtstr, istep)
    use iso_fortran_env, only : output_unit ! access computing environment
    IMPLICIT none
    CHARACTER(LEN=*), INTENT(in) :: fmtstr
    INTEGER,          INTENT(in) :: istep

    WRITE(output_unit,FMT=fmtstr) istep

  END SUBROUTINE write_log_i

  !===================================================

  !> Write log entry with one real arg
  subroutine write_log_r(fmtstr, fvar)
    use iso_fortran_env, only : output_unit ! access computing environment
    implicit none
    character(len=*), intent(in) :: fmtstr
    real(wp),         intent(in) :: fvar

    write(output_unit,FMT=fmtstr) fvar

  end subroutine write_log_r

  !===================================================

  !> Write log entry with just a string
  subroutine write_log_a(fmtstr, msg)
    use iso_fortran_env, only : output_unit ! access computing environment
    implicit none
    character(len=*), intent(in) :: fmtstr
    character(len=*), intent(in) :: msg

    write(output_unit,FMT=fmtstr) msg

  end subroutine write_log_a

  !===================================================

end module gocean_mod
