!> Module containing basic utilities
module gocean_mod
  use kind_params_mod
  implicit none

  private

  !> Interface to logging routines
  interface model_write_log
     module procedure write_log_a, write_log_ir, &
                      write_log_i, write_log_r
  end interface

  public gocean_initialise, gocean_finalise, gocean_stop
  public model_write_log

contains

  !===================================================

  !> Initialise the GOcean environment
  subroutine gocean_initialise()
#if _OPENACC
    use openacc
#endif
    use parallel_mod, only: parallel_init
    implicit none

    call parallel_init()

#if _OPENACC
    call acc_init(acc_device_nvidia)
#endif
  end subroutine gocean_initialise

  !===================================================

  !> Clean-up the GOcean environment
  subroutine gocean_finalise()
    use parallel_mod, only: parallel_finalise

    call parallel_finalise()

  end subroutine gocean_finalise

  !===================================================

  !> Stop the model run. Passes message on down to parallel_abort().
  !! @param[in] msg Message to print - reason we're stopping
  subroutine gocean_stop(msg)
    use parallel_mod, only: parallel_abort
    implicit none
    character(len=*), intent(in) :: msg

    call parallel_abort(msg)

  end subroutine gocean_stop

  !===================================================

  !> Write log entry with one integer and one real arg
  SUBROUTINE write_log_ir(fmtstr, istep, fvar, all_ranks)
    use iso_fortran_env, only : output_unit ! access computing environment
    use parallel_mod, only: on_master
    IMPLICIT none
    CHARACTER(LEN=*), INTENT(in) :: fmtstr
    INTEGER,          INTENT(in) :: istep
    REAL(go_wp),      INTENT(in) :: fvar
    LOGICAL, OPTIONAL            :: all_ranks
    lOGICAL                      :: all_ranks_value

    if(present(all_ranks)) then
        all_ranks_value = all_ranks
    else
        all_ranks_value = .false.
    endif

    if(all_ranks_value .or. on_master()) then
      WRITE(output_unit,FMT=fmtstr) istep, fvar
    endif

  END SUBROUTINE write_log_ir

  !===================================================

  !> Write log entry with one integer arg
  SUBROUTINE write_log_i(fmtstr, istep, all_ranks)
    use iso_fortran_env, only : output_unit ! access computing environment
    use parallel_mod, only: on_master
    IMPLICIT none
    CHARACTER(LEN=*), INTENT(in) :: fmtstr
    INTEGER,          INTENT(in) :: istep
    LOGICAL, OPTIONAL            :: all_ranks
    lOGICAL                      :: all_ranks_value

    if(present(all_ranks)) then
        all_ranks_value = all_ranks
    else
        all_ranks_value = .false.
    endif

    if(all_ranks_value .or. on_master()) then
      WRITE(output_unit,FMT=fmtstr) istep
    endif

  END SUBROUTINE write_log_i

  !===================================================

  !> Write log entry with one real arg
  subroutine write_log_r(fmtstr, fvar, all_ranks)
    use iso_fortran_env, only : output_unit ! access computing environment
    use parallel_mod, only: on_master
    implicit none
    character(len=*), intent(in) :: fmtstr
    real(go_wp),      intent(in) :: fvar
    LOGICAL, OPTIONAL            :: all_ranks
    lOGICAL                      :: all_ranks_value

    if(present(all_ranks)) then
        all_ranks_value = all_ranks
    else
        all_ranks_value = .false.
    endif

    if(all_ranks_value .or. on_master()) then
      write(output_unit,FMT=fmtstr) fvar
    endif

  end subroutine write_log_r

  !===================================================

  !> Write log entry with just a string
  subroutine write_log_a(fmtstr, msg, all_ranks)
    use iso_fortran_env, only : output_unit ! access computing environment
    use parallel_mod, only: on_master
    implicit none
    character(len=*), intent(in) :: fmtstr
    character(len=*), intent(in) :: msg
    LOGICAL, OPTIONAL            :: all_ranks
    lOGICAL                      :: all_ranks_value

    if(present(all_ranks)) then
        all_ranks_value = all_ranks
    else
        all_ranks_value = .false.
    endif

    if(all_ranks_value .or. on_master()) then
      write(output_unit,FMT=fmtstr) msg
    endif

  end subroutine write_log_a

  !===================================================

end module gocean_mod
