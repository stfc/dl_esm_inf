
module IO

    ! Contains functionality to write and read variables

    interface PSyDataDeclareVariable
        module procedure DeclareScalarInteger
        !module procedure DeclareScalarReal
        !module procedure DeclareScalarDoublePrecision
        module procedure DeclareRealField
    end interface PSyDataDeclareVariable
    interface PSyDataWriteVariable
        module procedure WriteScalarInteger
        !module procedure WriteScalarReal
        !module procedure WriteScalarDoublePrecision
        module procedure WriteRealField
    end interface PSyDataWriteVariable

    integer, private :: ncid
    integer, dimension(:), allocatable :: var_id
    integer :: next_var_index
Contains

    ! -------------------------------------------------------------------------
    function CheckError(retval) 
        use netcdf
        implicit none
        integer :: CheckError
        integer, intent(in) :: retval
        if (retval /= nf90_noerr) then
            print *,"NetCDF Error:"
            print *,trim(nf90_strerror(retval))
        endif
        CheckError = retval
    end function CheckError

    ! -------------------------------------------------------------------------
    subroutine PsyDataPreStart(module_name, kernel_name, num_vars)
        use netcdf, only : nf90_create, NF90_CLOBBER
        implicit none
        character(*), intent(in) :: module_name, kernel_name
        integer, intent(in)      :: num_vars
        integer :: retval

        retval = CheckError(nf90_create(module_name//kernel_name//".nc", &
                                        NF90_CLOBBER, ncid))
        allocate(var_id(num_vars))
        next_var_index = 1
    end subroutine PSyDataPreStart

    ! -------------------------------------------------------------------------
    subroutine PSyDataPreEnd(module_name, kernel_name)
        use netcdf, only : nf90_close
        implicit none
        character(*), intent(in) :: module_name, kernel_name
        integer :: retval
        retval = CheckError(nf90_close(ncid))
    end subroutine PSyDataPreEnd

    ! -------------------------------------------------------------------------
    subroutine CheckLastDeclaration()
        use netcdf
        implicit none
        integer :: retval
        ! Test if this was the last declaration
        if (next_var_index > size(var_id,1)) then
            next_var_index = 1
            retval = CheckError(nf90_enddef(ncid))
        endif   
    end subroutine CheckLastDeclaration

    ! -------------------------------------------------------------------------
    subroutine DeclareScalarInteger(name, value)
        use netcdf
        implicit none
        character(*), intent(in) :: name
        integer, intent(in) :: value
        integer :: retval
        retval = CheckError(nf90_def_var(ncid, name, Nf90_INT,     &
                                         var_id(next_var_index)))
        next_var_index = next_var_index + 1
        call CheckLastDeclaration()
    end subroutine DeclareScalarInteger

    ! -------------------------------------------------------------------------
    subroutine WriteScalarInteger(name, value)
        use netcdf
        implicit none
        character(*), intent(in) :: name
        integer, intent(in) :: value
        integer :: retval
        retval = CheckError(nf90_put_var(ncid, var_id(next_var_index), value))
        next_var_index = next_var_index + 1
    end subroutine WriteScalarInteger

    ! -------------------------------------------------------------------------
    subroutine WriteScalarReal(name, value)
        implicit none
        character(*), intent(in) :: name
        real, intent(in) :: value
    end subroutine WriteScalarReal

    ! -------------------------------------------------------------------------
    subroutine WriteScalarDoublePrecision(name, value)
        implicit none
        character(*), intent(in) :: name
        double precision, intent(in) :: value
    end subroutine WriteScalarDoublePrecision
    
    ! -------------------------------------------------------------------------
    subroutine DeclareRealField(name, value)
        use netcdf
        use field_mod, only : r2d_field
        implicit none
        character(*), intent(in) :: name
        type(r2d_field), intent(in) :: value
        integer :: x_dimid, y_dimid, retval
        integer, dimension(2) :: dimids

        retval = CheckError(nf90_def_dim(ncid, name//"dim1",  &
                                         value%grid%nx, x_dimid))
        retval = CheckError(nf90_def_dim(ncid, name//"dim2",  &
                                         value%grid%ny, y_dimid))
         dimids =  (/ x_dimid, y_dimid /)

        retval = CheckError(nf90_def_var(ncid, name, Nf90_REAL8,     &
                                         dimids, var_id(next_var_index)))

        next_var_index = next_var_index + 1
        call CheckLastDeclaration()
    end subroutine DeclareRealField

    ! -------------------------------------------------------------------------
    subroutine WriteRealField(name, value)
        use netcdf
        use field_mod, only : r2d_field
        implicit none
        character(*), intent(in) :: name
        type(r2d_field), intent(in) :: value
        integer :: retval
        print *,"XXX", shape(value%data)
        retval = CheckError(nf90_put_var(ncid, var_id(next_var_index), value%data(:,:)))

        next_var_index = next_var_index + 1
        print *,"Writing field ", name
    end subroutine WriteRealField

    ! -------------------------------------------------------------------------

end module IO