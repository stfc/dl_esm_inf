Program test
    USE field_mod
    USE kind_params_mod
    USE grid_mod
    use subdomain_mod, only : decomposition_type, decompose
      use parallel_mod, only : parallel_init

      USE io, ONLY: PsyDataPreStart, PSyDataDeclareVariable, PSyDataWriteVariable,  &
                    PSyDataPreEnd
      TYPE(r2d_field) :: cu_fld, p_fld, u_fld, cv_fld, v_fld, unew_fld, uold_fld
      TYPE(grid_type), target :: grid
      type(decomposition_type) :: decomp

      INTEGER j
      INTEGER i
      INTEGER istop, jstop

      call parallel_init()
      grid = grid_type(GO_ARAKAWA_C,(/GO_BC_PERIODIC,GO_BC_PERIODIC,GO_BC_NONE/), &
                       GO_OFFSET_SW)
      decomp = decompose(3, 3, 1, 1, 1, halo_width=1)

      ! Grid init adds a halo region of 2 automatically
      call grid_init(grid, decomp, 1.0_8, 1.0_8)

      cu_fld = r2d_field(grid, GO_T_POINTS); cu_fld%data = 0
      p_fld = r2d_field(grid, GO_T_POINTS); p_fld%data = 0
      v_fld = r2d_field(grid, GO_T_POINTS); v_fld%data = 0
      !
      ! Look-up loop bounds
      istop = cu_fld%grid%subdomain%internal%xstop
      jstop = cu_fld%grid%subdomain%internal%ystop
      !
      DO j=1,jstop
        DO i=1,istop
          cu_fld%data(i, j) = i+j
          p_fld%data(i, j) = i*i+j*j
          p_fld%data(i, j) = i-j
        END DO 
      END DO 
      !
      ! ExtractStart
      call PSyDataPreStart("tmp", "function", 5)
      CALL PSyDataDeclareVariable("cu_fld", cu_fld)
      CALL PSyDataDeclareVariable("p_fld", p_fld)
      CALL PSyDataDeclareVariable("v_fld", v_fld)
      CALL PSyDataDeclareVariable("i", i)
      CALL PSyDataDeclareVariable("j", j)
      CALL PSyDataWriteVariable("cu_fld", cu_fld)
      CALL PSyDataWriteVariable("p_fld", p_fld)
      CALL PSyDataWriteVariable("v_fld", v_fld)
      CALL PSyDataWriteVariable("i", i)
      CALL PSyDataWriteVariable("j", j)
      call PSyDataPreEnd("tmp", "function")
      !
      !$omp parallel do default(shared), private(i,j), schedule(static)
      DO j=2,jstop+1
        DO i=2,istop
        END DO 
      END DO 
      !$omp end parallel do
      !
      ! call PSyDataPostStart("tmp", "function", 3)
      ! CALL PSyDataDeclareVariable("cv_fld", cv_fld)
      ! CALL PSyDataDeclareVariable("i", i)
      ! CALL PSyDataDeclareVariable("j", j)
      ! CALL PSyDataWriteVariable("cv_fld", cv_fld)
      ! CALL PSyDataWriteVariable("i", i)
      ! CALL PSyDataWriteVariable("j", j)
      ! call PsyDataPostEnd("tmp", "function")

      ! ExtractEnd
      !
      DO j=1,jstop+1
        DO i=1,istop+1
        END DO 
      END DO 
end program test