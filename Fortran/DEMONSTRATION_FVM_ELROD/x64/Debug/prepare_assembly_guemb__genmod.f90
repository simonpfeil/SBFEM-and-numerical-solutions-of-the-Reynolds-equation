        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 16:35:38 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PREPARE_ASSEMBLY_GUEMB__genmod
          INTERFACE 
            SUBROUTINE PREPARE_ASSEMBLY_GUEMB(N_X,DOF_Y,DOF,N_NZ,       &
     &ROWINDEX_VEC,COEFINDEX_VEC,COLUMNS_VEC)
              INTEGER(KIND=4), INTENT(IN) :: N_NZ
              INTEGER(KIND=4), INTENT(IN) :: DOF
              INTEGER(KIND=4), INTENT(IN) :: N_X
              INTEGER(KIND=4), INTENT(IN) :: DOF_Y
              INTEGER(KIND=4), INTENT(OUT) :: ROWINDEX_VEC(DOF+1)
              INTEGER(KIND=4), INTENT(OUT) :: COEFINDEX_VEC(N_NZ)
              INTEGER(KIND=4), INTENT(OUT) :: COLUMNS_VEC(N_NZ)
            END SUBROUTINE PREPARE_ASSEMBLY_GUEMB
          END INTERFACE 
        END MODULE PREPARE_ASSEMBLY_GUEMB__genmod
