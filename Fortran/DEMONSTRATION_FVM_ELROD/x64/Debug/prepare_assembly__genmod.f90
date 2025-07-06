        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 19 11:01:18 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PREPARE_ASSEMBLY__genmod
          INTERFACE 
            SUBROUTINE PREPARE_ASSEMBLY(N_X,DOF_Y,DOF,N_NZ,ROWINDEX_VEC,&
     &COEFINDEX_VEC,COLUMNS_VEC,MD_VEC)
              INTEGER(KIND=4), INTENT(IN) :: N_NZ
              INTEGER(KIND=4), INTENT(IN) :: DOF
              INTEGER(KIND=4), INTENT(IN) :: N_X
              INTEGER(KIND=4), INTENT(IN) :: DOF_Y
              INTEGER(KIND=4), INTENT(IN) :: ROWINDEX_VEC(DOF+1)
              INTEGER(KIND=4), INTENT(OUT) :: COEFINDEX_VEC(N_NZ)
              INTEGER(KIND=4), INTENT(OUT) :: COLUMNS_VEC(N_NZ)
              INTEGER(KIND=4), INTENT(OUT) :: MD_VEC(DOF)
            END SUBROUTINE PREPARE_ASSEMBLY
          END INTERFACE 
        END MODULE PREPARE_ASSEMBLY__genmod
