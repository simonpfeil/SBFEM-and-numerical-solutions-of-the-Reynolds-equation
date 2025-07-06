        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 30 18:40:06 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FVM_ASSEMBLY_NOZEROS__genmod
          INTERFACE 
            SUBROUTINE FVM_ASSEMBLY_NOZEROS(N_X,DOF_Y,DOF,N_NZ,N,       &
     &ROWINDEX_VEC,VALUES_VEC,COLUMNS_VEC,MD_VEC,G_VEC,NNN_VEC,NNE_VEC, &
     &NNS_VEC,NNW_VEC,COEF0_VEC,COEFN_VEC,COEFE_VEC,COEFS_VEC,COEFW_VEC,&
     &SGN_U)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: N_NZ
              INTEGER(KIND=4), INTENT(IN) :: DOF
              INTEGER(KIND=4), INTENT(IN) :: N_X
              INTEGER(KIND=4), INTENT(IN) :: DOF_Y
              INTEGER(KIND=4), INTENT(OUT) :: ROWINDEX_VEC(DOF+1)
              REAL(KIND=8), INTENT(OUT) :: VALUES_VEC(N_NZ)
              INTEGER(KIND=4), INTENT(OUT) :: COLUMNS_VEC(N_NZ)
              INTEGER(KIND=4), INTENT(OUT) :: MD_VEC(DOF)
              INTEGER(KIND=4), INTENT(IN) :: G_VEC(N)
              INTEGER(KIND=4), INTENT(IN) :: NNN_VEC(DOF)
              INTEGER(KIND=4), INTENT(IN) :: NNE_VEC(DOF)
              INTEGER(KIND=4), INTENT(IN) :: NNS_VEC(DOF)
              INTEGER(KIND=4), INTENT(IN) :: NNW_VEC(DOF)
              REAL(KIND=8), INTENT(IN) :: COEF0_VEC(DOF)
              REAL(KIND=8), INTENT(IN) :: COEFN_VEC(DOF)
              REAL(KIND=8), INTENT(IN) :: COEFE_VEC(DOF)
              REAL(KIND=8), INTENT(IN) :: COEFS_VEC(DOF)
              REAL(KIND=8), INTENT(IN) :: COEFW_VEC(DOF)
              REAL(KIND=8), INTENT(IN) :: SGN_U
            END SUBROUTINE FVM_ASSEMBLY_NOZEROS
          END INTERFACE 
        END MODULE FVM_ASSEMBLY_NOZEROS__genmod
