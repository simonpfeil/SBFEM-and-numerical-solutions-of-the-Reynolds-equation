        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 19 11:01:18 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ASSEMBLY_STRCTSYM__genmod
          INTERFACE 
            SUBROUTINE ASSEMBLY_STRCTSYM(N_X,DOF_Y,DOF,N_NZ,N,          &
     &ROWINDEX_VEC,VALUES_VEC,COLUMNS_VEC,MD_VEC,G_VEC,NNN_VEC,NNE_VEC, &
     &NNS_VEC,NNW_VEC,COEF0_VEC,COEFN_VEC,COEFE_VEC,COEFS_VEC,COEFW_VEC,&
     &NN0_VEC)
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
              INTEGER(KIND=4), INTENT(IN) :: NN0_VEC(DOF)
            END SUBROUTINE ASSEMBLY_STRCTSYM
          END INTERFACE 
        END MODULE ASSEMBLY_STRCTSYM__genmod
