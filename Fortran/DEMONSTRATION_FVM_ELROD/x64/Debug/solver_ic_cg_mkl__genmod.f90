        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 16:36:05 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVER_IC_CG_MKL__genmod
          INTERFACE 
            SUBROUTINE SOLVER_IC_CG_MKL(N_IN,NZE_IN,AC,SPALTENINDEX,    &
     &ANF_ZEILE,RHS,SOL,ITER_OUT,AEPS,MAXI)
              INTEGER(KIND=4) :: NZE_IN
              INTEGER(KIND=4) :: N_IN
              REAL(KIND=8) :: AC(NZE_IN)
              INTEGER(KIND=4) :: SPALTENINDEX(NZE_IN)
              INTEGER(KIND=4) :: ANF_ZEILE(N_IN+1)
              REAL(KIND=8) :: RHS(N_IN)
              REAL(KIND=8) :: SOL(N_IN)
              INTEGER(KIND=4) :: ITER_OUT
              REAL(KIND=8) :: AEPS
              INTEGER(KIND=4) :: MAXI
            END SUBROUTINE SOLVER_IC_CG_MKL
          END INTERFACE 
        END MODULE SOLVER_IC_CG_MKL__genmod
