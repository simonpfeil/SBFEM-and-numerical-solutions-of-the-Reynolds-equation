        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul  3 16:35:38 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FVM_ELROD__genmod
          INTERFACE 
            SUBROUTINE FVM_ELROD(D_B,L_B,C,GROOVES,X_OS,L_X_OS,L_Y_OS,  &
     &P_OS,AC_VEC,T,ANGLE_SHELL,OMEGA_SHELL,OMEGA_SHAFT,DIS_H_SHELL,    &
     &DIS_V_SHELL,VEL_H_SHELL,VEL_V_SHELL,DIS_H_SHAFT,DIS_V_SHAFT,      &
     &VEL_H_SHAFT,VEL_V_SHAFT,TILT_H_SHELL,TILT_V_SHELL,TILT_DOT_H_SHELL&
     &,TILT_DOT_V_SHELL,TILT_H_SHAFT,TILT_V_SHAFT,TILT_DOT_H_SHAFT,     &
     &TILT_DOT_V_SHAFT,N_X,ITER_MAX,MU_VEC,GUEMBEL,N_Y,QUASISTATIC,SYMBC&
     &,F_H,F_V,M_H,M_V,M_FR,V_OIL,V_DOT_BB,PI_MAT,P_REF,CONVERGENT,ITER,&
     &ITER_SOL,PTS_VEC,TOL,ITER_MAX_SOLVER,PM)
              INTEGER(KIND=4), INTENT(IN) :: N_Y
              INTEGER(KIND=4), INTENT(IN) :: N_X
              INTEGER(KIND=4), INTENT(IN) :: GROOVES
              REAL(KIND=8), INTENT(IN) :: D_B
              REAL(KIND=8), INTENT(IN) :: L_B
              REAL(KIND=8), INTENT(IN) :: C
              REAL(KIND=8), INTENT(IN) :: X_OS
              REAL(KIND=8), INTENT(IN) :: L_X_OS
              REAL(KIND=8), INTENT(IN) :: L_Y_OS
              REAL(KIND=8), INTENT(IN) :: P_OS
              REAL(KIND=8), INTENT(IN) :: AC_VEC(N_X*N_Y)
              REAL(KIND=8), INTENT(IN) :: T
              REAL(KIND=8), INTENT(IN) :: ANGLE_SHELL
              REAL(KIND=8), INTENT(IN) :: OMEGA_SHELL
              REAL(KIND=8), INTENT(IN) :: OMEGA_SHAFT
              REAL(KIND=8), INTENT(IN) :: DIS_H_SHELL
              REAL(KIND=8), INTENT(IN) :: DIS_V_SHELL
              REAL(KIND=8), INTENT(IN) :: VEL_H_SHELL
              REAL(KIND=8), INTENT(IN) :: VEL_V_SHELL
              REAL(KIND=8), INTENT(IN) :: DIS_H_SHAFT
              REAL(KIND=8), INTENT(IN) :: DIS_V_SHAFT
              REAL(KIND=8), INTENT(IN) :: VEL_H_SHAFT
              REAL(KIND=8), INTENT(IN) :: VEL_V_SHAFT
              REAL(KIND=8), INTENT(IN) :: TILT_H_SHELL
              REAL(KIND=8), INTENT(IN) :: TILT_V_SHELL
              REAL(KIND=8), INTENT(IN) :: TILT_DOT_H_SHELL
              REAL(KIND=8), INTENT(IN) :: TILT_DOT_V_SHELL
              REAL(KIND=8), INTENT(IN) :: TILT_H_SHAFT
              REAL(KIND=8), INTENT(IN) :: TILT_V_SHAFT
              REAL(KIND=8), INTENT(IN) :: TILT_DOT_H_SHAFT
              REAL(KIND=8), INTENT(IN) :: TILT_DOT_V_SHAFT
              INTEGER(KIND=4), INTENT(IN) :: ITER_MAX
              REAL(KIND=8), INTENT(IN) :: MU_VEC(N_X*N_Y)
              INTEGER(KIND=4), INTENT(IN) :: GUEMBEL
              INTEGER(KIND=4), INTENT(IN) :: QUASISTATIC
              INTEGER(KIND=4), INTENT(IN) :: SYMBC
              REAL(KIND=8), INTENT(OUT) :: F_H
              REAL(KIND=8), INTENT(OUT) :: F_V
              REAL(KIND=8), INTENT(OUT) :: M_H
              REAL(KIND=8), INTENT(OUT) :: M_V
              REAL(KIND=8), INTENT(OUT) :: M_FR
              REAL(KIND=8), INTENT(OUT) :: V_OIL
              REAL(KIND=8), INTENT(OUT) :: V_DOT_BB
              REAL(KIND=8), INTENT(OUT) :: PI_MAT(N_X,N_Y)
              REAL(KIND=8), INTENT(OUT) :: P_REF
              INTEGER(KIND=4), INTENT(OUT) :: CONVERGENT
              INTEGER(KIND=4), INTENT(OUT) :: ITER
              INTEGER(KIND=4), INTENT(OUT) :: ITER_SOL
              REAL(KIND=8), INTENT(INOUT) :: PTS_VEC(N_X*N_Y+1)
              REAL(KIND=8), INTENT(IN) :: TOL
              INTEGER(KIND=4), INTENT(IN) :: ITER_MAX_SOLVER
              REAL(KIND=8), INTENT(IN) :: PM
            END SUBROUTINE FVM_ELROD
          END INTERFACE 
        END MODULE FVM_ELROD__genmod
