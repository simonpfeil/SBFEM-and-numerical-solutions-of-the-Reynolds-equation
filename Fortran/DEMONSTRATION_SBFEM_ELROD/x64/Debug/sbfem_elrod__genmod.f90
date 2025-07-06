        !COMPILER-GENERATED INTERFACE MODULE: Thu Aug  8 16:03:33 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SBFEM_ELROD__genmod
          INTERFACE 
            SUBROUTINE SBFEM_ELROD(GROOVES,N_X,ITER_MAX,QUASISTATIC,    &
     &GUEMBEL,N_Y,D_B,L_B,C,X_OS,L_X_OS,P_OS,T,ANGLE_SHELL,OMEGA_SHELL, &
     &DIS_H_SHELL,DIS_V_SHELL,VEL_H_SHELL,VEL_V_SHELL,OMEGA_SHAFT,      &
     &DIS_H_SHAFT,DIS_V_SHAFT,VEL_H_SHAFT,VEL_V_SHAFT,AC_VEC,MU_VEC,TAY,&
     &N_TAY,N_LD,RED,EPS_CONSTR,EPS_MAX,LD_VEC,VD_MAT,CONVERGENT,ITER,  &
     &M_FR,V_OIL,V_DOT_BB,P_REF,F_H,F_V,G_VEC,THETA_VEC,PI_MAT,PTS_VEC)
              INTEGER(KIND=4), INTENT(IN) :: N_LD
              INTEGER(KIND=4), INTENT(IN) :: N_Y
              INTEGER(KIND=4), INTENT(IN) :: N_X
              INTEGER(KIND=4), INTENT(IN) :: GROOVES
              INTEGER(KIND=4), INTENT(IN) :: ITER_MAX
              INTEGER(KIND=4), INTENT(IN) :: QUASISTATIC
              INTEGER(KIND=4), INTENT(IN) :: GUEMBEL
              REAL(KIND=8), INTENT(IN) :: D_B
              REAL(KIND=8), INTENT(IN) :: L_B
              REAL(KIND=8), INTENT(IN) :: C
              REAL(KIND=8), INTENT(IN) :: X_OS
              REAL(KIND=8), INTENT(IN) :: L_X_OS
              REAL(KIND=8), INTENT(IN) :: P_OS
              REAL(KIND=8), INTENT(IN) :: T
              REAL(KIND=8), INTENT(IN) :: ANGLE_SHELL
              REAL(KIND=8), INTENT(IN) :: OMEGA_SHELL
              REAL(KIND=8), INTENT(IN) :: DIS_H_SHELL
              REAL(KIND=8), INTENT(IN) :: DIS_V_SHELL
              REAL(KIND=8), INTENT(IN) :: VEL_H_SHELL
              REAL(KIND=8), INTENT(IN) :: VEL_V_SHELL
              REAL(KIND=8), INTENT(IN) :: OMEGA_SHAFT
              REAL(KIND=8), INTENT(IN) :: DIS_H_SHAFT
              REAL(KIND=8), INTENT(IN) :: DIS_V_SHAFT
              REAL(KIND=8), INTENT(IN) :: VEL_H_SHAFT
              REAL(KIND=8), INTENT(IN) :: VEL_V_SHAFT
              REAL(KIND=8), INTENT(IN) :: AC_VEC(N_X)
              REAL(KIND=8), INTENT(IN) :: MU_VEC(N_X)
              INTEGER(KIND=4), INTENT(IN) :: TAY
              INTEGER(KIND=4), INTENT(IN) :: N_TAY
              INTEGER(KIND=4), INTENT(IN) :: RED
              REAL(KIND=8), INTENT(IN) :: EPS_CONSTR
              REAL(KIND=8), INTENT(IN) :: EPS_MAX
              REAL(KIND=8), INTENT(IN) :: LD_VEC(N_LD)
              REAL(KIND=8), INTENT(IN) :: VD_MAT(N_X,N_LD)
              INTEGER(KIND=4), INTENT(OUT) :: CONVERGENT
              INTEGER(KIND=4), INTENT(OUT) :: ITER
              REAL(KIND=8), INTENT(OUT) :: M_FR
              REAL(KIND=8), INTENT(OUT) :: V_OIL
              REAL(KIND=8), INTENT(OUT) :: V_DOT_BB
              REAL(KIND=8), INTENT(OUT) :: P_REF
              REAL(KIND=8), INTENT(OUT) :: F_H
              REAL(KIND=8), INTENT(OUT) :: F_V
              INTEGER(KIND=4), INTENT(OUT) :: G_VEC(N_X)
              REAL(KIND=8), INTENT(OUT) :: THETA_VEC(N_X)
              REAL(KIND=8), INTENT(OUT) :: PI_MAT(N_X,N_Y)
              REAL(KIND=8), INTENT(INOUT) :: PTS_VEC(N_X+1)
            END SUBROUTINE SBFEM_ELROD
          END INTERFACE 
        END MODULE SBFEM_ELROD__genmod
