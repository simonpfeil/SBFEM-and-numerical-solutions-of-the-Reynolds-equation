SUBROUTINE SOLVER_IC_CG_MKL(n_in, nze_in, AC, SPALTENINDEX, ANF_ZEILE, RHS, SOL, iter_out, aeps, maxi)
    
  ! AUTHOR, AFFILIATION: Steffen Nitzschke, OvGU Magdeburg
  ! EDITS:
  ! 07.08.24 (Simon Pfeil, OvGU Magdeburg): n and nze are now passed to the program (via n_in and nze_in) rather than determined by SIZE()
  ! 07.08.24 (Simon Pfeil, OvGU Magdeburg): the arrays are now passed to the program with known shape instead of assumed shape
  ! 07.08.24 (Simon Pfeil, OvGU Magdeburg): the arguments that used to be optional are now assumed to always be present
  ! 07.08.24 (Simon Pfeil, OvGU Magdeburg): new definition of error in convergence criterion (for consistency with BiCGStab)
    
  USE MODULE_JPICC, ONLY: JPICC
  IMPLICIT NONE
  
  ! global variables
  INTEGER                                   :: n_in
  INTEGER                                   :: nze_in
  REAL(KIND=8),DIMENSION(nze_in)            :: AC 
  INTEGER,     DIMENSION(nze_in)            :: SPALTENINDEX
  INTEGER,     DIMENSION(n_in+1)            :: ANF_ZEILE
  REAL(KIND=8),DIMENSION(n_in)              :: RHS
  REAL(KIND=8),DIMENSION(n_in)              :: SOL
  INTEGER                                   :: iter_out
  REAL(KIND=8)                              :: aeps
  INTEGER                                   :: maxi
    
  ! local variables
  INTEGER                                   :: i
  INTEGER                                   :: j
  INTEGER                                   :: n
  INTEGER                                   :: nze
  INTEGER                                   :: iter
  INTEGER                                   :: retval
  INTEGER,DIMENSION(:),ALLOCATABLE          :: IA_IC
  INTEGER,DIMENSION(:),ALLOCATABLE          :: IE_IC
  INTEGER,DIMENSION(:),ALLOCATABLE          :: SI_IC
  INTEGER,DIMENSION(:),ALLOCATABLE          :: ITCOL
  INTEGER,DIMENSION(:),ALLOCATABLE          :: IFIRST
  INTEGER,DIMENSION(:),ALLOCATABLE          :: LIST
  
  REAL(KIND=8)                              :: one
  REAL(KIND=8)                              :: zero
  REAL(KIND=8)                              :: ka,kb,rtr0,rtr1,btb,btab,fehler
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE     :: TEMP
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE     :: R,S,V,H
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE     :: A_IC
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE     :: HD_IC
  
  n = n_in
  nze = nze_in
  
  ! Inits
  one    = 1.0d0
  zero   = 0d0
  iter   = 0
  fehler = 1.D30  
  
  ALLOCATE(A_IC(nze))
  ALLOCATE(SI_IC(nze))
  ALLOCATE(IA_IC(n))
  ALLOCATE(IE_IC(n))
  ALLOCATE(HD_IC(n))
  ALLOCATE(ITCOL(n))
  ALLOCATE(IFIRST(n))
  ALLOCATE(LIST(n))
  ALLOCATE(TEMP(n))
  ALLOCATE(R(n))
  ALLOCATE(S(n))
  ALLOCATE(V(n))
  ALLOCATE(H(n))
  TEMP = zero
  
  ! get preconditioner: 
  ! 
  ! JPICC requires lower triangular matrix with column-wise storage of the off-diagonal entries and separate storage of the main diagonal,
  ! but the matrix is given in upper triangular form with row-wise storage of all entries. Below, the format is converted accordingly.
  !
  ! isolate main diagonal
  A_IC = AC
  HD_IC= AC(ANF_ZEILE(1:n))
  ! prevent access to these elements in compact storate within JPICC
  IA_IC = ANF_ZEILE(1:n)+1
  IE_IC = ANF_ZEILE(2:n+1)-1
  IA_IC(n)=IA_IC(n)-1
  IE_IC(n)=IE_IC(n)-1
  SI_IC = SPALTENINDEX
  ! perform Incomplete Cholesky
  retval= JPICC(N,HD_IC,A_IC,IA_IC,IE_IC,SI_IC,TEMP,ITCOL,IFIRST,LIST)
  IF (retval.EQ.0) THEN
    ! re-incorporate main diagonal into compact storage format after successful IC decomposition
    A_IC(ANF_ZEILE(1:n))= HD_IC
  ELSE
    WRITE(*,*) "Error in preconditioning in SOLVER_IC_CG_MKL.f90"
    STOP
  END IF
  DEALLOCATE(ITCOL,IFIRST,LIST,IA_IC,IE_IC)
  
  ! determine residuum of initial solution res=rhs-A*x0
  CALL MKL_DCSRSYMV('U', N, AC, ANF_ZEILE, SPALTENINDEX, SOL, TEMP)
  R = RHS - TEMP
  
  S    = PREC_SSOR(R)
  V    = S
  rtr0 = DOT_PRODUCT(S,R)
  rtr1 = rtr0  
  
  ! norm of the right-hand side
  btb = DOT_PRODUCT(PREC_SSOR(RHS), RHS)  
  
  DO WHILE (fehler > aeps .AND. iter <= maxi)
    iter = iter + 1
    ! A*s=v
    ! ---------------------
    ! Matrix multiplied by vector: call mkl_dcsrsymv(uplo, n, a, ia, ja, x, y)
    ! ---------------------    
    !  Computes matrix-vector product y= A*x of a sparse symmetrical matrix stored in the 
    !  CSR format (3-array variation) with one-based indexing.
    ! 'U' - In: uplo-Parameter: 'U' means UpperTriangle, 'L' means LowerTriangle
    !  N  - In: number of rows in matrix
    !  A  - In: matrix entries (nonzeros)
    !  IA - In: index of the respective first entry of the row
    !           (dim. is N+1, IA(N+1)=ANZ_NNE+1 ) --> pointer to main diagonal elements
    !  JA - In: columns of the matrix entries
    CALL MKL_DCSRSYMV('U', N, AC, ANF_ZEILE, SPALTENINDEX, V, H)
    btab= DOT_PRODUCT(V,H)
    ka  = rtr1 / btab 
    SOL = SOL + ka * V
    
    R   = R   - ka * H
    S   = PREC_SSOR(R)
    rtr1= DOT_PRODUCT(R,S)
    
    kb  = rtr1 / rtr0
    rtr0= rtr1
    V   = S + kb * V
    !fehler = SQRT(abs(rtr1/btb)) 
    fehler = SQRT(SUM(R**2)/SUM(RHS**2))
    !write(*,*) iter, fehler
  ENDDO
  
  DEALLOCATE(TEMP)
  DEALLOCATE(A_IC)
  DEALLOCATE(R)
  DEALLOCATE(S)
  DEALLOCATE(V)
  DEALLOCATE(H)
  
  iter_out = iter
  
  CONTAINS
  
  ! --------------------
  ! Preconditioning
  ! --------------------  
  FUNCTION PREC_SSOR(Y) RESULT(X)
    IMPLICIT NONE
    
    REAL(KIND=8),DIMENSION(:)      :: Y
    REAL(KIND=8),DIMENSION(SIZE(Y)):: X
    
    CHARACTER(LEN=1), DIMENSION(6) :: MATDES
    
    ! solution of M*rho=r with M=L*L^T
    ! --> solution can be split in two steps:
    MATDES(1)='T'   ! Triangular
    MATDES(2)='U'   ! Upper Triangular
    !MATDES(3)='U'   ! (U)nit/(N)on-Unit Main Diagonal
    MATDES(3)='N'   ! (U)nit/(N)on-Unit Main Diagonal
    MATDES(4)='F'   ! Indizierung der NNe faengt mit "1" an (statt mit "0")
    ! step 1: L^T*y=r
    CALL MKL_DCSRSV('T', n, ONE, MATDES,     &
                    A_IC, SPALTENINDEX, ANF_ZEILE, ANF_ZEILE(2), Y, TEMP)
    ! step 2: L*rho=y
    CALL MKL_DCSRSV('N', n, ONE, MATDES,     &
                    A_IC, SPALTENINDEX, ANF_ZEILE, ANF_ZEILE(2), TEMP, X)    
    
  END FUNCTION
  
END SUBROUTINE SOLVER_IC_CG_MKL