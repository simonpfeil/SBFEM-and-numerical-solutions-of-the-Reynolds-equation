! FVM_ELROD.f90 
!
! FUNCTIONS:
! FVM_ELROD - program for solving the Reynolds equation
! PREPARE_ASSEMBLY_ELROD - used by FVM_ELROD to prepare the sparse matrix assembly under Elrod conditions
! PREPARE_ASSEMBLY_GUEMB - used by FVM_ELROD to prepare the sparse matrix assembly under Guembel conditions
!
! AUTHOR, AFFILIATION, DATE: Simon Pfeil, OvGU Magdeburg (Germany), 30.06.2025
!
! **********************************************************************************************************
!
! PROGRAM: FVM_ELROD
!
! PURPOSE: This program solves the Reynolds equation numerically under consideration of transient, mass-
! conserving cavitation. For clarification of the input and output variables and how to use this program,
! check the attached file DEMONSTRATION_FVM_ELROD.f90.
!
! **********************************************************************************************************
!
! Check the attached script DEMONSTRATION_FVM_ELROD.f90 for clarifications.
!
! **********************************************************************************************************



SUBROUTINE FVM_ELROD(d_b, l_b, c, grooves, X_os, L_X_os, l_y_os, p_os, ac_vec, t, angle_shell, &
  omega_shell, omega_shaft, dis_h_shell, dis_v_shell, vel_h_shell, vel_v_shell, dis_h_shaft, &
  dis_v_shaft, vel_h_shaft, vel_v_shaft, tilt_h_shell, tilt_v_shell, tilt_dot_h_shell, &
  tilt_dot_v_shell, tilt_h_shaft, tilt_v_shaft, tilt_dot_h_shaft, tilt_dot_v_shaft, n_x, &
  iter_max, mu_vec, guembel, n_y, quasistatic, symBC, F_h, F_v, M_h, M_v, M_fr, V_oil, &
  V_dot_bb, Pi_mat, p_ref, convergent, iter, iter_sol, pts_vec, tol, iter_max_solver, pm)
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Modules and variable types
  ! --------------------------------------------------------------------------------------------------------
  
  ! modules
  USE MODULE_SPARSKIT2_BCGSTAB
  
  ! no implicit variables
  IMPLICIT NONE
  
  ! input variables
  REAL(KIND=8),INTENT(IN)                               :: d_b, l_b
  INTEGER,INTENT(IN)                                    :: n_x, n_y, grooves, iter_max_solver, quasistatic
  INTEGER,INTENT(IN)                                    :: iter_max, symBC, guembel
  REAL(KIND=8),INTENT(IN)                               :: c, X_os, L_X_os, l_y_os, p_os, t, pm, tol
  REAL(KIND=8),INTENT(IN)                               :: angle_shell, omega_shell, omega_shaft
  REAL(KIND=8),INTENT(IN)                               :: dis_h_shell, dis_v_shell, vel_h_shell
  REAL(KIND=8),INTENT(IN)                               :: vel_v_shell, dis_h_shaft, dis_v_shaft
  REAL(KIND=8),INTENT(IN)                               :: vel_h_shaft, vel_v_shaft
  REAL(KIND=8),INTENT(IN)                               :: tilt_h_shell, tilt_v_shell, tilt_dot_h_shell
  REAL(KIND=8),INTENT(IN)                               :: tilt_dot_v_shell, tilt_h_shaft, tilt_v_shaft
  REAL(KIND=8),INTENT(IN)                               :: tilt_dot_h_shaft, tilt_dot_v_shaft
  REAL(KIND=8),DIMENSION(n_x*n_y),INTENT(IN)            :: mu_vec, ac_vec
  
  ! output variables
  INTEGER,INTENT(OUT)                                   :: convergent, iter, iter_sol
  REAL(KIND=8),INTENT(OUT)                              :: F_h, F_v, M_h, M_v, M_fr, V_oil, V_dot_bb, p_ref
  REAL(KIND=8),DIMENSION(n_x,n_y),INTENT(OUT)           :: Pi_mat
  
  ! input/output variables
  REAL(KIND=8),DIMENSION(n_x*n_y+1),INTENT(INOUT)       :: pts_vec
  
  ! solver-specific local variables
  INTEGER                                               :: error, dof, lfil, nnz_ilu, lws, liws, iter_out
  INTEGER,DIMENSION(16)                                 :: ipar_vec
  INTEGER,DIMENSION(:),ALLOCATABLE                      :: rowIndex_vec, columns_vec
  INTEGER,DIMENSION(:),ALLOCATABLE                      :: j_lu_vec, j_u_vec, iws_vec
  REAL(KIND=8)                                          :: tol_fil
  REAL(KIND=8),DIMENSION(16)                            :: fpar_vec
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE                 :: values_vec, R_vec, Pi_dof_vec
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE                 :: initial_guess_vec, a_lu_vec, ws_vec
  
  ! other local variables
  INTEGER                                               :: n, dof_y, j_start, i, j, k, l, n_os_x, n_i
  INTEGER                                               :: start_os, longgroove, end_os, n_os_y, n_nz
  INTEGER                                               :: restore13, g_os, s1, s2
  INTEGER,DIMENSION(grooves)                            :: start_os_vec, end_os_vec
  INTEGER,DIMENSION(n_x*n_y)                            :: g_pts_vec, g_vec, g_old_vec, g_old_old_vec
  INTEGER,DIMENSION(:),ALLOCATABLE                      :: nn0_vec, nnN_vec, nnE_vec, nnS_vec, nnW_vec
  INTEGER,DIMENSION(:),ALLOCATABLE                      :: os_x_vec, os_y_vec, os_vec, coefIndex_vec, md_vec
  REAL(KIND=8)                                          :: u, dis_h, dis_v, vel_h, vel_v, q, X_att, q_dot
  REAL(KIND=8)                                          :: X_att_dot, tilt_h, tilt_v, tilt_dot_h, tilt_dot_v
  REAL(KIND=8)                                          :: tilt, tilt_dot, X_tilt, X_tilt_dot, omega
  REAL(KIND=8)                                          :: r_b, sgn_u, mu_ref, epsil, epsil_dot, pi, L_X
  REAL(KIND=8)                                          :: L_Y, a, Delta_T, t_pts, pf, A_i, Pi_os
  REAL(KIND=8)                                          :: F_1, F_2, M_1, M_2
  REAL(KIND=8)                                          :: zero, q_squared, tilt_squared
  REAL(KIND=8),DIMENSION(grooves)                       :: X_os_vec
  REAL(KIND=8),DIMENSION(n_x*n_y)                       :: X_vec, Y_vec, Pi_pts_vec, cos_vec, sin_vec, H_vec
  REAL(KIND=8),DIMENSION(n_x*n_y)                       :: dHdT_vec, H3_vec, H3_over_mu_vec, Pi_vec
  REAL(KIND=8),DIMENSION(n_x*n_y)                       :: theta_vec, p_vec, A_vec
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE                 :: aN_vec, aE_vec, aS_vec, aW_vec, coef_vec
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Switch to reference frame of shell
  ! --------------------------------------------------------------------------------------------------------
  
  dis_h = dis_h_shaft-dis_h_shell                                                                           ! horizontal displacement of the shaft relative to the shell, still in the reference frame of the inertial system
  dis_v = dis_v_shaft-dis_v_shell                                                                           ! vertical displacement of the shaft relative to the shell, still in the reference frame of the inertial system
  q_squared = dis_h**2+dis_v**2
  zero = EPSILON(q_squared)
  IF ( q_squared .LT. 5.0d0*zero ) THEN
    q_squared = 5.0d0*zero
    q = SQRT(q_squared)
    dis_v = -q
    dis_h = 0.0d0
  ELSE
    q = SQRT(q_squared)                                                                                     ! absolute eccentricity
  END IF  
  X_att = ATAN2(dis_v,dis_h)                                                                                ! attitude angle in the reference frame of the inertial system
  X_att = X_att-angle_shell                                                                                 ! attitude angle in the shell-fixed reference frame in which the Reynolds equation is solved
  tilt_h = tilt_h_shaft-tilt_h_shell                                                                        ! tilting angle of the shaft relative to the shell around the horizontal axis, still in the reference frame of the inertial system
  tilt_v = tilt_v_shaft-tilt_v_shell                                                                        ! tilting angle of the shaft relative to the shell around the vertical axis, still in the reference frame of the inertial system
  tilt_squared = tilt_h**2+tilt_v**2
  IF ( tilt_squared .LT. 5.0d0*zero ) THEN
    tilt_squared = 5.0d0*zero
    tilt = SQRT(tilt_squared)
    tilt_v = -tilt
    tilt_h = 0.0d0
  ELSE
    tilt = SQRT(tilt_squared)                                                                               ! tilting angle
  END IF  
  X_tilt = ATAN2(tilt_v,tilt_h)                                                                             ! angle describing the resulting tilting axis in the reference frame of the inertial system   
  X_tilt = X_tilt-angle_shell                                                                               ! angle describing the tilting axis in the shell-fixed reference frame in which the Reynolds equation is solved
  IF ( quasistatic .EQ. 0 ) THEN
    omega = omega_shaft-omega_shell                                                                         ! rotational velocity of the shaft relative to the shell
    vel_h = vel_h_shaft-vel_h_shell                                                                         ! horizontal velocity of the shaft relative to the shell, still in the reference frame of the inertial system
    vel_v = vel_v_shaft-vel_v_shell                                                                         ! vertical velocity of the shaft relative to the shell, still in the reference frame of the inertial system
    q_dot = (vel_v*dis_v+vel_h*dis_h)/q                                                                     ! rate of change of absolute eccentricity
    X_att_dot = (vel_v*dis_h-vel_h*dis_v)/q**2                                                              ! rate of change of attitude angle in the reference frame of the inertial system
    tilt_dot_h = tilt_dot_h_shaft-tilt_dot_h_shell                                                          ! rate of change of the tilting angle of the shaft relative to the shell around the horizontal axis, still in the reference frame of the inertial system
    tilt_dot_v = tilt_dot_v_shaft-tilt_dot_v_shell                                                          ! rate of change of the tilting angle of the shaft relative to the shell around the vertical axis, still in the reference frame of the inertial system
    tilt_dot = (tilt_dot_v*tilt_v+tilt_dot_h*tilt_h)/tilt                                                   ! rate of change of the absolute tilting angle
    X_tilt_dot = (tilt_dot_v*tilt_h-tilt_dot_h*tilt_v)/tilt**2                                              ! rate of change of the angle describing the resulting tilting axis in the reference frame of the inertial system
    X_att_dot = X_att_dot-omega_shell                                                                       ! rate of change of attitude angle in the shell-fixed reference frame in which the Reynolds equation is solved
    X_tilt_dot = X_tilt_dot-omega_shell                                                                     ! rate of change of the angle describing the tilting axis in the shell-fixed reference frame in which the Reynolds equation is solved
  ELSE
    omega = omega_shaft                                                                                     ! quasistatic case: rotational velocity of shell is assumed to be zero
    q_dot = 0                                                                                               ! quasistatic case: rate of change of eccentricity is assumed to be zero
    X_att_dot = 0                                                                                           ! quasistatic case: rate of change of attitude angle is assumed to be zero
    tilt_dot = 0                                                                                            ! quasistatic case: rate of change of tilting angle is assumed to be zero
    X_tilt_dot = 0                                                                                          ! quasistatic case: rate of change of the angle describing the tilting axis is assumed to be zero
  END IF
  u = omega*(d_b/2.0d0)                                                                                     ! circumferential surface velocity of the shaft in the shell-fixed reference frame in which the Reynolds equation is solved (the surface velocity of the shell is zero in this reference frame) [m/s]
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Analyze grid
  ! --------------------------------------------------------------------------------------------------------
  
  pi = 3.14159265359d0                                                                                      ! define pi
  r_b = d_b/2                                                                                               ! bearing radius [m]
  L_X = 2*pi/n_x                                                                                            ! angular circumferential side length of control volume [rad]
  L_Y = (l_b/(n_y-1))/r_b                                                                                   ! nondimensionalized axial side length of control volume [-]
  IF ( symBC .EQ. 0 ) THEN                                                                                  ! if no symmetric BC is used
    dof_y = n_y-2                                                                                           ! number of DOFs in the axial direction (nummber of nodes without counting the nodes at the bearing boundaries)
    j_start = 2                                                                                             ! axial node number corresponding to the first DOF in the axial direction
  ELSE                                                                                                      ! if a symmetric BC is used
    IF ( MODULO(n_y,2) .EQ. 0 ) THEN                                                                        ! if the node number is even (no node is located at the axial center of the bearing Y=0)
      dof_y = n_y/2-1                                                                                       ! number of DOFs in the axial direction (one half of the bearing, not counting the node at the bearing boundary)
    ELSE                                                                                                    ! if the node number is uneven (a node is located at the axial center of the bearing Y=0)
      dof_y = (n_y-1)/2                                                                                     ! number of DOFs in the axial direction (one half of the bearing, not counting the node at the bearing boundary)
    END IF
    j_start = n_y-dof_y                                                                                     ! axial node number corresponding to the first DOF in the axial direction
  END IF
  n = n_x*n_y                                                                                               ! total number of nodes (whole bearing, nodes with Dirichlet BCs are included) [-]
  dof = dof_y*n_x                                                                                           ! total number of DOFs, including the nodes with oil supply BCs (because the penalty method is used) but not the nodes at the bearing boundaries
  n_nz = n_x*(5*(dof_y-2)+8)                                                                                ! number of nonzero matrix entries (or rather, number of stored matrix entries, in this case)
  IF ( guembel .EQ. 1 ) THEN
    n_nz = (n_nz-dof)/2+dof                                                                                 ! under Guembel conditions, the matrix is symmetric, meaning that fewer entries need to be stored
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Settings for the BiCGStab solver
  ! --------------------------------------------------------------------------------------------------------
  
  tol_fil = 1.0d-4                                                                                          ! fill-in threshold (smaller values are ignored) for ilu
  lfil = 5                                                                                                  ! max allowed number of fill-ins per row in ilu
  lws = 8*dof                                                                                               ! size of workspace array
  liws = 2*dof                                                                                              ! size of integer workspace array
  ipar_vec = 0
  ipar_vec(2) = 1                                                                                           ! status of the preconditioning: (0,1,2,3) == ( no prec., left prec. only, right prec. only, both prec.)
  ipar_vec(3) = 2                                                                                           ! use residual-based convergence criterion: || residual || <= rtol * || rhs || + atol
  ipar_vec(4) = lws                                                                                         ! size of workspace array
  ipar_vec(6) = iter_max_solver                                                                             ! maximum number of matrix-vector multiplications
  fpar_vec = 0.0d0
  fpar_vec(1) = tol                                                                                         ! relative tolerance
  fpar_vec(2) = 0.0d0                                                                                       ! absolute tolerance (if we set this to zero, a pure relative tolerance will be used))
  nnz_ilu = n_nz + 2*dof*lfil                                                                               ! upper estimate of the overall number of nonzeros after fill-in 
      
  
  ! --------------------------------------------------------------------------------------------------------
  ! Some allocations
  ! --------------------------------------------------------------------------------------------------------
  
  ALLOCATE(nn0_vec(dof))                                                                                    ! this array will later contain the numbers of the nodes that are DOFs
  ALLOCATE(nnN_vec(dof))                                                                                    ! for every DOF, this array will contain the number of the neighboring node N in the positive Y-direction
  ALLOCATE(nnE_vec(dof))                                                                                    ! for every DOF, this array will contain the number of the neighboring node E in the positive X-direction
  ALLOCATE(nnS_vec(dof))                                                                                    ! for every DOF, this array will contain the number of the neighboring node S in the negative Y-direction
  ALLOCATE(nnW_vec(dof))                                                                                    ! for every DOF, this array will contain the number of the neighboring node W in the negative X-direction
  ALLOCATE(aN_vec(dof))                                                                                     ! array for storing the dimensionless conductivities in the positive Y-direction (interaction with DOF N) for every DOF
  ALLOCATE(aE_vec(dof))                                                                                     ! array for storing the dimensionless conductivities in the positive X-direction (interaction with DOF E) for every DOF
  ALLOCATE(aS_vec(dof))                                                                                     ! array for storing the dimensionless conductivities in the negative Y-direction (interaction with DOF S) for every DOF
  ALLOCATE(aW_vec(dof))                                                                                     ! array for storing the dimensionless conductivities in the negative X-direction (interaction with DOF W) for every DOF
  ALLOCATE(coef_vec(5*dof))                                                                                 ! array for storing all matrix entries (before they are ordered accoring to the csr format)
  ALLOCATE(coefIndex_vec(n_nz))                                                                             ! this array will be used to organize the order in which the matrix entries are stored in values_vec (i.e. how their order is changed when they are copied from coef_vec to values_vec) to satisfy the csr format
  ALLOCATE(values_vec(n_nz))                                                                                ! array for storing all matrix entries
  ALLOCATE(columns_vec(n_nz))                                                                               ! array for storing the columns of all matrix entries
  ALLOCATE(rowIndex_vec(dof+1))                                                                             ! array for storing the indices of values_vec corresponding to the matrix entries that appear first in their respective matrix rows
  ALLOCATE(md_vec(dof))                                                                                     ! array for storing the indices of values_vec corresponding to the main diagonal entries (we will need this for the oil supply BCs)
  ALLOCATE(R_vec(dof))                                                                                      ! array for storing the RHS vector
  ALLOCATE(initial_guess_vec(dof))                                                                          ! array for storing the initial guess for the iterative solver
  ALLOCATE(Pi_dof_vec(dof))                                                                                 ! array for storing the solution vector (pressure-like function at the DOFs)
  ALLOCATE(a_lu_vec(nnz_ilu))                                                                               ! for ILU preconditioner
  ALLOCATE(j_lu_vec(nnz_ilu))                                                                               ! for ILU preconditioner
  ALLOCATE(j_u_vec(dof))                                                                                    ! for ILU preconditioner
  ALLOCATE(ws_vec(lws))                                                                                     ! workspace array for solver and ILU
  ALLOCATE(iws_vec(liws))                                                                                   ! workspace array for ILU
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Define some arrays organizing the node numbers and the nodal coordinate values
  ! --------------------------------------------------------------------------------------------------------
  
  nn0_vec = (/(i,i=(j_start-1)*n_x+1,(n_y-1)*n_x,1)/)                                                       ! array containing the numbers of the nodes which are DOFs
  nnN_vec = nn0_vec+n_x                                                                                     ! array containing, for every DOF, the number of the neighboring node N in the positive Y-direction
  nnE_vec = nn0_vec+1                                                                                       ! array containing, for every DOF, the number of the neighboring node E in the positive X-direction
  nnS_vec = nn0_vec-n_x                                                                                     ! array containing, for every DOF, the number of the neighboring node S in the negative Y-direction
  nnW_vec = nn0_vec-1                                                                                       ! array containing, for every DOF, the number of the neighboring node W in the negative X-direction
  IF ( ( symBC .EQ. 1 ) .AND. ( MODULO(n_y,2) .EQ. 0 ) ) THEN                                               ! if a symmetric BC is used and an even axial number of nodes has been prescribed for the entire bearing length (meaning that no nodes exist at the symmetric boundary)
    nnS_vec(1:n_x) = nn0_vec(1:n_x)                                                                         ! this manipulation ensures that, at the symmetric boundary, several nodal values (switch function, gap function, and viscosity) at S will be replaced by those at 0
  ELSEIF ( ( symBC .EQ. 1 ) .AND. ( MODULO(n_y,2) .EQ. 1 ) ) THEN                                           ! if a symmetric BC is used and an uneven axial number of nodes has been prescribed for the entire bearing length (meaning that nodes exist at the symmetric boundary)
    nnS_vec(1:n_x) = nnN_vec(1:n_x)                                                                         ! this manipulation ensures that, at the symmetric boundary, several nodal values (switch function, gap function, and viscosity) at S will be replaced by those at N
  END IF
  FORALL (j=1:dof_y) nnE_vec(j*n_x) = nnE_vec(j*n_x) - n_x                                                  ! correct some node numbers in nnE_vec according to the periodicity of the computational domain
  FORALL (j=0:(dof_y-1)) nnW_vec(j*n_x+1) = nnW_vec(j*n_x+1) + n_x                                          ! correct some node numbers in nnW_vec according to the periodicity of the computational domain
  X_vec(1:n_x) = ((/(i,i=0,n_x-1,1)/)*L_X)                                                                  ! nodal X-positions for one row of nodes in the circumferential direction
  FORALL (i=1:(n_y-1)) X_vec((i*n_x+1):((i+1)*n_x)) = X_vec(1:n_x)                                          ! definition of X_vec: this array contains all nodal X-positions (angular circumferential positions)
  FORALL (j=0:(n_y-1)) Y_vec((j*n_x+1):((j+1)*n_x)) = ((j*1.0d0)/(n_y-1)-0.5)*(l_b/r_b)                     ! definition of Y_vec: this array contains all nodal Y-positions (nondimensionalized axial positions)
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Analyze oil supply
  ! --------------------------------------------------------------------------------------------------------
  
  X_os_vec = (/(i,i=0,grooves-1,1)/)*(2*pi/grooves) + X_os                                                  ! X-values (circumferential coordinate values) describing the positions of the oil supply grooves [rad]
  start_os_vec = ANINT((X_os_vec-L_X_os/2)/L_X+1)                                                           ! circumferential node/DOF numbers where the oil supply grooves begin
  end_os_vec = ANINT((X_os_vec+L_X_os/2)/L_X+1)                                                             ! circumferential node/DOF numbers where the oil supply grooves end
  n_os_x = SUM(end_os_vec-start_os_vec+1)                                                                   ! circumferential number of nodes/DOFs where oil supply BCs will be defined
  ALLOCATE(os_x_vec(n_os_x))                                                                                ! this array serves for storing the circumferential node/DOF numbers where oil supply BCs will be defined
  j = 0                                                                                                     ! j serves to track the number of already considered nodes
  DO i = 1, grooves                                                                                         ! loop through all node numbers
    n_i = end_os_vec(i)-start_os_vec(i)+1                                                                   ! circumferential number of nodes in current groove
    os_x_vec((j+1):(j+n_i)) = (/(k,k=start_os_vec(i),end_os_vec(i),1)/)                                     ! save node numbers
    j = j+n_i                                                                                               ! update number of considered nodes
  END DO
  DO i = 1, n_os_x                                                                                          ! loop through all circumferential node/DOF numbers where oil supply BCs will be defined
    os_x_vec(i) = MODULO(os_x_vec(i)-1,n_x) + 1                                                             ! correct the circumferential node/DOF number if it's currently outside the range 1...nx
  END DO
  IF ( l_y_os .LE. 0.0d0 ) THEN                                                                             ! if the axial side length of the groove is less than or equal to zero
    a = 1.0d-10                                                                                             ! the variable a is now used for storing the axial side length of the groove (because l_y_os is not editable) and this length is set to a value slightly larger than zero (exactly zero causes difficulties if the axial number of nodes is even)
  ELSE                                                                                                      ! if the axial side length of the groove is positive
    a = l_y_os                                                                                              ! the axial side length of the groove is now stored in the variable a
  END IF
  start_os = ANINT( ((l_b-a)/2)/l_b*(n_y-1)+1 )                                                             ! axial node number where the groove begins
  IF ( start_os .LT. 2 ) THEN                                                                               ! if the groove begins at the bearing boundary (start_os = 1) or, illegaly, extends further than the bearing (start_os < 1)
    longgroove = 1                                                                                          ! set longgroove = 1 as a flag which will later allow to recover this information
    start_os = 2                                                                                            ! shift the axial starting index of the groove one node away from the bearing boundary (because the nodes located exactly at the bearing boundary have already been excluded from the DOFs)
  ELSE
    longgroove = 0
  END IF
  end_os = n_y-start_os+1                                                                                   ! axial node number where the groove ends
  IF ( symBC .EQ. 1 ) THEN                                                                                  ! if a symmetric BC is used
    n_os_y = end_os-j_start+1                                                                               ! axial number of DOFs where oil supply BCs will be defined
    ALLOCATE(os_y_vec(n_os_y))
    os_y_vec = (/(i,i=1,end_os-j_start+1,1)/)                                                               ! list of axial DOF numbers where oil supply BCs will be defined (not equal to the corresponding node numbers, due to the axial shift according to j_start)
  ELSE                                                                                                      ! if no symmetric BC is used
    n_os_y = end_os-start_os+1                                                                              ! axial number of DOFs where oil supply BCs will be defined
    ALLOCATE(os_y_vec(n_os_y))
    os_y_vec = (/(i,i=start_os-j_start+1,end_os-j_start+1,1)/)                                              ! list of axial DOF numbers where oil supply BCs will be defined (not equal to the corresponding node numbers, due to the axial shift according to j_start)
  END IF
  ALLOCATE(os_vec(n_os_x*n_os_y))
  FORALL (j=1:n_os_y) os_vec((j-1)*n_os_x+1:j*n_os_x) = os_x_vec + (os_y_vec(j)-1)*n_x                      ! definition of os_vec: this array contains the DOF numbers where oil supply BCs will be defined (now using a global numbering scheme that considers both directions)
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Nondimensionalization
  ! --------------------------------------------------------------------------------------------------------
  
  sgn_u = SIGN(1.0d0,u)                                                                                     ! sign of circumferential surface velocity [-]
  mu_ref = SUM(mu_vec)/n                                                                                    ! reference viscosity [Pa*s]
  p_ref = ABS(u)*mu_ref*r_b/(2*c**2)                                                                        ! reference pressure for nondimensionalization [Pa]
  epsil = q/c                                                                                               ! relative eccentricity [-]
  epsil_dot = q_dot/c                                                                                       ! rate of change of the relative eccentricity [1/s]
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Analyze boundary value prescribed in oil supply groove
  ! --------------------------------------------------------------------------------------------------------
  
  IF ( p_os .GE. 0.0d0 ) THEN                                                                               ! "IF a non-negative supply pressure is prescribed, then ..."
    Pi_os = p_os/p_ref                                                                                      ! nondimensionalize the supply pressure and save as pressure-like function
    g_os = 1                                                                                                ! set corresponding switch function to 1
  ELSE                                                                                                      ! otherwise, the prescribed value p_os will be interpreted as a film fraction with an offset of -1, i.e., as a negative pressure-like function
    Pi_os = p_os                                                                                            ! save prescribed value as pressure-like function
    g_os = 0                                                                                                ! set corresponding switch function to 0
  END IF
  IF ( ( guembel .EQ. 1 ) .AND. ( p_os .LT. 0.0d0 ) ) THEN                                                  ! if Guembel is used instead of Elrod, we cannot prescribe a film fraction as BC (which would be done via a negative p_os), so ...
    Pi_os = 0.0d0                                                                                           ! ... make sure that the oil supply pressure is at least zero
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Analyze data from previous time step
  ! --------------------------------------------------------------------------------------------------------
  
  Pi_pts_vec = pts_vec(1:n)                                                                                 ! extract the previous solution (pressure-like function of the previous time step) from pts_vec
  Pi_dof_vec(:) = Pi_pts_vec(nn0_vec)                                                                       ! solution from previous time step will be provided as initial guess for the iterative solver
  Pi_vec = 0.0d0                                                                                            ! initialize solution as zero; note that Pi_vec is not the solution vector passed to the solver (that would be Pi_dof_vec), so we don't need to worry about providing an initial guess here
  g_pts_vec = 1                                                                                             ! ... initialize the switch function from the previous time step as one, but ...
  IF ( guembel .EQ. 0 ) THEN                                                                                ! ... if Elrod cavitation is assumed, then ...
    WHERE(Pi_pts_vec .LT. 0.0d0) g_pts_vec = 0                                                              ! ... set this switch function to 0 at all nodes that, according to the previous solution, were part of the pressure zone
  END IF
  g_vec = g_pts_vec                                                                                         ! use the switch function of the previous time step as initial guess for the current switch function
  g_vec(os_vec) = g_os                                                                                      ! ensure that the switch functions in the oil supply groove are consistent with the prescribed oil supply pressure or oil supply film fraction
  g_old_vec = 0                                                                                             ! switch function of previous iteration (does not exist yet)
  IF ( quasistatic .EQ. 1 ) THEN                                                                            ! IF a quasistatic simulation is desired, then ...
    Delta_T = ABS(u)/(2*(d_b/2))*1.0d16                                                                     ! ... set the dimensionless time increment to an extremely large value, suppressing transient effects
  ELSE                                                                                                      ! IF a transient simulation is desired, then ...
    t_pts = pts_vec(n+1)                                                                                    ! ... extract the time of the previous time step from pts_vec
    Delta_T = ABS(u)/(2*(d_b/2))*(t-t_pts)                                                                  ! ... compute the dimensionless time increment
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Analyze gap function and its rate of change
  ! --------------------------------------------------------------------------------------------------------
  
  cos_vec(1:n_x) = COS(X_vec(1:n_x)-X_att)                                                                  ! cos(X-X_att) for one row of nodes in the circumferential direction, where X is the respective circumferential nodal position
  sin_vec(1:n_x) = SIN(X_vec(1:n_x)-X_att)                                                                  ! sin(X-X_att) for one row of nodes in the circumferential direction, where X is the respective circumferential nodal position
  FORALL (i=1:(n_y-1)) cos_vec((i*n_x+1):((i+1)*n_x)) = cos_vec(1:n_x)                                      ! cos(X-X_att) for all nodes
  FORALL (i=1:(n_y-1)) sin_vec((i*n_x+1):((i+1)*n_x)) = sin_vec(1:n_x)                                      ! sin(X-X_att) for all nodes
  H_vec = 1 - epsil*cos_vec    + ac_vec/c                                                                   ! nondimensionalized gap function without shaft tilting
  IF ( MINVAL(H_vec) .LE. 0.0d0 ) THEN                                                                      ! if the minimum gap width isn't positive
    PRINT *, 'warning: mininmum gap width is not positive'                                                  ! show a warning
  END IF
  dHdT_vec = (d_b/ABS(u))*(-epsil_dot*cos_vec-epsil*X_att_dot*sin_vec)                                      ! nondimensionalized time derivative of the nondimensionalized gap function without shaft tilting
  IF ( symBC .EQ. 0 ) THEN                                                                                  ! if no symmetric BC is used (i.e., if the consideration of shaft tilting is enabled)
    cos_vec(1:n_x) = COS(X_vec(1:n_x)-X_tilt)                                                               ! cos(X-X_tilt) for one row of nodes in the circumferential direction, where X is the respective circumferential nodal position
    sin_vec(1:n_x) = SIN(X_vec(1:n_x)-X_tilt)                                                               ! sin(X-X_tilt) for one row of nodes in the circumferential direction, where X is the respective circumferential nodal position
    FORALL (i=1:(n_y-1)) cos_vec((i*n_x+1):((i+1)*n_x)) = cos_vec(1:n_x)                                    ! cos(X-X_tilt) for all nodes
    FORALL (i=1:(n_y-1)) sin_vec((i*n_x+1):((i+1)*n_x)) = sin_vec(1:n_x)                                    ! sin(X-X_tilt) for all nodes
    H_vec = H_vec - (tilt*r_b/c)*Y_vec*sin_vec                                                              ! nondimensionalized gap function: consider shaft tilting 
    dHdT_vec = dHdT_vec + &                                                                                 ! nondimensionalized time derivative of the nondimensionalized gap function: consider shaft tilting 
      ((d_b/ABS(u))*(r_b/c))*Y_vec*(-tilt_dot*sin_vec+(tilt*X_tilt_dot)*cos_vec)
  END IF
  H3_vec = H_vec**3                                                                                         ! cubed nondimensionalized gap function at all nodes
  H3_over_mu_vec = (H3_vec/mu_vec)*mu_ref                                                                   ! cubed nondimensionalized gap function divided by nondimensionalized viscosity at all nodes
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Compute dimensionless conductivities
  ! --------------------------------------------------------------------------------------------------------
  
  aN_vec = (L_X/(24*L_Y))*(H3_over_mu_vec(nn0_vec)+H3_over_mu_vec(nnN_vec))                                 ! dimensionless conductivities in the positive Y-direction (interaction with DOF N) for every DOF
  aE_vec = (L_Y/(24*L_X))*(H3_over_mu_vec(nn0_vec)+H3_over_mu_vec(nnE_vec))                                 ! dimensionless conductivities in the positive X-direction (interaction with DOF E) for every DOF
  aS_vec = (L_X/(24*L_Y))*(H3_over_mu_vec(nn0_vec)+H3_over_mu_vec(nnS_vec))                                 ! dimensionless conductivities in the negative Y-direction (interaction with DOF S) for every DOF
  aW_vec = (L_Y/(24*L_X))*(H3_over_mu_vec(nn0_vec)+H3_over_mu_vec(nnW_vec))                                 ! dimensionless conductivities in the negative X-direction (interaction with DOF W) for every DOF  
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Some preparations for the assembly
  ! --------------------------------------------------------------------------------------------------------
  
  IF ( guembel .EQ. 0 ) THEN
    CALL PREPARE_ASSEMBLY_ELROD(n_x, dof_y, dof, n_nz, md_vec, rowIndex_vec, coefIndex_vec, columns_vec)    ! prepare sparse matrix assembly under Elrod conditions (this subroutine is defined below)
  ELSE
    CALL PREPARE_ASSEMBLY_GUEMB(n_x, dof_y, dof, n_nz, rowIndex_vec, coefIndex_vec, columns_vec)            ! prepare sparse matrix assembly under Guembel conditions, i.e., only for upper triangular matrix (this subroutine is defined below)
    md_vec = rowIndex_vec(1:dof)                                                                            ! this array points to the main diagonal entries (first entry of every row, because only the upper triangle will be stored)
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Execute iterative scheme
  ! --------------------------------------------------------------------------------------------------------
  
  iter_sol = 0                                                                                              ! for counting the number of iterations required by the solver (i.e., by BiCGStab or CG), accumulated over all calls of the solver
  convergent = 0                                                                                            ! this flag will be set to 1 once the solution has converged 
  
  DO iter = 1, iter_max                                                                                     ! fixed-point iteration for linearization of the nonlinear boundary value problem (BVP); iter counts the iterations
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Compute coefficients on the left-hand side of the discretized Reynolds equation
    ! ------------------------------------------------------------------------------------------------------
    
    coef_vec(1:dof) = (aE_vec+aW_vec+aN_vec+aS_vec)*g_vec(nn0_vec) + &                                      ! in this space of the array coef_vec, we store the coefficient corresponding to the DOF 0 for every instance of the discretized Reynolds equation
      L_Y*H_vec(nn0_vec)*(1-g_vec(nn0_vec)) + &    
      (L_X*L_Y)*(dHdT_vec(nn0_vec)+H_vec(nn0_vec)/Delta_T)*(1-g_vec(nn0_vec))
    coef_vec((2*dof+1):(3*dof)) = -aE_vec*g_vec(nnE_vec) + &                                                ! in this space of the array coef_vec, we store the coefficient corresponding to the DOF E for every instance of the discretized Reynolds equation
      ((sgn_u-1)*(L_Y/2))*H_vec(nnE_vec)*(1-g_vec(nnE_vec))
    coef_vec((4*dof+1):(5*dof)) = -aW_vec*g_vec(nnW_vec) + &                                                ! in this space of the array coef_vec, we store the coefficient corresponding to the DOF W for every instance of the discretized Reynolds equation
      ((-sgn_u-1)*(L_Y/2))*H_vec(nnW_vec)*(1-g_vec(nnW_vec))
    coef_vec((dof+1):(2*dof)) = -aN_vec*g_vec(nnN_vec)                                                      ! in this space of the array coef_vec, we store the coefficient corresponding to the DOF N for every instance of the discretized Reynolds equation
    coef_vec((3*dof+1):(4*dof)) = -aS_vec*g_vec(nnS_vec)                                                    ! in this space of the array coef_vec, we store the coefficient corresponding to the DOF S for every instance of the discretized Reynolds equation
    IF ( ( symBC .EQ. 1 ) .AND. ( MODULO(n_y,2) .EQ. 0 ) ) THEN                                             ! if a symmetric BC is used and an even axial number of nodes has been prescribed for the entire bearing length (meaning that no nodes exist at the symmetric boundary)
      coef_vec(1:n_x) = coef_vec(1:n_x) + coef_vec((3*dof+1):(3*dof+n_x))                                   ! at the symmetric boundary, the DOF S will be substituted by the DOF 0; to facilitate this substitution, we add the entries in coefS_vec(1:n_x) to coef0_vec(1:n_x); after that, anything stored in coefS_vec(1:n_x) is obsolete
    ELSEIF ( ( symBC .EQ. 1 ) .AND. ( MODULO(n_y,2) .EQ. 1 ) ) THEN                                         ! if a symmetric BC is used and an uneven axial number of nodes has been prescribed for the entire bearing length (meaning that nodes exist at the symmetric boundary)
      coef_vec((dof+1):(dof+n_x)) = coef_vec((dof+1):(dof+n_x)) + coef_vec((3*dof+1):(3*dof+n_x))           ! at the symmetric boundary, the DOF S will be substituted by the DOF N; to facilitate this substitution, we add the entries in coefS_vec(1:n_x) to coefN_vec(1:n_x); after that, anything stored in coefS_vec(1:n_x) is obsolete
    END IF
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Assembly
    ! ------------------------------------------------------------------------------------------------------
    
    R_vec(:) = (-sgn_u*L_Y/2)*(H_vec(nnE_vec)-H_vec(nnW_vec)) - (L_X*L_Y)*dHdT_vec(nn0_vec) + &             ! array containing the RHS vector
      (L_X*L_Y/Delta_T)*H_vec(nn0_vec)*(1-g_pts_vec(nn0_vec))*Pi_pts_vec(nn0_vec)
    values_vec = coef_vec(coefIndex_vec)                                                                    ! array containing the matrix entries (the assembly is that simple because of the preperations made by PREPARE_ASSEMBLY earlier)
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Some modifications due to boundary conditions
    ! ------------------------------------------------------------------------------------------------------
    
    IF ( ( symBC .EQ. 1 ) .AND. ( MODULO(n_y,2) .EQ. 1 ) ) THEN                                             ! if a symmetric BC is used and an uneven axial number of nodes has been prescribed for the entire bearing length (meaning that nodes exist at the symmetric boundary)
      R_vec(1:n_x) = R_vec(1:n_x)/2                                                                         ! for symmetric BC: the control volumes at the symmetric boundary are only half the regular size, which is hereby considered in the RHS vector
      values_vec(1:(rowIndex_vec(n_x+1)-1)) = values_vec(1:(rowIndex_vec(n_x+1)-1))/2                       ! for symmetric BC: the control volumes at the symmetric boundary are only half the regular size, which is hereby considered in the matrix entries
    END IF
    pf = pm*MAXVAL(ABS(values_vec(md_vec)))                                                                 ! for oil supply BCs: penalty factor equal to 100 times the maximum absolute main diagonal entry
    R_vec(os_vec) = pf*Pi_os                                                                                ! for oil supply BCs: modify RHS vector (penalty method)
    values_vec(md_vec(os_vec)) = values_vec(md_vec(os_vec)) + pf                                            ! for oil supply BCs: modify matrix entries (penalty method)
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Solve system of equations
    ! ------------------------------------------------------------------------------------------------------
    
    IF ( guembel .EQ. 0 ) THEN
      error = 0
      CALL ILUT(dof, values_vec, columns_vec, rowIndex_vec, lfil, tol_fil, a_lu_vec, j_lu_vec, j_u_vec, &   ! ILU (incomplete lower-upper) factorization, which will be used as preconditioner
        nnz_ilu, ws_vec(1:dof+1), iws_vec, error)
      IF ( error .NE. 0 ) THEN                                                                              ! if an error has occurred
        PRINT *, 'ILU error'
        PRINT *, error                                                                                      ! print error type
      END IF
      initial_guess_vec = Pi_dof_vec                                                                        ! choose initial guess according to previous solution
      CALL RUNRC(dof, R_vec, Pi_dof_vec, ipar_vec, fpar_vec, ws_vec, initial_guess_vec, values_vec, &       ! solve system of equations by BiCGStab iteration
        columns_vec, rowIndex_vec, a_lu_vec, j_lu_vec, j_u_vec, BCGSTAB)  
      iter_out = ipar_vec(7)                                                                                ! number of iterations required by solver
    ELSE
      CALL SOLVER_IC_CG_MKL(dof, n_nz, values_vec, columns_vec, rowIndex_vec, R_vec, Pi_dof_vec, &          ! under Guembel conditions, use ICCG instead of BiCGStab
        iter_out, tol, iter_max_solver)
    END IF
    Pi_vec(nn0_vec) = Pi_dof_vec(:)                                                                         ! save solution
    iter_sol = iter_sol + iter_out                                                                          ! count the number of iterations required by the solver (i.e., by BiCGStab or CG), accumulated over all calls of the solver
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Update switch function and check convergence
    ! ------------------------------------------------------------------------------------------------------
    
    g_old_old_vec = g_old_vec
    g_old_vec = g_vec
    g_vec = 0                                                                                               ! the new switch function (which will be used in the next iteration) is computed: it is first initialized as zero, then ...
    WHERE(Pi_vec .GE. 0.0d0) g_vec = 1                                                                      ! ... it is set to 1 where the axial average of the pressure-like function Pi_bar is non-negative
    IF (MAXVAL(ABS(g_vec(nn0_vec)-g_old_vec(nn0_vec))) .EQ. 0) THEN                                         ! if the new switch function is equal to the old one, then
      convergent = 1                                                                                        ! flag the solution as converged
    END IF
    IF ( (iter .GE. 2) .AND. (MAXVAL(ABS(g_vec(nn0_vec)-g_old_old_vec(nn0_vec))) .EQ. 0) ) THEN             ! if the same switch function has been computed in this iteration as two iterations before (indicating that the solution oscillates between two alternating configurations - although I doubt that this even happens under the current model assumptions), then
      convergent = 1                                                                                        ! flag the solution as converged anyway (the solution is as good as it gets)
    END IF
    IF ( (convergent .EQ. 1) .OR. (iter .EQ. iter_max) .OR. (guembel .EQ. 1)  ) THEN                        ! if this is the last iteration (solution converged or max allowed number of iterations reached or Guembel conditions assumed), then
      EXIT                                                                                                  ! exit the iterative scheme
    END IF
    
    
  END DO
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Postprocessing
  ! --------------------------------------------------------------------------------------------------------
  
  IF ( longgroove .EQ. 1 ) THEN                                                                             ! if the oil supply groove extends across the entire bearing length (simplification)
    Pi_vec(os_x_vec) = Pi_os                                                                                ! replace atmospheric pressure with supply pressure at the corresponding nodes at the first bearing boundary
    Pi_vec(os_x_vec+(n_y-1)*n_x) = Pi_os                                                                    ! replace atmospheric pressure with supply pressure at the corresponding nodes at the second bearing boundary
  END IF
  IF ( symBC .EQ. 1 ) THEN                                                                                  ! if a symmetric BC was used
    FORALL (j=1:(j_start-1)) Pi_vec((1+(j-1)*n_x):(j*n_x)) = Pi_vec((1+(n_y-j)*n_x):((n_y-j+1)*n_x))        ! the solution is mirrored to the other half of the bearing
  END IF
  FORALL (j=1:n_y) Pi_mat(:,j) = Pi_vec((1+(j-1)*n_x):(j*n_x))                                              ! the solution (pressure-like function) is stored in a matrix for output (for visualization purposes)
  g_vec = 0
  WHERE(Pi_vec .GE. 0.0d0) g_vec = 1                                                                        ! the switch function is computed again for all nodes and stored in an array
  theta_vec = (1-g_vec)*Pi_vec + 1                                                                          ! the film fraction is computed at all nodes and stored in an array
  IF ( guembel .EQ. 1 ) THEN                                                                                ! if Guembel is used instead of Elrod, ...
    theta_vec = 1                                                                                           ! ... set film fractions to 1
  END IF
  p_vec = p_ref*Pi_vec*g_vec                                                                                ! the physical (as opposed to nondimensionalized) pressure is computed at all nodes and stored in an array
  cos_vec(1:n_x) = COS(X_vec(1:n_x))                                                                        ! cos(X) for one row of nodes in the circumferential direction, where X is the respective circumferential nodal position
  sin_vec(1:n_x) = SIN(X_vec(1:n_x))                                                                        ! sin(X) for one row of nodes in the circumferential direction, where X is the respective circumferential nodal position
  FORALL (i=1:(n_y-1)) cos_vec((i*n_x+1):((i+1)*n_x)) = cos_vec(1:n_x)                                      ! cos(X) for all nodes
  FORALL (i=1:(n_y-1)) sin_vec((i*n_x+1):((i+1)*n_x)) = sin_vec(1:n_x)                                      ! sin(X) for all nodes
  L_X = L_X*r_b                                                                                             ! the circumferential control volume side length is converted from [rad] to [m]
  L_Y = L_Y*r_b                                                                                             ! the axial control volume side length is converted from [-] to [m]
  A_i = L_X*L_Y                                                                                             ! control volume area [m^2]
  A_vec = A_i                                                                                               ! in order to weight every nodal solution by the control volume area (for integration purposes, see below), an array with one entry per node is defined containing the control volume area
  A_vec(1:n_x) = A_i/2                                                                                      ! half control volumes are present at the first bearing boundary
  A_vec((1+(n_y-1)*n_x):(n_y*n_x)) = A_i/2                                                                  ! half control volumes are present at the second bearing boundary
  F_1 = SUM(cos_vec*p_vec*A_vec)                                                                            ! component of the hydrodynamic force (computed by integration of the pressure, i.e., summation of the nodal forces) [N]
  F_2 = SUM(sin_vec*p_vec*A_vec)                                                                            ! component of the hydrodynamic force (computed by integration of the pressure, i.e., summation of the nodal forces) [N]
  M_1 = SUM(sin_vec*p_vec*A_vec*Y_vec)*r_b                                                                  ! component of the hydrodynamic moment [Nm]
  M_2 = -SUM(cos_vec*p_vec*A_vec*Y_vec)*r_b                                                                 ! component of the hydrodynamic moment [Nm]
  F_h = F_1*COS(angle_shell) - F_2*SIN(angle_shell);                                                        ! horizontal hydrodynamic force acting on the shell [N]
  F_v = F_1*SIN(angle_shell) + F_2*COS(angle_shell);                                                        ! vertical hydrodynamic force acting on the shell [N]
  M_h = M_1*COS(angle_shell) - M_2*SIN(angle_shell);                                                        ! horizontal hydrodynamic force acting on the shell [Nm]
  M_v = M_1*SIN(angle_shell) + M_2*COS(angle_shell);                                                        ! vertical hydrodynamic force acting on the shell [Nm]
  M_fr = 0.0d0                                                                                              ! the friction moment M_fr will be computed below; it is initialized as zero [Nm]
  DO j = 1, n_y                                                                                             ! loop through the nodes in the axial direction (to evaluate the respective contributions of the nodal pressure gradients to the friction moment)
    M_fr = M_fr - (p_vec(2+(j-1)*n_x)-p_vec(j*n_x)) * H_vec(1+(j-1)*n_x) * A_vec(1+(j-1)*n_x)               ! contribution of the first node at the current axial position (the pressure gradient is evaluated under consideration of the periodicity of the domain)
    M_fr = M_fr - (p_vec(1+(j-1)*n_x)-p_vec(j*n_x-1)) * H_vec(j*n_x) * A_vec(j*n_x)                         ! contribution of all nodes at the current axial position except for the first one and the last one
    M_fr = M_fr - SUM( (p_vec((3+(j-1)*n_x):(j*n_x))-p_vec((1+(j-1)*n_x):(j*n_x-2))) * &                    ! contribution of the last node at the current axial position (the pressure gradient is evaluated under consideration of the periodicity of the domain)
      H_vec((2+(j-1)*n_x):(j*n_x-1)) * A_vec((2+(j-1)*n_x):(j*n_x-1)) )
  END DO
  M_fr = M_fr*(c*r_b/(4*L_X))                                                                               ! some factors have been left out in the summands of M_fr in the loop above; these are now considered
  M_fr = M_fr + SUM( (theta_vec*mu_vec*A_vec)/H_vec ) * (u*r_b/c)                                           ! contribution of the nodal film fractions to the friction moment
  V_oil = SUM(theta_vec*H_vec*A_vec)*c                                                                      ! oil volume in the fluid film under consideration of the film fraction [m^3]
  V_dot_bb = SUM( (p_vec(1:n_x)-p_vec((n_x+1):(2*n_x)))*(H3_vec(1:n_x)/mu_vec(1:n_x)) + &                   ! oil volume flow through the bearing boundaries [m^3/s]
    (p_vec((n-n_x+1):n)-p_vec((n-2*n_x+1):(n-n_x)))*(H3_vec((n-n_x+1):n)/mu_vec((n-n_x+1):n)) ) &
    * (c**3*L_X/(12*L_Y))
  pts_vec(1:n) = Pi_vec                                                                                     ! the solution (pressure-like function) is saved in the array pts_vec so that it will be available at the next time step
  pts_vec(n+1) = t                                                                                          ! the time at the current time step is saved in the array pts_vec so that it will be available at the next time step
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Release memory
  ! --------------------------------------------------------------------------------------------------------
  
  DEALLOCATE(nn0_vec)
  DEALLOCATE(nnN_vec)
  DEALLOCATE(nnE_vec)
  DEALLOCATE(nnS_vec)
  DEALLOCATE(nnW_vec)
  DEALLOCATE(aN_vec)
  DEALLOCATE(aE_vec)
  DEALLOCATE(aS_vec)
  DEALLOCATE(aW_vec)
  DEALLOCATE(coef_vec)
  DEALLOCATE(coefIndex_vec)
  DEALLOCATE(values_vec)
  DEALLOCATE(columns_vec)
  DEALLOCATE(rowIndex_vec)
  DEALLOCATE(md_vec)
  DEALLOCATE(R_vec)
  DEALLOCATE(initial_guess_vec)
  DEALLOCATE(Pi_dof_vec)
  DEALLOCATE(a_lu_vec)
  DEALLOCATE(j_lu_vec)
  DEALLOCATE(j_u_vec)
  DEALLOCATE(ws_vec)
  DEALLOCATE(iws_vec)
  DEALLOCATE(os_x_vec)
  DEALLOCATE(os_y_vec)
  DEALLOCATE(os_vec)
  
  
END SUBROUTINE FVM_ELROD



SUBROUTINE PREPARE_ASSEMBLY_ELROD(n_x, dof_y, dof, n_nz, md_vec, rowIndex_vec, coefIndex_vec, columns_vec)
  
  ! no implicit variables
  IMPLICIT NONE
  
  ! input variables
  INTEGER,INTENT(IN)                                    :: n_x, dof_y, dof, n_nz
  
  ! output variables
  INTEGER,DIMENSION(dof),INTENT(OUT)                    :: md_vec
  INTEGER,DIMENSION(dof+1),INTENT(OUT)                  :: rowIndex_vec
  INTEGER,DIMENSION(n_nz),INTENT(OUT)                    :: coefIndex_vec, columns_vec
  
  ! local variables
  INTEGER                                                :: i, j, k, s1, s2
  
  k = 0
  rowIndex_vec(1) = 1
  ! bottom left corner of the domain:
  k = k+1                                                                                                   ! global DOF number
  s1 = rowIndex_vec(k)                                                                                      ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  s2 = s1+3                                                                                                 ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  rowIndex_vec(k+1) = s2+1                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
  coefIndex_vec(s1:s2) = (/k,k+2*dof,k+4*dof,k+dof/)                                                        ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
  columns_vec(s1:s2) = (/k,k+1,k-1+n_x,k+n_x/)                                                              ! matrix column numbers corresponding to the entries of the k-th matrix row
  md_vec(k) = s1                                                                                            ! array index of values_vec where the main diagonal entry will be stored (this information is important for the oil supply BCs)
  ! bottom boundary of the domain except for corners:
  DO i = 2, (n_x-1)                                                                                         ! loop through circumferential DOF numbers
    k = k+1                                                                                                 ! global DOF number
    s1 = rowIndex_vec(k)                                                                                    ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    s2 = s1+3                                                                                               ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    rowIndex_vec(k+1) = s2+1                                                                                ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
    coefIndex_vec(s1:s2) = (/k+4*dof,k,k+2*dof,k+dof/)                                                      ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
    columns_vec(s1:s2) = (/k-1,k,k+1,k+n_x/)                                                                ! matrix column numbers corresponding to the entries of the k-th matrix row
    md_vec(k) = s1 + 1                                                                                      ! array index of values_vec where the main diagonal entry will be stored (this information is important for the oil supply BCs)
  END DO
  ! bottom right corner of the domain:
  k = k+1                                                                                                   ! global DOF number
  s1 = rowIndex_vec(k)                                                                                      ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  s2 = s1+3                                                                                                 ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  rowIndex_vec(k+1) = s2+1                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
  coefIndex_vec(s1:s2) = (/k+2*dof,k+4*dof,k,k+dof/)                                                        ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
  columns_vec(s1:s2) = (/k+1-n_x,k-1,k,k+n_x/)                                                              ! matrix column numbers corresponding to the entries of the k-th matrix row
  md_vec(k) = s1 + 2                                                                                        ! array index of values_vec where the main diagonal entry will be stored (this information is important for the oil supply BCs)
  DO j = 2, (dof_y-1)                                                                                       ! loop through axial DOF numbers
    ! left boundary of the domain except for corners:
    k = k+1                                                                                                 ! global DOF number
    s1 = rowIndex_vec(k)                                                                                    ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    s2 = s1+4                                                                                               ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    rowIndex_vec(k+1) = s2+1                                                                                ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
    coefIndex_vec(s1:s2) = (/k+3*dof,k,k+2*dof,k+4*dof,k+dof/)                                              ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
    columns_vec(s1:s2) = (/k-n_x,k,k+1,k-1+n_x,k+n_x/)                                                      ! matrix column numbers corresponding to the entries of the k-th matrix row
    md_vec(k) = s1 + 1                                                                                      ! array index of values_vec where the main diagonal entry will be stored (this information is important for the oil supply BCs)
    ! all except boundaries:
    DO i = 2, (n_x-1)                                                                                       ! loop through circumferential DOF numbers
      k = k+1                                                                                               ! global DOF number
      s1 = rowIndex_vec(k)                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
      s2 = s1+4                                                                                             ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
      rowIndex_vec(k+1) = s2+1                                                                              ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
      coefIndex_vec(s1:s2) = (/k+3*dof,k+4*dof,k,k+2*dof,k+dof/)                                            ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
      columns_vec(s1:s2) = (/k-n_x,k-1,k,k+1,k+n_x/)                                                        ! matrix column numbers corresponding to the entries of the k-th matrix row
      md_vec(k) = s1 + 2                                                                                    ! array index of values_vec where the main diagonal entry will be stored (this information is important for the oil supply BCs)
    END DO
    ! right boundary of the domain except for corners:
    k = k+1                                                                                                 ! global DOF number
    s1 = rowIndex_vec(k)                                                                                    ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    s2 = s1+4                                                                                               ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    rowIndex_vec(k+1) = s2+1                                                                                ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
    coefIndex_vec(s1:s2) = (/k+3*dof,k+2*dof,k+4*dof,k,k+dof/)                                              ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
    columns_vec(s1:s2) = (/k-n_x,k+1-n_x,k-1,k,k+n_x/)                                                      ! matrix column numbers corresponding to the entries of the k-th matrix row
    md_vec(k) = s1 + 3                                                                                      ! array index of values_vec where the main diagonal entry will be stored (this information is important for the oil supply BCs)
  END DO
  ! top left corner of the domain:
  k = k+1                                                                                                   ! global DOF number
  s1 = rowIndex_vec(k)                                                                                      ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  s2 = s1+3                                                                                                 ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  rowIndex_vec(k+1) = s2+1                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
  coefIndex_vec(s1:s2) = (/k+3*dof,k,k+2*dof,k+4*dof/)                                                      ! array indices of coef_vec where the entries of the k-th matrix row will be found  (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
  columns_vec(s1:s2) = (/k-n_x,k,k+1,k-1+n_x/)                                                              ! matrix column numbers corresponding to the entries of the k-th matrix row
  md_vec(k) = s1 + 1                                                                                        ! array index of values_vec where the main diagonal entry will be stored (this information is important for the oil supply BCs)
  ! top boundary of the domain except for corners:
  DO i = 2, (n_x-1)                                                                                         ! loop through circumferential DOF numbers
    k = k+1                                                                                                 ! global DOF number
    s1 = rowIndex_vec(k)                                                                                    ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    s2 = s1+3                                                                                               ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    rowIndex_vec(k+1) = s2+1                                                                                ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
    coefIndex_vec(s1:s2) = (/k+3*dof,k+4*dof,k,k+2*dof/)                                                    ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
    columns_vec(s1:s2) = (/k-n_x,k-1,k,k+1/)                                                                ! matrix column numbers corresponding to the entries of the k-th matrix row
    md_vec(k) = s1 + 2                                                                                      ! array index of values_vec where the main diagonal entry will be stored (this information is important for the oil supply BCs)
  END DO
  ! top right corner of the domain:
  k = k+1                                                                                                   ! global DOF number
  s1 = rowIndex_vec(k)                                                                                      ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  s2 = s1+3                                                                                                 ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  rowIndex_vec(k+1) = s2+1                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
  coefIndex_vec(s1:s2) = (/k+3*dof,k+2*dof,k+4*dof,k/)                                                      ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
  columns_vec(s1:s2) = (/k-n_x,k+1-n_x,k-1,k/)                                                              ! matrix column numbers corresponding to the entries of the k-th matrix row
  md_vec(k) = s1 + 3                                                                                        ! array index of values_vec where the main diagonal entry will be stored (this information is important for the oil supply BCs)
  
END SUBROUTINE PREPARE_ASSEMBLY_ELROD



SUBROUTINE PREPARE_ASSEMBLY_GUEMB(n_x, dof_y, dof, n_nz, rowIndex_vec, coefIndex_vec, columns_vec)
  
  ! no implicit variables
  IMPLICIT NONE
  
  ! input variables
  INTEGER,INTENT(IN)                                    :: n_x, dof_y, dof, n_nz
  
  ! output variables
  INTEGER,DIMENSION(dof+1),INTENT(OUT)                  :: rowIndex_vec
  INTEGER,DIMENSION(n_nz),INTENT(OUT)                    :: coefIndex_vec, columns_vec
  
  ! local variables
  INTEGER                                                :: i, j, k, s1, s2
  
  k = 0
  rowIndex_vec(1) = 1
  ! bottom left corner of the domain:
  k = k+1                                                                                                   ! global DOF number
  s1 = rowIndex_vec(k)                                                                                      ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  s2 = s1+3                                                                                                 ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  rowIndex_vec(k+1) = s2+1                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
  coefIndex_vec(s1:s2) = (/k,k+2*dof,k+4*dof,k+dof/)                                                        ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
  columns_vec(s1:s2) = (/k,k+1,k-1+n_x,k+n_x/)                                                              ! matrix column numbers corresponding to the entries of the k-th matrix row
  ! bottom boundary of the domain except for corners:
  DO i = 2, (n_x-1)                                                                                         ! loop through circumferential DOF numbers
    k = k+1                                                                                                 ! global DOF number
    s1 = rowIndex_vec(k)                                                                                    ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    s2 = s1+2                                                                                               ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    rowIndex_vec(k+1) = s2+1                                                                                ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
    coefIndex_vec(s1:s2) = (/k,k+2*dof,k+dof/)                                                              ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
    columns_vec(s1:s2) = (/k,k+1,k+n_x/)                                                                    ! matrix column numbers corresponding to the entries of the k-th matrix row
  END DO
  ! bottom right corner of the domain:
  k = k+1                                                                                                   ! global DOF number
  s1 = rowIndex_vec(k)                                                                                      ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  s2 = s1+1                                                                                                 ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  rowIndex_vec(k+1) = s2+1                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
  coefIndex_vec(s1:s2) = (/k,k+dof/)                                                                        ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
  columns_vec(s1:s2) = (/k,k+n_x/)                                                                          ! matrix column numbers corresponding to the entries of the k-th matrix row
  DO j = 2, (dof_y-1)                                                                                       ! loop through axial DOF numbers
    ! left boundary of the domain except for corners:
    k = k+1                                                                                                 ! global DOF number
    s1 = rowIndex_vec(k)                                                                                    ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    s2 = s1+3                                                                                               ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    rowIndex_vec(k+1) = s2+1                                                                                ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
    coefIndex_vec(s1:s2) = (/k,k+2*dof,k+4*dof,k+dof/)                                                      ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
    columns_vec(s1:s2) = (/k,k+1,k-1+n_x,k+n_x/)                                                            ! matrix column numbers corresponding to the entries of the k-th matrix row
    ! all except boundaries:
    DO i = 2, (n_x-1)                                                                                       ! loop through circumferential DOF numbers
      k = k+1                                                                                               ! global DOF number
      s1 = rowIndex_vec(k)                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
      s2 = s1+2                                                                                             ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
      rowIndex_vec(k+1) = s2+1                                                                              ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
      coefIndex_vec(s1:s2) = (/k,k+2*dof,k+dof/)                                                            ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
      columns_vec(s1:s2) = (/k,k+1,k+n_x/)                                                                  ! matrix column numbers corresponding to the entries of the k-th matrix row
    END DO
    ! right boundary of the domain except for corners:
    k = k+1                                                                                                 ! global DOF number
    s1 = rowIndex_vec(k)                                                                                    ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    s2 = s1+1                                                                                               ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    rowIndex_vec(k+1) = s2+1                                                                                ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
    coefIndex_vec(s1:s2) = (/k,k+dof/)                                                                      ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
    columns_vec(s1:s2) = (/k,k+n_x/)                                                                        ! matrix column numbers corresponding to the entries of the k-th matrix row
  END DO
  ! top left corner of the domain:
  k = k+1                                                                                                   ! global DOF number
  s1 = rowIndex_vec(k)                                                                                      ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  s2 = s1+2                                                                                                 ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  rowIndex_vec(k+1) = s2+1                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
  coefIndex_vec(s1:s2) = (/k,k+2*dof,k+4*dof/)                                                              ! array indices of coef_vec where the entries of the k-th matrix row will be found  (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
  columns_vec(s1:s2) = (/k,k+1,k-1+n_x/)                                                                    ! matrix column numbers corresponding to the entries of the k-th matrix row
  ! top boundary of the domain except for corners:
  DO i = 2, (n_x-1)                                                                                         ! loop through circumferential DOF numbers
    k = k+1                                                                                                 ! global DOF number
    s1 = rowIndex_vec(k)                                                                                    ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    s2 = s1+1                                                                                               ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
    rowIndex_vec(k+1) = s2+1                                                                                ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
    coefIndex_vec(s1:s2) = (/k,k+2*dof/)                                                                    ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
    columns_vec(s1:s2) = (/k,k+1/)                                                                          ! matrix column numbers corresponding to the entries of the k-th matrix row
  END DO
  ! top right corner of the domain:
  k = k+1                                                                                                   ! global DOF number
  s1 = rowIndex_vec(k)                                                                                      ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  s2 = s1                                                                                                   ! last array index of columns_vec and coefIndex_vec (and later values_vec) where information about the k-th matrix row will be stored
  rowIndex_vec(k+1) = s2+1                                                                                  ! first array index of columns_vec and coefIndex_vec (and later values_vec) where information about the (k+1)-th matrix row will be stored
  coefIndex_vec(s1:s2) = (/k/)                                                                              ! array indices of coef_vec where the entries of the k-th matrix row will be found (coef_vec will be defined later; it's an array containing all coefficients on the left-hand side of the discretized Reynolds equation); the indices k, k+dof, k+2*dof, k+3*dof, and k+4*dof point to the coefficients corresponding to the DOFs 0, N, E, S, and W in this order
  columns_vec(s1:s2) = (/k/)                                                                                ! matrix column numbers corresponding to the entries of the k-th matrix row
  
END SUBROUTINE PREPARE_ASSEMBLY_GUEMB
