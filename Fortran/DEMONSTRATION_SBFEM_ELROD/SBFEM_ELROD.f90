! SBFEM_ELROD.f90 
!
! FUNCTIONS:
! SBFEM_ELROD
!
! **********************************************************************************************************
!
! PROGRAM: SBFEM_ELROD
!
! PURPOSE: see below
!
! AUTHOR, AFFILIATION, DATE: Simon Pfeil, OvGU Magdeburg (Germany), 30.06.2025
!
! **********************************************************************************************************
!
! COMMENTS
!
! - This algorithm solves the Reynolds equation for hydroynamic journal bearings.
!
! - The computation is based on the semi-analytical "scaled boundary finite element method" (SBFEM).
!
! - A transient, nonlinear cavitation model is used (optionally, the simple Guembel approach can be used 
!   instead).
!
! - Shaft tilting is not considered.
!
! - The oil supply groove is modeled with simplified geometry: its axial length is assumed to be equal to
!   the bearing length.
!
! - This program is designed to be incorporated into time integration schemes for rotordynamic (or 
!   multibody) simulations; however, a quasi-static solution is also possible.
!
! - The input and output variables are explained below.
!
! - See the attached PDF for clarification of the coordinate systems and node numbers. x (or the 
!   nondimensionalized X) and y (or the nondimensionalized Y) are the coordinates of the lubrication gap 
!   and are used by the Reynolds equation. This coordinate system, as well as the computational grid, are 
!   fixed in the reference frame of the shell. The transformations between the reference frame of the 
!   shell and the inertial system are performed internally by this program. The kinematic variables of the 
!   shell and the shaft, when handed to this program, should be formulated from the perspective of the 
!   inertial system. The hydrodynamic forces computed by this program are automatically transformed into 
!   the inertial system before output. In contrast, those input and output variables that describe one- or 
!   two-dimensional fields (e.g., viscosity distribution, additional contour, and pressure-like function) 
!   never use the inertial system, as they are expressed by arrays whoose entries represent the values at 
!   the nodes (which are always fixed at the shell). The circumferential position of an oil supply groove 
!   X_os is also always expressed in the reference frame of the shell; otherwise, a rotating shell would 
!   require us to keep updating this input variable.
!
! - More thorough discussions of the computational method, of the cavitation model, of the boundary 
!   conditions (BCs), and of further assumptions will be given in the thesis "Simulating Hydrodynamic 
!   Bearings with the Scaled Boundary Finite Element Method" by Simon Pfeil, but this thesis isn't 
!   publically available yet
!
! - This program relies on BLAS and LAPACK routines, which are available, for example, through the Intel
!   oneAPI Math Kernel Library (oneMKL); this library is included in the Intel oneAPI Base Toolkit.
!   The code can be compiled via Visual Studio.
!
! - In the cavitation zone, the sparsity pattern of the coefficient matrix E2 depends on the flow direction
!   (due to the upwind scheme), which is determined by the sign of the surface velocity u. Parts of the
!   algorithm exploit this sparsity pattern for the sake of efficiency and are quite complicated, but they
!   would become even more complicated if we wanted to qualify them for both positive and negative values 
!   of u, so I restricted them to the positive case. Nonetheless, this program accepts negative values 
!   of u and produces the correct results, as it internally rotates the coordinate system of the Reynolds 
!   equation by 180 degrees (flipping the direction of the circumferential coordinate as well as the axial
!   coordinate) so that u becomes positive again. This transformation is applied to the corresponding
!   input variables (and later reversed for the output variables) internally. We don't need to worry
!   about this when applying this program, but we need to keep it in mind if we want to analyze this 
!   program or check intermediate results.
!
! - At some point in the algorithm, the node numbers are shifted depending on the position of the oil
!   supply groove, which ensures a tridiagonal E2-matrix. This affects the order in which the nodal values
!   of all discretized functions are stored. Since this operation is performed (and, at the end, reversed) 
!   internally, we don't need to worry about it when using this program. If, however, we want to analyze 
!   the algorithm step by step or check intermediate results, we need to keep this in mind.
!
! **********************************************************************************************************
!
! INPUT VARIABLES
!
! - grooves = number of oil supply grooves; the grooves are assumed to be distributed equidistantly along 
!   the circumference; at least 1 groove must be present under Elrod conditions (guembel = 0), and at least 
!   0 grooves must be present under Guembel conditions (guembel = 1) [-]
!
! - n_x = number of nodes for circumferential equidistant grid; the smallest node number 1 is located at
!   the angular position X=0 and the largest node number n_x at X=2*pi*(1-1/n_x); the node at X=2*pi is
!   not counted because the periodic BC merges this node with node 1 (see attached PDF file) [-]
!
! - iter_max = maximum allowed number of iterations for the cavitation model (usually, less than 10 
!   iterations are enough, but in some cases, more may be necessary; I usually choose iter_max = n_x) [-]
!
! - quasistatic = a flag that specifies whether a quasi-static solution is desired (1 = yes, 0 = no); 
!   if yes, then all velocities and angular velocities except for omega_shaft are set to zero and the film 
!   fraction is assumed to be constant over time [-]
!
! - guembel = determines whether Guembel conditions should be assumed, bypassing the Elrod cavitation 
!   model; guembel = 0 --> use Elrod, guembel = 1 --> use Guembel [-]
! 
! - n_y = number of points in the axial direction for evaluation of the solution (pressure-like function)
!   on a two-dimensional grid (optional postprocessing step for creation of a surface plot); this step is
!   skipped if n_y = 1 [-]
!
! - d_b = bearing diameter [m]
!
! - l_b = bearing length [m]
!
! - c = radial clearance (shell radius minus shaft radius) [m]
!
! - X_os = angular circumferential position of the center of one of the oil supply grooves (must not be 
!   smaller than 0 or greater than 2*pi, see attached PDF); since the coordinate X is fixed at the shell, 
!   a rotation of the shell does not change the value of X_os [rad]
!
! - L_X_os = angular circumferential side length of the oil supply groove (see attached PDF file) [rad]
!
! - p_os (if positive or zero) = pressure prescribed in the oil supply groove (zero corresponds to
!   atmospheric pressure) [Pa]
!
! - p_os (if negative) = film fraction prescribed in the oil supply groove minus 1 (for example, a film 
!   fraction of 0.7 is prescribed by setting p_os = 0.7-1 = -0.3) [-]
!
! - t = physical time at the current time step [s]
!
! - angle_shell = rotation angle of the shell from the perspective of the inertial system (see attached 
!   PDF); if the shell rotates, knowing this angle is crucial in order for this program to perform the 
!   transformations between the inertial system and the reference frame of the shell [rad]
!
! - omega_shell, omega_shaft = angular velocities of the shell and of the shaft, respectively, around the 
!   longitudinal axis, measured by a non-rotating observer [rad/s]
!
! - dis_h_shell, dis_h_shaft = horizontal displacements of the shell and the shaft, respectively, from the 
!   perspective of the inertial system (see attached PDF) [m]
!
! - dis_v_shell, dis_v_shaft = vertical displacements of the shell and the shaft, respectively, from the 
!   perspective of the inertial system (see attached PDF) [m]
!
! - vel_h_shell, vel_h_shaft = horizontal velocities of the shell and the shaft, respectively, from the 
!   perspective of the inertial system (see attached PDF) [m/s]
!
! - vel_v_shell, vel_v_shaft = vertical velocities of the shell and the shaft, respectively, from the 
!   perspective of the inertial system (see attached PDF) [m/s]
!
! - ac_vec = column vector of length n_x containing an additional contour along the circumferential (but 
!   not axial) direction in case the shell is not cylindrical; the entries of this vector are simply added 
!   to the nodal gap widths; since the nodes are fixed at the shell, a rotation of the shell does not
!   require this additional contour to be updated; note that an additional contour of the shaft is not 
!   allowed in the current version of this program [m]
!
! - mu_vec = column vector of length n_x containing the dynamic oil viscosities at the nodes, from node 1 
!   to node n_x; since the nodes are fixed at the possibly rotating shell, the viscosity distribution is 
!   always expressed in the shell's reference frame [Pa*s]
!
! - pts_vec = array of length n_x+1 where the first n_x entries are the nodal pressure-like functions
!   of the previous time step (averaged in the axial direction) and the entry n_x+1 states the time of the
!   previous time step; this variable is used for transferring data between time steps, as required by the 
!   transient cavitation model; note that pts_vec is also an output variable; simply take the pts_vec that 
!   was generated as output at the previous valid time step and provide it as input variable at the 
!   current time step [-]/[s]
!
! - tay = determines whether the eigenvalue problem is solved by a Taylor approximation (0 = no, 1 = yes); 
!   the provided database of Taylor coefficients assumes that the additional contour ac_vec is zero, that 
!   the viscosity distribution described by mu_vec is spatially constant, and that n_x = 100, grooves = 0, 
!   and guembel = 1, so we must choose our parameters accordingly when using tay = 1; the script
!   demonstration_sbfem_elrod.m explains how to use this database [-]
!
! - n_tay = order of the Taylor approximation (i.e. polynomial degree of the Taylor series); only used in 
!   case tay = 1; the provided database of Taylor coefficients assumes that n_tay = 2 [-]
!
! - n_Ld = length of the vector containing the Taylor coefficients for the eigenvalues [-]
!
! - red = number of considered modes of the pressure field due to modal reduction; only used in case 
!   tay = 1; if we use the provided database of Taylor coefficients, this database also determines what 
!   number of modes to choose for a given eccentricity (see demonstration_sbfem_elrod.m) [-]
! 
! - eps_constr = relative eccentricity where the Taylor series that will be employed was constructed; only 
!   used if tay = 1 [-] 
!
! - eps_max = maximum relative eccentricity at which Taylor approximations for solving the eigenvalue 
!   problem are available; if the current relative eccentricity exceeds epsilon_max, the eigenvalue problem 
!   is solved by an eigensolver even if tay = 1 [-]
!
! - Ld_vec = vector containing the eigenvalue derivatives (coefficients of the Taylor series to be used); 
!   only used if tay = 1 (otherwise, simply define Ld_vec as an empty variable); the script 
!   demonstration_sbfem_elrod.m demonstrates how to obtain Ld_vec from the database provided together with 
!   this program [-]
!
! - Vd_mat = matrix containing the eigenvector derivatives (coefficients of the Taylor series to be used); 
!   only used if tay = 1 (otherwise, simply define Vd_mat as an empty variable); the script 
!   demonstration_sbfem_elrod.m demonstrates how to obtain Vd_mat from the database provided together with 
!   this program [-]
!
! **********************************************************************************************************
!
! OUTPUT VARIABLES
!
! - convergent = a flag that indicates whether or not the solution has converged (1 = yes, 0 = no) [-]
!
! - iter = number of iterations [-]
!
! - M_fr = oil friction moment acting on the shell (see attached PDF); the same moment acts on the shaft
!   in the opposite direction [Nm]
!
! - V_oil = total oil volume in the bearing (equal to the volume of the lubrication gap minus the gas
!   volume in the cavitation zone) [m^3]
!
! - V_dot_bb = oil volume flow through the bearing boundaries (a negative sign indicates a flow direction
!   out of the bearing) [m^3/s]
!
! - p_ref = reference pressure due to nondimensionalization [Pa]
!
! - F_h = hydrodynamic force acting on the shell in the horizontal direction of the inertial system
!   (see attached PDF); the same force acts on the shaft in the opposite direction [N]
!
! - F_v = hydrodynamic force acting on the shell in the vertical direction of the inertial system (see 
!   attached PDF); the same force acts on the shaft in the opposite direction [N]
!
! - g_vec = column vector of length n_x containing the switch functions g at the circumferential 
!   positions of the nodes; g = 1 is the pressure zone and g = 0 is the cavitation zone; when analyzing 
!   g_vec for a system where the shell rotates, we need to keep in mind that the nodes are fixed at the 
!   shell (and so are the entries of g_vec) [-]
!
! - theta_vec = column vector of length n_x containing the film fractions theta at the circumferential 
!   positions of the nodes; theta = 1 indicates a fully oil-filled gap and theta < 1 indicates a cavitated 
!   fluid film; when analyzing theta_vec for a system where the shell rotates, we need to keep in mind 
!   that the nodes are fixed at the shell (and so are the entries of theta_vec); under Guembel conditions
!   (guembel = 1), theta_vec is quite meaningless [-]
!
! - Pi_mat = matrix of size n_x*n_y stating the two-dimensional distribution of the computed pressure-like 
!   function Pi; from this, the pressure distribution can be derived (p = p_ref*Pi*g), where g is the
!   switch function (see above) and p_ref is the reference pressure which was used for the 
!   nondimensionalization (see below); Pi_mat is only computed if n_y > 1; each row of Pi_mat represents a 
!   circumferential position where a node is located (if the shell rotates, we need to keep in mind that
!   the nodes are fixed at the shell), while each column represents an axial position [-]
!
! - pts_vec = array of length n_x+1 where the first n_x entries are the axially averaged nodal
!   pressure-like functions and the entry n_x+1 contains the current physical time [-]/[s]
!
! **********************************************************************************************************


  SUBROUTINE SBFEM_ELROD(grooves, n_x, iter_max, quasistatic, guembel, n_y, d_b, l_b, c, X_os, &
    L_X_os, p_os, t, angle_shell, omega_shell, dis_h_shell, dis_v_shell, vel_h_shell, vel_v_shell, &
    omega_shaft, dis_h_shaft, dis_v_shaft, vel_h_shaft, vel_v_shaft, ac_vec, mu_vec, tay, n_tay, &
    n_Ld, red, eps_constr, eps_max, Ld_vec, Vd_mat, convergent, iter, M_fr, V_oil, V_dot_bb, p_ref, &
    F_h, F_v, g_vec, theta_vec, Pi_mat, pts_vec)
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Variable types
  ! --------------------------------------------------------------------------------------------------------
  
  ! no implicit variables
  IMPLICIT NONE
  
  ! input variables
  INTEGER,INTENT(IN)                            :: n_x, n_y, grooves, iter_max, quasistatic
  INTEGER,INTENT(IN)                            :: guembel, tay, n_tay, n_Ld, red
  REAL(KIND=8),INTENT(IN)                       :: d_b, l_b, c, X_os, L_X_os, p_os, t, angle_shell
  REAL(KIND=8),INTENT(IN)                       :: omega_shell, dis_h_shell, dis_v_shell, vel_h_shell
  REAL(KIND=8),INTENT(IN)                       :: vel_v_shell, omega_shaft, dis_h_shaft, dis_v_shaft
  REAL(KIND=8),INTENT(IN)                       :: vel_h_shaft, vel_v_shaft, eps_constr, eps_max
  REAL(KIND=8),DIMENSION(n_x),INTENT(IN)        :: ac_vec, mu_vec
  REAL(KIND=8),DIMENSION(n_Ld),INTENT(IN)       :: Ld_vec
  REAL(KIND=8),DIMENSION(n_x,n_Ld),INTENT(IN)   :: Vd_mat
  
  ! output variables
  INTEGER,INTENT(OUT)                           :: convergent, iter
  INTEGER,DIMENSION(n_x),INTENT(OUT)            :: g_vec
  REAL(KIND=8),INTENT(OUT)                      :: F_h, F_v, M_fr, V_oil, V_dot_bb, p_ref
  REAL(KIND=8),DIMENSION(n_x),INTENT(OUT)       :: theta_vec
  REAL(KIND=8),DIMENSION(n_x,n_y),INTENT(OUT)   :: Pi_mat
  
  ! input/output variables
  REAL(KIND=8),DIMENSION(n_x+1),INTENT(INOUT)   :: pts_vec
  
  ! local
  INTEGER                                       :: i, j, k, l, n_os, dof, start_os1, n_p, n_c, liwork
  INTEGER                                       :: lwork, info, n1, n2, indx1, indx2, cav1, cav2, n_ab
  INTEGER                                       :: n_pz, n_cpz, incr, g_os, n_mod, di_floor
  INTEGER,DIMENSION(grooves)                    :: start_os_vec, end_os_vec
  INTEGER,DIMENSION(n_x)                        :: index_0_vec, index_E_vec, index_W_vec, shift_vec
  INTEGER,DIMENSION(n_x)                        :: allnodes_vec, g_pts_vec, g_old_vec, g_old_old_vec
  INTEGER,DIMENSION(n_x)                        :: states_vec
  INTEGER,DIMENSION(:),ALLOCATABLE              :: nodes_dof_vec, nodes_os_vec, nodes_p_vec, nodes_c_vec 
  INTEGER,DIMENSION(:),ALLOCATABLE              :: start_p_vec, end_p_vec, indx_sp_vec, indx_ep_vec
  INTEGER,DIMENSION(:),ALLOCATABLE              :: iwork_vec
  REAL(KIND=8)                                  :: pi, L_X, L_Y, epsil, epsil_dot, sgn_u, alpha, beta
  REAL(KIND=8)                                  :: mu_ref, Delta_T, t_pts, entry_a, entry_b, entry_c
  REAL(KIND=8)                                  :: Pi_os, u_pos, X_os_loc, F_1, F_2, u, q, X_att, q_dot
  REAL(KIND=8)                                  :: X_att_dot, omega, dis_h, dis_v, vel_h, vel_v
  REAL(KIND=8)                                  :: X_att_trafo, Delta_eps, fac, di, fl, fr
  REAL(KIND=8)                                  :: zero, q_squared
  REAL(KIND=8),DIMENSION(grooves)               :: X_os_vec
  REAL(KIND=8),DIMENSION(n_x)                   :: X_vec, Pi_bar_pts_vec, H_vec, dHdT_vec, H3_over_mu_vec
  REAL(KIND=8),DIMENSION(n_x)                   :: b_E_vec, b_W_vec, R_vec, diag_E2_vec, Ldiag_E2_vec
  REAL(KIND=8),DIMENSION(n_x)                   :: Rdiag_E2_vec, Pi_con_vec, Pi_bar_vec, dpdy_bb_vec
  REAL(KIND=8),DIMENSION(n_x)                   :: diag_E0_vec, mu_rel_vec, ac_rel_vec
  REAL(KIND=8),DIMENSION(n_x,n_x)               :: E2_mat
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE         :: diag_B_vec, Ldiag_B_vec, trafo_vec, work_vec, lambda_vec
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE         :: C_vec, C_i_vec, xi_vec, Vec1, Vec2, Vec3, lambda_i_vec
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE         :: trafo_i_vec, diag_B_i_vec, Ldiag_B_i_vec
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE       :: Phi_mat, Phi_i_mat, Mat1, Mat2, cosh_mat, Pip_mat
  
  ! external subroutines
  EXTERNAL                                      :: DSTEDC, DGEMM, DGEMV, DSYEVD
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Some settings for the DGEMM and DGEMV routines
  ! --------------------------------------------------------------------------------------------------------
  
  incr = 1
  alpha = 1.0d0
  beta = 0.0d0
  
  
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
  IF ( quasistatic .EQ. 0 ) THEN
    omega = omega_shaft-omega_shell                                                                         ! rotational velocity of the shaft relative to the shell
    vel_h = vel_h_shaft-vel_h_shell                                                                         ! horizontal velocity of the shaft relative to the shell, still in the reference frame of the inertial system
    vel_v = vel_v_shaft-vel_v_shell                                                                         ! vertical velocity of the shaft relative to the shell, still in the reference frame of the inertial system
    q_dot = (vel_v*dis_v+vel_h*dis_h)/q                                                                     ! rate of change of absolute eccentricity
    X_att_dot = (vel_v*dis_h-vel_h*dis_v)/q**2                                                              ! rate of change of attitude angle in the reference frame of the inertial system
    X_att_dot = X_att_dot-omega_shell                                                                       ! rate of change of attitude angle in the shell-fixed reference frame in which the Reynolds equation is solved
  ELSE
    omega = omega_shaft                                                                                     ! quasistatic case: rotational velocity of shell is assumed to be zero
    q_dot = 0                                                                                               ! quasistatic case: rate of change of eccentricity is assumed to be zero
    X_att_dot = 0                                                                                           ! quasistatic case: rate of change of attitude angle is assumed to be zero
  END IF
  u = omega*(d_b/2.0d0)                                                                                     ! circumferential surface velocity of the shaft in the shell-fixed reference frame in which the Reynolds equation is solved (the surface velocity of the shell is zero in this reference frame) [m/s]
  IF ( tay .EQ. 1 ) THEN                                                                                    ! if a Taylor-series solution of the eigenvalue problem is desired, we must shift the computational grid (and the coordinate system) so that node 1 is located at the position of the minimum gap; after this shift, the attitude angle is zero
    X_att_trafo = X_att                                                                                     ! since the attitude angle will be set to zero in the line below, we must now save the actual attitude angle to allow transforming the solution back at the end
    X_att = 0                                                                                               ! the attitude angle must be set to zero for now, but the correct attitude angle will be considered later via a simple transformation
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Nondimensionalization 1
  ! --------------------------------------------------------------------------------------------------------
  
  epsil = q/c                                                                                               ! relative eccentricity [-]
  epsil_dot = q_dot/c                                                                                       ! rate of change of relative eccentricity [1/s]
  mu_ref = SUM(mu_vec)/n_x                                                                                  ! reference viscosity for nondimensionalization [Pa*s]
  mu_rel_vec = mu_vec/mu_ref                                                                                ! relative (i.e. nondimensionalized) viscosities [-]
  ac_rel_vec = ac_vec/c                                                                                     ! relative (i.e. nondimensionalized) additional contour [-]
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Rotate coordinate system of Reynolds equation by 180 degrees if u is negative
  ! --------------------------------------------------------------------------------------------------------
  
  IF ( u .LT. 0.0d0 ) THEN                                                                                  ! "IF the circumferential surface velocity is negative, then ..."
    u_pos = -u                                                                                              ! create copy of surface velocity and change sign
    X_os_loc = -X_os                                                                                        ! create local copy of position of oil supply groove and change sign
    X_att = -X_att                                                                                          ! change sign of attitude angle
    X_att_dot = -X_att_dot                                                                                  ! change sign of rate of change of attitude angle
    X_att_trafo = -X_att_trafo                                                                              ! change sign of attitude angle
    index_0_vec(1) = 1                                                                                      ! this array will describe how the rotation of the coordinate system changes the node numbers
    index_0_vec(2:n_x) = (/(i,i=n_x,2,-1)/)                                                                 ! this array will describe how the rotation of the coordinate system changes the node numbers
    mu_rel_vec = mu_rel_vec(index_0_vec)                                                                    ! adjust node numbers (i.e. order of the entries) in the array of nodal viscosities
    ac_rel_vec = ac_rel_vec(index_0_vec)                                                                    ! adjust node numbers (i.e. order of the entries) in the array representing the additional contour
    pts_vec(1:n_x) = pts_vec(index_0_vec)                                                                   ! adjust node numbers (i.e. order of the entries) in the array of nodal solutions of the previous time step
  ELSE                                                                                                      ! otherwise (if u is positive), nothing needs to be changed, but ...
    u_pos = u                                                                                               ! no signs need to be changed, but the copies of these variables still need to be defined for consistency ...
    X_os_loc = X_os                                                                                         ! ...
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Nondimensionalization 2
  ! --------------------------------------------------------------------------------------------------------
  
  sgn_u = SIGN(1.0d0,u_pos)                                                                                 ! sign of the tangential surface velocity; actually, the variable sgn_u is obsolete (a negative sign doesn't occur anyway) in the current version of the algorithm (due to the rotation of the coordinate system, see above), but I left it like this [-]
  p_ref = ABS(u_pos)*mu_ref*(d_b/2)/(2*c**2)                                                                ! reference pressure for nondimensionalization; actually, the "ABS" is obsolete (a negative sign doesn't occur anyway) in the current version of the algorithm (due to the rotation of the coordinate system, see above), but I left it like this [-]
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
  ! Analyze grid
  ! --------------------------------------------------------------------------------------------------------
  
  pi = 3.14159265359d0                                                                                      ! define pi
  L_X = 2*pi/n_x                                                                                            ! dimensionless/angular circumferential length of a sector [rad]
  L_Y = l_b/d_b                                                                                             ! dimensionless axial length of one half of the bearing [-]
  index_0_vec = (/(i,i=1,n_x,1)/)                                                                           ! list of all node numbers from 1 to n_x in ascending order [-]
  index_E_vec(1:(n_x-1)) = index_0_vec(2:n_x)                                                               ! number of the neighboring node E (right neighbor) for every node [-]
  index_E_vec(n_x) = 1                                                                                      ! consider periodicity of domain
  index_W_vec(2:n_x) = index_0_vec(1:(n_x-1))                                                               ! number of the neighboring node W (left neighbor) for every node [-]
  index_W_vec(1) = n_x                                                                                      ! consider periodicity of domain
  X_vec = (index_0_vec-1)*L_X                                                                               ! X-values (circumferential coordinate values) describing the positions of all nodes [rad]
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Analyze locations of oil supply grooves
  ! --------------------------------------------------------------------------------------------------------
  
  X_os_vec = (/(i,i=0,grooves-1,1)/)*(2*pi/grooves) + X_os_loc                                              ! X-values (circumferential coordinate values) describing the positions of the oil supply grooves [rad]
  start_os_vec = ANINT((X_os_vec-L_X_os/2)/L_X+1)                                                           ! node numbers where the oil supply grooves begin
  end_os_vec = ANINT((X_os_vec+L_X_os/2)/L_X+1)                                                             ! node numbers where the oil supply grooves end
  n_os = SUM(end_os_vec-start_os_vec+1)                                                                     ! number of nodes with oil supply BCs
  dof = n_x-n_os                                                                                            ! number of degrees of freedom (DOFs)
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Shift node numbers to ensure tridiagonal E2-matrix
  ! --------------------------------------------------------------------------------------------------------
  
  IF ( grooves .GE. 1 ) THEN
    start_os1 = start_os_vec(1)                                                                             ! node number where the first oil supply groove begins
  ELSE
    start_os1 = 1;                                                                                          ! without oil supply BCs, there is no use in shifting the node numbers (a tridiagonal matrix is not achievable); setting start_os1 to 1 prevents this shift
  END IF
  start_os_vec = start_os_vec-start_os1+1                                                                   ! the node numbers are shifted so that the first oil supply groove begins at node 1
  end_os_vec = end_os_vec-start_os1+1                                                                       ! the node numbers are shifted so that the first oil supply groove begins at node 1
  start_os1 = MODULO(start_os1-1,n_x)+1                                                                     ! node number where the first oil supply groove began before shifting the node numbers (MODULO because only values from 1 to n_x are allowed this time)
  shift_vec(1:(n_x-start_os1+1)) = (/(i,i=start_os1,n_x,1)/)                                                ! this array describes how every node number is changed
  shift_vec((n_x-start_os1+2):n_x) = (/(i,i=1,start_os1-1,1)/)                                              ! this array describes how every node number is changed
  X_vec = X_vec(shift_vec)                                                                                  ! manipulate array of circumferential nodal positions according to shift of node numbers
  pts_vec(1:n_x) = pts_vec(shift_vec)                                                                       ! manipulate array of nodal solutions of previous time step according to shift of node numbers
  mu_rel_vec = mu_rel_vec(shift_vec)                                                                        ! manipulate array of nondimensionalized nodal viscosities according to shift of node numbers
  ac_rel_vec = ac_rel_vec(shift_vec)                                                                        ! manipulate array of nondimensionalized additional contour according to shift of node numbers
  
  
  ! ------------------------------------------------------------------------------------------------------
  ! Allocate arrays with size related to dof
  ! ------------------------------------------------------------------------------------------------------
  
  ALLOCATE(nodes_dof_vec(dof))                                                                              ! this array will be used for storing the node numbers corresponding to nodes with DOFs (i.e. nodes without oil supply BCs)
  ALLOCATE(nodes_os_vec(n_os))                                                                              ! this array will be used for storing the node numbers corresponding to nodes with oil supply BCs
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Save node numbers with and without oil supply BCs separately
  ! --------------------------------------------------------------------------------------------------------
  
  allnodes_vec = index_0_vec                                                                                ! allnodes_vec now contains a list of all node numbers
  DO i = 1, grooves                                                                                         ! loop through all grooves (if more than one groove is present)
    allnodes_vec(start_os_vec(i):end_os_vec(i)) = 0                                                         ! the entries of allnodes_vec corresponding to nodes with oil supply BCs are set to zero
  END DO  
  j = 1
  k = 1
  DO i = 1, n_x                                                                                             ! loop through all node numbers
    IF (allnodes_vec(i) .NE. 0) THEN
      nodes_dof_vec(j) = i                                                                                  ! if the current node has a DOF (has no oil supply BC), its number is saved in nodes_dof_vec
      j = j+1
    ELSE
      nodes_os_vec(k) = i                                                                                   ! if the current node has an oil supply BC, its number is saved in nodes_os_vec
      k = k+1
    END IF
  END DO
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Analyze data from previous time step
  ! --------------------------------------------------------------------------------------------------------
  
  Pi_bar_pts_vec = pts_vec(1:n_x)                                                                           ! extract the previous solution (pressure-like function of the previous time step) from pts_vec
  g_pts_vec = 1                                                                                             ! initialize the switch function from the previous time step as one, but ...
  IF ( guembel .EQ. 0 ) THEN                                                                                ! ... if Elrod cavitation is assumed, then ...
    WHERE(Pi_bar_pts_vec .LT. 0.0d0) g_pts_vec = 0                                                          ! ... set this switch function to 1 at all nodes that, according to the previous solution, were part of the pressure zone
  END IF
  g_vec = g_pts_vec                                                                                         ! use the switch function of the previous time step as initial guess for the current switch function
  g_vec(nodes_os_vec) = g_os                                                                                ! ensure that the switch functions in the oil supply groove are consistent with the prescribed oil supply pressure or oil supply film fraction
  g_old_vec = 0                                                                                             ! switch function of previous iteration (does not exist yet)
  IF ( quasistatic .EQ. 1 ) THEN                                                                            ! IF a quasistatic simulation is desired, then ...
    Delta_T = ABS(u_pos)/(2*(d_b/2))*1.0d16                                                                 ! ... set the dimensionless time increment to an extremely large value, suppressing transient effects
  ELSE                                                                                                      ! IF a transient simulation is desired, then ...
    t_pts = pts_vec(n_x+1)                                                                                  ! ... extract the time of the previous time step from pts_vec
    Delta_T = ABS(u_pos)/(2*(d_b/2))*(t-t_pts)                                                              ! ... compute the dimensionless time increment
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Analyze gap function
  ! --------------------------------------------------------------------------------------------------------
  
  H_vec = 1-epsil*COS(X_vec-X_att)+ac_rel_vec                                                               ! nondimensionalized nodal gap functions [-]
  dHdT_vec = (d_b/ABS(u_pos))*(-epsil_dot*COS(X_vec-X_att)-epsil*X_att_dot*SIN(X_vec-X_att))                ! nondimensionalized time derivatives of the nondimensionalized nodal gap functions [-]
  H3_over_mu_vec = (H_vec**3/mu_rel_vec)                                                                    ! cubed nondimensionalized gap function divided by nondimensionalized oil viscosity [-]
  b_E_vec = (H3_over_mu_vec+H3_over_mu_vec(index_E_vec))/(24*L_X)                                           ! coefficient b_E for describing the interaction between a node 0 and its neighbor E, for all nodes [-]
  b_W_vec = (H3_over_mu_vec+H3_over_mu_vec(index_W_vec))/(24*L_X)                                           ! coefficient b_W for describing the interaction between a node 0 and its neighbor W, for all nodes [-]
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Fixed-point iteration
  ! --------------------------------------------------------------------------------------------------------
  
  convergent = 0                                                                                            ! this variable will be set to 1 once the solution has converged
  
  DO iter = 1, iter_max                                                                                     ! fixed-point iteration for linearization of the nonlinear boundary value problem (BVP); iter counts the iterations
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Assembly and oil supply BCs
    ! ------------------------------------------------------------------------------------------------------
    
    diag_E0_vec = (L_X/(12*L_Y**2))*H3_over_mu_vec*g_vec                                                    ! diagonal of the diagonal matrix E0pp
    diag_E2_vec = (b_E_vec+b_W_vec)*g_vec + H_vec*(1-g_vec) + L_X*((dHdT_vec + H_vec/Delta_T)*(1-g_vec))    ! main diagonal of the coefficient matrix E2
    Ldiag_E2_vec = -b_W_vec*g_vec(index_W_vec)+(-sgn_u-1)*0.5d0*(H_vec(index_W_vec)*(1-g_vec(index_W_vec))) ! left/lower diagonal of the coefficient matrix E2
    Rdiag_E2_vec = -b_E_vec*g_vec(index_E_vec)+(sgn_u-1)*0.5d0*(H_vec(index_E_vec)*(1-g_vec(index_E_vec)))  ! right/upper diagonal of the coefficient matrix E2
    E2_mat = 0.0d0                                                                                          ! initialize the entries of E2 as zero
    DO i = 1, n_x                                                                                           ! loop through the nodes and fill ...
      E2_mat(i,index_W_vec(i)) = Ldiag_E2_vec(i)                                                            ! ... the left/lower diagonal of E2
      E2_mat(i,i) = diag_E2_vec(i)                                                                          ! ... the main diagonal of E2
      E2_mat(i,index_E_vec(i)) = Rdiag_E2_vec(i)                                                            ! ... the right/upper diagonal of E2
    END DO                                                                                                  ! note that E2 is originally a cyclic tridiagonal matrix but is rendered truly tridiagonal by the oil supply BCs
    R_vec = (sgn_u*0.5d0)*(H_vec(index_E_vec)-H_vec(index_W_vec)) + L_X*dHdT_vec - &                        ! right-hand side (RHS) vector of the SBFEM equation
      (L_X/Delta_T)*(H_vec*(1-g_pts_vec)*Pi_bar_pts_vec)
    R_vec = R_vec + SUM(E2_mat(:,nodes_os_vec),2)*Pi_os                                                     ! oil supply BCs: consideration of the influence of the oil supply pressure on the neighboring DOFs via the RHS vector (i.e., shift BCs to the right-hand side)
    DO i = 1, n_os                                                                                          ! oil supply BCs: loop through nodes with oil supply BCs and ...
      j = nodes_os_vec(i)                                                                                   ! ... determine current node number, ...
      E2_mat(j,:) = 0.0d0                                                                                   ! ... set matrix row of E2 to zero (elimination technique, but we leave the zeros in the matrix because we don't want to reorganize or even reallocate the whole thing)
      E2_mat(:,j) = 0.0d0                                                                                   ! ... set matrix column of E2 to zero (elimination technique, but we leave the zeros in the matrix because we don't want to reorganize or even reallocate the whole thing)
      diag_E0_vec(j) = 0.0d0                                                                                ! ... set diagonal entry of diagonal matrix E0 to zero (elimination technique, but we leave the zeros in the matrix/vector because we don't want to reorganize or even reallocate the whole thing)
      R_vec(j) = 0.0d0                                                                                      ! ... set entry of R to zero (elimination technique, but we leave the zeros in the vector because we don't want to reorganize or even reallocate the whole thing)
    END DO
    n_p = SUM(g_vec(nodes_dof_vec))                                                                         ! number of DOFs in the pressure zone
    n_c = dof-n_p                                                                                           ! number of DOFs in the cavitation zone
    states_vec = g_vec                                                                                      ! for every node, this array will contain the information as to whether the node is part of the oil supply (-1), is a DOF in the cavitation zone (0), or is a DOF in the pressure zone (1); initialize as equal to switch function
    states_vec(nodes_os_vec) = -1                                                                           ! indicate nodes with oil supply BCs by a value of -1
    IF ( ( tay .EQ. 0 ) .OR. ( epsil > eps_max ) ) THEN                                                     ! if no Taylor approximations are used (i.e., if the eigenvalue problem will be solved by an eigensolver), we will not perform a modal reduction, so ...
        n_mod = n_p                                                                                         ! ... the number of modes is equal to the number of nodes in the pressure regime
    ELSE                                                                                                    ! otherwise (i.e., if Taylor approximations are used), the computational effort is so small that the number of considered modes actually makes a difference, so ...
        n_mod = red                                                                                         ! ... we will perform a modal reduction (thus, we use the number of modes prescribed by the variable red)
    END IF
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Allocate arrays whoose size depends on cavitation states and initialize some arrays
    ! ------------------------------------------------------------------------------------------------------
    
    ALLOCATE(nodes_p_vec(n_p))                                                                              ! array for storing the node numbers corresponding to the DOFs in the pressure zone
    ALLOCATE(nodes_c_vec(n_c))                                                                              ! array for storing the node numbers corresponding to the DOFs in the cavitation zone
    ALLOCATE(Phi_mat(n_p,n_mod))                                                                            ! array for storing the modal matrix (matrix of eigenvectors) Phi
    ALLOCATE(lambda_vec(n_mod))                                                                             ! array for storing the eigenvalues
    ALLOCATE(C_vec(n_mod))                                                                                  ! array for storing the integration constants
    Phi_mat = 0.0d0                                                                                         ! initialize as zero (Phi might be block-diagonal, and the entries outside these blocks will not be edited but must be zero)
    Pi_con_vec = Pi_os                                                                                      ! the array Pi_con_vec for storing the constant component of the solution (constant in the axial direction) globally (both zones + oil supply) is initialized with the dimensionless supply pressure
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Analyze which nodes are located in which flow regime
    ! ------------------------------------------------------------------------------------------------------
    
    j = 1
    k = 1
    DO i = 1, dof                                                                                           ! loop through all DOFs
      l = nodes_dof_vec(i)
      IF (g_vec(l) .EQ. 1) THEN                                                                             ! IF the current DOF is in the pressure zone, then ...
        nodes_p_vec(j) = l                                                                                  ! ... save the node number in nodes_p_vec
        j = j+1
      ELSE                                                                                                  ! otherwise (cavitation zone), ...
        nodes_c_vec(k) = l                                                                                  ! ... save the node number in nodes_c_vec
        k = k+1
      END IF
    END DO
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Determine number of pressure zones
    ! ------------------------------------------------------------------------------------------------------
    
    IF ( n_p .GE. 1 ) THEN                                                                                  ! "IF at least 1 DOF is in a pressure zone, then ..."
      n_pz = 1                                                                                              ! this variable will count the number of pressure zones; initialize as 1
    ELSE
      n_pz = 0                                                                                              ! there is no pressure zone
    END IF
    DO i = 2, n_p                                                                                           ! loop through nodes in the pressure zones
      IF ( nodes_p_vec(i)-nodes_p_vec(i-1) .GT. 1 ) THEN                                                    ! "IF there is a jump in node number, then ..."
        n_pz = n_pz + 1                                                                                     ! ... increase the number of encountered pressure zone by 1
      END IF
    END DO
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Allocate arrays whoose size depends on the number of pressure zones
    ! ------------------------------------------------------------------------------------------------------
    
    ALLOCATE(start_p_vec(n_pz))                                                                             ! allocate array for storing the first node number of every pressure zone
    ALLOCATE(end_p_vec(n_pz))                                                                               ! allocate array for storing the last node number of every pressure zone
    ALLOCATE(indx_sp_vec(n_pz))                                                                             ! allocate array for storing the indices of nodes_p_vec corresponding to the nodes in start_p_vec
    ALLOCATE(indx_ep_vec(n_pz))                                                                             ! allocate array for storing the indices of nodes_p_vec corresponding to the nodes in end_p_vec
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Analyze which nodes are in which pressure zone
    ! ------------------------------------------------------------------------------------------------------
    
    IF ( n_pz .GE. 1 ) THEN                                                                                 ! "IF at least 1 pressure zone exists, then ..."
      indx_sp_vec(1) = 1                                                                                    ! store index of nodes_p_vec corresponding to first node of first pressure zone
      indx_ep_vec(n_pz) = n_p                                                                               ! store index of nodes_p_vec corresponding to last node of last pressure zone
    END IF
    j = 1                                                                                                   ! variable for tracking the number of the pressure zone
    DO i = 2, n_p                                                                                           ! loop through nodes in the pressure zones
      IF ( nodes_p_vec(i)-nodes_p_vec(i-1) .GT. 1 ) THEN                                                    ! "IF there is a jump in node number, then ..."
        indx_sp_vec(j+1) = i                                                                                ! store index of nodes_p_vec corresponding to first node of the current pressure zone
        indx_ep_vec(j) = i-1                                                                                ! store index of nodes_p_vec corresponding to last node of the previous pressure zone
        j = j + 1                                                                                           ! increase number of encountered pressure zones by one
      END IF
    END DO
    start_p_vec = nodes_p_vec(indx_sp_vec)                                                                  ! store numbers of the first nodes of all pressure zones
    end_p_vec = nodes_p_vec(indx_ep_vec)                                                                    ! store numbers of the last nodes of all pressure zones
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Loop through pressure zones
    ! ------------------------------------------------------------------------------------------------------
    
    DO i = 1, n_pz                                                                                          ! loop through pressure zones
      
      
      ! ----------------------------------------------------------------------------------------------------
      ! Allocate arrays with size dependent on size of current pressure zone
      ! ----------------------------------------------------------------------------------------------------
      
      n1 = start_p_vec(i)                                                                                   ! first node number in the current pressure zone
      n2 = end_p_vec(i)                                                                                     ! last node number in the current pressure zone
      n_cpz = n2-n1+1                                                                                       ! number of nodes in the current pressure zone
      ALLOCATE(trafo_i_vec(n_cpz))                                                                          ! array for storing the diagonal of the diagonal matrix which transforms the eigenvalue problem
      IF ( grooves .GE. 1 ) THEN                                                                            ! we can only have more than one pressure zone if there is at least one groove (note that the Elrod solution is only allowed with at least one groove)
        ALLOCATE(diag_B_i_vec(n_cpz))                                                                       ! array for storing the main diagonal of the matrix of the transformed eigenvalue problem
        ALLOCATE(Ldiag_B_i_vec(n_cpz-1))                                                                    ! array for storing the left/lower diagonal of the matrix of the transformed eigenvalue problem
        ALLOCATE(Phi_i_mat(n_cpz,n_cpz))                                                                    ! 2D array for storing the eigenvectors
        ALLOCATE(lambda_i_vec(n_cpz))                                                                       ! array for storing the eigenvalues of the current pressure zone
        ALLOCATE(Vec1(n_cpz))                                                                               ! array for temporary storage of intermediate results
        ALLOCATE(Vec2(n_cpz))                                                                               ! array for temporary storage of intermediate results
      END IF
      
      
      ! ----------------------------------------------------------------------------------------------------
      ! Transformation of eigenvalue problem from generalized to standard for current pressure zone
      ! ----------------------------------------------------------------------------------------------------
      
      IF ( ( tay .EQ. 0 ) .OR. ( epsil > eps_max ) ) THEN                                                   ! if no Taylor approximations are used (i.e., if the eigenvalue problem will be solved by an eigensolver)
        trafo_i_vec = diag_E0_vec(n1:n2)**(-0.5d0)                                                          ! diagonal of the diagonal matrix that transforms the eigenvalue problem
        IF ( grooves .EQ. 0 ) THEN                                                                          ! if no grooves exist, there is only one pressure zone and it is periodic; with this periodicity, the eigenvalue problem has cyclic tridiagonal form (instead of tridiagonal form), which requires a different treatment
          E2_mat(1,1) = E2_mat(1,1)/diag_E0_vec(1)                                                          ! top left entry of B, overwriting E2
          E2_mat(n_x,1) = E2_mat(n_x,1)*trafo_i_vec(1)*trafo_i_vec(n_x)                                     ! bottom left entry of B, overwriting E2; the top right entry is not required because the eigenvalue problem is symmetric
          DO k = 2, n_x    
            E2_mat(k,k) = E2_mat(k,k)/diag_E0_vec(k)                                                        ! main diagonal entries of B, overwriting E2
            E2_mat(k,k-1) = E2_mat(k,k-1)*trafo_i_vec(k)*trafo_i_vec(k-1)                                   ! left diagonal entries of B, overwriting E2; the right diagonal is not required because the eigenvalue problem is symmetric
          END DO
        ELSE                                                                                                ! if there is at least one cavitation zone or oil supply groove, the current pressure zone is not periodic; this allows a tridiagonal formulation of the eigenvalue problem
          diag_B_i_vec(1) = E2_mat(n1,n1)                                                                   ! start extracting the diagonal entries from E2_mat and storing them in diag_B_i_mat; first entry corresponding to the current pressure zone
          DO k = 1, (n_cpz-1)                                                                               ! loop through all nodes in the current pressure zone in order to ...
            diag_B_i_vec(k+1) = E2_mat(n1+k,n1+k)                                                           ! ... extract the main diagonal of E2pp (only the part corresponding to the current pressure zone) and store it in diag_B_i_vec
            Ldiag_B_i_vec(k) = E2_mat(n1+k,n1+k-1)                                                          ! ... extract the left/lower diagonal of E2pp (only the part corresponding to the current pressure zone) and store it in Ldiag_B_i_vec
          END DO
          diag_B_i_vec = diag_B_i_vec/diag_E0_vec(n1:n2)                                                    ! transformation of eigenvalue problem: main diagonal of B
          Ldiag_B_i_vec = Ldiag_B_i_vec*trafo_i_vec(1:(n_cpz-1))*trafo_i_vec(2:n_cpz)                       ! transformation of eigenvalue problem: left/lower diagonal of B
        END IF
      END IF
      
      
      ! ----------------------------------------------------------------------------------------------------
      ! Solution of eigenvalue problem for current pressure zone
      ! ----------------------------------------------------------------------------------------------------
      
      IF ( ( tay .EQ. 0 ) .OR. ( epsil > eps_max ) ) THEN                                                   ! if no Taylor approximations are used (i.e., if the eigenvalue problem will be solved by an eigensolver)
        ALLOCATE(iwork_vec(1))                                                                              ! allocate workspace array for workspace query
        ALLOCATE(work_vec(1))                                                                               ! allocate workspace array for workspace query
        lwork = -1                                                                                          ! set flag for workspace query
        liwork = -1                                                                                         ! set flag for workspace query
        IF ( grooves .EQ. 0 )  THEN                                                                         ! if no grooves exist, there is only one pressure zone and it is periodic; with this periodicity, the eigenvalue problem has cyclic tridiagonal form (instead of tridiagonal form), which requires a different treatment
          CALL DSYEVD('V','L',n_x,E2_mat,n_x,lambda_vec,work_vec,lwork,iwork_vec,liwork,info)               ! perform workspace query
          lwork = INT(work_vec(1))                                                                          ! workspace size according to workspace query
          liwork = iwork_vec(1)                                                                             ! workspace size according to workspace query
          DEALLOCATE(iwork_vec)                                                                             ! deallocate workspace array
          DEALLOCATE(work_vec)                                                                              ! deallocate workspace array
          ALLOCATE(iwork_vec(liwork))                                                                       ! allocate workspace array with size determined by workspace query
          ALLOCATE(work_vec(lwork))                                                                         ! allocate workspace array with size determined by workspace query
          CALL DSYEVD('V','L',n_x,E2_mat,n_x,lambda_vec,work_vec,lwork,iwork_vec,liwork,info)               ! solve eigenvalue problem; the solver stores the squared eigenvalues and the eigenvectors in lambda_vec and in E2_mat, respectively
          DEALLOCATE(iwork_vec)                                                                             ! deallocate workspace array
          DEALLOCATE(work_vec)                                                                              ! deallocate workspace array
          lambda_vec(1) = 0.0d0                                                                             ! the numerically-zero eigenvalue is set to exactly zero
          FORALL (k=1:n_x) Phi_mat(k,:) = E2_mat(k,:)*trafo_i_vec(k)                                        ! the eigenvectors (currently stored in E2_mat) are transformed back and stored in Phi_mat
          lambda_vec = SQRT(lambda_vec)                                                                     ! eigenvalues
        ELSE                                                                                                ! if there is at least one cavitation zone or oil supply groove, the current pressure zone is not periodic; this allows a tridiagonal formulation of the eigenvalue problem
          CALL DSTEDC('I',n_cpz,diag_B_i_vec,Ldiag_B_i_vec,Phi_i_mat,n_cpz,work_vec,lwork,iwork_vec, &      ! perform workspace query
            liwork, info)
          lwork = INT(work_vec(1))                                                                          ! workspace size according to workspace query
          liwork = iwork_vec(1)                                                                             ! workspace size according to workspace query
          DEALLOCATE(iwork_vec)                                                                             ! deallocate workspace array
          DEALLOCATE(work_vec)                                                                              ! deallocate workspace array
          ALLOCATE(iwork_vec(liwork))                                                                       ! allocate workspace array with size determined by workspace query
          ALLOCATE(work_vec(lwork))                                                                         ! allocate workspace array with size determined by workspace query
          CALL DSTEDC('I',n_cpz,diag_B_i_vec,Ldiag_B_i_vec,Phi_i_mat,n_cpz,work_vec,lwork,iwork_vec, &      ! solve eigenvalue problem; the solver stores the squared eigenvalues and the eigenvectors in diag_B_i_vec and in Phi_i_mat, respectively
            liwork, info)
          DEALLOCATE(iwork_vec)                                                                             ! deallocate workspace array
          DEALLOCATE(work_vec)                                                                              ! deallocate workspace array
          FORALL (k=1:n_cpz) Phi_i_mat(k,:) = Phi_i_mat(k,:)*trafo_i_vec(k)                                 ! transform eigenvectors back
          lambda_i_vec = SQRT(diag_B_i_vec)                                                                 ! eigenvalues of current pressure zone
          indx1 = indx_sp_vec(i)                                                                            ! index of nodes_p_vec corresponding to the first node of the current pressure zone
          indx2 = indx_ep_vec(i)                                                                            ! index of nodes_p_vec corresponding to the last node of the current pressure zone
          lambda_vec(indx1:indx2) = lambda_i_vec                                                            ! store eigenvalues of current pressure zone globally
          Phi_mat(indx1:indx2,indx1:indx2) = Phi_i_mat                                                      ! store eigenvectors of current pressure zone in the overall modal matrix Phi_mat; if there is more than one pressure zone, Phi_mat is block diagonal and Phi_i_mat makes up one block
        END IF
      ELSE
        Delta_eps = epsil-eps_constr                                                                        ! difference in relative eccentricity between the current shaft position and the one where the Taylor series coefficients were derived
        fac = 1.0d0
        Phi_mat = Vd_mat(:,1:n_mod)
        lambda_vec = Ld_vec(1:n_mod)
        DO j = 1, n_tay                                                                                     ! loop through all powers involved in the Taylor series
          fac = fac*Delta_eps/j                                                                             ! factor to multiply coefficient by
          Phi_mat = Phi_mat + Vd_mat(:,(1+j*n_mod):((j+1)*n_mod))*fac                                       ! add contribution of Taylor series coefficient corresponding to current power to the matrix of eigenvectors
          lambda_vec = lambda_vec + Ld_vec((1+j*n_mod):((j+1)*n_mod))*fac                                   ! add contribution of Taylor series coefficient corresponding to current power to the vector of eigenvalues
        END DO
        Phi_mat(:,1) = 1.0d0                                                                                ! first eigenvector (corresponding to the zero-eigenvalue) does not need to be approximated, as this eigenvector is theoretically constant and equal to [1 1 1 ...]^T (apart from the normalization)
        FORALL (j=1:n_mod) Phi_mat(:,j) = Phi_mat(:,j)/SQRT(SUM(Phi_mat(:,j)**2*diag_E0_vec))               ! normalization
        lambda_vec(1) = 0.0d0                                                                               ! the numerically-zero eigenvalue is set to exactly zero
        lambda_vec = SQRT(lambda_vec)                                                                       ! eigenvalues
      END IF
      
      
      ! ----------------------------------------------------------------------------------------------------
      ! Analyze interaction between current pressure zone and neighboring cavitation zones
      ! ----------------------------------------------------------------------------------------------------
      
      IF ( grooves .EQ. 0 ) THEN                                                                            ! if no grooves exist, this algorithm does not permit a simulation with cavitation (only Guembel), meaning that no cavitation zone is present right now --> "IF no grooves exist, then ..."
        cav1 = 0                                                                                            ! ... set cav1 = 0, which means that this pressure zone is not preceded by a cavitation zone
      ELSEIF ( states_vec(n1-1) .EQ. 0 ) THEN                                                               ! "otherwise, and IF this pressure zone is directly preceded by a cavitation zone, then ..."
        cav1 = 1                                                                                            ! ... save this information, setting the flag cav1 = 1
      ELSE                                                                                                  ! "otherwise, ..." (in this case, there must be an oil supply in front of this pressure zone)
        cav1 = 0                                                                                            ! ... set cav1 = 0, indicating that no cavitation zone directly precedes this pressure zone
      END IF
      IF ( n2 .EQ. n_x ) THEN                                                                               ! "IF this pressure zone reaches the end of the computational domain, then ..."
        cav2 = 0                                                                                            ! ... set the flag cav2 = 0, indicating that this pressure zone is not directly succeeded by a cavitation zone (due to the shift of node numbers conducted earlier, the end of the domain always marks the onset of an oil supply groove)
      ELSEIF ( states_vec(n2+1) .EQ. 0 ) THEN                                                               ! "IF this pressure zone is directly succeeded by a cavitation zone, then ..." (the alternative is that behind this pressure zone there is an oil supply groove)
        cav2 = 1                                                                                            ! ... save this information setting the flag cav2 = 1
      ELSE                                                                                                  ! otherwise,
        cav2 = 0                                                                                            ! ... set cav2 = 0, indicating that no cavitation zone directly succeeds this pressure zone
      END IF
      IF ( cav1 .EQ. 1 ) THEN                                                                               ! "IF this pressure zone is directly preceded by a cavitation zone, then ..."
        entry_c = E2_mat(n1,n1-1)                                                                           ! entry of E2pc coupling current pressure zone with preceding cavitation zone; regarding the notation in the comments: E2 is formally subdivided into E2pp and E2cc, which represent the influences of the flow regimes on themselves, as well as E2pc and E2cp, which represent the influences of the flow regimes on each other; similarly, R and Pi_con are formally subdivided into Rp and Rc as well as Pi_con_p and Pi_con_c, respectively
        entry_a = E2_mat(n1-1,n1)                                                                           ! entry of E2cp coupling current pressure zone with preceding cavitation zone
      END IF
      IF ( cav2 .EQ. 1 ) THEN                                                                               ! "IF this pressure zone is directly succeeded by a cavitation zone, then ..."
        entry_b = E2_mat(n2+1,n2)                                                                           ! entry of E2cp coupling current pressure zone with succeeding cavitation zone
      END IF                                                                                                ! from here on, it gets a bit abstract; what happens in the following is basically an efficient way of handling the series of matrix multiplications for the construction of the equation system for Pi_con_c (check the referenced thesis) under exploitation of the extreme sparsity of E2cp and E2pc as well as - if 2 or more pressure zones exist - the block-diagonal nature of the global variants of Phi and E2pp^-1
      IF ( grooves .GE. 1 ) THEN                                                                            ! "IF oil supply grooves are present, ..." (otherwise, the steps below aren't necessary, since only Guembel is used for this case and we need a different formulation anyway, due to the zero eigenvalue resulting from the absence of oil supply BCs)
        Vec1 = R_vec(n1:n2)                                                                                 ! Vec1 = Rp (actually, only the part of Rp corresponding to the current pressure zone, that is, Rp_i, but let's generally omit this _i when discussing the variables here in the comments)
        CALL DGEMV('T',n_cpz,n_cpz,alpha,Phi_i_mat,n_cpz,Vec1,incr,beta,Vec2,incr)                          ! Vec2 = Phi^T * Vec1 = Phi^T * Rp
        Vec2 = Vec2/diag_B_i_vec                                                                            ! Vec2 = Vec2 * Lambda^-2 = Lambda^-2 * Vec2 = Lambda^-2 * Phi^T * Rp; note that diag_B_i_vec currently containes the squared eigenvalues
        CALL DGEMV('N',n_cpz,n_cpz,alpha,Phi_i_mat,n_cpz,Vec2,incr,beta,Vec1,incr)                          ! Vec1 = Phi * Vec2 = Phi * Lambda^-2 * Phi^T * Rp = E2pp^-1 * Rp; this will be needed within the next IF-clause; we are already computing it here (outside the IF-clause) because it will be needed also at another point of the algorithm, namely for the computation of Pi_con_p; thus, we store it for later (see next line)
        Pi_con_vec(n1:n2) = -Vec1                                                                           ! we abuse the array Pi_con_vec to store -Vec1 because E2pp^-1 * Rp will be needed again much later (when Vec1 is already deleted), but with a negative sign; Pi_con_vec(n1:n2) did not contain any important values yet
      END IF
      n_ab = cav1+cav2                                                                                      ! determines how many of the considered entries of E2cp (namely entry_a and entry_b) actually exist (0, 1, or 2)
      IF ( n_ab .GE. 1 ) THEN                                                                               ! "IF the current pressure zone interacts with at least one cavitation zone, then ..."
        ALLOCATE(Mat1(n_ab,n_cpz))                                                                          ! array for temporary storage of intermediate results
        ALLOCATE(Mat2(n_ab,n_cpz))                                                                          ! array for temporary storage of intermediate results
        ALLOCATE(Vec3(n_ab))                                                                                ! array for temporary storage of intermediate results
        IF ( cav1 .EQ. 1 ) THEN                                                                             ! "IF this pressure zone is directly preceded by a cavitation zone, then ..."
          Mat1(1,:) = entry_a*Phi_i_mat(1,:)                                                                ! Mat1 = E2cp * Phi, but note that E2cp * Phi may have a lot of rows containing only zeros and that these are excluded; this line computes the nonzero row of E2cp * Phi corresponding to the interaction with the preceding cavitation zone
        END IF
        IF ( cav2 .EQ. 1 ) THEN                                                                             ! "IF this pressure zone is directly succeeded by a cavitation zone, then ..."
          Mat1(n_ab,:) = entry_b*Phi_i_mat(n_cpz,:)                                                         ! Mat1 = E2cp * Phi, but note that E2cp * Phi may have a lot of rows containing only zeros and that these are excluded; this line computes the nonzero row of E2cp * Phi corresponding to the interaction with the succeeding cavitation zone
        END IF
        Vec2 = TANH(lambda_i_vec)/lambda_i_vec                                                              ! Vec2 = diagonal of Lambda^-1 * tanh(Lambda)
        !FORALL (j=1:n_cpz) Mat1(:,j) = Mat1(:,j)*Vec2(j)                                                   ! replaced by the line below (which only works because Fortran doesn't care whether an array represents a row or a column)
        FORALL (j=1:n_ab) Mat1(j,:) = Mat1(j,:)*Vec2                                                        ! Mat1 = Mat1 * diag(Vec2) = Mat1 * Lambda^-1 * tanh(Lambda) = E2cp * Phi * Lambda^-1 * tanh(Lambda), but only the (maximum two) nonzero rows
        CALL DGEMM('N','T',n_ab,n_cpz,n_cpz,alpha,Mat1,n_ab,Phi_i_mat,n_cpz,beta,Mat2,n_ab)                 ! Mat2 = Mat1 * Phi^T = E2cp * Phi * Lambda^-1 * tanh(Lambda) * Phi^T, but only the nonzero rows
        !FORALL (j=1:n_cpz) Mat2(:,j) = -Mat2(:,j)*diag_E0_vec(n1-1+j)                                      ! replaced by the line below (which only works because Fortran doesn't care whether an array represents a row or a column)
        FORALL (j=1:n_ab) Mat2(j,:) = -Mat2(j,:)*diag_E0_vec(n1:n2)                                         ! Mat2 = - Mat2 * E0 = - E2cp * Phi * Lambda^-1 * tanh(Lambda) * Phi^T * E0 = - E2cp * Phi * Lambda^-1 * tanh(Lambda) * Phi^-1, but only the nonzero rows
        IF ( cav1 .EQ. 1 ) THEN                                                                             ! "IF this pressure zone is directly preceded by a cavitation zone, then ..."
          Mat2(1,1) = Mat2(1,1) + entry_a                                                                   ! Mat2 = E2cp + Mat2 = E2cp - E2cp * Phi * Lambda^-1 * tanh(Lambda) * Phi^-1, but only the nonzero rows; this line considers the part of the first term E2cp corresponding to the preceding cavitation zone
        END IF
        IF ( cav2 .EQ. 1 ) THEN                                                                             ! "IF this pressure zone is directly succeeded by a cavitation zone, then ..."
          Mat2(n_ab,n_cpz) = Mat2(n_ab,n_cpz) + entry_b                                                     ! Mat2 = E2cp + Mat2 = E2cp - E2cp * Phi * Lambda^-1 * tanh(Lambda) * Phi^-1, but only the nonzero rows; this line considers the part of the first term E2cp corresponding to the succeeding cavitation zone
        END IF
        CALL DGEMV('N',n_ab,n_cpz,alpha,Mat2,n_ab,Vec1,incr,beta,Vec3,incr)                                 ! Vec3 = Mat2 * Vec1 = ( E2cp - E2cp * Phi * Lambda^-1 * tanh(Lambda) * Phi^-1 ) * E2pp^-1 * Rp, but only the (maximum two) nonzero entries; Vec3 now describes how the current pressure zone affects the neighboring cavitation zones through the RHS
        IF ( cav1 .EQ. 1 ) THEN                                                                             ! "IF this pressure zone is directly preceded by a cavitation zone, then ..."
          R_vec(n1-1) = R_vec(n1-1) - Vec3(1)                                                               ! manipulate RHS vector in preceding cavitation zone to consider influence of pressure zone; we may interpret this as a static condensation of the overall equation system (for Pi_con and C) into a smaller equation system for the cavitation regime (for Pi_con_c); to avoid allocating a new vector, we simply edit parts of the existing one; after this manipulation, all further evaluation of R_vec must be handled with caution
        END IF
        IF ( cav2 .EQ. 1 ) THEN                                                                             ! "IF this pressure zone is directly succeeded by a cavitation zone, then ..."
          R_vec(n2+1) = R_vec(n2+1) - Vec3(n_ab)                                                            ! manipulate RHS vector in succeeding cavitation zone to consider influence of pressure zone; we may interpret this as a static condensation of the overall equation system (for Pi_con and C) into a smaller equation system for the cavitation regime (for Pi_con_c); to avoid allocating a new vector, we simply edit parts of the existing one; after this manipulation, all further evaluation of R_vec must be handled with caution
        END IF
        IF ( cav1 .EQ. 1 ) THEN                                                                             ! "IF this pressure zone is directly preceded by a cavitation zone, then ..."
          !FORALL (j=1:n_cpz) Vec1(j) = entry_c*Phi_i_mat(1,j)                                              ! replaced by the line below (which only works because Fortran doesn't care whether an array represents a row or a column)
          Vec1 = entry_c*Phi_i_mat(1,:)                                                                     ! Vec1 = Phi^T * E2pc; only one column of the matrix resulting from this multiplication is nonzero, we store this nonzero column as a vector
          Vec1 = Vec1/diag_B_i_vec                                                                          ! Vec1 = Vec1 * Lambda^-2 = Lambda^-2 * Vec1 = Lambda^-2 * Phi^T * E2pc, but we only consider the nonzero column; note that diag_B_i_vec currently containes the squared eigenvalues
          CALL DGEMV('N',n_cpz,n_cpz,alpha,Phi_i_mat,n_cpz,Vec1,incr,beta,Vec2,incr)                        ! Vec2 = Phi * Vec1 = Phi * Lambda^-2 * Phi^T * E2pc = E2pp^-1 * E2pc, but we only consider the nonzero column
          diag_E2_vec(n1:n2) = Vec2                                                                         ! we abuse the array diag_E2_vec to store Vec2 because E2pp^-1 * E2pc will be needed again much later (when Vec2 is already deleted); the values that were contained in diag_E2_vec(n1:n2) before aren't needed anymore
          CALL DGEMV('N',n_ab,n_cpz,alpha,Mat2,n_ab,Vec2,incr,beta,Vec3,incr)                               ! Vec3 = Mat2 * Vec2 = ( E2cp - E2cp * Phi * Lambda^-1 * tanh(Lambda) * Phi^-1 ) * E2pp^-1 * E2pc, but only the (maximum two) nonzero entries; the first entry of Vec3 now describes an influence of the current pressure zone on the preceding cavitation zone, the second entry (if existent) describes how the preceding and succeeding cavitation zones are coupled with each other through this pressure zone
          E2_mat(n1-1,n1-1) = E2_mat(n1-1,n1-1) - Vec3(1)                                                   ! manipulate E2-matrix in preceding cavitation zone to consider influence of pressure zone; we may interpret this as a static condensation of the overall equation system (for Pi_con and C) into a smaller equation system for the cavitation regime (for Pi_con_c); to avoid allocating a new matrix, we simply edit parts of the existing one; after this manipulation, all further evaluation of E2_mat must be handled with caution
          IF ( cav2 .EQ. 1 ) THEN                                                                           ! "IF this pressure zone is directly succeeded by a cavitation zone (in addition to being preceded by one), then ..."
            E2_mat(n2+1,n1-1) = - Vec3(2)                                                                   ! add coupling term between preceding and succeeding cavitation zones to E2-matrix; we may interpret this as a static condensation of the overall equation system (for Pi_con and C) into a smaller equation system for the cavitation regime (for Pi_con_c); to avoid allocating a new matrix, we simply edit parts of the existing one; after this manipulation, all further evaluation of E2_mat must be handled with caution
          END IF
        END IF
      END IF
      
      
      ! ----------------------------------------------------------------------------------------------------
      ! Deallocate arrays
      ! ----------------------------------------------------------------------------------------------------
      
      DEALLOCATE(trafo_i_vec)                                                                               ! size depends on size of current pressure zone
      IF ( grooves .GE. 1 ) THEN                                                                            ! we can only have more than one pressure zone if there is at least one groove (note that the Elrod solution is only allowed with at least one groove)
        DEALLOCATE(diag_B_i_vec)                                                                            ! size depends on size of current pressure zone
        DEALLOCATE(Ldiag_B_i_vec)                                                                           ! size depends on size of current pressure zone
        DEALLOCATE(Phi_i_mat)                                                                               ! size depends on size of current pressure zone
        DEALLOCATE(lambda_i_vec)                                                                            ! size depends on size of current pressure zone
        DEALLOCATE(Vec1)                                                                                    ! size depends on size of current pressure zone
        DEALLOCATE(Vec2)                                                                                    ! size depends on size of current pressure zone
      END IF
      IF ( n_ab .GE. 1 ) THEN                                                                               ! the arrays below only exist if the current pressure zone interacts with at least one cavitation zone
        DEALLOCATE(Mat1)                                                                                    ! size depends on interaction of current pressure zone with cavitation zones and on size of current pressure zone
        DEALLOCATE(Mat2)                                                                                    ! size depends on interaction of current pressure zone with cavitation zones and on size of current pressure zone
        DEALLOCATE(Vec3)                                                                                    ! size depends on interaction of current pressure zone with cavitation zones
      END IF
      
      
    END DO
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Solve bidiagonal equation system in the cavitation regime
    ! ------------------------------------------------------------------------------------------------------
    
    IF ( n_c .GE. 1 ) THEN                                                                                  ! "IF there is cavitation at all, then ..."
      j = nodes_c_vec(1)                                                                                    ! determine number of first node located in a cavitation zone
      Pi_con_vec(j) = R_vec(j)/E2_mat(j,j)                                                                  ! compute solution (pressure-like function) on this node
    END IF
    DO i = 2, n_c                                                                                           ! loop through all nodes with cavitation
      j = nodes_c_vec(i)                                                                                    ! determine number of current node
      k = nodes_c_vec(i-1)                                                                                  ! determine number of previous node with cavitation
      Pi_con_vec(j) = (R_vec(j)-E2_mat(j,k)*Pi_con_vec(k))/E2_mat(j,j)                                      ! compute solution (pressure-like function) on current node, exploiting the bidiagonal form; the part of the E2 matrix corresponding to the cavitation zone is bidiagonal (each node is affected only by the previous node) due to the optimal upwind scheme
    END DO
    Pi_con_vec(nodes_c_vec) = -Pi_con_vec(nodes_c_vec)                                                      ! the sign of the matrix on the LHS of the equation system was inverted; this is now compensated via the sign of the solution vector
    Pi_bar_vec = Pi_con_vec                                                                                 ! Pi_con_vec and Pi_bar_vec will describe the axially constant component of the solution and the axially averaged solution, respectively, for all nodes; right now, these two are still identical because the solution in the pressure regime is not included yet (the SBFEM model assumes the solution in the cavitation regime - but not the one in the pressure regime - to be constant in the axial direction)
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Loop through pressure zones
    ! ------------------------------------------------------------------------------------------------------
    
    DO i = 1, n_pz                                                                                          ! loop through pressure zones
      
      
      ! ----------------------------------------------------------------------------------------------------
      ! Allocate arrays with size dependent on size of current pressure zone
      ! ----------------------------------------------------------------------------------------------------
      
      n1 = start_p_vec(i)                                                                                   ! first node number in the current pressure zone
      n2 = end_p_vec(i)                                                                                     ! last node number in the current pressure zone
      n_cpz = n2-n1+1                                                                                       ! number of nodes in the current pressure zone
      IF ( grooves .GE. 1 ) THEN                                                                            ! we can only have more than one pressure zone if there is at least one groove (note that the Elrod solution is only allowed with at least one groove)
        ALLOCATE(Phi_i_mat(n_cpz,n_cpz))                                                                    ! 2D array for storing the eigenvectors of the current pressure zone
        ALLOCATE(C_i_vec(n_cpz))                                                                            ! array for storing the integration constants of the current pressure zone
        ALLOCATE(Vec2(n_cpz))                                                                               ! array for temporary storage of some intermediate results
        ALLOCATE(Vec1(n_cpz))                                                                               ! array for temporary storage of some intermediate results
      ELSE
        ALLOCATE(Vec1(n_mod))                                                                               ! array for temporary storage of some intermediate results
      END IF
      
      
      ! ----------------------------------------------------------------------------------------------------
      ! Compute solution in current pressure zone
      ! ----------------------------------------------------------------------------------------------------
      
      IF ( grooves .EQ. 0 ) THEN                                                                            ! if no grooves exist, a zero eigenvalue is present, requiring a slightly different formulation (at the same time, we will exploit the fact that the solution without grooves is always conducted under Guembel conditions)
        CALL DGEMV('T',n_x,n_mod,alpha,Phi_mat,n_x,R_vec,incr,beta,C_vec,incr)                              ! C = Phi^T * R (this is not really C yet, it's an intermediate step)
        C_vec(2:n_mod) = C_vec(2:n_mod)/(lambda_vec(2:n_mod)**2)                                            ! C = Lambda^-2 * Phi^T * R, using the intermediate result from the line above; note that we are excluding the mode that is in the null space (the mode corresponding to the zero eigenvalue)
        C_vec(1) = 0.0d0                                                                                    ! no contribution of the mode that is in the null space (this mode expresses a constant pressure distribution), assuming homogeneous Dirichlet BCs at the bearing boundary
        alpha = -1.0d0                                                                                      ! set alpha to -1, enforcing the negative sign in the multiplication below
        CALL DGEMV('N',n_x,n_mod,alpha,Phi_mat,n_x,C_vec,incr,beta,Pi_con_vec,incr)                         ! Pi_con = - Phi * C
        alpha = 1.0d0                                                                                       ! reset alpha to its usual value of 1
        Vec1(2:n_mod) = TANH(lambda_vec(2:n_mod))/lambda_vec(2:n_mod)*C_vec(2:n_mod)                        ! Vec1 = tanh(Lambda) * Lambda^-1 * C
        Vec1(1) = 0.0d0                                                                                     ! no contribution of the mode that is in the null space
        CALL DGEMV('N',n_x,n_mod,alpha,Phi_mat,n_x,Vec1,incr,beta,Pi_bar_vec,incr)                          ! Pi_bar = Phi * tanh(Lambda) * Lambda^-1 * C (this is not really Pi_bar yet, it's an intermediate step)
        Pi_bar_vec = Pi_bar_vec + Pi_con_vec                                                                ! Pi_bar = Phi * tanh(Lambda) * Lambda^-1 * C + Pi_con, using the intermediate result from the line above
        Vec1 = TANH(lambda_vec)*lambda_vec*C_vec                                                            ! Vec1 = tanh(Lambda) * Lambda * C
        alpha = 2*p_ref/l_b                                                                                 ! set additional factor alpha for multiplication below to 2*p_ref/l_b
        CALL DGEMV('N',n_x,n_mod,alpha,Phi_mat,n_x,Vec1,incr,beta,dpdy_bb_vec,incr)                         ! dpdy_bb = alpha * Phi * Vec1 = (2*p_ref/l_b) * Phi * tanh(Lambda) * Lambda * C (computation of axial pressure gradient at bearing boundary)
        alpha = 1.0d0                                                                                       ! reset alpha to usual value of 1 for future calls of DGEMV
      ELSE                                                                                                  ! otherwise (i.e., if oil supply grooves are present), we use an algorithm that assumes positive eigenvalues and is able to handle Elrod cavitation
        IF ( states_vec(n1-1) .EQ. 0 ) THEN                                                                 ! "IF this pressure zone is directly preceded by a cavitation zone, then ..."
          Pi_con_vec(n1:n2) = Pi_con_vec(n1:n2) - diag_E2_vec(n1:n2)*Pi_con_vec(n1-1)                       ! constant component of the solution in the current pressure zone Pi_con_p; the first summand -E2pp^-1 * Rp was computed earlier and already stored in Pi_con_vec(n1:n2); the second summand -E2pp^-1 * E2pc * Pi_con_c (influence of preceding cavitation zone) is considered now: the only nonzero column of E2pp^-1 * E2pc was stored in diag_E2_vec(n1:n2) earlier and the only relevant entry of Pi_con_c is given by Pi_con_vec(n1-1)
        END IF
        indx1 = indx_sp_vec(i)                                                                              ! index of nodes_p_vec corresponding to the first node of the current pressure zone
        indx2 = indx_ep_vec(i)                                                                              ! index of nodes_p_vec corresponding to the last node of the current pressure zone
        Phi_i_mat = Phi_mat(indx1:indx2,indx1:indx2)                                                        ! 2D array containing the eigenvectors of current pressure zone
        Vec1 = diag_E0_vec(n1:n2)*Pi_con_vec(n1:n2)                                                         ! Vec1 = E0 * Pi_con_p
        alpha = -1.0d0                                                                                      ! set additional factor alpha for multiplication below to -1
        CALL DGEMV('T',n_cpz,n_cpz,alpha,Phi_i_mat,n_cpz,Vec1,incr,beta,C_i_vec,incr)                       ! C = alpha * Phi^T * Vec1 = (-1) * Phi^T * E0 * Pi_con_p = - Phi^-1 * Pi_con_p (computation of the integration constants)
        alpha = 1.0d0                                                                                       ! reset alpha to usual value of 1 for future calls of DGEMV
        C_vec(indx1:indx2) = C_i_vec                                                                        ! store integration constants
        Vec1 = TANH(lambda_vec(indx1:indx2))/lambda_vec(indx1:indx2)*C_i_vec                                ! Vec1 = tanh(Lambda) * Lambda^-1 * C
        CALL DGEMV('N',n_cpz,n_cpz,alpha,Phi_i_mat,n_cpz,Vec1,incr,beta,Vec2,incr)                          ! Vec2 = Phi * Vec1 = Phi * tanh(Lambda) * Lambda^-1 * C (computation of the axial average of the axially variable component of the solution)
        Pi_bar_vec(n1:n2) = Pi_con_vec(n1:n2) + Vec2                                                        ! axial average of the solution in the current pressure zone: axially constant component + axial average of the axially non-constant component
        Vec1 = Vec1*lambda_vec(indx1:indx2)**2                                                              ! Vec1 = Vec1 * Lambda^2 = tanh(Lambda) * Lambda * C
        alpha = 2*p_ref/l_b                                                                                 ! set additional factor alpha for multiplication below to 2*p_ref/l_b
        CALL DGEMV('N',n_cpz,n_cpz,alpha,Phi_i_mat,n_cpz,Vec1,incr,beta,Vec2,incr)                          ! Vec2 = alpha * Phi * Vec1 = (2*p_ref/l_b) * Phi * tanh(Lambda) * Lambda * C (computation of axial pressure gradient at bearing boundary)
        alpha = 1.0d0                                                                                       ! reset alpha to usual value of 1 for future calls of DGEMV
        dpdy_bb_vec(n1:n2) = Vec2                                                                           ! save axial gradient at bearing boundary for current pressure zone
      END IF
      
      
      ! ----------------------------------------------------------------------------------------------------
      ! Deallocate arrays
      ! ----------------------------------------------------------------------------------------------------
      
      IF ( grooves .GE. 1 ) THEN                                                                            ! we can only have more than one pressure zone if there is at least one groove (note that the Elrod solution is only allowed with at least one groove)
        DEALLOCATE(Phi_i_mat)                                                                               ! size depends on size of current pressure zone
        DEALLOCATE(C_i_vec)                                                                                 ! size depends on size of current pressure zone
        DEALLOCATE(Vec2)                                                                                    ! size depends on size of current pressure zone
      END IF
      DEALLOCATE(Vec1)                                                                                      ! size depends on size of current pressure zone
      
      
    END DO
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Compute switch function and check convergence
    ! ------------------------------------------------------------------------------------------------------
    
    g_old_old_vec = g_old_vec                                                                               ! the switch function that was used in the previous iterations is saved
    g_old_vec = g_vec                                                                                       ! the switch function that was used in the current iteration is saved
    g_vec = 0                                                                                               ! the new switch function (which will be used in the next iteration) is computed: it is first initialized as zero, then ...
    WHERE(Pi_bar_vec .GE. 0.0d0) g_vec = 1                                                                  ! ... it is set to 1 where the axial average of the pressure-like function Pi_bar is non-negative
    IF (MAXVAL(ABS(g_vec-g_old_vec)) .EQ. 0) THEN                                                           ! "IF the new switch function is equal to the old one, then ..."
      convergent = 1                                                                                        ! the solution is flagged as convergent
    END IF
    IF ( (iter .GE. 2) .AND. (MAXVAL(ABS(g_vec-g_old_old_vec)) .EQ. 0) ) THEN                               ! "IF the same switch function has been computed in this iteration as two iterations before (indicating that the solution oscillates between two alternating configurations - although I doubt that this even happens under the current model assumptions), then ..."
      convergent = 1                                                                                        ! the solution is flagged as convergent anyway (the solution is as good as it gets)
    END IF
    
    
    ! ------------------------------------------------------------------------------------------------------
    ! Exit loop if this was the last iteration, deallocate arrays for next iteration otherwise
    ! ------------------------------------------------------------------------------------------------------
    
    IF ( (convergent .EQ. 1) .OR. (iter .EQ. iter_max) .OR. (guembel .EQ. 1) ) THEN                         ! "IF this is the last iteration (solution convergent or max allowed number of iterations reached or Guembel conditions assumed), then ..."
      EXIT                                                                                                  ! exit the loop
    ELSE                                                                                                    ! "otherwise" (meaning that there will be a next iteration; in that case, all arrays that may need to be allocated again in the next iteration with different sizes need to be deallocated now)
      DEALLOCATE(start_p_vec)                                                                               ! size depends on number of pressure zones
      DEALLOCATE(end_p_vec)                                                                                 ! size depends on number of pressure zones
      DEALLOCATE(indx_sp_vec)                                                                               ! size depends on number of pressure zones
      DEALLOCATE(indx_ep_vec)                                                                               ! size depends on number of pressure zones
      DEALLOCATE(nodes_p_vec)                                                                               ! size depends on number of DOFs per flow regime
      DEALLOCATE(nodes_c_vec)                                                                               ! size depends on number of DOFs per flow regime
      DEALLOCATE(Phi_mat)                                                                                   ! size depends on number of DOFs per flow regime
      DEALLOCATE(lambda_vec)                                                                                ! size depends on number of DOFs per flow regime
      DEALLOCATE(C_vec)                                                                                     ! size depends on number of DOFs per flow regime
    END IF
    
    
  END DO
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Save some data for the next time step, considering original node numbers
  ! --------------------------------------------------------------------------------------------------------
  
  IF ( tay .EQ. 1 ) THEN                                                                                    ! if the Taylor approximations were used for the eigenvalue problem, the computational grid was shifted (setting X_att=0); this now requires transforming some output variables back, which is done below for pts_vec
    di = MODULO((X_att_trafo/(2*pi)),1.0d0)*n_x                                                             ! number of elements between node 1 and the circumferential position defined by the attitude angle
    fl = di-FLOOR(di)                                                                                       ! weight for left node in interpolation
    fr = 1.0d0-fl                                                                                           ! weight for right node in interpolation
    di_floor = INT(FLOOR(di))                                                                               ! number of elements between node 1 and the circumferential position defined by the attitude angle, rounded down to an integer
    pts_vec(1:n_x) = (/&                                                                                    ! interpolation of pts_vec, so that the result is valid for the actual attitude angle (which was stored in X_att_trafo before setting X_att=0)
      fr*(/Pi_bar_vec((n_x-di_floor+1):n_x),Pi_bar_vec(1)/)+fl*Pi_bar_vec((n_x-di_floor):n_x),&
      fr*Pi_bar_vec(2:(n_x-di_floor))+fl*Pi_bar_vec(1:(n_x-di_floor-1))/)
  ELSE                                                                                                      ! if no Taylor approximations were used (implying that there is at least one oil supply groove), the node numbers were shiftet to ensure a tridiagonal eigenvalue problem; this shift must now be reversed for several output variables (below, this is done for pts_vec)
    pts_vec(shift_vec) = Pi_bar_vec                                                                         ! the nodal pressure-like functions, averaged in the axial direction, are stored in pts_vec in order to be available at the next time step; moreover, the shift of node numbers performed earlier is reversed    
  END IF
  pts_vec(n_x+1) = t                                                                                        ! the current time is stored in pts_vec in order to be available at the next time step
  

  ! --------------------------------------------------------------------------------------------------------
  ! Compute film fractions, averaged physical pressures, and circumferential gradients
  ! --------------------------------------------------------------------------------------------------------
  
  theta_vec = (1-g_vec)*Pi_bar_vec + 1.0d0                                                                  ! array containing the nodal film fractions
  Pi_bar_vec = p_ref*(Pi_bar_vec*g_vec)                                                                     ! the values in Pi_bar_vec aren't needed anymore; the array is overwritten and will now contain the nodal physical pressures averaged in the axial direction [Pa]
  L_X = L_X*d_b/2                                                                                           ! conversion of nondimensionalized circumferential sector length to physical circumferential sector length (i.e. in [m]) 
  IF ( guembel .EQ. 1 ) THEN                                                                                ! if Guembel is used instead of Elrod, ...
    theta_vec = 1                                                                                           ! ... set film fractions to 1 and ...
    dpdy_bb_vec = dpdy_bb_vec*g_vec                                                                         ! ... apply Guembel condition to axial pressure gradients at bearing boundary (zero in the cavitation zone)
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Compute bearing forces, friction moment, oil volume, and flow through the bearing boundary
  ! --------------------------------------------------------------------------------------------------------
  
  F_1 = SUM(COS(X_vec)*Pi_bar_vec)*l_b*L_X                                                                  ! component of hydrodynamic force (note that Pi_bar_vec does not contain pressure-like functions anymore, but physical pressures, averaged in the axial direction)
  F_2 = SUM(SIN(X_vec)*Pi_bar_vec)*l_b*L_X                                                                  ! component of hydrodynamic force
  IF ( tay .EQ. 1 ) THEN                                                                                    ! if the eigenvalue problem was solved by Taylor approximations, the attitude angle was manipulated (set to zero); this is now compensated via transforming the forces, see below
    F_h = F_1*cos(X_att_trafo) - F_2*sin(X_att_trafo)                                                       ! transform first force component (F_h is only used for temporary storage)
    F_v = F_1*sin(X_att_trafo) + F_2*cos(X_att_trafo)                                                       ! transform second force component (F_v is only used for temporary storage)
    F_1 = F_h                                                                                               ! save transformed force component (this force component now corresponds to the correct attitude angle)
    F_2 = F_v                                                                                               ! save transformed force component (this force component now corresponds to the correct attitude angle)
  END IF
  M_fr = (d_b/2*l_b*L_X)*((-0.25d0*(c/L_X))*SUM(H_vec*(Pi_bar_vec(index_E_vec)-Pi_bar_vec(index_W_vec)))&   ! computation of the oil friction moment; note that Pi_bar_vec does not contain pressure-like functions anymore, but physical pressures (averaged in the axial direction)
    +(u_pos*mu_ref/c)*SUM(theta_vec*(mu_rel_vec/H_vec)))
  V_oil = (l_b*L_X*c)*SUM(theta_vec*H_vec)                                                                  ! computation of the oil volume in the bearing, considering cavitation
  dpdy_bb_vec(nodes_c_vec) = 0.0d0                                                                          ! the axial gradients in the cavitation zone are zero
  dpdy_bb_vec(nodes_os_vec) = 0.0d0                                                                         ! the axial gradients of the supply pressure are zero
  V_dot_bb = (L_X/(6.0d0*mu_ref)*c**3)*SUM(H_vec**3*(dpdy_bb_vec/mu_rel_vec))                               ! computation of the oil volume flow through the bearing boundaries (a negative value indicates that the oil is leaving the bearing)
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Evaluation of the solution (pressure-like function) on a twodimensional grid, if desired
  ! --------------------------------------------------------------------------------------------------------
  
  IF (n_y .GE. 2) THEN                                                                                      ! "IF the number of desired data points in the axial direction is at least 2, then ..."
    j = ANINT(n_y/2.0d0+0.75)                                                                               ! axial starting index for evaluation of the solution (as one half of the bearing is initially left out because of the symmetry)
    l = n_y-j+1                                                                                             ! number of data points in the axial direction per half of the bearing
    FORALL (i=j:n_y) Pi_mat(:,i) = Pi_con_vec                                                               ! the array Pi_mat (which will contain the pressure-like function evaluated on a twodimensional grid) is initialized with the axially constant component of the solution (the non-constant component is added if a pressure zone is present)
    IF (n_p .GE. 1) THEN                                                                                    ! "IF a pressure zone is present, then ..."
      ALLOCATE(xi_vec(l))                                                                                   ! this array will be used for storing the xi-positions (xi is a dimensionless axial coordinate used by the SBFEM model) where the solution will be evaluated
      ALLOCATE(cosh_mat(n_mod,l))                                                                           ! this array will be used for storing the values of cosh(lambda_i*xi_j) for all i-j-combinations
      ALLOCATE(Pip_mat(n_p,l))                                                                              ! the component of the solution that is variable in the axial direction (pressure zone only) will be stored in this array before being added to Pi_mat
      FORALL (i=j:n_y) xi_vec(i-j+1) = (i-1.0d0)/(n_y-1.0d0)*2.0d0-1.0d0                                    ! computation of the xi-values
      FORALL (i=1:n_mod) Phi_mat(:,i) = Phi_mat(:,i)*(C_vec(i)/COSH(lambda_vec(i)))                         ! the array Phi_mat is repurposed for storing the matrix Phi*diag(C)*cosh(Lambda)^-1
      FORALL (i=1:l) cosh_mat(:,i) = COSH(lambda_vec*xi_vec(i))                                             ! cosh_mat now contains a matrix that provides the values of cosh(lambda_i*xi_j) for all i-j-combinations
      CALL DGEMM('N','N',n_p,l,n_mod,alpha,Phi_mat,n_p,cosh_mat,n_mod,beta,Pip_mat,n_p)                     ! the component of the solution that is variable in the axial direction is evaluated at all xi-positions and stored in Pip_mat; to this end, the matrix Phi*diag(C)*cosh(Lambda)^-1 (given by Phi_mat) is multiplied by the matrix of cosh(lambda_i*xi_j)-values (given by cosh_mat)
      Pi_mat(nodes_p_vec,j:n_y) = Pi_mat(nodes_p_vec,j:n_y) + Pip_mat(:,:)                                  ! the component of the solution that is variable in the axial direction (given by Pip_mat) is added to the overall solution Pi_mat
      DEALLOCATE(cosh_mat)
      DEALLOCATE(Pip_mat)
      DEALLOCATE(xi_vec)
    END IF
    IF ( tay .EQ. 1 ) THEN                                                                                  ! if the Taylor approximations were used for the eigenvalue problem, the computational grid was shifted (setting X_att=0); this now requires transforming some output variables back, which is done below for Pi_mat
      FORALL (i=j:n_y) Pi_mat(1:n_x,i) = (/&                                                                ! interpolation of Pi_mat, so that the result is valid for the actual attitude angle (which was stored in X_att_trafo before setting X_att=0)
        fr*(/Pi_mat((n_x-di_floor+1):n_x,i),Pi_mat(1,i)/)+fl*Pi_mat((n_x-di_floor):n_x,i),&
        fr*Pi_mat(2:(n_x-di_floor),i)+fl*Pi_mat(1:(n_x-di_floor-1),i)/)
    ELSE                                                                                                    ! if no Taylor approximations were used (implying that there is at least one oil supply groove), the node numbers were shiftet to ensure a tridiagonal eigenvalue problem; this shift must now be reversed for several output variables (below, this is done for Pi_mat)
      Pi_mat(shift_vec,j:n_y) = Pi_mat(:,j:n_y)                                                             ! the shift of node numbers performed earlier is now reversed
    END IF
    Pi_mat(:,l:1:-1) = Pi_mat(:,j:n_y)                                                                      ! the computed solution is copied to the other half of the bearing
  ELSE                                                                                                      ! "IF the number of desired data points in the axial direction is less than 2, then ..." (which is interpreted as a flag indicating that the evaluation of the solution on a twodimensional grid should be skipped)
    Pi_mat = 0.0                                                                                            ! Pi_mat is set to zero                                            
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Shift node numbers back
  ! --------------------------------------------------------------------------------------------------------
  
  IF ( tay .EQ. 1 ) THEN                                                                                    ! if the Taylor approximations were used for the eigenvalue problem, the computational grid was shifted (setting X_att=0); the axial average of the pressure-like function (contained in pts_vec) was already transformed back via an interpolation, now g_vec and theta_vec need to be adjusted
    g_vec = 0
    WHERE( pts_vec(1:n_x) .GE. 0.0d0 ) g_vec = 1                                                            ! adjust switch function to interpolated solution
    theta_vec = (1-g_vec)*pts_vec(1:n_x) + 1.0d0                                                            ! adjust film fraction to interpolated
  ELSE                                                                                                        ! if no Taylor approximations were used (implying that there is at least one oil supply groove), the node numbers were shiftet to ensure a tridiagonal eigenvalue problem; this shift must now be reversed for several output variables (below, this is done for Pi_mat)
    g_vec(shift_vec) = g_vec                                                                                ! prepare array of nodal switch functions for output: the shift of node numbers performed earlier is reversed
    theta_vec(shift_vec) = theta_vec                                                                        ! prepare array of nodal film fractions for output: the shift of node numbers performed earlier is reversed
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Reverse rotation of coordinate system
  ! --------------------------------------------------------------------------------------------------------
  
  IF ( u .LT. 0.0d0 ) THEN                                                                                  ! "IF the circumferential surface velocity is negative, then ..."
    F_2 = -F_2                                                                                              ! transform hydrodynamic force F2 back (F1 is not affected by the rotation)
    M_fr = -M_fr                                                                                            ! transform hydrodynamic friction moment back
    index_0_vec(1) = 1                                                                                      ! this array will describe how the rotation of the coordinate system changes the node numbers
    index_0_vec(2:n_x) = (/(i,i=n_x,2,-1)/)                                                                 ! this array will describe how the rotation of the coordinate system changes the node numbers
    g_vec = g_vec(index_0_vec)                                                                              ! change node numbers (i.e. order of the entries) back to normal in the output variable g_vec
    theta_vec = theta_vec(index_0_vec)                                                                      ! change node numbers (i.e. order of the entries) back to normal in the output variable theta_vec
    pts_vec(1:n_x) = pts_vec(index_0_vec)                                                                   ! change node numbers (i.e. order of the entries) back to normal in the output variable pts_vec
    Pi_mat = Pi_mat(index_0_vec,:)                                                                          ! change node numbers (i.e. order of the rows) back to normal in the output variable Pi_mat
  END IF
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Transform hydrodynamic forces
  ! --------------------------------------------------------------------------------------------------------
  
  F_h = F_1*COS(angle_shell) - F_2*SIN(angle_shell)                                                         ! transformation into inertial system: horizontal hydrodynamic force acting on the shell [N]
  F_v = F_1*SIN(angle_shell) + F_2*COS(angle_shell)                                                         ! transformation into inertial system: vertical hydrodynamic force acting on the shell [N]
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Deallocation
  ! --------------------------------------------------------------------------------------------------------
  
  DEALLOCATE(nodes_dof_vec)                                                                                 ! size depends on oil supply BCs
  DEALLOCATE(nodes_os_vec)                                                                                  ! size depends on oil supply BCs
  DEALLOCATE(start_p_vec)                                                                                   ! size depends on number of pressure zones
  DEALLOCATE(end_p_vec)                                                                                     ! size depends on number of pressure zones
  DEALLOCATE(indx_sp_vec)                                                                                   ! size depends on number of pressure zones
  DEALLOCATE(indx_ep_vec)                                                                                   ! size depends on number of pressure zones
  DEALLOCATE(nodes_p_vec)                                                                                   ! size depends on number of DOFs per flow regime
  DEALLOCATE(nodes_c_vec)                                                                                   ! size depends on number of DOFs per flow regime
  DEALLOCATE(Phi_mat)                                                                                       ! size depends on number of DOFs per flow regime
  DEALLOCATE(lambda_vec)                                                                                    ! size depends on number of DOFs per flow regime
  DEALLOCATE(C_vec)                                                                                         ! size depends on number of DOFs per flow regime
    
  
END SUBROUTINE SBFEM_ELROD
