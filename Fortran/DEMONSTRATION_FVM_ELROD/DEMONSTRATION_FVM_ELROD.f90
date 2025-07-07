! DEMONSTRATION_FVM_ELROD.f90 
!
! FUNCTIONS:
! DEMONSTRATION_FVM_ELROD
!
! AUTHOR, AFFILIATION, DATE: Simon Pfeil, OvGU Magdeburg (Germany), 30.06.2025
!
! **********************************************************************************************************
!
! PROGRAM: DEMONSTRATION_FVM_ELROD
!
! PURPOSE: This script demonstrates how to call the program FVM_ELROD (contained in the file FVM_ELROD.f90) 
!          for solving the Reynolds equation numerically.
!
! **********************************************************************************************************
!
! COMMENTS
!
! - The program FVM_ELROD solves the Reynolds equation by means of a numerical method which will be referred 
!   to as FVM (Finite Volume Method), although a discretized form of the Reynolds equation exactly identical 
!   to the one being used could also be derived based on an FEM or FDM formalism. The program is desinged to
!   be incorporated into time integration schemes, but it can also be used for quasistatic single calls (the 
!   script at hand will primarily do the latter, which should suffice for demonstration purposes).
!
! - A transient, mass-conserving cavitation model is used which is based on the assumptions proposed by
!   Kumar and Booker (1991) and is usually classified as a later variant of the well-known Elrod algorithm. 
!   Since this renders the problem nonlinear, the program FVM_ELROD may internally perform multiple 
!   iterations. To handle the transient cavitation term in the Reynolds equation, the program needs to 
!   transfer data to itself across time steps, for which it uses the input/output variable pts_vec.
!    Optionally, the cavitation model can be bypassed, meaning that Guembel conditions are assumed instead.
!
! - The input and output variables of FVM_ELROD will all be clarified throughout the script at hand.
!
! - See the attached PDF for clarification of the coordinate systems. x (or the nondimensionalized X) and y 
!   (or the nondimensionalized Y) are the coordinates of the lubrication gap and are used by the Reynolds 
!   equation. This coordinate system, as well as the computational grid, are fixed in the reference frame of 
!   the shell. The transformations between the reference frame of the shell and the inertial system are 
!   performed internally by this program. The kinematic variables of the shell and the shaft, when handed to 
!   FVM_ELROD, should be formulated from the perspective of the inertial system. The hydrodynamic forces 
!   and moments computed by FVM_ELROD are automatically transformed into the inertial system before output. 
!   In contrast, those input and output variables that describe one- or two-dimensional fields (e.g., 
!   viscosity distribution, additional contour, and pressure-like function) never use the inertial system, 
!   as they are expressed by arrays whoose entries represent the values at the nodes (which are always fixed 
!   at the shell). The circumferential position of an oil supply groove X_os is also always expressed in 
!   the reference frame of the shell; otherwise, a rotating shell would require us to keep updating this 
!   input variable.
!
! - More thorough discussions of the computational method, of the cavitation model, of the boundary 
!   conditions (BCs), and of further assumptions will be given in the thesis "Simulating Hydrodynamic 
!   Bearings with the Scaled Boundary Finite Element Method" by Simon Pfeil, but this thesis isn't 
!   publically available yet
!
! - FVM_ELROD uses BiCGStab [1] to solve the unsymmetric equation system. Under Guembel conditions, the 
!   equation system is symmetric, so ICCG is used for its solution instead. The ICCG subprogram is given 
!   by the file SOLVER_IC_CG_MKL.f90 (provided by Steffen Nitzschke, OvGU Magdeburg), which relies on BLAS 
!   routines (available via, e.g., Intel's Math Kernel Library) as well as on an IC factorization [2].
!   The subprograms [1] and [2] can be obtained from 
!   [1] Yousef Saad: SPARSKIT2. University of Minnesota, Department of Computer Science and Engineering, 
!       200 Union Street S.E., Minneapolis, MN 55455 USA, [saad -at- umn -dot- edu],
!       https://www-users.cse.umn.edu/~saad/software/SPARSKIT/
!   [2] Mark T. Jones, Paul E. Plassmann: Algorithm 740: Fortran subroutines to compute improved 
!       incomplete Cholesky factorizations. ACM Transactions on Mathematical Software, vol. 21,
!       no. 1, March 1995, p. 18-19. Association for Computing Machinery, New York, NY, USA, 0098-3500,
!       https://doi.org/10.1145/200979.200986
!    Please follow the intructions given by the comments in the files MODULE_SPARSKIT2_BCGSTAB.f90 and
!   MODULE_JPICC.f90. Note that it is also possible to modify the FVM algorithm to always use BiCGStab
!   (as opposed to switching between BiCGStab and ICCG depending on the cavitation assumptions).
!
! - The node numbering scheme is illustrated below. Note that the nodes are counted without including the
!   periodic nodes at X=2pi. The first note in the circumferential direction is located at X=0 and the last
!   one at X=2pi*(1-1/n_x), where n_x is the circumferential number of nodes. The first and last nodes in
!   the axial direction are located at the bearing boundaries. As mentioned above, the computational
!   grid is fixed in the reference frame of the shell.
!
! **********************************************************************************************************
!
! NODE NUMBERING SCHEME
!
!
!   axial direction (y [m] or Y [-]) with n_y nodes
!
!   ^
!   |
!   |       ...                                 ...      n_y*n_x
!   |
!   |       ...
!           2*n_x+1  2*n_x+2  ...
!           n_x+1    n_x+2    n_x+3    n_x+4    ...      2*n_x
!           1        2        3        4        ...      n_x
!
!           -------->  circumferential direction (x [m] or X [rad]) with n_x nodes
!
!
! **********************************************************************************************************



PROGRAM DEMONSTRATION_FVM_ELROD

  ! --------------------------------------------------------------------------------------------------------
  ! Modules and variable types
  ! --------------------------------------------------------------------------------------------------------
  
  ! no implicit variables
  IMPLICIT NONE
  
  ! input variables for FVM_ELROD: bearing properties, kinematic variables, discretization, ...
  REAL(KIND=8)                                     :: d_b, l_b
  INTEGER                                          :: n_x, n_y, grooves, iter_max_solver, quasistatic
  INTEGER                                          :: iter_max, symBC, guembel
  REAL(KIND=8)                                     :: c, X_os, L_X_os, l_y_os, p_os, t, pm, tol
  REAL(KIND=8)                                     :: angle_shell, omega_shell, omega_shaft
  REAL(KIND=8)                                     :: dis_h_shell, dis_v_shell, vel_h_shell
  REAL(KIND=8)                                     :: vel_v_shell, dis_h_shaft, dis_v_shaft
  REAL(KIND=8)                                     :: vel_h_shaft, vel_v_shaft
  REAL(KIND=8)                                     :: tilt_h_shell, tilt_v_shell, tilt_dot_h_shell
  REAL(KIND=8)                                     :: tilt_dot_v_shell, tilt_h_shaft, tilt_v_shaft
  REAL(KIND=8)                                     :: tilt_dot_h_shaft, tilt_dot_v_shaft
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE            :: mu_vec, ac_vec  
  
  ! output variables of FVM_ELROD, e.g., simulation results
  INTEGER                                          :: convergent, iter, iter_sol
  REAL(KIND=8)                                     :: F_h, F_v, M_h, M_v, M_fr, V_oil, V_dot_bb, p_ref
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE          :: Pi_mat
  
  ! input/output variable used by the program FVM_ELROD to communicate with itself accross time steps
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE            :: pts_vec
  
  ! variables not required for using FVM_ELROD (only required by the example script at hand)
  INTEGER                                          :: i
  REAL(KIND=8)                                     :: pi
  INTEGER,DIMENSION(:,:),ALLOCATABLE               :: g_mat
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE          :: theta_mat, p_mat
  
  
  ! --------------------------------------------------------------------------------------------------------
  ! Main script
  ! --------------------------------------------------------------------------------------------------------
  
  ! First things first, let's define pi.
  
  pi = 3.14159265359d0                                                                                      ! define pi [-]
  
  ! Now, let's define the parameters that determine the number of DOFs. These will typically stay 
  ! constant throughout a sequence of calls of the Reynolds equation (for example, throughout a time 
  ! integration).
  
  d_b = 0.1d0                                                                                               ! bearing diameter [m]
  l_b = 0.08d0                                                                                              ! bearing length [m]
  n_x = 100                                                                                                 ! circumferential number of nodes  (not counting the periodic node at X=2*pi) [-]
  n_y = ANINT((l_b/(pi*d_b))*n_x+1)                                                                         ! axial number of nodes across the entire bearing length (including the nodes at the two bearing boundaries); we can set n_y = ANINT((l_b/(pi*d_b))*n_x+1) to ensure an equally fine discretization in both directions, but that requires l_b and d_b to be defined first [-]
  symBC = 1                                                                                                 ! flag for the usage of a symmetric boundary condition: 0 = no, 1 = yes; setting symBC = 1 reduces the computational effort but also renders the model incapable of considering shaft tilting [-]
  
  ! As the total number of nodes (n_x*n_y) is known now, we can allocate pts_vec. In a time integration, 
  ! the program FVM_ELROD uses this array to communicate with itself across time steps. At every call, 
  ! FVM_ELROD reads the data contained in pts_vec and stores new data in pts_vec. When calling FVM_ELROD, 
  ! pts_vec must contain the values that were stored in pts_vec by FVM_ELROD at the previous valid time 
  ! step. At the first time step, we use pts_vec to prescribe an initial condition instead. Setting 
  ! pts_vec = 0.0d0 and then pts_vec(n_x*n_y+1) = -1.0d-5 will lead to the following initial condition: 
  ! a film fraction of 1 (i.e. only oil, no gas cavities) is assumed at all nodes at t=-1.0d-5 (exactly at 
  ! t=0 isn't allowed, so we define it shortly before that). If Guembel conditions or quasistatic 
  ! conditions are chosen, FVM_ELROD only uses the values provided by pts_vec to choose an initial guess 
  ! for the iterative solver of the linear/linarized system of equations and, in case of quasistatic Elrod 
  ! simulations, also to determine the initial cavitation states. For transient Elrod simulations, these 
  ! initial guesses are also used, but more importantly, pts_vec allows FVM_ELROD to formulate the 
  ! transient cavitation term (which requires knowing the film fractions of the previous time step).
  
  ALLOCATE(pts_vec(n_x*n_y+1))
  pts_vec = 0.0d0
  pts_vec(n_x*n_y+1) = -1.0d-5
  
  ! Also, we can now allocate the arrays that will represent the nodal viscosities mu_vec and the 
  ! additional contour ac_vec (these are values which are added to the nodal gap widths to model 
  ! bearings where the surfaces aren't cylindrical), which are input variables for the FVM solution. 
  ! We also allocate a 2D-array Pi_mat which the FVM algorithm will use for output of the solution 
  ! (the computed pressure-like function) in a form that can be used to create surface plots. The 
  ! arrays mu_vec and ac_vec assume the node numbering scheme illustrated above.
  
  ALLOCATE(mu_vec(n_x*n_y))                                                                                    
  ALLOCATE(ac_vec(n_x*n_y))                                                                                    
  ALLOCATE(Pi_mat(n_x,n_y))
  
  ! Now, let's define some parameters and call the FVM algorithm for solving the Reynolds equation
  ! (check the attached PDF for clarification of the kinematic variables).
  
  c = 0.00015d0                                                                                             ! radial clearance [m]
  grooves = 1                                                                                               ! number of oil supply grooves (evenly distributed across the circumference), under Elrod conditions, there must be at least 1 groove [-]
  X_os = 0.25d0*2.0d0*pi                                                                                    ! angular circumferential position of one of the oil supply grooves in the reference frame of the shell [rad]
  L_X_os = 15.0d0/360.0d0*2.0d0*pi                                                                          ! angular circumferential side length of the oil supply groove(s) [rad]
  l_y_os = 0.06d0                                                                                           ! axial side length of the oil supply groove [m]
  p_os = 70000.0d0                                                                                          ! boundary value prescribed in the oil supply groove; either a non-negative pressure in [Pa] (zero corresponds to atmospheric pressure) or a negative value stating a film fraction minus 1 [-] (for example, a film fraction of 0.7 is prescribed by p_os = 0.7-1 = -0.3)
  t = 0.0d0                                                                                                 ! time corresponding to the current time step (irrelevant for this quasistatic single call) [s]
  angle_shell = 0.0d0                                                                                       ! rotation angle of the shell in the reference frame of the inertial system [rad]
  omega_shell = 0.0d0                                                                                       ! angular velocity of the shaft in the reference frame of the inertial system [rad/s] (omega_shaft must not be equal to omega_shell)
  dis_h_shell = 0.0d0                                                                                       ! horizontal displacement of the shell in the reference frame of the inertial system [m]
  dis_v_shell = 0.0d0                                                                                       ! vertical displacement of the shell in the reference frame of the inertial system [m]
  vel_h_shell = 0.0d0                                                                                       ! horizontal velocity of the shell in the reference frame of the inertial system [m/s]
  vel_v_shell = 0.0d0                                                                                       ! vertical velocity of the shell in the reference frame of the inertial system [m/s]
  omega_shaft = (3000/60)*2*pi                                                                              ! angular velocity of the shaft in the reference frame of the inertial system [rad/s] (omega_shaft must not be equal to omega_shell)
  dis_h_shaft = c/2                                                                                         ! horizontal displacement of the shaft in the reference frame of the inertial system [m]
  dis_v_shaft = 0.0d0                                                                                       ! vertical displacement of the shaft in the reference frame of the inertial system [m]
  vel_h_shaft = 0.0d0                                                                                       ! horizontal velocity of the shaft in the reference frame of the inertial system [m/s]
  vel_v_shaft = 0.0d0                                                                                       ! vertical velocity of the shaft in the reference frame of the inertial system [m/s]
  tilt_h_shell = 0.0d0                                                                                      ! tilting angle of shell around horizontal axis [rad]
  tilt_v_shell = 0.0d0                                                                                      ! tilting angle of shell around vertical axis [rad]
  tilt_dot_h_shell = 0.0d0                                                                                  ! rate of change of tilting angle of shell around horizontal axis [rad/s]
  tilt_dot_v_shell = 0.0d0                                                                                  ! rate of change of tilting angle of shell around vertical axis [rad/s]
  tilt_h_shaft = 0.0d0                                                                                      ! tilting angle of shaft around horizontal axis [rad]
  tilt_v_shaft = 0.0d0                                                                                      ! tilting angle of shaft around vertical axis [rad]
  tilt_dot_h_shaft = 0.0d0                                                                                  ! rate of change of tilting angle of shaft around horizontal axis [rad/s]
  tilt_dot_v_shaft = 0.0d0                                                                                  ! rate of change of tilting angle of shaft around vertical axis [rad/s]
  iter_max = n_x                                                                                            ! max. allowed number of iterations [-]
  ac_vec = 0.0d0                                                                                            ! additional contour of the shell defined at the nodes, following the node numbering scheme illustrated above; positive values increase the gap width, negative values reduce the gap width, the shell is cylindrical if all entries are zero; since the nodes are fixed at the reference frame of the shell, a rotation of the shell does not require any adjustment of ac_vec [m] 
  mu_vec = 0.01d0                                                                                           ! oil viscosities prescribed at the nodes, following the node numbering scheme illustrated above; since the nodes are fixed at the possibly rotating shell, this viscosity distribution is always expressed in the shell's reference frame [Pa*s]
  quasistatic = 1                                                                                           ! quasistatic = 1 means that a quasistatic simulation is requested (all time derivatives in the Reynolds equation are set to zero); when incorporating the SBFEM solution into a time integration scheme, set quasistatic = 0 [-]
  tol = 1.0d-8                                                                                              ! tolerance for BiCGStab (i.e. for iterative solution of equation system) [-]
  iter_max_solver = 10000                                                                                   ! maximum number of iterations for solution of system of equations (if an iterative solver is used) [-]
  pm = 100.0d0                                                                                              ! determines the penalty factor for the boundary conditions: pf=pm*max(diag(K)), where K is the system matrix [-]
  guembel = 0                                                                                               ! use Guembel instead of Elrod cavitation? 0 = no (i.e. use Elrod), 1 = yes (i.e. use Guembel)
  
  CALL FVM_ELROD(d_b, l_b, c, grooves, X_os, L_X_os, l_y_os, p_os, ac_vec, t, angle_shell, &
    omega_shell, omega_shaft, dis_h_shell, dis_v_shell, vel_h_shell, vel_v_shell, dis_h_shaft, &
    dis_v_shaft, vel_h_shaft, vel_v_shaft, tilt_h_shell, tilt_v_shell, tilt_dot_h_shell, &
    tilt_dot_v_shell, tilt_h_shaft, tilt_v_shaft, tilt_dot_h_shaft, tilt_dot_v_shaft, n_x, &
    iter_max, mu_vec, guembel, n_y, quasistatic, symBC, F_h, F_v, M_h, M_v, M_fr, V_oil, &
    V_dot_bb, Pi_mat, p_ref, convergent, iter, iter_sol, pts_vec, tol, iter_max_solver, pm)
  
  ! Now we have solved the Reynolds equation. To check the results, we can print some information.
  
  PRINT *, 'iterations = '
  PRINT *, iter
  PRINT *, 'hydrodynamic forces F_h and F_v (N) = '
  PRINT *, (/F_h,F_v/)                                                                                      ! hydrodynamic forces acting on the shell (and, with opposite sign, on the shaft) in the horizontal and vertical directions of the inertial system; check the attached PDF for clarification of the directions
  PRINT *, 'hydrodynamic moments M_h, M_v, and M_fr (Nm) = '
  PRINT *, (/M_h,M_v,M_fr/)                                                                                 ! hydrodynamic moments due to shaft tilting and hydrodynamic friction moment acting on the shell (and, with opposite sign, on the shaft); check the attached PDF for clarification of the directions
  PRINT *, 'oil volume (m^3) and flow through bearing boundaries (m^3/s) = '
  PRINT *, (/V_oil,V_dot_bb/)                                                                               ! oil volume in the bearing (considering cavitation) and oil volume flow through the bearing boundaries (this volume flow should be negative, indicating that oil is leaving the bearing)
  
  ! If we want to evaluate the nodal cavitation states, the nodal pressures, or the nodal film fractions
  ! accoring to the solution of the Reynolds equation that we performed, we can derive that information
  ! from the from the pressure-like function Pi stored in the 2D-array Pi_mat. We can define a switch
  ! function g indicating the cavitation states (g = 0 for Pi < 0 --> cavitation zone, g = 1 for 
  ! Pi >= 0 --> pressure zone). We can then compute the film fraction theta with theta = (1-g)*Pi+1 and
  ! the pressure p with p = g*Pi*p_ref. Here, p_ref is a reference pressure which the FVM program
  ! provides as output. Note that p = 0 corresponds to atmospheric pressure. I wrote the code below in
  ! order to export these two-dimensional fields (and then visualize them using the MATLAB script
  ! matlab_plot_results_2.m) for debugging and verification purposes. The file format used here for these
  ! exports (text files) is of course relatively expensive; if exporting the results is still desired in
  ! the long run, switching to binary file format is recommended.
  
  ALLOCATE(g_mat(n_x,n_y))
  ALLOCATE(theta_mat(n_x,n_y))
  ALLOCATE(p_mat(n_x,n_y))
  g_mat = 1
  WHERE(Pi_mat .LT. 0.0d0) g_mat = 0
  theta_mat = (1-g_mat)*Pi_mat+1
  p_mat = g_mat*Pi_mat*p_ref
  OPEN (UNIT = 2, FILE = "results2D_g.txt", RECL=(n_y*24))                                                  ! for 2D distribution of the the switch function: g
  OPEN (UNIT = 3, FILE = "results2D_theta.txt", RECL=(n_y*24))                                              ! for 2D distribution of the film fraction: theta
  OPEN (UNIT = 4, FILE = "results2D_Pi.txt", RECL=(n_y*24))                                                 ! for 2D distribution of the pressure-like function: Pi
  OPEN (UNIT = 5, FILE = "results2D_p.txt", RECL=(n_y*24))                                                  ! for 2D distribution of the pressure: p
  DO i = 1, n_x
    WRITE(2,*) g_mat(i,:)
    WRITE(3,*) theta_mat(i,:)
    WRITE(4,*) Pi_mat(i,:)
    WRITE(5,*) p_mat(i,:)
  END DO
  CLOSE(2)
  CLOSE(3)
  CLOSE(4)
  CLOSE(5)
  DEALLOCATE(g_mat)
  DEALLOCATE(theta_mat)
  DEALLOCATE(p_mat)
  
  ! Let's deallocate the remaining allocatable arrays
  
  DEALLOCATE(mu_vec)
  DEALLOCATE(ac_vec)
  DEALLOCATE(pts_vec)
  DEALLOCATE(Pi_mat)
  
END PROGRAM DEMONSTRATION_FVM_ELROD
