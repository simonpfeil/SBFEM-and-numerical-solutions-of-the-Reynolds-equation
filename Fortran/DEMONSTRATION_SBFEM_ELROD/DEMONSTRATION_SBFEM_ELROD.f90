! DEMONSTRATION_SBFEM_ELROD.f90 
!
! FUNCTIONS:
! DEMONSTRATION_SBFEM_ELROD
!
! **********************************************************************************************************
!
! PROGRAM: DEMONSTRATION_SBFEM_ELROD
!
! PURPOSE: This script demonstrates how to use the attached program SBFEM_ELROD (contained in the file 
! SBFEM_ELROD.f90 for solving the Reynolds equation semi-analytically. For clarification of the input and 
! output variables, the solution technique, the assumptions, and the coordinate system, check the 
! information provided in SBFEM_ELROD.f90 and the attached PDF.
!
! AUTHOR, AFFILIATION, DATE: Simon Pfeil, OvGU Magdeburg (Germany), 30.06.2025
!
! **********************************************************************************************************

PROGRAM DEMONSTRATION_SBFEM_ELROD
  
  ! no implicit variables
  IMPLICIT NONE
  
  ! input variables for SBFEM_ELROD: bearing properties, kinematic variables, discretization, etc.
  INTEGER										:: n_x, n_y, grooves, iter_max, quasistatic
  INTEGER										:: guembel
  REAL(KIND=8) 									:: d_b, l_b, c, X_os, L_X_os, p_os, t, angle_shell
  REAL(KIND=8) 									:: omega_shell, dis_h_shell, dis_v_shell, vel_h_shell
  REAL(KIND=8) 									:: vel_v_shell, omega_shaft, dis_h_shaft, dis_v_shaft
  REAL(KIND=8) 									:: vel_h_shaft, vel_v_shaft
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE			:: ac_vec, mu_vec
  
  ! more input variables for SBFEM_ELROD: data for Taylor approximations
  INTEGER										:: tay, n_tay, n_Ld, red
  REAL(KIND=8)                                  :: eps_constr, eps_max
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE			:: Ld_vec
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE       :: Vd_mat
  
  ! output variables of SBFEM_ELROD, e.g., simulation results
  INTEGER										:: convergent, iter
  INTEGER,DIMENSION(:),ALLOCATABLE				:: g_vec
  REAL(KIND=8) 									:: F_h, F_v, M_fr, V_oil, V_dot_bb, p_ref
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE			:: theta_vec
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE		:: Pi_mat
  
  ! input/output variable used by the program SBFEM_ELROD to communicate with itself accross time steps
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE			:: pts_vec
  
  ! variables used by the code at hand but not by SBFEM_ELROD
  INTEGER										:: i, n_points, n_Ld_all, ind
  REAL(KIND=8)                                  :: pi, gamma_ref, epsil
  REAL(KIND=8),DIMENSION(6)			            :: param_vec
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE			:: constr_vec, switch_vec, red_vec, ind_vec
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE         :: Ld_allpoints_vec
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE       :: Vd_allpoints_mat
  
  ! Let's first define the parameters that should stay constant throughout all calls of the Reynolds 
  ! equation within a time integration. Although we are not performing a time integration in the exemplary 
  ! code at hand, I prefer already making this distinction between the parameters here.
  
  pi = 3.14159265359d0																						! define pi [-]
  d_b = 0.1d0																								! bearing diameter [m]
  l_b = 0.08d0																								! bearing length [m]
  grooves = 1																								! number of oil supply grooves; the grooves are assumed to be distributed equidistantly along the circumference; at least 1 groove must be present under Elrod conditions (guembel = 0), and at least 0 grooves must be present under Guembel conditions (guembel = 1) [-]
  X_os = 0.25d0*2.0d0*pi																					! angular circumferential position of one of the oil supply grooves in the reference frame of the shell [rad]
  L_X_os = 15.0d0/360.0d0*2.0d0*pi																			! angular circumferential width of the oil supply groove(s) [rad]
  quasistatic = 1																							! quasistatic = 1 means that a quasistatic simulation is requested (all velocities and angular velocities except for omega_shaft are internally set to zero and the film fraction is assumed to be constant over time); when incorporating the SBFEM solution into a time integration scheme, set quasistatic = 0 [-]
  guembel = 0																								! use Guembel instead of Elrod cavitation? 0 = no (i.e. use Elrod), 1 = yes (i.e. use Guembel)
  n_x = 100																									! circumferential number of nodes  (not counting the periodic node at X=2*pi) [-]
  n_y = 26 ! 1																									! axial number of data points for evaluation of the solution on a twodimensional grid (purely for visualization purposes); set n_y=1 to skip this optional postprocessing step and save a lot of computational effort [-]  
  iter_max = 100																							! max. allowed number of iterations [-]
  tay = 0																									! use Taylor approximations for eigenvalues and eigenvectors? 0 = no, 1 = yes; these approximations are only possible if the additional contour ac_vec is zero, if the viscosity distribution mu_vec is spatially constant, and if n_x = 100, grooves = 0, and guembel = 1, which means that choosing tay = 1 will require us to align our parameters with these assumptions [-]
  
  ! If we want the SBFEM algorithm to solve its eigenvalue problem by Taylor approximations (tay = 1) 
  ! instead of using an eigensolver, we need to load the database of Taylor coefficients. The Taylor 
  ! coefficients in the given database were computed under certain assumtions, including Guembel
  ! conditions, no oil supply grooves, a pre-defined circumferential number of nodes (n_x = 100), and
  ! a pre-defined order of the Taylor series (n_tay = 2). In the following, if tay = 1, the database
  ! is loaded and the parameters are adjusted to these assumptions.
  
  IF ( tay .EQ. 1 ) THEN                                                                                    ! if Taylor approximations will be used for solving the eigenvalue problem
    OPEN(UNIT=1,FILE='tay_parameters.bin',STATUS='OLD',ACCESS='STREAM',FORM='UNFORMATTED')                  ! open the file that states the parameters which were chosen when computing the Taylor coefficients and ...
    READ(UNIT=1) param_vec                                                                                  ! ... copy the data to the array param_vec. Below, this data is then extracted from the array; this data includes ...
    n_x = INT(param_vec(1))                                                                                 ! ... the circumferential number of nodes (which may overwrite the one chosen above), ...
    n_tay = INT(param_vec(2))                                                                               ! ... the order of the Taylor series, ...
    gamma_ref = param_vec(3)                                                                                ! ... the original slenderness ratio (required for adjusting the eigenvalues to the current slenderness ratio, which is achived via a simple scaling operation), ...
    n_points = INT(param_vec(4))                                                                            ! ... the number of points where Taylor series were constructed, ...
    n_Ld_all = INT(param_vec(5))                                                                            ! ... the size of the vector containing the eigenvalue derivatives (i.e. containing the Taylor coefficients for computing the eigenvalues), and ...
    eps_max = param_vec(6)                                                                                  ! ... the maximum relative eccentricity for which Taylor approximations are available in the given database
    guembel = 1                                                                                             ! the Taylor approximations require Guembel conditions, so even if Elrod has been chosen (guembel = 0), we resort to assuming Guembel conditions now
    grooves = 0                                                                                             ! the Taylor approximations require a bearing without oil supply grooves, so if a nonzero number of grooves has been chosen, we change it back to zero
    ALLOCATE(Ld_allpoints_vec(n_Ld_all))                                                                    ! this array will contain the eigenvalue derivatives (i.e. the Taylor coefficients for computing the eigenvalues)
    ALLOCATE(Vd_allpoints_mat(n_x,n_Ld_all))                                                                ! this 2D array will contain the eigenvector derivatives (i.e. the Taylor coefficients for computing the eigenvectors)
    ALLOCATE(constr_vec(n_points))                                                                          ! this array will state the points (i.e. the relative eccentricities) where Taylor series were constructed
    ALLOCATE(switch_vec(n_points))                                                                          ! this array will state the points (i.e. the relative eccentricities) where the algorithm is supposed to switch between two neighboring Taylor series
    ALLOCATE(red_vec(n_points))                                                                             ! this array will state, for every Taylor series, the number of considered modes of the pressure field
    ALLOCATE(ind_vec(n_points))                                                                             ! this array will state, for every Taylor series, where within Ld_allpoints_vec and Vd_allpoints_mat we can find the corresponding Taylor coefficients
    CLOSE(1)                                                                                                ! close file opened above; below, the other database files are opened and their content is copied into the arrays allocated above ...
    OPEN (UNIT=1,FILE='tay_constrpoints.bin',STATUS='OLD',ACCESS='STREAM',FORM='UNFORMATTED')               ! ...
    READ(UNIT=1) constr_vec
    CLOSE(1)
    OPEN (UNIT=1,FILE='tay_switchpoints.bin',STATUS='OLD',ACCESS='STREAM',FORM='UNFORMATTED')
    READ(UNIT=1) switch_vec
    CLOSE(1)
    OPEN (UNIT=1,FILE='tay_reductions.bin',STATUS='OLD',ACCESS='STREAM',FORM='UNFORMATTED')
    READ(UNIT=1) red_vec
    CLOSE(1)
    OPEN (UNIT=1,FILE='tay_indices.bin',STATUS='OLD',ACCESS='STREAM',FORM='UNFORMATTED')
    READ(UNIT=1) ind_vec
    CLOSE(1)
    OPEN (UNIT=1,FILE='tay_eigenvalues.bin',STATUS='OLD',ACCESS='STREAM',FORM='UNFORMATTED')
    READ(UNIT=1) Ld_allpoints_vec
    CLOSE(1)
    OPEN (UNIT=1,FILE='tay_eigenvectors.bin',STATUS='OLD',ACCESS='STREAM',FORM='UNFORMATTED')
    READ(UNIT=1) Vd_allpoints_mat
    CLOSE(1)
    Ld_allpoints_vec = Ld_allpoints_vec*(((l_b/d_b)/gamma_ref)**2)                                          ! scale eigenvalues and their derivatives according to the defined bearing dimensions
  ELSE
    n_tay = 0                                                                                               ! no Taylor approximations are used, so the order of Taylor series will be ignored, we may set it to zero (or to any other number)
    eps_max = 0.0d0                                                                                         ! no Taylor approximations are used, so the maximum eccentricity for which Taylor approximations are available will be ignored, we may set it to zero (or to any other number)
  END IF
  
  ! Since n_x and n_y are defined by now, we can perform some array allocations that require these 
  ! parameters. The concrete values contained in these arrays will be defined later and, in a time 
  ! integration, may change from one call of the Reynolds equation to another.
  
  ALLOCATE(pts_vec(n_x+1))                                                                                  ! in- and output variable of SBFEM_ELROD; the program uses this array to communicate with itself across neighboring time steps, which is necessary for formulating the transient cavitation term (rate of change of the film fraction)
  ALLOCATE(ac_vec(n_x))                                                                                     ! input variable of SBFEM_ELROD; additional contour of the shell
  ALLOCATE(mu_vec(n_x))                                                                                     ! input variable of SBFEM_ELROD; viscosity distribution
  ALLOCATE(g_vec(n_x))                                                                                      ! output variable of SBFEM_ELROD; switch function
  ALLOCATE(theta_vec(n_x))                                                                                  ! output variable of SBFEM_ELROD; film fraction
  ALLOCATE(Pi_mat(n_x,n_y))                                                                                 ! output variable of SBFEM_ELROD; pressure-like function
  
  ! If we were to perform a time integration, the time integration would start here. The variables
  ! defined in the following may vary between calls of the Reynolds equation.
  
  ac_vec = 0.0d0																							! additional contour of the shell defined point-wise at the circumferential positions where the nodes are located; positive values increase the gap width, negative values reduce the gap width, the shell is cylindrical if all entries are zero; since the nodes are fixed at the reference frame of the shell, a rotation of the shell does not require any adjustment of ac_vec [m]
  mu_vec = 0.01d0																							! oil viscosity distribution defined point-wise at the circumferential positions where the nodes are located; since the nodes are fixed at the reference frame of the possibly rotating shell, this viscosity distribution is always expressed in the shell's reference frame [Pa*s]
  t = 0.0d0																									! time corresponding to the current time step (irrelevant for this quasistatic single call) [s]
  c = 0.00015d0																								! radial clearance [m]
  p_os = 70000.0d0																							! boundary value prescribed in the oil supply groove; either a non-negative pressure in [Pa] (zero corresponds to atmospheric pressure) or a negative value stating a film fraction minus 1 [-] (for example, a film fraction of 0.7 is prescribed by p_os = 0.7-1 = -0.3)
  angle_shell = 0.0d0                                                      									! rotation angle of the shell in the reference frame of the inertial system [rad]
  omega_shell = 0.0d0                                                      									! angular velocity of the shaft in the reference frame of the inertial system [rad/s] (omega_shaft must not be equal to omega_shell)
  dis_h_shell = 0.0d0                                                      									! horizontal displacement of the shell in the reference frame of the inertial system [m]
  dis_v_shell = 0.0d0                                                      									! vertical displacement of the shell in the reference frame of the inertial system [m]
  vel_h_shell = 0.0d0                                                      									! horizontal velocity of the shell in the reference frame of the inertial system [m/s]
  vel_v_shell = 0.0d0                                                      									! vertical velocity of the shell in the reference frame of the inertial system [m/s]
  omega_shaft = (3000/60)*2*pi                                                                              ! angular velocity of the shaft in the reference frame of the inertial system [rad/s] (omega_shaft must not be equal to omega_shell)
  dis_h_shaft = c/2                                                        									! horizontal displacement of the shaft in the reference frame of the inertial system [m]
  dis_v_shaft = 0.0d0                                                      									! vertical displacement of the shaft in the reference frame of the inertial system [m]
  vel_h_shaft = 0.0d0                                                      									! horizontal velocity of the shaft in the reference frame of the inertial system [m/s]
  vel_v_shaft = 0.0d0                                                      									! vertical velocity of the shaft in the reference frame of the inertial system [m/s]
  
  ! The next variable, pts_vec, requires special attention. In a time integration, the program SBFEM_ELROD 
  ! uses this array to communicate with itself across time steps. At every call, SBFEM_ELROD reads the data 
  ! contained in pts_vec and stores new data in pts_vec. When calling SBFEM_ELROD, pts_vec must contain the 
  ! values that were stored in pts_vec by SBFEM_ELROD at the previous valid time step. At the first time 
  ! step, we use pts_vec to prescribe an initial condition instead. Setting pts_vec = 0.0d0 and then 
  ! pts_vec(n_x+1) = -1.0d-5 will lead to the following initial condition: a film fraction of 1 (i.e. only 
  ! oil, no gas cavities) is assumed at all nodes at t=-1.0d-5 (exactly at t=0 isn't allowed, so we define 
  ! it shortly before that). If Guembel conditions are assumed, SBFEM_ELROD ignores pts_vec. If a 
  ! quasistatic solution with Elrod cavitation is performed, SBFEM_ELROD uses pts_vec only to choose an 
  ! initial guess for the cavitation states. For transient simulations with Elrod, this initial guess is 
  ! also used, but more importantly, pts_vec allows fvm_elrod to formulate the transient cavitation term 
  ! (which requires knowing the film fractions of the previous time step).
  
  pts_vec = 0.0d0
  pts_vec(n_x+1) = -1.0d-5
  
  ! Now the parameters are defined, but there is still one thing to take care of before solving the 
  ! Reynolds equation. As mentioned above, the SBFEM algorithm offer the option to tackle its eigenvalue 
  ! problem by Taylor approximations instead of using an eigensolver. This requires setting the variable 
  ! tay = 1 (see above). If we have chosen this option, we now need to extract the Taylor coefficients 
  ! that are suitable for the current shaft position (for the current relative eccentricity) from the 
  ! database (this is done in the following). This step must be executed before every call of the Reynolds 
  ! equation. The Taylor coefficients in the given database were computed under certain assumtions, 
  ! including the simplifications that the additional contour ac_vec is zero and that the viscosity 
  ! distribution described by mu_vec is spatially constant. If tay = 1, these conditions must be enforced 
  ! (this is also done in the following).
  
  IF ( tay .EQ. 1 ) THEN                                                                                    ! if Taylor approximations will be used for solving the eigenvalue problem
    epsil = SQRT((dis_h_shaft-dis_h_shell)**2+(dis_v_shaft-dis_v_shell)**2)/c								! relative eccentricity
    i = MINLOC(ABS(switch_vec-epsil),1)																		! switch_vec contains the epsilon-positions where we switch from one Taylor series to another; this command line determines which of these switching points our current epsilon is closest to
    IF ( epsil .LT. switch_vec(i) ) THEN																	! depending on whether we are above or below this switching point, ...
      i = i-1                                                                                               ! ... we choose the proper Taylor series (the i-th Taylor series in our database)
    END IF
    eps_constr = constr_vec(i)                                                                              ! point where Taylor series was constructed
    red = red_vec(i)                                                                                        ! number of modes to consider
    n_Ld = red*(n_tay+1)                                                                                    ! number of entries in the array of eigenvalue derivatives Ld_vec that we will define for the current Taylor series
    ind = ind_vec(i)                                                                                        ! ind will help us find the requested Taylor coefficients in the arrays Vd_mat and Ld_vec
    ALLOCATE(Ld_vec(n_Ld))                                                                                  ! when calling the SBFEM solution, this array should contain the eigenvalue derivatives for the concrete Taylor series to be used at this call (or contain nothing, if tay = 0)
    ALLOCATE(Vd_mat(n_x,n_Ld))                                                                              ! when calling the SBFEM solution, this 2D array should contain the eigenvector derivatives for the concrete Taylor series to be used at this call (or contain nothing, if tay = 0)
    Ld_vec = Ld_allpoints_vec(ind:(ind-1+n_Ld))                                                             ! vector of eigenvalue derivatives (Taylor coefficients for approximating the eigenvalues) for the concrete Taylor series that will be used at the next call of the SBFEM algorithm
    Vd_mat = Vd_allpoints_mat(:,ind:(ind-1+n_Ld))                                                           ! matrix of eigenvector derivatives (Taylor coefficients for approximating the eigenvectors) for the concrete Taylor series that will be used at the next call of the SBFEM algorithm
    ac_vec = 0.0d0                                                                                          ! the Taylor approximations assume that there is no additional contour, so we set ac_vec to zero
    mu_vec = (SUM(mu_vec)/SIZE(mu_vec))                                                                     ! the Taylor approximations assume that the viscosity distribution described by mu_vec is spatially constant, so we replace the viscosities by their average
  ELSE
    eps_constr = 0.0d0                                                                                      ! no Taylor approximations are used, so the point of construction will be ignored, we may set it to zero (or to any other number)
    red = 0                                                                                                 ! no Taylor approximations are used, so the number of modes will be ignored, we may set it to zero (or to any other number)
    n_Ld = 0                                                                                                ! no Taylor approximations are used, so let's set n_Ld to zero, preventing the arrays Ld_vec and Vd_mat from occupying memory they don't require
    ALLOCATE(Ld_vec(n_Ld))                                                                                  ! when calling the SBFEM solution, this array should contain the eigenvalue derivatives for the concrete Taylor series to be used at this call (or contain nothing, if tay = 0)
    ALLOCATE(Vd_mat(n_x,n_Ld))                                                                              ! when calling the SBFEM solution, this 2D array should contain the eigenvector derivatives for the concrete Taylor series to be used at this call (or contain nothing, if tay = 0)
  END IF
  
  ! We call the SBFEM algorithm to solve the Reynolds equation.
  
  CALL SBFEM_ELROD(grooves, n_x, iter_max, quasistatic, guembel, n_y, d_b, l_b, c, X_os, &
	L_X_os, p_os, t, angle_shell, omega_shell, dis_h_shell, dis_v_shell, vel_h_shell, vel_v_shell, &
	omega_shaft, dis_h_shaft, dis_v_shaft, vel_h_shaft, vel_v_shaft, ac_vec, mu_vec, tay, n_tay, &
    n_Ld, red, eps_constr, eps_max, Ld_vec, Vd_mat, convergent, iter, M_fr, V_oil, V_dot_bb, p_ref, &
    F_h, F_v, g_vec, theta_vec, Pi_mat, pts_vec)
  
  ! Since the arrays Ld_vec and Vd_mat are allocated before every call of SBFEM_ELROD (possibly with a 
  ! different size each time), we should deallocate them after every call.
  
  DEALLOCATE(Ld_vec)
  DEALLOCATE(Vd_mat)
  
  ! Let's print some results to the screen.
  
  PRINT *, 'iterations = '
  PRINT *, iter
  PRINT *, 'hydrodynamic forces F_h and F_v (N) = '
  PRINT *, (/F_h,F_v/)																						! hydrodynamic forces acting on the shell (and, in the opposite direction, on the shaft); check the attached PDF for clarification of the coordinate system 
  PRINT *, 'hydrodynamic friction moment M_fr (Nm) = '
  PRINT *, (/M_fr/) 																						! hydrodynamic friction moment acting on the shell (and, in the opposite direction, on the shaft); check the attached PDF for clarification of the coordinate system 
  PRINT *, 'oil volume (m^3) and flow through bearing boundaries (m^3/s) = '
  PRINT *, (/V_oil,V_dot_bb/)																				! oil volume in the bearing (considering cavitation) and oil volume flow through the bearing boundaries (this volume flow should be negative, indicating that oil is leaving the bearing)
  
  ! The variable n_y defines the number of data points in the axial direction for evaluating the pressure-
  ! like function Pi (a two-dimensional field) on a discrete grid and saving it to Pi_mat. I usually set
  ! n_y to 1, which prompts SBFEM_ELROD to skip this optional postprocessing step, saving a lot of 
  ! computational effort (note that the hydrodynamic forces etc. are not affected by n_y, as their
  ! computation relies on a continuous, analytical formulation in the axial direction). However, if we want
  ! to evaluate these 2D results for some reason, we must set n_y to a reasonably large value. For example,
  ! you can execute this script again with n_y equal to, e.g., 20 (you need to scroll upwards a little to
  ! find the line where n_y is defined). This will cause the pressure-like function to be evaluated at 20
  ! data points in the axial direction. From that, we can also derive the 2D pressure distribution. I wrote
  ! the code below to check whether these 2D results are computed correctly (i.e., for debugging). The
  ! pressure-like function and the pressure (2D results) are written to text files, the switch function
  ! and the film fraction (1D results) as well. We can later use the attached Matlab script 
  ! matlab_plot_results.m to visualize the results, if Matlab is available. Of course, writing text files 
  ! is expensive; if exporting the results to files is still desired in the long term, switching to binary 
  ! file format is recommended. 
    
  OPEN (UNIT = 2, FILE = "results1D_g.txt", RECL=(24))                                                      ! for 1D distribution of the the switch function
  OPEN (UNIT = 3, FILE = "results1D_theta.txt", RECL=(24))                                                  ! for 1D distribution of the film fraction
  OPEN (UNIT = 4, FILE = "results2D_Pi.txt", RECL=(n_y*24))													! for 2D distribution of the pressure-like function
  OPEN (UNIT = 5, FILE = "results2D_p.txt", RECL=(n_y*24))													! for 2D distribution of the pressure
  DO i = 1, n_x
	WRITE(2,*) g_vec(i)
	WRITE(3,*) theta_vec(i)
	WRITE(4,*) Pi_mat(i,:)
	WRITE(5,*) Pi_mat(i,:)*(g_vec(i)*p_ref)
  END DO
  CLOSE(2)
  CLOSE(3)
  CLOSE(4)
  CLOSE(5)
  
  ! Final step after our last solution (or only solution) of the Reynolds equation: deallocate arrays.
  
  DEALLOCATE(ac_vec)
  DEALLOCATE(mu_vec)
  DEALLOCATE(g_vec)
  DEALLOCATE(theta_vec)
  DEALLOCATE(Pi_mat)
  DEALLOCATE(pts_vec)
  IF ( tay .EQ. 1 ) THEN
    DEALLOCATE(constr_vec)
    DEALLOCATE(switch_vec)
    DEALLOCATE(red_vec)
    DEALLOCATE(ind_vec)
    DEALLOCATE(Ld_allpoints_vec)
    DEALLOCATE(Vd_allpoints_mat)
  END IF
  
END PROGRAM DEMONSTRATION_SBFEM_ELROD

INCLUDE 'SBFEM_ELROD.f90'