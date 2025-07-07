function [F_h,F_v,M_fr,V_oil,V_dot_bb,pts_vec,g_vec,theta_vec,Pi_mat,...
    p_ref,convergent,iter] ...
    = sbfem_elrod(d_b,l_b,c,grooves,X_os,L_X_os,p_os,t,angle_shell,...
    omega_shell,dis_h_shell,dis_v_shell,vel_h_shell,vel_v_shell,...
    omega_shaft,dis_h_shaft,dis_v_shaft,vel_h_shaft,vel_v_shaft,...
    n_x,iter_max,ac_vec,mu_vec,pts_vec,guembel,quasistatic,n_y,...
    tay,n_tay,epsilon_max,epsilon_constr,red,Vd_mat,Ld_vec)

% PROGRAM FOR SBFEM SOLUTION OF THE REYNOLDS EQUATION WITH CAVITATION

% Author: Simon Pfeil
% Affiliation: OvGU Magdeburg, Germany (until 30.06.2025)
% Date: 30.06.2025

% --- comments ------------------------------------------------------------
%
% - This algorithm solves the Reynolds equation for hydroynamic journal 
%   bearings.
%
% - The computation is based on the semi-analytical "scaled boundary finite 
%   element method" (SBFEM).
%
% - A transient, nonlinear cavitation model is used (optionally, the simple
%   Guembel approach can be used instead).
%
% - Shaft tilting is not considered.
%
% - The oil supply groove is modeled with simplified geometry: its axial 
%   length is assumed to be equal to the bearing length.
%
% - Given a plain bearing without oil supply grooves and under Guembel
%   conditions, the efficiency of this program can be enhanced by means of
%   eigenvalue problem derivatives, meaning that an eigenvalue problem is
%   solved by Taylor approximations (instead of by an eigensolver). A set 
%   of Taylor coefficients (as well as some further data) that can be used
%   for this purpose is provided via the files 1_locations.mat,
%   2_reductions.mat, and 3_coefficients.mat, assuming that the relative
%   eccentricity does not exceed 0.99 (if it does, the algorithm will
%   automatically resort to using an eigensolver instead of the Taylor
%   approximation).
%
% - This program is designed to be incorporated into time integration 
%   schemes for rotordynamic (or multibody) simulations; however, a quasi-
%   static solution is also possible.
%
% - The input and output variables are explained below.
%
% - See the attached PDF for clarification of the coordinate systems and
%   node numbers. x (or the nondimensionalized X) and y (or the 
%   nondimensionalized Y) are the coordinates of the lubrication gap and 
%   are used by the Reynolds equation. This coordinate system, as well as 
%   the computational grid, are fixed in the reference frame of the shell. 
%   The transformations between the reference frame of the shell and the
%   inertial system are performed internally by this program. The kinematic 
%   variables of the shell and the shaft, when handed to this program, 
%   should be formulated from the perspective of the inertial system. The 
%   hydrodynamic forces computed by this program are automatically 
%   transformed into the inertial system for output. In contrast, those 
%   input and output variables that describe one- or two-dimensional fields 
%   (e.g., viscosity distribution, additional contour, and pressure-like 
%   function) don't use the inertial system, as they are expressed by 
%   arrays whoose entries represent the values at the nodes (which are 
%   always fixed at the shell). The circumferential position of the oil 
%   supply groove X_os is also always expressed in the reference frame of 
%   the shell; otherwise, a rotating shell would require us to keep 
%   updating this input variable.
%
% - More thorough discussions of the computational method, of the 
%   cavitation model, of the boundary conditions (BCs), and of further 
%   assumptions will be given in the thesis "Simulating Hydrodynamic 
%   Bearings with the Scaled Boundary Finite Element Method" by Simon 
%   Pfeil, but this thesis isn't publically available yet
%
% - Some steps of the algorithm are numbered; the assigned numbers 
%   loosely correlate to the overview given in the above mentioned thesis 
%   (the overview given for the case with nonlinear cavitation, not the
%   one that assumes Guembel conditions).
%
% - At some point in the algorithm, the node numbers are shifted depending 
%   on the position of the oil supply groove, which ensures a tridiagonal 
%   E2-matrix. This affects the order in which the nodal values of all 
%   discretized functions are stored. Since this operation is performed 
%   (and, at the end, reversed) internally, we don't need to worry about 
%   it when using this program. However, if we want to analyze the 
%   algorithm or check intermediate results, we need to keep this in mind.

% --- input variables -----------------------------------------------------
%
% - d_b = bearing diameter [m]
%
% - l_b = bearing length [m]
%
% - c = radial clearance (shell radius minus shaft radius) [m]
%
% - grooves = number of oil supply grooves; the grooves are assumed to be 
%   distributed equidistantly along the circumference; at least 1 groove 
%   must be present under Elrod conditions (guembel = 0), and at least 0
%   grooves must be present under Guembel conditions (guembel = 1) [-]
%
% - X_os = angular circumferential position of the center of one of the
%   oil supply grooves (must not be smaller than 0 or greater than 2*pi,
%   see attached PDF); since the coordinate X is fixed at the shell, a 
%   rotation of the shell does not change the value of X_os [rad]
%
% - L_X_os = angular circumferential side length of the oil supply groove
%   (see attached PDF file) [rad]
%
% - p_os (if positive or zero) = pressure prescribed in the oil supply 
%   groove (zero corresponds to atmospheric pressure) [Pa]
%
% - p_os (if negative) = film fraction prescribed in the oil supply groove
%   minus 1 (for example, a film fraction of 0.7 is prescribed by setting 
%   p_os = 0.7-1 = -0.3) [-]
%
% - t = physical time at the current time step [s]
%
% - angle_shell = rotation angle of the shell from the perspective of the
%   inertial system (see attached PDF); if the shell rotates, knowing this
%   angle is crucial in order for this program to perform the 
%   transformations between the inertial system and the reference frame of 
%   the shell [rad]
%
% - omega_shell, omega_shaft = angular velocities of the shell and of the 
%   shaft, respectively, around the longitudinal axis, measured by a non-
%   rotating observer [rad/s]
%
% - dis_h_shell, dis_h_shaft = horizontal displacements of the shell and 
%   the shaft, respectively, from the perspective of the inertial system 
%   (see attached PDF) [m]
%
% - dis_v_shell, dis_v_shaft = vertical displacements of the shell and the
%   shaft, respectively, from the perspective of the inertial system (see
%   attached PDF) [m]
%
% - vel_h_shell, vel_h_shaft = horizontal velocities of the shell and the
%   shaft, respectively, from the perspective of the inertial system (see
%   attached PDF) [m/s]
%
% - vel_v_shell, vel_v_shaft = vertical velocities of the shell and the
%   shaft, respectively, from the perspective of the inertial system (see
%   attached PDF) [m/s]
%
% - n_x = number of nodes for circumferential equidistant grid; the  
%   smallest node number 1 is located at the angular position X=0 and the  
%   largest node number n_x at X=2*pi*(1-1/n_x); the node at X=2*pi is not 
%   counted because the periodic BC merges this node with node 1 (see 
%   attached PDF file) [-]
%
% - iter_max = maximum allowed number of iterations for the cavitation 
%   model (usually, less than 10 iterations are enough, but in some cases,
%   more may be necessary; I usually choose iter_max = n_x) [-]
%
% - ac_vec = column vector of length n_x containing an additional contour 
%   along the circumferential (but not axial) direction in case the shell
%   is not cylindrical; the entries of this vector are simply added to the 
%   nodal gap widths; since the nodes are fixed at the shell, a rotation 
%   of the shell does not require this additional contour to be updated; 
%   note that an additional contour of the shaft is not allowed in the 
%   current version of this program [m]
%
% - mu_vec = column vector of length n_x containing the dynamic oil
%   viscosities at the nodes, from node 1 to node n_x; since the 
%   nodes are fixed at the possibly rotating shell, the viscosity 
%   distribution is always expressed in the shell's reference frame [Pa*s]
%
% - pts_vec = column vector of length n_x+1 where the first n_x entries are 
%   the nodal pressure-like functions of the previous time step (averaged
%   in the axial direction) and the entry n_x+1 states the time of the 
%   previous time step; this variable is used for transferring data 
%   between time steps, as required by the transient cavitation model; note 
%   that pts_vec is also an output variable; simply take the pts_vec that 
%   was generated as output at the previous time step and provide it as 
%   input variable at the current time step [-]/[s]
% 
% - quasistatic = a flag that specifies whether a quasi-static solution is
%   desired (1 = yes, 0 = no); if yes, then all velocities and angular 
%   velocities except for omega_shaft are set to zero and the film fraction 
%   is assumed to be constant over time [-]
%
% - n_y = number of points in the axial direction for evaluation of the
%   solution (pressure-like function) on a two-dimensional grid (optional
%   postprocessing step for creation of a surface plot); this step is 
%   skipped if n_y = 1 [-]
%
% - guembel = determines whether Guembel conditions should be assumed,
%   bypassing the Elrod cavitation model; guembel = 0 --> use Elrod, 
%   guembel = 1 --> use Guembel [-]
%
% - tay = determines whether the eigenvalue problem is solved by a Taylor
%   approximation (0 = no, 1 = yes); the provided database of Taylor
%   coefficients assumes that the additional contour ac_vec is zero, that 
%   the viscosity distribution described by mu_vec is spatially constant, 
%   and that n_x = 100, grooves = 0, and guembel = 1, so we must choose
%   our parameters accordingly when using tay = 1; the script
%   demonstration_sbfem_elrod.m explains how to use this database [-]
% 
% - n_tay = order of the Taylor approximation (i.e. polynomial degree of 
%   the Taylor series); only used in case tay = 1; the provided database 
%   of Taylor coefficients assumes that n_tay = 2 [-]
% 
% - epsilon_max = maximum relative eccentricity at which Taylor
%   approximations for solving the eigenvalue problem are available;
%   if the current relative eccentricity exceeds epsilon_max, the
%   eigenvalue problem is solved by an eigensolver even if tay = 1 [-]
% 
% - epsilon_constr = relative eccentricity where the Taylor series that
%   will be employed was constructed; only used if tay = 1 [-] 
% 
% - red = number of considered modes of the pressure field due to modal
%   reduction; only used in case tay = 1; if we use the provided database
%   of Taylor coefficients, this database also determines what number of 
%   modes to choose for a given eccentricity (see
%   demonstration_sbfem_elrod.m) [-]
% 
% - Vd_mat = matrix containing the eigenvector derivatives (coefficients 
%   of the Taylor series to be used); only used if tay = 1 (otherwise,
%   simply define Vd_mat as an empty variable); the script 
%   demonstration_sbfem_elrod.m demonstrates how to obtain Vd_mat from the
%   database provided together with this program [-]
%
% - Ld_vec = vector containing the eigenvalue derivatives (coefficients 
%   of the Taylor series to be used); only used if tay = 1 (otherwise,
%   simply define Ld_vec as an empty variable); the script 
%   demonstration_sbfem_elrod.m demonstrates how to obtain Ld_vec from the
%   database provided together with this program [-]

% --- output variables ----------------------------------------------------
%
% - F_h = hydrodynamic force acting on the shell in the horizontal 
%   direction of the inertial system (see attached PDF); the same force 
%   acts on the shaft in the opposite direction [N]
%
% - F_v = hydrodynamic force acting on the shell in the vertical
%   direction of the inertial system (see attached PDF); the same force 
%   acts on the shaft in the opposite direction [N]
%
% - M_fr = oil friction moment acting on the shell (see attached PDF); the 
%   same moment acts on the shaft in the opposite direction [Nm]
%
% - V_oil = total oil volume in the bearing (equal to the volume of the
%   lubrication gap minus the gas volume in the cavitation zone) [m^3]
%
% - V_dot_bb = oil volume flow through the bearing boundaries (a negative 
%   sign indicates a flow direction out of the bearing) [m^3/s]
%
% - pts_vec = column vector of length n_x+1 where the first n_x entries
%   are the axially averaged nodal pressure-like functions and the entry 
%   n_x+1 contains the current physical time [-]/[s]
%
% - g_vec = column vector of length n_x containing the switch functions g 
%   at the circumferential positions of the nodes; g = 1 is the pressure 
%   zone and g = 0 is the cavitation zone; when analyzing g_vec for a 
%   system where the shell rotates, we need to keep in mind that the nodes 
%   are fixed at the shell (and so are the entries of g_vec) [-]
%
% - theta_vec = column vector of length n_x containing the film fractions 
%   theta at the circumferential positions of the nodes; theta = 1 
%   indicates a fully oil-filled gap and theta < 1 indicates a cavitated 
%   fluid film; when analyzing theta_vec for a system where the shell 
%   rotates, we need to keep in mind that the nodes are fixed at the shell 
%   (and so are the entries of theta_vec) [-]
%
% - Pi_mat = matrix of size n_x*n_y stating the two-dimensional 
%   distribution of the computed pressure-like function Pi; from this, the 
%   pressure distribution can be derived (p = p_ref*Pi*g), where g is the
%   switch function (see above) and p_ref is the reference pressure which
%   was used for the nondimensionalization (see below); Pi_mat is only 
%   computed if n_y > 1; each row of Pi_mat represents a circumferential
%   position where a node is located (if the shell rotates, we need to 
%   keep in mind that the nodes are fixed at the shell), while each column 
%   represents an axial position [-]
%
% - p_ref = reference pressure due to nondimensionalization [Pa]
%
% - convergent = a flag that states whether or not the solution has 
%   converged (1 = yes, 0 = no) [-]
%
% - iter = number of iterations [-]



% ----------------------------------------------------------------------- %
% --- EXPRESS INPUTS W.R.T. COORDINATE SYSTEM FOR REYNOLDS {x,y} -------- %
% ----------------------------------------------------------------------- %


dis_h = dis_h_shaft-dis_h_shell;                                           % horizontal displacement of the shaft relative to the shell, still in the reference frame of the inertial system
dis_v = dis_v_shaft-dis_v_shell;                                           % vertical displacement of the shaft relative to the shell, still in the reference frame of the inertial system
q = max(sqrt(dis_h^2+dis_v^2),eps);                                        % absolute eccentricity
X_att = atan2(dis_v,dis_h);                                                % attitude angle in the reference frame of the inertial system
X_att = X_att-angle_shell;                                                 % attitude angle in the shell-fixed reference frame in which the Reynolds equation is solved

if quasistatic == 0
    omega = omega_shaft-omega_shell;                                       % rotational velocity of the shaft relative to the shell
    vel_h = vel_h_shaft-vel_h_shell;                                       % horizontal velocity of the shaft relative to the shell, still in the reference frame of the inertial system
    vel_v = vel_v_shaft-vel_v_shell;                                       % vertical velocity of the shaft relative to the shell, still in the reference frame of the inertial system
    q_dot = (vel_v*dis_v+vel_h*dis_h)/q;                                   % rate of change of absolute eccentricity
    X_att_dot = (vel_v*dis_h-vel_h*dis_v)/q^2;                             % rate of change of attitude angle in the reference frame of the inertial system
    X_att_dot = X_att_dot-omega_shell;                                     % rate of change of attitude angle in the shell-fixed reference frame in which the Reynolds equation is solved
else
    omega = omega_shaft;                                                   % quasistatic case: rotational velocity of shell is assumed to be zero
    q_dot = 0;                                                             % quasistatic case: rate of change of eccentricity is assumed to be zero
    X_att_dot = 0;                                                         % quasistatic case: rate of change of attitude angle is assumed to be zero
end

if tay == 1                                                                % if a Taylor-series solution of the eigenvalue problem is desired, we must shift the computational grid (and the coordinate system) so that node 1 is located at the position of the minimum gap; after this shift, the attitude angle is zero
    X_att_trafo = X_att;                                                   % since the attitude angle will be set to zero in the line below, we must now save the actual attitude angle to allow transforming the solution back at the end
    X_att = 0;                                                             % the attitude angle must be set to zero for now, but the correct attitude angle will be considered later via a simple transformation
end



% ----------------------------------------------------------------------- %
% --- NONDIMENSIONALIZATION --------------------------------------------- %
% ----------------------------------------------------------------------- %


u = omega*(d_b/2);                                                         % circumferential surface velocity of the shaft in the shell-fixed reference frame in which the Reynolds equation is solved (the surface velocity of the shell is zero in this reference frame) [m/s]
sgn_u = sign(u);                                                           % sign of circumferential surface velocity
mu_ref = mean(mu_vec);                                                     % reference viscosity (order of magnitude should be similar to the actual viscosities) [Pa*s]
p_ref = abs(u)*mu_ref*(d_b/2)/(2*c^2);                                     % reference pressure for nondimensionalization [Pa]
mu_rel_vec = mu_vec/mu_ref;                                                % dimensionless viscosities at all nodes
epsilon = q/c;                                                             % relative eccentricity [-]
epsilon_dot = q_dot/c;                                                     % rate of change of the relative eccentricity [1/s]



% ----------------------------------------------------------------------- %
% --- ANALIZE GRID ------------------------------------------------------ %
% ----------------------------------------------------------------------- %


L_X = 2*pi/n_x;                                                            % dimensionless circumferential side length of a sector
L_Y = l_b/d_b;                                                             % dimensionless axial side length of a sector (equal to slenderness ratio)

X_vec = linspace(0,2*pi-L_X,n_x)';                                         % vector containing the X-positions of all nodes

index_0_vec = linspace(1,n_x,n_x)';                                        % vector containing all node numbers
index_E_vec = vertcat(linspace(2,n_x,n_x-1)',1);                           % vector containing the number of the neighboring node E for every node 0
index_W_vec = vertcat(n_x,linspace(1,n_x-1,n_x-1)');                       % vector containing the number of the neighboring node W for every node 0



% ----------------------------------------------------------------------- %
% --- ANALYZE LOCATIONS OF THE OIL SUPPLY GROOVES ----------------------- %
% ----------------------------------------------------------------------- %


X_os_vec = linspace(0,grooves-1,grooves)'*(2*pi/grooves) + X_os;           % vector containing the angular circumferential positions of the centers of all oil supply grooves

start_os_vec = round((X_os_vec-L_X_os/2)/L_X+1);                    	   % vector containing the node numbers where the oil supply grooves begin (if they are outside the range 1...n_x, this will be corrected later)
end_os_vec = round((X_os_vec+L_X_os/2)/L_X+1);                      	   % vector containing the node numbers where the oil supply grooves end (if they are outside the range 1...n_x, this will be corrected later)

dof = n_x-sum(end_os_vec-start_os_vec+1);                                  % dof = number of degrees of freedom after incorporation of the oil supply BCs



% ----------------------------------------------------------------------- %
% --- SHIFT NODE NUMBERS TO ENSURE TRIDIAGONAL E2-MATRIX AFTER BCs ------ %
% ----------------------------------------------------------------------- %


if grooves >= 1
    start_os1 = start_os_vec(1,1);                                         % original node number at the onset of the first defined oil supply groove; this node will be defined as the new node number 1
else
    start_os1 = 1;                                                         % without oil supply BCs, there is no use in shifting the node numbers (a tridiagonal matrix is not achievable); setting start_os1 to 1 prevents this shift
end

start_os_vec = start_os_vec-start_os1+1;                                   % the node numbers given in the vector that describes where the oil supply grooves begin are shifted so that the first groove begins at number 1 ...
end_os_vec = end_os_vec-start_os1+1;                                       % ... and the node numbers given in the vector that describes where the oil supply grooves end undergo the same shift

start_os1 = mod(start_os1-1,n_x)+1;                                        % the original number of the node that will be defined as number 1, given by start_os1, will be required again below; this time, it is not allowed to be outside the range 1...n_x (if that is the case, start_os1 is now shifted by a multiple of n_x into this range)

shift_vec = vertcat(...                                                    % this vector describes how each node number is changed; for example, if the 12th vector entry is 20, the number 12 will be assigned to the node that originally had the number 20
    linspace(start_os1,n_x,n_x-start_os1+1)',...
    linspace(1,start_os1-1,start_os1-1)');

X_vec = X_vec(shift_vec,1);                                                % adjust the order of the vector entries to the new node numbers
mu_vec = mu_vec(shift_vec,1);                                              % adjust the order of the vector entries to the new node numbers
ac_vec = ac_vec(shift_vec,1);                                              % adjust the order of the vector entries to the new node numbers
pts_vec(1:n_x,1) = pts_vec(shift_vec,1);                                   % adjust the order of the vector entries to the new node numbers



% ----------------------------------------------------------------------- %
% --- SOME PREPARATIONS FOR THE OIL SUPPLY BCs -------------------------- %
% ----------------------------------------------------------------------- %


nodes_dof_vec = index_0_vec;                                               % see comment below

for k = 1:grooves
    nodes_dof_vec(start_os_vec(k,1):end_os_vec(k,1),1) = 0;                % this vector has one entry per node; an entry equal to 1 will indicate that the corresponding node belongs to the oil supply, the other entries remain 0
end

nodes_dof_vec = nonzeros(nodes_dof_vec);                                   % vector containing all node numbers that are not part of the oil supply; these are the DOFs that remain after incorporation of the oil supply BCs



% ----------------------------------------------------------------------- %
% --- ANALYZE BOUNDARY VALUE PRESCRIBED IN OIL SUPPLY GROOVE ------------ %
% ----------------------------------------------------------------------- %


if p_os >= 0                                                               % "IF a non-negative supply pressure is prescribed, then ..."
    Pi_os = p_os/p_ref;                                                    % nondimensionalize the supply pressure and save as pressure-like function
    g_os = 1;                                                              % set corresponding switch function to 1
else                                                                       % otherwise, the prescribed value p_os will be interpreted as a film fraction with an offset of -1, i.e., as a negative pressure-like function
    Pi_os = p_os;                                                          % save prescribed value as pressure-like function
    g_os = 0;                                                              % set corresponding switch function to 0
end

if ( guembel == 1 && p_os < 0 ) || grooves == 0                            % if Guembel is used instead of Elrod, we cannot prescribe a film fraction as BC (which would be done via a negative p_os), so p_os represents a pressure and can't be negative; moreover, we set p_os to zero in case no oil supply is modeled (as one out of several steps necessary to bypass this BC)
    Pi_os = 0;                                                             % set oil supply pressure to zero
end



% ----------------------------------------------------------------------- %
% --- READ DATA FROM PREVIOUS TIME STEP --------------------------------- %
% ----------------------------------------------------------------------- %


Pi_bar_pts_vec = pts_vec(1:n_x,1);                                         % vector of nodal pressure-like functions, averaged in the axial direction, at previous time step

if guembel == 0
    g_pts_vec = sign(sign(Pi_bar_pts_vec)+1);                              % vector of nodal switch functions at previous time step
else
    g_pts_vec = ones(n_x,1);                                               % Guembel --> set swtich function to 1
end

g_vec = ones(n_x,1)*g_os;                                                  % initialize switch function, step 1: make it consistent with supply BC
g_vec(nodes_dof_vec,1) = g_pts_vec(nodes_dof_vec,1);                       % initialize switch function, step 2: use switch function from previous time step as initial guess at all nodes without oil supply BCs
g_vec_old = [];                                                            % switch function from previous iteration (does not exist yet)

if quasistatic == 1                                                        % --> a quasi-static solution is desired
    Delta_T = abs(u)/(2*(d_b/2))*1e16;                                     % set time increment to a very large number for quasi-static solutions, suppressing the time derivative of the film fraction
else                                                                       % --> a transient solution is desired, which takes into account the results (concretely, the film fraction) of a previous time step
    t_pts = pts_vec(n_x+1,1);                                              % time at previous time step
    Delta_T = abs(u)/(2*(d_b/2))*(t-t_pts);                                % dimensionless time increment
end



% ----------------------------------------------------------------------- %
% --- ANALYZE GAP FUNCTION AND ITS RATE OF CHANGE ----------------------- %
% ----------------------------------------------------------------------- %


H_vec = 1-epsilon*cos(X_vec-X_att)+ac_vec/c;                               % nondimensionalized gap function at every node
H3_vec = H_vec.^3;                                                         % cubed nondimensionalized gap function at every node
dHdT_vec = (d_b/abs(u))*(-epsilon_dot*cos(X_vec-X_att)-...                 % nondimensionalized time derivative of the nondimensionalized gap function
    epsilon*X_att_dot*sin(X_vec-X_att));



% ----------------------------------------------------------------------- %
% --- DEFINE SOME VECTORS IN PREPARATION OF THE ASSEMBLY ---------------- %
% ----------------------------------------------------------------------- %


H_E_vec = H_vec(index_E_vec,1);                                            % vector containing the dimensionless gap function at neighboring node E for each node 0
H_W_vec = H_vec(index_W_vec,1);                                            % vector containing the dimensionless gap function at neighboring node W for each node 0
H3_over_mu_rel_vec = H3_vec./mu_rel_vec;                                   % vector containing the cubed gap function divided by the relative viscosity at each node
b_E_vec = (1/(24*L_X))*...                                                 % vector containing the b_E-coefficients
    (H3_over_mu_rel_vec+H3_over_mu_rel_vec(index_E_vec,1));
b_W_vec = (1/(24*L_X))*...                                                 % vector containing the b_W-coefficients
    (H3_over_mu_rel_vec+H3_over_mu_rel_vec(index_W_vec,1));



% ----------------------------------------------------------------------- %
% --- FIXED-POINT ITERATION --------------------------------------------- %
% ----------------------------------------------------------------------- %


convergent = 0;                                                            % flag for convergence; will be set to 1 once the solution has converged

for iter = 1:iter_max                                                      % loop for iterative solution of nonlinear BVP
	
	
	
    % ------------------------------------------------------------------- %
    % --- (1;2) ASSEMBLY UNDER CONSIDERATION OF OIL SUPPLY BCs ---------- %
    % ------------------------------------------------------------------- %
	
	
    g_E_vec = g_vec(index_E_vec,1);                                        % vector containing the switch function at neighboring node E for each node 0
    g_W_vec = g_vec(index_W_vec,1);                                        % vector containing the switch function at neighboring node W for each node 0
	
    R_vec = (sgn_u*0.5)*(H_E_vec-H_W_vec) + L_X*dHdT_vec ...               % RHS vector; note that (1-g_pts_vec).*Pi_bar_pts_vec describes the nodal film fractions (minus 1) at the previous time step
        - (L_X/Delta_T)*(H_vec.*(1-g_pts_vec).*Pi_bar_pts_vec);
	
    diag_E2_vec = (b_E_vec+b_W_vec).*g_vec+H_vec.*(1-g_vec)...             % main diagonal of E2 before BCs
        +L_X*((dHdT_vec+H_vec/Delta_T).*(1-g_vec)); 
    Ldiag_E2_vec = ...                                                     % left/lower diagonal of E2 before BCs
        -b_W_vec.*g_W_vec+(-sgn_u-1)*0.5*(H_W_vec.*(1-g_W_vec));  
    Rdiag_E2_vec = ...                                                     % right/upper diagonal of E2 before BCs
        -b_E_vec.*g_E_vec+(sgn_u-1)*0.5*(H_E_vec.*(1-g_E_vec));   
	
    R_vec(end_os_vec+1,1) = ...                                            % consideration of the effect of the oil supply pressure on the nodes next to the oil supply groove via the RHS vector
        R_vec(end_os_vec+1,1) + Pi_os*Ldiag_E2_vec(end_os_vec+1,1);
    R_vec(vertcat(n_x,start_os_vec(2:grooves,1)-1),1) = ...
        R_vec(vertcat(n_x,start_os_vec(2:grooves,1)-1),1) + ...
        Pi_os*Rdiag_E2_vec(vertcat(n_x,start_os_vec(2:grooves,1)-1),1);
	
    if grooves > 1
        Ldiag_E2_vec(start_os_vec(2:grooves,1),1) = 0;                     % replace by zero (elimination due to oil supply BCs)
        Rdiag_E2_vec(start_os_vec(2:grooves,1)-1,1) = 0;                   % replace by zero (elimination due to oil supply BCs)
    end
	
    R_vec_BC = R_vec(nodes_dof_vec,1);                                     % elimination due to oil supply BCs
    diag_E2_vec_BC = diag_E2_vec(nodes_dof_vec,1);                         % elimination due to oil supply BCs
    Ldiag_E2_vec_BC = Ldiag_E2_vec(nodes_dof_vec(1:(dof-1),1)+1,1);        % elimination due to oil supply BCs
    Rdiag_E2_vec_BC = Rdiag_E2_vec(nodes_dof_vec(1:(dof-1),1),1);          % elimination due to oil supply BCs
	
    indices_Ldiag_vec = ...                                                % indices corresponding to the left/lower diagonal of the matrix (using one instead of two indices per matrix entry)
        dof*(index_0_vec(1:(dof-1),1)-1)+index_0_vec(2:dof,1);
	
    E2_mat_BC = diag(diag_E2_vec_BC);                                      % define E2-matrix after enforcement of BCs, first only the main diagonal, ...
    E2_mat_BC(indices_Ldiag_vec) = Ldiag_E2_vec_BC;                        % ... then add left/lower diagonal, ...
    E2_mat_BC(indices_Ldiag_vec+(dof-1)) = Rdiag_E2_vec_BC;                % ... then add right/upper diagonal
	
    if grooves == 0                                                        % for the special case without oil supply BCs, we need to ...
        E2_mat_BC(1,dof) = Ldiag_E2_vec(1,1);                              % ... add the E2-entry that couples the first DOF with the last DOF (these DOFs have been decoupled because there is usually an oil supply groove between them; this decoupling is reversed here)
        E2_mat_BC(dof,1) = Rdiag_E2_vec(dof,1);                            % ... add the E2-entry that couples the last DOF with the first DOF (these DOFs have been decoupled because there is usually an oil supply groove between them; this decoupling is reversed here)
    end
    
    g_vec_BC = g_vec(nodes_dof_vec,1);                                     % switch function vector only for the nodes with DOFs (without oil supply BCs)
	
    diag_E0_vec_BC = (L_X/(12*L_Y^2))*...                                  % vector representing the diagonal of E0 (diagonal matrix) after enforcement of BCs
        (H3_over_mu_rel_vec(nodes_dof_vec,1).*g_vec_BC);
	
	
	
    % ------------------------------------------------------------------- %
    % --- (3) ANALYZE WHICH DOFs ARE LOCATED IN WHICH FLOW REGIME ------- %
    % ------------------------------------------------------------------- %
	
	
    n_p = sum(g_vec_BC);                                                   % how many DOFs are located in the pressure zone
    n_c = dof-n_p;                                                         % how many DOFs are located in the cavitation zone
	
    nodes_p_vec = nonzeros((1:dof)'.*g_vec_BC);                            % all DOFs located in the pressure zone (now assuming that the DOFs are numbered from 1 to dof)
    nodes_c_vec = nonzeros((1:dof)'.*(1-g_vec_BC));                        % all DOFs located in the cavitation zone (now assuming that the DOFs are numbered from 1 to dof)

    if tay == 0 || epsilon > epsilon_max                                   % if the eigenvalue problem is tackled with an eigensolver (as opposed to using Taylor approximations)
        n_mod = n_p;                                                       % we use the full set of modes (i.e., the number of modes n_mod is equal to the number of nodes in the pressure zone) because the effort for evaluating unnecessary modes is negligible compared to the effort of the eigensolver (also note that eig is more efficient than eigs for our matrix size)
    else                                                                   % if the eigenvalue problem is solved by Taylor approximations (as opposed to using an eigensolver)
        n_mod = red;                                                       % a modal reduction will be used (the number of modes n_mod is prescribed by red, with red < n_p) because the algorithm with Taylor approximations, while fast, is sensitive to the number of considered modes (and the database doesn't provide data for more modes than prescribed by red)
    end
    
    
	
    % ------------------------------------------------------------------- %
    % --- (4) TRANSFORMATION OF EIGENVALUE PROBLEM ---------------------- %
    % ------------------------------------------------------------------- %
	
	
    if tay == 0 || epsilon > epsilon_max                                   % if no Taylor approximations are used (i.e., if the eigenvalue problem will be tackled by an eigensolver)
        trafo_vec = diag_E0_vec_BC(nodes_p_vec,1).^(-0.5);                 % diagonal of the diagonal matrix which transforms the eigenvalue problem from generalized to standard
        B_mat = trafo_vec.*E2_mat_BC(nodes_p_vec,nodes_p_vec).*...         % matrix B of standard eigenvalue problem
            trafo_vec';
        B_mat = (B_mat+B_mat')/2;                                          % in case rounding errors render this matrix unsymmetric, make it symmetric
    end

	
	
    % ------------------------------------------------------------------- %
    % --- (5;6;7;8) SOLUTION OF EIGENVALUE PROBLEM ---------------------- %
    % ------------------------------------------------------------------- %
	
	
    if tay == 0 || epsilon > epsilon_max                                   % if no Taylor approximations are used (i.e., if the eigenvalue problem will be tackled by an eigensolver)
        [Phi_p_mat,lambda_square_mat,~] = eig(B_mat);                      % solution of the tridiagonal symmetric eigenvalue problem --> the matrix Phi_p_mat contains the eigenvectors as columns and the diagonal matrix lambda_square_mat contains the squares of the eigenvalues
        lambda_square_vec = diag(lambda_square_mat);                       % vector containing the squares of the eigenvalues
        if grooves == 0 && n_c == 0                                        % if no oil supply BCs are modeled and all nodes are in the pressure zone, then the smallest eigenvalue is numerically zero; a negative value should be avoided, so we ...
            [lambda_square_vec,sort_vec] = ...                             % ... sort the eigenvalues in ascending order, ...
                sort(lambda_square_vec,'ascend');
            lambda_square_vec(1,1) = 0;                                    % ... set the first eigenvalue to exactly zero, and ...
            Phi_p_mat = Phi_p_mat(:,sort_vec);                             % ... sort the eigenvectors for consistency
        end
        % Phi_p_mat = Phi_p_mat.*(diag(Phi_p_mat'*Phi_p_mat).^(-0.5))';      % normalization of the eigenvectors (not necessary because the chosen eigensolver already normalizes them as desired)
        Phi_p_mat = trafo_vec.*Phi_p_mat;                                  % transformation of the eigenvectors
    else                                                                   % otherwise (i.e., if the eigenvalue problem is solved by Taylor approximations)
        Delta_eps = epsilon-epsilon_constr;                                % difference in relative eccentricity between the current shaft position and the one where the Taylor series coefficients were derived
        Phi_p_mat = zeros(n_x,n_mod);                                      % initialization of the matrix that will later contain the eigenvectors
        lambda_square_vec = zeros(n_mod,1);                                % initialization of the vector that will later contain the squared eigenvalues
        for j = 0:n_tay                                                    % loop through all powers involved in the Taylor series
            fac = Delta_eps^j/factorial(j);                                % factor to multiply coefficient by
            Phi_p_mat = Phi_p_mat + ...                                    % add contribution of Taylor series coefficient corresponding to current power to the matrix of eigenvectors
                Vd_mat(:,(1+j*n_mod):((j+1)*n_mod))*fac;
            lambda_square_vec = lambda_square_vec + ...                    % add contribution of Taylor series coefficient corresponding to current power to the vector of eigenvalues
                Ld_vec(1,(1+j*n_mod):((j+1)*n_mod))'*fac;
        end
    	Phi_p_mat(:,1) = 1;                                                % first eigenvector (corresponding to the zero-eigenvalue) does not need to be approximated, as this eigenvector is theoretically constant and equal to [1 1 1 ...]^T (apart from the normalization)
        Phi_p_mat = Phi_p_mat.*(diag(Phi_p_mat'*...                        % normalization
            (diag_E0_vec_BC.*Phi_p_mat)).^(-0.5))';
        lambda_square_vec(1,1) = 0;                                        % set the first eigenvalue to exactly zero
    end
    lambda_vec = sqrt(lambda_square_vec);                                  % vector containing the eigenvalues
    tanh_lambda_vec = tanh(lambda_vec);                                    % vector containing the tanh of every eigenvalue
    
	
	
    % ------------------------------------------------------------------- %
    % --- (9;10;11) LINEAR EQUATION SYSTEMS ----------------------------- %
    % ------------------------------------------------------------------- %
	
	
    if n_c >= 1                                                            % if a cavitation zone exists, follow the general algorithm
        
        Phi_p_inv_mat = (Phi_p_mat').*(diag_E0_vec_BC(nodes_p_vec,1)');    % inverse of the eigenvector matrix
        
        E2_pp_inv_mat_BC = Phi_p_mat*((1./lambda_square_vec).*...          % inverse of the E2-matrix in the pressure zone
            (Phi_p_mat'));
        
        A_mat = ( E2_mat_BC(nodes_c_vec,nodes_p_vec) - ...                 % this matrix A represents the bulky expression that appears on both sides of the equation system for Pi_con_c_vec
            E2_mat_BC(nodes_c_vec,nodes_p_vec)*Phi_p_mat*...
            ((tanh_lambda_vec./lambda_vec).*Phi_p_inv_mat) ) * ...
            E2_pp_inv_mat_BC;
        
        D_mat = A_mat*E2_mat_BC(nodes_p_vec,nodes_c_vec) - ...             % matrix on the left-hand side of the equation system for Pi_con_c_vec
            E2_mat_BC(nodes_c_vec,nodes_c_vec);
        
        rhs_vec = R_vec_BC(nodes_c_vec,1) - A_mat*R_vec_BC(nodes_p_vec,1); % vector on the right-hand side of the equation system for Pi_con_c_vec
        
        Pi_con_c_vec = zeros(n_c,1);                                       % the solution in the cavitation zone will be stored in this vector
        
        if sgn_u == 1
            Pi_con_c_vec(1,1) = rhs_vec(1,1)/D_mat(1,1);                   % computation of the first vector entry of Pi_con_c_vec (the solution of the equation system is trivial because the upwind scheme renders the matrix D bidiagonal)
            for i = 2:n_c
                Pi_con_c_vec(i,1) = (rhs_vec(i,1)-D_mat(i,i-1)*...         % computation of the other vector entries of Pi_con_c_vec (the solution of the equation system is trivial because the upwind scheme renders the matrix D bidiagonal)
                    Pi_con_c_vec(i-1,1))/D_mat(i,i);
            end
        else
            Pi_con_c_vec(n_c,1) = rhs_vec(n_c,1)/D_mat(n_c,n_c);           % computation of the last vector entry of Pi_con_c_vec (the solution of the equation system is trivial because the upwind scheme renders the matrix D bidiagonal)
            for i = (n_c-1):(-1):1
                Pi_con_c_vec(i,1) = (rhs_vec(i,1)-D_mat(i,i+1)*...         % computation of the other vector entries of Pi_con_c_vec (the solution of the equation system is trivial because the upwind scheme renders the matrix D bidiagonal)
                    Pi_con_c_vec(i+1,1))/D_mat(i,i);
            end
        end
        
        Pi_con_p_vec = - ( E2_pp_inv_mat_BC * ( ...                        % compute axially constant component of the solution in the pressure zone
            E2_mat_BC(nodes_p_vec,nodes_c_vec)*Pi_con_c_vec + ...
            R_vec_BC(nodes_p_vec,1) ) );
        
        C_vec = - ( Phi_p_inv_mat*Pi_con_p_vec );                          % compute integration constants
	    
    else                                                                   % else (i.e., if all nodes are in the pressure zone), use a simpler algorithm
        
        if grooves == 0                                                    % if no oil supply groove is present, then ...
            C_vec = vertcat(0,1./lambda_square_vec(2:n_mod,1)).*...        % ... compute integration constants, excluding the singular mode
                ((Phi_p_mat')*R_vec_BC);
        else                                                               % otherwise, ...
            C_vec = (1./lambda_square_vec).*((Phi_p_mat')*R_vec_BC);       % ... compute integration constants (no singular mode exists)
        end
        
        Pi_con_p_vec = - ( Phi_p_mat*C_vec );                              % compute axially constant component of the solution in the pressure zone
        Pi_con_c_vec = zeros(0,1);                                         % the solution in the cavitation zone does not exist
        
    end
    
	
    
    % ------------------------------------------------------------------- %
    % --- (12) EXPRESS RESULTS GLOBALLY (BOTH ZONES + OIL SUPPLY) ------- %
    % ------------------------------------------------------------------- %
	
	
    Pi_con_vec = Pi_os*ones(n_x,1);                                        % define one global vector for the axially constant component of the solution (including cavitation zone, pressure zone, and oil supply); first, fill all entries with the oil supply pressure, then ...
    Pi_con_vec(nodes_dof_vec(nodes_p_vec,1),1) = Pi_con_p_vec;             % ... overwrite the oil supply pressure by the computed solution at the DOFs in the cavitation zone and ...
    Pi_con_vec(nodes_dof_vec(nodes_c_vec,1),1) = Pi_con_c_vec;             % ... overwrite the oil supply pressure by the computed solution at the DOFs in the pressure zone
	
    Phi_mat = zeros(n_x,n_mod);
    Phi_mat(nodes_dof_vec(nodes_p_vec,1),:) = Phi_p_mat;                   % matrix containing the eigenvectors in the pressure zone and zeros everywhere else
	
	
	
    % ------------------------------------------------------------------- %
    % --- (13;14;15) AXIAL AVERAGE OF SOLUTION AND SWITCH FUNCTION ------ %
    % ------------------------------------------------------------------- %
    
	
    Pi_bar_vec = Pi_con_vec;                                               % vector for storing the axial average of the pressure-like function Pi at all nodes
    if n_p >= 1
        if grooves == 0 && n_c == 0                                        % if no oil supply BCs are modeled and all nodes are in the pressure zone, ...
            Pi_bar_vec = Pi_bar_vec + Phi_mat*(vertcat(1,...               % ... perform integration while accounting for singular mode
                tanh_lambda_vec(2:n_mod,1)./...
                lambda_vec(2:n_mod,1)).*C_vec);
        else                                                               % otherwise, ...
            Pi_bar_vec = Pi_bar_vec + ...                                  % ... perform integration (no singular mode exists)
                Phi_mat*(tanh_lambda_vec./lambda_vec.*C_vec);
        end
    end
    
    g_vec_old_old = g_vec_old;                                             % save old switch function as old old switch function
    g_vec_old = g_vec;                                                     % save current switch function as old switch function
    g_vec = sign(sign(Pi_bar_vec)+1);                                      % compute new switch function and save as current switch function
	
    if max(abs(g_vec-g_vec_old)) == 0                                      % if the solution is convergent, then ...
        convergent = 1;                                                    % ... flag it as convergent and ...
        break                                                              % ... end the iterative scheme
    end
	
    if iter >= 2 && max(abs(g_vec-g_vec_old_old)) == 0                     % if the switch function oscillates between two configurations (actually, I think this doesn't happen with the current SBFEM model), then ...
        convergent = 1;                                                    % ... consider the solution convergent anyway and ...
        break                                                              % ... end the iterative scheme
    end
    
    if guembel == 1                                                        % if Guembel is used instead of Elrod ...
        break                                                              % ... quit iterative scheme after first iteration
    end
	
end



% ----------------------------------------------------------------------- %
% --- (16;17) FILM FRACTIONS, PRESSURES, AXIAL GRADIENTS ---------------- %
% ----------------------------------------------------------------------- %


theta_vec = (1-g_vec).*Pi_bar_vec+1;                                       % vector of nodal film fractions
p_bar_vec = p_ref*(Pi_bar_vec.*g_vec);                                     % vector of nodal physical (as opposed to dimensionless) pressures, averaged in the axial direction
dpdy_bb_vec = (Phi_mat*(tanh_lambda_vec.*lambda_vec.*C_vec))*...           % vector of axial pressure gradients at the bearing boundary at the nodes
    (2*p_ref/l_b);

if guembel == 1                                                            % if Guembel is used instead of Elrod, ...
    theta_vec = ones(n_x,1);                                               % ... set film fractions to 1 and ...
    dpdy_bb_vec = dpdy_bb_vec.*g_vec;                                      % ... apply Guembel condition to axial pressure gradients at bearing boundary (zero in the cavitation zone)
end



% ----------------------------------------------------------------------- %
% --- (18) FORCES, FRICTION MOMENT, OIL VOLUME, VOLUME FLOW ------------- %
% ----------------------------------------------------------------------- %


l_x = L_X*d_b/2;                                                           % physical (as opposed to nondimensionalized) circumferential sector side length

F_1 = (cos(X_vec')*p_bar_vec)*l_b*l_x;                                     % component of hydrodynamic force
F_2 = (sin(X_vec')*p_bar_vec)*l_b*l_x;                                     % component of hydrodynamic force

if tay == 1                                                                % if the eigenvalue problem was solved by Taylor approximations, the attitude angle was manipulated (set to zero); this is now compensated via transforming the forces, see below
    F_h = F_1*cos(X_att_trafo) - F_2*sin(X_att_trafo);                     % transform first force component (F_h is only used for temporary storage)
    F_v = F_1*sin(X_att_trafo) + F_2*cos(X_att_trafo);                     % transform second force component (F_v is only used for temporary storage)
    F_1 = F_h;                                                             % save transformed force component (this force component now corresponds to the correct attitude angle)
    F_2 = F_v;                                                             % save transformed force component (this force component now corresponds to the correct attitude angle)
end

F_h = F_1*cos(angle_shell) - F_2*sin(angle_shell);                         % transformation into inertial system --> horizontal hydrodynamic force acting on the shell [N]
F_v = F_1*sin(angle_shell) + F_2*cos(angle_shell);                         % transformation into inertial system --> vertical hydrodynamic force acting on the shell [N]

dpdx_vec = (vertcat(p_bar_vec(2:n_x,1),p_bar_vec(1,1))-...                 % vector of axially averaged circumferential pressure gradients at the nodes
    vertcat(p_bar_vec(n_x,1),p_bar_vec(1:(n_x-1),1)))/(2*l_x);

M_fr = (d_b/2*l_b*l_x) * (-0.5*c*(H_vec'*dpdx_vec)+...                     % oil friction moment acting on the shell
    (u/c)*(theta_vec'*(mu_vec./H_vec)));

V_oil = (l_b*l_x*c)*(theta_vec'*H_vec);                                    % oil volume in the bearing

V_dot_bb = 0;
if n_p >= 1
    V_dot_bb = V_dot_bb + (l_x/6*c^3)*(H3_vec'*(dpdy_bb_vec./mu_vec));     % oil volume flow through the bearing boundaries
end



% ----------------------------------------------------------------------- %
% --- EVALUATE SOLUTION ON TWO-DIMENSIONAL GRID, IF DESIRED ------------- %
% ----------------------------------------------------------------------- %


if tay == 1                                                                % if the eigenvalue problem was solved by Taylor approximations, the attitude angle was manipulated (set to zero); the lines below contain some preperations for interpolating the solution in a way that compensates this manipulation
    di = mod((X_att_trafo/(2*pi)),1)*n_x;                                  % number of elements between node 1 and the circumferential position defined by the attitude angle
    fl = di-floor(di);                                                     % weight for left node in interpolation
    fr = 1-fl;                                                             % weight for right node in interpolation
    di = floor(di);                                                        % number of elements between node 1 and the circumferential position defined by the attitude angle, rounded down to an integer
end

if n_y > 1
    Pi_mat = zeros(n_x,n_y);                                               % this matrix will contain the two-dimensional distribution of the computed pressure-like function Pi
    j = round(n_y/2+0.75);                                                 % axial starting index for evaluation of the solution (as one half of the bearing is initially ignored due to the symmetry)
    l = n_y-j+1;                                                           % number of data points in the axial direction per half of the bearing
    xi_vec = (linspace(j,n_y,l)-1)*(2/(n_y-1))-1;                          % vector containing the axial positions where the solution will be evaluated (dimensionless axial coordinate xi)
    Pi_mat(:,j:n_y) = repmat(Pi_con_vec,[1,l]);                            % consider axially constant component of the solution
    if n_p >= 1
        Pi_mat(:,j:n_y) = Pi_mat(:,j:n_y) + ...
            (Phi_mat.*(C_vec./cosh(lambda_vec))')*cosh(lambda_vec*xi_vec); % consider axially variable component of the solution
    end
    if tay == 0                                                            % if the eigenvalue problem was solved by an eigensolver (i.e. not by a Taylor approximation), oil supply BCs might be present, which may have caused the algorithm to shift the node numbers in order to ensure a tridiagonal eigenvalue problem; this shift is now reversed in Pi_mat
        Pi_mat(shift_vec,j:n_y) = Pi_mat(:,j:n_y);                         % shift node numbers back to normal
    else                                                                   % if the eigenvalue problem was solved by a Taylor approximation (i.e. not by an eigensolver), the attitude angle was manipulated (set to zero); compensating this manipulation requires, amongst other steps, the interpolation performed below
        Pi_mat(:,j:n_y) = vertcat( ...                                     % interpolate pressure-like function to compensate manipulation of attitude angle
            fr*Pi_mat([(n_x-di+1):n_x,1],j:n_y) + ...
            fl*Pi_mat((n_x-di):n_x,j:n_y), ...
            fr*Pi_mat(2:(n_x-di),j:n_y) + ...
            fl*Pi_mat(1:(n_x-di-1),j:n_y) );
    end
    Pi_mat(:,1:l) = fliplr(Pi_mat(:,j:n_y));                               % copy solution to the other half of the bearing
else
    Pi_mat = [];                                                           % if no evaluation of the solution on a twodimensional grid is desired, Pi_mat is defined as an empty variable
end



% ----------------------------------------------------------------------- %
% --- SAVE SOME RESULTS AND SHIFT NODE NUMBERS BACK TO NORMAL ----------- %
% ----------------------------------------------------------------------- %


if tay == 0                                                                % if the eigenvalue problem was solved by an eigensolver (i.e. not by a Taylor approximation), oil supply BCs might be present, which may have caused the algorithm to shift the node numbers in order to ensure a tridiagonal eigenvalue problem; this shift is now reversed in some vectors
    
    pts_vec(shift_vec,1) = Pi_bar_vec;                                     % save axially averaged pressure-like function while shifting the node numbers back to normal (this pressure-like function will be considered at the next time step as previous solution)
    g_vec(shift_vec,1) = g_vec;                                            % vector of switch functions: shift node numbers back to normal
    theta_vec(shift_vec,1) = theta_vec;                                    % vector of film fractions: shift node numbers back to normal
    
else                                                                       % if the eigenvalue problem was solved by a Taylor approximation (i.e. not by an eigensolver), the attitude angle was manipulated (set to zero); compensating this manipulation requires, amongst other steps, the interpolations performed below
    
    pts_vec(1:n_x,1) = vertcat( ...                                        % interpolate pressure-like function (averaged in the axial direction) to compensate manipulation of attitude angle
        fr*Pi_bar_vec([(n_x-di+1):n_x,1],1) + ...
        fl*Pi_bar_vec((n_x-di):n_x,1), ...
        fr*Pi_bar_vec(2:(n_x-di),1) + ...
        fl*Pi_bar_vec(1:(n_x-di-1),1) );
    g_vec = sign(sign(pts_vec(1:n_x,1))+1);                                % compute switch function again for interpolated solution
    theta_vec = (1-g_vec).*pts_vec(1:n_x,1)+1;                             % compute film fraction again for interpolated solution
    
end

pts_vec(n_x+1,1) = t;                                                      % save current time (will be used at the next time step for formulating the rate of change of the film fraction)



end
