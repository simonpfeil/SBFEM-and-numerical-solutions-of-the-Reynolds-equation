function [F_h,F_v,M_h,M_v,M_fr,V_oil,V_dot_bb,Pi_mat,pts_vec,p_ref,...
    convergent,iter] = fvm_elrod(d_b,l_b,c,grooves,X_os,L_X_os,l_y_os,...
    p_os,ac_vec,t,angle_shell,omega_shell,dis_h_shell,dis_v_shell,...
    vel_h_shell,vel_v_shell,omega_shaft,dis_h_shaft,dis_v_shaft,...
    vel_h_shaft,vel_v_shaft,tilt_h_shell,tilt_v_shell,tilt_dot_h_shell,...
    tilt_dot_v_shell,tilt_h_shaft,tilt_v_shaft,tilt_dot_h_shaft,...
    tilt_dot_v_shaft,n_x,iter_max,mu_vec,n_y,pts_vec,quasistatic,...
    symBC,guembel)

% PROGRAM FOR NUMERICAL SOLUTION OF THE REYNOLDS EQUATION WITH CAVITATION

% Author: Simon Pfeil
% Affiliation: OvGU Magdeburg, Germany (until 30.06.2025)
% Date: 30.06.2025

% --- comments ------------------------------------------------------------
%
% - This algorithm solves the Reynolds equation for hydroynamic journal 
%   bearings.
%
% - The discretized form of the Reynolds equation used in this computation
%   can be derived by the finite element method (FEM), the finite volume 
%   method (FVM), or the finite difference method (FDM). To keep the 
%   comments and explanations brief, we will only use the FVM terminology.
%
% - A transient, nonlinear cavitation model is used (optionally, the simple
%   Guembel approach can be used instead).
%
% - This program is designed to be incorporated into time integration 
%   schemes for rotordynamic (or multibody) simulations; however, a quasi-
%   static solution is also possible.
%
% - The input and output variables are explained below.
%
% - Some input variables (ac_vec and mu_vec) can only be defined with
%   knowledge of the node numbering scheme (unless a uniform value in both 
%   directions is desired). The node numbering scheme is illustrated below
%   the discussion of input and output variables.
%
% - See the attached PDF for clarification of the coordinate systems. 
%   x (or the nondimensionalized X) and y (or the nondimensionalized Y) are 
%   the coordinates of the lubrication gap and are used by the Reynolds 
%   equation. This coordinate system, as well as the computational grid, 
%   are fixed in the reference frame of the shell. The transformations 
%   between the reference frame of the shell and the inertial system are 
%   performed internally by this program. The kinematic variables of the 
%   shell and the shaft, when handed to this program, should be formulated 
%   from the perspective of the inertial system. The hydrodynamic forces 
%   computed by this program are automatically transformed into the 
%   inertial system before output. In contrast, those input and output 
%   variables that describe one- or two-dimensional fields (e.g., viscosity 
%   distribution, additional contour, and pressure-like function) never use 
%   the inertial system, as they are expressed by arrays whoose entries 
%   represent the values at the nodes (which are always fixed at the 
%   shell). The circumferential position of the oil supply groove X_os is 
%   also always expressed in the reference frame of the shell; otherwise, 
%   a rotating shell would require us to keep updating this input variable.
%
% - More thorough discussions of the computational method, of the 
%   cavitation model, of the boundary conditions (BCs), and of further 
%   assumptions will be given in the thesis "Simulating Hydrodynamic 
%   Bearings with the Scaled Boundary Finite Element Method" by Simon 
%   Pfeil, but this thesis isn't publically available yet

% --- input variables -----------------------------------------------------
%
% - d_b = bearing diameter [m]
%
% - l_b = bearing length [m]
%
% - c = radial clearance (shell radius minus shaft radius) [m]
%
% - grooves = number of oil supply grooves (at least 1, due to the 
%   consideration of cavitation), which are assumed to be distributed 
%   equidistantly over the circumference [-]
%
% - X_os = angular circumferential position of the center of one of the
%   oil supply grooves (must not be smaller than 0 or greater than 2*pi,
%   see attached PDF); since the coordinate X is fixed at the shell, a 
%   rotation of the shell does not change the value of X_os [rad]
%
% - L_X_os = angular circumferential side length of the oil supply groove
%   (see attached PDF file) [rad]
%
% - l_y_os = axial side length of the oil supply groove [m]
%
% - p_os (if positive or zero) = pressure prescribed in the oil supply 
%   groove (zero corresponds to atmospheric pressure) [Pa]
%
% - p_os (if negative) = film fraction prescribed in the oil supply groove
%   minus 1 (for example, a film fraction of 0.7 is prescribed by setting 
%   p_os = 0.7-1 = -0.3) [-]
%
% - ac_vec = column vector of length n_x*n_y containing an additional 
%   contour in case the shell is not cylindrical; the entries of this 
%   vector are simply added to the nodal gap widths (consider the given
%   node numbering scheme); since the nodes are fixed at the shell, a
%   rotation of the shell does not require this additional contour to be 
%   updated; note that an additional contour of the shaft is not allowed 
%   in the current version of this program [m]
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
% - tilt_h_shell, tilt_h_shaft = tilting angles of the shell and the
%   shaft, respectively, around a horizontal axis (see attached PDF) [rad]
%
% - tilt_v_shell, tilt_v_shaft = tilting angles of the shell and the
%   shaft, respectively, around a vertical axis (see attached PDF) [rad]
%
% - tilt_dot_h_shell, tilt_dot_h_shaft = rates of change of the tilting 
%   angles of the shell and the shaft, respectively, around a horizontal
%   axis (see attached PDF) [rad/s]
%
% - tilt_dot_v_shell, tilt_dot_v_shaft = rates of change of the tilting 
%   angles of the shell and the shaft, respectively, around a vertical 
%   axis (see attached PDF) [rad/s]
%
% - n_x = number of nodes in the circumferential direction; the smallest
%   circumferential node number 1 is located at the angular position X = 0 
%   and the largest one n_x is located at X = 2*pi*(1-1/n_x); the node at 
%   X = 2*pi is not included because the periodic BC merges this node with 
%   the one at X = 0 (see attached PDF file and illustration of the node
%   numbering scheme below) [-]
%
% - iter_max = maximum allowed number of iterations for the cavitation 
%   model (usually, less than 10 iterations are enough, but in some cases,
%   more may be necessary; I usually choose iter_max = n_x) [-]
%
% - mu_vec = column vector of length n_x*n_y containing the dynamic oil
%   viscosities at the nodes (consider the given node numbering scheme); 
%   since the nodes are fixed at the possibly rotating shell, the 
%   viscosity distribution is always expressed in the shell's reference 
%   frame [Pa*s]
%
% - n_y = axial number of nodes (see illustration of the node numbering 
%   scheme below); the smallest and largest axial node numbers, 1 and n_y, 
%   respectively, are located at the two bearing boundaries; ideally, n_y 
%   should be defined as a function of n_x and of the bearing dimensions in 
%   such a way that the control volume side lengths are almost equal in 
%   both directions [-]
%
% - pts_vec = column vector of length n_x*n_y+1 where the first n_x*n_y 
%   entries are the nodal pressure-like functions of the previous time step 
%   and the entry n_x*n_y+1 states the time of the previous time step; this 
%   variable is used for transferring data between time steps, as required 
%   by the transient cavitation model; note that pts_vec is also an output
%   variable; simply take the pts_vec that was generated as output at the 
%   previous valid time step and provide it as input variable at the 
%   current time step [-]/[s]
%
% - quasistatic = a flag that specifies whether a quasi-static solution is
%   desired (1 = yes, 0 = no), i.e., whether time derivatives terms in the 
%   Reynolds equation should be ignored [-]
%
% - symBC = a flag which states whether or not a symmetric boundary
%   condition at the axial center of the bearing y = 0 should be used: 
%   0 = no, 1 = yes; setting symBC = 1 reduces the computational cost but 
%   also renders the model incapable of considering shaft tilting [-]
%
% - guembel = determines whether Guembel conditions should be assumed,
%   bypassing the Elrod cavitation model; guembel = 0 --> use Elrod, 
%   guembel = 1 --> use Guembel [-]

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
% - M_h = hydrodynamic moment acting on the shell around a horizontal 
%   axis (see attached PDF); the same moment acts on the shaft in the 
%   opposite direction [Nm]
%
% - M_v = hydrodynamic moment acting on the shell around a vertical
%   axis (see attached PDF); the same moment acts on the shaft in the 
%   opposite direction [Nm]
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
% - Pi_mat = matrix of size n_x*n_y stating the computed pressure-like 
%   function Pi at all nodes; from this, the pressure distribution can be 
%   derived (p = p_ref*Pi*g), where g is the switch function (see above) 
%   and p_ref is the reference pressure which was used for the 
%   nondimensionalization (see below); since the nodes are fixed at the 
%   possibly rotating shell, the two-dimensional field represented by 
%   Pi_mat is always expressed in the shell's reference frame [-]
%
% - pts_vec = column vector of length n_x*n_y+1 where the first n_x*n_y 
%   entries contain the newly computed nodal pressure-like functions Pi 
%   and the entry n_x*n_y+1 contains the current time t [-]/[s]
%
% - p_ref = reference pressure due to nondimensionalization [Pa]
%
% - convergent = a flag that indicates whether or not the solution has 
%   converged (1 = yes, 0 = no) [-]
%
% - iter = number of iterations [-]

% --- node numbering scheme -----------------------------------------------
%
%
%   axial direction (y or Y) with n_y nodes
%
%   ^
%   |
%   |       ...                                 ...      n_y*n_x
%   |
%   |       ...
%           2*n_x+1  2*n_x+2  ...
%           n_x+1    n_x+2    n_x+3    n_x+4    ...      2*n_x
%           1        2        3        4        ...      n_x
%
%           -------->  circumferential direction (x or X) with n_x nodes
%
%
% -------------------------------------------------------------------------



% ----------------------------------------------------------------------- %
% --- EXPRESS INPUTS W.R.T. COORDINATE SYSTEM FOR REYNOLDS {x,y} -------- %
% ----------------------------------------------------------------------- %


dis_h = dis_h_shaft-dis_h_shell;                                           % horizontal displacement of the shaft relative to the shell, still in the reference frame of the inertial system
dis_v = dis_v_shaft-dis_v_shell;                                           % vertical displacement of the shaft relative to the shell, still in the reference frame of the inertial system
q = max(sqrt(dis_h^2+dis_v^2),eps);                                        % absolute eccentricity
X_att = atan2(dis_v,dis_h);                                                % attitude angle in the reference frame of the inertial system

tilt_h = tilt_h_shaft-tilt_h_shell;                                        % tilting angle of the shaft relative to the shell around the horizontal axis, still in the reference frame of the inertial system
tilt_v = tilt_v_shaft-tilt_v_shell;                                        % tilting angle of the shaft relative to the shell around the vertical axis, still in the reference frame of the inertial system
tilt = max(sqrt(tilt_h^2+tilt_v^2),eps);                                   % absolute tilting angle
X_tilt = atan2(tilt_v,tilt_h);                                             % angle describing the resulting tilting axis in the reference frame of the inertial system

X_att = X_att-angle_shell;                                                 % attitude angle in the shell-fixed reference frame in which the Reynolds equation is solved
X_tilt = X_tilt-angle_shell;                                               % angle describing the tilting axis in the shell-fixed reference frame in which the Reynolds equation is solved

if quasistatic == 0
    
    omega = omega_shaft-omega_shell;                                       % rotational velocity of the shaft relative to the shell
    
    vel_h = vel_h_shaft-vel_h_shell;                                       % horizontal velocity of the shaft relative to the shell, still in the reference frame of the inertial system
    vel_v = vel_v_shaft-vel_v_shell;                                       % vertical velocity of the shaft relative to the shell, still in the reference frame of the inertial system
    q_dot = (vel_v*dis_v+vel_h*dis_h)/q;                                   % rate of change of absolute eccentricity
    X_att_dot = (vel_v*dis_h-vel_h*dis_v)/q^2;                             % rate of change of attitude angle in the reference frame of the inertial system
    
    tilt_dot_h = tilt_dot_h_shaft-tilt_dot_h_shell;                        % rate of change of the tilting angle of the shaft relative to the shell around the horizontal axis, still in the reference frame of the inertial system
    tilt_dot_v = tilt_dot_v_shaft-tilt_dot_v_shell;                        % rate of change of the tilting angle of the shaft relative to the shell around the vertical axis, still in the reference frame of the inertial system
    tilt_dot = (tilt_dot_v*tilt_v+tilt_dot_h*tilt_h)/tilt;                 % rate of change of the absolute tilting angle
    X_tilt_dot = (tilt_dot_v*tilt_h-tilt_dot_h*tilt_v)/tilt^2;             % rate of change of the angle describing the resulting tilting axis in the reference frame of the inertial system
    
    X_att_dot = X_att_dot-omega_shell;                                     % rate of change of attitude angle in the shell-fixed reference frame in which the Reynolds equation is solved
    X_tilt_dot = X_tilt_dot-omega_shell;                                   % rate of change of the angle describing the tilting axis in the shell-fixed reference frame in which the Reynolds equation is solved
    
else
    
    omega = omega_shaft;                                                   % quasistatic case: rotational velocity of shell is assumed to be zero
    q_dot = 0;                                                             % quasistatic case: rate of change of eccentricity is assumed to be zero
    X_att_dot = 0;                                                         % quasistatic case: rate of change of attitude angle is assumed to be zero
    tilt_dot = 0;                                                          % quasistatic case: rate of change of tilting angle is assumed to be zero
    X_tilt_dot = 0;                                                        % quasistatic case: rate of change of the angle describing the tilting axis is assumed to be zero
    
end



% ----------------------------------------------------------------------- %
% --- NONDIMENSIONALIZATION --------------------------------------------- %
% ----------------------------------------------------------------------- %


u = omega*(d_b/2);                                                         % circumferential surface velocity [m/s]
r_b = d_b/2;                                                               % bearing radius [m]
sgn_u = sign(u);                                                           % sign of circumferential surface velocity [-]
mu_ref = mean(mu_vec);                                                     % reference viscosity [Pa*s]
p_ref = abs(u)*mu_ref*r_b/(2*c^2);                                         % reference pressure for nondimensionalization [Pa]
epsilon = q/c;                                                             % relative eccentricity [-]
epsilon_dot = q_dot/c;                                                     % rate of change of the relative eccentricity [1/s]



% ----------------------------------------------------------------------- %
% --- ANALYZE BOUNDARY VALUE PRESCRIBED IN OIL SUPPLY GROOVE ------------ %
% ----------------------------------------------------------------------- %


if p_os >= 0               												   % "IF a non-negative supply pressure is prescribed, then ..."
	Pi_os = p_os/p_ref;													   % nondimensionalize the supply pressure and save as pressure-like function
	g_os = 1;															   % set corresponding switch function to 1
else																	   % otherwise, the prescribed value p_os will be interpreted as a film fraction with an offset of -1, i.e., as a negative pressure-like function
    Pi_os = p_os;														   % save prescribed value as pressure-like function
    g_os = 0;															   % set corresponding switch function to 0
end

if guembel == 1 && p_os < 0                                                % if Guembel is used instead of Elrod, we cannot prescribe a film fraction as BC (which would be done via a negative p_os), so ...
    Pi_os = 0;                                                             % ... make sure that the oil supply pressure is at least zero
end



% ----------------------------------------------------------------------- %
% --- ANALIZE GRID ------------------------------------------------------ %
% ----------------------------------------------------------------------- %


L_X = 2*pi/n_x;                                                            % angular circumferential side length of control volume [rad]
L_Y = (l_b/(n_y-1))/r_b;                                                   % nondimensionalized axial side length of control volume [-]

if symBC == 0                                                              % if no symmetric BC is used
    dof_y = n_y-2;                                                         % number of DOFs in the axial direction (nummber of nodes without counting the nodes at the bearing boundaries)
    j_start = 2;                                                           % axial node number corresponding to the first DOF in the axial direction
else                                                                       % if a symmetric BC is used
    if mod(n_y,2) == 0                                                     % if the node number is even (no node is located at the axial center of the bearing Y=0)
        dof_y = n_y/2-1;                                                   % number of DOFs in the axial direction (one half of the bearing, not counting the node at the bearing boundary)
    else                                                                   % if the node number is uneven (a node is located at the axial center of the bearing Y=0)
        dof_y = (n_y-1)/2;                                                 % number of DOFs in the axial direction (one half of the bearing, not counting the node at the bearing boundary)
    end
    j_start = n_y-dof_y;                                                   % axial node number corresponding to the first DOF in the axial direction
end

n = n_y*n_x;                                                               % total number of nodes (whole bearing, nodes with Dirichlet BCs are included)
dof = dof_y*n_x;                                                           % total number of DOFs, including the nodes with oil supply BCs (because the penalty method is used) but not the nodes at the bearing boundaries

nn0_vec = linspace((j_start-1)*n_x+1,(n_y-1)*n_x,dof_y*n_x)';              % vector containing the numbers of the nodes which are DOFs
nnN_vec = nn0_vec+n_x;                                                     % vector containing, for every DOF, the number of the neighboring node N in the positive Y-direction
nnE_vec = nn0_vec+1;                                                       % vector containing, for every DOF, the number of the neighboring node E in the positive X-direction
nnS_vec = nn0_vec-n_x;                                                     % vector containing, for every DOF, the number of the neighboring node S in the negative Y-direction
nnW_vec = nn0_vec-1;                                                       % vector containing, for every DOF, the number of the neighboring node W in the negative X-direction

if symBC == 1                                                              % if a symmetric BC is used
    if mod(n_y,2) == 0                                                     % if an even axial number of nodes has been prescribed for the entire bearing length (meaning that no nodes exist at the symmetric boundary)
        nnS_vec(1:n_x,1) = nn0_vec(1:n_x,1);                               % symmetric BC
    else                                                                   % if an uneven axial number of nodes has been prescribed for the entire bearing length (meaning that nodes exist at the symmetric boundary)
        nnS_vec(1:n_x,1) = nnN_vec(1:n_x,1);                               % symmetric BC
    end
end

nnE_vec((1:dof_y)*n_x,1) = nnE_vec((1:dof_y)*n_x,1) - n_x;                 % periodic BC
nnW_vec((0:(dof_y-1))*n_x+1,1) = nnW_vec((0:(dof_y-1))*n_x+1,1) + n_x;     % periodic BC

if symBC == 1                                                              % if a symmetric BC is used
    rows_vec = vertcat(nn0_vec,nn0_vec(1:(dof-n_x),1),nn0_vec,...          % this vector will facilitate the matrix assembly by providing the number of the matrix row for every matrix entry; note that the DOF N is removed (replaced by P=0, i.e., atmospheric pressure) when the discrete Reynolds equation is evaluated at a node next to the bearing boundary
        nn0_vec,nn0_vec);
    cols_vec = vertcat(nn0_vec,nnN_vec(1:(dof-n_x),1),nnS_vec,...          % this vector will facilitate the matrix assembly by providing the number of the matrix column for every matrix entry; note that the DOF N is removed (replaced by P=0, i.e., atmospheric pressure) when the discrete Reynolds equation is evaluated at a node next to the bearing boundary
        nnE_vec,nnW_vec);
    if mod(n_y,2) == 1                                                     % if the axial number of nodes is uneven (meaning that nodes exist at the center of the bearing Y=0)
        halfcv_vec = horzcat(1:n_x,(dof+1):(dof+n_x),(2*dof-n_x+1):...     % some entries of the stiffness matrix will be divided by 2 due to half-sized control volumes (these occur at the center of the bearing Y=0 if a symmetric BC is used and the axial number of nodes is uneven); the vector halfcv_vec will be used to address these entries
            (2*dof),(3*dof-n_x+1):(3*dof),(4*dof-n_x+1):(4*dof));
    end
else                                                                       % if no symmetric BC is used
    rows_vec = vertcat(nn0_vec,nn0_vec(1:(dof-n_x),1),...                  % this vector will facilitate the matrix assembly by providing the number of the matrix row for every matrix entry; note that the DOF N or S, respectively, is removed (replaced by P=0, i.e., atmospheric pressure) when the discrete Reynolds equation is evaluated at a node next to the bearing boundary
        nn0_vec((n_x+1):dof,1),nn0_vec,nn0_vec);
    cols_vec = vertcat(nn0_vec,nnN_vec(1:(dof-n_x),1),...                  % this vector will facilitate the matrix assembly by providing the number of the matrix column for every matrix entry; note that the DOF N or S, respectively, is removed (replaced by P=0, i.e., atmospheric pressure) when the discrete Reynolds equation is evaluated at a node next to the bearing boundary
        nnS_vec((n_x+1):dof,1),nnE_vec,nnW_vec);
end

rows_vec = rows_vec - (j_start-1)*n_x;                                     % shift from node numbers to DOF numbers
cols_vec = cols_vec - (j_start-1)*n_x;                                     % shift from node numbers to DOF numbers

X_vec = repmat(linspace(0,2*pi*(1-1/n_x),n_x)',[n_y,1]);                   % X-positions of all nodes
Y_vec = kron(linspace(-1,1,n_y)'*l_b/d_b,ones(n_x,1));                     % Y-positions of all nodes



% ----------------------------------------------------------------------- %
% --- ANALYZE OIL SUPPLY BCs -------------------------------------------- %
% ----------------------------------------------------------------------- %


X_os_vec = linspace(0,grooves-1,grooves)'*(2*pi/grooves) + X_os;           % vector containing the angular circumferential positions of the centers of all oil supply grooves
start_os_vec = round((X_os_vec-L_X_os/2)/L_X+1);                    	   % vector containing the circumferential node/DOF numbers where the oil supply grooves begin (if they are outside the range 1...nx, this will be corrected later)
end_os_vec = round((X_os_vec+L_X_os/2)/L_X+1);                      	   % vector containing the circumferential node/DOF numbers where the oil supply grooves end (if they are outside the range 1...nx, this will be corrected later)

n_os_x = sum(end_os_vec-start_os_vec+1);                                   % circumferential number of nodes/DOFs where oil supply BCs will be defined
os_x_vec = zeros(n_os_x,1);                                                % this vector serves for storing the circumferential node/DOF numbers where oil supply BCs will be defined
j = 0;                                                                     % j serves to track the number of already considered nodes
for i = 1:grooves                                                          % loop through all grooves
    n_i = end_os_vec(i,1)-start_os_vec(i,1)+1;                             % circumferential number of nodes in current groove
    os_x_vec((j+1):(j+n_i),1) = start_os_vec(i,1):end_os_vec(i,1);         % save node/DOF numbers
    j = j+n_i;                                                             % update number of considered nodes
end
for i = 1:n_os_x                                                           % loop through all circumferential node/DOF numbers where oil supply BCs will be defined
    os_x_vec(i,1) = mod(os_x_vec(i,1)-1,n_x)+1;                            % correct the circumferential node/DOF number if it's currently outside the range 1...nx
end

if l_y_os == 0                                                             % if the axial side length of the groove is exactly equal to zero
    l_y_os = 1e-10;                                                        % set the axial side length of the groove to a value slightly larger than zero (exactly zero causes difficulties if the axial number of nodes is even)
end

start_os = round( ((l_b-l_y_os)/2)/l_b*(n_y-1)+1 );                        % axial node number where the groove begins

if start_os < 2                                                            % if the groove begins at the bearing boundary (start_os = 1) or, illegaly, extends further than the bearing (start_os < 1)
    longgroove = 1;                                                        % set longgroove = 1 as a flag which will later allow to recover this information
    start_os = 2;                                                          % shift the axial starting index of the groove one node away from the bearing boundary (because the nodes located exactly at the bearing boundary have already been excluded from the DOFs)
else
    longgroove = 0;
end

end_os = n_y-start_os+1;                                                   % axial node number where the groove ends

if symBC == 1                                                              % if a symmetric BC is used
    os_y_vec = (1:(end_os-j_start+1))';                                    % list of axial DOF numbers where oil supply BCs will be defined (not equal to the corresponding node numbers, due to the axial shift according to j_start)
else                                                                       % if no symmetric BC is used
    os_y_vec = (start_os:end_os)'-j_start+1;                               % list of axial DOF numbers where oil supply BCs will be defined (not equal to the corresponding node numbers, due to the axial shift according to j_start)
end

n_os_y = length(os_y_vec);                                                 % axial number of DOFs where oil supply BCs will be defined

os_vec = repmat(os_x_vec,[n_os_y,1])+...                                   % DOF numbers where oil supply BCs will be defined (now using a global numbering scheme that considers both directions)
    kron((os_y_vec-1)*n_x,ones(n_os_x,1));



% ----------------------------------------------------------------------- %
% --- READ DATA FROM PREVIOUS TIME STEP --------------------------------- %
% ----------------------------------------------------------------------- %


if guembel == 0
    g_pts_vec = sign(sign(pts_vec(1:n,1))+1);                              % switch function from the previous time step
else
    g_pts_vec = ones(n,1);                                                 % Guembel --> set switch function to 1
end

g_vec = g_pts_vec;                                                         % previous switch function is used as initial guess for current switch function
g_vec(os_vec) = g_os;   												   % ensure that the switch functions in the oil supply groove are consistent with the prescribed oil supply pressure or oil supply film fraction
g_vec_old = [];                                                            % vector for storing the switch function of the previous iteration within the current time step (does not exist yet)

if quasistatic == 1                                                        % if a quasistatic simulation is requested (meaning that all time derivatives are set to zero)
    Delta_T = (abs(u)/(2*r_b))*1e16;                                       % set the nondimensionless time increment to an incredibly large value
else                                                                       % if a transient simulation is performed
    t_pts = pts_vec(n+1,1);                                                % time at the previous time step
    Delta_T = (abs(u)/(2*r_b))*(t-t_pts);                                  % nondimensionalized time increment
end



% ----------------------------------------------------------------------- %
% --- ANALYZE GAP FUNCTION AND ITS RATE OF CHANGE ----------------------- %
% ----------------------------------------------------------------------- %


cos_vec = repmat(cos(X_vec(1:n_x,1)-X_att),[n_y,1]);                       % cos(X-X_att) at every node, where X is the angular circumferential position of the node and X_att is the attitude angle
sin_vec = repmat(sin(X_vec(1:n_x,1)-X_att),[n_y,1]);                       % sin(X-X_att) at every node, where X is the angular circumferential position of the node and X_att is the attitude angle

H_vec = 1 + ac_vec/c - epsilon*cos_vec;                                    % nondimensionalized gap function without shaft tilting [-]

dHdT_vec = (d_b/abs(u))*...                                                % nondimensionalized time derivative of the nondimensionalized gap function without shaft tilting [-]
    (-epsilon_dot*cos_vec-(epsilon*X_att_dot)*sin_vec);

if symBC == 0                                                              % if no symmetric BC is used (i.e., if the consideration of shaft tilting is enabled)
    
    cos_vec = repmat(cos(X_vec(1:n_x,1)-X_tilt),[n_y,1]);                  % cos(X-X_tilt) at every node, where X is the angular circumferential position of the node and X_tilt is an angle describing the tilting direction of the shaft
    sin_vec = repmat(sin(X_vec(1:n_x,1)-X_tilt),[n_y,1]);                  % sin(X-X_tilt) at every node, where X is the angular circumferential position of the node and X_tilt is an angle describing the tilting direction of the shaft
    
    H_vec = H_vec - (tilt*r_b/c)*Y_vec.*sin_vec;                           % nondimesionalized gap function: consideration of shaft tilting
    
    dHdT_vec = dHdT_vec + ((d_b/abs(u))*(r_b/c))*Y_vec.*...                % nondimensionalized time derivative of the nondimensionalized gap function: consideration of shaft tilting
        (-tilt_dot*sin_vec+(tilt*X_tilt_dot)*cos_vec);
    
end

H3_vec = H_vec.^3;                                                         % vector of cubed nondimensionalized nodal gap functions
H3_over_mu_vec = (H3_vec./mu_vec)*mu_ref;                                  % vector of cubed nondimensionalized nodal gap functions divided by nondimensionalized nodal viscosities



% ----------------------------------------------------------------------- %
% --- DIMENSIONLESS CONDUCTIVITIES BETWEEN NODES ------------------------ %
% ----------------------------------------------------------------------- %


aN_vec = (L_X/(24*L_Y)) * ...                                              % vector of dimensionless conductivities in the positive Y-direction (interaction with DOF N) for every DOF
    (H3_over_mu_vec(nn0_vec,1)+H3_over_mu_vec(nnN_vec,1));
aE_vec = (L_Y/(24*L_X)) * ...                                              % vector of dimensionless conductivities in the positive X-direction (interaction with DOF E) for every DOF
    (H3_over_mu_vec(nn0_vec,1)+H3_over_mu_vec(nnE_vec,1));
aS_vec = (L_X/(24*L_Y)) * ...                                              % vector of dimensionless conductivities in the negative Y-direction (interaction with DOF S) for every DOF
    (H3_over_mu_vec(nn0_vec,1)+H3_over_mu_vec(nnS_vec,1));
aW_vec = (L_Y/(24*L_X)) * ...                                              % vector of dimensionless conductivities in the negative X-direction (interaction with DOF W) for every DOF
    (H3_over_mu_vec(nn0_vec,1)+H3_over_mu_vec(nnW_vec,1));



% ----------------------------------------------------------------------- %
% --- FIXED-POINT ITERATION --------------------------------------------- %
% ----------------------------------------------------------------------- %


convergent = 0;                                                            % this will be set to 1 once the solution has converged

for iter = 1:iter_max                                                      % loop for iterative solution of nonlinear BVP
    
    
    
    % ------------------------------------------------------------------- %
    % --- ASSEMBLY ------------------------------------------------------ %
    % ------------------------------------------------------------------- %
    
    
    R_vec = (-sgn_u*L_Y/2)*(H_vec(nnE_vec,1)-H_vec(nnW_vec,1)) ...         % RHS vector
        -(L_X*L_Y)*dHdT_vec(nn0_vec,1) + ...
        (L_X*L_Y/Delta_T)*H_vec(nn0_vec,1).*(1-g_pts_vec(nn0_vec,1)).*...
        pts_vec(nn0_vec,1);
    
    if symBC == 1                                                          % if a symmetric BC is used
        K_vec = vertcat( ...                                               % the assembly of a full matrix K is avoided; all matrix entries are stored in a vector
            (aN_vec+aE_vec+aS_vec+aW_vec).*g_vec(nn0_vec,1)+...
            L_Y*H_vec(nn0_vec,1).*(1-g_vec(nn0_vec,1))+...
            (L_X*L_Y)*(dHdT_vec(nn0_vec,1)+H_vec(nn0_vec,1)/Delta_T).*...
            (1-g_vec(nn0_vec,1)), ...
            -aN_vec(1:(dof-n_x),1).*g_vec(nnN_vec(1:(dof-n_x),1),1), ...   % the DOF N vanishes when located at the bearing boundary
            -aS_vec.*g_vec(nnS_vec,1), ...
            -aE_vec.*g_vec(nnE_vec,1)+...
            ((sgn_u-1)*L_Y/2).*H_vec(nnE_vec,1).*(1-g_vec(nnE_vec,1)), ...
            -aW_vec.*g_vec(nnW_vec,1)+...
            ((-sgn_u-1)*L_Y/2).*H_vec(nnW_vec,1).*(1-g_vec(nnW_vec,1)) );
        if mod(n_y,2) == 1                                                 % if the axial number of nodes is uneven (if nodes exist at the symmetric boundary)
            R_vec(1:n_x,1) = R_vec(1:n_x,1)/2;                             % some equations contained in the equation system are divided by 2 because they correspond to half-sized control volumes
            K_vec(halfcv_vec,1) = K_vec(halfcv_vec,1)/2;                   % some equations contained in the equation system are divided by 2 because they correspond to half-sized control volumes
        end
    else                                                                   % if no symmetric BC is used
        K_vec = vertcat( ...                                               % the assembly of a full matrix K is avoided; all matrix entries are stored in a vector
            (aN_vec+aE_vec+aS_vec+aW_vec).*g_vec(nn0_vec,1)+...
            L_Y*H_vec(nn0_vec,1).*(1-g_vec(nn0_vec,1))+...
            (L_X*L_Y)*(dHdT_vec(nn0_vec,1)+H_vec(nn0_vec,1)/Delta_T).*...
            (1-g_vec(nn0_vec,1)), ...
            -aN_vec(1:(dof-n_x),1).*g_vec(nnN_vec(1:(dof-n_x),1),1), ...   % the DOF N vanishes when located at the second bearing boundary
            -aS_vec((n_x+1):dof,1).*g_vec(nnS_vec((n_x+1):dof,1),1), ...   % the DOF S vanishes when located at the first bearing boundary
            -aE_vec.*g_vec(nnE_vec,1)+...
            ((sgn_u-1)*L_Y/2).*H_vec(nnE_vec,1).*(1-g_vec(nnE_vec,1)), ...
            -aW_vec.*g_vec(nnW_vec,1)+...
            ((-sgn_u-1)*L_Y/2).*H_vec(nnW_vec,1).*(1-g_vec(nnW_vec,1)) );
    end
    
    pf = 100*max(abs(K_vec(1:dof,1)));                                     % penalty factor for enforcing the oil supply BCs
    R_vec(os_vec,1) = pf*Pi_os;                                            % oil supply BCs applied to R_vec (RHS vector): penalty method (H.R.Schwarz p.162)
    K_vec(os_vec,1) = K_vec(os_vec,1) + pf;                                % oil supply BCs applied to K_vec (vector containing the matrix entries of K): penalty method (H.R.Schwarz p.162)
    
    K_mat = sparse(rows_vec,cols_vec,K_vec,dof,dof);                       % construction of the matrix K in sparse format
    
	
    
    % ------------------------------------------------------------------- %
    % --- SOLVE EQUATION SYSTEM ----------------------------------------- %
    % ------------------------------------------------------------------- %
    
	
    Pi_vec = zeros(n,1);                                                   % this vector will be used for storing the solution (pressure-like function Pi); its length is equal to the number of nodes
    Pi_vec(nn0_vec,1) = K_mat\R_vec;                                       % the equation system is solved and the solution is stored in Pi_vec; the entries of Pi_vec corresponding to nodes that aren't part of the computational domain (because they are located at a bearing boundary or beyond the symmetric boundary, if existent) are ignored for now
    
    
    
    % ------------------------------------------------------------------- %
    % --- UPDATE SWITCH FUNCTION AND CHECK CONVERGENCE ------------------ %
    % ------------------------------------------------------------------- %
    
    
    g_vec_old_old = g_vec_old;                                             % save old switch function as old old switch function
    g_vec_old = g_vec;                                                     % save current switch function as old switch function
    g_vec = sign(sign(Pi_vec)+1);                                          % compute new switch function and save as current switch function
	
    if max(abs(g_vec(nn0_vec,1)-g_vec_old(nn0_vec,1))) == 0                % if the solution is convergent
        convergent = 1;                                                    % flag as convergent
        break                                                              % end the iterative scheme
    end
	
    if iter >= 2 && ...                                                    % if the switch function oscillates between two configurations (actually, I think this doesn't happen with the current numerical model)
            max(abs(g_vec(nn0_vec,1)-g_vec_old_old(nn0_vec,1))) == 0
        convergent = 1;                                                    % consider the solution convergent anyway
        break                                                              % end the iterative scheme
    end
    
    if guembel == 1                                                        % if Guembel is used instead of Elrod ...
        break                                                              % ... quit iterative scheme after first iteration
    end
    
end



% ----------------------------------------------------------------------- %
% --- POSTPROCESSING ---------------------------------------------------- %
% ----------------------------------------------------------------------- %


if longgroove == 1                                                         % if the groove extends across the whole axial bearing length
    Pi_vec(os_x_vec,1) = Pi_os;                                            % change atmospheric pressure to supply pressure at the nodes that are part of the oil supply groove but were treated as part of the first bearing boundary (this is only neccessary if symBC = 0; if symBC = 1, it doesn't do any harm though)
    Pi_vec(os_x_vec+(n_y-1)*n_x,1) = Pi_os;                                % change atmospheric pressure to supply pressure at the nodes that are part of the oil supply groove but were treated as part of the second bearing boundary
end

Pi_mat = reshape(Pi_vec,[n_x,n_y]);                                        % express solution as matrix

if symBC == 1                                                              % if symmetric BC was used 
    Pi_mat(:,2:(dof_y+1)) = fliplr(Pi_mat(:,j_start:(n_y-1)));             % add other half of the bearing
    Pi_vec = reshape(Pi_mat,[n,1]);                                        % define solution as vector again, now containing both halves of the bearing
end

g_vec = sign(sign(Pi_vec)+1);                                              % vector of nodal switch functions [-]
theta_vec = (1-g_vec).*Pi_vec + 1;                                         % vector of nodal film fractions [-]
p_vec = p_ref*(Pi_vec.*g_vec);                                             % vector of nodal pressures (zero is atmospheric pressure) [Pa]

if guembel == 1                                                            % if Guembel is used instead of Elrod, ...
    theta_vec = ones(n,1);                                                 % ... set film fractions to 1 and ...
end

cos_vec = repmat(cos(X_vec(1:n_x,1)),[n_y,1]);                             % vector containing the cosines of the angular circumferential positions of the nodes
sin_vec = repmat(sin(X_vec(1:n_x,1)),[n_y,1]);                             % vector containing the sines of the angular circumferential positions of the nodes

l_x = L_X*r_b;                                                             % physical (as opposed to nondimensionalized) control volume side length in the circumferential direction [m]
l_y = L_Y*r_b;                                                             % physical (as opposed to nondimensionalized) control volume side length in the axial direction [m]
A_i = l_x*l_y;                                                             % physical (as opposed to nondimensionalized) area of a regular-sized control volume [m^2]
A_vec = ones(n,1)*A_i;                                                     % vector for storing the control volume areas
A_vec(1:n_x,1) = A_i/2;                                                    % half-sized control volumes at the first bearing boundary
A_vec((1:n_x)+(n_y-1)*n_x,1) = A_i/2;                                      % half-sized control volumes at the second bearing boundary

F_1 = (cos_vec.*p_vec)'*A_vec;                                             % component of hydrodynamic force
F_2 = (sin_vec.*p_vec)'*A_vec;                                             % component of hydrodynamic force

F_h = F_1*cos(angle_shell) - F_2*sin(angle_shell);                         % transformation into inertial system --> horizontal hydrodynamic force acting on the shell [N]
F_v = F_1*sin(angle_shell) + F_2*cos(angle_shell);                         % transformation into inertial system --> vertical hydrodynamic force acting on the shell [N]

M_1 = ((sin_vec.*p_vec.*Y_vec)'*A_vec)*r_b;                                % component of hydrodynamic moment (the factor r_b converts the nondimensionalized axial coordinate values to physical ones)
M_2 = -((cos_vec.*p_vec.*Y_vec)'*A_vec)*r_b;                               % component of hydrodynamic moment (the factor r_b converts the nondimensionalized axial coordinate values to physical ones)

M_h = M_1*cos(angle_shell) - M_2*sin(angle_shell);                         % transformation into inertial system --> hydrodynamic moment around horizontal axis acting on the shell [Nm]
M_v = M_1*sin(angle_shell) + M_2*cos(angle_shell);                         % transformation into inertial system --> hydrodynamic moment around vertical axis acting on the shell [Nm]

nnE_full_vec = (1:n)'+1;                                                   % vector containing the number of the neighboring node E for every node 0
nnW_full_vec = nnE_full_vec-2;                                             % vector containing the number of the neighboring node W for every node 0
nnE_full_vec((1:n_y)*n_x,1) = nnE_full_vec((1:n_y)*n_x,1)-n_x;             % periodic BC
nnW_full_vec((0:(n_y-1))*n_x+1,1) = nnW_full_vec((0:(n_y-1))*n_x+1,1)+n_x; % periodic BC

M_fr = -((H_vec.*(p_vec(nnE_full_vec,1)-p_vec(nnW_full_vec,1)))'*A_vec)... % friction moment acting on the shell [Nm]
    *(c*r_b/(4*l_x)) ...
    + ((theta_vec.*mu_vec./H_vec)'*A_vec)*(u*r_b/c);

V_oil = c*((theta_vec.*H_vec)'*A_vec);                                     % oil volume (not equal to the overall fluid volume if cavitation occurs) [m^3]

V_dot_bb = (c^3*l_x/(12*l_y)) * ( ...                                      % oil volume flow through the bearing boundaries [m^3/s]
    (p_vec(1:n_x,1)-p_vec((n_x+1):(2*n_x),1))'*...
    (H3_vec(1:n_x,1)./mu_vec(1:n_x,1)) + ...
    (p_vec((n-n_x+1):n,1)-p_vec((n-2*n_x+1):(n-n_x),1))'*...
    (H3_vec((n-n_x+1):n,1)./mu_vec((n-n_x+1):n,1)) );



% ----------------------------------------------------------------------- %
% --- SAVE DATA FOR NEXT TIME STEP -------------------------------------- %
% ----------------------------------------------------------------------- %


pts_vec(1:n,1) = Pi_vec;                                                   % make sure that the solution (pressure-like function) of the current time step is known at the next time step
pts_vec(n+1,1) = t;                                                        % make sure that the time of the current time step is known at the next time step



end