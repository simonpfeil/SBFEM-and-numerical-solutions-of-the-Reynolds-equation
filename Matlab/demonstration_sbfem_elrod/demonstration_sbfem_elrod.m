% SINGLE CALL OF THE REYNOLDS EQUATION

% Author: Simon Pfeil
% Affiliation: OvGU Magdeburg, Germany (until 30.06.2025)
% Date: 30.06.2025

% This script performs a single call of the Reynolds equation for
% demonstration and testing purposes. For clarification of the input and 
% output variables, the solution technique, the assumptions, and the 
% coordinate system, check the information provided in 'sbfem_elrod.m' and
% the attached PDF.

clear variables
dbstop if error
% close all
clc

% Let's first define the parameters that should stay constant throughout 
% all calls of the Reynolds equation within a time integration. Although 
% we are not performing a time integration in the exemplary code at hand, 
% I prefer already making this distinction between the parameters here.

d_b = 0.1;                                                                 % bearing diameter [m]
l_b = 0.08;                                                                % bearing length [m]
grooves = 1;                                                               % number of oil supply grooves; the grooves are assumed to be distributed equidistantly along the circumference; at least 1 groove must be present under Elrod conditions (guembel = 0), and at least 0 grooves must be present under Guembel conditions (guembel = 1) [-]
X_os = 0.25*2*pi;                                                          % angular circumferential position of one of the oil supply grooves in the reference frame of the shell [rad]
L_X_os = 15/360*2*pi;                                                      % angular circumferential width of the oil supply groove(s) [rad]
quasistatic = 1;                                                           % quasistatic = 1 means that a quasistatic simulation is requested (all velocities and angular velocities except for omega_shaft are internally set to zero and the film fraction is assumed to be constant over time); when incorporating the SBFEM solution into a time integration scheme, set quasistatic = 0 [-]
guembel = 0;                                                               % use Guembel instead of Elrod cavitation? 0 = no (i.e. use Elrod), 1 = yes (i.e. use Guembel)
n_x = 100;                                                                 % circumferential number of nodes (not counting the periodic node at X=2*pi) [-]
n_y = round((l_b/(pi*d_b))*n_x+1);                                         % only for visualization of the solution: defines the number of points in the axial direction (set to 1 to skip this optional postprocessing step and save computational effort) [-]
iter_max = 100;                                                            % max. allowed number of iterations [-]
tay = 0;                                                                   % use Taylor approximations for eigenvalues and eigenvectors? 0 = no, 1 = yes; these approximations are only possible if the additional contour ac_vec is zero, if the viscosity distribution mu_vec is spatially constant, and if n_x = 100, grooves = 0, and guembel = 1, which means that choosing tay = 1 will cause the program to manipulate these parameters if necessary [-]

% If we want the SBFEM algorithm to solve its eigenvalue problem by Taylor 
% approximations (tay = 1) instead of using an eigensolver, we need to 
% load the database of Taylor coefficients. The Taylor coefficients in the 
% given database were computed under certain assumtions, including Guembel
% conditions, no oil supply grooves, a pre-defined circumferential number 
% of nodes (n_x = 100), and a pre-defined order of the Taylor series 
% (n_tay = 2). In the following, if tay = 1, the database is loaded and 
% the parameters are adjusted to these assumptions.

if tay == 1                                                                % if Taylor approximations will be used for solving the eigenvalue problem
    fileID1 = fopen('tay_parameters.bin');                                 % open the file that states the parameters which were chosen when computing the Taylor coefficients and ...
    parameters_vec = fread(fileID1,[6,1],'double');                        % ... save the data in a vector parameters_vec. This data includes ...
    n_x = parameters_vec(1,1);                                             % ... the circumferential number of nodes (which may overwrite the one chosen above), ...
    n_tay = parameters_vec(2,1);                                           % ... the order of the Taylor series, ...
    gamma_ref = parameters_vec(3,1);                                       % ... the original slenderness ratio (required for adjusting the eigenvalues to the current slenderness ratio, which is achived via a simple scaling operation), ...
    n_points = parameters_vec(4,1);                                        % ... the number of points where Taylor series were constructed, ...
    n_Ld_all = parameters_vec(5,1);                                        % ... the size of the vector containing the eigenvalue derivatives (i.e. containing the Taylor coefficients for computing the eigenvalues), and ...
    epsilon_max = parameters_vec(6,1);                                     % ... the maximum relative eccentricity for which Taylor approximations are available in the given database
    guembel = 1;                                                           % the Taylor approximations require Guembel conditions, so even if Elrod has been chosen (guembel = 0), we resort to assuming Guembel conditions now
    grooves = 0;                                                           % the Taylor approximations require a bearing without oil supply grooves, so if a nonzero number of grooves has been chosen, we change it back to zero
    fileID2 = fopen('tay_constrpoints.bin');                               % open the file that states the points (i.e. the relative eccentricities) where Taylor series were constructed
    fileID3 = fopen('tay_switchpoints.bin');                               % open the file that states the points (i.e. the relative eccentricities) where the algorithm is supposed to switch between two neighboring Taylor series
    fileID4 = fopen('tay_reductions.bin');                                 % open the file that states, for every Taylor series, the number of considered modes of the pressure field
    fileID5 = fopen('tay_eigenvalues.bin');                                % open the file that contains the eigenvalue derivatives (i.e. the Taylor coefficients for computing the eigenvalues)
    fileID6 = fopen('tay_eigenvectors.bin');                               % open the file that contains the eigenvector derivatives (i.e. the Taylor coefficients for computing the eigenvectors)
    fileID7 = fopen('tay_indices.bin');                                    % open the file that states, for every Taylor series, where within the data we can find the Taylor coefficients for this series
    constr_vec = fread(fileID2,[n_points,1],'double');                     % vector stating the points (i.e. the relative eccentricities) where Taylor series were constructed
    switch_vec = fread(fileID3,[n_points,1],'double');                     % vector stating the points (i.e. the relative eccentricities) where the algorithm is supposed to switch between two neighboring Taylor series
    red_vec = fread(fileID4,[n_points,1],'double');                        % vector stating, for every Taylor series, the number of considered modes of the pressure field
    Ld_allpoints_vec = fread(fileID5,[1,n_Ld_all],'double');               % vector containing the eigenvalue derivatives (i.e. the Taylor coefficients for computing the eigenvalues)
    Ld_allpoints_vec = Ld_allpoints_vec*((l_b/d_b)/gamma_ref)^2;           % scale eigenvalues and their derivatives according to the defined bearing dimensions
    Vd_allpoints_mat = fread(fileID6,[n_x,n_Ld_all],'double');             % matrix containing the eigenvector derivatives (i.e. the Taylor coefficients for computing the eigenvectors)
    ind_vec = fread(fileID7,[n_points,1],'double');                        % vector stating, for every Taylor series, where within Ld_allpoints_vec and Vd_allpoints_mat we can find the corresponding Taylor coefficients
    fclose(fileID1);                                                       % close database files ...
    fclose(fileID2);                                                       % ...
    fclose(fileID3);
    fclose(fileID4);
    fclose(fileID5);
    fclose(fileID6);
    fclose(fileID7);
else
    n_tay = 0;                                                             % no Taylor approximations are used, so the order of Taylor series will be ignored, we may set it to zero (or to any other number)
    epsilon_max = 0;                                                       % no Taylor approximations are used, so the maximum eccentricity for which Taylor approximations are available will be ignored, we may set it to zero (or to any other number)
end

% If we were to perform a time integration, the time integration would 
% start here. The variables defined in the following may vary between 
% calls of the Reynolds equation.

ac_vec = zeros(n_x,1);                                                     % additional contour of the shell defined point-wise at the circumferential positions where the nodes are located; positive values increase the gap width, negative values reduce the gap width, the shell is cylindrical if all entries are zero; since the nodes are fixed at the reference frame of the shell, a rotation of the shell does not require any adjustment of ac_vec [m] 
mu_vec = 0.01*ones(n_x,1);                                                 % oil viscosity distribution defined point-wise at the circumferential positions where the nodes are located; since the nodes are fixed at the reference frame of the possibly rotating shell, this viscosity distribution is always expressed in the shell's reference frame [Pa*s]
t = 0;                                                                     % time corresponding to the current time step (irrelevant for this quasistatic single call) [s]
c = 150e-6;                                                                % radial clearance [m]
p_os = 70000;                                                              % boundary value prescribed in the oil supply groove; either a non-negative pressure in [Pa] (zero corresponds to atmospheric pressure) or a negative value stating a film fraction minus 1 [-] (for example, a film fraction of 0.7 is prescribed by p_os = 0.7-1 = -0.3)
angle_shell = 0;                                                           % rotation angle of the shell in the reference frame of the inertial system [rad]
omega_shell = 0;                                                           % angular velocity of the shaft in the reference frame of the inertial system [rad/s] (omega_shaft must not be equal to omega_shell)
dis_h_shell = 0;                                                           % horizontal displacement of the shell in the reference frame of the inertial system [m]
dis_v_shell = 0;                                                           % vertical displacement of the shell in the reference frame of the inertial system [m]
vel_h_shell = 0;                                                           % horizontal velocity of the shell in the reference frame of the inertial system [m/s]
vel_v_shell = 0;                                                           % vertical velocity of the shell in the reference frame of the inertial system [m/s]
omega_shaft = (3000/60)*2*pi;                                              % angular velocity of the shaft in the reference frame of the inertial system [rad/s] (omega_shaft must not be equal to omega_shell)
dis_h_shaft = c/2;                                                         % horizontal displacement of the shaft in the reference frame of the inertial system [m]
dis_v_shaft = 0;                                                           % vertical displacement of the shaft in the reference frame of the inertial system [m]
vel_h_shaft = 0;                                                           % horizontal velocity of the shaft in the reference frame of the inertial system [m/s]
vel_v_shaft = 0;                                                           % vertical velocity of the shaft in the reference frame of the inertial system [m/s]

% The next variable, pts_vec, requires special attention. In a time 
% integration, the program sbfem_elrod uses this array to communicate with 
% itself across time steps. At every call, sbfem_elrod reads the data 
% contained in pts_vec and stores new data in pts_vec. When calling 
% sbfem_elrod, pts_vec must contain the values that were stored in pts_vec 
% by sbfem_elrod at the previous valid time step. At the first time step, 
% we use pts_vec to prescribe an initial condition instead. Setting 
% pts_vec = zeros(n_x+1,1) and then pts_vec(n_x+1,1) = -1e-5 will lead to 
% the following initial condition: a film fraction of 1 (i.e. only oil, no 
% gas cavities) is assumed at all nodes at t=-1e-5 (exactly at t=0 isn't 
% allowed, so we define it shortly before that). If Guembel conditions are 
% assumed, sbfem_elrod ignores pts_vec. If a quasistatic solution with 
% Elrod cavitation is performed, sbfem_elrod uses pts_vec only
% to choose an initial guess for the cavitation states. For transient 
% simulations with Elrod, this initial guess is also used, but more
% importantly, pts_vec allows sbfem_elrod to formulate the transient 
% cavitation term (which requires knowing the film fractions of the 
% previous time step).

pts_vec = zeros(n_x+1,1);
pts_vec(n_x+1,1) = -1e-5;

% Now the parameters are defined, but there is still one thing to take 
% care of before solving the Reynolds equation. As mentioned above, the 
% SBFEM algorithm offer the option to tackle its eigenvalue problem by 
% Taylor approximations instead of using an eigensolver. This requires 
% setting the variable  tay = 1 (see above). If we have chosen this 
% option, we now need to extract the Taylor coefficients that are suitable 
% for the current shaft position (for the current relative eccentricity) 
% from the database (this is done in the following). This step must be 
% executed before every call of the Reynolds  equation. The Taylor 
% coefficients in the given database were computed under certain 
% assumtions, including the simplifications that the additional contour 
% ac_vec is zero and that the viscosity distribution described by mu_vec 
% is spatially constant. If tay = 1, these conditions must be enforced 
% (this is also done in the following).

if tay == 1
    epsilon = sqrt((dis_h_shaft-dis_h_shell)^2+...                         % relative eccentricity
        (dis_v_shaft-dis_v_shell)^2)/c;
    [~,k] = min(abs(switch_vec-epsilon));                                  % switch_vec contains the epsilon-positions where we switch from one Taylor series to another; this command line determines which of these switching points our current epsilon is closest to
    if epsilon < switch_vec(k,1)                                           % depending on whether we are above or below this switching point, ...
        k = k-1;                                                           % ... we choose the proper Taylor series (the k-th Taylor series in our database)
    end
    epsilon_constr = constr_vec(k,1);                                      % point where Taylor series was constructed
    red = red_vec(k,1);                                                    % number of modes to consider
    ind = ind_vec(k,1);                                                    % ind will help us find the requested Taylor coefficients in the arrays Vd_mat and Ld_vec
    Vd_mat = Vd_allpoints_mat(:,ind:(ind-1+red*(n_tay+1)));                % matrix of eigenvector derivatives (Taylor coefficients for approximating the eigenvectors) for the concrete Taylor series that will be used at the next call of the SBFEM algorithm
    Ld_vec = Ld_allpoints_vec(1,ind:(ind-1+red*(n_tay+1)));                % vector of eigenvalue derivatives (Taylor coefficients for approximating the eigenvalues) for the concrete Taylor series that will be used at the next call of the SBFEM algorithm
    ac_vec = zeros(n_x,1);                                                 % the Taylor approximations assume that there is no additional contour, so we set ac_vec to zero
    mu_vec = mean(mu_vec)*ones(n_x,1);                                     % the Taylor approximations assume that the viscosity distribution described by mu_vec is spatially constant, so we replace the viscosities by their average
else
    epsilon_constr = 0;                                                    % no Taylor approximations are used, so the point of construction will be ignored, we may set it to zero (or to any other number)
    red = 0;                                                               % no Taylor approximations are used, so the number of modes will be ignored, we may set it to zero (or to any other number)
    Vd_mat = 0;                                                            % no Taylor approximations are used, so the variable representing the eigenvector derivatives will be ignored, we may set it to zero (or to any other number)
    Ld_vec = 0;                                                            % no Taylor approximations are used, so the variable representing the eigenvalue derivatives will be ignored, we may set it to zero (or to any other number)
end

% We call the SBFEM algorithm to solve the Reynolds equation.

[F_h,F_v,M_fr,V_oil,V_dot_bb,pts_vec,g_vec,theta_vec,Pi_mat,...
    p_ref,convergent,iter] ...
    = sbfem_elrod(d_b,l_b,c,grooves,X_os,L_X_os,p_os,t,angle_shell,...
    omega_shell,dis_h_shell,dis_v_shell,vel_h_shell,vel_v_shell,...
    omega_shaft,dis_h_shaft,dis_v_shaft,vel_h_shaft,vel_v_shaft,...
    n_x,iter_max,ac_vec,mu_vec,pts_vec,guembel,quasistatic,n_y,...
    tay,n_tay,epsilon_max,epsilon_constr,red,Vd_mat,Ld_vec);

% If we want, we can visualize some results.

if n_y > 1 
    n_x = length(g_vec);                                                   % circumferential number of nodes
    p_mat = (p_ref*g_vec).*Pi_mat;                                         % pressures (zero is atmospheric pressure) [Pa]
    g_mat = repmat(g_vec,[1,n_y]);                                         % switch functions (0 = cavitation zone, 1 = pressure zone) [-]
    theta_mat = repmat(theta_vec,[1,n_y]);                                 % film fractions [-]
    X_vec = linspace(0,1-1/n_x,n_x)'*360;                                  % vector of circumferential angular coordinate values X of all nodes
    xi_vec = linspace(-1,1,n_y);                                           % vector of nondimensionalized (from -1 to 1) axial coordinate values xi of the axial data points
    X_mat = repmat(X_vec,[1,n_y]);                                         % matrix containing the X-coordinates of all nodes and data points
    xi_mat = repmat(xi_vec,[n_x,1]);                                       % matrix containing the xi-coordinates of all nodes and data points
    % figure
    % surf(X_mat,xi_mat,Pi_mat)                                              % surface plot of pressure-like function
    % xlim([0,360])
    % xticks([0,90,180,270,360])
    % xlabel('\itX\rm (deg)')
    % ylabel('\it\xi\rm ()')
    % zlabel('\it\Pi\rm ()')
    figure
    surf(X_mat,xi_mat,p_mat)                                               % surface plot of pressure
    xlim([0,360])
    xticks([0,90,180,270,360])
    xlabel('\itX\rm (deg)')
    ylabel('\it\xi\rm ()')
    zlabel('\itp\rm (Pa)')
    % figure
    % surf(X_mat,xi_mat,g_mat)                                               % surface plot of switch function
    % xlim([0,360])
    % xticks([0,90,180,270,360])
    % xlabel('\itX\rm (deg)')
    % ylabel('\it\xi\rm ()')
    % zlabel('\itg\rm ()')
    % figure
    % surf(X_mat,xi_mat,theta_mat)                                           % surface plot of film fraction
    % xlim([0,360])
    % zlim([0,1])
    % xticks([0,90,180,270,360])
    % xlabel('\itX\rm (deg)')
    % ylabel('\it\xi\rm ()')
    % zlabel('\it\vartheta\rm ()')
end

clearvars -except F_h F_v M_fr V_oil V_dot_bb pts_vec g_vec theta_vec ...
    Pi_mat p_ref convergent iter
