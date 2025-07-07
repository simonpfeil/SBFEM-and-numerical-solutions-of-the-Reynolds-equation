% SINGLE CALL OF THE REYNOLDS EQUATION

% Author: Simon Pfeil
% Affiliation: OvGU Magdeburg, Germany (until 30.06.2025)
% Date: 30.06.2025

% This script performs a single call of the Reynolds equation for
% demonstration and testing purposes. For clarification of the input and 
% output variables, the solution technique, the assumptions, and the 
% coordinate system, check the information provided in 'fvm_elrod.m' and
% the attached PDF.

clear variables
dbstop if error
% close all
clc

% ----------------------------------------------------------------------- %
% --- DEFINITION OF PHYSICAL AND NUMERICAL INPUT VARIABLES -------------- %
% ----------------------------------------------------------------------- %

d_b = 0.1;                                                                 % bearing diameter [m]
l_b = 0.08;                                                                % bearing length [m]
c = 150e-6;                                                                % radial clearance [m]
grooves = 1;                                                               % number of oil supply grooves (evenly distributed across the circumference), 1 minimum [-]
X_os = 0.25*2*pi;                                                          % angular circumferential position of one of the oil supply grooves in the reference frame of the shell [rad]
L_X_os = 15/360*2*pi;                                                      % angular circumferential side length of the oil supply groove(s) [rad]
l_y_os = 0.06;                                                             % axial side length of the oil supply groove(s) [m]
p_os = 70000;                                                              % boundary value prescribed in the oil supply groove; either a non-negative pressure in [Pa] (zero corresponds to atmospheric pressure) or a negative value stating a film fraction minus 1 [-] (for example, a film fraction of 0.7 is prescribed by p_os = 0.7-1 = -0.3)
t = 0;                                                                     % time corresponding to the current time step (irrelevant for this quasistatic single call) [s]
angle_shell = 0;                                                           % rotation angle of the shell in the reference frame of the inertial system [rad]
omega_shell = 0;                                                           % angular velocity of the shaft in the reference frame of the inertial system [rad/s] (omega_shaft must not be equal to omega_shell)
dis_h_shell = 0;                                                           % horizontal displacement of the shell in the reference frame of the inertial system [m]
dis_v_shell = 0;                                                           % vertical displacement of the shell in the reference frame of the inertial system [m]
vel_h_shell = 0;                                                           % horizontal velocity of the shell in the reference frame of the inertial system [m/s]
vel_v_shell = 0;                                                           % vertical velocity of the shell in the reference frame of the inertial system [m/s]
omega_shaft = 3000/60*2*pi;                                                % angular velocity of the shaft in the reference frame of the inertial system [rad/s] (omega_shaft must not be equal to omega_shell)
dis_h_shaft = c/2;                                                         % horizontal displacement of the shaft in the reference frame of the inertial system [m]
dis_v_shaft = 0;                                                           % vertical displacement of the shaft in the reference frame of the inertial system [m]
vel_h_shaft = 0;                                                           % horizontal velocity of the shaft in the reference frame of the inertial system [m/s]
vel_v_shaft = 0;                                                           % vertical velocity of the shaft in the reference frame of the inertial system [m/s]
tilt_h_shell = 0;                                                          % tilting angle of shell around horizontal axis [rad]
tilt_v_shell = 0;                                                          % tilting angle of shell around vertical axis [rad]
tilt_dot_h_shell = 0;                                                      % rate of change of tilting angle of shell around horizontal axis [rad/s]
tilt_dot_v_shell = 0;                                                      % rate of change of tilting angle of shell around vertical axis [rad/s]
tilt_h_shaft = 0;                                                          % tilting angle of shaft around horizontal axis [rad]
tilt_v_shaft = 0;                                                          % tilting angle of shaft around vertical axis [rad]
tilt_dot_h_shaft = 0;                                                      % rate of change of tilting angle of shaft around horizontal axis [rad/s]
tilt_dot_v_shaft = 0;                                                      % rate of change of tilting angle of shaft around vertical axis [rad/s]
symBC = 1;                                                                 % flag for the usage of a symmetric boundary condition: 0 = no, 1 = yes; setting symBC = 1 reduces the computational effort but also renders the model incapable of considering shaft tilting [-]
n_x = 100;                                                                 % circumferential number of nodes (not counting the periodic node at X=2*pi) [-]
iter_max = 100;                                                            % max. allowed number of iterations [-]
quasistatic = 1;                                                           % quasistatic = 1 means that a quasistatic simulation is requested (all velocities and angular velocities except for omega_shaft are internally set to zero and the film fraction is assumed to be constant over time); when incorporating the SBFEM solution into a time integration scheme, set quasistatic = 0 [-]
n_y = round((l_b/(pi*d_b))*n_x+1);                                         % axial number of nodes (whole bearing, including the nodes at the bearing boundaries), at least 5 [-]
ac_vec = zeros(n_x*n_y,1);                                                 % additional contour of the shell defined at the nodes, following the overall node numbering scheme; positive values increase the gap width, negative values reduce the gap width, the shell is cylindrical if all entries are zero; since the nodes are fixed at the reference frame of the shell, a rotation of the shell does not require any adjustment of ac_vec [m] 
mu_vec = 0.01*ones(n_x*n_y,1);                                             % oil viscosities prescribed at the nodes, following the overall node numbering scheme; since the nodes are fixed at the possibly rotating shell, this viscosity distribution is always expressed in the shell's reference frame [Pa*s]
guembel = 0;                                                               % use Guembel instead of Elrod cavitation? 0 = no (i.e. use Elrod), 1 = yes (i.e. use Guembel)

% The next variable, pts_vec, requires special attention. In a time 
% integration, the program fvm_elrod uses this array to communicate with 
% itself across time steps. At every call, fvm_elrod reads the data 
% contained in pts_vec and stores new data in pts_vec. When calling 
% fvm_elrod, pts_vec must contain the values that were stored in pts_vec 
% by fvm_elrod at the previous valid time step. At the first time step, 
% we use pts_vec to prescribe an initial condition instead. Setting 
% pts_vec = zeros(n_x*n_y+1,1) and then pts_vec(n_x*n_y+1,1) = -1e-5 will 
% lead to the following initial condition: a film fraction of 1 (i.e. only 
% oil, no gas cavities) is assumed at all nodes at t=-1e-5 (exactly at t=0
% isn't allowed, so we define it shortly before that). If Guembel 
% conditions are assumed, fvm_elrod pts_vec is ignored. If a quasistatic
% solution with Elrod cavitation is performed, fvm_elrod uses pts_vec only
% to choose an initial guess for the cavitation states. For transient 
% simulations with Elrod, this initial guess is also used, but more
% importantly, pts_vec allows fvm_elrod to formulate the transient 
% cavitation term (which requires knowing the film fractions of the 
% previous time step).

pts_vec = zeros(n_x*n_y+1,1);
pts_vec(n_x*n_y+1,1) = -1e-5;

% ----------------------------------------------------------------------- %
% --- CALL FVM ALGORITHM FOR SOLVING THE REYNOLDS EQUATION -------------- %
% ----------------------------------------------------------------------- %

[F_h,F_v,M_h,M_v,M_fr,V_oil,V_dot_bb,Pi_mat,pts_vec,p_ref,...
    convergent,iter] = fvm_elrod(d_b,l_b,c,grooves,X_os,L_X_os,l_y_os,...
    p_os,ac_vec,t,angle_shell,omega_shell,dis_h_shell,dis_v_shell,...
    vel_h_shell,vel_v_shell,omega_shaft,dis_h_shaft,dis_v_shaft,...
    vel_h_shaft,vel_v_shaft,tilt_h_shell,tilt_v_shell,tilt_dot_h_shell,...
    tilt_dot_v_shell,tilt_h_shaft,tilt_v_shaft,tilt_dot_h_shaft,...
    tilt_dot_v_shaft,n_x,iter_max,mu_vec,n_y,pts_vec,quasistatic,...
    symBC,guembel);

% ----------------------------------------------------------------------- %
% --- VISUALIZATION OF THE RESULTS -------------------------------------- %
% ----------------------------------------------------------------------- %

X_vec = linspace(0,1-1/n_x,n_x)'*360;
xi_vec = linspace(-1,1,n_y);
X_mat = repmat(X_vec,[1,n_y]);
xi_mat = repmat(xi_vec,[n_x,1]);
g_mat = sign(sign(Pi_mat)+1);                                              % nodal switch functions (0 = cavitation zone, 1 = pressure zone) [-]
theta_mat = (1-g_mat).*Pi_mat+1;                                           % nodal film fractions [-]
p_mat = (Pi_mat.*g_mat)*p_ref;                                             % nodal pressures, where zero is atmospheric pressure [Pa]
% figure
% surf(X_mat,xi_mat,Pi_mat)                                                  % surface plot of pressure-like function
% xlim([0,360])
% xticks([0,90,180,270,360])
% xlabel('\itX\rm (deg)')
% ylabel('\it\xi\rm ()')
% zlabel('\it\Pi\rm ()')
figure
surf(X_mat,xi_mat,p_mat)                                                   % surface plot of pressure
xlim([0,360])
xticks([0,90,180,270,360])
xlabel('\itX\rm (deg)')
ylabel('\it\xi\rm ()')
zlabel('\itp\rm (Pa)')
% figure
% surf(X_mat,xi_mat,g_mat)                                                   % surface plot of switch function
% xlim([0,360])
% xticks([0,90,180,270,360])
% xlabel('\itX\rm (deg)')
% ylabel('\it\xi\rm ()')
% zlabel('\itg\rm ()')
% figure
% surf(X_mat,xi_mat,theta_mat)                                               % surface plot of film fraction
% xlim([0,360])
% zlim([0,1])
% xticks([0,90,180,270,360])
% xlabel('\itX\rm (deg)')
% ylabel('\it\xi\rm ()')
% zlabel('\it\vartheta\rm ()')

clearvars -except F_h F_v M_h M_v M_fr V_oil V_dot_bb Pi_mat pts_vec ...
    p_ref convergent iter
