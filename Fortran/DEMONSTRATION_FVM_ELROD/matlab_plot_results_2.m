clear variables
dbstop if error
% close all
clc

% adjust this to your settings
n_x = 100;
n_y = 26;

% this often doesn't work well, better double-check
results1D_g = textread('results2D_g.txt');
results1D_theta = textread('results2D_theta.txt');
results2D_Pi = textread('results2D_Pi.txt');
results2D_p = textread('results2D_p.txt');

g_mat = results1D_g(1:n_x,1:n_y);
theta_mat = results1D_theta(1:n_x,1:n_y);
Pi_mat = results2D_Pi(1:n_x,1:n_y);
p_mat = results2D_p(1:n_x,1:n_y);

X_vec = linspace(0,1-1/n_x,n_x)'*360;
X_mat = repmat(X_vec,[1,n_y]);
xi_vec = linspace(-1,1,n_y);
xi_mat = repmat(xi_vec,[n_x,1]);

figure
surf(X_mat(1:n_x,1:n_y),xi_mat(1:n_x,1:n_y),g_mat(1:n_x,1:n_y))            % plot switch function
xlim([0,360])
xticks([0,90,180,270,360])
xlabel('\itX\rm (deg)')
ylabel('\it\xi\rm ()')
zlabel('\itg\rm ()')

figure
surf(X_mat(1:n_x,1:n_y),xi_mat(1:n_x,1:n_y),theta_mat(1:n_x,1:n_y))        % plot film fraction
xlim([0,360])
xticks([0,90,180,270,360])
xlabel('\itX\rm (deg)')
ylabel('\it\xi\rm ()')
zlabel('\it\vartheta\rm ()')

figure
surf(X_mat(1:n_x,1:n_y),xi_mat(1:n_x,1:n_y),Pi_mat(1:n_x,1:n_y))           % plot pressure-like function
xlim([0,360])
xticks([0,90,180,270,360])
xlabel('\itX\rm (deg)')
ylabel('\it\xi\rm ()')
zlabel('\it\Pi\rm ()')

figure
surf(X_mat(1:n_x,1:n_y),xi_mat(1:n_x,1:n_y),p_mat(1:n_x,1:n_y))            % plot physical pressure
xlim([0,360])
xticks([0,90,180,270,360])
xlabel('\itX\rm (deg)')
ylabel('\it\xi\rm ()')
zlabel('\itp\rm (Pa)')
