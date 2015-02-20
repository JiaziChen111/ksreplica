%% Housekeeping 
clear 
close all
clc

%% Parameter (see Heer et al)
aalpha = 0.36;
eeta = 2.0;
bbeta = 0.995;
xxi = 0.25;
ddelta = 0.005;

w = 4.7;
ttau = 0.015; % as an initial guess
b = xxi*(1-ttau)*w;

P = [0.5,0.5;0.0435,0.9565];
u = 1; e = 2;

%% Accuracy control
tol = 1e-3;
n_coarse = 301;
n_fine = 603;
n_K = 8;

%% Main Boday
% Step 0
fine_grid = linspace(-2,3000,n_fine);
coarse_grid = linspace(-2,3000,n_coarse);
K_grid = linspace(140,340,n_K);
[a_mesh,K_mesh] = meshgrid(coarse_grid,K_grid);

% Step 1 Initiaze capital distribution
[~,i_right] = min(abs(fine_grid-300));
f_old = zeros(1,n_fine);
f_old(1:i_right) = 1/i_right;
F_old = cumsum(f_old);
K_old= sum(f_old.*fine_grid);

ggamma0 = 0.09; ggamma1 = 0.98;

% Step 2 - 8
err = 10;
existence = 0;
iter = 0;
while err > tol
    % Step 3

    % Step 4
    initialize_valuefunction;
    get_decision;
    
    % Step 5, 7.22
    simulate_forward;
    
    % Regress to update ggamma
    
    
    % Maybe normalization here to guarantee sums to one
    
    % Step 6-7 Update new guess
    K_new = sum(a_dist_new)/n_hh;
    ttau_new = b*(1-N)/(w*N+r*K);
    
    % Error and Update
    iter = iter+1;
    disp(iter);
    err = norm([K;ttau]-[K_new;ttau_new],Inf);
    K = K_new;
    ttau = ttau_new;
    disp(err);
    disp(K_new);
    disp('==========');
end
