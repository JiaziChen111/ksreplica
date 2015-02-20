%% Housekeeping 
clear 
close all
clc

%% Parameter (see Heer et al)
aalpha = 0.36;
eeta = 2.0;
bbeta = 0.995;
P = [0.5,0.5;0.0435,0.9565];
u = 1; e = 2;
b = 1.199;
ddelta = 0.005;

%% Accuracy control
tol = 1e-5;
n_coarse = 250;
n_fine = 100;

%% Main Boday
% Step 1
N = 0.92; % stationary employment rate
kss = (aalpha/(1/bbeta-1+ddelta))^(1/(1-aalpha))*N;
fine_grid = linspace(-2,1500,n_fine);
coarse_grid = linspace(-2,1500,n_coarse);

% Step 2 - 8
err = 10;
K = kss; ttau = 0.017;
existence = 0;
iter = 0;
while err > tol
    % Step 3
    w = (1-aalpha)*(K/N)^aalpha;
    r = aalpha*(K/N)^(aalpha-1)-ddelta;
    
    % Step 4
    get_decision;
    
    % Step 5, 7.22
    algor_724;
    
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
