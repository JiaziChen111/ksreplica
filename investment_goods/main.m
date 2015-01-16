%% Housekeeping
clear; close all; clc;
deep_para; % load parameters

%% Accuracy control 
nk = 50; % number of grid points on capital stock
nx = 7; % number of grid points on idiosyncractic prod.
nK = 50; % ... on aggregate capital level
m = 3; % support is m s.d. away from mean
nq = 30;
tol = 1e-5;

%% Grids
% ssigmagrid = [0.8;1.2];
% Pssigma = 0.5*ones(2,2);
[X,PX] = tauchen(nx,0,rrhox,ssigmax,m);
X = exp(X); % Convert back to level
K = linspace(0.5*10,2*10,nk)'; % calibrated from ig_calibration_10_4
w = 3.63; % wage, taken as given, calibrated from ig_calibration_10_4

%% Pretend everything else is solved
q_old = ones(nK,1);
q_grid = linspace(0.1,10,nq);
ttheta_old = 0.1*ones(nK,1);
mmu_old = 1./(1+ttheta_old.^(-aalpha0)).^(1/aalpha0);
pphi_KK = 0.99; pphi_KC = 0.1; % Aggregate Law of motion for aggregate capital

%% Initialize value functions
W_old = ones(nk,nx,nK); % value of matching with investment goods producer after paying the search cost
W_new = W_old;
U_old = ones(nk,nx,nK); % value of not going to search 
U_new = U_old;
V_old = ones(nk,nx,nK); % maximized value after discrete choice
V_new = ones(nk,nx,nK); % 

profit = zeros(nk,nx); % Compute operation profit resulting from labor choice
for i_k = 1:nk
    for i_x = 1:nx
        profit(i_k,i_x) = K(i_k)*X(i_x)^(1/aalpha)*w^((aalpha-1)/aalpha)*(1-aalpha)^(1/aalpha)*(aalpha/(1-aalpha));
    end
end

err = 10;
iter = 0;
while err > tol
    for i_K = 1:nK
        q = q_old(i_K);
        ttheta = ttheta_old(i_K);
        mmu = mmu_old(i_K);
        [Kplus,i_Kplus] = min(abs(K-exp(pphi_KC+pphi_KK*log(K(i_k)))));
        EV = V_old(:,:,i_K)*PX';
        
        for i_x = 1:nx
            for i_k = 1:nk
                [~,i_k_noinvest] = min(abs(K-(1-ddelta)*K(i_k)));
                U_new(i_k,i_x,i_K) = profit(i_k,i_x) + bbeta*EV(i_k_noinvest,i_x);
                moving_part = bbeta*EV(:,i_x) - q*K; % Only two terms in W is moving with k_t+1
                [~,i_k_investopt] = max(moving_part);
                W_new(i_k,i_x,i_K) = profit(i_k,i_x) -mmu*q*(K(i_k_investopt)-(1-ddelta)*K(i_k)) - ttau + bbeta*(mmu*EV(i_k_investopt,i_x)+(1-mmu)*EV(i_k_noinvest,i_x));
            end
        end
        
        V_new(:,:,i_K) = max(W_new(:,:,i_K),U_new(:,:,i_K));
    end
    err = norm([V_old(:);W_old(:);U_old(:)]-[V_new(:);W_new(:);U_new(:)],Inf);
    V_old = V_new;
    W_old = W_new;
    U_old = U_new;
    iter = iter + 1;
    disp(iter);
    disp(err);
end

%% When converged, find policy functions
kopt = zeros(nk,nx,nK);
koptind = zeros(nk,nx,nK);
active = zeros(nk,nx,nK);

for i_K = 1:nK
    q = q_old(i_K);
    ttheta = ttheta_old(i_K);
    mmu = mmu_old(i_K);
    [Kplus,i_Kplus] = min(abs(K-exp(pphi_KC+pphi_KK*log(K(i_k)))));
    EV = V_old(:,:,i_K)*PX';
    
    for i_x = 1:nx
        for i_k = 1:nk
            [~,i_k_noinvest] = min(abs(K-(1-ddelta)*K(i_k)));
            U_new(i_k,i_x,i_K) = profit(i_k,i_x) + bbeta*EV(i_k_noinvest,i_x);
            moving_part = bbeta*EV(:,i_x) - q*K; % Only two terms in W is moving with k_t+1
            [~,i_k_investopt] = max(moving_part);
            W_new(i_k,i_x,i_K) = profit(i_k,i_x) -mmu*q*(K(i_k_investopt)-(1-ddelta)*K(i_k)) - ttau + bbeta*(mmu*EV(i_k_investopt,i_x)+(1-mmu)*EV(i_k_noinvest,i_x));
            if W_new(i_k,i_x,i_K) > U_new(i_k,i_x,i_K)
                kopt(i_k,i_x,i_K) = K(i_k_investopt);
                koptind(i_k,i_x,i_K) = i_k_investopt;
                active(i_k,i_x,i_K) = 1;
            else
                kopt(i_k,i_x,i_K) = K(i_k_noinvest);
                koptind(i_k,i_x,i_K) = i_k_noinvest;
                active(i_k,i_x,i_K) = 0;
            end
        end
    end    
end

%% Once we have the policy, find the optimal q policy assuming buyers believe "same old" happens tomorow and so on
revenue = zeros(nq,1);
for i_K = 1:nK
    for i_q = 1:nq
        q = q_grid(i_q);
        ttheta = ttheta_old(i_K);
        mmu = mmu_old(i_K);
        [Kplus,i_Kplus] = min(abs(K-exp(pphi_KC+pphi_KK*log(K(i_k)))));
        EV = V_old(:,:,i_K)*PX';
        
        % only q changed compared to question. Assuming uniform
        % distribution
        revenue(i_q) = 0;
        for i_x = 1:nx
            for i_k = 1:nk
                [~,i_k_noinvest] = min(abs(K-(1-ddelta)*K(i_k)));
                U_new(i_k,i_x,i_K) = profit(i_k,i_x) + bbeta*EV(i_k_noinvest,i_x);
                moving_part = bbeta*EV(:,i_x) - q*K; % Only two terms in W is moving with k_t+1
                [~,i_k_investopt] = max(moving_part);
                W_new(i_k,i_x,i_K) = profit(i_k,i_x) -mmu*q*(K(i_k_investopt)-(1-ddelta)*K(i_k)) - ttau + bbeta*(mmu*EV(i_k_investopt,i_x)+(1-mmu)*EV(i_k_noinvest,i_x));
                if W_new(i_k,i_x,i_K) > U_new(i_k,i_x,i_K)
                    kopt(i_k,i_x,i_K) = K(i_k_investopt);
                    koptind(i_k,i_x,i_K) = i_k_investopt;
                    active(i_k,i_x,i_K) = 1;
                else
                    kopt(i_k,i_x,i_K) = K(i_k_noinvest);
                    koptind(i_k,i_x,i_K) = i_k_noinvest;
                    active(i_k,i_x,i_K) = 0;
                end
                
                % Compute revenue 
                revenue(i_q) = revenue(i_q) + active(i_k,i_x,i_K)*q*(kopt(i_k,i_x,i_K)-(1-ddelta)*K(i_k));
            end
        end

    end
    [~,i_q_opt] = max(revenue);
    q_old(i_K) = q_grid(i_q_opt);
end

save main.mat

%% Given individual policies, simulate a large panel to update aggregate law of motion
T = 1000;
N = 1000;
kss = 10;
Ksim = kss*ones(1,T);

dist_k_now = kss*ones(N,1);
dist_k_tmr = dist_k_now;
[~,i_kss] = min(abs(kss-K));
dist_kind_now = i_kss*ones(N,1);
dist_kind_tmr = dist_kind_now;
dist_xind_now = ceil(nx/2)*ones(N,1);
dist_xind_tmr = ceil(nx/2)*ones(N,1);
x_cdf = cumsum(PX,2);

for t = 1:T
    [~,i_K] = min(abs(Ksim(t)-K));
    % According to policy functions, find the next period state (x,k)
    uu = rand(N,1);
    for i = 1:N
        i_x = dist_xind_now(i);
        dist_xind_tmr(i) = find(x_cdf(i_x,:)>=uu(i),1,'first');
        i_k = dist_kind_now(i);
        dist_kind_tmr(i) = koptind(i_k,i_x,i_K);
        dist_k_tmr(i) = K(dist_kind_tmr(i));
    end
    if t < T
        Ksim(t+1) = sum(dist_k_tmr,1)/N;
    end
    
    % Find the optimal q this period
    
    
    
    % Get ready for next period
    dist_k_now = dist_k_tmr;
    dist_xind_now = dist_xind_tmr;
    dist_kind_now = dist_kind_tmr;
    disp(t);
end

% Regress to get coefficients
X = [ones(T-1,1) log(Ksim(1:T-1))'];
Y = [log(Ksim(2:T)')];
bbeta = (X'*X)\(X'*Y);
pphi_KC = bbeta(1); pphi_KK = bbeta(2);

