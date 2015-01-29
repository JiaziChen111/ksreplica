%% Housekeeping
clear; close all; clc;
deep_para; % load parameters

%% Accuracy control
nk = 40; % number of grid points on capital stock
nx = 7; % number of grid points on idiosyncractic prod.
nK = 40; % ... on aggregate capital level
nq = 10;
m = 3; % support is m s.d. away from mean
tol = 1e-2;
damp = 0.5;

%% Grids
% ssigmagrid = [0.8;1.2];
% Pssigma = 0.5*ones(2,2);
[X,PX] = tauchen(nx,0,rrhox,ssigmax,m);
X = exp(X);% Convert back to level
uncond_X = PX^3000;

% Capital stuff
K = linspace(5,40,nk)'; % calibrated from ig_calibration_10_4
noinvest_ind = ones(nk,1); % for each k, the index of tmr k if no invest
for i_k = 1:nk
    [~,noinvest_ind(i_k)] = min(abs(K-(1-ddelta)*K(i_k)));
end
w = 3.63; % wage, taken as given, calibrated from ig_calibration_10_4
inv_mat = repmat(K',nk,1)-(1-ddelta)*repmat(K,1,nk);
pos_inv = inv_mat>0;
neg_inv = inv_mat<=0;

q_grid = linspace(0.1,1.1,nq); % grid for current q
pphi_qC = log(0.1);
pphi_qK = 0; % aggregate prediction rule for q
q_old = ones(nK,1);

pphi_tthetaC = log(0.1); % aggregate prediction rule for ttheta
pphi_tthetaK = 0;
ttheta_old = ones(nK,1);
for i_K = 1:nK
    ttheta_old(i_K) = exp(pphi_tthetaC+pphi_tthetaK*log(K(i_K)));
end
mmu_old = 1./(1+ttheta_old.^(-aalpha0)).^(1/aalpha0);

pphi_KK = 0.99; pphi_KC = 0.1; % Aggregate Law of motion for aggregate capital

%% Initialize value functions
W_old = ones(nk,nx,nK,nq); % value of matching with investment goods producer after paying the search cost
W_new = W_old;
U_old = ones(nk,nx,nK); % value of not going to search, not related to current q
U_new = U_old;
V_old = ones(nk,nx,nK,nq); % maximized value after discrete choice
V_new = ones(nk,nx,nK,nq); %
max_movingpart = zeros(nk,nx); % k,kplus
koptind_active = zeros(nk,nx,nK,nq);

profit = zeros(nk,nx); % Compute operation profit resulting from labor choice
L = zeros(nk,nx);
% for i_k = 1:nk
%     for i_x = 1:nx
%         profit(i_k,i_x) = K(i_k)*X(i_x)^(1/aalpha)*w^((aalpha-1)/aalpha)*(1-aalpha)^(1/aalpha)*(aalpha/(1-aalpha));
%     end
% end
nu = 0.5; % y = x*k^aalpha*L^nu
for i_k = 1:nk
    for i_x = 1:nx
        L(i_k,i_x) = (X(i_x)*nu*K(i_k)^aalpha/w)^(1/(1-nu));
        profit(i_k,i_x) = X(i_x)*K(i_k)^aalpha*L(i_k,i_x)^nu-L(i_k,i_x)*w;
    end
end

% Prepare for Simulation stuff
T = 1000;
N = 1000;
burnin = 100;
kss = 10;
Ksim = kss*ones(1,T);
qsim = mean(q_grid)*ones(1,T);
x_cdf = cumsum(PX,2);
activesim = zeros(N,1);
uncond_draw = zeros(N,1);
uu = rand(N,1);
uncond_xind = zeros(N,1); uncond_x = zeros(N,1);
for i = 1:N
    uncond_xind(i) = find(cumsum(uncond_X(1,:),2)>=uu(i),1,'first');
    uncond_x(i) = X(uncond_xind(i));
end
revenue = zeros(nq,1);
uu = rand(N,T);

tic
%% Main Body of KS iter
diff = 10;
outer_iter = 0;
while diff > tol
    % Start VFI here
    err = 10;
    iter = 0;
    while err > tol
        for i_K = 1:nK
            % Predict future aggregate variable
            [Kplus,i_Kplus] = min(abs(K-exp(pphi_KC+pphi_KK*log(K(i_K)))));
            [qplus,i_qplus] = min(abs(K-exp(pphi_qC+pphi_qK*log(K(i_K)))));
            
            ttheta = ttheta_old(i_K);
            mmu = mmu_old(i_K);
            EV = V_old(:,:,i_K,i_qplus)*PX';
            
            % U has no choice, so no need to find max in a loop.
            U_new(:,:,i_K) = profit + bbeta*EV(noinvest_ind,:);
            
            for i_q = 1:nq
                q = q_grid(i_q);
                %==== Not Vectorized ================================%
                %             for i_x = 1:nx
                %                 for i_k = 1:nk
                %                     % moving_part = bbeta*EV(:,i_x) - q*K; % Only two terms in W is moving with k_t+1
                %                     moving_part = bbeta*EV(:,i_x) -q*(K-(1-ddelta)*K(i_k)).*(K>(1-ddelta)*K(i_k))+0.0*q*(K-(1-ddelta)*K(i_k)).*(K<(1-ddelta)*K(i_k));
                %                     [~,i_k_investopt] = max(moving_part);
                %                     W_new(i_k,i_x,i_K) = profit(i_k,i_x) + mmu*moving_part(i_k_investopt) - ttau + bbeta*(1-mmu)*EV(noinvest_ind(i_k),i_x);
                %                 end
                %             end
                %====================================================%
                %==== Vectorized ====================================%
                for i_x = 1:nx
                    [max_movingpart(:,i_x),koptind_active(:,i_x,i_K,i_q)] = max(bbeta*repmat(EV(:,i_x)',nk,1) - q*inv_mat.*pos_inv - 0.1*q*inv_mat.*neg_inv,[],2);
                end
                W_new(:,:,i_K,i_q) = profit + mmu*max_movingpart - ttau + bbeta*(1-mmu)*EV(noinvest_ind,:);
                
                %====================================================%
                V_new(:,:,i_K,i_q) = max(W_new(:,:,i_K,i_q),U_new(:,:,i_K));
            end
        end
        
        err = norm([V_old(:);W_old(:);U_old(:)]-[V_new(:);W_new(:);U_new(:)],Inf);
        V_old = V_new;
        W_old = W_new;
        U_old = U_new;
        iter = iter + 1;
        % disp(iter);
        % disp(err);
    end
    
    % When converged, find policy functions
    active = W_new > repmat(U_new,1,1,1,nq);
    koptind = koptind_active.*active + repmat(noinvest_ind,1,nx,nK,nq).*(1-active);
    K_tile = repmat(K,1,nx,nK,nq);
    kopt = K_tile(koptind);
    
    
    %% Given individual policies, simulate a large panel to update aggregate law of motion
    tthetasim = exp(pphi_tthetaC+pphi_tthetaK*log(kss))*ones(1,T);
    dist_k_now = kss*ones(N,1);
    dist_k_tmr = dist_k_now;
    [~,i_kss] = min(abs(kss-K));
    dist_kind_now = i_kss*ones(N,1);
    dist_kind_tmr = dist_kind_now;
    dist_xind_now = uncond_xind;
    dist_xind_tmr = dist_xind_now;

    for t = 1:T
        [~,i_K] = min(abs(Ksim(t)-K));
        
        % According to policy functions, find the optimal q
        for i_q = 1:nq
            revenue(i_q) = 0;
            for i = 1:N
                i_x = dist_xind_now(i);
                i_k = dist_kind_now(i);
                dist_kind_tmr(i) = koptind(i_k,i_x,i_K,i_q);
                dist_k_tmr(i) = K(dist_kind_tmr(i));
                activesim(i) = active(i_k,i_x,i_K,i_q);
                
                % Find revenue from this guy
                inv = dist_k_tmr(i) - (1-ddelta)*dist_k_now(i);
                revenue(i_q) = revenue(i_q) + activesim(i)*q_grid(i_q)*inv/N;
            end
        end
        [~,i_qmax] = max(revenue);
        qsim(t) = q_grid(i_qmax);
        
        % Evolution under the argmax q
        for i = 1:N
            i_x = dist_xind_now(i);
            dist_xind_tmr(i) = find(x_cdf(i_x,:)>=uu(i,t),1,'first');
            i_k = dist_kind_now(i);
            dist_kind_tmr(i) = koptind(i_k,i_x,i_K,i_qmax);
            dist_k_tmr(i) = K(dist_kind_tmr(i));
            activesim(i) = active(i_k,i_x,i_K,i_qmax);
        end
                
        if t < T
            Ksim(t+1) = sum(dist_k_tmr,1)/N;
        end
        
        % Find the ttheta, the measure of active firms as a function of agg. K
        tthetasim(t) = sum(activesim)/N;
        
        % Get ready for next period
        dist_k_now = dist_k_tmr;
        dist_xind_now = dist_xind_tmr;
        dist_kind_now = dist_kind_tmr;
        % disp(t);
        % disp(tthetasim(t));
    end
    
    % update kss
    kss = mean(Ksim(burnin+1:T));
    
    % Regress to get coefficients of K law
    X = [ones(T-burnin-1,1) log(Ksim(burnin+1:T-1))'];
    Y = log(Ksim(2+burnin:T)');
    bbeta_K = (X'*X)\(X'*Y);
    pphi_KC_new = damp*bbeta_K(1)+(1-damp)*pphi_KC; pphi_KK_new = damp*bbeta_K(2)+(1-damp)*pphi_KK;
    
    % Regress to get ttheta law
    Y = log(tthetasim(2+burnin:T)+1e-10)';
    bbeta_ttheta = (X'*X)\(X'*Y);
    pphi_tthetaC_new = damp*bbeta_ttheta(1)+(1-damp)*pphi_tthetaC; pphi_tthetaK_new = damp*bbeta_ttheta(2)+(1-damp)*pphi_tthetaK;
    
    % Regress to get q law
    Y = log(qsim(2+burnin:T))';
    bbeta_q = (X'*X)\(X'*Y);
    pphi_qC_new = damp*bbeta_q(1)+(1-damp)*pphi_qC; pphi_qK_new = damp*bbeta_q(2)+(1-damp)*pphi_qK;
    
    diff = norm([pphi_KC,pphi_KK,pphi_tthetaC,pphi_tthetaK,pphi_qC,pphi_qK]-[pphi_KC_new,pphi_KK_new,pphi_tthetaC_new,pphi_tthetaK_new,pphi_qC_new,pphi_qK_new],Inf);
    
    % Update mmu_old as well
    pphi_tthetaC = pphi_tthetaC_new; pphi_tthetaK = pphi_tthetaK_new;
    pphi_KC = pphi_KC_new; pphi_KK = pphi_KK_new;
    pphi_qC = pphi_qC_new; pphi_qK = pphi_qK_new;
    for i_K = 1:nK
        ttheta_old(i_K) = exp(pphi_tthetaC+pphi_tthetaK*log(K(i_K)));
    end
    mmu_old = 1./(1+ttheta_old.^(-aalpha0)).^(1/aalpha0);
    
    outer_iter = outer_iter + 1;
    
    disp('============');
    disp(outer_iter);
    disp(diff);
    disp('============');

end
toc

save main.mat
