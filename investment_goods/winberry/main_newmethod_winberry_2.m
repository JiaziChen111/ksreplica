%% Housekeeping
% question 1: Old question...do we really need q in the state? Can we just plug q = phi*K
%             into that?
%          2: It seems the starting point of ttheta is too low, maybe this prevents the convergence? Can we
%          restrict it in [0.5,1.0], which is sensible?
%          3: It seems there lacks decipline on ttheta, such as
%          ttheta=1/sum(active). Is that a problem?

clear; close all; clc;
deep_para; % load parameters
resale = 0.0;

%% Accuracy control
nk = 50; % number of grid points on capital stock
nfine = 200; % number of finer grid points for simulation
nx = 7; % number of grid points on idiosyncractic prod.
ns = 2*nx; % number of joint exo state
nK = 10;
nssigmak = 5; % ... on aggregate capital level
nq = 10;
m = 2.5; % support is m s.d. away from mean
tol = 1e-3;
damp = 0.5;

%% Grids
% ssigmagrid = [0.8;1.2];
% Pssigma = 0.5*ones(2,2);
[X,PX_low] = tauchen(nx,0,rrhox,ssigmax_low,m);
PX_high = tauchen_givengrid(0,rrhox,ssigmax_high,X);
uncond_X_low = PX_low^20000;
uncond_X_high = PX_high^20000;
X = exp(X);% Convert back to level
P = zeros(2*nx,2*nx);
low = 1;
high = 2;
for i = 1:2*nx
    [i_x,i_ssigmax] = ind2sub([nx 2],i);
    for j = 1:2*nx
        [j_x,j_ssigmax] = ind2sub([nx 2],j);
        P(i,j) = Pssigmax(i_ssigmax,j_ssigmax)*( (j_ssigmax==low)*PX_low(i_x,j_x)+ (j_ssigmax==high)*PX_high(i_x,j_x) );
    end
end
uncond_S = P^3000;

% Capital stuff
max_k = 2000;
min_k = 5;
% k_grid = zeros(1,nk)';
% for i_k = 1:nk
%     k_grid(end-i_k+1) = (1-ddelta)^(i_k-1)*max_k;
% end
k_grid = linspace(min_k,max_k,nk)'; % calibrated from ig_calibration_10_4


fine_grid = linspace(k_grid(1),k_grid(nk),nfine)';
noinvest_ind = ones(nk,1); % for each k, the index of tmr k if no invest
for i_k = 1:nk
    [~,noinvest_ind(i_k)] = min(abs(k_grid-(1-ddelta)*k_grid(i_k)));
end
w = 1; % Normalized to one
inv_mat = repmat(k_grid',nk,1)-(1-ddelta)*repmat(k_grid,1,nk);
pos_inv = inv_mat>0;
neg_inv = inv_mat<=0;
% K_tile = repmat(k_grid,1,ns,nK,nssigmak,nq);
min_k = max_k*(1-ddelta)^nk-1;
K_grid = linspace(min_k,max_k,nK)'; % Aggregate capital grid

% Variance grid
ssigmak_grid = linspace(sqrt(1/nfine),1500,nssigmak)';

q_grid = linspace(0.005,1.5,nq); % grid for current q

% Rule of q as function of aggregate states
pphi_qC = log(mean(q_grid));
pphi_qK = 0; % aggregate prediction rule for q
pphi_qssigmak = 0;

pphi_tthetaC = log(0.5); % aggregate prediction rule for ttheta
pphi_tthetaK = 0;
pphi_tthetassigmak = 0;

pphi_KK = 0.99; pphi_KC = log(mean(K_grid)); pphi_Kssigmak = 0.01;% Aggregate Law of motion for aggregate capital
pphi_ssigmakK = 0.1; pphi_ssigmakC = 0.1; pphi_ssigmakssigmak = 0.01;% Aggregate Law of motion for aggregate capital

%% Initialize value functions
W_old = ones(nk,ns,nK,nssigmak,nq); % value of matching with investment goods producer after paying the search cost
W_new = W_old;
U_old = ones(nk,ns,nK,nssigmak); % value of not going to search, not related to current q
U_new = U_old;
V_old = ones(nk,ns,nK,nssigmak,nq); % maximized value after discrete choice
V_new = ones(nk,ns,nK,nssigmak,nq); %
max_movingpart = zeros(nk,ns); % k,kplus
koptind_active = zeros(nk,ns,nK,nssigmak,nq);
kopt_fine = zeros(nfine,ns,nK,nssigmak,nq);
W_new_fine = zeros(nfine,ns,nK,nssigmak,nq);
U_new_fine = zeros(nfine,ns,nK,nssigmak);

profit = zeros(nk,ns); % Compute operation profit resulting from labor choice
L = zeros(nk,ns);
% for i_k = 1:nk
%     for i_x = 1:nx
%         profit(i_k,i_x) = K(i_k)*X(i_x)^(1/aalpha)*w^((aalpha-1)/aalpha)*(1-aalpha)^(1/aalpha)*(aalpha/(1-aalpha));
%     end
% end
nu = 0.5; % y = x*k^aalpha*L^nu
for i_k = 1:nk
    for i_s = 1:ns
        [i_x,~] = ind2sub([nx 2],i_s);
        L(i_k,i_s) = (X(i_x)*nu*k_grid(i_k)^aalpha/w)^(1/(1-nu));
        profit(i_k,i_s) = X(i_x)*k_grid(i_k)^aalpha*L(i_k,i_s)^nu-L(i_k,i_s)*w;
    end
end

% Prepare for Simulation stuff
T = 500;
N = 1000;
burnin = 100;
kss = 10;
Ksim = kss*ones(1,T);
qsim = mean(q_grid)*ones(1,T);
ssigmaksim = Ksim;
x_cdf_low = cumsum(PX_low,2);
x_cdf_high = cumsum(PX_high,2);
activesim = zeros(N,1);
uncond_draw = zeros(N,1);
uu = rand(N,1);
uncond_xind = zeros(N,1); uncond_x = zeros(N,1);
for i = 1:N
    uncond_xind(i) = find(cumsum(uncond_X_low(1,:),2)>=uu(i),1,'first');
    uncond_x(i) = X(uncond_xind(i));
end
revenue = zeros(nq,T);
uu = rand(N,T);
% Initialize distribution
dist_k = zeros(nfine,nx,T);
for i_x = 1:nx
    dist_k(floor(nfine/2),i_x,1) = uncond_X_low(1,i_x); % Economy starts at low volatility state
end

% New Discrete Simulation 
rng('default');
rng(2015);
markov_shock = rand(1,T);
ssigmaxsim = ones(1,T);
ssigmax_cdf = cumsum(Pssigmax,2);
for t = 1:T-1
    ssigmaxsim(t+1) = find(ssigmax_cdf(ssigmaxsim(t),:) >= markov_shock(t),1,'first');
end

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
            for i_ssigmak = 1:nssigmak
                % Predict future aggregate variable
                [Kplus,i_Kplus] = min(abs(K_grid-exp(pphi_KC+pphi_KK*log(K_grid(i_K))+pphi_Kssigmak*log(ssigmak_grid(i_ssigmak))) ));
                [qplus,i_qplus] = min(abs(q_grid-exp(pphi_qC+pphi_qK*log(K_grid(i_K))+pphi_qssigmak*log(ssigmak_grid(i_ssigmak)) ) ));
                [ssigmakplus,i_ssigmakplus] = min(abs(ssigmak_grid- exp( pphi_ssigmakC+pphi_ssigmakK*log(K_grid(i_K))+pphi_ssigmakssigmak*log(ssigmak_grid(i_ssigmak)) ) ));
                
                % Current aggregate variables
                ttheta = exp( pphi_tthetaC+pphi_tthetaK*log(K_grid(i_K))+pphi_tthetassigmak*log(ssigmak_grid(i_ssigmak)) );
                mmu = 1./(1+ttheta.^(-aalpha0)).^(1/aalpha0);
                
                % Find fucking expected value tomorrow
                EV = V_old(:,:,i_Kplus,i_ssigmakplus,i_qplus)*P';
                
                % U has no choice, so no need to find max in a loop.
                U_new(:,:,i_K,i_ssigmak) = profit + bbeta*EV(noinvest_ind,:);
                
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
                    for i_s = 1:ns
                        [max_movingpart(:,i_s),koptind_active(:,i_s,i_K,i_ssigmak,i_q)] = max(bbeta*repmat(EV(:,i_s)',nk,1) - q*inv_mat.*pos_inv - resale*q*inv_mat.*neg_inv,[],2);
                    end
                    W_new(:,:,i_K,i_ssigmak,i_q) = profit + mmu*max_movingpart - ttau + bbeta*(1-mmu)*EV(noinvest_ind,:);
                    
                    %====================================================%
                    V_new(:,:,i_K,i_ssigmak,i_q) = max(W_new(:,:,i_K,i_ssigmak,i_q),U_new(:,:,i_K,i_ssigmak));
                end
            end
        end        
            
        err = norm([V_old(:);W_old(:);U_old(:)]-[V_new(:);W_new(:);U_new(:)],Inf);
        V_old = V_new;
        W_old = W_new;
        U_old = U_new;
        iter = iter + 1;
        if mod(iter,100) == 0
            disp_text = sprintf('KS Iter = %d, KS err = %d, Current VFI Iter = %d, err = %d',outer_iter,diff,iter,err);
            disp(disp_text);
        end
    end
    
    % When converged, find policy functions
    active = W_new > repmat(U_new,1,1,1,1,nq);
    koptind = koptind_active.*active + repmat(noinvest_ind,1,ns,nK,nssigmak,nq).*(1-active);
    kopt = k_grid(koptind);
    
    % Interpolate on finer grid
    for i_s = 1:ns
        for i_K = 1:nK
            for i_ssigmak = 1:nssigmak
                U_new_fine(:,i_s,i_K,i_ssigmak) = interp1(k_grid,U_new(:,i_s,i_K,i_ssigmak),fine_grid,'linear')';
                for i_q = 1:nq
                    kopt_fine(:,i_s,i_K,i_ssigmak,i_q) = interp1(k_grid,kopt(:,i_s,i_K,i_ssigmak,i_q),fine_grid,'linear')';
                    W_new_fine(:,i_s,i_K,i_ssigmak,i_q) = interp1(k_grid,W_new(:,i_s,i_K,i_ssigmak,i_q),fine_grid,'linear')';
                end
            end
        end
    end
    active_fine = W_new_fine > repmat(U_new_fine,1,1,1,1,nq);
        
    %% Given individual policies, simulate a large panel to update aggregate law of motion
    tthetasim = exp(pphi_tthetaC+pphi_tthetaK*log(kss)++pphi_tthetassigmak*log(ssigmak_grid(1)))*ones(1,T);
    
    % I forgot to zero out after each outer iteration.
    dist_k(:,:,2:T) = zeros(nfine,nx,T-1);
    for t = 1:T
        % Find aggregate stuff today
        Ksim(t) = sum(vec(dist_k(:,:,t).*repmat(fine_grid,1,nx)));
        ssigmaksim(t) = sum(vec(dist_k(:,:,t).*(repmat(fine_grid,1,nx)-Ksim(t)).^2 ));
        
        [~,i_K] = min(abs(Ksim(t)-K_grid));
        [~,i_ssigmak] = min(abs(ssigmaksim(t)-ssigmak_grid));
        
        whichxind = (1:nx)+nx*(ssigmaxsim(t)-1);
        if ssigmaxsim(t) == low
            whichprob = PX_low;
        elseif ssigmaxsim(t) == high
            whichprob = PX_high;
        end
        
        % According to policy functions, find the optimal q
        for i_q = 1:nq
            tot_revenue_grid = q_grid(i_q)*dist_k(:,:,t).*(kopt_fine(:,whichxind,i_K,i_ssigmak,i_q)-(1-ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichxind,i_K,i_ssigmak,i_q));
            tot_revenue_grid(tot_revenue_grid<0) = resale*tot_revenue_grid(tot_revenue_grid<0);
            revenue(i_q,t) = sum(tot_revenue_grid(:));
        end
        [~,i_qmax] = max(revenue(:,t));
        qsim(t) = q_grid(i_qmax);
        
        % Evolution under the argmax q
        if (t<=T-1)
            for i_k = 1:nfine
                for i_x = 1:nx
                    i_s = sub2ind([nx 2],i_x,ssigmaxsim(t));
                    kplus = kopt_fine(i_k,i_s,i_K,i_ssigmak,i_qmax);
                    
                    % Assign mass to tomorrow's distribution
                    for i_xplus = 1:nx
                        if (kplus>fine_grid(1) && kplus<fine_grid(nfine))
                            lower_ind = find(fine_grid<=kplus,1,'last');
                            upper_ind = lower_ind + 1;
                            denom = fine_grid(upper_ind)-fine_grid(lower_ind);
                            dist_k(lower_ind,i_xplus,t+1) = dist_k(lower_ind,i_xplus,t+1) + whichprob(i_x,i_xplus)*dist_k(i_k,i_x,t)*(kplus-fine_grid(lower_ind))/denom;
                            dist_k(upper_ind,i_xplus,t+1) = dist_k(upper_ind,i_xplus,t+1) + whichprob(i_x,i_xplus)*dist_k(i_k,i_x,t)*(fine_grid(upper_ind)-kplus)/denom;
                        elseif (kplus<=fine_grid(1))
                            dist_k(1,i_xplus,t+1) = dist_k(1,i_xplus,t+1) + whichprob(i_x,i_xplus)*dist_k(i_k,i_x,t);
                        elseif (kplus>=fine_grid(nfine))
                            dist_k(nfine,i_xplus,t+1) = dist_k(nfine,i_xplus,t+1) + whichprob(i_x,i_xplus)*dist_k(i_k,i_x,t);
                        end
                    end
                    
                end
            end
        end
        
        % Importantly, fine the mass of active firms
        tthetasim(t) = sum(vec(dist_k(:,:,t).*(active_fine(:,whichxind,i_K,i_ssigmak,i_qmax)) ));
    end
    
    % update kss
    kss = mean(Ksim(burnin+1:T));
    
    % Regress to get coefficients of K law
    X = [ones(T-burnin-1,1) log(Ksim(burnin+1:T-1))' log(ssigmaksim(burnin+1:T-1))'];
    Y = log(Ksim(2+burnin:T)');
    bbeta_K = (X'*X)\(X'*Y);
    e = Y-X*bbeta_K;
    ytilde = Y-mean(Y);
    Rsq_K = 1-(e'*e)/(ytilde'*ytilde);
    pphi_KC_new = damp*bbeta_K(1)+(1-damp)*pphi_KC; pphi_KK_new = damp*bbeta_K(2)+(1-damp)*pphi_KK;
    pphi_Kssigmak_new = damp*bbeta_K(3)+(1-damp)*pphi_Kssigmak;
    
    % Regress to get ttheta law
    Y = log(tthetasim(1+burnin:T-1))';
    bbeta_ttheta = (X'*X)\(X'*Y);
    e = Y-X*bbeta_ttheta;
    ytilde = Y-mean(Y);
    Rsq_ttheta = 1-(e'*e)/(ytilde'*ytilde);      
    pphi_tthetaC_new = damp*bbeta_ttheta(1)+(1-damp)*pphi_tthetaC; pphi_tthetaK_new = damp*bbeta_ttheta(2)+(1-damp)*pphi_tthetaK;
    pphi_tthetassigmak_new = damp*bbeta_ttheta(3)+(1-damp)*pphi_tthetassigmak;
    
    % Regress to get q law
    Y = log(qsim(1+burnin:T-1))';
    bbeta_q = (X'*X)\(X'*Y);
    e = Y-X*bbeta_q;
    ytilde = Y-mean(Y);
    Rsq_q = 1-(e'*e)/(ytilde'*ytilde);      
    pphi_qC_new = damp*bbeta_q(1)+(1-damp)*pphi_qC; pphi_qK_new = damp*bbeta_q(2)+(1-damp)*pphi_qK;
    pphi_qssigmak_new = damp*bbeta_q(3)+(1-damp)*pphi_qssigmak;

    diff = norm([pphi_KC,pphi_KK,pphi_Kssigmak,pphi_tthetaC,pphi_tthetaK,pphi_tthetassigmak,pphi_qC,pphi_qK,pphi_qssigmak]-[pphi_KC_new,pphi_KK_new,pphi_Kssigmak_new,pphi_tthetaC_new,pphi_tthetaK_new,pphi_tthetassigmak_new,pphi_qC_new,pphi_qK_new,pphi_qssigmak_new],Inf);
    % Update mmu_old as well
    pphi_tthetaC = pphi_tthetaC_new; pphi_tthetaK = pphi_tthetaK_new; pphi_tthetassigmak = pphi_tthetassigmak_new;
    pphi_KC = pphi_KC_new; pphi_KK = pphi_KK_new; pphi_Kssigmak = pphi_Kssigmak_new;
    pphi_qC = pphi_qC_new; pphi_qK = pphi_qK_new; pphi_qssigmak = pphi_qssigmak_new;
            
    outer_iter = outer_iter + 1;
    disp_text = sprintf('Rsq_K = %d, Rsq_ttheta = %d, Rsq_q = %d',Rsq_K,Rsq_ttheta,Rsq_q);
    disp(disp_text);
    disp_text = sprintf('log(q) = %d + %d * log(K) + %d * log(ssigmak)',pphi_qC,pphi_qK,pphi_qssigmak);
    disp(disp_text);
    disp_text = sprintf('log(ttheta) = %d + %d * log(K)+%d * log(ssigmak)',pphi_tthetaC,pphi_tthetaK,pphi_tthetassigmak);
    disp(disp_text);
    disp_text = sprintf('log(Kplus) = %d + %d * log(K) + %d * log(ssigmak)',pphi_KC,pphi_KK,pphi_Kssigmak);
    disp(disp_text);

end
toc

save main.mat
