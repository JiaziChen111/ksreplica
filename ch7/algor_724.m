n_hh = 1e5;
periods = 1000;
%shocks_e = rand(n_hh,periods);
%shocks_u = rand(n_hh,periods);
epsilon_dist_old = zeros(n_hh,1);
a_dist_old = kss*ones(n_hh,1);
epsilon_dist_new = epsilon_dist_old;
a_dist_new = a_dist_old;
epsilon_dist_old(1:n_hh*0.08,1) = 1;
epsilon_dist_old(n_hh*0.08+1:end,1) = 2;

moments_diff = 10;
toltol = 1e-3;
T_sim = 0;
while moments_diff>toltol
    % Simulate employment status
    shocks = rand(n_hh,1);
    u_to_u = (shocks < P(u,u)) .* (epsilon_dist_old == u);
    u_to_e = (shocks > P(u,u)) .* (epsilon_dist_old == u);
    e_to_u = (shocks < P(e,u)) .* (epsilon_dist_old == e);
    e_to_e = (shocks > P(e,u)) .* (epsilon_dist_old == e);
    epsilon_dist_new(logical(u_to_u+e_to_u)) = u;
    epsilon_dist_new(logical(u_to_e+e_to_e)) = e;
    
    % Maybe rescaling to maintain 0.08 unemployment rate
    discrepancy = sum(epsilon_dist_new==u) - 0.08*n_hh;
    if discrepancy > 0
        become_e = find(epsilon_dist_new==u,abs(discrepancy),'last');
        epsilon_dist_new(become_e) = e;
    elseif discrepancy < 0
        become_u = find(epsilon_dist_new==e,abs(discrepancy),'last');
        epsilon_dist_new(become_u) = u;
    end
    
    % Determine next period's asset
    a_dist_new(epsilon_dist_old==u) = interp1(coarse_grid,coarse_grid(ind_u),a_dist_old(epsilon_dist_old==u));
    a_dist_new(epsilon_dist_old==e) = interp1(coarse_grid,coarse_grid(ind_e),a_dist_old(epsilon_dist_old==e));
    
    % Compare moments
    stat_new = [mean(a_dist_new);std(a_dist_new)];
    stat_old = [mean(a_dist_old);std(a_dist_old)];
    moments_diff = norm(stat_new-stat_old,2)/norm(stat_old,2);
                
    % Update
    epsilon_dist_old = epsilon_dist_new;
    a_dist_old = a_dist_new;
    T_sim = T_sim + 1;
end;
