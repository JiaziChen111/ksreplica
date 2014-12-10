% Algorithm 7.22 to compute the stationary distribution

% Interpolate the decision on finer grid
aprime_e = interp1(coarse_grid,coarse_grid(ind_e),fine_grid);
aprime_u = interp1(coarse_grid,coarse_grid(ind_u),fine_grid);

% Fit interpolated decision on the grid
[~,optind_e] = min(abs(repmat(aprime_e,n_fine,1)-repmat(fine_grid',1,n_fine)));
[~,optind_u] = min(abs(repmat(aprime_u,n_fine,1)-repmat(fine_grid',1,n_fine)));
aprime_e = fine_grid(optind_e);
aprime_u = fine_grid(optind_u);
for i_aprime = 1:n_fine
	a_inverse_e(i_aprime) = find(optind_e==i_aprime,1,'last');
	a_inverse_u(i_aprime) = find(optind_u==i_aprime,1,'last');
end

% Find the stationary distribution
[~,i_kss] = min(abs(fine_grid-kss));
F_old = zeros(2,n_fine);
F_old(u,i_kss:end) = 0.08;
F_old(e,i_kss:end) = 0.92;
F_new = F_old;
dist_err = 10; dist_tol = 1e-5;
while dist_err > dist_tol
    for i_aprime = 1:n_fine
        F_new(u,i_aprime) = P(u,u)*F_old(u,a_inverse_u(i_aprime)) + P(e,u)*F_old(e,a_inverse_e(i_aprime)); 
        F_new(e,i_aprime) = P(u,e)*F_old(u,a_inverse_u(i_aprime)) + P(e,e)*F_old(e,a_inverse_e(i_aprime)); 
    end
    dist_err = norm(F_old(:)-F_new(:),Inf);
    F_old = F_new;
end
