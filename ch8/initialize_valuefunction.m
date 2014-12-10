% Initial Guess for Value functions
ttau = 0.015;
V_e_old = zeros(n_coarse,n_K);
V_u_old = zeros(n_coarse,n_K);
for j_K = 1:n_K
    r = aalpha*(K_grid(j_K)/0.92)^(aalpha-1) - ddelta;
    w = (1-aalpha)*(K_grid(j_K)/0.92)^(aalpha);
    b = xxi*(1-ttau)*(1-aalpha)*(K_grid(j_K)/0.92)^(-aalpha);
    for i_a = 1:n_coarse
        V_e_old(i_a,j_K) = ((1-ttau)*(r*coarse_grid(i_a)+w))^(1-eeta)/(1-eeta)/(1-bbeta);
        V_u_old(i_a,j_K) = ((1-ttau)*(r*coarse_grid(i_a))+b)^(1-eeta)/(1-eeta)/(1-bbeta);
    end
end
V_e_new = V_e_old;
V_u_new = V_u_old;