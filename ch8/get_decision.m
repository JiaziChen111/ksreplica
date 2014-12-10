% Accuracy Control
temp_err = 10;
temp_tol = 1e-5;
while temp_err>temp_tol
    for j_K = 1:n_K
        r = aalpha*(K_grid(j_K)/0.92)^(aalpha-1) - ddelta;
        w = (1-aalpha)*(K_grid(j_K)/0.92)^(aalpha);
        ttau = xxi*(1-aalpha)/(0.92/0.08+xxi*(1-aalpha));
        b = ttau*0.92/0.08*(K_grid(j_K)/0.92)^(aalpha);        
        Kprime = exp(ggamma0+ggamma1*reallog(K_grid(j_K)));
        ve1 = interp2(a_mesh,K_mesh,V_e_old',coarse_grid,Kprime*ones(1,n_coarse),'spline');
        vu1 = interp2(a_mesh,K_mesh,V_u_old',coarse_grid,Kprime*ones(1,n_coarse),'spline');
        for i_a = 1:n_coarse
            a = coarse_grid(i_a);
            max_a_e = (1+(1-ttau)*r)*a+(1-ttau)*w-1e-4;
            max_a_u = (1+(1-ttau)*r)*a+b-1e-4;
            maximand_e = @(aprime) rhs_e(aprime,ttau,r,a,w,coarse_grid,ve1,vu1,eeta,bbeta,P);
            maximand_u = @(aprime) rhs_u(aprime,ttau,r,a,b,coarse_grid,ve1,vu1,eeta,bbeta,P);
            [V_e_new(i_a,j_K),aopt(i_a,j_K)] = gss(-2,max_a_e,maximand_e);
            [V_u_new(i_a,j_K),aopt(i_a,j_K)] = gss(-2,max_a_u,maximand_u);
        end
    end

    temp_err = norm([V_e_new;V_u_new]-[V_e_old;V_e_old],Inf)
    V_e_old = V_e_new;
    V_u_old = V_u_new;
end
disp('Individual Policies Found.');

