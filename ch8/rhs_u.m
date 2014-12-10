function value = rhs_u(aprime,ttau,r,a,b,coarse_grid,ve1,vu1,eeta,bbeta,P)
    u = 1; e = 2;
    c = (1+(1-ttau)*r)*a+b - aprime;
    V_e_interp = interp1(coarse_grid,ve1,aprime,'spline');
    V_u_interp = interp1(coarse_grid,vu1,aprime,'spline');
    value = (c)^(1-eeta)/(1-eeta) + bbeta*(P(u,e)*V_e_interp+P(u,u)*V_u_interp);
end