function value = rhs_e(aprime,ttau,r,a,w,coarse_grid,ve1,vu1,eeta,bbeta,P)
    u = 1; e = 2;
    c = (1+(1-ttau)*r)*a+(1-ttau)*w - aprime;
    V_e_interp = interp1(coarse_grid,ve1,aprime,'spline');
    V_u_interp = interp1(coarse_grid,vu1,aprime,'spline');
    value = (c)^(1-eeta)/(1-eeta) + bbeta*(P(e,e)*V_e_interp+P(e,u)*V_u_interp);
end