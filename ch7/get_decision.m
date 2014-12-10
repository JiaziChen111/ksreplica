%Initialization
temp_err = 10;
temp_tol = 1e-5;

if existence ~= 1
    V_old = zeros(2,n_coarse);
    V_new = V_old;
    exisitence = 1;
end

cc_e = (1+(1-ttau)*r)*repmat(coarse_grid,n_coarse,1)-repmat(coarse_grid',1,n_coarse)+(1-ttau)*w;
cc_u = (1+(1-ttau)*r)*repmat(coarse_grid,n_coarse,1)-repmat(coarse_grid',1,n_coarse)+b;

cc_e(cc_e<=0) = 51709394; 
cc_u(cc_u<=0) = 51709394; 

util_e = (cc_e~=51709394).*((cc_e).^(1-eeta))/(1-eeta) + (cc_e==51709394).*(-9e10);
util_u = (cc_u~=51709394).*((cc_u).^(1-eeta))/(1-eeta) + (cc_u==51709394).*(-9e10);

while temp_err>temp_tol
    rhs_e = util_e + bbeta*(P(e,e)*repmat(V_old(e,:)',1,n_coarse)+P(e,u)*repmat(V_old(u,:)',1,n_coarse));
    rhs_u = util_u + bbeta*(P(u,e)*repmat(V_old(e,:)',1,n_coarse)+P(u,u)*repmat(V_old(u,:)',1,n_coarse));
    
    [V_new(e,:),ind_e] = max(rhs_e);
    [V_new(u,:),ind_u] = max(rhs_u);

    temp_err = norm(abs(V_new(:)-V_old(:)),Inf);
    V_old = V_new;
end
disp('Individual Policies Found.');

