% Accuracy Control
temp_err = 10;
temp_tol = 1e-5;
[a_mesh,K_mesh] = meshgrid(coarse_grid,K_grid);
while temp_err>temp_tol
    % Find aggregate implied variables
    r = aalpha*(K_grid/0.92).^(aalpha-1) - ddelta;
    w = (1-aalpha)*(K_grid/0.92).^(aalpha);
    ttau = xxi*(1-aalpha)/(0.92/0.08+xxi*(1-aalpha));
    b = ttau*0.92/0.08*(K_grid/0.92).^(aalpha);
    Kprime = exp(ggamma0+ggamma1*reallog(K_grid));
    % Process Kprime to make within bound
    for j_K = 1:n_K
        if (Kprime(j_K)<min(K_grid))
            Kprime(j_K) = min(K_grid);
        elseif (Kprime(j_K)>max(K_grid))
            Kprime(j_K) = max(K_grid);
        end
    end
    
    
    cc_e = zeros(n_coarse,n_K,n_coarse);
    cc_u = cc_e;
    
    ind_e = zeros(n_coarse,n_K);
    ind_u = zeros(n_coarse,n_K);
    
    % In aprime,a format aprime on the first dimension
    for j_K = 1:n_K
        cc_e(:,j_K,:) = (1+(1-ttau)*r(j_K))*repmat(coarse_grid,n_coarse,1)-repmat(coarse_grid',1,n_coarse)+(1-ttau)*w(j_K);
        cc_u(:,j_K,:) = (1+(1-ttau)*r(j_K))*repmat(coarse_grid,n_coarse,1)-repmat(coarse_grid',1,n_coarse)+b(j_K);
    end

    cc_e(cc_e<=0) = 51709394;
    cc_u(cc_u<=0) = 51709394;
    
    util_e = (cc_e~=51709394).*((cc_e).^(1-eeta))/(1-eeta) + (cc_e==51709394).*(-9e10);
    util_u = (cc_u~=51709394).*((cc_u).^(1-eeta))/(1-eeta) + (cc_u==51709394).*(-9e10);
    
    for j_K = 1:n_K
        ve1 = qinterp2(a_mesh,K_mesh,V_e_old',coarse_grid,Kprime(j_K)*ones(1,n_coarse));
        vu1 = qinterp2(a_mesh,K_mesh,V_u_old',coarse_grid,Kprime(j_K)*ones(1,n_coarse));
        
        rhs_e = squeeze(util_e(:,j_K,:)) + bbeta*(P(e,e)*repmat(ve1',1,n_coarse)+P(e,u)*repmat(vu1',1,n_coarse));
        rhs_u = squeeze(util_u(:,j_K,:)) + bbeta*(P(u,e)*repmat(ve1',1,n_coarse)+P(u,u)*repmat(vu1',1,n_coarse));
        
        [junkV_e,junkind_e] = max(rhs_e);
        [junkV_u,junkind_u] = max(rhs_u);
        V_u_new(:,j_K) = junkV_u'; ind_u(:,j_K) = junkind_u';
        V_e_new(:,j_K) = junkV_u'; ind_e(:,j_K) = junkind_e';
    end
    
    temp_err = norm([V_e_new(:)-V_e_old(:);V_e_new(:)-V_e_old(:)],2);
    V_e_old = V_e_new;
    V_u_old = V_u_new;
    
%     for j_K = 1:n_K
% 
%         
%         

%         for i_a = 1:n_coarse
%             a = coarse_grid(i_a);
%             max_a_e = (1+(1-ttau)*r)*a+(1-ttau)*w-1e-4;
%             max_a_u = (1+(1-ttau)*r)*a+b-1e-4;
%             maximand_e = @(aprime) rhs_e(aprime,ttau,r,a,w,coarse_grid,ve1,vu1,eeta,bbeta,P);
%             maximand_u = @(aprime) rhs_u(aprime,ttau,r,a,b,coarse_grid,ve1,vu1,eeta,bbeta,P);
%             [V_e_new(i_a,j_K),aopt_e(i_a,j_K)] = gss(-2,max_a_e,maximand_e);
%             [V_u_new(i_a,j_K),aopt_u(i_a,j_K)] = gss(-2,max_a_u,maximand_u);
%         end
%     end
end
disp('Individual Policies Found.');

