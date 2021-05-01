function F_tost_vec = F_tost_varinf_vec(t_tost,theta_mat,k_inc,gamma,k_I)

    % Cumulative TOST distribution function for our mechanistic approach at
    % the scalar time t_tost, for each set of estimated parameters given by
    % the rows of theta_mat.
    
    % Vectors of values of estimated parameters
    k_E_vec = exp(theta_mat(:,1));
    k_P_vec = k_inc - k_E_vec;
    mu_vec = exp(theta_mat(:,2));
    alpha_vec = exp(theta_mat(:,3));
    
    % Values of TOST to integrate over
    t_min = t_tost-100;
    dt = 0.1;
    t_vec = (t_min+dt/2):dt:t_tost;
    
    % Grids of pairs of values of parameter values and the TOST
    [k_P_grid,t_grid] = ndgrid(k_P_vec,t_vec);
    [mu_grid,~] = ndgrid(mu_vec,t_vec);
    [alpha_grid,~] = ndgrid(alpha_vec,t_vec);
    
    % Evaluate TOST density function on the grid
    C_grid = k_inc*gamma*mu_grid./(alpha_grid.*k_P_grid.*mu_grid+k_inc*gamma);
    
    fm = alpha_grid.*C_grid.*(1-gamcdf(-t_grid,k_P_grid,1/(k_inc*gamma)));
    fp = C_grid.*(1-gamcdf(t_grid,k_I,1./(k_I*mu_grid)));
    
    indm = (t_grid<0);
    ind0 = (t_grid==0);
    indp = logical(1-(indm+ind0));
    
    f_tost_grid = zeros(size(t_grid));
    f_tost_grid(indm) = fm(indm);
    f_tost_grid(indp) = fp(indp);
    f_tost_grid(ind0) = fp(ind0);
    
    % Intergrate density function over time values to obtain cumulative
    % distribution function at the input time t_tost (for each parameter
    % set)
    F_tost_vec = sum(f_tost_grid,2)*dt;
end