function b_cond = b_cond_form_mechanistic(t_tost,t_inc1,params)
    
    % Infectiousness of a host at time t_tost since symptom onset,
    % conditional on incubation period t_inc1, under our mechanistic
    % approach with parameters given by params.
    
    if length(t_inc1)==1
        t_inc1 = repmat(t_inc1,size(t_tost));
    elseif length(t_tost)==1
        t_tost = repmat(t_tost,size(t_inc1));
    end

    gamma = params(1); mu = params(2);
    k_inc = params(3); k_E = params(4); k_I = params(5);
    alpha = params(6);
    k_P = k_inc-k_E;

    C = k_inc*gamma*mu/(alpha*k_P*mu+k_inc*gamma);
    
    fm = alpha*C*(1-betacdf(-t_tost./t_inc1,k_P,k_E));
    fp = C*(1-gamcdf(t_tost,k_I,1/(k_I*mu)));

    indm = (t_tost<0);
    ind0 = (t_tost==0);
    indp = logical(1-(indm+ind0));

    b_cond = zeros(size(t_tost));
    b_cond(indm) = fm(indm);
    b_cond(indp) = fp(indp);
    b_cond(ind0) = fp(ind0);

end