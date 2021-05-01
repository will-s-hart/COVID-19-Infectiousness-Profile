function f_gen_cheb = get_gen_dist_mechanistic(params)

    % Calculate the generation time distribution for our mechanistic
    % approach with parameters given by params.
    
    % This function requires the Chebfun package to execute (freely
    % available at https://www.chebfun.org/download/).

    gamma = params(1); mu = params(2);
    k_inc = params(3); k_E = params(4); k_I = params(5);
    alpha = params(6);
    k_P = k_inc-k_E;

    C = k_inc*gamma*mu/(alpha*k_P*mu+k_inc*gamma);

    tmax = 500;

    f_E = @(t)gampdf(t,k_E,1/(k_inc*gamma));
    f_E_cheb = chebfun(f_E,[0,10,50,tmax]);
    f_P = @(t)gampdf(t,k_P,1/(k_inc*gamma));
    f_P_cheb = chebfun(f_P,[0,10,50,tmax]);
    F_P = @(t)gamcdf(t,k_P,1/(k_inc*gamma));
    F_P_cheb = chebfun(F_P,[0,10,50,tmax]);
    F_I = @(t)gamcdf(t,k_I,1/(k_I*mu));
    F_I_cheb = chebfun(F_I,[0,10,50,tmax]);

    f1 = conv(f_P_cheb,F_I_cheb,'same','old');
    f_star_cheb = alpha*C*(1-F_P_cheb)+C*(F_P_cheb-f1);
    f_star_cheb = chebfun(f_star_cheb,[0,10,50,tmax]);

    f_gen_cheb = conv(f_E_cheb,f_star_cheb,'same','old');

end