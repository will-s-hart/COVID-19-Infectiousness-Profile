function params = get_params_constinf(theta,k_inc,gamma,k_I)

    % Recover the full vector of parameters for our mechanistic approach,
    % params = [gamma,mu,k_inc,k_E,k_I,alpha], from the vector of
    % parameters updated in the MCMC model fitting procedure for the
    % constant infectiousness model, theta = [log(k_E),log(mu)].

    alpha = 1;

    k_E = exp(theta(1));
    mu = exp(theta(2));

    params = [gamma,mu,k_inc,k_E,k_I,alpha];

end