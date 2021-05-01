function params = get_params_ferretti(theta,m_inc)

    % Recover the full vector of parameters for the Ferretti model, params
    % = [mu,sigma,alpha,m_inc], from the vector of parameters updated in
    % the MCMC model fitting procedure, theta = [mu,log(sigma),log(alpha)].

    mu = theta(1);
    sigma = exp(theta(2));
    alpha = exp(theta(3));

    params = [mu,sigma,alpha,m_inc];

end