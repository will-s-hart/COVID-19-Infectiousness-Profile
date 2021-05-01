function b_cond = b_cond_form_ferretti(t_tost,t_inc1,params)

    % Infectiousness of a host at time t_tost since symptom onset,
    % conditional on incubation period t_inc1, under the Ferretti model
    % with parameters given by params.

    if length(t_inc1)==1
        t_inc1 = repmat(t_inc1,size(t_tost));
    elseif length(t_tost)==1
        t_tost = repmat(t_tost,size(t_inc1));
    end

    mu = params(1);
    sigma = params(2);
    alpha = params(3);
    m_inc = params(4);

    C = alpha/(sigma*(1-(1+exp((m_inc+mu)/sigma))^(-alpha))); %constant of proportionality to give correct scaling for infectiousness

    f_m = C*exp(alpha*((t_tost./(t_inc1/m_inc))-mu)/sigma)./(1+exp(((t_tost./(t_inc1/m_inc))-mu)/sigma)).^(alpha+1);
    f_p = C*exp(-(t_tost-mu)/sigma)./(1+exp(-(t_tost-mu)/sigma)).^(alpha+1);

    ind_m = (t_tost<0);
    ind_p = logical(1-ind_m);

    b_cond = zeros(size(t_tost));
    b_cond(ind_m) = f_m(ind_m);
    b_cond(ind_p) = f_p(ind_p);

    b_cond(t_tost<-t_inc1) = 0; %truncate distribution to avoid non-zero infectiousness before infection

end