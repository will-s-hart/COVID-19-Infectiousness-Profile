function F_tost_vec = F_tost_indep_vec(t_tost,gen_shape_vec,gen_scale_vec,inc_shape,inc_scale)

    % Numerically calculate the cumulative TOST distribution function for
    % the independent transmission and symptoms model at scalar time
    % t_tost, for each set of generation time parameters given by the
    % entries in gen_shape_vec and gen_scale_vec, by integrating the
    % generation time distribution weighted by the proportion of hosts who
    % have developed symptoms.

    % Grid of times since infection to integrate over
    tau_max = 50;
    dt = 0.1;
    tau_vec = ((dt/2):dt:(tau_max-dt/2))';

    % Grids of pairs of values of the shape/scale parameter of the generation
    % time distribution, and of the time since infection
    [gen_shape_grid,tau_grid] = ndgrid(gen_shape_vec,tau_vec);
    [gen_scale_grid,~] = ndgrid(gen_scale_vec,tau_vec);

    % Evaluate the density of the generation time and the cumulative
    % distribution of the incubation period, at each grid point
    F_inc_vec = gamcdf(tau_vec-t_tost,inc_shape,inc_scale);
    f_gen_grid = gampdf(tau_grid,gen_shape_grid,gen_scale_grid);

    % Numerically integrate over the generation time distribution, weighted by
    % the proportion of hosts who have developed symptoms, to calculate the
    % TOST density for each parameter set.
    F_tost_vec = 1-f_gen_grid*F_inc_vec*dt;

end