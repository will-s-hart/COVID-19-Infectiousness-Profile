function l_vec = log_likelihood_indiv_indep(f_inc,f_gen,data_struct_augmented)

    % Calculate the contribution from each transmission pair to the
    % log-likelihood of the augmented data, data_struct_augmented, for the
    % independent transmission and symptoms model with generation time
    % distribution f_gen and incubation period f_inc.
    
    % Extract cell arrays containing vectors of possible infection and
    % symptom onset times for each source and recipient
    t_s1_data = data_struct_augmented.t_s1_data;
    t_s2_data = data_struct_augmented.t_s2_data;
    t_i1_data = data_struct_augmented.t_i1_data;
    t_i2_data = data_struct_augmented.t_i2_data;
    
    % Time step for numerical intergration
    dt = data_struct_augmented.dt;

    % Bounds for incubation period and generation time
    t_inc_min = data_struct_augmented.t_inc_min;
    t_inc_max = data_struct_augmented.t_inc_max;
    t_gen_min = data_struct_augmented.t_gen_min;
    t_gen_max = data_struct_augmented.t_gen_max;
    
    % Vectors of discretised incubation period and generation time values
    t_inc_vec_all = ((t_inc_min+(dt/2)):dt:(t_inc_max-(dt/2)))';
    t_gen_vec_all = (t_gen_min:dt:t_gen_max)';
    
    % Evaluate incubation period and generation time distributions on the
    % vectors of possible values
    f_inc_vec_all = f_inc(t_inc_vec_all);
    f_gen_vec_all = f_gen(t_gen_vec_all);
    f_gen_vec_all(isinf(f_gen_vec_all)) = 0; %take left limit when distribution undefined at 0 (if shape parameter less than 1; right limit is infinity)

    % Create interpolation functions to enable the incubation period and
    % generation time distributions to be evaluated efficiently on
    % different subsets of the grid.
    f_inc_interp = griddedInterpolant(t_inc_vec_all,f_inc_vec_all);
    f_gen_interp = griddedInterpolant(t_gen_vec_all,f_gen_vec_all);
    
    % Evaluate the likelihood for each individual transmission pair
    t_s1_data = num2cell(t_s1_data); t_s2_data = num2cell(t_s2_data);
    L_vec = cellfun(@L_indiv,t_s1_data,t_s2_data,t_i1_data,t_i2_data);
    
    % Calculate vector of contributions to the log-likelihood from each
    % transmission pair
    l_vec = log(L_vec);

        function L = L_indiv(t_s1,t_s2,t_i1_vec,t_i2_vec)
            
            % Calculate the likelihood for a single transmission pair, for
            % known onset times t_s1 and t_s2, and vectors of possible
            % infection times t_i1_vec and t_i2_vec.
            
            % Evaluate incubation period distribution at possible
            % incubation period values for the source and recipient
            f_inc1_vec = f_inc_interp(t_s1-t_i1_vec);
            f_inc2_vec = f_inc_interp(t_s2-t_i2_vec);

            % Evaluate generation time distribution at possible values
            f_gen_grid = f_gen_interp(t_i2_vec-t_i1_vec');

            % Calculate likelihood by discretising integrals over possible
            % infection times            
            L = f_inc2_vec'*f_gen_grid*f_inc1_vec*dt*dt;
        end
end