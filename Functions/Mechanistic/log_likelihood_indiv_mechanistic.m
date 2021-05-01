function l_vec = log_likelihood_indiv_mechanistic(f_inc,b_cond,data_struct_augmented)

    % Calculate the contribution from each transmission pair to the
    % log-likelihood of the augmented data, data_struct_augmented, for our
    % mechanistic approach with conditional infectiousness b_cond and
    % incubation period distribution f_inc.

    % Extract cell arrays containing vectors of possible infection and
    % symptom onset times for each source and recipient
    t_s1_data = data_struct_augmented.t_s1_data;
    t_s2_data = data_struct_augmented.t_s2_data;
    t_i1_data = data_struct_augmented.t_i1_data;
    t_i2_data = data_struct_augmented.t_i2_data;
    
    % Time step for numerical intergration
    dt = data_struct_augmented.dt;
    
    % Bounds for incubation period and TOST
    t_inc_min = data_struct_augmented.t_inc_min;
    t_inc_max = data_struct_augmented.t_inc_max;
    t_tost_min = data_struct_augmented.t_tost_min;
    t_tost_max = data_struct_augmented.t_tost_max;
    
    % Vectors of discretised incubation period and TOST values, and grids
    % containing each pair of values in these vectors
    t_inc_vec_all = ((t_inc_min+(dt/2)):dt:(t_inc_max-(dt/2)))';
    t_tost_vec_all = ((t_tost_min+(dt/2)):dt:(t_tost_max-(dt/2)))';
    [t_tost_grid_all,t_inc_grid_all] = ndgrid(t_tost_vec_all,t_inc_vec_all);
    
    % Evaluate incubation period distribution on the vector of possible
    % values
    f_inc_vec_all = f_inc(t_inc_vec_all);
    
    % Evaluate conditional infectiousness on the grid of possible TOST and
    % incubation period values
    b_cond_grid_all = b_cond(t_tost_grid_all,t_inc_grid_all);
    
    % Create interpolation functions to enable the incubation
    % period distribution and conditional infectiousness to be evaluated
    % efficiently on different subsets of the grid.
    f_inc_interp = griddedInterpolant(t_inc_vec_all,f_inc_vec_all);
    b_cond_interp = griddedInterpolant(t_tost_grid_all,t_inc_grid_all,b_cond_grid_all);
    
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

            % Evaluate conditional infectiousness at possible values of the
            % TOST and of the source's incubation period
            [t_tost_grid,t_inc1_grid] = ndgrid(t_i2_vec-t_s1,t_s1-t_i1_vec);
            b_cond_mat = b_cond_interp(t_tost_grid,t_inc1_grid);
            
            % Calculate likelihood by discretising integrals over possible
            % infection times
            L = f_inc2_vec'*b_cond_mat*f_inc1_vec*dt*dt;
        end
end