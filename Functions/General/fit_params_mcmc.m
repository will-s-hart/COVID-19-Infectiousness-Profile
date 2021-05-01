function [theta_mat,ll_vec,acceptance_rate] = fit_params_mcmc(ll_indiv_form,data_struct_observed,no_steps,theta_prop_sd,theta_init)
    
    % Use data augmentation MCMC to fit the parameters, theta, of the
    % model of infectiousness under consideration, to the observed
    % transmission pair data contained in data_struct.
    
    % Number of fitted parameters
    no_params_fitted = length(theta_init);
    
    % Create matrices containing each possible symptom onset time for each
    % source and recipient
    t_s1_mat = cell2mat((data_struct_observed.t_s1_data)')';
    t_s2_mat = cell2mat((data_struct_observed.t_s2_data)')';
    
    % Number of transmission pairs
    no_pairs = size(t_s1_mat,1);
    
    % Length of grid of possible onset times
    t_s_grid_length = size(t_s1_mat,2);
    
    % Initialise vector of parameters, theta
    theta = theta_init;
    
    % Initialise vectors t_s1 and t_s2 containing onset times for each
    % source and recipient
    t_s1_inds = randi(t_s_grid_length,no_pairs,1);
    t_s2_inds = randi(t_s_grid_length,no_pairs,1);

    t_s1 = t_s1_mat(sub2ind([no_pairs,t_s_grid_length],(1:no_pairs)',t_s1_inds));
    t_s2 = t_s2_mat(sub2ind([no_pairs,t_s_grid_length],(1:no_pairs)',t_s2_inds));
        
    % Structure array containing initial augmented data
    data_struct_augmented = rmfield(data_struct_observed,{'t_s1_data','t_s2_data'});
    data_struct_augmented.t_s1_data = t_s1; data_struct_augmented.t_s2_data = t_s2;
    
    % Calculate initial likelihood (ll_indiv is a vector of the likelihood
    % contributions from each transmission pair)
    ll_indiv = ll_indiv_form(theta,data_struct_augmented);
    ll = sum(ll_indiv);
    
    % Create vectors to hold output of fitting procedure
    theta_mat = zeros(no_steps,no_params_fitted); %parameters at each step
    acceptance_vec = zeros(no_steps,1); %to record which steps are accepted
    ll_vec = zeros(no_steps,1); %likelihood at each step
    
    % Break steps into 100 groups to record progress
    step_no_mat = reshape(1:no_steps,no_steps/100,100);
    
    % Run chain
    for j = 1:size(step_no_mat,2)
        for i = 1:size(step_no_mat,1)

            step_no = step_no_mat(i,j);

            if mod(step_no,2)>0.5
                
                % For odd step numbers, update model parameters using
                % independent normal proposal distributions
                theta_prop = theta + randn(1,no_params_fitted).*theta_prop_sd;

                % Calculate the likelihood using the proposed parameters and
                % onset times
                ll_indiv_prop = ll_indiv_form(theta_prop,data_struct_augmented);
                ll_prop = sum(ll_indiv_prop);
                
                % Ratio between proposed and previous likelihoods
                a = exp(ll_prop-ll);

                % Accept the proposed parameters and onset times with
                % probability a
                if rand < a

                    % Update the current parameters
                    theta = theta_prop;
                    
                    ll_indiv = ll_indiv_prop;
                    ll = ll_prop;

                    acceptance_vec(step_no) = 1;
                end
            else
                
                % For even step numbers, resample the symptom onset times
                % of each host
                t_s1_inds_prop = randi(t_s_grid_length,no_pairs,1);
                t_s2_inds_prop = randi(t_s_grid_length,no_pairs,1);
             
                t_s1_prop = t_s1_mat(sub2ind([no_pairs,t_s_grid_length],(1:no_pairs)',t_s1_inds_prop));
                t_s2_prop = t_s2_mat(sub2ind([no_pairs,t_s_grid_length],(1:no_pairs)',t_s2_inds_prop));
                
                % Populate the structure array data_struct_augmented_prop
                % with the proposed onset times
                data_struct_augmented_prop = data_struct_augmented;
                data_struct_augmented_prop.t_s1_data = t_s1_prop;
                data_struct_augmented_prop.t_s2_data = t_s2_prop;

                % Calculate the likelihood contributions from each
                % transmission pair using the proposed parameters and onset
                % times
                ll_indiv_prop = ll_indiv_form(theta,data_struct_augmented_prop);
                
                % Ratio between each proposed and previous likelihood
                % contribution
                la_vec = ll_indiv_prop-ll_indiv;
                
                % Accept individual changes with probabilities given by
                % entries in exp(la_vec)
                accept_inds = log(rand(no_pairs,1))<la_vec;
                
                % Update order of transmission
                t_s1_inds(accept_inds) = t_s1_inds_prop(accept_inds);
                t_s2_inds(accept_inds) = t_s2_inds_prop(accept_inds);
                t_s1(accept_inds) = t_s1_prop(accept_inds);
                t_s2(accept_inds) = t_s2_prop(accept_inds);
                
                data_struct_augmented.t_s1_data = t_s1;
                data_struct_augmented.t_s2_data = t_s2;
                
                ll_indiv(accept_inds) = ll_indiv_prop(accept_inds);
                ll = sum(ll_indiv);
                
                acceptance_vec(step_no) = mean(accept_inds);
            end
            
            % Record parameter values and likelihood
            theta_mat(step_no,:) = theta;
            ll_vec(step_no) = ll;
        end
        
        % Display progress when an integer percentage of steps has been
        % completed
        fprintf('%d%% complete\n',100*step_no/no_steps);
    end
    
    % Calculate acceptance rates: overall, and when either the parameters
    % or the onset times were updated.
    acceptance_rate_overall = mean(acceptance_vec)
    acceptance_rate_theta = mean(acceptance_vec(1:2:end))
    acceptance_rate_data = mean(acceptance_vec(2:2:end))

    acceptance_rate.overall = acceptance_rate_overall;
    acceptance_rate.theta = acceptance_rate_theta;
    acceptance_rate.data = acceptance_rate_data;
end