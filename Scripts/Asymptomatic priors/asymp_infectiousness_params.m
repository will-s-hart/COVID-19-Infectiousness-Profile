% Obtain the parameters of a lognormal distribution for the relative
% infectiousness of asymptomatic infected individuals, so that this
% distribution is consistent with the confidence bounds in Buitrago-Garcia
% et al., "﻿Occurrence and transmission potential of asymptomatic and
% presymptomatic SARS- CoV-2 infections: A living systematic review and
% meta-analysis", PLOS Med. (2020).

clear all; close all; clc;

% Target 95% confidence bounds
target_c025 = 0.1;
target_c975 = 1.27;

% Vectors of possible parameters, and grid of parameter pairs
mu_vec = linspace(-2,0,101);
sigma_vec = linspace(0.01,1,100);
[mu_grid,sigma_grid] = ndgrid(mu_vec,sigma_vec);

% Calculate actual cumulative density at target confidence bounds for each
% grid point
F1_grid = logncdf(target_c025,mu_grid,sigma_grid);
F2_grid = logncdf(target_c975,mu_grid,sigma_grid);

% Find parameters giving closest to desired confidence bounds
diff_grid = (F1_grid-0.025).^2 + (F2_grid-0.975).^2;
[diff_best,posn_best] = min(diff_grid(:));
mu_best = mu_grid(posn_best);
sigma_best = sigma_grid(posn_best);

% Save best parameters
save('../../Results/asymp_infectiousness_params','mu_best','sigma_best')