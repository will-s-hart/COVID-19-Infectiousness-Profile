% Obtain the parameters of a beta distribution for the proportion of
% asymptomatic cases, so that this distribution is consistent with the
% confidence bounds in Buitrago-Garcia et al., "﻿Occurrence and
% transmission potential of asymptomatic and presymptomatic SARS- CoV-2
% infections: A living systematic review and meta-analysis", PLOS Med.
% (2020).

clear all; close all; clc;

% Target 95% confidence bounds
target_c025 = 0.26;
target_c975 = 0.37;

% Vectors of possible parameters, and grid of parameter pairs
a_vec = linspace(1,200,200);
b_vec = linspace(1,200,200);
[a_grid,b_grid] = ndgrid(a_vec,b_vec);

% Calculate actual cumulative density at target confidence bounds for each
% grid point
F1_grid = betacdf(target_c025,a_grid,b_grid);
F2_grid = betacdf(target_c975,a_grid,b_grid);

% Find parameters giving closest to desired confidence bounds
diff_grid = (F1_grid-0.025).^2 + (F2_grid-0.975).^2;
[diff_best,posn_best] = min(diff_grid(:));
a_best = a_grid(posn_best);
b_best = b_grid(posn_best);

% Save best parameters
save('../../Results/asymp_proportion_params','a_best','b_best')