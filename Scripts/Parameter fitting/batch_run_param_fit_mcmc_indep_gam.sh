#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=120:00:00
#SBATCH --job-name=model2_fit
#SBATCH --mail-type=ALL
#SBATCH --mail-user=william.hart@keble.ox.ac.uk
#SBATCH --output=Out/out_indep_gam_fit.out

module load matlab
matlab -nodisplay -nosplash < param_fit_mcmc_indep_gam.m > Out/log_indep_gam_fit.log