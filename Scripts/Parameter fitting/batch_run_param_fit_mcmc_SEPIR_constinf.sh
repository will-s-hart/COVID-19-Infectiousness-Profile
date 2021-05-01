#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=120:00:00
#SBATCH --job-name=model3_fit
#SBATCH --mail-type=ALL
#SBATCH --mail-user=william.hart@keble.ox.ac.uk
#SBATCH --output=Out/out_SEPIR_constinf_fit.out

module load matlab
matlab -nodisplay -nosplash < param_fit_mcmc_SEPIR_constinf.m > Out/log_SEPIR_constinf_fit.log