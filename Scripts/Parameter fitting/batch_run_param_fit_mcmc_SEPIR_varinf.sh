#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=120:00:00
#SBATCH --job-name=model4_fit
#SBATCH --mail-type=ALL
#SBATCH --mail-user=william.hart@keble.ox.ac.uk
#SBATCH --output=Out/out_SEPIR_varinf_fit.out

module load matlab
matlab -nodisplay -nosplash < param_fit_mcmc_SEPIR_varinf.m > Out/log_SEPIR_varinf_fit.log