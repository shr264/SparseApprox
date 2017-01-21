#!/bin/bash
#SBATCH --job-name=admm
#SBATCH -o=cscs_%A_%a.out
#SBATCH -e=cscs_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shr264@ufl.edu
#SBATCH --account=statistics
#SBATCH --qos=statistics 
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH -t 40:00:00
#SBATCH --array 1-9
[[ -d $SLURM_SUBMIT_DIR ]] & cd $SLURM_SUBMIT_DIR
date;hostname;pwd

module load R/3.2.0

Rscript admm2_${SLURM_ARRAY_TASK_ID}.R

find . -group statistics -exec chgrp presnell {} \;

date
