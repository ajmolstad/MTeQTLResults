#!/bin/bash

#SBATCH -o Results/Run_%a.Rout
#SBATCH --array=1-7500
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 96:00:00

export OMP_NUM_THREADS=1

module load R/3.6

R CMD BATCH --vanilla Simulations_Main.R  Results/Run_${SLURM_ARRAY_TASK_ID}.Rout
