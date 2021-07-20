#!/bin/bash

# SLURM batch file to run the R script `simstudy-good_leverage.R` on an HPC system.

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G

## Load R/4
module load intel
module load r/4

## Ensure OpenMP threads are "close"
export OMP_PROC_BIND="CLOSE"

Rscript --no-save "simstudy-good_leverage.R" \
  --ncores ${SLURM_CPUS_PER_TASK} \
  --job ${SLURM_ARRAY_TASK_ID}
