#!/bin/bash

# SLURM batch file to run the R script `simstudy-scenario_3.R` on an HPC system.

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G

## Load R/4
module load gnu10 r

## Ensure OpenMP threads are "close"
export OMP_PROC_BIND="CLOSE"

Rscript --no-save "simstudy-scenario_3.R" \
  --ncores ${SLURM_CPUS_PER_TASK} \
  --job ${SLURM_ARRAY_TASK_ID}
