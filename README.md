# Numerical Studies for Paper "Robust Variable Selection and Estimation Via Adaptive Elastic Net S-Estimators for Linear Regression"

[![DOI](https://zenodo.org/badge/387899778.svg)](https://zenodo.org/badge/latestdoi/387899778)

This repository contains all necessary code and data to reproduce the results in "Robust Variable Selection and Estimation Via Adaptive Elastic Net S-Estimators for Linear Regression".

Re-computing all results takes about 33,000 CPU hours when using an R build linked to an optimized BLAS/LAPACK library (e.g., Intel® MKL).
Since it is not feasible to reproduce all results without access to a high performance cluster (HPC), the repository also contains R data files with results for individual replications (called "jobs").
This allows re-computing only some of these jobs and still generate all results in the manuscript.

## Dependencies

The results have been generated with R version 4.1.0, but the R code is compatible with R ≥ 4.0.
This repository contains all relevant code to set-up all other dependencies automatically via the [renv](https://rstudio.github.io/renv/) R package.
To install all necessary dependencies, start R within this repository and run the following line of code:

```r
renv::restore()
```

## How to run the code

The project is organized as several stand-alone R scripts and R markdown files to compute the estimates and create figures/tables, respectively.
All scripts are designed to be called from the command line with *Rscript* and require one or more command-line arguments.
For demonstration purposes, this repository also includes SLURM scripts which can be used as reference for reproducing the results on a HPC which uses the SLURM workload manager.

The code is divided into two groups: scripts for the [real data analysis](#real-data-analysis) are prefixed by `real_data_example`, and scripts for the [simulation study](#simulation-study) start with `simstudy`.

### Real data analysis

The real study consists of the following R scripts and R markdown files:

* *real_data_analysis-full_data.R* … computes estimates on the full glass vessel data set ([details](#computing-esimates-on-the-full-data-set)).
* *real_data_analysis-prediction_error.R* … estimates the prediction performance via nested CV ([details](#estimating-the-prediction-accuracy)).
* *real_data_analysis.Rmd* … generates a report of the results, including the figures and tables for the manuscript.

In addition, the folder *results/real_da* contains the R data files as generated by the R scripts *real_data_analysis-full_data.R* and *real_data_analysis-prediction_error.R*.

#### Computing esimates on the full data set

The R script *real_data_analysis-full_data.R* computes estimates using all considered methods on the full glass vessel data set.
It uses 10 replications of 6-fold CV to select the hyper-parameters for each method.
To run the R script using 4 CPUs, for example, the script is called from the command line as follows:

```sh
Rscript real_data_analysis-full_data.R --ncores 4
```

The script takes ~6 CPU hours to compute and accepts the arguments `--ncores` (the number of CPUs to use in parallel) and `--results-dir` (the path where result file *full_estimates.rda* will be stored).
By default, the result file is saved under *results/real_da*, and the folder is already populated with the result file computed for the manuscript: *results/real_da/full_estimates.rda*.

#### Estimating the prediction accuracy

Prediction errors for **a single CV fold** are estimated via the R script *real_data_analysis-prediction_error.R*.
In addition to the arguments `--ncores` and `--results-dir`, this script requires the user to specify the `--job` argument.
This `--job` argument takes an integer value which defines which CV fold is to be performed.
One job corresponds to a specific seed for the random number generator and the CV fold.
The seed is determined as `1 + floor((job - 1) / 6)` and the fold as `1 + (job - 1) %% 6`.

For example, to compute the estimates for job 47 (i.e., `seed = 8`, `fold = 5`) using 4 CPUs call the R script from the command line as follows:

```sh
Rscript real_data_analysis-prediction_error.R --ncores 4 --job 47 --
```

(Please note the two trailing dashes (`--`), which are required to delineate the end of the job specifier.)

Assuming the `--result-dir` argument is *results/real_da*, the script saves the result file under *results/real_da/cv/{seed}-{fold}.rds*.
This repository contains all 300 result files computed for the manuscript (*0001-01.rds* to *0050-06.rds*).

For the manuscript, a total of 50 replications of 6 fold CV were performed (jobs 1 to 300).
Each job takes ~3 CPU hours.
To replicate all results in the manuscript, the R script must be called for each job:

```sh
#!/bin/bash

for ((j = 1; j <= 300; j++)); do
  Rscript real_data_analysis-prediction_error.R --ncores 4 --job $j --
done
```

On HPC systems using the SLURM workload manager, the results can be replicated with

```sh
sbatch \
  --array=1-300 \
  --time=2:30:00 \
  --output="logs/real_da/prediction_error-%a.out" \
  --error="logs/real_da/prediction_error-%a.err"  \
  --open-mode=append \
  real_data_analysis-slurm_cv.sh
```

The script *real_data_analysis-prediction_error-varbdp.R* estimates the prediction error of adaptive PENSE for varying breakdown points.
It can be run like the script *real_data_analysis-prediction_error.R*, except that it is separated into only 50 tasks, one per seed.
A single tasks requires at most 35 CPU hours.


#### Re-creating the figures and tables

The file `real_data_analysis.Rmd` re-creates all plots and other results reported in the manuscript.
To compile the R markdown file to an HTML file, start a new R session and run

```r
rmarkdown::render("real_data_analysis.Rmd")
```

or, when using the RStudio IDE, press the "Knit" button.

This R markdown file reads the result files from the path *results/real_da*.

### Simulation study

The simulation study consists of two main scenarios and the results concerning good leverage points.
Results for the three simulation scenarios are computed by separate R scripts: *simstudy-scenario_1.R*, *simstudy-scenario_2.R* and *simstudy-scenario_3.R*.
They are almost identical, except for the differences in how the data is generated.
The instructions below only reference the script for scenario 1, and results for scenarios 2 and 3 can be re-created by the exact same steps (substituting all occurrences of `scenario_1` with `scenario_2` or `scenario_3`).

The script *simstudy-scenario_1.R* computes all estimates for all settings for a single seed (i.e., a single replication).
To compute the estimates for the 15th replication using 4 CPUs, for example, the script is called on the command line as follows:

```sh
Rscript simstudy-scenario_1.R --ncores 4 --job 15 --
```

One job requires approximately 480 CPU hours for scenario 1.
Scenario 2 takes less than 60 CPU hours, scenario 3 less than 20 CPU hours.
By default, the script uses the result path *results/simstudy*, and in this case the result files are stored under *results/simstudy/scenario_01/{job_setting_id}.rds*, where *{job_setting_id}* is a unique identifier for the job and setting.

On HPC systems using the SLURM workload manager, the results for can be replicated with

```sh
sbatch \
  --array=1-50 \
  --time=24:00:00 \
  --output="logs/simstudy/scenario_1-%a.out" \
  --error="logs/simstudy/scenario_1-%a.err" \
  --open-mode=append \
  simstudy-scenario_1-slurm.sh
```

#### Good leverage points

The results pertaining to the behavior of estimators under good leverage points are re-created with the script *simstudy-good_leverage.R*.
It can be run just like the other two scripts for the simulation study with jobs 1 -- 50.
For example, to re-compute results for the 22nd replication, the script would be run from the command line as follows:

```sh
Rscript simstudy-good_leverage.R --ncores 4 --job 22 --
```

This job takes up to 36 CPU hours.
Result files are saved under *results/simstudy/good_leverage*.

#### Estimation accuracy

The results pertaining to increasing estimation accuracy with increasing $n$ can be reproduced with the script *simstudy-scenario_3.R*.
It can be run just like the other scripts for the simulation study with jobs 1 -- 50.
For example, to re-compute results for the 22nd replication, the script would be run from the command line as follows:

```sh
Rscript simstudy-scenario_3.R --ncores 4 --job 22 --
```

This job takes up to 55 CPU hours.
Result files are saved under *results/simstudy/scenario_03*.

#### Re-creating the figures and tables

The R markdown file *simstudy.Rmd* reads in the result files from *results/simstudy/scenario_01*, *results/simstudy/scenario_02*, and *results/simstudy/good_leverage* and re-creates the figures and other information reported in the manuscript.

To compile the R markdown file to an HTML file, start a new R session and run

```r
rmarkdown::render("real_data_analysis.Rmd")
```

or, when using the RStudio IDE, press the "Knit" button.

## License

The code and derived files in this repository are © David Kepplinger (2021) and published under the [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)](https://creativecommons.org/licenses/by-sa/4.0/) license.

The file *renv/active.R* is © RStudio, PBC (2021) and distributed under the MIT License

See file [*LICENSE*](./LICENSE) for details.
