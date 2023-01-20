##
## Compute all estimates on the full glass vessel data set.
##
library(argparser)
library(magrittr)

args <- arg_parser('Compute all estimates on the full glass vessel data set.') %>%
  add_argument('--results-dir', type = 'character', help = 'Result path',
               default = file.path('results', 'real_da')) %>%
  add_argument('--ncores', default = 4L, help = 'Number of CPUs') %>%
  parse_args()

## Cache path settings
CACHE_PATH <- file.path(args$results_dir, 'cache')

## General settings
CV_K <- 6     # 6-fold cross-validation
CV_REPL <- 20 # 20 replications of cross-validation
ALPHA_SEQUENCE <- c(0.5, 0.75, 1)  # sequence of alpha parameters for EN-type estimators
ZETA_SEQUENCE <- c(1, 2)  # sequence of zeta parameters for adaptive estimators
BASE_SEED <- 12345  # base seed for all computations
NUMERIC_EPS <- 1e-6  # numerical tolerance level
PENALTY_LEVELS <- 50  # number of penalty levels to consider

## Settings for PENSE and adaptive PENSE
PENSE_INITIAL_PENALTY_LEVELS <- 10
PENSE_RETAIN_INITIAL_CANDIDATES <- 25
PENSE_BDP <- 0.28  # ~50 obs

## Determine the parallelization (threading or multiple processes)
if (args$ncores > 1L) {
  args$total_ncores <- args$ncores
  if (pense:::.k_multithreading_support) {
    args$cluster <- NULL
  } else {
    args$cluster <- parallel::makeCluster(args$ncores)
    args$ncores <- 1L
  }
} else {
  args$total_ncores <- 1L
  args$ncores <- 1L
  args$cluster <- NULL
}

## Load the utilities
source('utilities.R')

## Load the data set from other R packages
requireNamespace('chemometrics')
requireNamespace('cellWise')
glass_data_env <- new.env(parent = emptyenv())
data(glass, package = 'chemometrics', envir = glass_data_env)
data(data_glass, package = 'cellWise', envir = glass_data_env)

glass_y_orig <- glass_data_env$glass[, 'P2O5']
glass_y <- log(glass_y_orig)
glass_x <- as.matrix(glass_data_env$data_glass[, 15:500])

# Create results and cache path if necessary
if (!dir.exists(args$results_dir)) {
  dir.create(args$results_dir, recursive = TRUE, mode = '0700')
}
if (!dir.exists(CACHE_PATH)) {
  dir.create(CACHE_PATH, recursive = TRUE, mode = '0700')
}

estimates <- list()

## Compute PENSE estimates
estimates$pense <- compute_adapense_cv(
  y = glass_y, x = glass_x,
  nlambda = PENALTY_LEVELS,
  nlambda_enpy = PENSE_INITIAL_PENALTY_LEVELS,
  lambda_min_ratio = 2e-2,
  seed = BASE_SEED,
  alpha = ALPHA_SEQUENCE,
  cv_repl = CV_REPL,
  cv_k = CV_K,
  bdp = PENSE_BDP,
  fit_all = c('min', 'se'),
  retain_initial = PENSE_RETAIN_INITIAL_CANDIDATES,
  ncores = args$ncores,
  cl = args$cluster,
  eps = NUMERIC_EPS,
  en_algo_opts = en_lars_options(),
  cache_path = CACHE_PATH,
  penalty_loadings = NULL)

## Compute adaptive PENSE estimates
estimates$adapense <- compute_adapense_cv_zeta(
  y = glass_y, x = glass_x,
  nlambda = PENALTY_LEVELS,
  nlambda_enpy = PENSE_INITIAL_PENALTY_LEVELS,
  lambda_min_ratio = c(2e-2, 2e-3),
  lambda_min_ratio_prelim = 1e-1,
  seed = BASE_SEED,
  alpha = ALPHA_SEQUENCE,
  zeta_seq = ZETA_SEQUENCE,
  cv_repl = CV_REPL,
  cv_k = CV_K,
  bdp = PENSE_BDP,
  fit_all = c('min', 'se'),
  retain_initial = PENSE_RETAIN_INITIAL_CANDIDATES,
  ncores = args$ncores,
  cl = args$cluster,
  eps = NUMERIC_EPS,
  en_algo_opts = en_lars_options(),
  cache_path = CACHE_PATH)

## Compute MM estimates from the scale estimated via PENSE-Ridge
estimates$mm <- compute_adammest_cv(
  y = glass_y, x = glass_x,
  pense_ridge = estimates$adapense$preliminary$cv_min,
  nlambda = PENALTY_LEVELS,
  lambda_min_ratio = 1e-1,
  seed = BASE_SEED,
  alpha = ALPHA_SEQUENCE,
  cv_repl = CV_REPL,
  cv_k = CV_K,
  mscale_bdp = PENSE_BDP,
  fit_all = c('min', 'se'),
  ncores = args$ncores,
  cl = args$cluster,
  eps = NUMERIC_EPS,
  en_algo_opts = en_lars_options(),
  cache_path = CACHE_PATH)

## Compute adaptive MM estimates
estimates$adamm <- compute_adammest_cv_zeta(
  y = glass_y, x = glass_x,
  pense_ridge = estimates$adapense$preliminary$cv_min,
  nlambda = PENALTY_LEVELS,
  lambda_min_ratio = c(5e-2, 2e-3),
  seed = BASE_SEED,
  alpha = ALPHA_SEQUENCE,
  zeta_seq = ZETA_SEQUENCE,
  cv_repl = CV_REPL,
  cv_k = CV_K,
  mscale_bdp = PENSE_BDP,
  fit_all = c('min', 'se'),
  ncores = args$ncores,
  cl = args$cluster,
  eps = NUMERIC_EPS,
  en_algo_opts = en_lars_options(),
  cache_path = CACHE_PATH)

## Compute ILAMM estimates
estimates$ilamm <- compute_ilamm(
  y = glass_y, x = glass_x,
  nlambda = PENALTY_LEVELS,
  seed = BASE_SEED,
  cv_k = CV_K,
  cv_repl = CV_REPL,
  ncores = args$total_ncores,
  cache_path = CACHE_PATH)

## Compute ILAMM (SCAD) estimates
estimates$ilamm_scad <- compute_ilamm(
  y = glass_y, x = glass_x,
  nlambda = PENALTY_LEVELS,
  seed = BASE_SEED,
  cv_k = CV_K,
  cv_repl = CV_REPL,
  penalty = 'SCAD',
  ncores = args$total_ncores,
  cache_path = CACHE_PATH)

## Compute LS-EN estimates
estimates$en <- comp_glmnet(
  y = glass_y, x = glass_x,
  nlambda = PENALTY_LEVELS,
  alpha = ALPHA_SEQUENCE,
  seed = BASE_SEED,
  cv_k = CV_K,
  cv_repl = CV_REPL,
  cache_path = CACHE_PATH)

## Compute adaptive LS-EN estimates
estimates$adaen <- comp_glmnet_zeta(
  y = glass_y, x = glass_x,
  nlambda = PENALTY_LEVELS,
  alpha = ALPHA_SEQUENCE,
  zeta_seq = ZETA_SEQUENCE,
  seed = BASE_SEED,
  cv_k = CV_K,
  cv_repl = CV_REPL,
  cache_path = CACHE_PATH)

if (!is.null(args$cluster)) {
  parallel::stopCluster(args$cluster)
}

## Save results for later use
saveRDS(estimates, file = file.path(args$results_dir, 'full_estimates.rds'))
