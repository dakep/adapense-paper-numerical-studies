##
## Simulation study -- Scenario 3
##
suppressPackageStartupMessages({
  library(rlang)
  library(argparser)
  library(magrittr)
  library(stringr)
  library(digest)
  library(pense)
})

args <- arg_parser('Run the simulation study for scenario 3.') %>%
  add_argument('--results-dir', type = 'character', help = 'Result path',
               default = file.path('results', 'simstudy', 'scenario_03')) %>%
  add_argument('--ncores', default = 4L, help = 'Number of CPUs to utilize.') %>%
  add_argument('--job', type = 'integer', help = 'Job number, i.e., the seed for the RNG.') %>%
  parse_args()

## Source other scripts
source('simstudy-generate_data.R')
source('utilities.R')

## Simulation settings to consider
SIM_P <- 2^6 # number of predictors
SIM_N <- c(10, 15, 20, 25)^2 # number of observations
SIM_RESID_DIST <- 'stable(alpha = 1.33, beta = 0)'  # error distributions
SIM_CONT_PROP <- 0 # contamination proportion

## Cache path settings
CACHE_PATH <- file.path(args$results_dir, 'cache')

## Estimator settings
CV_K <- 5L # use 5-fold cross-validation
CV_REPL <- 20L # replicate CV 10 times
ALPHA_SEQUENCE <- 0.75 # alpha parameters for EN-type estimators
ZETA_SEQUENCE <- 1  # sequence of zeta parameters for adaptive estimators
NUMERIC_EPS <- 1e-6
PENALTY_LEVELS <- 50  # number of penalty levels to consider
PENSE_INITIAL_PENALTY_LEVELS <- 10 # number of penalty levels where initial estimates are computed
PENSE_RETAIN_INITIAL_CANDIDATES <- 25  # use only the best 25 initial estimates
PENSE_BDP <- 1/3  # desired breakdown point

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


#' Generate Data for Scenario 3
#'
#' @param p number of predictors.
#' @param resid_dist distribution of the residuals.
#' @param contamination_percentage percentage of observations affected by outliers in the
#'   response and bad leverage points.
#' @param test_n,test_seed number of test observations and the seed for the RNG to
#'   generate these test observations.
generate_data_scenario_3 <- function (n, test_n = 1000, test_seed = 7357) {
  settings <- list(
    n = n,
    p = SIM_P,
    s = floor(0.9 * sqrt(SIM_P)),
    cont = SIM_CONT_PROP,
    k_vert = -1,
    bad_lev_p = log2(SIM_P),
    bad_lev_p_relevant = 0,
    bad_lev_pos = 2,
    good_lev_p = NA_integer_,
    good_lev_prop = 0.2,
    good_lev_pos = 8,
    pve = 0.25,
    resid_dist = SIM_RESID_DIST,
    pred_cor_type = 'grouped'
  )

  settings$good_lev_p <- with(settings, ceiling(0.5 * min(n, (p - s))))
  generate_data(settings, test_n, test_seed)
}

## Ensure that paths exist or create them if necessary
if (!dir.exists(args$results_dir)) {
  dir.create(args$results_dir, recursive = TRUE, mode = '0700')
}

## Use the LARS algorithm for all settings
en_algo_opts <- en_lars_options()

## Run all combinations for the given seed
for (n in SIM_N) {
  job_id <- sprintf('s=%05d-n=%05d', args$job, n)
  job_hash <- digest(job_id, serialize = TRUE)
  job_fname <- paste(str_replace_all(job_id, '[^a-zA-Z0-9-_]+', '_'),
                     str_sub(job_hash, end = 7),
                     sep = '-')
  job_file <- file.path(args$results_dir, paste(job_fname, 'rds', sep = '.'))

  job_cache_path <- file.path(CACHE_PATH, job_fname)
  if (!dir.exists(job_cache_path)) {
    dir.create(job_cache_path, recursive = TRUE, mode = '0700')
  }

  print_log("Computing estimates for job with id %s.", job_id)
  # Generate data
  set.seed(args$job)
  simdat <- generate_data_scenario_3(n)

  # Build final result object
  cv_results <- list(
    seed = args$job,
    settings = simdat$settings,
    oracle = c(sd = sd(simdat$residuals),
               tau_size = tau_size(simdat$residuals),
               mscale = mscale(simdat$residuals, bdp = PENSE_BDP)),
    true_beta = simdat$beta,
    true_sd = simdat$sd_err,
    cont_pred_good = simdat$cont_pred_good,
    cont_pred_bad = simdat$cont_pred_bad,
    estimates = list()
  )

  ## Compute the intercept-only model
  cv_results$estimates$intonly <- compute_intonly(simdat$x, simdat$y, PENSE_BDP)

  ## Compute PENSE estimates
  cv_results$estimates$pense <- compute_adapense_cv(
    y = simdat$y, x = simdat$x,
    nlambda = PENALTY_LEVELS,
    nlambda_enpy = PENSE_INITIAL_PENALTY_LEVELS,
    seed = args$job,
    alpha = ALPHA_SEQUENCE,
    cv_repl = CV_REPL,
    cv_k = CV_K,
    bdp = PENSE_BDP,
    fit_all = c('min', 'se'),
    retain_initial = PENSE_RETAIN_INITIAL_CANDIDATES,
    ncores = args$ncores,
    cl = args$cluster,
    eps = NUMERIC_EPS,
    en_algo_opts = en_algo_opts,
    cache_path = job_cache_path,
    penalty_loadings = NULL,
    log_indent = 1)

  ## Compute adaptive PENSE estimates
  cv_results$estimates$adapense <- compute_adapense_cv_zeta(
    y = simdat$y, x = simdat$x,
    nlambda = PENALTY_LEVELS,
    nlambda_enpy = PENSE_INITIAL_PENALTY_LEVELS,
    seed = args$job,
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
    en_algo_opts = en_algo_opts,
    cache_path = job_cache_path,
    log_indent = 1)

  ## Compute MM estimates from the scale estimated via PENSE-Ridge
  cv_results$estimates$mm <- compute_adammest_cv(
    y = simdat$y, x = simdat$x,
    pense_ridge = cv_results$estimates$adapense$preliminary$cv_min,
    nlambda = PENALTY_LEVELS,
    seed = args$job,
    alpha = ALPHA_SEQUENCE,
    cv_repl = CV_REPL,
    cv_k = CV_K,
    fit_all = c('min', 'se'),
    mscale_bdp = PENSE_BDP,
    ncores = args$ncores,
    cl = args$cluster,
    eps = NUMERIC_EPS,
    en_algo_opts = en_algo_opts,
    cache_path = job_cache_path,
    log_indent = 1)

  ## Compute adaptive MM estimates
  cv_results$estimates$adamm <- compute_adammest_cv_zeta(
    y = simdat$y, x = simdat$x,
    pense_ridge = cv_results$estimates$adapense$preliminary$cv_min,
    nlambda = PENALTY_LEVELS,
    seed = args$job,
    alpha = ALPHA_SEQUENCE,
    zeta_seq = ZETA_SEQUENCE,
    cv_repl = CV_REPL,
    cv_k = CV_K,
    fit_all = c('min', 'se'),
    mscale_bdp = PENSE_BDP,
    ncores = args$ncores,
    cl = args$cluster,
    eps = NUMERIC_EPS,
    en_algo_opts = en_algo_opts,
    cache_path = job_cache_path,
    log_indent = 1)

  ## Compute ILAMM estimates
  cv_results$estimates$ilamm <- compute_ilamm(
    y = simdat$y, x = simdat$x,
    nlambda = PENALTY_LEVELS,
    seed = args$job,
    cv_k = CV_K,
    cv_repl = CV_REPL,
    ncores = args$total_ncores,
    cache_path = job_cache_path,
    log_indent = 1)

  ## Compute ILAMM (SCAD) estimates
  cv_results$estimates$ilamm_scad <- compute_ilamm(
    y = simdat$y, x = simdat$x,
    nlambda = PENALTY_LEVELS,
    seed = args$job,
    cv_k = CV_K,
    cv_repl = CV_REPL,
    cache_path = job_cache_path,
    penalty = 'SCAD',
    ncores = args$total_ncores,
    log_indent = 1)

  ## Compute LS-EN estimates
  cv_results$estimates$en <- comp_glmnet(
    y = simdat$y, x = simdat$x,
    nlambda = PENALTY_LEVELS,
    alpha = ALPHA_SEQUENCE,
    seed = args$job,
    cv_k = CV_K,
    cv_repl = CV_REPL,
    cache_path = job_cache_path,
    log_indent = 1)

  ## Compute adaptive LS-EN estimates
  cv_results$estimates$adaen <- comp_glmnet_zeta(
    y = simdat$y, x = simdat$x,
    nlambda = PENALTY_LEVELS,
    alpha = ALPHA_SEQUENCE,
    zeta_seq = ZETA_SEQUENCE,
    seed = args$job,
    cv_k = CV_K,
    cv_repl = CV_REPL,
    cache_path = job_cache_path,
    log_indent = 1)

  cv_results$estimates <- evaluate_estimate(cv_results$estimates,
                                            test_x = simdat$test_data$x,
                                            test_y = simdat$test_data$y)

  print_log("Saving results for job %s to '%s'.", job_id, job_file)
  saveRDS(cv_results, file = job_file)
}
