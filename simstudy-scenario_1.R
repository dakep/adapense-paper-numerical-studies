##
## Simulation study -- Scenario 1
##
suppressPackageStartupMessages({
  library(rlang)
  library(argparser)
  library(magrittr)
  library(stringr)
  library(digest)
  library(pense)
})

args <- arg_parser('Run the simulation study for scenario 1.') %>%
  add_argument('--results-dir', type = 'character', help = 'Result path',
               default = file.path('results', 'simstudy', 'scenario_01')) %>%
  add_argument('--ncores', default = 4L, help = 'Number of CPUs to utilize.') %>%
  add_argument('--job', type = 'integer', help = 'Job number, i.e., the seed for the RNG.') %>%
  parse_args()

## Source other scripts
source('simstudy-generate_data.R')
source('utilities.R')

## Simulation settings to consider
SIM_P <- 2^c(5, 7, 9) # number of predictors
SIM_RESID_DIST <- c('norm', 'cauchy', 'stable(alpha = 1.33, beta = 0)')  # error distributions
SIM_CONT_PROP <- c(0, 0.1) # contamination proportions

## Cache path settings
CACHE_PATH <- file.path(args$results_dir, 'cache')

## Estimator settings
CV_K <- 5L # use 5-fold cross-validation
CV_REPL <- 10L # replicate CV 10 times
ALPHA_SEQUENCE <- c(0.5, 0.75, 1) # alpha parameters for EN-type estimators
ZETA_SEQUENCE <- c(1, 2)  # sequence of zeta parameters for adaptive estimators
NUMERIC_EPS <- 1e-6
PENALTY_LEVELS <- 50  # number of penalty levels to consider
PENSE_INITIAL_PENALTY_LEVELS <- 10 # number of penalty levels where initial estimates are computed
PENSE_RETAIN_INITIAL_CANDIDATES <- 25  # use only the best 25 initial estimates
PENSE_BDP <- 1/3  # desired breakdown point

## Determine the parallelization (threading or multiple processes)
if (args$ncores > 1L) {
  if (pense:::.k_multithreading_support) {
    args$cluster <- NULL
  } else {
    args$cluster <- parallel::makeCluster(args$ncores)
    args$ncores <- 1L
  }
} else {
  args$ncores <- 1L
  args$cluster <- NULL
}

#' Generate Data for Scenario 1
#'
#' @param p number of predictors.
#' @param resid_dist distribution of the residuals.
#' @param contamination_percentage percentage of observations affected by outliers in the
#'   response and bad leverage points.
#' @param test_n,test_seed number of test observations and the seed for the RNG to
#'   generate these test observations.
generate_data_scenario_1 <- function (p, resid_dist, contamination_percentage,
                                      test_n = 1000, test_seed = 7357) {
  settings <- list(
    n = 200,
    p = p,
    s = floor(0.9 * sqrt(p)),
    cont = contamination_percentage,
    k_vert = -1,
    bad_lev_p = log2(p),
    bad_lev_p_relevant = 0,
    bad_lev_pos = 2,
    good_lev_p = NA_integer_,
    good_lev_prop = 0.2,
    good_lev_pos = 8,
    pve = 0.25,
    resid_dist = resid_dist,
    pred_cor_type = 'grouped'
  )

  settings$good_lev_p <- with(settings, ceiling(0.5 * min(n, (p - s))))
  generate_data(settings, test_n, test_seed)
}

## Ensure that paths exist or create them if necessary
if (!dir.exists(args$results_dir)) {
  dir.create(args$results_dir, recursive = TRUE, mode = '0700')
}

## Run all combinations for the given seed
for (p in SIM_P) {
  for (resid_dist in SIM_RESID_DIST) {
    for (cont in SIM_CONT_PROP) {
      job_id <- sprintf('s=%05d-p=%05d-rd=%s-c=%.2f', args$job, p, resid_dist, cont)
      job_hash <- digest(job_id, serialize = TRUE)
      job_fname <- paste(str_replace_all(job_id, '[^a-zA-Z0-9-_]+', '_'),
                         str_sub(job_hash, end = 7),
                         sep = '-')
      job_file <- file.path(args$results_dir, paste(job_fname, 'rds', sep = '.'))

      job_cache_path <- file.path(CACHE_PATH, job_fname)
      if (!dir.exists(job_cache_path)) {
        dir.create(job_cache_path, recursive = TRUE, mode = '0700')
      }

      en_algo_opts <- if (isTRUE(p > 400)) {
        en_dal_options()
      } else {
        en_lars_options()
      }

      print_log("Computing estimates for job with id %s.", job_id)
      # Generate data
      set.seed(args$job)
      simdat <- generate_data_scenario_1(p, resid_dist, cont)

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
        cache_path = job_cache_path,
        log_indent = 1)

      ## Compute ILAMM estimates
      cv_results$estimates$ilamm <- compute_ilamm(
        y = simdat$y, x = simdat$x,
        nlambda = PENALTY_LEVELS,
        seed = args$job,
        cv_k = CV_K,
        cv_repl = CV_REPL,
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
  }
}
