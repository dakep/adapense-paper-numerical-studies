##
## Estimate prediction accuracy of adaptive PENSE with varying breakdown point
## for the glass vessel data set via nested CV.
##
suppressPackageStartupMessages({
  library(rlang)
  library(argparser)
  library(magrittr)
})

args <- arg_parser(paste('Estimate prediction accuracy of adaptive PENSE with',
                         'varying breakdown point for the glass vessel data set',
                         'via nested CV.')) %>%
  add_argument('--results-dir', type = 'character', help = 'Result path',
               default = file.path('results', 'real_da')) %>%
  add_argument('--ncores', default = 4L, help = 'Number of CPUs to utilize.') %>%
  add_argument('--job', type = 'integer', help = 'Job number(s).', nargs = Inf) %>%
  parse_args()

## Cache path settings
CACHE_PATH <- file.path(args$results_dir, 'cache')
JOB_RESULTS_PATH <- file.path(args$results_dir, 'cv')

## General settings
CV_K <- 6     # 6-fold cross-validation
CV_REPL <- 20 # 20 replications of cross-validation
ALPHA_SEQUENCE <- c(0.5, 0.75, 1)  # sequence of alpha parameters for EN-type estimators
ZETA_SEQUENCE <- c(1, 2)  # sequence of zeta parameters for adaptive estimators
NUMERIC_EPS <- 1e-6  # numerical tolerance level
PENALTY_LEVELS <- 50  # number of penalty levels to consider

## Settings for PENSE and adaptive PENSE
PENSE_INITIAL_PENALTY_LEVELS <- 10
PENSE_RETAIN_INITIAL_CANDIDATES <- 25
PENSE_BDP <- c(0.2, 0.28, 0.33, 0.4)  # 30-60 potential outliers

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

## Load the utilities
source('utilities.R')

## Load the data set from other R packages
requireNamespace('chemometrics', quietly = TRUE)
requireNamespace('cellWise', quietly = TRUE)
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
if (!dir.exists(JOB_RESULTS_PATH)) {
  dir.create(JOB_RESULTS_PATH, recursive = TRUE, mode = '0700')
}
if (!dir.exists(CACHE_PATH)) {
  dir.create(CACHE_PATH, recursive = TRUE, mode = '0700')
}

if (anyNA(args$job)) {
  args$job <- 1L
}

for (job in unique(args$job)) {
  ## Determine seed and CV fold
  job_seed <- 1L + (job - 1L) %/% CV_K
  job_fold <- 1L + (job - 1L) %% CV_K

  ## Split data into training and test set
  set.seed(job_seed)
  cv_folds <- split(seq_along(glass_y), sample(rep_len(seq_len(CV_K), length(glass_y))))

  test_glass_y <- glass_y[cv_folds[[job_fold]]]
  test_glass_x <- glass_x[cv_folds[[job_fold]], ]

  train_glass_y <- glass_y[-cv_folds[[job_fold]]]
  train_glass_x <- glass_x[-cv_folds[[job_fold]], ]

  for (bdp in PENSE_BDP) {
    job_cache_path <- file.path(CACHE_PATH, 'cv',
                                sprintf('%04d-%02d-%04d', job_seed, job_fold, 1000 * bdp))
    if (!dir.exists(job_cache_path)) {
      dir.create(job_cache_path, recursive = TRUE, mode = '0700')
    }

    print_log('Computing prediction errors for bdp %f in job %d (seed: %03d, fold: %d)',
              bdp, job, job_seed, job_fold)
    ## Set up the result structure
    cv_results <- list(
      job = list(
        id = job,
        bdp = bdp,
        seed = job_seed,
        fold = job_fold,
        test_fold = cv_folds[[job_fold]]
      ),
      estimates = list()
    )
    ## Compute adaptive PENSE estimates
    cv_results$estimates$adapense <- compute_adapense_cv_zeta(
      y = train_glass_y, x = train_glass_x,
      nlambda = PENALTY_LEVELS,
      nlambda_enpy = PENSE_INITIAL_PENALTY_LEVELS,
      lambda_min_ratio = c(5e-2, 4e-4),
      lambda_min_ratio_prelim = 1e-2,
      seed = job_seed,
      alpha = ALPHA_SEQUENCE,
      zeta_seq = ZETA_SEQUENCE,
      cv_repl = CV_REPL,
      cv_k = CV_K,
      bdp = bdp,
      fit_all = c('min', 'se'),
      retain_initial = PENSE_RETAIN_INITIAL_CANDIDATES,
      ncores = args$ncores,
      cl = args$cluster,
      eps = NUMERIC_EPS,
      cache_path = job_cache_path,
      en_algo_opts = en_lars_options(),
      log_indent = 1)

    ## Evaluate all estimates on the test set
    cv_results$estimates <- evaluate_estimate(cv_results$estimates,
                                              test_x = test_glass_x,
                                              test_y = test_glass_y)
    ## Save results
    saveRDS(cv_results, file = file.path(JOB_RESULTS_PATH,
                                         sprintf('%04d-%02d-%04d.rds',
                                                 job_seed, job_fold, bdp)))
  }
}
