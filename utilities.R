#' Compute PENSE and adaptive PENSE estimates for a sequence of *alpha* values **and** cache
#' individual results.
#'
#' _NOTE_: this is also possible with the [pense_cv()] package directly, but here we add
#' caching to avoid having to re-compute all results!
#'
#' The same seed is used for all cross-validations (i.e., for the different *alpha* values).
#' Only the estimate with the best prediction performance (among all *alpha* values) is returned.
#'
#' @param ... passed on to `pense::pense_cv()`.
#' @param penalty_loadings passed on to `pense::pense_cv()`
#' @param cache_path path for caching intermediate results.
#' @param alpha a set of different *alpha* values.
#' @param seed seed for the individual cross-validations (the same seed is used for all
#' *alpha* values).
#' @param retain_initial passed on to `pense::enpy_options()`.
#' @param log_indent indentation of log messages.
#' @param en_algo_opts specify the EN algorithm options.
#' @param cv_k,cv_repl,sparse,max_solutions,enpy_opts,algorithm_opts further arguments to
#'    `pense::pense_cv()`.
compute_adapense_cv <- function (..., cv_k, seed, alpha, retain_initial, cache_path,
                                 cv_repl, penalty_loadings = NULL, sparse = FALSE,
                                 max_solutions = 1, en_algo_opts = NULL, enpy_opts, algorithm_opts,
                                 zeta = NA_real_, log_indent = 0L) {
  if (!require(pense, quietly = TRUE)) {
    stop('`pense` package not available')
  }
  if (!require('Matrix', quietly = TRUE)) {
    stop('`Matrix` package not available')
  }

  est_name <- if (is.null(penalty_loadings)) { 'PENSE' } else { 'adaptive PENSE' }

  if (is.null(en_algo_opts)) {
    en_algo_opts <- en_cd_options()
  }
  if (missing(algorithm_opts)) {
    algorithm_opts <- mm_algorithm_options(en_algorithm_opts = en_algo_opts)
  }
  if (missing(enpy_opts)) {
    enpy_opts <- enpy_options(en_algorithm_opts = en_algo_opts, retain_max = retain_initial)
  }

  # Run for every `alpha`
  all_alpha_results <- lapply(alpha, function (alpha) {
    print_log("Computing %s estimate for alpha=%0.2f.", est_name, alpha,
              .indent = log_indent)

    try_catch({
      cache <- load_cache('%010d,%010d,%010d,%.2f|%s', seed, cv_k, cv_repl, alpha,
                          hash_num_vec(penalty_loadings),
                          .prefix = 'adapense_', .cache_path = cache_path,
                          .indent = log_indent + 1L)

      if (!is.null(cache$object)) {
        cache$object
      } else {
        # Set the seed so that we get the same CV splits for every alpha
        set.seed(seed)
        alpha_start <- proc.time()[["elapsed"]]
        estimates <- pense_cv(..., alpha = alpha, penalty_loadings = penalty_loadings,
                              sparse = sparse, max_solutions = max_solutions, cv_k = cv_k,
                              cv_repl = cv_repl, enpy_opts = enpy_opts,
                              algorithm_opts = algorithm_opts)
        alpha_end <- proc.time()[["elapsed"]]

        pense_alpha_res <- extract_pense_cv_coef(estimates,
                                                 lambda = if (isTRUE(cv_repl > 1)) c('min', 'se')
                                                 else 'min',
                                                 zeta = zeta)
        pense_alpha_res$duration <- alpha_end - alpha_start
        print_log("Saving results to cache file '%s'.", cache$file, .indent = log_indent + 1L)
        saveRDS(pense_alpha_res, file = cache$file)

        pense_alpha_res
      }
    })
  })
  duration_all_alpha <- sum(vapply(all_alpha_results, FUN.VALUE = numeric(1),
                                   FUN = `[[`, 'duration'))

  # Select the best alpha value
  try_catch({
    adapense_ests <- list(cv_min = get_best_estimate(all_alpha_results, 'min'),
                          cv_se = get_best_estimate(all_alpha_results, 'se'))

    if (isFALSE(sparse)) {
      # Make sparse coefficient vectors (if necessary -- sparse is only a suggestion)
      adapense_ests <- lapply(adapense_ests, function (est) {
        if (!is.null(est) && is.numeric(est$beta)) {
          nnz_ind <- which(abs(est$beta) > 1e-12)
          est$beta <- sparseVector(est$beta[nnz_ind], i = nnz_ind, length = length(est$beta))
        }
        return(est)
      })
    }

    adapense_ests$cv_min$duration_all_alpha <- duration_all_alpha
    adapense_ests$cv_se$duration_all_alpha <- duration_all_alpha

    adapense_ests
  })
}

#' Compute regularized MM and adaptive MM estimates for a sequence of *alpha* values and cache
#' intermediate results.
#'
#' The same seed is used for all cross-validations (i.e., for the different *alpha* values).
#' Only the estimate with the best prediction performance (among all *alpha* values) is returned.
#'
#' @param ... passed to `pense::regmest_cv()`.
#' @param pense_ridge the PENSE ridge estimate
#' @param penalty_loadings passed on to `pense::regmest_cv()`
#' @param cache_path path for caching intermediate results.
#' @param alpha a set of different *alpha* values.
#' @param seed seed for the individual cross-validations (the same seed is used for all
#'  *alpha* values).
#' @param cv_k,cv_repl,sparse,max_solutions,algorithm_opts further arguments to `pense::regmest_cv()`.
compute_adammest_cv <- function (..., pense_ridge, bdp, seed, alpha, cache_path,
                                 cv_k, cv_repl, penalty_loadings = NULL, sparse = FALSE,
                                 max_solutions = 1, en_algo_opts = NULL, algorithm_opts,
                                 zeta = NA_real_,
                                 log_indent = 0L) {
  if (!require(pense, quietly = TRUE)) {
    stop('`pense` package not available')
  }
  if (!require('Matrix', quietly = TRUE)) {
    stop('`Matrix` package not available')
  }

  est_name <- if (is.null(penalty_loadings)) { 'MM' } else { 'adaptive MM' }

  if (is.null(en_algo_opts)) {
    en_algo_opts <- en_cd_options()
  }

  if (missing(algorithm_opts)) {
    algorithm_opts <- mm_algorithm_options(en_algorithm_opts = en_algo_opts)
  }

  scale <- try_catch({
    with(pense_ridge,
         sqrt(2 * (objf_value - lambda * (alpha * sum(abs(std_beta)) +
                                            0.5 * (1 - alpha) * sum(std_beta^2)))))
  })

  ridge_starting_point <- try_catch({
    as_starting_point(structure(list(
      alpha = 0,
      estimates = list(pense_ridge[c('intercept', 'beta', 'alpha', 'lambda')])),
      class = "pense_fit"))
  })

  # Run for every `alpha`
  all_alpha_results <- lapply(alpha, function (alpha) {
    print_log("Computing %s estimate for alpha=%0.2f.", est_name, alpha,
              .indent = log_indent)

    try_catch({
      cache <- load_cache('%010d,%010d,%010d,%.2f|%s', seed, cv_k, cv_repl,
                          alpha, hash_num_vec(penalty_loadings),
                          .prefix = 'adammest_', .cache_path = cache_path,
                          .indent = log_indent + 1L)

      if (!is.null(cache$object)) {
        cache$object
      } else {
        # Set the seed so that we get the same CV splits for every alpha
        set.seed(seed)
        alpha_start <- proc.time()[["elapsed"]]
        estimates <- regmest_cv(..., penalty_loadings = penalty_loadings, scale = scale,
                                starting_points = ridge_starting_point,
                                alpha = alpha, sparse = sparse, max_solutions = max_solutions,
                                cv_k = cv_k, cv_repl = cv_repl,
                                algorithm_opts = algorithm_opts)
        alpha_end <- proc.time()[["elapsed"]]

        mmest_alpha_res <- extract_pense_cv_coef(estimates,
                                                 lambda = if (isTRUE(cv_repl > 1)) c('min', 'se')
                                                 else 'min',
                                                 zeta = zeta)
        mmest_alpha_res$duration <- alpha_end - alpha_start
        print_log("Saving results to cache file '%s'.", cache$file, .indent = log_indent + 1L)
        saveRDS(mmest_alpha_res, file = cache$file)

        mmest_alpha_res
      }
    })
  })
  duration_all_alpha <- sum(vapply(all_alpha_results, FUN.VALUE = numeric(1),
                                   FUN = `[[`, 'duration'))

  try_catch({
    # Select the best alpha value
    adamm_ests <- list(cv_min = get_best_estimate(all_alpha_results, 'min'),
                       cv_se = get_best_estimate(all_alpha_results, 'se'))

    if (isFALSE(sparse)) {
      # Make sparse coefficient vectors
      adamm_ests <- lapply(adamm_ests, function (est) {
        if (!is.null(est) && is.numeric(est$beta)) {
          nnz_ind <- which(abs(est$beta) > 1e-12)
          est$beta <- sparseVector(est$beta[nnz_ind], i = nnz_ind, length = length(est$beta))
        }
        return(est)
      })
    }
    adamm_ests$cv_min$duration_all_alpha <- duration_all_alpha
    adamm_ests$cv_se$duration_all_alpha <- duration_all_alpha

    adamm_ests
  })
}

#' Compute adaptive PENSE estimates from a PENSE-Ridge preliminary estimate for a sequence of
#' *alpha* and *zeta* values.
#'
#' This is a two-step procedure, leveraging the function `compute_adapense_cv()` in both steps.
#'
#' 1. Compute a preliminary PENSE-Ridge estimate.
#' 2. Compute the adaptive PENSE estimates.
#'
#' @param lambda_min_ratio can be a vector as long as `zeta_seq`, in which case each
#'   `zeta` will use a different ratio.
#' @param ...,alpha,fit_all,en_algo_opts passed on to `compute_adapense_cv()`.
#' @param zeta_seq a set of different *zeta* values.
#' @param penalty_loadings ignored (they are computed automatically from the preliminary
#' PENSE-Ridge estimate).
compute_adapense_cv_zeta <- function (..., alpha, zeta_seq,
                                      lambda_min_ratio, lambda_min_ratio_prelim = lambda_min_ratio,
                                      fit_all = c('min', 'se'),
                                      penalty_loadings, en_algo_opts = NULL,
                                      log_indent = 0L) {
  if (!missing(penalty_loadings)) {
    warning("`penalty_loadings` cannot be specified as they will be computed from a
            preliminary PENSE-Ridge estimate.")
  }

  print_log('Computing a prelminary PENSE-Ridge estimate.', .indent = log_indent)

  lambda_min_ratio <- if (!missing(lambda_min_ratio)) {
    rep_len(lambda_min_ratio, length(zeta_seq))
  } else {
    NULL
  }

  if (missing(lambda_min_ratio_prelim)) {
    lambda_min_ratio_prelim <- NULL
  }

  prelim_pense <- try_catch({
    compute_adapense_cv(..., alpha = 0, fit_all = c('min', 'se'),
                        lambda_min_ratio = lambda_min_ratio_prelim[[1L]],
                        en_algo_opts = en_lars_options(),
                        log_indent = log_indent + 1L)
  })

  zeta_results <- lapply(seq_along(zeta_seq), function (zi) {
    zeta <- zeta_seq[[zi]]
    print_log("Computing adaptive PENSE estimates with zeta=%0.1f.", zeta,
              .indent = log_indent)
    try_catch({
      pen_loadings <- abs(prelim_pense$cv_min$beta)^(-zeta)
      adapr <- compute_adapense_cv(penalty_loadings = pen_loadings,
                                   alpha = alpha, fit_all = fit_all,
                                   lambda_min_ratio = lambda_min_ratio[[zi]],
                                   en_algo_opts = en_algo_opts, zeta = zeta,
                                   log_indent = log_indent + 1L, ...)
      lapply(adapr, function (x) {
        x$zeta <- zeta
        x$duration_prelim <- prelim_pense$cv_min$duration_all_alpha
        return(x)
      })
    })
  })
  duration_all_zeta <- sum(vapply(zeta_results, FUN.VALUE = numeric(1),
                                  FUN = function (x) x$cv_min$duration_all_alpha))

  try_catch({
    res <- list(cv_min = get_best_estimate(zeta_results, 'cv_min'),
                cv_se = get_best_estimate(zeta_results, 'cv_se'),
                preliminary = prelim_pense)
    res$cv_min$duration_all_zeta <- duration_all_zeta
    res$cv_se$duration_all_zeta <- duration_all_zeta

    res
  })
}


#' Compute adaptive MM estimates from a PENSE-Ridge preliminary estimate for a sequence
#' of *alpha* and *zeta* values.
#'
#' @param lambda_min_ratio can be a vector as long as `zeta_seq`, in which case each
#'   `zeta` will use a different ratio.
#' @param ...,pense_ridge passed on to `compute_adamm_cv()`.
#' @param zeta_seq a set of different *zeta* values.
#' @param penalty_loadings ignored (they are computed automatically from the preliminary
#'   PENSE-Ridge estimate).
compute_adammest_cv_zeta <- function (..., lambda_min_ratio, pense_ridge, zeta_seq,
                                      penalty_loadings, log_indent = 0L) {
  if (!missing(penalty_loadings)) {
    warning("`penalty_loadings` cannot be specified as they will be computed from a preliminary
            PENSE-Ridge estimate.")
  }

  lambda_min_ratio <- if (!missing(lambda_min_ratio)) {
    rep_len(lambda_min_ratio, length(zeta_seq))
  } else {
    NULL
  }

  zeta_results <- lapply(seq_along(zeta_seq), function (zi) {
    zeta <- zeta_seq[[zi]]
    print_log("Computing adaptive MM estimates with zeta=%0.1f.", zeta,
              .indent = log_indent)
    try_catch({
      pen_loadings <- abs(pense_ridge$beta)^(-zeta)
      adamm <- compute_adammest_cv(penalty_loadings = pen_loadings, pense_ridge = pense_ridge,
                                   lambda_min_ratio = lambda_min_ratio[[zi]],
                                   zeta = zeta, log_indent = log_indent + 1L, ...)
      lapply(adamm, function (x) {
        x$zeta <- zeta
        x$duration_prelim <- pense_ridge$duration_all_alpha
        return(x)
      })
    })
  })
  duration_all_zeta <- sum(vapply(zeta_results, FUN.VALUE = numeric(1),
                                  FUN = function (x) x$cv_min$duration_all_alpha))

  try_catch({
    res <- list(cv_min = get_best_estimate(zeta_results, 'cv_min'),
                cv_se = get_best_estimate(zeta_results, 'cv_se'))

    res$cv_min$duration_all_zeta <- duration_all_zeta
    res$cv_se$duration_all_zeta <- duration_all_zeta

    res
  })


}

#' Compute I-LAMM Estimates
#'
#' @param x,y regression data.
#' @param nlambda number of lambda values in the grid.
#' @param ntau number of tau parameters to try.
#' @param cache_path path for caching intermediate results.
#' @param cv_k number of CV folds and train-test splits, the size of the test set is
#'             roughly *n / nrepl*.
#' @param penalty the penalty function to use for `ILAMM::cvNcvxHuberReg()`.
compute_ilamm <- function (x, y, nlambda, cv_k, seed, cache_path, cv_repl = 1,
                           penalty = 'Lasso', log_indent = 0L, ncores = 1) {
  if (!require('Matrix', quietly = TRUE)) {
    stop('`Matrix` package not available')
  }
  if (!requireNamespace('ILAMM', quietly = TRUE)) {
    stop('`ILAMM` package not available')
  }

  if (isTRUE(ncores > 1L)) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    lapply <- function (X, FUN, ...) {
      parallel::clusterApply(cl, x = X, fun = FUN, ...)
    }
  }

  print_log("Computing I-LAMM estimate (%s).", penalty, .indent = log_indent)

  try_catch({
    cache <- load_cache('%010d,%010d,%010d,%s', seed, cv_k, cv_repl, penalty,
                        .prefix = 'ilamm_', .cache_path = cache_path,
                        .indent = log_indent + 1L)

    if (!is.null(cache$object)) {
      cache$object
    } else {
      set.seed(seed)
      start <- proc.time()[["elapsed"]]
      ilr <- ILAMM::cvNcvxHuberReg(x, y, nlambda = nlambda, penalty = penalty, nfolds = cv_k,
                                   intercept = TRUE)

      if (cv_repl > 1L) {
        mae <- lapply(seq_len(cv_repl - 1), function (., .x, .y, nlambda,
                                               penalty, cv_k) {
          perm <- sample.int(length(.y))
          ILAMM::cvNcvxHuberReg(.x[perm, ], .y[perm], nlambda = nlambda,
                                penalty = penalty,
                                nfolds = cv_k, intercept = TRUE)$mae
        }, .x = x, .y = y, nlambda = nlambda, penalty = penalty, cv_k = cv_k)

        mae <- array(c(ilr$mae, unlist(mae, recursive = FALSE, use.names = FALSE)),
                     dim = c(dim(ilr$mae), cv_repl))

        mae_stat <- apply(mae, c(1, 2), function (x) {
          c(mean = mean(x), sd = sd(x))
        })
        best_hp <- arrayInd(which.min(mae_stat['mean', , , drop = TRUE]), dim(mae_stat)[-1L])
        repl_ilr <- ILAMM::ncvxHuberReg(x, y, intercept = TRUE,
                                        penalty = penalty,
                                        lambda = ilr$lambdaSeq[[best_hp[[1]]]],
                                        tau = ilr$tauSeq[[best_hp[[2]]]])
        ilr <- repl_ilr
      }
      end <- proc.time()[["elapsed"]]

      nnz_ind <- which(abs(ilr$beta[-1]) > 1e-12)

      est <- with(ilr, list(
        intercept = beta[[1]],
        beta = sparseVector(beta[1 + nnz_ind], nnz_ind, length(beta) - 1),
        pred_err = min(mae / length(y)),
        duration = end - start
      ))

      print_log("Saving results to cache file '%s'.", cache$file, .indent = log_indent + 1L)
      saveRDS(est, file = cache$file)

      est
    }
  })
}

#' Compute the Classical EN estimator
#'
#' @param bdp desired breakdown-point for the *M-scale of the residuals*.
#' @param nrepl number of CV-folds
comp_glmnet <- function (x, y, alpha, nlambda, cv_k, seed, cache_path, cv_repl = 1,
                         penalty_loadings = NULL, zeta = NA_real_, log_indent = 0L) {
  if (!require('Matrix', quietly = TRUE)) {
    stop('`Matrix` package not available')
  }
  if (!requireNamespace('glmnet', quietly = TRUE)) {
    stop('`glmnet` package not available')
  }

  est_name <- if (is.null(penalty_loadings)) { 'LS-EN' } else { 'adaptive LS-EN' }

  print_log("Computing classical %s estimates.", est_name, .indent = log_indent)

  en_res <- try_catch({
    cache <- load_cache('%010d,%010d,%010d|%s|%s', seed, cv_k, cv_repl,
                        hash_num_vec(alpha, 2),
                        hash_num_vec(penalty_loadings),
                        .prefix = 'lsen_', .cache_path = cache_path,
                        .indent = log_indent + 1L)

    if (!is.null(cache$object)) {
      cache$object
    } else {
      pen_fact <- if (!is.null(penalty_loadings)) {
        penalty_loadings
      } else {
        rep.int(1, ncol(x))
      }

      en_res <- lapply(alpha, function (alpha) {
        print_log("Computing %s estimate for alpha=%.2f and zeta=%f.", est_name, alpha,
                  zeta, .indent = log_indent + 1L)

        set.seed(seed)
        start <- proc.time()[["elapsed"]]
        cv_res <- glmnet::cv.glmnet(x = x, y = y, nlambda = nlambda, alpha = alpha,
                                    type.measure = 'mae', nfolds = cv_k,
                                    penalty.factor = pen_fact,
                                    parallel = FALSE)

        # glmnet may not compute the solution for all requested values
        cv_res$cvm <- c(rep_len(cv_res$cvm[[1L]], nlambda - length(cv_res$cvm)),
                        cv_res$cvm)
        cv_res$lambda <- c(rep_len(cv_res$lambda[[1L]], nlambda - length(cv_res$lambda)),
                           cv_res$lambda)

        if (cv_repl > 1) {
          cvm <- vapply(seq_len(cv_repl - 1), FUN.VALUE = cv_res$cvm, FUN = function (.) {
            repl_res_cvm <- glmnet::cv.glmnet(x = x, y = y, nlambda = nlambda, alpha = alpha,
                                              type.measure = 'mae', nfolds = cv_k,
                                              parallel = FALSE)$cvm

            c(rep_len(repl_res_cvm[[1L]], nlambda - length(repl_res_cvm)), repl_res_cvm)
          })
          cvm <- cbind(cv_res$cvm, cvm)
          cv_res$cvm <- rowMeans(cvm)
          cv_res$cvsd <- apply(cvm, 1, sd)
          cv_res$cvup <- cv_res$cvm + cv_res$cvsd
          cv_res$cvlo <- cv_res$cvm - cv_res$cvsd
        }
        end <- proc.time()[["elapsed"]]

        glmnet_coefs_min <- coef(cv_res, 'lambda.min')
        glmnet_coefs_se <- coef(cv_res, 'lambda.1se')
        list(min = list(intercept = glmnet_coefs_min@x[[1]],
                        beta = sparseVector(glmnet_coefs_min@x[-1], glmnet_coefs_min@i[-1],
                                            length = glmnet_coefs_min@Dim[[1]] - 1L),
                        alpha = alpha,
                        zeta = zeta,
                        pred_err = cv_res$cvm[[cv_res$index['min', 1L]]]),
             se = list(intercept = glmnet_coefs_se@x[[1]],
                       beta = sparseVector(glmnet_coefs_se@x[-1], glmnet_coefs_se@i[-1],
                                           length = glmnet_coefs_se@Dim[[1]] - 1L),
                       alpha = alpha,
                       zeta = zeta,
                       pred_err = cv_res$cvm[[cv_res$index['1se', 1L]]]),
             duration = end - start)
      })

      print_log("Saving results to cache file '%s'.", cache$file, .indent = log_indent + 1L)
      saveRDS(en_res, file = cache$file)

      en_res
    }
  })
  duration_all_alpha <- sum(vapply(en_res, FUN.VALUE = numeric(1),
                                   FUN = `[[`, 'duration'))

  try_catch({
    res <- list(cv_min = get_best_estimate(en_res, 'min'),
                cv_se = get_best_estimate(en_res, 'se'))

    res$cv_min$duration_all_alpha <- duration_all_alpha
    res$cv_se$duration_all_alpha <- duration_all_alpha

    res
  })
}

#' Compute the Classical EN estimator
#'
#' @param bdp desired breakdown-point for the *M-scale of the residuals*.
#' @param nrepl number of CV-folds
comp_glmnet_zeta <- function (x, ..., alpha, zeta_seq, zeta, penalty_loadings,
                              log_indent = 0L) {
  if (!missing(penalty_loadings)) {
    warning("`penalty_loadings` cannot be specified as they will be computed from a
            preliminary PENSE-Ridge estimate.")
  }
  if (!missing(zeta)) {
    warning("`zeta` cannot be specified. Use `zeta_seq`.")
  }

  print_log('Computing a prelminary LS-Ridge estimate.', .indent = log_indent)

  prelim_en <- try_catch(comp_glmnet(x = x, ..., alpha = 0, log_indent = log_indent + 1L,
                                     penalty_loadings = NULL))

  x_sd <- apply(x, 2, sd)

  zeta_results <- lapply(zeta_seq, function (zeta) {
    print_log("Computing adaptive EN estimates with zeta=%0.1f.", zeta, .indent = log_indent)

    try_catch({
      pen_loadings <- as.numeric(abs(prelim_en$cv_min$beta * x_sd)^(-zeta))
      if (any(!is.finite(pen_loadings))) {
        pen_loadings[!is.finite(pen_loadings)] <- max(pen_loadings[is.finite(pen_loadings)])
      }

      adaen <- comp_glmnet(x = x, penalty_loadings = pen_loadings,
                           alpha = alpha, zeta = zeta, log_indent = log_indent + 1L,
                           ...)
      lapply(adaen, function (x) {
        x$zeta <- zeta
        x$duration_preliminary <- prelim_en$cv_min$duration_all_alpha
        return(x)
      })
    })
  })
  duration_all_zeta <- sum(vapply(zeta_results, FUN.VALUE = numeric(1),
                                  FUN = function (x) x$cv_min$duration_all_alpha))

  try_catch({
    res <- list(cv_min = get_best_estimate(zeta_results, 'cv_min'),
                cv_se = get_best_estimate(zeta_results, 'cv_se'),
                preliminary = prelim_en)

    res$cv_min$duration_all_zeta <- duration_all_zeta
    res$cv_se$duration_all_zeta <- duration_all_zeta

    res
  })
}

#' Compute the intercept-only model using the M-location estimate
#'
#'
#' @param x,y regression data.
#' @param bdp breakdown point of the estimator
compute_intonly <- function (x, y, bdp) {
  if (!require(pense, quietly = TRUE)) {
    stop('`pense` package not available')
  }

  try_catch(list(
    intercept = mlocscale(y, bdp = bdp)[['location']],
    beta = sparseVector(numeric(0L), integer(0L), ncol(x))
  ))
}

#' Recursively evaluate all estimates on test data.
#' Computes the prediction errors and the tau-size of the prediction errors.
#'
#' @param estimate a (possibly nested) list of estimate object. An estimate object is defined
#'                 by having fields `intercept` and `beta`.
#' @param test_x a matrix of predictor values from the test set.
#' @param test_x a vector of response values from the test set.
#' @return an updated list of estimate objects. Each estimate will have fields
#'         `test_pred_err` and `test_pred_err_tau_size` added to it.
evaluate_estimate <- function (estimate, test_x, test_y) {
  if (is.list(estimate)) {
    if (!is.null(estimate$beta)) {
      try_catch({
        estimate$test_pred <- with(estimate, as.numeric(test_x %*% beta) + intercept)
        estimate$test_pred_err <- test_y - estimate$test_pred
        estimate$test_pred_err_tau_size <- tau_size(estimate$test_pred_err)
        estimate
      })
    }
    return(lapply(estimate, evaluate_estimate, test_x = test_x, test_y = test_y))
  }
  # Don't do anything
  estimate
}

#' Print a log message, if .verbose = TRUE
print_log <- function (..., .indent = 0L) {
  indent_spaces <- tryCatch(strrep('  ', max(as.integer(.indent), 0L)),
                            warning = function (w) '',
                            error = function (e) '')
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S | "), indent_spaces,
      sprintf(...), '\n', sep = '')
}

#' Extract coefficients from the `pense_cv()` result.
#'
#' @param ... additional arguments passed on to `coef()` and `prediction_performance()`.
extract_pense_cv_coef <- function (pense_cv_res, lambda, zeta = NA_real_, ...) {
  coefs <- lapply(lambda, function (l) {
    coef_list <- coef(pense_cv_res, lambda = l, concat = FALSE, sparse = TRUE)
    c(coef_list,
      list(pred_err = with(pense_cv_res$cvres, cvavg[which.min(abs(lambda - coef_list$lambda))]),
           cv_curve = pense_cv_res$cvres,
           zeta = zeta,
           bdp = pense_cv_res$bdp))
  })
  names(coefs) <- as.character(lambda)
  coefs
}

#' Get the best estimate based on the value of `what` for the specified estimate level
get_best_estimate <- function (x, level, what = 'pred_err', p = 0L) {
  best_ind <- which.min(vapply(x, FUN.VALUE = numeric(1L), function (y) {
    max(y[[level]][[what]], 0)
  }))
  x[[best_ind]][[level]]
}


#' Load a result object from a cache file.
#'
#' The cache name is created by hashing the given format string,
#' applied to the additional arguments, and concatenating the hash with the given
#' pre- and suffix under the provided cache path.
#'
#' @param format_str,... the format string and optional arguments passed on to [sprintf()].
#' @param .cache_path path to the cache files.
#' @param .prefix,.suffix pre- and suffix applied to the cache name.
#' @param .indent the indentation for log messages
#' @return a list with components `file`, `hash` and (if available) `object`.
load_cache <- function (format_str, ..., .prefix = '', .suffix = '', .cache_path, .indent = 0L) {
  if (!requireNamespace('digest', quietly = TRUE)) {
    stop("Package `digest` is required")
  }
  hash <- digest::digest(sprintf(format_str, ...), algo = 'md5', serialize = FALSE)
  cache_name <- file.path(.cache_path, paste(.prefix, hash, .suffix, '.rds', sep = ''))
  if (file.exists(cache_name)) {
    print_log("Fetching object from cache '%s'", cache_name, .indent = .indent)
    list(file = cache_name, hash = hash, object = readRDS(cache_name))
  } else {
    list(file = cache_name, hash = hash)
  }
}

hash_num_vec <- function (x, min_digits = 3) {
  if (!requireNamespace('digest', quietly = TRUE)) {
    stop("Package `digest` is required")
  }
  if (is.null(x)) {
    'NULL'
  } else {
    digits <- max(min_digits, abs(floor(log10(abs(median(x)))) - min_digits))
    if (!is.finite(digits)) {
      digits <- min_digits
    }
    digest::digest(paste(sprintf(sprintf('%%.%df', digits), x), collapse = ','),
                   serialize = FALSE, algo = 'md5')
  }
}

try_catch <- function (expr, env = parent.frame()) {
  if (!requireNamespace('rlang', quietly = TRUE)) {
    stop("Package `rlang` is required.")
  }
  if (!requireNamespace('methods', quietly = TRUE)) {
    stop("Package `methods` is required.")
  }
  rlang::try_fetch(eval(expr, envir = env),
                   error = function (cnd) {
                     rlang::try_fetch(
                       rlang::abort("Error:", parent = cnd),
                       error = function (cnd) {
                         print(cnd, simplify = 'branch')
                         NULL
                       })
                     NULL
                   })
}
