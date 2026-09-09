## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  SEQUENTIAL MEAN-THEN-VARIANCE FIT ##################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Cluster variables in the mean, then in the residual covariance
#'
#' Fits the two model families one after the other: a mean-block model
#' ([NormalBlockMeanBase]) groups the variables by how they respond to the
#' covariates, then a variance-block model ([NormalBlockVarBase]) groups the
#' residuals of that fit by how they co-vary. The two answer different
#' questions and generally return unrelated partitions, so running both is
#' often more informative than choosing one.
#'
#' The second stage uses an intercept-only design on purpose: the covariate
#' effects have already been removed by the first stage.
#'
#' This is a heuristic two-stage estimator, not a joint model. On simulated
#' data carrying two genuinely distinct structures it recovers both exactly
#' (see `inst/mean_block_analyses/sequential_mean_then_variance.R`).
#'
#' @param data a [NormalBlockData] object
#' @param blocks_mean number of clusters for the mean-block stage: an integer,
#' a vector of integers to explore, or a p x q indicator matrix
#' @param blocks_var idem for the variance-block stage, run on the residuals
#' @param crit criterion used to pick a model when a range is explored,
#' "ICL" (the default) or "BIC"
#' @param zero_inflation whether Y carries structural zeros. Both stages are
#' then zero-inflated, sharing the same mask: the second one would otherwise
#' re-derive it from residuals, which are never exactly zero.
#' @param control_mean control list for the mean-block stage, see [NB_control()]
#' @param control_var control list for the variance-block stage
#' @return an object of class `normal_block_sequential`, a list with the fitted
#' `mean` and `var` models and the residual matrix `residuals` handed from one
#' stage to the other.
#' @examples
#' ex   <- generate_normal_block_mean_data(n = 80, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' fit  <- normal_block_sequential(data, blocks_mean = 3, blocks_var = 2)
#' fit
#' @importFrom stats fitted
#' @export
normal_block_sequential <- function(data, blocks_mean, blocks_var,
                                    crit = c("ICL", "BIC"),
                                    zero_inflation = FALSE,
                                    control_mean = NB_control(verbose = FALSE),
                                    control_var  = NB_control(verbose = FALSE)) {
  stopifnot("data must be a NormalBlockData object" = inherits(data, "NormalBlockData"))
  crit <- match.arg(crit)

  pick <- function(fit) if (isNBcollection(fit)) fit$get_best_model(crit) else fit

  fit_mean  <- pick(normal_block(data, blocks_mean, model = "mean",
                                 zero_inflation = zero_inflation, control = control_mean))
  residuals <- data$Y * matrix(data$Y_scale, data$n, data$p, byrow = TRUE) - fitted(fit_mean)

  intercept <- matrix(1, data$n, 1, dimnames = list(NULL, "(Intercept)"))
  data_var  <- NormalBlockData$new(residuals, intercept, zeros = data$zeros)
  fit_var   <- pick(normal_block(data_var, blocks_var,
                                 zero_inflation = zero_inflation, control = control_var))

  structure(list(mean = fit_mean, var = fit_var, residuals = residuals),
            class = "normal_block_sequential")
}

#' @title Print a Sequential Mean-then-Variance Fit
#' @description Reports both stages and, when \pkg{aricode} is available, how
#' related the two partitions are.
#' @param x an object of class `normal_block_sequential`
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns `x`.
#' @export
print.normal_block_sequential <- function(x, ...) {
  cat("A sequential mean-then-variance normal-block fit\n")
  cat("===========================================================================\n")
  cat("  mean-block stage    :", x$mean$q, "clusters --", x$mean$who_am_I, "\n")
  cat("  variance-block stage:", x$var$q, "clusters --", x$var$who_am_I, "\n")
  if (requireNamespace("aricode", quietly = TRUE))
    cat("  ARI between the two partitions:",
        round(aricode::ARI(x$mean$clustering, x$var$clustering), 3),
        "\n  (near 0 means the two structures are unrelated, as is usual)\n")
  cat("===========================================================================\n")
  cat("* Useful fields\n")
  cat("    $mean, $var (both ordinary fitted models), $residuals \n")
  invisible(x)
}
