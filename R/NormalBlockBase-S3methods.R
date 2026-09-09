## ===================================================================================
##
## PUBLIC S3 METHODS FOR NormalBlockBase (all families)
##
## ===================================================================================

#' @title Check if an Object is a Normal-Block Model
#' @description Checks if a model is of class [NormalBlockBase()], i.e. any
#' normal-block model (variance-block or mean-block family).
#' @param object An R object.
#' @return A boolean telling whether object inherits from the NormalBlockBase class.
#' @export
isNB <- function(object) {inherits(object, "NormalBlockBase")}

#' @title Extract Model Coefficients
#' @description Extract coefficients from a normal-block model.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return A matrix of coefficients extracted from the model.
#' @export
coef.NormalBlockBase <- function(object, ...){
  stopifnot(isNB(object))
  object$model_par$B
}

#' @title Extract the Covariance Matrix
#' @description Extract the covariance matrix `Omega^-1`: between latent
#' blocks (q x q) for the variance-block models, between variables (p x p)
#' for the mean-block models.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return The covariance matrix, of size q x q or p x p depending on the model.
#' @importFrom stats sigma
#' @export
sigma.NormalBlockBase <- function(object, ...){
  stopifnot(isNB(object))
  chol2inv(chol(object$model_par$Omega))
}

#' @title Extract Fitted Values
#' @description Extract fitted values from a normal-block model.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return A matrix of fitted values extracted from the object.
#' @export
fitted.NormalBlockBase <- function(object, ...){
  stopifnot(isNB(object))
  object$fitted
}


#' @title Predict Method for Variance-Block Models
#' @description Predicts observations Y for new covariates X. Specific to the
#' variance-block family: the mean-block models have their own formula
#' (mu = C B' X).
#' @param object An object of class NormalBlockVarBase.
#' @param new_X New set of covariates.
#' @param ... not used, only here for S3 compatibility
#' @return A n*p prediction matrix for new observations
#' @export
predict.NormalBlockBase <- function(object, new_X, ...){
  stopifnot(isNB(object))
  object$predict(new_X)
}

#' @title Extract Log-Likelihood of a Normal-Block Model
#' @description Returns the (variational) log-likelihood of a fitted
#' normal-block model as a `"logLik"` object, compatible with [stats::AIC()]
#' and [stats::BIC()].
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return An object of class `"logLik"`. The numeric value is the
#' log-likelihood or its variational lower bound (ELBO). Attributes `df` and
#' `nobs` hold the number of parameters and observations.
#' @importFrom stats logLik
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' logLik(model)
logLik.NormalBlockBase <- function(object, ...) {
  stopifnot(isNB(object))
  structure(object$loglik, class = "logLik", df = object$nb_param, nobs = object$n)
}

#' @title Bayesian Information Criterion for a Normal-Block Model
#' @description Extracts the (variational) BIC of a fitted normal-block
#' model, computed as `deviance + log(n) * nb_param` (lower is better).
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return A scalar: the (variational) BIC.
#' @importFrom stats BIC
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' BIC(model)
BIC.NormalBlockBase <- function(object, ...) {
  stopifnot(isNB(object))
  object$BIC
}

#' @title Print a Normal-Block Model
#' @description Print a short summary of a fitted normal-block model: model
#' type, goodness-of-fit criteria, and the useful fields/methods to explore
#' it further.
#' @param x An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns `x`.
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' print(model)
print.NormalBlockBase <- function(x, ...) {
  stopifnot(isNB(x))
  x$print()
  invisible(x)
}

#' @title Summarize a Normal-Block Model
#' @description Summarizes a fitted normal-block model: model type,
#' goodness-of-fit criteria, cluster sizes, and the density of the inferred
#' network between blocks.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return An object of class `summary.NormalBlockBase` (a list with the
#' model's `who_am_I`, `criteria`, `cluster_sizes` and network `density`),
#' printed with a dedicated [print.summary.NormalBlockBase()] method.
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' summary(model)
summary.NormalBlockBase <- function(object, ...) {
  stopifnot(isNB(object))
  ## the network lives on Omega: q x q (variance-block) or p x p (mean-block)
  size <- ncol(object$model_par$Omega)
  res <- list(
    who_am_I      = object$who_am_I,
    criteria      = object$criteria,
    cluster_sizes = object$cluster_sizes,
    n_edges       = object$n_edges,
    density       = if (size > 1) object$n_edges / choose(size, 2) else NA
  )
  class(res) <- "summary.NormalBlockBase"
  res
}

#' @title Print a Normal-Block Model Summary
#' @description Print method for objects returned by [summary.NormalBlockBase()].
#' @param x An object of class `summary.NormalBlockBase`.
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns `x`.
#' @export
print.summary.NormalBlockBase <- function(x, ...) {
  cat("A", x$who_am_I, ".\n")
  cat("===========================================================================\n")
  print(as.data.frame(round(x$criteria, digits = 3), row.names = ""))
  cat("===========================================================================\n")
  cat("* Cluster sizes:", paste(x$cluster_sizes, collapse = ", "), "\n")
  cat("* Network: ", x$n_edges, " edge(s), density = ",
      if (is.na(x$density)) NA else round(x$density, 3), "\n", sep = "")
  invisible(x)
}

#' @title Plot a Normal-Block Model
#' @description Plots the evolution of the objective (log-likelihood or
#' ELBO) across the (V)EM iterations of the last call to `optimize()`, see
#' `$plot_loglik()`.
#' @param x An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns the [ggplot2::ggplot] object; called for its
#' side effect of plotting.
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' plot(model)
plot.NormalBlockBase <- function(x, ...) {
  stopifnot(isNB(x))
  p <- x$plot()
  print(p)
  invisible(p)
}
