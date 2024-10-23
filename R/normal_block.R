#' Normal-block model
#'
#' Fit a normal-block model with a variational algorithm
#' @param Y observations data, n observations of p categories, dim n*p.
#' @param X covariates data, n observations of d covariates, dim n*d.
#' @param blocks either a integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * Q matrix (for indicating block membership when its known)
#' @param sparsity boolean to say whether the model should be sparsified (several penalties tested)
#' OR float to run model with a single penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param noise_cov character the type of covariance for the noise: either diagonal of spherical
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood / elbo threshold under which optimization stops
#' @param verbose telling if information should be printed during optimization
#' @param control a list-like structure for detailed control on parameters
#' @param optimize boolean stating whether the model should just be instantiated (F) or optimized too (T)
#' @return an R6 object with class [`NB`] or [`NB_unknown`] or [`NB_unknown_ZI`]
#' @examples
#' data("example_data")
#' Y <- example_data$Y
#' X <- example_data$X
#' my_normal_block <- normal_block(Y, X, blocks = 1:6)
#' \dontrun{
#' my_normal_block$plot_criterion("loglik")
#' my_normal_block$plot_criterion("BIC")
#' my_normal_block$plot_criterion("ICL")
#' Y_hat <- my_normal_block$get_best_model()$fitted
#' plot(Y, Y_hat); abline(0,1)
#' }
#' @export
normal_block <- function(Y, X, blocks,
                         sparsity = FALSE,
                         zero_inflation = FALSE,
                         noise_cov = c("diagonal","spherical"),
                         niter = 100, threshold = 1e-4,
                         verbose=TRUE, control = NULL,
                         optimize = TRUE) {

  ## Recovering the requested model from the function arguments
  stopifnot(is.numeric(blocks) | is.matrix(blocks))
  noise_cov <- match.arg(noise_cov)
  block_class <- ifelse(is.matrix(blocks), "fixed_blocks",
                        ifelse(length(blocks) > 1, "unknown", "fixed_Q"))
  zi_class <- ifelse(zero_inflation, "_zi",  "")
  # sparse_class <- ifelse(sparsity, "_sparse", "") # a enlever a terme

### FIX until all models have their spherical variant & their sparse variant
  noise_cov <- ifelse((!zero_inflation & block_class == "fixed_blocks") |(typeof(sparsity) == "logical" & sparsity), paste0("_",noise_cov), "")
  sparse_class <- ifelse(typeof(sparsity) == "logical" & sparsity, "_sparse", "")
  ## Instantiating model
  myClass <- eval(str2lang(paste0("NB_", block_class, zi_class, noise_cov, sparse_class)))
  if (!is.null(control)) {
    model <- myClass$new(Y, X, blocks, sparsity, control = control)
  }else{
    if(typeof(sparsity) == "logical"){model <- myClass$new(Y, X, blocks)
    }else{model <- myClass$new(Y, X, blocks, sparsity)}
  }

  ## Estimation/optimization
  if(optimize){
    if(verbose)
      cat("Fitting a", sub('.', '', noise_cov),
          sub('.', '',sparse_class),
          ifelse(zero_inflation, "zero-inflated",  ""),
          "normal-block model with", block_class, "...\n")
    model$optimize(niter, threshold)
  }else{
    if(verbose){
      cat("Instantiating a", sub('.', '', noise_cov),
          sub('.', '',sparse_class),
          ifelse(zero_inflation, "zero-inflated",  ""),
          "normal-block model with", block_class, "...\n")
    }
  }
  ## Finishing
  if(verbose) cat("\n DONE\n")
  model
}
