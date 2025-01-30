#' Creates appropriate new normal block model depending on the parametrization
#' @param blocks either a integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * Q matrix (for indicating block membership when its known)
#' @param sparsity boolean to say whether the model should be sparsified (several penalties tested)
#' OR float to run model with a single penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param noise_cov character the type of covariance for the noise: either diagonal of spherical
#' @param control a list-like structure for detailed control on parameters should be
#' generated with normal_block_param() for collections of sparse models
#' @param Y response matrix
#' @param X design matrix
#' #' @examples
#' myModel <- get_model(blocks = 2)
#' @export
get_model <- function(Y, X, blocks, sparsity = FALSE,
                      zero_inflation = FALSE,
                      noise_cov = c("diagonal","spherical"),
                      control= NULL){

  noise_cov <- match.arg(noise_cov)
  block_class <- ifelse(is.matrix(blocks), "fixed_blocks",
                        ifelse(length(blocks) > 1, "unknown", "fixed_Q"))
  zi_class <- ifelse(zero_inflation, "_zi",  "")
  noise_cov_class <- ifelse((block_class == "fixed_blocks" | block_class == "fixed_Q"), paste0("_",noise_cov), "")
  sparse_class <- ifelse(typeof(sparsity) == "logical" & sparsity, "_sparse", "")
  ## Instantiating model
  if(sparse_class != "_sparse"){
    myClass <- eval(str2lang(paste0("NB_", block_class, zi_class, noise_cov_class, sparse_class)))
    if (!is.null(control)) {
      model <- myClass$new(Y, X, blocks, sparsity, control = control)
    }else{
      if(typeof(sparsity) == "logical"){model <- myClass$new(Y, X, blocks)
      }else{model <- myClass$new(Y, X, blocks, sparsity)}
    }
  }else{
    myClass <- eval(str2lang(paste0("NB", ifelse(block_class == "unknown", "_unknown", ""), sparse_class)))
    if(is.null(control)){control <- normal_block_param()}
    model   <- myClass$new(Y, X, blocks = blocks,
                           zero_inflation = zero_inflation,
                           noise_cov = noise_cov,
                           control = control)
  }
  model
}



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
#' @param control a list-like structure for detailed control on parameters should be
#' generated with normal_block_param() for collections of sparse models
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
                         control = normal_block_param()) {
  ## Recovering the requested model from the function arguments
  stopifnot(is.numeric(blocks) | is.matrix(blocks))
  stopifnot(is.null(control$sparsity_weights) | is.matrix(control$sparsity_weights))
  if(!is.null(control$sparsity_weights)) stopifnot(isSymmetric(control$sparsity_weights))
  noise_cov <- match.arg(noise_cov)

  model <- get_model(Y, X, blocks, sparsity = sparsity,
                    zero_inflation = zero_inflation,
                    noise_cov = noise_cov,
                    control = control)
  ## Estimation/optimization
  if(control$verbose) cat("Fitting a", model$who_am_I, "\n")

  model$optimize(niter, threshold)

  ## Finishing
  if(control$verbose) cat("\nDONE\n")

  model
}

#' normal_block_param
#'
#' @param sparsity_weights weights with which penalty should be applied in case
#' sparsity is required, non-0 values on the diagonal mean diagonal shall be
#' penalized too (default is non-penalized diagonal and 1s off-diagonal)
#' @param penalties list of penalties the user wants to test, other parameters
#' are only used if penalties is not specified
#' @param n_penalties number of penalties to test.
#' @param min_ratio ratio between max penalty (0 edge penalty) and min penalty to test
#' @param fixed_tau whether tau should be fixed at clustering_init during optimization
#' @param verbose telling if information should be printed during optimization
#' @param clustering_init proposal of initial value for tau , for when fixed_tau = TRUE
#' useful for calls to fixed_Q models in stability_selection
#' Generates control parameters for NB sparse models
#' @export
normal_block_param <- function(
    sparsity_weights = NULL,
    penalties        = NULL,
    n_penalties      = 30,
    min_ratio        = 0.01,
    fixed_tau        = FALSE,
    clustering_init  = NULL,
    verbose          = TRUE) {
  structure(list(sparsity_weights  = sparsity_weights ,
                 penalties         = penalties        ,
                 n_penalties       = n_penalties      ,
                 min_ratio         = min_ratio        ,
                 fixed_tau         = fixed_tau        ,
                 clustering_init   = clustering_init   ,
                 verbose           = verbose))
}
