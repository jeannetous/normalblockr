#' Creates appropriate new normal block model depending on the parametrization
#' @param blocks either a integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * Q matrix (for indicating block membership when its known)
#' @param sparsity boolean to say whether the model should be sparsified (several penalties tested)
#' OR float to run model with a single penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param noise_cov character the type of covariance for the noise: either diagonal of spherical
#' @param control a list-like structure for detailed control on parameters should be
#' generated with either NB_param() or NB_sparse_param() for collections of sparse models
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
    if(is.null(control)){control <- NB_sparse_param()}
    model   <- myClass$new(Y, X, blocks = blocks, zero_inflation = zero_inflation,
                           noise_cov = noise_cov,
                           control = control, verbose=TRUE)
  }
  return(model)
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
#' @param verbose telling if information should be printed during optimization
#' @param control a list-like structure for detailed control on parameters should be
#' generated with either NB_param() or NB_sparse_param() for collections of sparse models
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
                         verbose=TRUE, control = NULL) {
  ## Recovering the requested model from the function arguments
  stopifnot(is.numeric(blocks) | is.matrix(blocks))
  stopifnot(is.null(control$sparsity_weights) | is.matrix(control$sparsity_weights))
  if(!is.null(control$sparsity_weights)) stopifnot(isSymmetric(control$sparsity_weights))

  noise_cov <- match.arg(noise_cov)
  block_class <- ifelse(is.matrix(blocks), "fixed_blocks",
                        ifelse(length(blocks) > 1, "unknown", "fixed_Q"))
  sparse_class <- ifelse(typeof(sparsity) == "logical" & sparsity, "_sparse", "")


  model <- get_model(Y, X, blocks, sparsity = sparsity,
                    zero_inflation = zero_inflation,
                    noise_cov = noise_cov,
                    control = control)
  ## Estimation/optimization
  if(verbose){
    cat("Fitting a", noise_cov,
        sub('.', '',sparse_class),
        ifelse(zero_inflation, "zero-inflated",  ""),
        "normal-block model with", block_class,
        ifelse(block_class == "unknown", " blocks", ""), "...\n")
  }
  model$optimize(niter, threshold)

  ## Finishing
  if(verbose) cat("\n DONE\n")
  model
}
