#' Normal-block model
#'
#' Fit a normal-block model with a variational algorithm
#' @param Y observations data, n observations of p categories, dim n*p.
#' @param X covariates data, n observations of d covariates, dim n*d.
#' @param nb_blocks number of blocks (if fixed) or list of possible number of blocks otherwise
#' @param blocks, if known, under the form of a p * Q matrix
#' @param sparsity sparsity factor to apply to blocks covariance matrix (or possibly list of such values with same length as nb_blocks)
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood / elbo threshold under which optimization stops
#' @return an R6 object with class [`NB`] or [`NB_unknown`] or [`NB_unknown_ZI`]
#' @export
normal_block <- function(Y, X, nb_blocks = NULL, blocks = NULL, sparsity = 0,
                         zero_inflation = FALSE, niter = 100, threshold = 1e-4) {
  if(!is.null(blocks)){
    if(zero_inflation){
      model <- NB_fixed_blocks_zi$new(Y, X, blocks, sparsity = sparsity)
      model$optimize(niter, threshold)
      return(model)
    }else{
      model <- NB_fixed_blocks$new(Y, X, blocks, sparsity = sparsity)
      model$optimize(niter, threshold)
      return(model)
    }
  }

  if(is.null(nb_blocks)){nb_blocks = 1:ncol(Y)}
  if(length(nb_blocks) == 1){
    if(zero_inflation){
      model <- NB_fixed_Q_zi$new(Y, X, nb_blocks, sparsity = sparsity)
      model$optimize(niter, threshold)
      return(model)
    }else{
      model <- NB_fixed_Q$new(Y, X, nb_blocks, sparsity = sparsity)
      model$optimize(niter, threshold)
      return(model)
    }
  }

  if(zero_inflation){
    models <- NB_unknown_zi$new(Y, X, nb_blocks, sparsity = sparsity)
    models$optimize(niter, threshold)
    return(models)
  }else{
    models <- NB_unknown$new(Y, X, nb_blocks, sparsity = sparsity)
    models$optimize(niter, threshold)
    return(models)
  }
}
