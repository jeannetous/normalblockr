#' Normal-block model
#'
#' Fit a normal-block model with a variational or heuristic algorithm
#' @param data contains the matrix of responses (Y, n x p) and the design matrix (X, n x d)."
#' @param blocks either a integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * Q matrix (for indicating block membership when its known)
#' @param sparsity either TRUE to run the optimization for different sparsity penalty values
#' OR float to run model with a single sparsity penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param control a list-like structure for detailed control on parameters should be
#' generated with NB_control().
#' @return an R6 object with one of the NB classes
#' @examples
#' ## Normal Data
#' ex_data <- generate_normal_block_data(n=50, p=50, d=1, Q=3)
#' data <- NBData$new(ex_data$Y, ex_data$X)
#' my_normal_block <- normal_block(data, blocks = 1:6)
#' \dontrun{
#' my_normal_block$plot(c("deviance", "BIC", "ICL"))
#' Y_hat <- my_normal_block$get_best_model()$fitted
#' plot(data$Y, Y_hat, log = "xy"); abline(0,1)
#' }
#' ## Normal Data with Zero Inflation
#' ex_data_zi <- generate_normal_block_data(n=50, p=50, d=1, Q=3, kappa = rep(0.5,50))
#' zidata <- NBData$new(ex_data_zi$Y, ex_data_zi$X)
#' my_normal_block <- normal_block(zidata, blocks = 1:6, zero_inflation = TRUE)
#'
#' @export
normal_block <- function(data,
                         blocks,
                         sparsity = 0,
                         zero_inflation = FALSE,
                         control = NB_control()) {

  ## Recovering the requested model from the function arguments
  stopifnot(is.numeric(blocks) | is.matrix(blocks))
  stopifnot(is.null(control$sparsity_weights) | is.matrix(control$sparsity_weights))
  if (!is.null(control$sparsity_weights)) stopifnot(isSymmetric(control$sparsity_weights))
  if (is.list(control$clustering_init)) stopifnot(length(control$clustering_init) == length(blocks))

  model <- get_model(data, blocks, sparsity = sparsity,
                     zero_inflation = zero_inflation,
                     control = control)

  ## Estimation/optimization
  if (control$verbose) cat("Fitting a", model$who_am_I, "\n")

  model$optimize(control)

  ## Finishing
  if (control$verbose) cat("\nDONE\n")

  model
}

#' NB_control
#'
#' Control the model settings and various optimization parameters
#'
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood / elbo threshold under which optimization stops
#' @param sparsity_weights weights with which the penalty should be applied in case
#' sparsity is required, non-0 values on the diagonal mean diagonal shall be
#' penalized too (default is non-penalized diagonal and 1s off-diagonal)
#' @param sparsity_penalties list of penalties the user wants to test, other parameters
#' are only used if penalties is not specified
#' @param n_sparsity_penalties number of penalties to test.
#' @param min_ratio ratio for sparsity between max penalty (0 edge penalty) and min penalty to test
#' @param fixed_tau whether tau should be fixed at clustering_init during optimization
#' useful for calls to fixed_Q models in stability_selection
#' @param clustering_init proposal of initial value for clustering, when Q is
#' unknown, can be a list with one clustering for each Q value
#' @param verbose telling if information should be printed during optimization
#' @param noise_covariance variance can be variable specific ("diagonal", the default) or common ("spherical")
#' @param heuristic weather to use heuristic approach or not. Default is FALSE
#' @param clustering_approx to use for clustering with heuristic inference method
#' @export
NB_control <- function(
    niter                = 100,
    threshold            = 1e-4,
    sparsity_weights     = NULL,
    sparsity_penalties   = NULL,
    n_sparsity_penalties = 30,
    min_ratio            = 0.01,
    fixed_tau            = FALSE,
    clustering_init      = NULL,
    verbose              = TRUE,
    heuristic            = FALSE,
    noise_covariance     = c("diagonal", "spherical"),
    clustering_approx    = c("ward2", "kmeans", "sbm")) {

  if (!is.null(sparsity_weights))
    stopifnot(all(is.matrix(sparsity_weights), isSymmetric(sparsity_weights)))

  structure(list(niter                = niter                ,
                 threshold            = threshold            ,
                 sparsity_weights     = sparsity_weights     ,
                 sparsity_penalties   = sparsity_penalties   ,
                 n_sparsity_penalties = n_sparsity_penalties ,
                 min_ratio            = min_ratio            ,
                 fixed_tau            = fixed_tau            ,
                 clustering_init      = clustering_init      ,
                 verbose              = verbose              ,
                 heuristic            = heuristic            ,
                 noise_covariance     = match.arg(noise_covariance),
                 clustering_approx    = match.arg(clustering_approx)))
}


#' Creates appropriate new normal block model depending on the parametrization
#' @param blocks either an integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * Q matrix (for indicating block membership when its known)
#' @param sparsity boolean to say whether the model should have a changing penalty
#' OR float to run model with a single penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param control a list-like structure for detailed control on parameters should be
#' generated with normal_block_control() for collections of sparse models
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @export
get_model <- function(data,
                      blocks,
                      sparsity = 0,
                      zero_inflation = FALSE,
                      control = NB_control()) {

  sparse_class  <- ifelse(typeof(sparsity) == "logical" &
                            sparsity,
                          "_changing_sparsity",
                          "")

  block_class   <- ifelse(
    sparse_class == "_changing_sparsity" &
      (is.matrix(blocks) | length(blocks) == 1),
    "",
    ifelse(
      is.matrix(blocks),
      "_fixed_blocks",
      ifelse(length(blocks) > 1, "_unknown_Q", "_fixed_Q")
    )
  )

  is_collection <- ifelse(sparse_class == "_changing_sparsity" | block_class == "_unknown_Q", TRUE, FALSE)

  class_name <- paste0("NB", block_class, sparse_class)
  if (!is_collection & zero_inflation) class_name <- paste0("ZI", class_name)

  myClass <- eval(str2lang(class_name))
  if (is_collection) {
    if (sparse_class == "_changing_sparsity") {
      model   <- myClass$new(data, blocks, zero_inflation, control = control)
    } else{
      model   <- myClass$new(data, blocks, zero_inflation, sparsity, control = control)
    }
  } else{
    model   <- myClass$new(data, blocks, sparsity, control = control)
  }
  model
}
