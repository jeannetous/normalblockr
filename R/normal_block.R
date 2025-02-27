#' NB_control
#'
#' @param sparsity_weights weights with which penalty should be applied in case
#' sparsity is required, non-0 values on the diagonal mean diagonal shall be
#' penalized too (default is non-penalized diagonal and 1s off-diagonal)
#' @param penalties list of penalties the user wants to test, other parameters
#' are only used if penalties is not specified
#' @param n_penalties number of penalties to test.
#' @param min_ratio ratio between max penalty (0 edge penalty) and min penalty to test
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
    sparsity_weights  = NULL,
    penalties         = NULL,
    n_penalties       = 30,
    min_ratio         = 0.01,
    fixed_tau         = FALSE,
    clustering_init   = NULL,
    verbose           = TRUE,
    heuristic         = FALSE,
    noise_covariance  = c("diagonal", "spherical"),
    clustering_approx = c("residuals", "covariance")) {
  structure(list(sparsity_weights = sparsity_weights ,
                 penalties        = penalties        ,
                 n_penalties      = n_penalties      ,
                 min_ratio        = min_ratio        ,
                 fixed_tau        = fixed_tau        ,
                 clustering_init  = clustering_init  ,
                 verbose          = verbose          ,
                 heuristic        = heuristic        ,
                 noise_covariance = match.arg(noise_covariance),
                 clustering_approx = match.arg(clustering_approx)))
}


#' Creates appropriate new normal block model depending on the parametrization
#' @param blocks either an integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * Q matrix (for indicating block membership when its known)
#' @param sparsity boolean to say whether the model should have a changing penalty
#' OR float to run model with a single penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param noise_cov character the type of covariance for the noise: either diagonal of spherical
#' @param control a list-like structure for detailed control on parameters should be
#' generated with normal_block_control() for collections of sparse models
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @export
get_model <- function(data, blocks, sparsity = 0,
                      zero_inflation = FALSE, control= NB_control()){
  block_class   <- ifelse(is.matrix(blocks), "_fixed_blocks",
                        ifelse(length(blocks) > 1, "_unknown_Q", "_fixed_Q"))
  sparse_class  <- ifelse(typeof(sparsity) == "logical" & sparsity, "_changing_sparsity", "")
  zi_class      <- ifelse(block_class != "unknown_Q" & sparse_class != "_changing_sparsity" &
                         zero_inflation, "_zi",  "")
  is_collection <- ifelse( sparse_class == "_changing_sparsity" | block_class == "_unknown_Q",
                           TRUE, FALSE)

  myClass <- eval(str2lang(paste0("NB", zi_class, block_class, sparse_class)))
  if(is_collection){
    model   <- myClass$new(data, blocks, zero_inflation, sparsity, control = control)
  }else{model   <- myClass$new(data, blocks, sparsity, control = control)}
  model
}

