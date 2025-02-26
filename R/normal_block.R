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
#' @param noise_variance variance can be variable specific ("diagonal", the default) or common ("spherical")
#' @param heuristic weither to use heuristical approach or not. Default is FALSE
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
