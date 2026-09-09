#' Normal-block model
#'
#' Fit a normal-block model with a variational or heuristic algorithm
#' @param data NormalBlockData object, contains the matrix of responses (Y, n x p) and the design matrix (X, n x d), must be created with NormalBlockData$new.
#' @param blocks either a integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * q matrix (for indicating block membership when its known)
#' @param sparsity either TRUE to run the optimization for different sparsity penalty values
#' OR float to run model with a single sparsity penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param model which model family to fit: "var" (the default) structures the
#' clustering in the latent covariance (see [NormalBlockVarBase]), "mean"
#' structures it in the mean (mu_i = C B' X_i, see [NormalBlockMeanBase]).
#' The mean-block family covers the same collections as the variance-block
#' one ([NormalBlockMeanCollectionClusters],
#' [NormalBlockMeanCollectionSparsity] and their crossing,
#' [NormalBlockMeanCollectionClustersSparsity]) -- each q simply fit
#' independently, without the variance-block family's SBM-path shortcut.
#' Zero-inflation is supported ([ZINormalBlockMeanKnownClusters],
#' [ZINormalBlockMeanUnknownClusters], and collections of those over a range
#' of q), with a diagonal or spherical Sigma only, hence not along a sparsity
#' path -- see [NormalBlockMeanBase].
#' @param control a list-like structure for detailed control on parameters should be
#' generated with NB_control().
#' @return an R6 object with one of the model classes (or a collection of model objects).
#' @examples
#' ## Normal Data
#' ex_data <- generate_normal_block_var_data(n=100, p=30, d=1, q=3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' my_normal_block <- normal_block(data, blocks = 1:6)
#' my_normal_block$plot(c("deviance", "BIC", "ICL"))
#' Y_hat <- my_normal_block$get_best_model()$fitted
#' plot(data$Y, Y_hat, log = "xy"); abline(0,1)
#' ## Normal Data with Zero Inflation
#' ex_data_zi <- generate_normal_block_var_data(n=50, p=50, d=1, q=3, kappa = rep(0.5,50))
#' zidata <- NormalBlockData$new(ex_data_zi$Y, ex_data_zi$X)
#' my_normal_block <- normal_block(zidata, blocks = 1:6, zero_inflation = TRUE)
#' ## Mean-Block model (clustering in the mean rather than the covariance)
#' ex_mean <- generate_normal_block_mean_data(n=50, p=20, d=1, q=3)
#' mean_data <- NormalBlockData$new(ex_mean$Y, ex_mean$X)
#' my_mean_block <- normal_block(mean_data, blocks = 3, model = "mean")
#'
#' @export
normal_block <- function(data,
                         blocks,
                         sparsity = 0,
                         zero_inflation = FALSE,
                         control = NB_control(),
                         model = c("var", "mean")) {
  ## Recovering the requested model from the function arguments
  model <- match.arg(model)
  stopifnot(is.numeric(blocks) | is.matrix(blocks))
  if (is.numeric(blocks) && !is.matrix(blocks) && length(blocks) > 1 && anyDuplicated(blocks))
    stop("`blocks` looks like a vector of cluster labels; a known clustering must be passed as a ", "p x q indicator matrix (e.g. model.matrix(~ 0 + factor(labels))). A numeric vector is ", "read as the range of cluster counts to explore, so its values must be distinct.", call. = FALSE)
  stopifnot(is.null(control$sparsity_weights) | is.matrix(control$sparsity_weights))
  if (!is.null(control$sparsity_weights)) stopifnot(isSymmetric(control$sparsity_weights))
  if (is.list(control$clustering_init)) stopifnot(length(control$clustering_init) == length(blocks))
  stopifnot(
    "clustering_init = 'best_of_inits' is not supported together with sparsity = TRUE (a sparsity path warm-starts a single clustering across all penalties)" =
      !(uses_best_of_inits(control) && isTRUE(sparsity))
  )

  fit <- get_model(data, blocks, sparsity = sparsity,
                   zero_inflation = zero_inflation,
                   model = model,
                   control = control)

  ## Estimation/optimization
  if (control$verbose) cat("Fitting a", fit$who_am_I, "\n")

  if (uses_best_of_inits(control) && !is.null(fit$best_of_inits)) {
    fit <- fit$best_of_inits(control = control)
  } else {
    fit$optimize(control)
  }

  ## Finishing
  if (control$verbose) cat("\nDONE\n")

  fit
}

#' NB_control
#'
#' Control the model settings and various optimization parameters
#'
#' @param niter number of iterations in model optimization
#' @param fixed_point_niter number of sweeps of the tau update for
#' Normal-Block-Mean with unknown clusters. Each sweep visits the rows of
#' tau sequentially and maximizes the ELBO exactly in each, so it can never
#' decrease the ELBO.
#' @param threshold loglikelihood / elbo threshold under which optimization stops
#' @param sparsity_weights weights with which the penalty should be applied in case
#' sparsity is required, non-0 values on the diagonal mean diagonal shall be
#' penalized too (default is non-penalized diagonal and 1s off-diagonal)
#' @param sparsity_penalties list of penalties the user wants to test, other parameters
#' are only used if penalties is not specified
#' @param n_sparsity_penalties number of penalties to test.
#' @param min_ratio ratio for sparsity between max penalty (0 edge penalty) and min penalty to test
#' @param fixed_tau whether tau should be fixed at clustering_init during optimization
#' useful for calls to fixed_q models in stability_selection
#' @param clustering_init how to obtain the initial clustering of the q
#' unknown blocks: a heuristic name ("ward2", "kmeans", "sbm" or
#' "spectral"), an actual clustering (a vector of labels or a p x q indicator
#' matrix, or a list of either per q for a collection), or "best_of_inits" to
#' try several heuristics per model and keep the best-ELBO fit (see
#' [NormalBlockVarBase]'s `best_of_inits()`; not supported with `sparsity =
#' TRUE`). Default `NULL`, resolved per model family at fit time: "ward2"
#' for variance-block models, "kmeans" for mean-block models ("ward2" was
#' benchmarked substantially worse there -- see [NormalBlockMeanBase]).
#' See `inst/methods_initialization_and_refine.md` for the
#' heuristics' rationale, why no single one dominates, and how this interacts
#' with `refine` (below).
#' @param verbose telling if information should be printed during optimization
#' @param noise_covariance shape of the residual covariance. Variance-block
#' models accept "diagonal" (variable-specific) or "spherical" (common);
#' mean-block models, whose Sigma is the full p x p residual covariance, also
#' accept "full". Default `NULL`, resolved per model family at fit time to
#' "diagonal" -- except for a mean-block model with `sparsity > 0`, which
#' implies "full" since a penalty on a diagonal precision matrix would have
#' nothing to act on (explicitly asking for both is an error). The mean-block
#' "diagonal"/"spherical" variants need no matrix inversion, hence no
#' `n > p` requirement, and they select the number of clusters markedly
#' better than a full Sigma once p approaches n: the p(p+1)/2 covariance
#' parameters otherwise drown the mean structure BIC/ICL are weighing. Use
#' "full" when the residual associations are themselves of interest.
#' @param heuristic whether to use the heuristic approach (moment-based, no (V)EM
#' recursion) instead of the full (V)EM. Default is FALSE. In heuristic mode, no
#' likelihood/ELBO is computed, so `entropy`, `loglik`, `BIC`, `ICL` and `EBIC`
#' are all `NA` on the resulting model.
#' @param refine for [NormalBlockVarCollectionClusters] only: whether
#' `optimize()` should automatically call `refine()` afterwards. Default
#' `FALSE` since it adds real cost; call `collection$refine()` directly at
#' any point afterwards for the same effect without setting this.
#' @return A named list of parameters to pass to [normal_block()]'s `control`
#' argument.
#' @export
NB_control <- function(
    niter                = 500,
    threshold            = 1e-4,
    fixed_point_niter    = 5,
    sparsity_weights     = NULL,
    sparsity_penalties   = NULL,
    n_sparsity_penalties = 30,
    min_ratio            = 0.01,
    fixed_tau            = FALSE,
    clustering_init      = NULL,
    verbose              = TRUE,
    heuristic            = FALSE,
    noise_covariance     = NULL,
    refine               = FALSE) {

  if (!is.null(sparsity_weights))
    stopifnot(all(is.matrix(sparsity_weights), isSymmetric(sparsity_weights)))
  if (!is.null(noise_covariance))
    stopifnot("noise_covariance must be one of 'diagonal', 'spherical' or (mean-block models only) 'full'" =
                length(noise_covariance) == 1 &&
                noise_covariance %in% c("diagonal", "spherical", "full"))
  if (is.character(clustering_init) && length(clustering_init) == 1) {
    stopifnot("clustering_init, when given as a single string, must name a known heuristic ('kmeans', 'ward2', 'sbm' or 'spectral'), or be 'best_of_inits' -- otherwise pass an actual clustering (a vector of labels, a p x q indicator matrix, or a list of either for a collection over several q values)" =
                clustering_init %in% c("kmeans", "ward2", "sbm", "spectral", "best_of_inits"))
  }

  structure(list(niter                = niter                ,
                 threshold            = threshold            ,
                 fixed_point_niter    = fixed_point_niter    ,
                 sparsity_weights     = sparsity_weights     ,
                 sparsity_penalties   = sparsity_penalties   ,
                 n_sparsity_penalties = n_sparsity_penalties ,
                 min_ratio            = min_ratio            ,
                 fixed_tau            = fixed_tau            ,
                 clustering_init      = clustering_init      ,
                 verbose              = verbose              ,
                 heuristic            = heuristic            ,
                 noise_covariance     = noise_covariance     ,
                 refine               = refine               ))
}


#' Create a Normal-Block Model Object
#'
#' Creates the appropriate normal-block model (or collection of models)
#' depending on the parametrization.
#' @param blocks either an integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * q matrix (for indicating block membership when its known)
#' @param sparsity boolean to say whether the model should have a changing penalty
#' OR float to run model with a single penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param model which model family to fit, "var" (the default) or "mean" --
#' see [normal_block()]
#' @param control a list-like structure for detailed control on parameters should be
#' generated with NB_control() for collections of sparse models
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @export
get_model <- function(data,
                      blocks,
                      sparsity = 0,
                      zero_inflation = FALSE,
                      control = NB_control(),
                      model = c("var", "mean")) {
  model <- match.arg(model)

  changing_sparsity <- isTRUE(sparsity)
  unknown_q_list    <- !is.matrix(blocks) && length(blocks) > 1
  is_collection     <- changing_sparsity || unknown_q_list

  if (model == "mean") {
    stopifnot(
      "zero-inflation is not available for a mean-block sparsity path: it needs a full Sigma, which zero-inflated mean-block models do not carry" =
        !(zero_inflation && changing_sparsity)
    )
    if (changing_sparsity) {
      return(if (unknown_q_list)
               NormalBlockMeanCollectionClustersSparsity$new(data, blocks, control = control)
             else
               NormalBlockMeanCollectionSparsity$new(data, blocks, control = control))
    }
    if (unknown_q_list) {
      return(NormalBlockMeanCollectionClusters$new(data, blocks, zero_inflation,
                                                   sparsity = sparsity, control = control))
    }
    class_name <- if (is.matrix(blocks)) "NormalBlockMeanKnownClusters" else "NormalBlockMeanUnknownClusters"
    if (zero_inflation) class_name <- paste0("ZI", class_name)
    myClass <- eval(str2lang(class_name))
    return(myClass$new(data, blocks, sparsity, control = control))
  }

  if (is_collection) {
    class_name <- if (changing_sparsity && unknown_q_list) {
      "NormalBlockVarCollectionClustersSparsity"
    } else if (changing_sparsity) {
      "NormalBlockVarCollectionSparsity"
    } else {
      "NormalBlockVarCollectionClusters"
    }
  } else {
    class_name <- if (is.matrix(blocks)) "NormalBlockVarKnownClusters" else "NormalBlockVarUnknownClusters"
    if (zero_inflation) class_name <- paste0("ZI", class_name)
  }

  myClass <- eval(str2lang(class_name))
  if (is_collection) {
    if (changing_sparsity) {
      model <- myClass$new(data, blocks, zero_inflation, control = control)
    } else {
      model <- myClass$new(data, blocks, zero_inflation, sparsity, control = control)
    }
  } else {
    model <- myClass$new(data, blocks, sparsity, control = control)
  }
  model
}
