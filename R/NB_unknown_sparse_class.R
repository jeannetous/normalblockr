## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_unknown #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param nb_blocks list of number of blocks values to be tested
#' @param penalties list of penalties to be tested
#' @param n_penalties number of penalty values to be tested (if penalties not provided, default is 30)
#' @param min_ratio ratio between min and max penalty to be tested  (if penalties not provided, default is 0.05)
#' @param models uderlying NB_unknown models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_unknown_sparse <- R6::R6Class(
  classname = "NB_unknown_sparse",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses
    Y  = NULL,
    #' @field X the design matrix
    X = NULL,
    #' @field nb_blocks list of number of blocks values to be tested
    nb_blocks = NULL,
    #' @field penalties list of penalties values for omegaQ sparsity
    penalties = NULL,
    #' @field n_penalties number of penalty values
    n_penalties = NULL,
    #' @field min_ratio ratio between min and max penalty to be tested
    min_ratio = NULL,
    #' @field models list of NB_unknown models corresponding to each nb_block value
    models = NULL,
    #' @field verbose say whether information should be given about the optimization
    verbose = NULL,
    #' @field latest_niter latest niter value used for optimization
    latest_niter = NULL,
    #' @field latest_threshold latest threshold value used for optimization
    latest_threshold = NULL,

    #' @description Create a new [`NB_unknown_sparse`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param Q number of blocks
    #' @param control structured list of specific parameters for sparsity
    #' @return A new [`nb_fixed_Q_sparse`] object
    initialize = function(Y, X, nb_blocks, control = NB_unknown_sparse_param(),
                          verbose=TRUE) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y, X and C must be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y <- Y
      self$X <- X
      self$nb_blocks <- nb_blocks
      penalties <- control$penalties
      n_penalties <- control$n_penalties
      min_ratio <- control$min_ratio
      if(!is.null(penalties)){
        self$n_penalties <- length(penalties)
        self$penalties   <- penalties[order(penalties)]
      }else{
        self$n_penalties <- n_penalties
        self$min_ratio   <- min_ratio
      }
      self$verbose <- verbose
    },


    #' @description optimizes an NB_unknown object for each value of Q
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    optimize = function(niter = 100, threshold = 1e-4) {
      self$latest_niter     <- niter
      self$latest_threshold <- threshold

      self$models <- furrr::future_map(seq_along(self$models), function(m) {
        model <- self$models[[m]]
        if(self$verbose) cat("\t penalty =", self$models[[m]]$penalty, "          \r")
        flush.console()
        model$optimize(niter, threshold)
        model
      }, .options = furrr_options(seed=TRUE))
    },

    #' @description returns a collection of NB_unknown models corresponding to given penalty
    #' or one single model if Q is also given
    #' @param penalty penalty asked by user
    #' @return A NB_unknown_sparse object with given value penalty
    get_model = function(Q, penalty = NA) {
      if(!(Q %in% self$nb_blocks)) {
        stop("No such model in the collection. Acceptable parameter values can be found via $nb_blocks")
      }
      Q_rank <- which(sort(self$nb_blocks) == Q)
      model  <- self$models[[Q_rank]]$clone()
      if(is.na(penalty)){
        return(model)
      }else{
        if(!(penalty %in% model$penalties)) {
          penalty <-  model$penalties[[which.min(abs(model$penalties - penalty))]]
          cat(paste0("No model with this penalty in the collection. Returning model with closest penalty : ", penalty,  " Collection penalty values can be found via $penalties \n"))
        }
        penalty_rank <- which(sort(model$penalties) == penalty)
        return(model$models[[penalty_rank]])
      }
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "BIC", "AIC" or "loglik" (-loglik so that criterion to be minimized)
    #' "loglik" is the default criterion
    #' @return a [`NB_unknown`] object
    get_best_model = function(crit = c("loglik", "BIC", "AIC", "ICL")) {
      crit <- match.arg(crit)
      stopifnot(!anyNA(self$criteria[[crit]]))
      id <- 1
      if (length(self$criteria[[crit]]) > 1) {
        id       <- which.min(self$criteria[[crit]])
        best_pen <- self$criteria$penalty[[id]]
        best_Q   <- self$criteria$Q[[id]]
      }
      model <- self$get_model(best_Q, best_pen)$clone()
      model
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples
    n = function() nrow(self$Y),
    #' @field p number of responses per sample
    p = function() ncol(self$Y),
    #' @field d number of variables (dimensions in X)
    d = function() ncol(self$X),
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC, AIC) for the collection of models
    criteria = function() crit <- purrr::map(self$models, "criteria") %>% purrr::reduce(rbind)
  )
)


#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param nb_blocks number of blocks
#' @param penalties list of penalties to be tested
#' @param n_penalties number of penalty values to be tested (if penalties not provided, default is 30)
#' @param min_ratio ratio between min and max penalty to be tested  (if penalties not provided, default is 0.05)
#' @param models uderlying NB_unknown models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_unknown_diagonal_sparse <- R6::R6Class(
  classname = "NB_unknown_diagonal_sparse",
  inherit = NB_unknown_sparse,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NB_unknown_diagonal_sparse`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param C group matrix.
    #' @param control structured list of specific parameters for sparsity
    #' @return A new [`NB_unknown_diagonal_sparse`] object
    initialize = function(Y, X, nb_blocks, control = NB_unknown_diagonal_sparse_param(),
                          verbose=TRUE) {
      super$initialize(Y, X, nb_blocks, control, verbose)
      # instantiates an NB_unknown_diagonal model for each penalty
      # For now NB_unknown_spherical not defined --> it's all NB_unknown_diagonal
      self$models <- map(self$nb_blocks[order(self$nb_blocks)],
                         function(Q) {
                           model <- NB_fixed_Q_diagonal_sparse$new(self$Y, self$X,
                                                                   Q, control,
                                                                   verbose = F)
                         })
    }
  )
)

#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param list of number of blocks values to be tested
#' @param penalties list of penalties to be tested
#' @param n_penalties number of penalty values to be tested (if penalties not provided, default is 30)
#' @param min_ratio ratio between min and max penalty to be tested  (if penalties not provided, default is 0.05)
#' @param models uderlying NB_unknown models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_unknown_spherical_sparse <- R6::R6Class(
  classname = "NB_unknown_spherical_sparse",
  inherit = NB_unknown_sparse,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NB_unknown_sparse`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param nb_blocks list of number of blocks values to be tested.
    #' @param control structured list of specific parameters to handle sparsity
    #' @return A new [`nb_fixed_Q_sparse`] object
    initialize = function(Y, X, nb_blocks, control, verbose=TRUE) {
      super$initialize(Y, X, nb_blocks, control, verbose)
      # instantiates an NB_unknown_spherical model for each Q in nb_blocks
      self$models <- map(self$nb_blocks[order(self$nb_blocks)],
                         function(Q) {
                           model <- NB_fixed_Q_spherical_sparse$new(self$Y, self$X,
                                                                   Q, self$penalties,
                                                                   self$n_penalties,
                                                                   self$min_ratio,
                                                                   verbose = F)
                         })
    }
  )
)


#' NB_unknown_sparse_param
#'
#' Generates control parameters for the NB_fixed_blocks_sparse class
NB_unknown_sparse_param <- function(penalties = NULL, n_penalties = 30,
                                    min_ratio = 0.05){
  structure(list(penalties         = penalties        ,
                 n_penalties       = n_penalties      ,
                 min_ratio         = min_ratio        ))
}


#' NB_unknown_diagonal_sparse_param
#'
#' Generates control parameters for the NB_fixed_blocks_sparse class
#' @export
NB_unknown_diagonal_sparse_param <- function(penalties = NULL, n_penalties = 30,
                                             min_ratio = 0.05){
  NB_unknown_sparse_param(penalties, n_penalties, min_ratio)
}

#' NB_unknown_diagonal_spherical_param
#'
#' Generates control parameters for the NB_fixed_blocks_sparse class
#' @export
NB_unknown_spherical_sparse_param <- function(penalties = NULL, n_penalties = 30,
                                              min_ratio = 0.05){
  NB_unknown_sparse_param(penalties, n_penalties, min_ratio)
}

