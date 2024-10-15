## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_unknown #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param penalties list of penalties to be tested
#' @param n_penalties number of penalty values to be tested (if penalties not provided, default is 30)
#' @param min_ratio ratio between min and max penalty to be tested  (if penalties not provided, default is 0.05)
#' @param models uderlying NB_fixed_blocks models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_fixed_blocks_sparse <- R6::R6Class(
  classname = "NB_fixed_blocks_sparse",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses
    Y  = NULL,
    #' @field X the design matrix
    X = NULL,
    #' @field C the group matrix
    C = NULL,
    #' @field penalties list of penalties values for omegaQ sparsity
    penalties = NULL,
    #' @field n_penalties number of penalty values
    n_penalties = NULL,
    #' @field min_ratio ratio between min and max penalty to be tested
    min_ratio = NULL,
    #' @field models list of NB_fixed_Q models corresponding to each nb_block value
    models = NULL,
    #' @field verbose say whether information should be given about the optimization
    verbose = NULL,

    #' @description Create a new [`NB_fixed_blocks_sparse`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param C group matrix.
    #' @param penalties list of penalties on the network density
    #' @param n_penalties number of penalty values
    #' @param min_ratio ratio between min and max penalty to be tested
    #' @return A new [`nb_fixed_blocks_sparse`] object
    initialize = function(Y, X, C, penalties = NULL,  n_penalties = 30,
                          min_ratio = 0.05, verbose=TRUE) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y, X and C must be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y <- Y
      self$X <- X
      self$C <- C
      if(!is.null(penalties)){
        self$n_penalties <- length(penalties)
        self$penalties <- penalties[order(penalties)]
      }else{
        self$n_penalties <- n_penalties
        self$min_ratio   <- min_ratio
        init_model <- normal_block(Y, X, C, verbose = FALSE)
        init_model$optimize(5)
        sigmaQ    <- solve(init_model$model_par$omegaQ)
        max_pen   <- max(abs(sigmaQ[upper.tri(sigmaQ, diag = FALSE)]))
        self$penalties <- 10^seq(log10(max_pen), log10(max_pen* self$min_ratio), len = self$n_penalties)
        }
      self$verbose <- verbose
    },


    #' @description optimizes an NB_fixed_Q object for each value of Q
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    optimize = function(niter = 100, threshold = 1e-4) {
      self$models <- furrr::future_map(seq_along(self$models), function(m) {
        model <- self$models[[m]]
        if(self$verbose) cat("\t penalty =", self$models[[m]]$sparsity, "          \r")
        flush.console()
        model$optimize(niter, threshold)
        if(m < self$n_penalties){
          self$models[[m + 1]]$update(B      = model$model_par$B,
                                      dm1    = model$model_par$dm1,
                                      omegaQ = model$model_par$omegaQ,
                                      gamma  = model$posterior_par$gamma,
                                      mu     = model$posterior_par$mu)
        }
        model
      }, .options = furrr_options(seed=TRUE))
    },

    #' @description returns the NB_fixed_Q model corresponding to given Q
    #' @param penalty penalty asked by user
    #' @return A NB_fixed_Q object with given value Q
    get_model = function(penalty) {
      if(!(penalty %in% self$penalties)) {
        closest_penalty <-  self$penalties[[which.min(abs(self$penalties - penalty))]]
        paste0("No model with this penalty in the collection. Returning model with closest penalty : ", closest_penalty,  " Collection penalty values can be found via $penalties")
      }
      penalty_rank <- which(sort(self$penalties) == penalty)
      return(self$models[[penalty_rank]])
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "BIC", "AIC" or "loglik" (-loglik so that criterion to be minimized)
    #' "loglik" is the default criterion
    #' @return a [`NB_fixed_Q`] object
    get_best_model = function(crit = c("loglik", "BIC", "AIC", "ICL")) {
      crit <- match.arg(crit)
      stopifnot(!anyNA(self$criteria[[crit]]))
      id <- 1
      if (length(self$criteria[[crit]]) > 1) {
        id <- which.min(self$criteria[[crit]])
      }
      model <- self$models[[id]]$clone()
      model
    },

    #' @param criterion criterion to plot
    #' @param type char for line type (see plot.default)
    #' @param log char for logarithmic axes (see plot.default)
    #' @param neg boolean plot negative criterion (useful when log="y")
    #' @description plots given criterion as a function of Q
    plot_criterion = function(criterion = "loglik", type = "b", log = "",
                              neg = FALSE) {
      neg <- ifelse(neg, -1, 1)
      y   <- self$criteria[[criterion]]
      plot(self$penalties, neg * y, type=type, log=log, ylab = paste0(ifelse(neg == -1, "-", ""), criterion),
           xlab = "penalties")
      axis(1, at = seq_along(y), labels = self$penalties)
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
    #' @field Q number of blocks
    Q = function() ncol(self$C),
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC, AIC) for the collection of models
    criteria = function() purrr::map(self$models, "criteria") %>% purrr::reduce(rbind)
  )
)


#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param penalties list of penalties to be tested
#' @param n_penalties number of penalty values to be tested (if penalties not provided, default is 30)
#' @param min_ratio ratio between min and max penalty to be tested  (if penalties not provided, default is 0.05)
#' @param models uderlying NB_fixed_Q models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_fixed_blocks_diagonal_sparse <- R6::R6Class(
  classname = "NB_fixed_blocks_diagonal_sparse",
  inherit = NB_fixed_blocks_sparse,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NB_fixed_blocks_diagonal_sparse`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param C group matrix.
    #' @param penalties list of penalties on the network density
    #' @param n_penalties number of penalty values
    #' @param min_ratio ratio between min and max penalty to be tested
    #' @return A new [`NB_fixed_blocks_diagonal_sparse`] object
    initialize = function(Y, X, C, penalties = NULL,  n_penalties = 30,
                          min_ratio = 0.05, verbose=TRUE) {
      super$initialize(Y, X, C, penalties, n_penalties,
                       min_ratio, verbose)
      # instantiates an NB_fixed_blocks_diagonal model for each Q in nb_blocks
      self$models <- map(self$penalties[order(self$penalties)],
                         function(penalty) {
                           model <- NB_fixed_blocks_diagonal$new(self$Y, self$X,
                                                                 self$C,
                                                                 sparsity = penalty)
                         })
    }
  )
)

#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param penalties list of penalties to be tested
#' @param n_penalties number of penalty values to be tested (if penalties not provided, default is 30)
#' @param min_ratio ratio between min and max penalty to be tested  (if penalties not provided, default is 0.05)
#' @param models uderlying NB_fixed_Q models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_fixed_blocks_spherical_sparse <- R6::R6Class(
  classname = "NB_fixed_blocks_spherical_sparse",
  inherit = NB_fixed_blocks_sparse,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NB_fixed_blocks_sparse`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param C group matrix.
    #' @param penalties list of penalties on the network density
    #' @param n_penalties number of penalty values
    #' @param min_ratio ratio between min and max penalty to be tested
    #' @return A new [`nb_fixed_blocks_sparse`] object
    initialize = function(Y, X, C, penalties = NULL,  n_penalties = 30,
                          min_ratio = 0.05, verbose=TRUE) {
      super$initialize(Y, X, C, penalties, n_penalties,
                       min_ratio, verbose)
      # instantiates an NB_fixed_blocks_spherical model for each Q in nb_blocks
      self$models <- map(self$penalties[order(self$penalties)],
                         function(penalty) {
                           model <- NB_fixed_blocks_spherical$new(self$Y, self$X,
                                                                 self$C, penalty)
                         })
    }
  )
)
