## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_unknown_sparse ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param blocks either group matrix C or number of blocks Q
#' @param zero_inflation boolean to specify whether data is zero-inflated
#' @param noise_cov character the type of covariance for the noise: either diagonal of spherical
#' @param control structured list of parameters to handle sparsity control
#' @param models underlying models for each penalty
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
    #' @field blocks list of number of blocks values to be tested.
    blocks = NULL,
    #' @field zero_inflation boolean stating whether model is zero_inflated
    zero_inflation = NULL,
    #' @field noise_cov character the type of covariance for the noise: either diagonal of spherical
    noise_cov = NULL,
    #' @field n_penalties number of penalty values
    n_penalties = NULL,
    #' @field min_ratio ratio between min and max penalty to be tested
    min_ratio = NULL,
    #' @field sparsity_weights weights with which penalty should be applied
    sparsity_weights = NULL,
    #' @field models list of NB_fixed_Q models corresponding to each nb_block value
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
    #' @param blocks group matrix or number of blocks.
    #' @param zero_inflation boolean to specify whether data is zero-inflated
    #' @param noise_cov character the type of covariance for the noise: either diagonal of spherical
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`nb_fixed_blocks_sparse`] object
    initialize = function(Y, X, blocks, zero_inflation = F,
                          noise_cov = "diagonal", control = NB_sparse_param(),
                          verbose=TRUE) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y and X must be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y      <- Y
      self$X      <- X
      self$blocks <- blocks

      self$zero_inflation   <- zero_inflation
      self$noise_cov        <- noise_cov
      self$n_penalties      <- control$n_penalties
      self$min_ratio        <- control$min_ratio
      self$sparsity_weights <- control$sparsity_weights
      self$verbose <- verbose
      self$models <- map(self$blocks[order(self$blocks)],
                         function(Q) {
                           model <- NB_sparse$new(self$Y, self$X, Q,
                                                  zero_inflation = self$zero_inflation,
                                                  noise_cov = self$noise_cov,
                                                  control = control, verbose = F)
                         })
    },


    #' @description optimizes an NB_fixed_Q object for each value of Q
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    optimize = function(niter = 100, threshold = 1e-4) {

      self$latest_niter     <- niter
      self$latest_threshold <- threshold

      self$models <- furrr::future_map(seq_along(self$models), function(m) {
        model <- self$models[[m]]
        if(self$verbose) cat("\t Q =", self$models[[m]]$Q, "          \r")
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
  ## PRIVATE MEMBERS ------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    stab_path = NULL # a field to store the stability path,
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
    nb_blocks = function() self$blocks,
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC, AIC) for the collection of models
    criteria = function(){
      crit <- purrr::map(self$models, "criteria") %>% purrr::reduce(rbind)
      crit},
    #' @field penalties list of penalties used for each Q
    penalties = function(){
      self$criteria %>%
        dplyr::group_by(Q) %>%
        dplyr::summarize(penalties = paste(round(penalty, 2), collapse = ", "))
    }

  )
)
