## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_unknown_zi #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param nb_blocks list of number of blocks values to be tested
#' @param models uderlying NB_fixed_Q_zi models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_unknown_zi <- R6::R6Class(
  classname = "NB_unknown_zi",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses
    Y  = NULL,
    #' @field X the matrix of covariates
    X = NULL,
    #' @field nb_blocks list of Q values to be tested for the number of blocks
    nb_blocks = NULL,
    #' @field sparsity penalty on the network density
    sparsity = NULL,
    #' @field models list of NB_fixed_Q_zi models corresponding to each nb_block value
    models = NULL,
    #' @field verbose say whether information should be given about the optimization
    verbose = NULL,

    #' @description Create a new [`NB_unknown_zi`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param sparsity penalty on the network density
    #' @param threshold loglikelihood threshold under which optimization stops
    #' @return A new [`nb_fixed`] object
    initialize = function(Y, X, nb_blocks, sparsity = 0, verbose = TRUE) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y, X and C must be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      if (length(nb_blocks) != length(unique(nb_blocks))) {
        stop("each nb_blocks value can only be present once in nb_blocks")
      }
      self$Y <- Y
      self$X <- X
      if (length(sparsity) == 1) sparsity <- rep(sparsity, length(nb_blocks))
      stopifnot(all.equal(length(sparsity), length(nb_blocks)))
      self$sparsity <- sparsity
      self$nb_blocks <- nb_blocks
      self$verbose   <- verbose

      # instantiates an NB_fixed_Q_zi model for each Q in nb_blocks
      self$models <- map2(order(self$nb_blocks), self$sparsity[order(self$nb_blocks)],
                                 function(block_rank, sparsity_sorted ) {
                                   model <- NB_fixed_Q_zi$new(self$Y, self$X,
                                                           nb_blocks[[block_rank]],
                                                           sparsity_sorted)
                                 })
    },

    #' @description optimizes an NB_fixed_Q_zi object for each value of Q
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    optimize = function(niter = 100, threshold = 1e-4) {
      self$models <- furrr::future_map(self$models, function(model) {
        if(self$verbose) cat("\tnumber of blocks =", model$Q, "          \r")
        model$optimize(niter, threshold)
        model
      }, .options = furrr_options(seed=TRUE))
    },

    #' @description returns the NB_fixed_Q_zi model corresponding to given Q
    #' @param Q number of blocks asked by user
    #' @return A NB_fixed_Q_zi object with given value Q
    get_model = function(Q) {
      if(!(Q %in% self$nb_blocks)) {
        stop("No such model in the collection. Acceptable parameter values can be found via $nb_blocks")
        }
      Q_rank <- which(sort(self$nb_blocks) == Q)
      return(self$models[[Q_rank]])
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "BIC", "AIC" or "loglik" (-loglik so that criterion to be minimized)
    #' "loglik" is the default criterion
    #' @return a [`NB_fixed_Q_zi`] object
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
      plot(seq_along(y), neg * y, type=type, log=log, ylab = criterion,
           xlab = "nb_blocks", xaxt = "n")
      axis(1, at = seq_along(y), labels = seq_along(y))
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
    criteria = function() purrr::map(self$models, "criteria") %>% purrr::reduce(rbind)
  )

)
