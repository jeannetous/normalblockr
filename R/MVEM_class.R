## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS MVEM    #######################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 generic class for Multivariate Model inferred with EM
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @export
MVEM <- R6::R6Class(
  classname = "MVEM",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses
    Y  = NULL,
    #' @field X the matrix of covariates
    X = NULL,

    #' @description Create a new [`NB`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param penalty penalty on the network density
    #' @return A new [`nb_fixed`] object
    initialize = function(Y, X,  penalty = 0) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y, X and C must be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y <- Y
      self$X <- X
      private$XtXm1   <- solve(crossprod(X, X))
      private$ll_list <- 0
    },

    #' @description
    #' Update a [`NB`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB`] object
    update = function(B = NA, dm1 = NA, ll_list = NA) {
      if (!anyNA(B))          private$B       <- B
      if (!anyNA(dm1))        private$dm1     <- dm1
      if (!anyNA(ll_list))    private$ll_list <- ll_list
    },

    #' @description calls EM optimization and updates relevant fields
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    #' @return optimizes the model and updates its parameters
    optimize = function(niter = 100, threshold = 1e-4) {
      optim_out <- private$EM_optimize(niter, threshold)
      do.call(self$update, optim_out)
    },

    #' @param type char for line type (see plot.default)
    #' @param log char for logarithmic axes (see plot.default)
    #' @param neg boolean plot negative log-likelihood (useful when log="y")
    #' @description plots log-likelihood values during model optimization
    plot_loglik = function(type = "b", log = "", neg = FALSE) {
      neg <- ifelse(neg, -1, 1)
      plot(seq_along(self$objective), neg * self$objective, type = type, log = log)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    XtXm1     = NA, # inverse of XtX, useful for EM calculations
    B         = NA, # regression matrix
    dm1       = NA, # diagonal vector of inverse variance matrix
    ll_list   = NA, # list of log-likelihood values during optimization

    EM_optimize = function(niter, threshold) {
      parameters <- do.call(private$EM_initialize, list())
      ll_list    <- do.call(private$compute_loglik, parameters)
      for (h in 2:niter) {
        parameters <- do.call(private$EM_step, parameters)
        ll_list    <- c(ll_list, do.call(private$compute_loglik, parameters))
        if (abs(ll_list[h] - ll_list[h - 1]) < threshold)
          break
      }
      c(parameters, list(ll_list = ll_list))
    },

    EM_step = function() {},
    EM_initialize = function() {},
    compute_loglik  = function() {}
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples
    n = function() nrow(self$Y),
    #' @field p number of responses per sample
    p = function() ncol(self$Y),
    #' @field d number of variables (dimensions in X)
    d = function() ncol(self$X),
    #' @field nb_param number of parameters in the model
    nb_param = function() as.integer(self$p * self$d + self$p),
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance)
    model_par  = function() list(B = private$B, dm1 = private$dm1),
    #' @field loglik (or its variational lower bound)
    loglik = function() private$ll_list[[length(private$ll_list)]] + self$penalty_term,
    #' @field penalty_term (for cases when a penalty is placed on the precision matrix)
    penalty_term = function() 0,
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() 0,
    #' @field deviance (or its variational lower bound)
    deviance = function() - 2 * self$loglik,
    #' @field BIC (or its variational lower bound)
    BIC = function() self$deviance + log(self$n) * self$nb_param,
    #' @field AIC (or its variational lower bound)
    AIC = function() self$deviance + 2 * self$nb_param,
    #' @field ICL variational lower bound of the ICL
    ICL        = function() self$BIC + 2 * self$entropy,
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      data.frame(nb_param = self$nb_param, loglik = self$loglik,
                 deviance = self$deviance, BIC = self$BIC,AIC = self$AIC, ICL = self$ICL)
    },
    #' @field objective evolution of the objective function during (V)EM algorithm
    objective = function() private$ll_list[-1]
  )
)
