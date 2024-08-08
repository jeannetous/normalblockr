## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood threshold under which optimization stops
#' @export
NB <- R6::R6Class(
  classname = "NB",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses
    Y  = NULL,
    #' @field X the matrix of covariates
    X = NULL,
    #' @field niter number of iterations in model optimization
    niter = NULL,
    #' @field threshold loglikelihood threshold under which optimization stops
    threshold = NULL,

    #' @description Create a new [`NB`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    #' @return A new [`nb_fixed`] object
    initialize = function(Y, X, niter = 50, threshold = 1e-4) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y, X and C must be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y <- Y
      self$X <- X
      self$niter <- niter
      self$threshold <- threshold
      private$n <- nrow(Y)
      private$p <- ncol(Y)
      private$d <- ncol(X)
      private$XtXm1   <- solve(crossprod(X, X))
      private$ll_list <- 0
    },

    #' @description
    #' Update a [`NB`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param omegaQ groups inverse variance matrix
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`zi_normal`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, ll_list=NA) {
      if (!anyNA(B))          private$B       <- B
      if (!anyNA(dm1))        private$dm1     <- dm1
      if (!anyNA(omegaQ))     private$omegaQ  <- omegaQ
      if (!anyNA(ll_list))    private$ll_list <- ll_list
    },

    #' @description calls EM optimization and updates relevant fields
    #' @return optimizes the model and updates its parameters
    optimize = function() {
      optim_out <- do.call(private$EM_optimize, list(Y = self$Y, X = self$X,
                                                     C = self$C,
                                                     niter = self$niter,
                                                     threshold = self$threshold))
      do.call(self$update, optim_out)
    },
    #' @description returns the model parameters B, dm1 and kappa
    #' @return A list containing the model parameters B, dm1, kappa
    get_model_parameters = function() {
      list("B" = private$B, "dm1" = private$dm1, "omegaQ" = private$omegaQ,
           "n" = private$n, "p" = private$p, "d" = private$d, "Q" = private$Q)
    },
    #' @description returns the model variables Y, X
    #' @return A list containing the model parameters Y, X
    get_model_variables = function() {
      list(Y = self$Y, X = self$X)
    },
    #' @description plots log-likelihood values during model optimization
    plot_loglik = function(){
      plot(1:length(private$ll_list), private$ll_list)
    },

    #' @description computes nparam, number of parameters
    nparam = function() {
      private$p * private$d + private$p + .5 * private$Q * (private$Q + 1)
    },

    #' @description computes BIC of the model
    BIC = function(){
      - 2 * private$ll_list[[length(private$ll_list)]] + log(private$n) * self$nparam()
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    n       = NULL, # number of samples
    p       = NULL, # number of responses
    Q       = NULL, # number of groups
    d       = NULL, # number of covariates
    XtXm1   = NULL, # inverse of XtX, useful for EM calculations
    B       = NA,   # regression matrix
    dm1     = NA,   # diagonal vector of inverse variance matrix
    omegaQ  = NA,   # groups variance matrix
    kappa   = NA,   # vector of zero-inflation probabilities
    rho     = NA,   # posterior probabilities of zero-inflation
    ll_list = NA,   # list of log-likelihood values during optimization

    EM_optimize = function(Y, X, C, niter, threshold) {
      variables  <- self$get_model_variables()
      parameters <- do.call(private$EM_initialize, variables)
      current    <- c(variables, parameters)
      ll_list    <- do.call(private$loglik, current)
      for (h in 2:niter) {
        parameters <- do.call(private$EM_step, current)
        current    <- c(variables, parameters)
        ll_list    <- c(ll_list, do.call(private$loglik, current))
        if (abs(ll_list[h] - ll_list[h - 1]) < threshold)
          break
      }
      c(parameters, list(ll_list = ll_list))
    },

    EM_step = function(){},
    EM_initialize = function(){},
    loglik  = function() {}
  )
)
