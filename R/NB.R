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
    #' @field C the matrix of species groups
    C = NULL,
    #' @field niter number of iterations in model optimization
    niter = NULL,
    #' @field threshold loglikelihood threshold under which optimization stops
    threshold = NULL,

    #' @description Create a new [`NB`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param C group matrix C_jq = 1 if species j belongs to group q
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    #' @return A new [`nb_fixed`] object
    initialize = function(Y, X, C, niter = 50, threshold = 1e-4) {
      if (!is.matrix(Y) || !is.matrix(X) || !is.matrix(C)) {
        stop("Y, X and C must be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y <- Y
      self$X <- X
      self$C <- C
      self$niter <- niter
      self$threshold <- threshold
      private$n <- nrow(Y)
      private$p <- ncol(Y)
      private$d <- ncol(X)
      private$Q <- ncol(C)
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
      return(list("B" = private$B, "dm1" = private$dm1,
                  "omegaQ" = private$omegaQ, "n" = private$n, "p" = private$p,
                  "d" = private$d, "Q" = private$Q))
    },
    #' @description plots log-likelihood values during model optimization
    plot_loglik = function(){
      plot(1:length(private$ll_list), private$ll_list)
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
    B       = NA,   # regression matrix
    dm1     = NA,   # diagonal vector of inverse variance matrix
    kappa   = NA,   # vector of zero-inflation probabilities
    rho     = NA,   # posterior probabilities of zero-inflation
    ll_list = NA,   # list of log-likelihood values during optimization

    EM_optimize = function(Y, X, C, niter, threshold) {
      n <- nrow(Y); p <- ncol(Y); d <- ncol(X) ; Q <- ncol(C)
      current_parameters <- private$EM_initialize(Y, X, C)
      ll_list            <- do.call(private$loglik, current_parameters)
      for (h in 2:niter) {
        current_parameters <- do.call(private$EM_step, current_parameters)
        ll_list            <- do.call(private$loglik, current_parameters)
        if (abs(ll_list[h] - ll_list[h - 1]) < threshold)
          break
      }
      current_parameters
    },
    EM_step = function(){},
    EM_initialize = function(){},
    loglik  = function() {}
  )
)
