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
    #' @field sparsity penalty on the network density
    sparsity = NULL,
    #' @field sparsity_weights distribution of sparsity on the lement of the network
    sparsity_weights = NULL,
    #' @field niter number of iterations in model optimization
    niter = NULL,
    #' @field threshold loglikelihood threshold under which optimization stops
    threshold = NULL,
    #' @field Q number of blocks
    Q = NULL,

    #' @description Create a new [`NB`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param sparsity penalty on the network density
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    #' @return A new [`nb_fixed`] object
    initialize = function(Y, X,  sparsity = 0, niter = 50, threshold = 1e-4) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y, X and C must be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y <- Y
      self$X <- X
      self$sparsity <- sparsity
      if (sparsity > 0) {
        sparsity_weights <- matrix(1, self$Q, self$Q)
        diag(sparsity_weights) <- 0
        self$sparsity_weights  <- sparsity_weights
      }
      self$niter <- niter
      self$threshold <- threshold
      private$XtXm1   <- solve(crossprod(X, X))
      private$ll_list <- 0
    },

    #' @description
    #' Update a [`NB`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param omegaQ groups inverse variance matrix
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, ll_list = NA) {
      if (!anyNA(B))          private$B       <- B
      if (!anyNA(dm1))        private$dm1     <- dm1
      if (!anyNA(omegaQ))     private$omegaQ  <- omegaQ
      if (!anyNA(ll_list))    private$ll_list <- ll_list
    },

    #' @description calls EM optimization and updates relevant fields
    #' @return optimizes the model and updates its parameters
    optimize = function() {
      optim_out <- private$EM_optimize(niter = self$niter,
                                       threshold = self$threshold)
      do.call(self$update, optim_out)
    },
    #' @description plots log-likelihood values during model optimization
    plot_loglik = function() {
      plot(seq_along(private$ll_list), private$ll_list)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    XtXm1   = NULL, # inverse of XtX, useful for EM calculations
    B       = NA,   # regression matrix
    dm1     = NA,   # diagonal vector of inverse variance matrix
    omegaQ  = NA,   # groups variance matrix
    ll_list = NA,   # list of log-likelihood values during optimization

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
    nb_param = function() as.integer(self$p * self$d + self$p + .5 * self$Q * (self$Q + 1)),
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), omegaQ (groups precision matrix))
    model_par  = function() list(B = private$B, dm1 = private$dm1, omegaQ = private$omegaQ),
    #' @field loglik (or its variational lower bound)
    loglik = function() private$ll_list[[length(private$ll_list)]],
    #' @field penalty (penalty on log-likelihood due to sparsity)
    penalty = function() - self$sparsity * sum(abs(self$sparsity_weights * private$omegaQ)),
    #' @field BIC (or its variational lower bound)
    BIC = function() - 2 * self$loglik + log(self$n) * self$nb_param,
    #' @field AIC (or its variational lower bound)
    AIC = function() - 2 * self$loglik + 2 * self$nb_param,
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() 0,
    #' @field ICL variational lower bound of the ICL
    ICL        = function() self$BIC - self$entropy,
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      data.frame(Q = self$Q, nb_param = self$nb_param, loglik = - self$loglik,
                 BIC = self$BIC,AIC = self$AIC, ICL = self$ICL)}
    )
)
