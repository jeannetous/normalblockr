## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS nb_fixed #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed groups
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood threshold under which optimization stops
#' @export
NB_fixed <- R6::R6Class(
  classname = "NB_fixed",

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

    #' @description Create a new [`zi_normal`] object.
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
    #' Update a [`nb_fixed`] object
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
      optim_out <- do.call(private$NB_fixed_EM, list(Y = self$Y, X = self$X,
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
    d       = NULL, # number of covariates
    Q       = NULL, # number of groups
    B       = NA,   # regression matrix
    dm1     = NA,   # diagonal vector of inverse variance matrix
    omegaQ  = NA,   # groups inverse variance matrix
    ll_list = NA,   # list of log-likelihood values during optimization

    NB_fixed_loglik  = function(Y, X, C, B, dm1, omegaQ, gamma, mu) {
      ## problem dimensions
      n   <- nrow(Y); p <- ncol(Y); d <- ncol(X)

      ## useful matrices
      Dm1 <- diag(dm1)
      R   <- Y - X %*% B
      J <- - .5 * n * (p + Q) * log(2 * pi) + .5 * n * sum(log(dm1))
      J <- J - .5 * sum(R %*% Dm1 %*% t(R)) + sum(R %*% Dm1 %*% C %*% t(mu))
      J <- J - .5 * n * sum(diag(t(C) %*% Dm1 %*% C %*% gamma))
      J <- J - .5 * sum(diag(mu %*% t(C) %*% Dm1 %*% C %*% t(mu)))
      J <- J + .5 * n * log(det(omegaQ))
      J <- J - .5 * n * sum(diag(omegaQ %*% gamma))
      J <- J - .5 * n * sum(diag(mu %*% omegaQ %*% t(mu)))
      J
    },


    NB_fixed_EM = function(Y, X, C, niter, threshold) {
      ## problem dimensions
      n <- nrow(Y); p <- ncol(Y); d <- ncol(X) ; Q <- ncol(C)

      ## Initialization
      B      <- solve(crossprod(X, X)) %*% t(X) %*% Y
      R      <- t(Y - X %*% B)
      Ddiag  <- check_one_boundary(check_zero_boundary(diag(cov(t(R)))))
      D      <- diag(Ddiag)
      dm1    <- 1 / Ddiag
      Dm1    <- diag(dm1)
      gamma  <- solve(t(C) %*% Dm1 %*% C)
      mu     <- t(gamma %*% t(C) %*% Dm1 %*% R)
      omegaQ <- solve(gamma + (1/n) * t(mu) %*% mu)

      ll_list <- private$NB_fixed_loglik(Y, X, C, B, dm1, omegaQ, gamma, mu)

      for (h in 2:niter) {
        ## E step
        gamma <- solve(omegaQ + t(C) %*% Dm1 %*% C)
        mu     <- t(gamma %*% t(C) %*% Dm1 %*% R)

        ## M step
        B     <- solve(crossprod(X, X)) %*% crossprod(X, (Y - mu %*% t(C)))
        Ddiag <- (1/n) * (t(Y - X%*% B) - C %*% t(mu))^2 %*% as.vector(rep(1, n)) + diag(C %*% gamma %*% t(C))
        dm1   <- as.vector(1/Ddiag)
        omegaQ <- solve(gamma + (1/n) * t(mu) %*% mu)

        loglik <- private$NB_fixed_loglik(Y, X, C, B, dm1, omegaQ, gamma, mu)
        ll_list <- c(ll_list, loglik)

        if (abs(ll_list[h] - ll_list[h - 1]) < threshold)
          break
      }
      list(B = B, dm1 = dm1, omegaQ = omegaQ,  ll_list = ll_list)
    }
  )
)
