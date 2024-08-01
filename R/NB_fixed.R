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
    #' @field thresold loglikelihood threshold under which optimization stops
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
      self$X <- C
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
    update = function(B = NA, dm1 = NA, omegaQ = NA, kappa = NA, rho = NA,
                      ll_list=NA) {
      if (!anyNA(B))          private$B       <- B
      if (!anyNA(dm1))        private$dm1     <- dm1
      if (!anyNA(omegaQ))     private$omegaQ     <- omegaQ
      if (!anyNA(kappa))      private$kappa   <- kappa
      if (!anyNA(rho))        private$rho     <- rho
      if (!anyNA(ll_list))    private$ll_list <- ll_list
    },

    #' @description calls EM optimization and updates relevant fields
    #' @return optimizes the model and updates its parameters
    optimize = function() {
      optim_out <- do.call(private$zi_normal_EM, list(Y = self$Y, X = self$X,
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
    B       = NA,   # regression matrix
    dm1     = NA,   # diagonal vector of inverse variance matrix
    omegaQ  = NA,   # groups inverse variance matrix
    ll_list = NA,   # list of log-likelihood values during optimization

    NB_fixed_loglik  = function(Y, X, C, B, dm1, omegaQ, gamma, mu) {
      ## problem dimensions
      n   <- nrow(Y); p <- ncol(Y); d <- ncol(X)

      ## useful matrices
      Dm1 <- diag(dm1)
      indicator <- matrix(sapply(Y, function(x) ifelse(x == 0, 1, 0)),
                          nrow = nrow(Y), ncol = ncol(Y))

      J <- sum(indicator * rho)
      J <- J - .5 * sum((1 - rho) * ((Y - X %*%B )^2 %*% Dm1))
      J <- J - .5 * sum((1 - rho) %*% (-log(diag(Dm1)) + log(2 * pi)))
      J <- J + sum((indicator * rho) %*% log(kappa))
      J <- J + sum((1 - indicator * rho) %*% log(1 - kappa))
      J
    },


    NB_fixed_EM = function(Y, X, C, niter, threshold) {
      ## problem dimensions
      n <- nrow(Y); p <- ncol(Y); d <- ncol(X)

      ## useful matrices
      indicator <- matrix(sapply(Y, function(x) ifelse(x == 0, 1, 0)),
                          nrow = nrow(Y), ncol = ncol(Y))
      ## Initialization
      rho   <- check_one_boundary(check_zero_boundary(indicator))
      kappa <- check_one_boundary(check_zero_boundary(colMeans(indicator)))
      B     <- solve(crossprod(X, X)) %*% t(X) %*% Y
      R     <- t(Y - X %*% B)
      Ddiag <- check_one_boundary(check_zero_boundary(diag(cov(t(R)))))
      D     <- diag(Ddiag)
      dm1   <- 1 / Ddiag
      Dm1   <- diag(dm1)

      ll_list <- private$NB_fixed_loglik(Y, X, C, B, dm1, omegaQ, gamma, mu)

      for (h in 2:niter) {
        ## E step
        kfactor <- diag((1 - kappa) / kappa)
        gamma   <- 1 / sqrt(2 * pi) * sqrt(Dm1) %*% exp(- .5 * Dm1 %*%
                                                          crossprod(B, t(X))^2)
        rho     <- check_one_boundary(check_zero_boundary(indicator * t((1 / (1 + kfactor %*% gamma)))))

        ## M step
        R     <- t(Y - X %*% B)
        Ddiag <- check_zero_boundary(rowSums(t(1 - rho) * R^2)  / rowSums(1 - t(rho)))
        dm1   <- 1 / Ddiag
        Dm1   <- diag(dm1)
        kappa <- colMeans(rho)
        B     <- private$zi_normal_optim_B(B, Y, X, dm1, rho, kappa)

        loglik <- private$NB_fixed_loglik(Y, X, C, B, dm1, omegaQ, gamma, mu)
        ll_list <- c(ll_list, loglik)

        if (abs(ll_list[h] - ll_list[h - 1]) < threshold)
          break
      }
      list(B = B, dm1 = dm1, omegaQ = omegaQ, kappa = kappa, rho = rho,
           ll_list = ll_list)
    }
  )
)
