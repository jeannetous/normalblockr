## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS zi_normal #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for zero-inflated normal model
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood threshold under which optimization stops
#' @export
zi_normal <- R6::R6Class(
  classname = "zi_normal",

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

    #' @description Create a new [`zi_normal`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    #' @return A new [`zi_normal`] object
    initialize = function(Y, X, niter = 50, threshold = 1e-4) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y and X must be matrices.")
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
    },

    #' @description
    #' Update a [`zi_normal`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of inverse variance matrix
    #' @param kappa vector of zero-inflation probabilities
    #' @param rho posterior probabilities of zero-inflation
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`zi_normal`] object
    update = function(B = NA, dm1 = NA, kappa = NA, rho = NA, ll_list=NA) {
      if (!anyNA(B))          private$B       <- B
      if (!anyNA(dm1))        private$dm1     <- dm1
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
                  "kappa" = private$kappa, "n" = private$n,
                  "p" = private$p, "d" = private$d))
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
    kappa   = NA,   # vector of zero-inflation probabilities
    rho     = NA,   # posterior probabilities of zero-inflation
    ll_list = NA,   # list of log-likelihood values during optimization

    zi_normal_loglik  = function(Y, X, B, dm1, rho, kappa) {
      ## problem dimensions
      n   <- nrow(Y); p <- ncol(Y); d <- ncol(X)

      ## useful matrices
      Dm1 <- diag(dm1)
      indicator <- Y == 0

      J <- sum(indicator * rho)
      J <- J - .5 * sum((1 - rho) * ((Y - X %*%B)^2 %*% Dm1))
      J <- J - .5 * sum((1 - rho) %*% (-log(diag(Dm1)) + log(2 * pi)))
      J <- J + sum((indicator * rho) %*% log(kappa))
      J <- J + sum((1 - indicator * rho) %*% log(1 - kappa))
      J
    },


    zi_normal_gradB   = function(Y, X, B, dm1, rho) {
      Dm1    <- diag(dm1)
      grad <- t(X)  %*% ((1 - rho) * (Y - X %*% B)) %*% Dm1
      grad
    },

    zi_normal_optim_B = function(B0, Y, X, dm1,  rho, kappa) {
      B0_vec <- as.vector(B0)
      res <- nloptr::nloptr(
        x0 = B0_vec,
        eval_f = function(B_vec) {
          B <- matrix(B_vec, nrow = ncol(X), ncol = ncol(Y))
          Dm1 <- diag(dm1)
          YmXB <- Y - X %*%B
          list("objective" = .5 * sum((1 - rho) * (YmXB^2 %*% Dm1)),
               "gradient"  = -crossprod(X, (1 - rho) * YmXB) %*% Dm1)
        },
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          xtol_rel = 1e-6,
          maxeval = 1000
        )
      )
      newB <- matrix(res$solution, nrow = ncol(X), ncol = ncol(Y))
      newB
    },

    zi_normal_EM = function(Y, X, niter, threshold) {
      ## problem dimensions
      n <- nrow(Y); p <- ncol(Y); d <- ncol(X)

      ## useful matrices
      indicator <- 1*(Y == 0)

      ## Initialization
      rho   <- check_one_boundary(check_zero_boundary(indicator))
      kappa <- check_one_boundary(check_zero_boundary(colMeans(indicator)))
      B     <- solve(crossprod(X, X)) %*% t(X) %*% Y
      R     <- t(Y - X %*% B)
      Ddiag <- check_one_boundary(check_zero_boundary(diag(cov(t(R)))))
      D     <- diag(Ddiag)
      dm1   <- 1 / Ddiag
      Dm1   <- diag(dm1)

      ll_list <- private$zi_normal_loglik(Y, X, B, dm1, rho, kappa)

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

        loglik <- private$zi_normal_loglik(Y, X, B, dm1, rho, kappa)
        ll_list <- c(ll_list, loglik)

        if (abs(ll_list[h] - ll_list[h - 1]) < threshold)
          break
      }
      list(B = B, dm1 = dm1, kappa = kappa, rho = rho, ll_list = ll_list)
    }
  )
)
