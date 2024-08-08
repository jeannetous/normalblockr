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
    #' @field zeros indicator matrix of zeros in Y
    zeros = NULL,
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
      self$zeros <- 1 * (Y == 0)
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
      optim_out <- private$zi_normal_EM(niter = self$niter, threshold = self$threshold)
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
      plot(1:length(private$ll_list), private$ll_list, type = 'b')
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

    zi_normal_loglik  = function(B, dm1, rho, kappa) {

      J <- sum(self$zeros * rho)
      J <- J - .5 * sum((1 - rho) * t(t(self$Y - self$X %*% B)^2 * dm1))
      J <- J - .5 * sum((1 - rho) %*% (-log(dm1) + log(2 * pi)))
      J <- J + sum((self$zeros * rho) %*% log(kappa))
      J <- J + sum((1 - self$zeros * rho) %*% log(1 - kappa))
      J
    },

    zi_normal_obj_grad_B = function(B_vec, dm1, rho) {
      YmXB <- self$Y - self$X %*% matrix(B_vec, nrow = private$d, ncol = private$p)

      ## A %*% diag(d) = t(t(A) * d)
      grad <- -t(t(crossprod(self$X, (1 - rho) * YmXB)) * dm1)
      obj <- .5 * sum((1 - rho) * t(t(YmXB^2) * dm1))

      res <- list("objective" = obj, "gradient"  = grad)
      res
    },

    zi_normal_optim_B = function(B0, dm1,  rho, kappa) {
      B0_vec <- as.vector(B0)
      res <- nloptr::nloptr(
        x0 = B0_vec,
        eval_f = private$zi_normal_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        dm1 = dm1,
        rho = rho
      )
      newB <- matrix(res$solution, nrow = private$d, ncol = private$p)
      newB
    },

    zi_normal_EM = function(niter, threshold) {

      ## Initialization
      rho   <- check_one_boundary(check_zero_boundary(self$zeros))
      kappa <- check_one_boundary(check_zero_boundary(colMeans(self$zeros)))
      B     <- solve(crossprod(self$X)) %*% crossprod(self$X, self$Y)
      dm1   <- 1 / check_one_boundary(check_zero_boundary(diag(cov(self$Y - self$X %*% B))))

      ll_list <- private$zi_normal_loglik(B, dm1, rho, kappa)

      for (h in 2:niter) {

        ## E step
        kfactor <- (1 - kappa) / kappa
        gamma   <- 1 / sqrt(2 * pi) * sqrt(dm1) * exp(- .5 * dm1 * t(self$X %*% B)^2)
        rho     <- check_one_boundary(check_zero_boundary(self$zeros * t((1 / (1 + kfactor * gamma)))))

        ## M step
        R     <- self$Y - self$X %*% B
        dm1   <- 1 / check_zero_boundary(colSums((1 - rho) * R^2)  / colSums(1 - rho))
        kappa <- colMeans(rho)
        B     <- private$zi_normal_optim_B(B, dm1, rho, kappa)

        ## Assessing convergence
        loglik <- private$zi_normal_loglik(B, dm1, rho, kappa)
        ll_list <- c(ll_list, loglik)
        if (abs(ll_list[h] - ll_list[h - 1]) < threshold)
          break

      }
      list(B = B, dm1 = dm1, kappa = kappa, rho = rho, ll_list = ll_list)
    }
  )
)
