## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_unknown #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed groups
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood threshold under which optimization stops
#' @export
NB_unknown <- R6::R6Class(
  classname = "NB_unknown",
  inherit = NB,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description
    #' Update a [`NB_unknown`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param omegaQ groups inverse variance matrix
    #' @param alpha vector of groups probabilities
    #' @param tau posterior probabilities for group affectation
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`zi_normal`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, alpha = NA, tau = NA,
                      M = NA, S = NA, ll_list=NA) {
      super$update(B, dm1, omegaQ, ll_list)
      if (!anyNA(alpha)) private$alpha <- alpha
      if (!anyNA(tau))   private$tau   <- tau
      if (!anyNA(M))     private$M     <- M
      if (!anyNA(S))     private$M     <- S
    },


    #' @description calls EM optimization and updates relevant fields
    #' @return optimizes the model and updates its parameters
    optimize = function() {
      optim_out <- do.call(private$NB_unknown_EM, list(Y = self$Y, X = self$X,
                                                     C = self$C,
                                                     niter = self$niter,
                                                     threshold = self$threshold))
      do.call(self$update, optim_out)
    },

    #' @description returns the model parameters B, dm1 and kappa
    #' @return A list containing the model parameters B, dm1, kappa
    get_model_parameters = function() {
      parameters       <- super$get_model_parameters()
      parameters$alpha <- private$alpha
      parameters$tau   <- private$tau
      parameters$M     <- private$M
      parameters$S     <- private$S
      return(parameters)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    alpha   = NULL, # vector of groups probabilities
    tau     = NULL, # posterior probabilities for group affectation
    M       = NULL, # variational mean for posterior distribution of W
    S       = NULL, # variational diagonal of variances for posterior distribution of W

    NB_unknown_loglik  = function(Y, X, C, B, dm1, omegaQ, gamma, mu) {
      ## problem dimensions
      n   <- nrow(Y); p <- ncol(Y); d <- ncol(X); Q <- ncol(C)

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


    NB_unknown_EM = function(Y, X, C, niter, threshold) {
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

      ll_list <- private$NB_unknown_loglik(Y, X, C, B, dm1, omegaQ, gamma, mu)

      for (h in 2:niter) {
        ## E step
        gamma <- solve(omegaQ + t(C) %*% Dm1 %*% C)
        mu     <- t(gamma %*% t(C) %*% Dm1 %*% R)

        ## M step
        B     <- solve(crossprod(X, X)) %*% crossprod(X, (Y - mu %*% t(C)))
        Ddiag <- (1/n) * (t(Y - X%*% B) - C %*% t(mu))^2 %*% as.vector(rep(1, n)) + diag(C %*% gamma %*% t(C))
        dm1   <- as.vector(1/Ddiag)
        omegaQ <- solve(gamma + (1/n) * t(mu) %*% mu)

        loglik <- private$NB_unknown_loglik(Y, X, C, B, dm1, omegaQ, gamma, mu)
        ll_list <- c(ll_list, loglik)

        if (abs(ll_list[h] - ll_list[h - 1]) < threshold)
          break
      }
      list(B = B, dm1 = dm1, omegaQ = omegaQ,  ll_list = ll_list)
    }
  )
)
