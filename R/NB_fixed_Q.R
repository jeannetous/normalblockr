## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed groups
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood threshold under which optimization stops
#' @export
NB_fixed_Q <- R6::R6Class(
  classname = "NB_fixed_Q",
  inherit = NB,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(


    #' @description Create a new [`NB`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param Q required number of groups
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    #' @return A new [`NB_fixed_Q`] object
    initialize = function(Y, X, Q, niter = 50, threshold = 1e-4) {
      super$initialize(Y, X, niter, threshold)
      private$Q <- Q
    },

    #' @description
    #' Update a [`NB_fixed_Q`] object
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

    #' @description returns the model parameters B, dm1 and kappa
    #' @return A list containing the model parameters B, dm1, kappa
    get_model_parameters = function() {
      parameters       <- super$get_model_parameters()
      parameters$alpha <- private$alpha
      parameters$tau   <- private$tau
      parameters$M     <- private$M
      parameters$S     <- private$S
      return(parameters)
    },
    #' @return nparam, number of parameters
    nparam = function() {
      number_parameters <- super$nparam()
      number_parameters + private$Q + private$p * private$Q + 2 * private$n * private$Q
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

    loglik  = function(Y, X, B, dm1, omegaQ, alpha, tau, M, S) {
      ## problem dimensions
      n   <- private$n; p <-  private$p; d <-  private$d; Q <-  private$Q

      R              <- t(Y - X %*% B)
      ones           <- as.vector(rep(1, n))
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      # expectation of log(p(Y | W, C))
      elbo <- - 0.5 * n * p * log(2*pi) + 0.5 * n * sum(log(dm1))
      elbo <- elbo - 0.5 * sum(dm1*((R^2 + (tau %*% t(M^2)) - 2*R* (tau %*% t(M))) %*% ones + n * (tau %*% S)))

      # expectation of log(p(W))
      elbo <- elbo - 0.5 * n * Q * log(2*pi) + 0.5 * n * log_det_omegaQ
      elbo <- elbo - 0.5 * sum(crossprod(ones, (M %*% omegaQ) * M))
      elbo <- elbo - 0.5 * n * S %*% diag(omegaQ)

      # expectation of log(p(C))
      elbo <- elbo + sum(crossprod(log(alpha), t(tau)))

      # Entropy term for W
      elbo <- elbo + 0.5 * n * Q * log(2*pi* exp(1)) + .5 * n * sum(log(S))

      # Entropy term for C
      elbo <- elbo - sum(xlogx(tau))

      elbo
    },

    EM_initialize = function(Y, X) {
      n   <- private$n; p <-  private$p; d <-  private$d; Q <-  private$Q
      B                 <- private$XtXm1%*% t(X) %*% Y
      R                 <- t(Y - X %*% B)
      clustering_kmeans <- kmeans(R, Q, nstart=30)
      cl                <- clustering_kmeans$cluster
      tau               <- as_indicator(cl)
      alpha             <- colMeans(tau)
      S                 <- rep(0.1, Q)
      M                 <- matrix(rep(0, n * Q), nrow=n)
      dm1               <- as.vector(rep(1, p))
      omegaQ            <- diag(rep(1, Q))
      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, tau = tau, M = M, S = S)
    },

    EM_step = function(Y, X, B, dm1, omegaQ, alpha, tau, M, S) {
      n    <- private$n ; p <- private$p
      R    <- t(Y - X %*% B)
      G    <- solve(diag(colSums(as.vector(dm1) * tau)) + omegaQ)
      ones <- as.vector(rep(1, n))

      # E step
      M         <- crossprod(as.vector(dm1) * R, tau) %*% G
      S         <- as.vector(1/(as.vector(dm1) %*% tau + t(diag(omegaQ))))
      pre_tau   <- -.5 * dm1 %*% t(ones) %*% M^2 -.5 * dm1 %*% t(n*S)
      pre_tau   <- pre_tau + dm1 * (R %*% M)  + outer(rep(1,p), log(alpha)) - 1
      tau       <- t(check_zero_boundary(check_one_boundary(apply(pre_tau, 1, softmax))))

      # M step
      omegaQ <- n * solve((t(M) %*% M) + n * diag(S))
      B <- private$XtXm1 %*% t(X) %*% (Y- M %*% t(tau))
      dm1 <- as.vector(n/((R^2 - 2*R*(tau %*% t(M)) + tau %*% t(M^2) +  tau %*% t(ones %*% t(S))) %*% ones))
      alpha <- colMeans(tau)

      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, tau = tau, M = M, S = S)
    }
  )
)
