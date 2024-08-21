## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed Q (number of groups)
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param Q number of blocks
#' @param sparsity to add on blocks precision matrix
#' @export
NB_fixed_Q <- R6::R6Class(
  classname = "NB_fixed_Q",
  inherit = NB,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NB_fixed_Q`] object.
    #' @param Q required number of groups
    #' @return A new [`NB_fixed_Q`] object
    initialize = function(Y, X, Q, sparsity = 0) {
      self$Q <- Q
      super$initialize(Y, X, sparsity)
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
    #' @return Update the current [`NB_fixed_Q`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, alpha = NA, tau = NA,
                      M = NA, S = NA, ll_list = NA) {
      super$update(B, dm1, omegaQ, ll_list)
      if (!anyNA(alpha)) private$alpha <- alpha
      if (!anyNA(tau))   private$tau   <- tau
      if (!anyNA(M))     private$M     <- M
      if (!anyNA(S))     private$S     <- S
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    alpha   = NA, # vector of groups probabilities
    M       = NA, # variational mean for posterior distribution of W
    S       = NA, # variational diagonal of variances for posterior distribution of W
    tau     = NA, # posterior probabilities for group affectation

    compute_loglik  = function(B, dm1, omegaQ, alpha, tau, M, S) {
      ## problem dimensions
      n   <- self$n; p <-  self$p; d <-  self$d; Q <-  self$Q

      R              <- t(self$Y - self$X %*% B)
      ones           <- as.vector(rep(1, n))
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      # expectation of log(p(Y | W, C))
      elbo <- - 0.5 * n * p * log(2 * pi) + 0.5 * n * sum(log(dm1))
      elbo <- elbo - 0.5 * sum(dm1 * ((R^2 + (tau %*% t(M^2)) - 2 * R * (tau %*% t(M))) %*% ones + n * (tau %*% S)))

      # expectation of log(p(W))
      elbo <- elbo - 0.5 * n * Q * log(2 * pi) + 0.5 * n * log_det_omegaQ
      elbo <- elbo - 0.5 * sum(crossprod(ones, (M %*% omegaQ) * M))
      elbo <- elbo - 0.5 * n * S %*% diag(omegaQ)

      # expectation of log(p(C))
      elbo <- elbo + sum(crossprod(log(alpha), t(tau)))

      # Entropy term for W
      elbo <- elbo + 0.5 * n * Q * log(2 * pi * exp(1)) + .5 * n * sum(log(S))

      # Entropy term for C
      elbo <- elbo - sum(xlogx(tau))

      if (self$sparsity == 0 ) {
        elbo
      }else {
        elbo - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
    },

    EM_initialize = function() {
      n   <- self$n; p <-  self$p; d <-  self$d; Q <-  self$Q
      B                 <- private$XtXm1 %*% t(self$X) %*% self$Y
      R                 <- t(self$Y - self$X %*% B)
      clustering_kmeans <- kmeans(R, Q, nstart = 30)
      cl                <- clustering_kmeans$cluster
      tau               <- as_indicator(cl)
      alpha             <- colMeans(tau)
      S                 <- rep(0.1, Q)
      M                 <- matrix(rep(0, n * Q), nrow = n)
      dm1               <- as.vector(rep(1, p))
      omegaQ            <- diag(rep(1, Q))
      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, tau = tau, M = M, S = S)
    },

    EM_step = function(B, dm1, omegaQ, alpha, tau, M, S) {
      n    <- self$n ; p <- self$p
      R    <- t(self$Y - self$X %*% B)
      G    <- solve(diag(colSums(as.vector(dm1) * tau)) + omegaQ)
      ones <- as.vector(rep(1, self$n))

      # E step
      M         <- crossprod(as.vector(dm1) * R, tau) %*% G
      S         <- as.vector(1 / (as.vector(dm1) %*% tau + t(diag(omegaQ))))
      eta       <- -.5 * dm1 %*% t(ones) %*% M^2 - .5 * dm1 %*% t(self$n * S)
      eta       <- eta + dm1 * (R %*% M)  + outer(rep(1, self$p), log(alpha)) - 1
      tau       <- t(check_zero_boundary(check_one_boundary(apply(eta, 1, softmax))))

      # M step
      if (self$sparsity == 0) {
        omegaQ <- self$n * solve((t(M) %*% M) + self$n * diag(S))
      }else {
        sigma_hat <- (1 / self$n) * (t(M) %*% M + self$n * diag(S))
        glasso_out <- glassoFast::glassoFast(sigma_hat, rho = self$sparsity * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) stop("GLasso fails")
        omegaQ <- Matrix::symmpart(glasso_out$wi)
      }
      B <- private$XtXm1 %*% t(self$X) %*% (self$Y - M %*% t(tau))
      dm1 <- as.vector(self$n / ((R^2 - 2 * R * (tau %*% t(M)) + tau %*% t(M^2) +  tau %*% t(ones %*% t(S))) %*% ones))
      alpha <- colMeans(tau)

      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, tau = tau, M = M, S = S)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function() super$nb_param + self$Q + self$p * self$Q + 2 * self$n * self$Q,
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), omegaQ (groups precision matrix))
    model_par  = function() {
      parameters       <- super$model_par
      parameters$alpha <- private$alpha
      parameters},
    #' @field var_par a list with the matrices of the variational parameters: M (means), S (variances), tau (posterior group probabilities)
    var_par    = function() list(M = private$M,  S = private$S, tau = private$tau),
    #' @field clustering a list of labels giving the clustering obtained in the model
    clustering = function() get_clusters(private$tau),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      ent <- 0.5 * self$n * self$Q * log(2 * pi* exp(1)) + .5 * self$n * sum(log(private$S))
      ent <- ent - sum(xlogx(private$tau))
      return(ent)
    }
    ),
)
