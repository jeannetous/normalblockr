## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_blocks #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed groups
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param sparsity to add on blocks precision matrix
NB_fixed_blocks <- R6::R6Class(
  classname = "NB_fixed_blocks",
  inherit = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field C the matrix of species groups
    C = NULL,

    #' @description Create a new [`NB_fixed_blocks`] object.
    #' @param C group matrix C_jq = 1 if species j belongs to group q
    #' @return A new [`NB_fixed_blocks`] object
    initialize = function(Y, X, C, sparsity = 0) {
      if (!is.matrix(C)) stop("C must be a matrix.")
      super$initialize(Y, X, ncol(C), sparsity)
      self$C <- C
    },

    #' @description
    #' Update a [`NB_fixed_blocks`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param omegaQ groups inverse variance matrix
    #' @param gamma  variance of  posterior distribution of W
    #' @param mu  mean for posterior distribution of W
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`zi_normal`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, gamma = NA, mu = NA,
                      ll_list = NA) {
      super$update(B, dm1, omegaQ, ll_list)
      if (!anyNA(gamma)) private$gamma <- gamma
      if (!anyNA(mu))   private$mu   <- mu
    }),


  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    gamma   = NA, # variance of  posterior distribution of W
    mu      = NA, #  mean for posterior distribution of W

    compute_loglik  = function(B, dm1, omegaQ, gamma, mu) {
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)
      log_det_gamma  <- as.numeric(determinant(gamma, logarithm = TRUE)$modulus)

      J <- -.5 * self$n * self$p * log(2 * pi * exp(1))
      J <- J + .5 * self$n * sum(log(dm1)) + .5 * self$n * log_det_omegaQ
      J <- J + .5 * self$n * log_det_gamma
      if (self$sparsity > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + self$n *self$Q / 2 - .5 * sum(diag(omegaQ %*% (self$n * gamma + t(mu) %*% mu)))
        J <- J - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
      J
    },

    EM_initialize = function() {
      B      <- private$XtXm1 %*% t(self$X) %*% self$Y
      R      <- self$Y - self$X %*% B
      dm1    <- 1 / apply(R, 2, var)
      gamma  <- diag(1/colSums(dm1 * self$C))
      mu     <- R %*% (dm1 * self$C) %*% gamma
      omegaQ <- solve(gamma + (1 / self$n) * crossprod(mu))
      list(B = B, dm1 = dm1, omegaQ = omegaQ, gamma = gamma, mu = mu)
    },

    EM_step = function(B, dm1, omegaQ, gamma, mu) {

      ## E step
      gamma <- solve(omegaQ + diag(colSums(dm1 * self$C), self$Q, self$Q))
      mu    <- (self$Y - self$X %*% B) %*% (dm1 * self$C) %*% gamma

      ## M step
      muCT   <- mu %*% t(self$C)
      B      <- private$XtXm1 %*% crossprod(self$X, self$Y - muCT)
      ddiag  <- colMeans((self$Y - self$X %*% B - muCT)^2) + self$C %*% diag(gamma)
      dm1    <- 1 / as.vector(ddiag)
      sigmaQ <- gamma + (1 / self$n) * crossprod(mu)

      if (self$sparsity == 0) {
        omegaQ <- solve(sigmaQ)
      } else {
        glasso_out <- glassoFast::glassoFast(sigmaQ, rho = self$sparsity * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) stop("GLasso fails")
        omegaQ <- Matrix::symmpart(glasso_out$wi)
      }
      list(B = B, dm1 = dm1, omegaQ = omegaQ, gamma = gamma, mu = mu)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par  = function() list(gamma = private$gamma, mu = private$mu),
    #' @field cond_par a list with the matrices of the conditional latent distirbution: mu (mean), Gamma (variance)
    cond_par    = function() list(mu = private$mu,  Gamma = private$gamma),
    #' @field clustering given as a list of labels
    clustering = function() get_clusters(self$C),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      log_det_Gamma <- as.numeric(determinant(private$gamma)$modulus)
      ent <- 0.5 * self$n * self$Q * log(2 * pi* exp(1)) + .5 * self$n * log_det_Gamma
      ent
    },
    #' @field fitted Y values predicted by the model
    fitted = function() self$X %*% private$B + private$mu %*% t(self$C)

  )
)
