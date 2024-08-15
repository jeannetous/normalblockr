## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_blocks #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed groups
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param sparsity to add on blocks precision matrix
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood threshold under which optimization stops
#' @export
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
    initialize = function(Y, X, C, sparsity = 0, niter = 50, threshold = 1e-4) {
      if (!is.matrix(C)) stop("C must be a matrix.")
      self$C <- C
      self$Q <- ncol(C)
      super$initialize(Y, X, sparsity,niter, threshold)
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
                      ll_list=NA) {
      super$update(B, dm1, omegaQ, ll_list)
      if (!anyNA(gamma)) private$gamma <- gamma
      if (!anyNA(mu))   private$mu   <- mu
    }),


  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    gamma   = NA, # variance of  posterior distribution of W
    mu      = NA,  #  mean for posterior distribution of W


    compute_loglik  = function(B, dm1, omegaQ, gamma, mu) {
      R   <- self$Y - self$X %*% B
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      J <- - .5 * self$n * (self$p + self$Q) * log(2 * pi) + .5 * self$n * sum(log(dm1))
      J <- J - .5 * sum(R %*% (dm1 * t(R))) + sum(R %*%  (dm1 * self$C) %*% t(mu))
      J <- J - .5 * self$n * sum(diag(t(self$C) %*% (dm1 * self$C) %*% gamma))
      J <- J - .5 * sum(diag(mu %*% t(self$C) %*%  (dm1 * self$C) %*% t(mu)))
      J <- J + .5 * self$n * log_det_omegaQ
      J <- J - .5 * self$n * sum(diag(omegaQ %*% gamma))
      J <- J - .5 * sum(diag(mu %*% omegaQ %*% t(mu)))
      if(self$sparsity == 0 ) {J
      }else{
        J - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
    },

    EM_initialize = function() {
      B      <- private$XtXm1 %*% t(self$X) %*% self$Y
      R      <- t(self$Y - self$X %*% B)
      dm1    <- 1 / check_one_boundary(check_zero_boundary(diag(cov(t(R)))))
      gamma  <- solve(t(self$C) %*% (dm1 * self$C))
      mu     <- t(gamma %*% t(self$C) %*% (dm1 * R))
      omegaQ <- solve(gamma + (1 / self$n) * t(mu) %*% mu)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, gamma = gamma, mu = mu)
    },

    EM_step = function(B, dm1, omegaQ, gamma, mu) {
      R      <- t(self$Y - self$X %*% B)

      ## E step
      gamma  <- solve(omegaQ + t(self$C) %*% (dm1 * self$C))
      mu     <- t(gamma %*% t(self$C) %*% (dm1 * R))

      ## M step
      B      <- private$XtXm1 %*% t(self$X) %*% (self$Y - mu %*% t(self$C))
      ddiag  <- (1 / self$n) * (diag(R %*% t(R)) - 2 * diag(R %*% mu %*% t(self$C)))
      ddiag  <- ddiag + (1 / self$n) * diag(self$C %*% t(mu) %*% mu %*% t(self$C))
      ddiag  <- ddiag + diag(self$C %*% gamma %*% t(self$C))
      dm1    <- as.vector(1 / ddiag)
      if (self$sparsity == 0 ) {
        omegaQ <- solve(gamma + (1 / self$n) * t(mu) %*% mu)
      }else {
        sigma_hat <- gamma + (1 / self$n) * t(mu) %*% mu
        glasso_out <- glassoFast::glassoFast(sigma_hat, rho = self$sparsity * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) break
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
    posterior_par  = function() list(gamma = private$gamma, mu = private$mu)
  )
)
