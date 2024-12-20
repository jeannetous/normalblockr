## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_blocks ##############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for normal-block model with fixed groups
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param penalty to add on blocks precision matrix for sparsity
#' @param control structured list of more specific parameters, to generate with NB_param()
NB_fixed_blocks <- R6::R6Class(
  classname = "NB_fixed_blocks",
  inherit = NB,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    C       = NA, # the matrix of species groups
    gamma   = NA, # variance of  posterior distribution of W
    mu      = NA, #  mean for posterior distribution of W

    compute_loglik  = function(B, dm1, omegaQ, gamma, mu) {
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)
      log_det_gamma  <- as.numeric(determinant(gamma, logarithm = TRUE)$modulus)

      J <- -.5 * self$n * self$p * log(2 * pi * exp(1))
      J <- J + .5 * self$n * sum(log(dm1)) + .5 * self$n * log_det_omegaQ
      J <- J + .5 * self$n * log_det_gamma
      if (self$penalty > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + self$n *self$Q / 2 - .5 * sum(diag(omegaQ %*% (self$n * gamma + t(mu) %*% mu)))
        J <- J - self$penalty * sum(abs(self$sparsity_weights * omegaQ))
      }
      J
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NB_fixed_blocks`] object.
    #' @param C group matrix C_jq = 1 if species j belongs to group q
    #' @return A new [`NB_fixed_blocks`] object
    initialize = function(Y, X, C, penalty = 0, control = NB_param()) {
      if (!is.matrix(C)) stop("C must be a matrix.")
      super$initialize(Y, X, ncol(C), penalty, control = control)
      private$C     <- C
      private$mu    <- matrix(0, self$n, self$Q)
      private$gamma <- diag(1, self$Q, self$Q)
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
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par  = function() list(gamma = private$gamma, mu = private$mu),
    #' @field clustering given as a list of labels
    clustering = function() get_clusters(private$C),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      log_det_Gamma <- as.numeric(determinant(private$gamma)$modulus)
      ent <- .5 * self$n * self$Q * log(2 * pi* exp(1)) + .5 * self$n * log_det_Gamma
      ent
    },
    #' @field fitted Y values predicted by the model
    fitted = function() self$X %*% private$B + private$mu %*% t(private$C)
  )
)

#' R6 class for normal-block model with fixed groups and diagonal residual covariance
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param penalty to add on blocks precision matrix for sparsity
NB_fixed_blocks_diagonal <- R6::R6Class(
  classname = "NB_fixed_blocks_diagonal",
  inherit = NB_fixed_blocks,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    EM_initialize = function() {
      if(any(sapply(self$model_par, function(x) any(is.na(x))))){
        B      <- private$XtXm1 %*% t(self$X) %*% self$Y
        dm1    <- as.vector(1/colMeans((self$Y - self$X %*% B)^2))
        omegaQ <- diag(colSums(dm1 * private$C), self$Q, self$Q)
        list(B = B, dm1 = dm1, omegaQ = omegaQ, gamma = private$gamma, mu = private$mu)
      }else{
        list(B = private$B, dm1 = private$dm1, omegaQ = private$omegaQ,
             gamma = private$gamma, mu = private$mu)
      }
    },

    EM_step = function(B, dm1, omegaQ, gamma, mu) {

      ## E step
      gamma <- solve(omegaQ + diag(colSums(dm1 * private$C), self$Q, self$Q))
      mu    <- (self$Y - self$X %*% B) %*% (dm1 * private$C) %*% gamma

      ## M step
      YmmuCT <- self$Y - mu %*% t(private$C)
      B      <- private$XtXm1 %*% crossprod(self$X, YmmuCT)
      ddiag  <- colMeans((YmmuCT - self$X %*% B)^2) + private$C %*% diag(gamma)
      dm1    <- 1 / as.vector(ddiag)
      omegaQ <- private$get_omegaQ(crossprod(mu)/self$n + gamma)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, gamma = gamma, mu = mu)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(){
      return("diagonal normal-block model with fixed blocks... \n")
    }
  )
)

#' R6 class for normal-block model with fixed groups and spherical residual covariance
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param penalty to add on blocks precision matrix for sparsity
NB_fixed_blocks_spherical <- R6::R6Class(
  classname = "NB_fixed_blocks_spherical",
  inherit = NB_fixed_blocks,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    EM_initialize = function() {
      B      <- private$XtXm1 %*% t(self$X) %*% self$Y
      dm1    <- rep(1/mean((self$Y - self$X %*% B)^2), self$p)
      omegaQ <- diag(colSums(dm1 * private$C), self$Q, self$Q)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, gamma = private$gamma, mu = private$mu)
    },

    EM_step = function(B, dm1, omegaQ, gamma, mu) {

      ## E step
      gamma <- solve(omegaQ + dm1[1] * diag(colSums(private$C), self$Q, self$Q))
      mu    <- (self$Y - self$X %*% B) %*% private$C %*% gamma * dm1[1]

      ## M step
      YmmuCT <- self$Y - mu %*% t(private$C)
      B      <- private$XtXm1 %*% crossprod(self$X, YmmuCT)
      sigma2 <- mean((YmmuCT - self$X %*% B)^2) + sum(colMeans(private$C) * diag(gamma))
      omegaQ <- private$get_omegaQ(crossprod(mu)/self$n + gamma)

      list(B = B, dm1 = rep(1/sigma2, self$p), omegaQ = omegaQ, gamma = gamma, mu = mu)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(){
      return("spherical normal-block model with fixed blocks")
    }
  )
)
