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
  inherit = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field C the matrix of species groups
    C = NULL,

    #' @description Create a new [`NB`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param C group matrix C_jq = 1 if species j belongs to group q
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    #' @return A new [`NB_unknown`] object
    initialize = function(Y, X, C, niter = 50, threshold = 1e-4) {
      super$initialize(Y, X, niter, threshold)
      if (!is.matrix(C)) stop("C must be a matrix.")
      self$C <- C
      private$Q <- ncol(C)
    },

    #' @description
    #' Update a [`NB_unknown`] object
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
    },

    #' @description returns the model parameters B, dm1, omegaQ, n, p, d, Q, gamma, mu
    #' @return A list containing the model parameters  B, dm1, omegaQ, n, p, d, Q, gamma, mu
    get_model_parameters = function() {
      parameters = super$get_model_parameters()
      c(parameters, list("gamma" = private$gamma, "mu" = private$mu))
    },
    #' @description returns the model variables Y, X, C
    #' @return A list containing the model parameters Y, X, C
    get_model_variables = function() {
      list(Y = self$Y, X = self$X, C = self$C)
    }),


  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    gamma   = NULL, # variance of  posterior distribution of W
    mu      = NULL,  #  mean for posterior distribution of W


    loglik  = function(Y, X, C, B, dm1, omegaQ, gamma, mu) {
      ## problem dimensions
      n   <- private$n; p <-  private$p; d <-  private$d; Q <-  private$Q

      Dm1 <- diag(dm1)
      R   <- Y - X %*% B
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      J <- - .5 * n * (p + Q) * log(2 * pi) + .5 * n * sum(log(dm1))
      J <- J - .5 * sum(R %*% Dm1 %*% t(R)) + sum(R %*% Dm1 %*% C %*% t(mu))
      J <- J - .5 * n * sum(diag(t(C) %*% Dm1 %*% C %*% gamma))
      J <- J - .5 * sum(diag(mu %*% t(C) %*% Dm1 %*% C %*% t(mu)))
      J <- J + .5 * n * log_det_omegaQ
      J <- J - .5 * n * sum(diag(omegaQ %*% gamma))
      J <- J - .5 * n * sum(diag(mu %*% omegaQ %*% t(mu)))
      J
    },

    EM_initialize = function(Y, X, C) {
      B      <- private$XtXm1%*% t(X) %*% Y
      R      <- t(Y - X %*% B)
      Ddiag  <- check_one_boundary(check_zero_boundary(diag(cov(t(R)))))
      D      <- diag(Ddiag)
      dm1    <- 1 / Ddiag
      Dm1    <- diag(dm1)
      gamma  <- solve(t(C) %*% Dm1 %*% C)
      mu     <- t(gamma %*% t(C) %*% Dm1 %*% R)
      omegaQ <- solve(gamma + (1/private$n) * t(mu) %*% mu)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, gamma = gamma, mu = mu)
    },

    EM_step = function(Y, X, C, B, dm1, omegaQ, gamma, mu) {
      Dm1    <- diag(dm1)
      R      <- t(Y - X %*% B)
      ## E step
      gamma  <- solve(omegaQ + t(C) %*% Dm1 %*% C)
      mu     <- t(gamma %*% t(C) %*% Dm1 %*% R)
      ## M step
      B      <- private$XtXm1 %*% t(X) %*% (Y - mu %*% t(C))
      Ddiag  <- (1/private$n) * (t(Y - X%*% B) - C %*% t(mu))^2 %*% as.vector(rep(1, private$n)) + diag(C %*% gamma %*% t(C))
      dm1    <- as.vector(1/Ddiag)
      omegaQ <- solve(gamma + (1/private$n) * t(mu) %*% mu)

      list(B = B, dm1 = dm1, omegaQ = omegaQ, gamma = gamma, mu = mu)
    }
  )
)
