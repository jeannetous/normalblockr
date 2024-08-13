## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  NB_fixed_blocks_zi #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed groups
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood threshold under which optimization stops
#' @export
NB_fixed_blocks_zi <- R6::R6Class(
  classname = "NB_fixed_blocks_zi",
  inherit = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field C the matrix of species groups
    C = NULL,

    #' @description Create a new [`NB_fixed_blocks_zi`] object.
    #' @param C group matrix C_jq = 1 if species j belongs to group q
    #' @return A new [`NB_fixed_blocks`] object
    initialize = function(Y, X, C, niter = 50, threshold = 1e-4) {
      if (!is.matrix(C)) stop("C must be a matrix.")
      self$C <- C
      self$Q <- ncol(C)
      super$initialize(Y = Y, X = X, niter = niter, threshold = threshold)
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
    update = function(B = NA, dm1 = NA, omegaQ = NA, kappa = NA,
                      M = NA, S = NA, rho = NA, ll_list=NA) {
      super$update(B, dm1, omegaQ, ll_list)
      if (!anyNA(kappa)) private$kappa <- kappa
      if (!anyNA(M))     private$M     <- M
      if (!anyNA(S))     private$S     <- S
      if (!anyNA(rho))   private$rho   <- rho
    }),


  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    M       = NA, # variational mean for posterior distribution of W
    S       = NA, # variational diagonal of variances for posterior distribution of W
    kappa   = NA,   # vector of zero-inflation probabilities
    rho     = NA,   # posterior probabilities of zero-inflation


    compute_loglik  = function(B, dm1, omegaQ, kappa, rho, M, S) {
      R              <- self$Y - self$X %*% B
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      elbo <- -.5 * sum((1 - rho) * log(2 * pi))
      ### ligne ci-dessous à récrire pour enlever le diag(dm1)
      elbo <- elbo + .5 * sum((1 - rho) %*% diag(log(dm1)))
      elbo <- elbo - .5 * sum(dm1 * ((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C))) )
      elbo <- elbo - .5 * self$n * self$Q * log(2 * pi) + .5 * self$n * log_det_omegaQ
      elbo <- elbo - .5 * sum((M %*% omegaQ) * M) - .5 * sum(S %*% diag(omegaQ))
      elbo <- elbo + sum(rho %*% log(kappa) + (1 - rho) %*% log(1 - kappa))
      elbo <- elbo + .5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * sum(log(S))
      elbo <- elbo - sum(rho * log(rho) + (1 - rho) * log(1 - rho))
      elbo
    },

    EM_initialize = function() {
      init_model <- zi_normal$new(self$Y, self$X)
      init_model$optimize()
      B          <- init_model$model_par$B
      dm1        <- init_model$model_par$dm1
      omegaQ     <- t(self$C) %*% diag(init_model$model_par$dm1) %*% self$C
      kappa      <- init_model$model_par$kappa ## mieux qu'une 0-initialisation ?
      rho        <- init_model$model_par$rho
      G          <- solve(diag(colSums(as.vector(dm1) * self$C)) + omegaQ)
      R          <- self$Y - self$X %*% B
      M          <- dm1 * R %*% self$C %*% G
      S          <- matrix(rep(0.1, self$n * self$Q), nrow=self$n)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, kappa = kappa, rho = rho, M = M,
           S = S)
    },

    zi_nb_fixed_blocks_obj_grad_M = function(M_vec, B, dm1, omegaQ, rho, kappa, S) {
      M    <- matrix(M_vec, nrow = self$n, ncol = self$Q)
      R    <- self$Y - X %*% B
      grad <- - (((1 - rho) * R) %*%  diag(dm1) %*% self$C - ((1 - rho) %*% diag(dm1) %*% self$C) * M - M %*% omegaQ)

      obj  <- - .5 * sum(((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C))) %*% diag(dm1))
      obj  <- obj - .5 * sum((M %*% omegaQ) * M)

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_blocks_nlopt_optim_M = function(M0, B, dm1, omegaQ, rho, kappa,
                                                S) {
      M0_vec <- as.vector(M0)
      res <- nloptr::nloptr(
        x0 = M0_vec,
        eval_f = private$zi_nb_fixed_blocks_obj_grad_M,
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        B      = B,
        dm1    = dm1,
        omegaQ = omegaQ,
        rho    = rho,
        kappa  = kappa,
        S      = S
      )
      newM <- matrix(res$solution, nrow = self$n, ncol = self$Q)
      newM
    },

    zi_nb_fixed_blocks_obj_grad_B = function(B_vec, dm1, omegaQ, rho, kappa,
                                             M, S) {
      B    <- matrix(B_vec, nrow = self$d, ncol = self$p)
      R    <- self$Y - X %*% B
      grad <- t(X) %*% ((1 - rho) * (R - M %*% t(C))) %*% diag(dm1)
      obj  <- - .5 * sum(dm1 * ((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C))) )

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_blocks_nlopt_optim_B = function(B0, dm1, omegaQ, rho, kappa, M,
                                                S) {
      B0_vec <- as.vector(B0)
      res <- nloptr::nloptr(
        x0 = B0_vec,
        eval_f = private$zi_nb_fixed_blocks_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        dm1    = dm1,
        omegaQ = omegaQ,
        rho    = rho,
        kappa  = kappa,
        M      = M,
        S      = S
      )
      newB <- matrix(res$solution, nrow = self$d, ncol = self$p)
      newB
    },

    EM_step = function(B, dm1, omegaQ, kappa, rho, M, S) {
      R   <- self$Y - self$X %*% B
      # rho_ind <- which(self$Y == 0, arr.ind=TRUE) ## à modifier

      # E step
      M <- private$zi_nb_fixed_blocks_nlopt_optim_M(M, B, dm1, omegaQ, rho,
                                                    kappa, S)
      S <-  1/sweep(((1 - rho) %*% (dm1 * C)), 2, diag(omegaQ), "+")

      eta <- as.vector(rep(1, n)) %*% t(as.vector(rep(1, p))) * log(2 * pi)
      eta <- eta - as.vector(rep(1, n)) %*% t(log(dm1))
      eta <- eta + t(t(((self$X %*% B)^2 +  2 * (self$X %*% B) * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C))) * dm1)
      rho <- 1 /(1 + exp(-.5 * eta) * as.vector(rep(1, n)) %*% t((1 - kappa)/kappa))
      rho <- check_one_boundary(check_zero_boundary((self$Y == 0)* rho))

      # M step
      B        <- private$zi_nb_fixed_blocks_nlopt_optim_B(B, dm1, omegaQ, rho,
                                                           kappa, M, S)
      omegaQ   <- self$n * solve((t(M) %*% M) + diag(self$n * diag(S)))
      kappa    <- check_boundaries(colMeans(rho))
      dd       <- (1/(t(1 - rho) %*% as.vector(rep(1, self$n)))) * t((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + ((M^2 + S) %*% t(self$C)))) %*% as.vector(rep(1, self$n))
      dm1      <- as.vector(1/dd)

      list(B = B, dm1 = dm1, omegaQ = omegaQ,  kappa = kappa, rho = rho, M = M,
           S = S)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with variational parameters
    posterior_par  = function() {list(M = private$M, S = private$S, rho = private$rho)},
    #' @field model_par a list with model parameters: B (covariates), dm1 (species variance), omegaQ (groups precision matrix), kappa (zero-inflation probabilities)
    model_par  = function() {
      par       = super$model_par
      par$kappa = private$kappa
      par}
  )
)

