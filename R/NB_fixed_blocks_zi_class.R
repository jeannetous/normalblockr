## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  NB_fixed_blocks_zi #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed groups
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param sparsity to add on blocks precision matrix
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
    #' @field zeros indicator matrix of zeros in Y
    zeros = NULL,

    #' @description Create a new [`NB_fixed_blocks_zi`] object.
    #' @param C group matrix C_jq = 1 if species j belongs to group q
    #' @return A new [`NB_fixed_blocks`] object
    initialize = function(Y, X, C, sparsity = 0) {
      if (!is.matrix(C)) stop("C must be a matrix.")
      super$initialize(Y = Y, X = X, ncol(C), sparsity = sparsity)
      self$C <- C
      self$zeros <- 1 * (Y == 0)
    },

    #' @description
    #' Update a [`NB_fixed_Q`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param omegaQ groups inverse variance matrix
    #' @param kappa zero-inflation probability
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal variances for posterior distribution of W
    #' @param rho variational parameter for zero-inflation probability
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB_fixed_Q`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, kappa = NA,
                      M = NA, S = NA, rho = NA, ll_list = NA) {
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
    M       = NA, # variational mean for W posterior distribution
    S       = NA, # variational diagonal variances for W posterior distribution
    kappa   = NA, # vector of zero-inflation probabilities
    rho     = NA, # posterior probabilities of zero-inflation

    compute_loglik_old  = function(B, dm1, omegaQ, kappa, M, S, rho) {
      R              <- self$Y - self$X %*% B
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      elbo <- -.5 * sum((1 - rho) * log(2 * pi)) + .5 * sum(log(dm1) * t(1 - rho))
      elbo <- elbo - .5 * sum((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C)) %*% diag(dm1))
      elbo <- elbo - .5 * self$n * self$Q * log(2 * pi) + .5 * self$n * log_det_omegaQ
      elbo <- elbo - .5 * sum((M %*% omegaQ) * M) - .5 * sum(S %*% diag(omegaQ))
      elbo <- elbo + sum(rho %*% log(kappa) + (1 - rho) %*% log(1 - kappa))
      elbo <- elbo + .5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * sum(log(S))
      elbo <- elbo - sum(rho * log(rho) + (1 - rho) * log(1 - rho))
      if (self$sparsity > 0) {
        elbo - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
      elbo
    },

    compute_loglik  = function(B, dm1, omegaQ, kappa, M, S, rho) {
      rho_bar        <- 1 - rho
      A              <- (self$Y - self$X %*% B - M %*% t(self$C))^2 + S %*% t(self$C)
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)
      J <- -.5 * sum(rho_bar %*% (log(2 * pi) - log(dm1)))
      J <- J - .5 * sum(rho_bar * A %*% diag(dm1))
      J <- J + .5 * self$n * log_det_omegaQ
      J <- J  + .5 * sum(log(S))
      J <- J + sum(rho %*% log(kappa) + rho_bar %*% log(1 - kappa))
      J <- J - sum(rho * log(rho)) - sum(rho_bar*log(rho_bar))
      if (self$sparsity > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + .5 * self$n *self$Q - .5 * sum(diag(omegaQ %*% (crossprod(M) + diag(colSums(S)))))
        J <- J - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
      J
    },

    EM_initialize = function() {
      init_model <- normal_zi$new(self$Y, self$X)
      init_model$optimize()
      B          <- init_model$model_par$B
      dm1        <- init_model$model_par$dm1
      omegaQ     <- t(self$C) %*% diag(dm1) %*% self$C
      kappa      <- init_model$model_par$kappa ## mieux qu'une 0-initialisation ?
      rho        <- init_model$model_par$rho
      G          <- solve(diag(colSums(dm1 * self$C)) + omegaQ)
      R          <- self$Y - self$X %*% B
      M          <- R %*% (dm1 * self$C) %*% G
      S          <- matrix(rep(0.1, self$n * self$Q), nrow = self$n)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, kappa = kappa, rho = rho, M = M,
           S = S)
    },

    zi_nb_fixed_blocks_obj_grad_M = function(M_vec, dm1_1mrho, omegaQ, R) {
      M    <- matrix(M_vec, nrow = self$n, ncol = self$Q)
      RmMC <- R - M %*% t(self$C)
      MO   <- M %*% omegaQ

      obj  <- -.5 * sum(dm1_1mrho * RmMC^2) - .5 * sum(MO * M)
      grad <- (dm1_1mrho * RmMC) %*% self$C - MO

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_blocks_nlopt_optim_M = function(M0, dm1, omegaQ, B, rho) {
      res <- nloptr::nloptr(
        x0 = as.vector(M0),
        eval_f = private$zi_nb_fixed_blocks_obj_grad_M,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        dm1_1mrho = t(dm1 * t(1 - rho)),
        omegaQ    = omegaQ,
        R         = self$Y - self$X %*% B
      )
      newM <- matrix(res$solution, nrow = self$n, ncol = self$Q)
      newM
    },

    zi_nb_fixed_blocks_obj_grad_B = function(B_vec, dm1_1mrho, MC) {
      R   <- self$Y - self$X %*% matrix(B_vec, nrow = self$d, ncol = self$p)
      RmMC <- R - MC

      obj  <- -.5 * sum(dm1_1mrho * RmMC^2)
      grad <- crossprod(self$X, dm1_1mrho * RmMC)

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_blocks_nlopt_optim_B = function(B0, dm1, omegaQ, M, rho) {
      res <- nloptr::nloptr(
        x0 = as.vector(B0),
        eval_f = private$zi_nb_fixed_blocks_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        dm1_1mrho = t(dm1 * t(1 - rho)),
        MC = M %*% t(self$C)
      )
      newB <- matrix(res$solution, nrow = self$d, ncol = self$p)
      newB
    },

    EM_step = function(B, dm1, omegaQ, kappa, M, S, rho) {
      R   <- self$Y - self$X %*% B
      ones <- as.vector(rep(1, self$n))

      # E step
      M <- private$zi_nb_fixed_blocks_nlopt_optim_M(M, dm1, omegaQ, B, rho)
      S <-  1 / sweep((1 - rho) %*% (dm1 * self$C), 2, diag(omegaQ), "+")
      A <- ((R - M %*% t(self$C))^2 + S %*% t(self$C))
      nu <- log(2 * pi) - outer(ones, log(dm1)) + A %*% diag(dm1)
      rho <- 1 / (1 + exp(-.5 * nu) * outer(ones, (1 - kappa) / kappa))
      rho <- check_one_boundary(check_zero_boundary(self$zeros * rho))

      # M step
      B <- private$zi_nb_fixed_blocks_nlopt_optim_B(B, dm1, omegaQ, M, rho)
      dm1   <- colSums(1 - rho) / colSums((1 - rho) * A)
      kappa <- colMeans(rho)
      sigmaQ <- (1 / self$n) * (crossprod(M) + diag(colSums(S)))

      if (self$sparsity == 0 ) {
        omegaQ <- solve(sigmaQ)
      } else {
        glasso_out <- glassoFast::glassoFast(sigmaQ, rho = self$sparsity * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) stop("GLasso fails")
        omegaQ <- Matrix::symmpart(glasso_out$wi)
      }

      list(B = B, dm1 = dm1, omegaQ = omegaQ,  kappa = kappa, rho = rho, M = M,
           S = S)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field var_par a list with variational parameters
    var_par  = function() list(M = private$M, S = private$S, rho = private$rho),
    #' @field nb_param number of parameters in the model
    nb_param = function() as.integer(super$nb_param + self$p),
    #' @field model_par a list with model parameters: B (covariates), dm1 (species variance), omegaQ (groups precision matrix), kappa (zero-inflation probabilities)
    model_par  = function() {
      par       <- super$model_par
      par$kappa <- private$kappa
      par
    },
    #' @field clustering given as a list of labels
    clustering = function() get_clusters(self$C),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      ent <- 0.5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * n * sum(log(private$S))
      ent <- ent - sum(private$rho * log(private$rho) + (1 - private$rho) * log(1 - private$rho))
      return(ent)
    }
  )
)
