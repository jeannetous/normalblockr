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
      self$C <- C
      self$Q <- ncol(C)
      self$zeros <- 1 * (Y == 0)
      super$initialize(Y = Y, X = X, sparsity = sparsity)
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


    compute_loglik  = function(B, dm1, omegaQ, kappa, M, S, rho) {
      R              <- self$Y - self$X %*% B
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      elbo <- -.5 * sum((1 - rho) * log(2 * pi))
      elbo <- elbo + .5 * sum(t(t(1 - rho) * log(dm1)))
      elbo <- elbo - .5 * sum(dm1 * ((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C))))
      elbo <- elbo - .5 * self$n * self$Q * log(2 * pi) + .5 * self$n * log_det_omegaQ
      elbo <- elbo - .5 * sum((M %*% omegaQ) * M) - .5 * sum(S %*% diag(omegaQ))
      elbo <- elbo + sum(rho %*% log(kappa) + (1 - rho) %*% log(1 - kappa))
      elbo <- elbo + .5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * sum(log(S))
      elbo <- elbo - sum(rho * log(rho) + (1 - rho) * log(1 - rho))
      if (self$sparsity == 0) {
        elbo
      }else {
        elbo - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
    },

    EM_initialize = function() {
      init_model <- normal_zi$new(self$Y, self$X)
      init_model$optimize()
      B          <- init_model$model_par$B
      dm1        <- init_model$model_par$dm1
      omegaQ     <- t(self$C) %*% diag(dm1) %*% self$C
      kappa      <- init_model$model_par$kappa ## mieux qu'une 0-initialisation ?
      rho        <- init_model$model_par$rho
      G          <- solve(diag(colSums(as.vector(dm1) * self$C)) + omegaQ)
      R          <- self$Y - self$X %*% B
      M          <- dm1 * R %*% self$C %*% G
      S          <- matrix(rep(0.1, self$n * self$Q), nrow = self$n)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, kappa = kappa, rho = rho, M = M,
           S = S)
    },

    zi_nb_fixed_blocks_obj_grad_M = function(M_vec, B, dm1, omegaQ, kappa, S, rho) {
      M    <- matrix(M_vec, nrow = self$n, ncol = self$Q)
      R    <- self$Y - self$X %*% B
      grad <- - ( t(t((1 - rho) * R) * dm1) %*% self$C - ( t(dm1 * t(1 - rho)) %*% self$C) * M - M %*% omegaQ)

      obj  <- - .5 * sum(t(dm1 * t((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C)))))
      obj  <- obj - .5 * sum((M %*% omegaQ) * M)

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_blocks_nlopt_optim_M = function(M0, B, dm1, omegaQ, kappa,
                                                S, rho) {
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
        kappa  = kappa,
        S      = S,
        rho    = rho
      )
      newM <- matrix(res$solution, nrow = self$n, ncol = self$Q)
      newM
    },

    zi_nb_fixed_blocks_obj_grad_B = function(B_vec, dm1, omegaQ, kappa,
                                             M, S, rho) {
      B    <- matrix(B_vec, nrow = self$d, ncol = self$p)
      R    <- self$Y - self$X %*% B
      grad <- t(self$X) %*% t(dm1 * t((1 - rho) * (R - M %*% t(self$C))))
      obj  <- - .5 * sum(dm1 * ((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C))))

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_blocks_nlopt_optim_B = function(B0, dm1, omegaQ, kappa, M,
                                                S, rho) {
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
        kappa  = kappa,
        M      = M,
        S      = S,
        rho    = rho
      )
      newB <- matrix(res$solution, nrow = self$d, ncol = self$p)
      newB
    },

    EM_step = function(B, dm1, omegaQ, kappa, M, S, rho) {
      R   <- self$Y - self$X %*% B
      # E step
      M <- private$zi_nb_fixed_blocks_nlopt_optim_M(M, B, dm1, omegaQ,
                                                    kappa, S,  rho)
      S <-  1 / sweep(((1 - rho) %*% (dm1 * self$C)), 2, diag(omegaQ), "+")

      nu <- as.vector(rep(1, self$n)) %*% t(as.vector(rep(1, self$p))) * log(2 * pi)
      nu <- nu - as.vector(rep(1, self$n)) %*% t(log(dm1))
      nu <- nu + t(t(((self$X %*% B)^2 +  2 * (self$X %*% B) * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C))) * dm1)
      rho <- 1 / (1 + exp(-.5 * nu) * as.vector(rep(1, self$n)) %*% t((1 - kappa) / kappa))
      rho <- check_one_boundary(check_zero_boundary(self$zeros * rho))

      # M step
      B        <- private$zi_nb_fixed_blocks_nlopt_optim_B(B, dm1, omegaQ,
                                                           kappa, M, S, rho)

      if (self$sparsity == 0 ) {
        omegaQ <- self$n * solve((t(M) %*% M) + diag(self$n * diag(S)))
      }else {
        sigma_hat <- (1 / self$n) * (t(M) %*% M) + diag(self$n * diag(S))
        glasso_out <- glassoFast::glassoFast(sigma_hat, rho = self$sparsity * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) stop("GLasso fails")
        omegaQ <- Matrix::symmpart(glasso_out$wi)
      }
      kappa    <- check_one_boundary(check_zero_boundary(colMeans(rho)))
      dd       <- (1 / (t(1 - rho) %*% as.vector(rep(1, self$n)))) * t((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + ((M^2 + S) %*% t(self$C)))) %*% as.vector(rep(1, self$n))
      dm1      <- as.vector(1 / dd)

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
    nb_param = function() as.integer(super$nb_param + 2 * self$n * self$Q + self$p * (self$n + 1)),
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
