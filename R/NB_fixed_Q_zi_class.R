## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  NB_fixed_Q_zi #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed blocks
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param Q number of blocks in the model
#' @param sparsity to add on blocks precision matrix
#' @export
NB_fixed_Q_zi <- R6::R6Class(
  classname = "NB_fixed_Q_zi",
  inherit = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field zeros indicator matrix of zeros in Y
    zeros = NULL,

    #' @description Create a new [`NB_fixed_Q_zi`] object.
    #' @param C block matrix C_jq = 1 if species j belongs to block q
    #' @return A new [`NB_fixed_blocks`] object
    initialize = function(Y, X, Q, sparsity = 0) {
      super$initialize(Y = Y, X = X, Q, sparsity = sparsity)
      self$zeros <- 1 * (Y == 0)
    },

    #' @description
    #' Update a [`NB_fixed_Q`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param omegaQ blocks inverse variance matrix
    #' @param alpha vector of blocks probabilities
    #' @param kappa zero-inflation probability
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param tau posterior probabilities for block affectation
    #' @param rho variational parameter for zero-inflation probability
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB_fixed_Q`] object
    update = function(B = NA, dm1 = NA, alpha = NA, omegaQ = NA, kappa = NA,
                      M = NA, S = NA, tau = NA, rho = NA, ll_list = NA) {
      super$update(B, dm1, omegaQ, ll_list)
      if (!anyNA(alpha)) private$alpha <- alpha
      if (!anyNA(kappa)) private$kappa <- kappa
      if (!anyNA(M))     private$M     <- M
      if (!anyNA(S))     private$S     <- S
      if (!anyNA(tau))   private$tau   <- tau
      if (!anyNA(rho))   private$rho   <- rho
    }),


  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    kappa   = NA,   # vector of zero-inflation probabilities
    alpha   = NA, # vector of groups probabilities
    M       = NA, # variational mean for posterior distribution of W
    S       = NA, # variational diagonal of variances for posterior distribution of W
    tau     = NA, # posterior probabilities for group affectation
    rho     = NA,   # posterior probabilities of zero-inflation

    compute_loglik  = function(B, dm1, omegaQ, alpha, kappa, M, S, tau, rho) {
      R              <- self$Y - self$X %*% B
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      elbo <- -.5 * sum((1 - rho) * log(2 * pi)) + .5 * sum(log(dm1) * t(1 - rho))
      elbo <- elbo - .5 * sum((1 - rho) * (R^2 - 2 * R * (M %*% t(tau)) + (M^2 + S) %*% t(tau)) %*% diag(dm1))
      elbo <- elbo - .5 * self$n * self$Q * log(2 * pi) + .5 * self$n * log_det_omegaQ
      elbo <- elbo - .5 * sum((M %*% omegaQ) * M) - .5 * sum(S %*% diag(omegaQ))
      elbo <- elbo + sum(tau %*% log(alpha))
      elbo <- elbo + sum(rho %*% log(kappa) + (1 - rho) %*% log(1 - kappa))
      elbo <- elbo + .5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * sum(log(S))
      elbo <- elbo - sum(tau * log(tau))
      elbo <- elbo - sum(rho * log(rho) + (1 - rho) * log(1 - rho))
      if (self$sparsity > 0) {
        elbo <- elbo - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
      elbo
    },

    EM_initialize = function() {
      init_model <- normal_zi$new(self$Y, self$X)
      init_model$optimize()
      B          <- init_model$model_par$B
      dm1        <- init_model$model_par$dm1
      R          <- self$Y - self$X %*% B
      R[self$Y == 0]  <- 0
      cl         <- kmeans(t(R), self$Q, nstart = 30)$cluster
      tau        <- check_one_boundary(check_zero_boundary(as_indicator(cl)))
      alpha      <- colMeans(tau)
      omegaQ     <- t(tau) %*% diag(dm1) %*% tau
      kappa      <- init_model$model_par$kappa ## mieux qu'une 0-initialisation ?
      rho        <- init_model$model_par$rho
      G          <- solve(diag(colSums(as.vector(dm1) * tau), self$Q, self$Q) + omegaQ)
      M          <- t(dm1 * t(R)) %*% tau %*% G
      S          <- matrix(rep(0.1, self$n * self$Q), nrow = self$n)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, kappa = kappa,
           M = M, S = S, tau = tau, rho = rho)
    },

    zi_nb_fixed_Q_obj_grad_M = function(M_vec, B, dm1, omegaQ, kappa, S,
                                             tau, rho) {
      M    <- matrix(M_vec, nrow = self$n, ncol = self$Q)
      R    <- self$Y - self$X %*% B
      MO   <- M %*% omegaQ
      dm1T <- dm1 * tau

      grad <- ((1 - rho) * R) %*% dm1T - ((1 - rho) %*% dm1T) * M - MO
      obj  <- sum((1 - rho) * (R * (M %*% t(dm1T)) - .5 * M^2 %*% t(dm1T))) - .5 * sum(MO * M)

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_Q_nlopt_optim_M = function(M0, B, dm1, omegaQ, kappa,
                                           S, tau, rho) {
      M0_vec <- as.vector(M0)
      res <- nloptr::nloptr(
        x0 = M0_vec,
        eval_f = private$zi_nb_fixed_Q_obj_grad_M,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        B      = B,
        dm1    = dm1,
        omegaQ = omegaQ,
        kappa  = kappa,
        S      = S,
        tau    = tau,
        rho    = rho
      )
      newM <- matrix(res$solution, nrow = self$n, ncol = self$Q)
      newM
    },

    zi_nb_fixed_Q_obj_grad_B = function(B_vec, dm1_1mrho, M, tau) {
      B    <- matrix(B_vec, nrow = self$d, ncol = self$p)
      R    <- self$Y - self$X %*% B
      MT   <- M %*% t(tau)

      grad <- crossprod(self$X, dm1_1mrho * (R - MT))
      obj  <- -.5 * sum(dm1_1mrho * (R^2 - 2 * R * MT))

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_Q_nlopt_optim_B = function(B0, dm1, omegaQ, kappa, M, S,
                                                tau, rho) {
      B0_vec <- as.vector(B0)
      dm1_1mrho <- t(dm1 * t(1 - rho))
      res <- nloptr::nloptr(
        x0 = B0_vec,
        eval_f = private$zi_nb_fixed_Q_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        dm1_1mrho = dm1_1mrho,
        M      = M,
        tau    = tau
      )
      newB <- matrix(res$solution, nrow = self$d, ncol = self$p)
      newB
    },

    EM_step = function(B, dm1, omegaQ, alpha, kappa, M, S, tau, rho) {
      R <- self$Y - self$X %*% B
      ones <- as.vector(rep(1, self$n))

      # E step
      M <- private$zi_nb_fixed_Q_nlopt_optim_M(M, B, dm1, omegaQ, kappa, S, tau, rho)
      S <-  1 / sweep((1 - rho) %*% (dm1 * tau), 2, diag(omegaQ), "+")

      nu <- tcrossprod(ones, as.vector(rep(1, self$p)) * log(2 * pi) - log(dm1))
      nu <- nu + (R^2 - 2 * R * (M %*% t(tau)) + (M^2 + S) %*% t(tau)) %*% diag(dm1)
      rho <- 1 / (1 + exp(-.5 * nu) * tcrossprod(ones, (1 - kappa) / kappa))
      rho <- check_one_boundary(check_zero_boundary((self$Y == 0) * rho))

      if (self$Q > 1) {
        eta <- -.5 * dm1 * t(1 - rho) %*% (M^2 + S)
        eta <- eta + dm1 * t((1 - rho) * R) %*% M  + outer(rep(1, self$p), log(alpha)) - 1
        tau <- t(check_zero_boundary(check_one_boundary(apply(eta, 1, softmax))))
      }

      # M step
      B   <- private$zi_nb_fixed_Q_nlopt_optim_B(B, dm1, omegaQ, kappa, M, S, tau, rho)
      nSigmaQ <- crossprod(M) + diag(colSums(S), self$Q, self$Q)
      if (self$sparsity == 0) {
        omegaQ <- self$n * solve(nSigmaQ)
      } else {
        glasso_out <- glassoFast::glassoFast(nSigmaQ/self$n, rho = self$sparsity * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) stop("GLasso fails")
        omegaQ <- Matrix::symmpart(glasso_out$wi)
      }
      alpha <- colMeans(tau)
      kappa <- check_one_boundary(check_zero_boundary(colMeans(rho)))
      dm1   <- colSums(1 - rho) / colSums((1 - rho) * (R^2 - 2 * R * (M %*% t(tau)) + (M^2 + S) %*% t(tau)))

      list(B = B, dm1 = dm1, alpha = alpha, omegaQ = omegaQ, kappa = kappa,
           M = M, S = S, tau = tau, rho = rho)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par a list with model parameters: B (covariates), dm1 (species variance), omegaQ (blocks precision matrix), kappa (zero-inflation probabilities)
    model_par  = function() {
      par       <- super$model_par
      par$kappa <- private$kappa
      par$alpha <- private$alpha
      par
    },
    #' @field var_par a list with variational parameters
    var_par  = function() {list(M = private$M, S = private$S,
                                     rho = private$rho, tau = private$tau)},
    #' @field clustering given as a list of labels
    clustering = function() get_clusters(private$tau),
    #' @field nb_param number of parameters in the model
    nb_param = function() {as.integer(super$nb_param + self$p + self$Q - 1)},
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      ent <- 0.5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * self$n * sum(log(private$S))
      ent <- ent - sum(xlogx(private$tau))
      ent <- ent - sum(private$rho * log(private$rho) + (1 - private$rho) * log(1 - private$rho))
      ent
    }
  )
)
