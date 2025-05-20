## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZINB_fixed_Q ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for zero-inflated normal-block model with fixed number of groups#'
#' @export
ZINB_fixed_Q <- R6::R6Class(
  classname = "ZINB_fixed_Q",
  inherit   = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field fixed_tau whether tau should be fixed at clustering_init during optimization, useful for stability selection
    fixed_tau = NULL,

    #' @description Create a new [`ZINB_fixed_Q`] object.
    #' @param data object of NBData class, with responses and design matrix
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param Q required number of groups
    #' @param control structured list of more specific parameters
    #' @return A new [`ZINB_fixed_Q`] object
    initialize = function(data, Q, sparsity = 0, control = NB_control()) {
      super$initialize(data, Q, sparsity, control)
      self$fixed_tau  <- control$fixed_tau
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS -------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    compute_loglik  = function(B, dm1, OmegaQ, alpha, kappa, M, S, C) {
      log_det_OmegaQ <- as.numeric(determinant(OmegaQ, logarithm = TRUE)$modulus)

      J <- -.5 * self$data$npY * log(2 * pi * exp(1)) + .5 * sum(self$data$nY * log(dm1))
      J <- J + .5 * self$n * log_det_OmegaQ + .5 * sum(log(S))
      J <- J +  sum(C %*% log(alpha)) - sum(xlogx(C))
      J <- J + private$ZI_cond_mean

      if (private$sparsity_ > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + .5 * self$n * self$Q - .5 * sum(diag(OmegaQ %*% (crossprod(M) + diag(colSums(S), nrow = self$Q, ncol = self$Q))))
        J <- J - private$sparsity_ * sum(abs(self$sparsity_weights * OmegaQ))
      }
      J
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference -----------------------

    get_heuristic_parameters = function() {
      zi_diag <- private$zi_diag_normal_inference()
      if (anyNA(private$C))
        private$C <- private$clustering_approx(zi_diag$R)
      private$C <- check_one_boundary(check_zero_boundary(private$C))
      SigmaQ <- private$heuristic_SigmaQ_from_Sigma(cov(zi_diag$R))
      OmegaQ <- private$get_OmegaQ(SigmaQ)
      list(B = zi_diag$B, dm1 = zi_diag$dm1, OmegaQ = OmegaQ,
           alpha = colMeans(private$C), kappa = private$kappa,
           C = private$C)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated inference ------------------------

    EM_initialize = function() {
      c(private$get_heuristic_parameters(),  list(
          M = matrix(rep(0, self$n * self$Q), nrow = self$n),
          S = matrix(rep(0.1, self$n * self$Q), nrow = self$n)
        )
      )
    },

    EM_step = function(B, dm1, OmegaQ, alpha, kappa, C, M, S) {

      R <- self$data$Y - self$data$X %*% B
      DM1 <- matrix(dm1, self$n, self$p, byrow = TRUE) * self$data$zeros_bar

      # E step
      M <- private$zi_NB_fixed_Q_nlopt_optim_M(M, DM1, R, OmegaQ, C)
      S <-  1 / sweep(DM1 %*% C, 2, diag(OmegaQ), "+")
      if (self$Q > 1 & !self$fixed_tau) {
        eta <- -.5 * crossprod(DM1, (M^2 + S))
        eta <- eta + crossprod(DM1 * R, M) + outer(rep(1, self$p), log(alpha)) - 1
        C <- t(check_zero_boundary(check_one_boundary(apply(eta, 1, softmax))))
      }
      A <- R^2 - 2 * R * tcrossprod(M,C) + tcrossprod(M^2 + S, C)

      # M step
      B   <- private$zi_NB_fixed_Q_nlopt_optim_B(B, DM1, M, C)
      dm1  <- switch(private$res_covariance,
                     "diagonal"  = self$data$nY / colSums(self$data$zeros_bar * A),
                     "spherical" = rep(self$data$npY / sum(self$data$zeros_bar * A), self$p))
      alpha <- colMeans(C)
      OmegaQ <- private$get_OmegaQ(crossprod(M)/self$n + diag(colMeans(S), self$Q, self$Q))

      list(B = B, dm1 = dm1, OmegaQ = OmegaQ, alpha = alpha, kappa = kappa,
           C = C, M = M, S = S)
    },

    zi_NB_fixed_Q_obj_grad_M = function(M_vec, DM1, DM1RC, DM1C, R, C, OmegaQ) {
      M    <- matrix(M_vec, nrow = self$n, ncol = self$Q)
      MO   <- M %*% OmegaQ
      grad <- DM1RC - DM1C * M - MO
      obj <- -.5 * ( sum(DM1 * (M^2 %*% t(C)) - 2*R * (M %*% t(C))) + sum(MO * M) )
      res  <- list("objective" = -obj, "gradient"  = -grad)
      res
    },

    zi_NB_fixed_Q_nlopt_optim_M = function(M0, DM1, R, OmegaQ, C) {
      M0_vec <- as.vector(M0)
      res <- nloptr::nloptr(
        x0 = M0_vec,
        eval_f = private$zi_NB_fixed_Q_obj_grad_M,
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          maxeval = 100
        ),
        DM1    = DM1,
        DM1RC  = (DM1 * R) %*% C,
        DM1C   = DM1 %*% C,
        R      = R,
        C      = C,
        OmegaQ = OmegaQ
      )
      newM <- matrix(res$solution, nrow = self$n, ncol = self$Q)
      newM
    },

    zi_NB_fixed_Q_obj_grad_B = function(B_vec, dm1_mat, MC) {
      R    <- self$data$Y - self$data$X %*% matrix(B_vec, nrow = self$d, ncol = self$p)
      grad <- crossprod(self$data$X, dm1_mat * (R - MC))
      obj  <- -.5 * sum(dm1_mat * (R^2 - 2 * R * MC))
      res  <- list("objective" = -obj, "gradient"  = -grad)
      res
    },

    zi_NB_fixed_Q_nlopt_optim_B = function(B0, dm1_mat, M, C) {
      res <- nloptr::nloptr(
        x0 = as.vector(B0),
        eval_f = private$zi_NB_fixed_Q_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          maxeval = 100
        ),
        dm1_mat = dm1_mat,
        MC = tcrossprod(M, C)
      )
      newB <- matrix(res$solution, nrow = self$d, ncol = self$p)
      newB
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function() super$nb_param + self$p + self$Q - 1, # adding kappa and alpha
    #' @field var_par a list with variational parameters
    var_par  = function() {list(M = private$M, S = private$S, tau = private$C)},
    #' @field model_par a list with model parameters: B (covariates), dm1 (species variance), OmegaQ (blocks precision matrix), kappa (zero-inflation probabilities)
    model_par  = function() {
      par       <- super$model_par
      par$kappa <- private$kappa
      par$alpha <- private$alpha
      par
    },
    #' @field entropy Entropy of the conditional distribution
    entropy    = function() {
      if (!private$approx) {
        res <- 0.5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * sum(log(private$S))
        res <- res - sum(xlogx(private$C))
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model Y values predicted by the model
    fitted = function(){
      if (private$approx) {
        res <- self$data$X %*% private$B
      } else {
        res <- self$data$X %*% private$B + private$M %*% t(private$C)
      }
      res <- res * self$data$zeros_bar
      res
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function(value)
    {paste("zero-inflated", private$res_covariance, "normal-block model with", self$Q, "unknown blocks")}
  )
)

