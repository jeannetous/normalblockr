## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_zi_fixed_blocks_fixed_sparsity ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param data object of normal_data class, with responses and design matrix
#' @param C clustering matrix C_jq = 1 if species j belongs to cluster q
#' @param sparsity to apply on variance matrix when calling GLASSO
#' @param control structured list of more specific parameters, to generate with NB_control
NB_zi_fixed_blocks <- R6::R6Class(
  classname = "NB_zi_fixed_blocks",
  inherit   = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field zeros indicator matrix of zeros in Y
    zeros = NULL,

    #' @description Create a new [`NB_zi_fixed_blocks_fixed_sparsity`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param C clustering matrix C_jq = 1 if species j belongs to cluster q
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NB_zi_fixed_blocks_fixed_sparsity`] object
    initialize = function(data, C, sparsity = 0, control = NB_control()) {
      stopifnot("C must be a matrix" = is.matrix(C))
      stopifnot("There cannot be empty clusters" = min(colSums(C)) > 0)
      super$initialize(data, ncol(C), sparsity, control = control)
      private$C  <- C
      self$zeros <- 1 * (data$Y == 0)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    compute_loglik  = function(B, OmegaQ, dm1 = NA, kappa = NA, M = NA, S = NA, rho = NA, R = NA) {
      rho_bar        <- 1 - rho
      A              <- (self$data$Y - self$data$X %*% B - M %*% t(private$C))^2 + S %*% t(private$C)
      log_det_OmegaQ <- as.numeric(determinant(OmegaQ, logarithm = TRUE)$modulus)
      J <- -.5 * sum(rho_bar %*% (log(2 * pi) - log(dm1)))
      J <- J - .5 * sum(rho_bar * A %*% diag(dm1))
      J <- J + .5 * self$n * log_det_OmegaQ
      J <- J  + .5 * sum(log(S))
      J <- J + sum(rho %*% log(kappa) + rho_bar %*% log(1 - kappa))
      J <- J - sum(rho * log(rho)) - sum(rho_bar*log(rho_bar))
      if (private$sparsity_ > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + .5 * self$n *self$Q - .5 * sum(diag(OmegaQ %*% (crossprod(M) + diag(colSums(S), self$Q, self$Q))))
        J <- J - private$sparsity_ * sum(abs(self$sparsity_weights * OmegaQ))
      }
      J
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated inference ------------------------

    EM_initialize = function() {
      init_model <- normal_diag_zi$new(self$data)
      init_model$optimize()
      B          <- init_model$model_par$B
      ddiag      <- 1/init_model$model_par$dm1
      dm1 <- switch(private$res_covariance,
                      "diagonal"  = 1 / as.vector(ddiag),
                      "spherical" = rep(1/mean(ddiag), self$p))
      OmegaQ     <- t(private$C) %*% diag(dm1) %*% private$C
      kappa      <- init_model$model_par$kappa ## mieux qu'une 0-initialisation ?
      rho        <- init_model$model_par$rho
      G          <- solve(diag(colSums(dm1 * private$C), self$Q, self$Q) + OmegaQ)
      R          <- self$data$Y - self$data$X %*% B
      M          <- R %*% (dm1 * private$C) %*% G
      S          <- matrix(rep(0.1, self$n * self$Q), nrow = self$n)
      list(B = B, dm1 = dm1, OmegaQ = OmegaQ, kappa = kappa, rho = rho, M = M,
           S = S)
    },

    EM_step = function(B, dm1, OmegaQ, kappa, M, S, rho) {
      R   <- self$data$Y - self$data$X %*% B
      ones <- as.vector(rep(1, self$n))

      # E step
      M <- private$zi_NB_fixed_blocks_nlopt_optim_M(M, dm1, OmegaQ, B, rho)
      S <-  1 / sweep((1 - rho) %*% (dm1 * private$C), 2, diag(OmegaQ), "+")
      A <- (R - tcrossprod(M, private$C))^2 + tcrossprod(S, private$C)
      nu <- log(2 * pi) - outer(ones, log(dm1)) + A %*% diag(dm1)
      rho <- 1 / (1 + exp(-.5 * nu) * outer(ones, (1 - kappa) / kappa))
      rho <- check_one_boundary(check_zero_boundary(self$zeros * rho))

      # M step
      B <- private$zi_NB_fixed_blocks_nlopt_optim_B(B, dm1, OmegaQ, M, rho)
      dm1  <- switch(private$res_covariance,
        "diagonal"  = colSums(1 - rho) / colSums((1 - rho) * A),
        "spherical" = rep(sum(1 - rho) / sum((1 - rho) * A), self$p))
      kappa <- colMeans(rho)
      OmegaQ <- private$get_OmegaQ(crossprod(M)/self$n + diag(colMeans(S), self$Q, self$Q))

      list(B = B, dm1 = dm1, OmegaQ = OmegaQ,  kappa = kappa, rho = rho, M = M,
           S = S)
    },

    zi_NB_fixed_blocks_obj_grad_M = function(M_vec, dm1_1mrho, OmegaQ, R) {
      M    <- matrix(M_vec, nrow = self$n, ncol = self$Q)
      RmMC <- R - M %*% t(private$C)
      MO   <- M %*% OmegaQ

      obj  <- -.5 * sum(dm1_1mrho * RmMC^2) - .5 * sum(MO * M)
      grad <- (dm1_1mrho * RmMC) %*% private$C - MO

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_NB_fixed_blocks_nlopt_optim_M = function(M0, dm1, OmegaQ, B, rho) {
      res <- nloptr::nloptr(
        x0 = as.vector(M0),
        eval_f = private$zi_NB_fixed_blocks_obj_grad_M,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        dm1_1mrho = t(dm1 * t(1 - rho)),
        OmegaQ    = OmegaQ,
        R         = self$data$Y - self$data$X %*% B
      )
      newM <- matrix(res$solution, nrow = self$n, ncol = self$Q)
      newM
    },

    zi_NB_fixed_blocks_obj_grad_B = function(B_vec, dm1_1mrho, MC) {
      R   <- self$data$Y - self$data$X %*% matrix(B_vec, nrow = self$d, ncol = self$p)
      RmMC <- R - MC

      obj  <- -.5 * sum(dm1_1mrho * RmMC^2)
      grad <- crossprod(self$data$X, dm1_1mrho * RmMC)

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_NB_fixed_blocks_nlopt_optim_B = function(B0, dm1, OmegaQ, M, rho) {
      res <- nloptr::nloptr(
        x0 = as.vector(B0),
        eval_f = private$zi_NB_fixed_blocks_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        dm1_1mrho = t(dm1 * t(1 - rho)),
        MC = tcrossprod(M, private$C)
      )
      newB <- matrix(res$solution, nrow = self$d, ncol = self$p)
      newB
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference -------------------------
    get_heuristic_parameters = function(){
      model  <- normal_diag_zi$new(self$data) ; model$optimize()
      B      <- model$model_par$B ; kappa <- model$model_par$kappa
      rho    <- model$model_par$rho
      R      <- self$data$Y - self$data$X %*% B ; R[rho > 0.7] <- 0
      Sigma  <- (t(R) %*% R) / model$n
      SigmaQ <- private$heuristic_SigmaQ_from_Sigma(Sigma)
      OmegaQ <- private$get_OmegaQ(SigmaQ)
      list(B = B, OmegaQ = OmegaQ, rho = rho, kappa = kappa)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function() super$nb_param + self$p, # adding kappa
    #' @field var_par a list with variational parameters
    var_par  = function() list(M = private$M, S = private$S, rho = private$rho),
    #' @field model_par a list with model parameters: B (covariates), dm1 (species variance), OmegaQ (groups precision matrix), kappa (zero-inflation probabilities)
    model_par  = function() {
      par       <- super$model_par
      par$kappa <- private$kappa
      par
    },
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      if (!private$approx){
        res <- 0.5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * sum(log(private$S))
        res <- res - sum(private$rho * log(private$rho) + (1 - private$rho) * log(1 - private$rho))
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model Y values predicted by the model
    fitted = function(){
      if (private$approx) {
        res <- self$data$X %*% private$B ; res[private$rho > 0.7] <- 0
      } else {
        res <- (1 - private$rho) * (self$data$X %*% private$B + private$M %*% t(private$C))
      }
      res
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function(value)
    {paste("zero-inflated", private$res_covariance, "normal-block model with fixed blocks")}
  )
)

