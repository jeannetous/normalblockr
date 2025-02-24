## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS normal_zi_fixed_sparsity #####################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param penalty to apply on variance matrix when calling GLASSO
#' @param control structured list of more specific parameters, to generate with normal_control
normal_zi_fixed_sparsity <- R6::R6Class(
  classname = "normal_zi_fixed_sparsity",
  inherit   = normal_fixed_sparsity,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @field zeros indicator matrix of zeros in Y
    zeros = NULL,

    #' @description Create a new [`normal_zi_fixed_sparsity`] object.
    #' @param data contains the matrix of responses (Y) and the design matrix (X).
    #' @param penalty to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with normal_control
    #' @return A new [`normal_zi_fixed_sparsity`] object
    initialize = function(data,  penalty = 0, control = normal_control()) {
      super$initialize(data,  penalty, control)
      self$zeros <- 1 * (data$Y == 0)
    },

    #' @description
    #' Update a [`normal_zi_fixed_sparsity`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of inverse variance matrix
    #' @param kappa vector of zero-inflation probabilities
    #' @param rho posterior probabilities of zero-inflation
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`normal_zi_fixed_sparsity`] object
    update = function(B = NA, dm1 = NA, kappa = NA, rho = NA, ll_list=NA) {
      super$update(B=B, dm1=dm1, ll_list=ll_list)
      if (!anyNA(kappa))      private$kappa   <- kappa
      if (!anyNA(rho))        private$rho     <- rho
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    kappa     = NA,   # vector of zero-inflation probabilities
    rho       = NA,   # posterior probabilities of zero-inflation

    compute_loglik  = function(B, dm1, rho, kappa) {
      rho_bar <- 1 - rho
      J <- - .5 * sum(rho_bar %*% (log(2 * pi) - log(dm1)))
      J <- J - .5 * sum(rho_bar * t(t(self$data$Y - self$data$X %*% B)^2 * dm1))
      J <- J + sum(rho %*% log(kappa)) + sum(rho_bar %*% log(1 - kappa))
      J <- J - sum(rho * log(rho)) - sum(rho_bar*log(rho_bar))
      J
    },

    EM_initialize = function() {
      rho   <- check_one_boundary(check_zero_boundary(self$zeros))
      kappa <- colMeans(rho)
      B     <- self$data$XtXm1 %*% crossprod(self$data$X, self$data$Y)
      dm1   <- 1 / check_one_boundary(check_zero_boundary(diag(cov(self$data$Y - self$data$X %*% B))))
      list(B = B, dm1 = dm1, kappa = kappa, rho = rho)
    },

    EM_step = function(B, dm1, kappa, rho) {

      ## E step
      rho <- 1/(1 + outer(rep(1, self$n), (1 - kappa) / kappa) * dnorm(0, self$data$X %*% B, sqrt(1/dm1)))
      rho <- check_one_boundary(check_zero_boundary(self$zeros * rho))

      ## M step
      B     <- private$normal_zi_optim_B(B, dm1, rho)
      dm1   <- colSums(1 - rho) / colSums((1 - rho) * (self$data$Y - self$data$X %*% B)^2)
      kappa <- colMeans(rho)

      list(B = B, dm1 = dm1, kappa = kappa, rho = rho)
    },

    normal_zi_obj_grad_B = function(B_vec, dm1_1mrho) {
      R <- self$data$Y - self$data$X %*% matrix(B_vec, nrow = self$d, ncol = self$p)
      grad <- crossprod(self$data$X, dm1_1mrho * R)
      obj <- - .5 * sum(dm1_1mrho * R^2)
      res <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    normal_zi_optim_B = function(B0, dm1, rho) {
      dm1_1mrho <- t(dm1 * t(1 - rho))
      res <- nloptr::nloptr(
        x0 = as.vector(B0),
        eval_f = private$normal_zi_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        dm1_1mrho = dm1_1mrho
      )
      newB <- matrix(res$solution, nrow = self$d, ncol = self$p)
      newB
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function() as.integer(super$nb_param + self$p),
    #' @field model_par a list with the matrices of the model parameters:
    #' B (regression coefficients), dm1 (species variance),
    #' kappa (zero-inflation probas)), rho (zero-inflation posterior proba)
    model_par  = function() {
      par <- super$model_par
      par$kappa <- private$kappa
      par$rho <- private$rho
      par
    },
    #' @field fitted Y values predicted by the model Y values predicted by the model
    fitted = function() (1 - private$rho) * (self$X %*% private$B),
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value){"zero-inflated diagonal normal model"}
  )
)
