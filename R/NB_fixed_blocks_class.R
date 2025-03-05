## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_blocks                ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param data object of normal_data class, with responses and design matrix
#' @param C clustering matrix C_jq = 1 if species j belongs to cluster q
#' @param sparsity sparsity penalty to apply on variance matrix when calling GLASSO
#' @param control structured list of more specific parameters, to generate with NB_control
NB_fixed_blocks <- R6::R6Class(
  classname = "NB_fixed_blocks",
  inherit   = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NB_fixed_blocks`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param C clustering matrix C_jq = 1 if species j belongs to cluster q
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NB_fixed_blocks`] object
    initialize = function(data, C, sparsity = 0, control = NB_control()) {
      stopifnot("C must be a matrix" = is.matrix(C))
      stopifnot("There cannot be empty clusters" = min(colSums(C)) > 0)
      super$initialize(data, ncol(C), sparsity, control)
      private$C <- C
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    compute_loglik  = function(B, OmegaQ, dm1 = NA, gamma = NA, mu = NA) {
      log_det_OmegaQ <- as.numeric(determinant(OmegaQ, logarithm = TRUE)$modulus)
      log_det_gamma  <- as.numeric(determinant(gamma, logarithm = TRUE)$modulus)

      J <- -.5 * self$n * self$p * log(2 * pi * exp(1))
      J <- J + .5 * self$n * sum(log(dm1)) + .5 * self$n * log_det_OmegaQ
      J <- J + .5 * self$n * log_det_gamma
      if (self$sparsity > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + self$n*self$Q / 2 - .5 * sum(diag(OmegaQ %*% (self$n * gamma + t(mu) %*% mu)))
        J <- J - self$sparsity * sum(abs(self$sparsity_weights * OmegaQ))
      }
      J
    },

    EM_initialize = function() {
      B     <- self$data$XtXm1 %*% self$data$XtY
      ddiag <- colMeans((self$data$Y - self$data$X %*% B)^2)
      dm1   <- switch(private$res_covariance,
        "diagonal"  = 1 / as.vector(ddiag),
        "spherical" = rep(1/mean(ddiag), self$p))
      OmegaQ <- diag(colSums(dm1 * private$C), self$Q, self$Q)
      mu    <- matrix(0, self$n, self$Q)
      gamma <- diag(1, self$Q, self$Q)
      list(B = B, dm1 = dm1, OmegaQ = OmegaQ, gamma = gamma, mu = mu)
    },

    EM_step = function(B, dm1, OmegaQ, gamma, mu) {
      ## E step
      gamma <- solve(OmegaQ + diag(colSums(dm1 * private$C), self$Q, self$Q))
      mu    <- (self$data$Y - self$data$X %*% B) %*% (dm1 * private$C) %*% gamma

      ## M step
      YmmuCT <- self$data$Y - mu %*% t(private$C)
      B      <- self$data$XtXm1 %*% crossprod(self$data$X, YmmuCT)
      ddiag  <- colMeans((YmmuCT - self$data$X %*% B)^2) + private$C %*% diag(gamma)
      dm1  <- switch(private$res_covariance,
        "diagonal"  = 1 / as.vector(ddiag),
        "spherical" = rep(1/mean(ddiag), self$p))
      OmegaQ <- private$get_OmegaQ(crossprod(mu)/self$n + gamma)
      list(B = B, OmegaQ = OmegaQ, dm1 = dm1, gamma = gamma, mu = mu)
    },

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      SigmaQ    <- private$heuristic_SigmaQ_from_Sigma(reg_res$Sigma)
      OmegaQ    <- private$get_OmegaQ(SigmaQ)
      list(B = reg_res$B, OmegaQ = OmegaQ)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par = function(value) list(gamma = private$gamma, mu = private$mu),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function(value) {
      if (!private$approx){
        res <- .5 * self$n * self$Q * log(2 * pi * exp(1)) +
          .5 * self$n * as.numeric(determinant(private$gamma)$modulus)
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model
    fitted = function(value){
      if (private$approx) {
        res <- self$data$X %*% private$B
      } else {
        res <- self$data$X %*% private$B + private$mu %*% t(private$C)
      }
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function(value)
      {paste(private$res_covariance, "normal-block model with fixed blocks")}
  )
)

