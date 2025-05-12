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
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    compute_loglik  = function(B, OmegaQ, dm1 = NA, kappa = NA, gamma = NA, mu = NA) {
      log_det_OmegaQ <- as.numeric(determinant(OmegaQ, logarithm = TRUE)$modulus)
      log_det_Gamma  <- gamma %>%
        map(determinant, logarithm = TRUE) %>%
        map("modulus") %>% map(as.numeric)

      J <- -.5 * self$data$npY * log(2 * pi * exp(1))
      J <- J + .5 * sum(self$data$nY * log(dm1)) + .5 * self$n * log_det_OmegaQ
      J <- J + .5 * sum(unlist(log_det_Gamma))
      J <- J + sum(self$data$zeros %*% log(kappa)) + sum(self$data$zeros_bar %*% log(1 - kappa))
      J <- J - self$data$n * sum(xlogx(kappa) - xlogx(1 - kappa))

      if (private$sparsity_ > 0) {
        gamma_bar <- reduce(gamma, `+`)
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + .5 * self$n * self$Q - .5 * sum(diag(OmegaQ %*% (gamma_bar + crossprod(mu))))
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
      kappa      <- init_model$model_par$kappa
      mu    <- matrix(0, self$n, self$Q)
      gamma <- rep(list(diag(1, self$Q, self$Q)), self$data$n)
      list(B = B, dm1 = dm1, OmegaQ = OmegaQ, kappa = kappa, gamma = gamma, mu = mu)
    },

    EM_step = function(B, dm1, OmegaQ, kappa, gamma = gamma, mu = mu) {

      ## Auxiliary variables
      R <- self$data$Y - self$data$X %*% B
      dm1_mat <- matrix(dm1, self$data$n, self$data$p, byrow = TRUE) * self$data$zeros_bar
      dm1C    <- dm1_mat %*% private$C

      # E step
      gamma <- apply(dm1C, 1, function(dm1C_) {
        solve(OmegaQ + diag(dm1C_, self$Q, self$Q))}, simplify = FALSE)
      Rdm1C <- (R * dm1_mat) %*% private$C
      mu    <- t(sapply(1:length(gamma), function(i) Rdm1C[i, ] %*% gamma[[i]]))

      # M step
      B <- private$zi_NB_fixed_blocks_optim_B(B, dm1_mat, mu)
      RmmuC <- R - mu %*% t(private$C)
      CgC   <- t(sapply(gamma, function(gamma_) diag(gamma_)[self$clustering]))
      A     <- RmmuC^2  + CgC

      dm1  <- switch(private$res_covariance,
        "diagonal"  = self$data$nY / colSums(self$data$zeros_bar * A),
        "spherical" = rep(self$data$npY / sum(self$data$zeros_bar * A), self$p))
      OmegaQ <- private$get_OmegaQ((crossprod(mu) + reduce(gamma, `+`))/self$n)

      list(B = B, dm1 = dm1, OmegaQ = OmegaQ,  kappa = kappa, gamma = gamma, mu = mu)
    },

    zi_NB_fixed_blocks_obj_grad_B = function(B_vec, dm1_mat, muC) {
      RmmuC <- self$data$Y - self$data$X %*% matrix(B_vec, nrow = self$d, ncol = self$p) - muC
      grad <- crossprod(self$data$X, dm1_mat * RmmuC)
      obj <- -.5 * sum(dm1_mat * RmmuC^2)
      res <- list("objective" = -obj, "gradient"  = -grad)
      res
    },

    zi_NB_fixed_blocks_optim_B = function(B0, dm1_mat, mu) {
      muC <- mu %*% t(private$C)
      res <- nloptr::nloptr(
        x0 = as.vector(B0),
        eval_f = private$zi_NB_fixed_blocks_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          maxeval = 100
        ),
        dm1_mat = dm1_mat,
        muC = muC
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
      Sigma  <- crossprod(R) / model$n
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
      if (!private$approx) {
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
        res <- self$data$X %*% private$B + private$mu %*% t(private$C)
        res <- sweep(res, MARGIN = 2, STATS = 1 - private$kappa, FUN = "*")
      }
      res
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function(value)
    {paste("zero-inflated", private$res_covariance, "normal-block model with fixed blocks")}
  )
)

