## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZINB_fixed_blocks_fixed_sparsity ############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @export
ZINB_fixed_blocks <- R6::R6Class(
  classname = "ZINB_fixed_blocks",
  inherit   = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`ZINB_fixed_blocks_fixed_sparsity`] object.
    #' @param data object of NBData class, with responses and design matrix
    #' @param C clustering matrix C_jq = 1 if species j belongs to cluster q
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`ZINB_fixed_blocks_fixed_sparsity`] object
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
        map("modulus") %>% map(as.numeric) %>% unlist()

      J <- -.5 * self$data$npY * log(2 * pi * exp(1)) + .5 * sum(self$data$nY * log(dm1))
      J <- J + .5 * (self$n * log_det_OmegaQ + sum(log_det_Gamma))
      J <- J +  private$ZI_cond_mean

      if (private$sparsity_ > 0) {
        gamma_bar <- reduce(gamma, `+`)
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + .5 * self$n * self$Q - .5 * sum(diag(OmegaQ %*% (gamma_bar + crossprod(mu))))
        J <- J - private$sparsity_ * sum(abs(self$sparsity_weights * OmegaQ))
      }
      J
    },

    get_heuristic_parameters = function(){
      zi_diag <- private$zi_diag_normal_inference()
      SigmaQ  <- private$heuristic_SigmaQ_from_Sigma(cov(zi_diag$R))
      OmegaQ  <- private$get_OmegaQ(SigmaQ)
      list(B = zi_diag$B, dm1 = zi_diag$dm1, OmegaQ = OmegaQ, kappa = zi_diag$kappa)
    },

    EM_initialize = function() {
      c(private$get_heuristic_parameters(),  list(
        gamma = rep(list(diag(1, self$Q, self$Q)), self$n),
        mu    = matrix(0, self$n, self$Q)
        )
      )
    },

    EM_step = function(B, dm1, OmegaQ, kappa, gamma, mu) {

      ## Auxiliary variables
      R <- self$data$Y - self$data$X %*% B
      dm1_mat <- matrix(dm1, self$n, self$p, byrow = TRUE) * self$data$zeros_bar
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
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par = function(value) list(gamma = private$gamma, mu = private$mu),
    #' @field entropy Entropy of the conditional distribution
    entropy    = function(value) {
      if (!private$approx){
        log_det_Gamma <- private$gamma %>%
          map(determinant, logarithm = TRUE) %>%
          map("modulus") %>% map(as.numeric) %>% unlist()
        res <- .5 * (self$n * self$Q * log(2*pi*exp(1)) + sum(log_det_Gamma))
      } else {res <- NA}
      res
    },
    #' @field nb_param number of parameters in the model
    nb_param = function() super$nb_param + self$p, # adding kappa
    #' @field model_par a list with model parameters: B (covariates), dm1 (species variance), OmegaQ (groups precision matrix), kappa (zero-inflation probabilities)
    model_par  = function() {
      par       <- super$model_par
      par$kappa <- private$kappa
      par
    },
    #' @field fitted Y values predicted by the model
    fitted = function(){
      if (private$approx) {
        res <- self$data$X %*% private$B
      } else {
        res <- self$data$X %*% private$B + private$mu %*% t(private$C)
      }
      res <- res * self$data$zeros_bar
      res
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function(value)
    {paste("zero-inflated", private$res_covariance, "normal-block model with fixed blocks")}
  )
)

