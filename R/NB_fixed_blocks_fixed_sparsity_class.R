## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_blocks_fixed_sparsity ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param data object of normal_data class, with responses and design matrix
#' @param C clustering matrix C_jq = 1 if species j belongs to cluster q
#' @param penalty to apply on variance matrix when calling GLASSO
#' @param control structured list of more specific parameters, to generate with NB_control
NB_fixed_blocks_fixed_sparsity <- R6::R6Class(
  classname = "NB_fixed_blocks_fixed_sparsity",
  inherit   = NB_fixed_sparsity,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NB_fixed_blocks_fixed_sparsity`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param C clustering matrix C_jq = 1 if species j belongs to cluster q
    #' @param penalty to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NB_fixed_blocks_fixed_sparsity`] object
    initialize = function(data, C, penalty = 0,
                          control = NB_control()) {
      if (!is.matrix(C)) stop("C must be a matrix.")
      if (min(colSums(C)) < 1) stop("There cannot be empty clusters.")
      super$initialize(data, ncol(C), penalty, control = control)
      private$C     <- C
      if(control$inference_method == "integrated"){
        private$mu    <- matrix(0, self$n, self$Q)
        private$gamma <- diag(1, self$Q, self$Q)
      }
    },

    #' @description
    #' Update a [`NB_fixed_blocks_fixed_sparsity`] object
    #' @param B regression matrix
    #' @param OmegaQ groups inverse variance matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param gamma  variance of  posterior distribution of W
    #' @param mu  mean for posterior distribution of W
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB_fixed_blocks_fixed_sparsity`] object
    update = function(B = NA, OmegaQ = NA,  dm1 = NA, gamma = NA, mu = NA,
                      ll_list = NA) {
      super$update(B, OmegaQ, dm1, ll_list)
      if (!anyNA(gamma)) private$gamma <- gamma
      if (!anyNA(mu))   private$mu   <- mu
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    gamma   = NA, # variance of  posterior distribution of latent variable
    mu      = NA, #  mean for posterior distribution of latent variable

    compute_loglik  = function(B, OmegaQ, dm1 = NA, gamma = NA, mu = NA) {
      if(self$inference_method == "integrated"){
        log_det_OmegaQ <- as.numeric(determinant(OmegaQ, logarithm = TRUE)$modulus)
        log_det_gamma  <- as.numeric(determinant(gamma, logarithm = TRUE)$modulus)

        J <- -.5 * self$n * self$p * log(2 * pi * exp(1))
        J <- J + .5 * self$n * sum(log(dm1)) + .5 * self$n * log_det_OmegaQ
        J <- J + .5 * self$n * log_det_gamma
        if (self$penalty > 0) {
          ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
          J <- J + self$n *self$Q / 2 - .5 * sum(diag(OmegaQ %*% (self$n * gamma + t(mu) %*% mu)))
          J <- J - self$penalty * sum(abs(self$sparsity_weights * OmegaQ))
        }
      }else{
        J <- private$heuristic_loglik(B, OmegaQ)
      }
      J
    },

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      SigmaQ    <- private$heuristic_SigmaQ_from_Sigma(reg_res$Sigma)
      OmegaQ    <- private$get_Omega(SigmaQ)
      list(B = reg_res$B, OmegaQ = OmegaQ)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par  = function() list(gamma = private$gamma, mu = private$mu),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      if(self$inference_method == "integrated"){
        log_det_Gamma <- as.numeric(determinant(private$gamma)$modulus)
        ent <- .5 * self$n * self$Q * log(2 * pi* exp(1)) + .5 * self$n * log_det_Gamma
        ent
      }else{NA}
    },
    #' @field fitted Y values predicted by the model
    fitted = function(){
      if(self$inference_method == "integrated"){self$data$X %*% private$B + private$mu %*% t(private$C)
      }else{self$data$X %*% private$B }
    }
  )
)


#' R6 class for normal-block model with fixed clusters and diagonal residual covariance
#' @param data object of normal_data class, with responses and design matrix
#' @param C clustering matrix C_jq = 1 if species j belongs to cluster q
#' @param penalty to apply on variance matrix when calling GLASSO
#' @param control structured list of more specific parameters, to generate with NB_control
NB_fixed_blocks_fixed_sparsity_diagonal <- R6::R6Class(
  classname = "NB_fixed_blocks_fixed_sparsity_diagonal",
  inherit = NB_fixed_blocks_fixed_sparsity,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    EM_initialize = function() {
      if(any(sapply(self$model_par, function(x) any(is.na(x))))){
        B      <- self$data$XtXm1 %*% t(self$data$X) %*% self$data$Y
        dm1    <- as.vector(1/colMeans((self$data$Y - self$data$X %*% B)^2))
        OmegaQ <- diag(colSums(dm1 * private$C), self$Q, self$Q)
        list(B = B, dm1 = dm1, OmegaQ = OmegaQ, gamma = private$gamma, mu = private$mu)
      }else{
        list(B = private$B, OmegaQ = private$OmegaQ, dm1 = private$dm1,
             gamma = private$gamma, mu = private$mu)
      }
    },

    EM_step = function(B, dm1, OmegaQ, gamma, mu) {
      ## E step
      gamma <- solve(OmegaQ + diag(colSums(dm1 * private$C), self$Q, self$Q))
      mu    <- (self$data$Y - self$data$X %*% B) %*% (dm1 * private$C) %*% gamma

      ## M step
      YmmuCT <- self$data$Y - mu %*% t(private$C)
      B      <- self$data$XtXm1 %*% crossprod(self$data$X, YmmuCT)
      ddiag  <- colMeans((YmmuCT - self$data$X %*% B)^2) + private$C %*% diag(gamma)
      dm1    <- 1 / as.vector(ddiag)
      OmegaQ <- private$get_Omega(crossprod(mu)/self$n + gamma)
      list(B = B, OmegaQ = OmegaQ, dm1 = dm1, gamma = gamma, mu = mu)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function(){
      if(self$inference_method == "integrated"){
        as.integer(super$nb_param + self$p)
      }else{as.integer(super$nb_param)}
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function(value) {"diagonal normal-block model with fixed blocks"}
  )
)

#' R6 class for normal-block model with fixed clusters and spherical residual covariance
#' @param data object of normal_data class, with responses and design matrix
#' @param C clustering matrix C_jq = 1 if species j belongs to cluster q
#' @param penalty to apply on variance matrix when calling GLASSO
#' @param control structured list of more specific parameters, to generate with NB_control
NB_fixed_blocks_fixed_sparsity_spherical <- R6::R6Class(
  classname = "NB_fixed_blocks_fixed_sparsity_spherical",
  inherit = NB_fixed_blocks_fixed_sparsity,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    EM_initialize = function() {
      B      <- self$data$XtXm1 %*% t(self$data$X) %*% self$data$Y
      dm1    <- rep(1/mean((self$data$Y - self$data$X %*% B)^2), self$p)
      OmegaQ <- diag(colSums(dm1 * private$C), self$Q, self$Q)
      list(B = B, dm1 = dm1, OmegaQ = OmegaQ, gamma = private$gamma, mu = private$mu)
    },

    EM_step = function(B, dm1, OmegaQ, gamma, mu) {

      ## E step
      gamma <- solve(OmegaQ + dm1[1] * diag(colSums(private$C), self$Q, self$Q))
      mu    <- (self$data$Y - self$data$X %*% B) %*% private$C %*% gamma * dm1[1]

      ## M step
      YmmuCT <- self$data$Y - mu %*% t(private$C)
      B      <- self$data$XtXm1 %*% crossprod(self$data$X, YmmuCT)
      sigma2 <- mean((YmmuCT - self$data$X %*% B)^2) + sum(colMeans(private$C) * diag(gamma))
      OmegaQ <- private$get_Omega(crossprod(mu)/self$n + gamma)

      list(B = B, dm1 = rep(1/sigma2, self$p), OmegaQ = OmegaQ, gamma = gamma, mu = mu)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function(){
      if(self$inference_method == "integrated"){
        as.integer(super$nb_param - self$p + 1)
      }else{as.integer(super$nb_param)}
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function(value) {"spherical normal-block model with fixed blocks"}
  )
)
