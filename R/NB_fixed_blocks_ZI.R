## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  NB_fixed_blocks_ZI #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed groups
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood threshold under which optimization stops
#' @export
NB_fixed_blocks_ZI <- R6::R6Class(
  classname = "NB_fixed_blocks_ZI",
  inherit = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field C the matrix of species groups
    C = NULL,

    #' @description Create a new [`NB_fixed_blocks_ZI`] object.
    #' @param C group matrix C_jq = 1 if species j belongs to group q
    #' @return A new [`NB_fixed_blocks`] object
    initialize = function(Y, X, C, sparsity = 0, niter = 50, threshold = 1e-4) {
      if (!is.matrix(C)) stop("C must be a matrix.")
      self$C <- C
      self$Q <- ncol(C)
      super$initialize(Y, X, sparsity,niter, threshold)
    },

    #' @description
    #' Update a [`NB_fixed_Q`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param omegaQ groups inverse variance matrix
    #' @param alpha vector of groups probabilities
    #' @param tau posterior probabilities for group affectation
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB_fixed_Q`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, kappa = NA,
                      M = NA, S = NA, rho = NA, ll_list=NA) {
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
    M       = NA, # variational mean for posterior distribution of W
    S       = NA, # variational diagonal of variances for posterior distribution of W
    kappa   = NA,   # vector of zero-inflation probabilities
    rho     = NA,   # posterior probabilities of zero-inflation


    compute_loglik  = function(B, dm1, omegaQ, kappa, rho, M, S) {
      R              <- self$Y - self$X %*% B
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      elbo <- -.5 * sum((1 - rho) * log(2 * pi))
      ### deux lignes ci-dessous à récrire pour enlever le diag(dm1)
      elbo <- elbo + .5 * sum((1 - rho) %*% diag(log(dm1)))
      elbo <- elbo - .5 * sum(((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(C))) %*% diag(dm1))
      elbo <- elbo - .5 * self$n * self$Q * log(2 * pi) + .5 * self$n * log_det_omegaQ
      elbo <- elbo - .5 * sum((M %*% omegaQ) * M) - .5 * sum(S %*% diag(omegaQ))
      elbo <- elbo + sum(rho %*% log(kappa) + (1 - rho) %*% log(1 - kappa))
      elbo <- elbo + .5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * sum(log(S))
      elbo <- elbo - sum(rho * log(rho) + (1 - rho) * log(1 - rho))
      elbo
    },

    EM_initialize = function() {
      init_model <- zi_normal$new(self$Y, self$X)
      init_model$optimize()
      B          <- init_model$model_par$B
      dm1        <- init_model$model_par$dm1
      omegaQ     <- t(self$C) %*% diag(init_model$model_par$dm1) %*% self$C
      kappa      <- init_model$model_par$kappa ## mieux qu'une 0-initialisation ?
      rho        <- init_model$model_par$rho
      G          <- solve(diag(colSums(as.vector(dm1) * self$C)) + omegaQ)
      R          <- self$Y - self$X %*% B
      M          <- dm1 * R %*% self$C %*% G
      S          <- matrix(rep(0.1, self$n * self$Q), nrow=self$n)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, kappa = kappa, rho = rho, M = M,
           S = S)
    },

    zi_nb_fixed_blocks_obj_grad_M = function(M_vec, B, dm1, omegaQ, rho, kappa, S) {
      M    <- matrix(M_vec, nrow = self$n, ncol = self$Q)
      R   <- Y - X %*% B
      grad <- ((1 - rho) * R) %*%  diag(dm1) %*% C - ((1 - rho) %*% diag(dm1) %*% C) * M - M %*% omegaQ

      obj  <- - .5 * sum(((1 - rho) * (R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(C))) %*% diag(dm1))
      obj  <- obj - .5 * sum((M %*% omegaQ) * M)

      res  <- list("objective" = obj, "gradient"  = grad)
      res
    },

    zi_nb_fixed_blocks_nlopt_optim_M = function(M0, B, dm1, omegaQ, rho, kappa,
                                                S) {
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
        rho    = rho,
        kappa  = kappa,
        S      = S
      )
      newM <- matrix(res$solution, nrow = self$n, ncol = self$Q)
      newM
    },

    EM_step = function(B, dm1, omegaQ, kappa, rho, M, S) {
      R      <- self$Y - self$X %*% B

      ## E step
      M <-  private$zi_nb_fixed_blocks_nlopt_optim_M(M, B, dm1, omegaQ, rho, kappa, S)
      S <-  1/sweep(((1 - rho) %*% (dm1 * self$C)), 2, diag(omegaQ), "+")

      a   <- R^2 - 2 * R * (M %*% t(self$C)) + (M^2 + S) %*% t(self$C)
      b   <- .5 * as.vector(rep(1, self$n)) %*% t(as.vector(rep(1, self$p)))
      eta <- b * log(2 * pi) + b %*% diag(log(dm1)) ### à récrire
      eta <- eta + .5 * a %*% diag(dm1) + as.vector(rep(1, self$n)) %*% t(log(kappa))
      eta <- eta - as.vector(rep(1, self$n)) %*% t(log(1 - kappa)) # ERREUR CALCUL ETA
      rho =  (exp(eta) / (1 + exp(eta)))

      ## M step
      B       <- private$XtXm1 %*% t(self$X) %*% (self$Y - M %*% t(self$C))  ###### CALCULS FAUX - ceux de NB FB sans ZI
      omegaQ  <- self$n * solve(t(M) %*% M + diag(t(S) %*% as.vector(rep(1, self$n))))
      dm1     <- as.vector((as.vector(rep(1, self$n)) %*% (1 - rho)) / (as.vector(rep(1, self$n)) %*% ((1 - rho) * (R^2 - 2 * R * M %*% t(C) + (M^2 + S) %*% t(C)))))
      kappa   <- check_one_boundary(check_zero_boundary(colMeans(rho)))

      list(B = B, dm1 = dm1, omegaQ = omegaQ,  kappa = kappa, rho = rho, M = M,
           S = S)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par  = function() {list(gamma = private$gamma, mu = private$mu)}
  )
)

