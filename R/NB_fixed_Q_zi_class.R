## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  NB_fixed_Q_zi #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed blocks
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param Q number of blocks in the model
#' @param sparsity to add on blocks precision matrix
#' @param verbose telling if information should be printed during optimization
NB_fixed_Q_zi <- R6::R6Class(
  classname = "NB_fixed_Q_zi",
  inherit = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field zeros indicator matrix of zeros in Y
    zeros = NULL,
    #' @field clustering_init model initial clustering
    clustering_init = NULL,

    #' @description Create a new [`NB_fixed_Q_zi`] object.
    #' @param C block matrix C_jq = 1 if species j belongs to block q
    #' @param clustering_init model initial clustering
    #' @return A new [`NB_fixed_blocks`] object
    initialize = function(Y, X, Q, sparsity = 0, clustering_init = NULL) {
      super$initialize(Y = Y, X = X, Q, sparsity = sparsity)
      self$zeros <- 1 * (Y == 0)
      if (!is.null(clustering_init)) {
        if(!is.vector(clustering_init) & !is.matrix(clustering_init)) stop("Labels must be encoded in list of labels or indicator matrix")
        if (is.vector(clustering_init)) {
          if (any(clustering_init < 1 | clustering_init > Q))
            stop("Cluster labels must be between 1 and Q")
          if (length(clustering_init) != ncol(Y))
            stop("Cluster labels must match the number of Y's columns")
        } else {
          if (nrow(clustering_init) != ncol(Y))
            stop("Cluster-indicating matrix must have as many rows as Y has columns")
          if (ncol(clustering_init) != Q)
            stop("Cluster-indicating matrix must have Q columns")
        }
      }
      self$clustering_init <- clustering_init
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
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)
      rho_bar <- 1 - rho
      R <- self$Y - self$X %*% B
      A <- R^2 - 2 * R * tcrossprod(M,tau) + tcrossprod(M^2 + S, tau)
      J <- -.5 * sum(rho_bar %*% (log(2 * pi) - log(dm1)))
      J <- J - .5 * sum(rho_bar * A %*% diag(dm1))
      J <- J + .5 * self$n * log_det_omegaQ + sum(tau %*% log(alpha))
      J <- J + sum(rho %*% log(kappa) + rho_bar %*% log(1 - kappa))
      J <- J - sum(rho * log(rho)) - sum(rho_bar*log(rho_bar))
      J <- J  + .5 * sum(log(S)) - sum(tau * log(tau))
      if (self$sparsity > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + .5 * self$n *self$Q - .5 * sum(diag(omegaQ %*% (crossprod(M) + diag(colSums(S), self$Q, self$Q))))
        J <- J - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
      J
    },

    EM_initialize = function() {
      init_model <- normal_zi$new(self$Y, self$X)
      init_model$optimize()
      B          <- init_model$model_par$B
      dm1        <- init_model$model_par$dm1
      R          <- self$Y - self$X %*% B
      R[self$Y == 0]  <- 0 # improve final value of objective
      # cl         <- kmeans(t(R), self$Q, nstart = 30)$cluster
      # tau        <- check_one_boundary(check_zero_boundary(as_indicator(cl)))
      if(is.null(self$clustering_init)){
        tau     <- as_indicator(kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster)
      }else{
        if(is.vector(self$clustering_init)){tau <- as_indicator(self$clustering_init)
        }else{ tau <- self$clustering_init}
      }
      tau     <- check_one_boundary(check_zero_boundary(tau))
      alpha      <- colMeans(tau)
      omegaQ     <- t(tau) %*% diag(dm1) %*% tau
      kappa      <- init_model$model_par$kappa ## mieux qu'une 0-initialisation ?
      rho        <- init_model$model_par$rho
      G          <- solve(diag(colSums(dm1 * tau), self$Q, self$Q) + omegaQ)
      M          <- R %*% (dm1 * tau) %*% G
      S          <- matrix(rep(0.1, self$n * self$Q), nrow = self$n)
      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, kappa = kappa,
           M = M, S = S, tau = tau, rho = rho)
    },

    zi_nb_fixed_Q_obj_grad_M = function(M_vec, R, dm1T, omegaQ, rho) {
      M    <- matrix(M_vec, nrow = self$n, ncol = self$Q)
      MO   <- M %*% omegaQ

      grad <- ((1 - rho) * R) %*% dm1T - ((1 - rho) %*% dm1T) * M - MO
      obj  <- sum((1 - rho) * (R * (M %*% t(dm1T)) - .5 * M^2 %*% t(dm1T))) - .5 * sum(MO * M)

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_Q_nlopt_optim_M = function(M0, B, dm1, omegaQ, tau, rho) {
      M0_vec <- as.vector(M0)
      res <- nloptr::nloptr(
        x0 = M0_vec,
        eval_f = private$zi_nb_fixed_Q_obj_grad_M,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        R = self$Y - self$X %*% B,
        dm1T    = dm1 * tau,
        omegaQ = omegaQ,
        rho    = rho
      )
      newM <- matrix(res$solution, nrow = self$n, ncol = self$Q)
      newM
    },

    zi_nb_fixed_Q_obj_grad_B = function(B_vec, dm1_1mrho, Mtau) {
      R    <- self$Y - self$X %*% matrix(B_vec, nrow = self$d, ncol = self$p)
      grad <- crossprod(self$X, dm1_1mrho * (R - Mtau))
      obj  <- -.5 * sum(dm1_1mrho * (R^2 - 2 * R * Mtau))

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_nb_fixed_Q_nlopt_optim_B = function(B0, dm1, omegaQ, M, tau, rho) {
      dm1_1mrho <- t(dm1 * t(1 - rho))
      res <- nloptr::nloptr(
        x0 = as.vector(B0),
        eval_f = private$zi_nb_fixed_Q_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        dm1_1mrho = t(dm1 * t(1 - rho)),
        Mtau = tcrossprod(M, tau)
      )
      newB <- matrix(res$solution, nrow = self$d, ncol = self$p)
      newB
    },

    EM_step = function(B, dm1, omegaQ, alpha, kappa, M, S, tau, rho) {
      R <- self$Y - self$X %*% B
      ones <- as.vector(rep(1, self$n))

      # E step
      M <- private$zi_nb_fixed_Q_nlopt_optim_M(M, B, dm1, omegaQ, tau, rho)
      S <-  1 / sweep((1 - rho) %*% (dm1 * tau), 2, diag(omegaQ), "+")
      if (self$Q > 1) {
        eta <- -.5 * dm1 * t(1 - rho) %*% (M^2 + S)
        eta <- eta + dm1 * t((1 - rho) * R) %*% M  + outer(rep(1, self$p), log(alpha)) - 1
        tau <- t(check_zero_boundary(check_one_boundary(apply(eta, 1, softmax))))
      }
      A <- R^2 - 2 * R * tcrossprod(M,tau) + tcrossprod(M^2 + S, tau)
      nu <- log(2 * pi) - outer(ones, log(dm1)) + A %*% diag(dm1)
      rho <- 1 / (1 + exp(-.5 * nu) * outer(ones, (1 - kappa) / kappa))
      rho <- check_one_boundary(check_zero_boundary(self$zeros * rho))

      # M step
      B   <- private$zi_nb_fixed_Q_nlopt_optim_B(B, dm1, omegaQ, M, tau, rho)
      dm1   <- colSums(1 - rho) / colSums((1 - rho) * A)
      alpha <- colMeans(tau)
      kappa <- colMeans(rho)
      omegaQ <- private$get_omegaQ((1 / self$n) * (crossprod(M) + diag(colSums(S), self$Q, self$Q)))

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
      ent <- 0.5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * sum(log(private$S))
      ent <- ent - sum(xlogx(private$tau))
      ent <- ent - sum(private$rho * log(private$rho) + (1 - private$rho) * log(1 - private$rho))
      ent
    },
    #' @field fitted Y values predicted by the model Y values predicted by the model
    fitted = function()(1 - private$rho) *(self$X %*% private$B + tcrossprod(private$M, private$tau))
  )
)
