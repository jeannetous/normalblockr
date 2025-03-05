## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_zi_fixed_Q ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param data object of normal_data class, with responses and design matrix
#' @param Q number of clusters
#' @param penalty to apply on variance matrix when calling GLASSO
#' @param control structured list of more specific parameters, to generate with NB_control
NB_zi_fixed_Q <- R6::R6Class(
  classname = "NB_zi_fixed_Q",
  inherit   = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field clustering_init model initial clustering
    clustering_init = NULL,
    #' @field fixed_tau whether tau should be fixed at clustering_init during optimization, useful for stability selection
    fixed_tau = NULL,
    #' @field zeros indicator matrix of zeros in Y
    zeros = NULL,

    #' @description Create a new [`NB_zi_fixed_Q`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param C block matrix C_jq = 1 if species j belongs to block q
    #' @param control structured list of more specific parameters
    #' @return A new [`NB_zi_fixed_Q`] object
    initialize = function(data, Q, penalty = 0, control = NB_control()) {
      if (Q > ncol(data$Y)) stop("There cannot be more blocks than there are entities to cluster.")
      self$fixed_tau  <- control$fixed_tau
      clustering_init <- control$clustering_init
      super$initialize(data, Q, penalty = penalty, control = control)
      self$zeros <- 1 * (data$Y == 0)
      if (!is.null(clustering_init)) {
        if (!is.vector(clustering_init) & !is.matrix(clustering_init)) stop("Labels must be encoded in list of labels or indicator matrix")
        if (is.vector(clustering_init)) {
          if (any(clustering_init < 1 | clustering_init > Q))
            stop("Cluster labels must be between 1 and Q")
          if (length(clustering_init) != ncol(Y))
            stop("Cluster labels must match the number of Y's columns")
          if (length(unique(clustering_init)) != Q)
            stop("The number of clusters in the initial clustering must be equal to Q.")
        } else {
          if (nrow(clustering_init) != ncol(Y))
            stop("Cluster-indicating matrix must have as many rows as Y has columns")
          if (ncol(clustering_init) != Q)
            stop("Cluster-indicating matrix must have Q columns")
          if ((min(colSums(clustering_init)) < 1) & !self$fixed_tau)
            stop("The number of clusters in the initial clustering must be equal to Q.")
        }
      }
      self$clustering_init <- clustering_init
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS -------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    compute_loglik  = function(B, dm1, OmegaQ, alpha, kappa, M, S, tau, rho) {
      log_det_OmegaQ <- as.numeric(determinant(OmegaQ, logarithm = TRUE)$modulus)
      rho_bar <- 1 - rho
      R <- self$data$Y - self$data$X %*% B
      A <- R^2 - 2 * R * tcrossprod(M,tau) + tcrossprod(M^2 + S, tau)
      J <- -.5 * sum(rho_bar %*% (log(2 * pi) - log(dm1)))
      J <- J - .5 * sum(rho_bar * A %*% diag(dm1))
      J <- J + .5 * self$n * log_det_OmegaQ + sum(tau %*% log(alpha))
      J <- J + sum(rho %*% log(kappa) + rho_bar %*% log(1 - kappa))
      J <- J - sum(rho * log(rho)) - sum(rho_bar*log(rho_bar))
      J <- J  + .5 * sum(log(S)) - sum(tau * log(tau))
      if (self$penalty > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat
        J <- J + .5 * self$n *self$Q - .5 * sum(diag(OmegaQ %*% (crossprod(M) + diag(colSums(S), self$Q, self$Q))))
        J <- J - self$penalty * sum(abs(self$sparsity_weights * OmegaQ))
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
      R          <- self$data$Y - self$data$X %*% B
      R[self$data$Y == 0]  <- 0 # improve final value of objective
      if(is.null(self$clustering_init)){
        clustering <- kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster
        if(length(unique(clustering)) < self$Q){
          # We try to ensure the optimization does not start with an empty cluster
          clustering <- cutree( ClustOfVar::hclustvar(t(R)), Q)
        }
        tau     <- as_indicator(clustering)
        if(min(colSums(tau)) < 0.5) warning("Initialization failed to place elements in each cluster")
      }else{
        if(is.vector(self$clustering_init)){tau <- as_indicator(self$clustering_init)
        }else{ tau <- self$clustering_init}
      }
      tau     <- check_one_boundary(check_zero_boundary(tau))
      alpha      <- colMeans(tau)
      OmegaQ     <- t(tau) %*% diag(dm1) %*% tau
      kappa      <- init_model$model_par$kappa ## mieux qu'une 0-initialisation ?
      rho        <- init_model$model_par$rho
      G          <- solve(diag(colSums(dm1 * tau), self$Q, self$Q) + OmegaQ)
      M          <- R %*% (dm1 * tau) %*% G
      S          <- matrix(rep(0.1, self$n * self$Q), nrow = self$n)
      list(B = B, dm1 = dm1, OmegaQ = OmegaQ, alpha = alpha, kappa = kappa,
           M = M, S = S, tau = tau, rho = rho)
    },

    EM_step = function(B, dm1, OmegaQ, alpha, kappa, M, S, tau, rho) {
      R <- self$data$Y - self$data$X %*% B
      ones <- as.vector(rep(1, self$n))

      # E step
      M <- private$zi_NB_fixed_Q_nlopt_optim_M(M, B, dm1, OmegaQ, tau, rho)
      S <-  1 / sweep((1 - rho) %*% (dm1 * tau), 2, diag(OmegaQ), "+")
      if (self$Q > 1 & !self$fixed_tau) {
        eta <- -.5 * dm1 * t(1 - rho) %*% (M^2 + S)
        eta <- eta + dm1 * t((1 - rho) * R) %*% M  + outer(rep(1, self$p), log(alpha)) - 1
        tau <- t(check_zero_boundary(check_one_boundary(apply(eta, 1, softmax))))
      }
      A <- R^2 - 2 * R * tcrossprod(M,tau) + tcrossprod(M^2 + S, tau)
      nu <- log(2 * pi) - outer(ones, log(dm1)) + A %*% diag(dm1)
      rho <- 1 / (1 + exp(-.5 * nu) * outer(ones, (1 - kappa) / kappa))
      rho <- check_one_boundary(check_zero_boundary(self$zeros * rho))

      # M step
      B   <- private$zi_NB_fixed_Q_nlopt_optim_B(B, dm1, OmegaQ, M, tau, rho)
      dm1  <- switch(private$res_covariance,
                     "diagonal"  = colSums(1 - rho) / colSums((1 - rho) * A),
                     "spherical" = rep(sum(1 - rho) / sum((1 - rho) * A), self$p))
      alpha <- colMeans(tau)
      kappa <- colMeans(rho)
      OmegaQ <- private$get_OmegaQ(crossprod(M)/self$n + diag(colMeans(S), self$Q, self$Q))

      list(B = B, dm1 = dm1, alpha = alpha, OmegaQ = OmegaQ, kappa = kappa,
           M = M, S = S, tau = tau, rho = rho)
    },

    zi_NB_fixed_Q_obj_grad_M = function(M_vec, R, dm1T, OmegaQ, rho) {
      M    <- matrix(M_vec, nrow = self$n, ncol = self$Q)
      MO   <- M %*% OmegaQ

      grad <- ((1 - rho) * R) %*% dm1T - ((1 - rho) %*% dm1T) * M - MO
      obj  <- sum((1 - rho) * (R * (M %*% t(dm1T)) - .5 * M^2 %*% t(dm1T))) - .5 * sum(MO * M)

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_NB_fixed_Q_nlopt_optim_M = function(M0, B, dm1, OmegaQ, tau, rho) {
      M0_vec <- as.vector(M0)
      res <- nloptr::nloptr(
        x0 = M0_vec,
        eval_f = private$zi_NB_fixed_Q_obj_grad_M,
        opts = list(
          algorithm = "NLOPT_LD_MMA",
          xtol_rel = 1e-6,
          maxeval = 1000
        ),
        R = self$data$Y - self$data$X %*% B,
        dm1T    = dm1 * tau,
        OmegaQ = OmegaQ,
        rho    = rho
      )
      newM <- matrix(res$solution, nrow = self$n, ncol = self$Q)
      newM
    },

    zi_NB_fixed_Q_obj_grad_B = function(B_vec, dm1_1mrho, Mtau) {
      R    <- self$data$Y - self$data$X %*% matrix(B_vec, nrow = self$d, ncol = self$p)
      grad <- crossprod(self$data$X, dm1_1mrho * (R - Mtau))
      obj  <- -.5 * sum(dm1_1mrho * (R^2 - 2 * R * Mtau))

      res  <- list("objective" = - obj, "gradient"  = - grad)
      res
    },

    zi_NB_fixed_Q_nlopt_optim_B = function(B0, dm1, OmegaQ, M, tau, rho) {
      dm1_1mrho <- t(dm1 * t(1 - rho))
      res <- nloptr::nloptr(
        x0 = as.vector(B0),
        eval_f = private$zi_NB_fixed_Q_obj_grad_B,
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

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference -----------------------
    get_heuristic_parameters = function() {
      model <- normal_diag_zi$new(self$data) ; model$optimize()
      B      <- model$model_par$B ; kappa <- model$model_par$kappa
      rho    <- model$model_par$rho
      R      <- self$data$Y - self$data$X %*% B ; R[self$rho > 0.7] <- 0
      Sigma  <- (t(R) %*% R) / model$n
      # C and tau play the same role here, but they're used in different
      # parts of the model: tau is used in $clustering bc we're in a fixed_Q context
      # and C is used in generic heuristic methods
      private$C   <- private$clustering_approx(R)
      private$tau <- private$clustering_approx(R)
      SigmaQ <- private$heuristic_SigmaQ_from_Sigma(Sigma)
      OmegaQ <- private$get_OmegaQ(SigmaQ)
      list("B" = B, "OmegaQ" = OmegaQ, rho = rho)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function() super$nb_param + self$p + self$Q - 1, # adding kappa and alpha
    #' @field var_par a list with variational parameters
    var_par  = function() {list(M = private$M, S = private$S,
                                rho = private$rho, tau = private$tau)},
    #' @field model_par a list with model parameters: B (covariates), dm1 (species variance), OmegaQ (blocks precision matrix), kappa (zero-inflation probabilities)
    model_par  = function() {
      par       <- super$model_par
      par$kappa <- private$kappa
      par$alpha <- private$alpha
      par
    },
    #' @field clustering a list of labels giving the clustering obtained in the model
    clustering = function() get_clusters(private$tau),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      if (!private$approx){
        res <- 0.5 * self$n * self$Q * log(2 * pi * exp(1)) + .5 * sum(log(private$S))
        res <- res - sum(private$rho * log(private$rho) + (1 - private$rho) * log(1 - private$rho))
        res <- res - sum(xlogx(private$tau))
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model Y values predicted by the model
    fitted = function(){
      if (private$approx) {
        res <- self$data$X %*% private$B ; res[private$rho > 0.7] <- 0
      } else {
        res <- (1 - private$rho) * (self$data$X %*% private$B + private$M %*% t(private$tau))
      }
      res
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function(value)
    {paste("zero-inflated", private$res_covariance, "normal-block model with", self$Q, "unknown blocks")}
  )
)

