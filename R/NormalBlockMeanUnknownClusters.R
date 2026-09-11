## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanUnknownClusters ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Mean-Block Model with Unknown Clustering
#'
#' R6 class for a Normal-Block-Mean model with a fixed number of clusters
#' (but unknown clustering), inferred by variational EM.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' model <- normal_block(data, blocks = 3, model = "mean")
#' model$clustering
#' @export
NormalBlockMeanUnknownClusters <- R6::R6Class(
  classname = "NormalBlockMeanUnknownClusters",
  inherit   = NormalBlockMeanBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockMeanUnknownClusters`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param q number of clusters
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockMeanUnknownClusters`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    optim_initialize = function() {
      if (private$warm_started) {
        list(B = private$B, Omega = private$Omega, tau = private$C)
      } else {
        private$get_heuristic_parameters()
      }
    },

    get_heuristic_parameters = function(){
      reg_res <- private$multivariate_normal_inference()
      Omega   <- private$omega_from_sigma(reg_res$Sigma)
      ## cluster on the fitted mean trajectory X %*% B_j (n x p): unlike the
      ## raw coefficient B_j (d x 1), it has enough rows for the shared
      ## cor()/cov()-based heuristics (ward2/sbm/spectral) even when d is small
      if (anyNA(private$C))
        private$C <- private$heuristic_clustering(self$data$X %*% reg_res$B)
      tau <- private$C
      ## soften the hard indicator away from 0/1: compute_loglik()'s entropy
      ## term (tau * log(tau)) is initialized on this tau
      tau <- check_one_boundary(check_zero_boundary(tau))
      tau <- tau / rowSums(tau)
      B   <- private$heuristic_cluster_B_from_variable_B(reg_res$B, tau)
      list(B = B, Omega = Omega, tau = tau)
    },


    B_estimator = function(Omega = private$Omega,
                           tau = private$C,
                           Psi = private$Psi){
      return(self$data$XtXm1 %*% self$data$XtY %*% Omega %*% tau %*% solve(Psi))
    },

    Lambda_estimator = function(B = private$B,
                                tau = private$C) {

    M      <- crossprod(B, self$data$XtX) %*% B
    diag_M <- diag(M)

    diag_val1 <- as.vector(tau %*% diag_M)
    term1 <- diag(diag_val1, nrow = self$data$p, ncol = self$data$p)

    diag_val2 <- diag(tau %*% M %*% t(tau))
    term2 <- diag(diag_val2, nrow = self$data$p, ncol = self$data$p)

    return((term1 - term2) / self$data$n)
    },

    Sigma_estimator = function(Omega = private$Omega, B = private$B,
                               tau = private$C, Lambda = private$Lambda){

      R_bar <- self$data$Y - self$data$X %*% B %*% t(tau)

      Sigma <- (crossprod(R_bar) / self$data$n) + Lambda
      return(Sigma)
    },

    Psi_estimator = function(Omega = private$Omega, tau = private$C) {
      diag_Omega <- diag(Omega)
      term1 <- crossprod(tau, Omega) %*% tau

      diag_product <- as.vector(crossprod(tau, diag_Omega))
      term2 <- diag(diag_product, nrow = self$q, ncol = self$q)

      term3 <- crossprod(tau, diag(diag_Omega, self$data$p)) %*% tau

      Psi <- term1 + term2 - term3

      Psi <- Psi + 1e-8 * diag(self$q)

      return(Psi)
    },

    Phi_estimator = function(Omega = private$Omega, tau = private$C,
                             Psi = private$Psi) {
      tau_Omega_tau <- crossprod(tau, Omega) %*% tau
      Phi <- Psi - tau_Omega_tau
      return(Phi)
    },

    alpha_estimator = function(tau = private$C) {
      return(colMeans(tau))
    },

    ## Sequential (Gauss-Seidel) sweep over the rows of tau: each row's
    ## softmax maximizes the ELBO exactly.
    tau_estimator = function(Omega = private$Omega,
                             B     = private$B,
                             alpha = private$alpha,
                             tau = private$C) {
      M      <- crossprod(B, self$data$XtX) %*% B
      diag_M <- diag(M)
      w      <- diag(Omega)
      ZtXB   <- crossprod(self$data$Y, self$data$X %*% B)
      log_alpha <- log(pmax(alpha, 1e-300))

      G <- ZtXB - tau %*% M ## kept in sync with tau, row by row
      for (j in seq_len(self$data$p)) {
        expo  <- log_alpha + as.vector(Omega[j, ] %*% G) +
          w[j] * as.vector(tau[j, ] %*% M) - 0.5 * w[j] * diag_M
        tau_j <- exp(expo - max(expo))
        tau_j <- tau_j / sum(tau_j)
        tau[j, ] <- tau_j
        G[j, ]   <- ZtXB[j, ] - tau_j %*% M
      }
      tau <- check_one_boundary(check_zero_boundary(tau))
      tau <- tau / rowSums(tau)
      return(tau)
    },


    compute_loglik = function(B = private$B, Omega = private$Omega,
                              tau = private$C, alpha = private$alpha,
                              Phi = private$Phi){

      R_bar <- self$data$Y - tcrossprod(self$data$X %*% B, tau)
      M     <- crossprod(B, self$data$XtX) %*% B
      l <- - (self$n * self$p / 2) * log(2 * pi) +
           (self$n / 2) * as.numeric(determinant(Omega, logarithm = TRUE)$modulus) +
           sum(tau * outer(rep(1, self$p), log(alpha))) - sum(tau * log(tau)) -
           0.5 * (sum(R_bar * (R_bar %*% Omega)) + sum(diag(Phi %*% M)))
      ## penalized objective (see the C++ core): $loglik adds the penalty back
      if (self$sparsity > 0) l <- l - self$sparsity * sum(abs(self$sparsity_weights * Omega))

      return(as.numeric(l))
    },

    ## Moment-based fit (NB_control(heuristic = TRUE)); see the known-clusters
    ## class for why this is just the heuristic initialization.
    heuristic_optimize = function(control) {
      private$niter <- 0
      init <- private$get_heuristic_parameters()
      list(B = init$B, Omega = init$Omega, C = init$tau,
           alpha = private$alpha_estimator(init$tau))
    },

    ## Runs the VEM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## NormalBlockMeanUnknownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- NormalBlockMeanUnknownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        B0 = init$B, Omega0 = init$Omega, tau0 = init$tau,
        sparsity = self$sparsity, sparsity_weights = self$sparsity_weights,
        noise_covariance = private$res_covariance,
        fixed_point_niter = control$fixed_point_niter,
        fixed_tau = isTRUE(control$fixed_tau),
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, Omega = res$Omega, C = res$C, alpha = res$alpha,
           Psi = res$Psi, Phi = res$Phi, Lambda = res$Lambda, ll_list = res$objective)
    },

    ## Reference R implementation of the same recursion, kept since the C++
    ## port is being validated against it (test-cpp-normal-block-mean.R).
    EM_optimize_R = function(control){
      init_params <- private$optim_initialize()
      B           <- init_params$B
      Omega       <- init_params$Omega
      tau         <- init_params$tau
      alpha       <- private$alpha_estimator(tau)
      Psi         <- private$Psi_estimator(Omega, tau)
      Phi         <- private$Phi_estimator(Omega, tau, Psi)
      Lambda      <- private$Lambda_estimator(B, tau)
      ll_prev     <- private$compute_loglik(B, Omega, tau, alpha, Phi)
      ll_list     <- c(ll_prev)
      warm        <- NULL

      for(i in 1:control$niter){
        # M-step
        Psi    <- private$Psi_estimator(Omega, tau)
        B      <- private$B_estimator(Omega, tau, Psi)
        Lambda <- private$Lambda_estimator(B, tau)
        Sigma  <- private$Sigma_estimator(Omega, B, tau, Lambda)
        og     <- private$omega_from_sigma_warm(Sigma, warm)
        Omega  <- og$Omega; warm <- og$warm
        alpha  <- private$alpha_estimator(tau)

        # VE-step
        if (!isTRUE(control$fixed_tau)) {
          for (k in 1:control$fixed_point_niter) {
            tau <- private$tau_estimator(Omega, B, alpha, tau)
          }
        }


        # elbo updating (Phi depends on the freshly updated Omega/tau)
        Phi        <- private$Phi_estimator(Omega, tau, private$Psi_estimator(Omega, tau))
        ll_current <- private$compute_loglik(B, Omega, tau, alpha, Phi)
        ll_list     <- c(ll_list, ll_current)
        if(abs(ll_current - ll_prev) < control$threshold){
          break
        }else{ll_prev <- ll_current}
      }
      private$niter <- i
      list(B = B, Omega = Omega, C = tau, alpha = alpha, Psi = Psi, Phi = Phi,
           Lambda = Lambda, ll_list = ll_list)
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field fitted Y values predicted by the model, in Y's original units
    fitted = function(){
      ## mean-block model: mu_i = C B' X_i, i.e. X B C' in matrix form
      private$rescale_to_original(self$data$X %*% private$B %*% t(private$C))
    },
    #' @field var_par a list with the variational parameter: tau (posterior group probabilities)
    var_par    = function() list(tau = private$C),
    #' @field nb_param number of parameters in the model
    nb_param = function() {as.integer(super$nb_param + self$q - 1)}, # adding alpha
    #' @field entropy Entropy of the conditional distribution
    #' The only latent variable is the clustering, so the entropy of the
    #' variational distribution reduces to -sum(tau * log(tau)).
    entropy    = function() -sum(xlogx(private$C)),
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function()
    {paste(private$res_covariance, "normal-block-mean model with", self$q, "unknown blocks")}
  )
)
