## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q                                  ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed Q (number of groups)
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param Q number of clusters
#' @param penalty to add on blocks precision matrix for sparsity
#' @param control structured list for specific parameters (including initial clustering proposal)
NB_fixed_Q <- R6::R6Class(
  classname = "NB_fixed_Q",
  inherit = NB,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field clustering_init model initial clustering
    clustering_init = NULL,
    #' @field fixed_tau whether tau should be fixed at clustering_init during optimization, useful for stability selection
    fixed_tau = NULL,

    #' @description Create a new [`NB_fixed_Q`] object.
    #' @param Q required number of groups
    #' @param control structured list for specific parameters
    #' @return A new [`NB_fixed_Q`] object
    initialize = function(data, Q, penalty = 0, control = NB_control()) {
      super$initialize(data, Q, penalty, control)
      stopifnot("There cannot be more blocks than there are entities to cluster" = Q <= ncol(self$data$Y))
      self$fixed_tau <- control$fixed_tau
      clustering_init <- control$clustering_init
      if (!is.null(clustering_init)) {
        if(!is.vector(clustering_init) & !is.matrix(clustering_init)) stop("Labels must be encoded in list of labels or indicator matrix")
        if (is.vector(clustering_init)) {
          if (any(clustering_init < 1 | clustering_init > Q))
            stop("Cluster labels must be between 1 and Q")
          if (length(clustering_init) != ncol(data$Y))
            stop("Cluster labels must match the number of Y's columns")
          if (length(unique(clustering_init)) != Q)
            stop("The number of clusters in the initial clustering must be equal to Q.")
        } else {
          if (nrow(clustering_init) != ncol(data$Y))
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
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    compute_loglik  = function(B, dm1, OmegaQ, alpha, tau, M, S) {
      log_det_OmegaQ <- as.numeric(determinant(OmegaQ, logarithm = TRUE)$modulus)

      J <- -.5 * self$n * self$p * log(2 * pi * exp(1)) + .5 * self$n * sum(log(dm1))
      J <- J  + .5 * self$n * log_det_OmegaQ
      J <- J + sum(tau %*% log(alpha))
      J <- J - sum(xlogx(tau)) + .5 * self$n * sum(log(S))

      if (self$penalty > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat and simplifies
        J <- J + self$n*self$Q / 2 - .5 * sum(diag(OmegaQ %*% (crossprod(M) + self$n * diag(S, self$Q, self$Q))))
        J <- J - self$penalty * sum(abs(self$penalty_weights * OmegaQ))
      }
      J
    },

    EM_initialize = function() {
      B <- self$data$XtXm1 %*% self$data$XtY
      R <- self$data$Y - self$data$X %*% B
      if (is.null(self$clustering_init)){
        clustering <- kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster
        if(length(unique(clustering)) < self$Q){
          # We try to ensure the optimization does not start with an empty cluster
          clustering <- cutree(ClustOfVar::hclustvar(t(R)), Q)
        }
        tau <- as_indicator(clustering)
        if (min(colSums(tau)) < 1) warning("Initialization failed to place elements in each cluster")
      } else {
        if (is.vector(self$clustering_init)){tau <- as_indicator(self$clustering_init)
        } else { tau <- self$clustering_init}
      }
      # self$clustering_init <- get_clusters(tau)
      tau     <- check_one_boundary(check_zero_boundary(tau))
      alpha   <- colMeans(tau)
      S       <- rep(0.1, self$Q)
      M       <- matrix(rep(0, self$n * self$Q), nrow = self$n)
      ddiag <- colMeans((self$data$Y - self$data$X %*% B)^2)
      dm1   <- switch(private$res_covariance,
                      "diagonal"  = 1 / as.vector(ddiag),
                      "spherical" = rep(1/mean(ddiag), self$p))
      OmegaQ  <- diag(rep(1, self$Q), self$Q, self$Q)
      self$clustering_init <- get_clusters(tau)
      list(B = B, OmegaQ = OmegaQ, dm1 = dm1,  alpha = alpha, tau = tau, M = M, S = S)
    },

    EM_step = function(B, OmegaQ, dm1, alpha, tau, M, S) {

      ## Auxiliary variables
      R     <- self$data$Y - self$data$X %*% B
      Gamma <- solve(OmegaQ + diag(colSums(dm1 * tau), self$Q, self$Q))

      # E step
      M <- R %*% (dm1 * tau) %*% Gamma
      S <- diag(Gamma)

      if (self$Q > 1 & !self$fixed_tau) {
        eta <- dm1 * (t(R) %*% M) - .5 * outer(dm1,  colSums(M^2) + self$n * S)
        eta <- eta + outer(rep(1, self$p), log(alpha)) - 1
        tau <- t(check_zero_boundary(check_one_boundary(apply(eta, 1, softmax))))
      }

      # M step
      MtauT <- M %*% t(tau)
      B     <- self$data$XtXm1 %*% t(self$data$X) %*% (self$data$Y - MtauT)
      ddiag  <- colMeans(R^2 - 2 * R * MtauT + (M^2 + outer(rep(1, self$n), S)) %*% t(tau))
      dm1  <- switch(private$res_covariance,
                     "diagonal"  = 1 / as.vector(ddiag),
                     "spherical" = rep(1/mean(ddiag), self$p))
      alpha <- colMeans(tau)
      OmegaQ <- private$get_OmegaQ(crossprod(M)/self$n +  diag(S, self$Q, self$Q))

      list(B = B, OmegaQ = OmegaQ, dm1 = dm1, alpha = alpha, tau = tau, M = M, S = S)
    },

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      private$C <- private$tau <- private$clustering_approx(reg_res$R)
      SigmaQ    <- private$heuristic_SigmaQ_from_Sigma(reg_res$Sigma)
      OmegaQ    <- private$get_OmegaQ(SigmaQ)
      list(B = reg_res$B, OmegaQ = OmegaQ)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), OmegaQ (groups precision matrix))
    model_par  = function() {
      parameters       <- super$model_par
      parameters$alpha <- private$alpha
      parameters},
    #' @field nb_param number of parameters in the model
    nb_param = function() {as.integer(super$nb_param + self$Q - 1)}, # adding alpha
    #' @field var_par a list with the matrices of the variational parameters: M (means), S (variances), tau (posterior group probabilities)
    var_par    = function() list(M = private$M,  S = private$S, tau = private$tau),
    #' @field clustering a list of labels giving the clustering obtained in the model
    clustering = function() get_clusters(private$tau),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      if (!private$approx){
        res <- .5 * self$n * self$Q * log(2 * pi * exp(1)) +
          .5 * self$n * sum(log(private$S)) - sum(xlogx(private$tau))
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model
    fitted = function(){
      if(private$approx) {
        res <- self$data$X %*% private$B
      } else {
        res <- self$data$X %*% private$B + private$M %*% t(private$tau)
      }
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value) {
      paste(private$res_covariance, "normal-block model with", self$Q, "unknown blocks")
    }
  )
)
