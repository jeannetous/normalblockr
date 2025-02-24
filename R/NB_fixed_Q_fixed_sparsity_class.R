## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q_fixed_sparsity ####################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed Q (number of groups)
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param Q number of blocks
#' @param penalty to add on blocks precision matrix for sparsity
#' @param control structured list for specific parameters (including initial clustering proposal)
NB_fixed_Q_fixed_sparsity <- R6::R6Class(
  classname = "NB_fixed_Q_fixed_sparsity",
  inherit = NB_fixed_sparsity,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field clustering_init model initial clustering
    clustering_init = NULL,
    #' @field clustering_method to use for clustering with heuristic inference method
    clustering_method = NULL,
    #' @field fixed_tau whether tau should be fixed at clustering_init during optimization, useful for stability selection
    fixed_tau = NULL,

    #' @description Create a new [`NB_fixed_Q_fixed_sparsity`] object.
    #' @param Q required number of groups
    #' @param control structured list for specific parameters
    #' @return A new [`NB_fixed_Q_fixed_sparsity`] object
    initialize = function(data, Q, penalty = 0, control = normal_control()) {
      super$initialize(data, Q, penalty, control)
      if (Q > ncol(self$data$Y)) stop("There cannot be more blocks than there are entities to cluster.")
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
      if(control$inference_method == "heuristic"){
        self$clustering_method <- control$clustering_method
      }
    },

    #' @description
    #' Update a [`NB_fixed_Q_fixed_sparsity`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param OmegaQ groups inverse variance matrix
    #' @param alpha vector of groups probabilities
    #' @param tau posterior probabilities for group affectation
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB_fixed_Q_fixed_sparsity`] object
    update = function(B = NA, OmegaQ = NA, dm1 = NA, alpha = NA, tau = NA,
                      M = NA, S = NA, ll_list = NA) {
      super$update(B, OmegaQ, dm1, ll_list)
      if (!anyNA(alpha)) private$alpha <- alpha
      if (!anyNA(tau))   private$tau   <- tau
      if (!anyNA(M))     private$M     <- M
      if (!anyNA(S))     private$S     <- S
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    alpha   = NA, # vector of groups probabilities
    M       = NA, # variational mean for posterior distribution of W
    S       = NA, # variational diagonal of variances for posterior distribution of W
    tau     = NA, # posterior probabilities for group affectation

    compute_loglik  = function(B, dm1, OmegaQ, alpha, tau, M, S) {
      if(self$inference_method == "integrated"){
        log_det_OmegaQ <- as.numeric(determinant(OmegaQ, logarithm = TRUE)$modulus)

        J <- -.5 * self$n * self$p * log(2 * pi * exp(1)) + .5 * self$n * sum(log(dm1))
        J <- J  + .5 * self$n * log_det_OmegaQ
        J <- J + sum(tau %*% log(alpha))
        J <- J - sum(xlogx(tau)) + .5 * self$n * sum(log(S))

        if (self$penalty > 0) {
          ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat and simplifies
          J <- J + self$n *self$Q / 2 - .5 * sum(diag(OmegaQ %*% (crossprod(M) + self$n * diag(S, self$Q, self$Q))))
          J <- J - self$penalty * sum(abs(self$sparsity_weights * OmegaQ))
        }
      }else{
        J <- private$heuristic_loglik(B, OmegaQ)
      }
      J
    },

    get_heuristic_parameters = function(){
      reg_res <- private$multivariate_normal_inference()
      private$C <- private$get_clustering(Sigma = reg_res$Sigma, R = reg_res$R)
      SigmaQ  <- private$heuristic_SigmaQ_from_Sigma(reg_res$Sigma)
      OmegaQ  <- private$get_Omega(SigmaQ)
      list(B = reg_res$B, OmegaQ = OmegaQ,dm1 = NA, alpha = NA, tau = private$C,
           M = NA, S = NA)
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
    nb_param = function() {as.integer(super$nb_param + self$Q - 1)},
    #' @field var_par a list with the matrices of the variational parameters: M (means), S (variances), tau (posterior group probabilities)
    var_par    = function() list(M = private$M,  S = private$S, tau = private$tau),
    #' @field clustering a list of labels giving the clustering obtained in the model
    clustering = function() get_clusters(private$tau),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() {
      ent <- .5 * self$n * self$Q * log(2 * pi* exp(1)) + .5 * self$n * sum(log(private$S))
      ent <- ent - sum(xlogx(private$tau))
      ent
    },
    #' @field fitted Y values predicted by the model
    fitted = function(){
      if(self$inference_method == "integrated"){self$data$X %*% private$B + tcrossprod(private$M, private$tau)
        }else{self$data$X %*% private$B}}
  ),
)


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q_fixed_sparsity_diagonal ###########
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' R6 class for normal-block model with fixed groups and diagonal residual covariance
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param Q number of blocks
#' @param penalty to add on blocks precision matrix for sparsity
NB_fixed_Q_fixed_sparsity_diagonal <- R6::R6Class(
  classname = "NB_fixed_Q_fixed_sparsity_diagonal",
  inherit = NB_fixed_Q_fixed_sparsity,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    EM_initialize = function() {
      if(any(sapply(self$model_par, function(x) any(is.na(x)))) | any(sapply(self$var_par, function(x) any(is.na(x))))){
        B       <- self$data$XtXm1 %*% t(self$data$X) %*% self$data$Y
        R       <- self$data$Y - self$data$X %*% B
        if(is.null(self$clustering_init)){
          clustering <- kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster
          if(length(unique(clustering)) < self$Q){
            # We try to ensure the optimization does not start with an empty cluster
            clustering <- cutree( ClustOfVar::hclustvar(t(R)), Q)
          }
          tau     <- as_indicator(clustering)
          if(min(colSums(tau)) < 1) warning("Initialization failed to place elements in each cluster")
        }else{
          if(is.vector(self$clustering_init)){tau <- as_indicator(self$clustering_init)
          }else{ tau <- self$clustering_init}
        }

        tau     <- check_one_boundary(check_zero_boundary(tau))
        alpha   <- colMeans(tau)
        S       <- rep(0.1, self$Q)
        M       <- matrix(rep(0, self$n * self$Q), nrow = self$n)
        # dm1     <- as.vector(rep(1, self$p))
        dm1     <- as.vector(1/colMeans((self$data$Y - self$data$X %*% B)^2))
        OmegaQ  <- diag(rep(1, self$Q), self$Q, self$Q)
      }else{
        B <- private$B ; dm1 <- private$dm1 ; OmegaQ <- private$OmegaQ
        alpha <- private$alpha
        tau <- private$tau ; M <- private$M ; S <- private$S
        parameters <-list(B = B, OmegaQ = OmegaQ, dm1 = dm1, alpha = alpha, tau = tau, M = M, S = S)
      }
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
      dm1   <- 1/colMeans(R^2 - 2 * R * MtauT + (M^2 + outer(rep(1, self$n), S)) %*% t(tau))
      alpha <- colMeans(tau)
      OmegaQ <- private$get_Omega(crossprod(M)/self$n +  diag(S, self$Q, self$Q))

      list(B = B, OmegaQ = OmegaQ, dm1 = dm1, alpha = alpha, tau = tau, M = M, S = S)
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value){paste("diagonal normal-block model with", self$Q, "unknown blocks")}
  )
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q_fixed_sparsity_spherical ##########
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for normal-block model with fixed groups and spherical residual covariance
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param Q number of blocks
#' @param penalty to add on blocks precision matrix for sparsity
NB_fixed_Q_fixed_sparsity_spherical <- R6::R6Class(
  classname = "NB_fixed_Q_fixed_sparsity_spherical",
  inherit = NB_fixed_Q_fixed_sparsity,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    EM_initialize = function() {
      if(any(sapply(self$model_par, function(x) any(is.na(x)))) | any(sapply(self$var_par, function(x) any(is.na(x))))){
        B       <- self$data$XtXm1 %*% t(self$data$X) %*% self$data$Y
        R       <- self$data$Y - self$data$X %*% B
        if(is.null(self$clustering_init)){
          clustering <- kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster
          if(length(unique(clustering)) < self$Q){
            # We try to ensure the optimization does not start with an empty cluster
            clustering <- cutree( ClustOfVar::hclustvar(t(R)), Q)
          }
          tau     <- as_indicator(clustering)
          if(min(colSums(tau)) < 1) warning("Initialization failed to place elements in each cluster")
        }else{
          if(is.vector(self$clustering_init)){tau <- as_indicator(self$clustering_init)
          }else{ tau <- self$clustering_init}
        }
        self$clustering_init <- get_clusters(tau)
        tau     <- check_one_boundary(check_zero_boundary(tau))
        alpha   <- colMeans(tau)
        S       <- rep(0.1, self$Q)
        M       <- matrix(rep(0, self$n * self$Q), nrow = self$n)
        dm1     <- rep(1/mean((self$data$Y - self$data$X %*% B)^2), self$p)
        OmegaQ  <- diag(rep(1, self$Q), self$Q, self$Q)
      }else{
        B <- private$B ; dm1 <- private$dm1 ; OmegaQ <- private$OmegaQ
        alpha <- private$alpha
        tau <- private$tau ; M <- private$M ; S <- private$S
      }

      list(B = B, OmegaQ = OmegaQ, dm1 = dm1, alpha = alpha, tau = tau, M = M, S = S)
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
      dm1_diag   <- 1/colMeans(R^2 - 2 * R * MtauT + (M^2 + outer(rep(1, self$n), S)) %*% t(tau))
      dm1   <- rep(mean(dm1_diag), self$p)
      alpha <- colMeans(tau)
      OmegaQ <- private$get_OmegaQ(crossprod(M)/self$n +  diag(S, self$Q, self$Q))

      list(B = B, OmegaQ = OmegaQ, dm1 = dm1, alpha = alpha, tau = tau, M = M, S = S)
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value){paste("spherical normal-block model with", self$Q, "unknown blocks")}
  )
)

