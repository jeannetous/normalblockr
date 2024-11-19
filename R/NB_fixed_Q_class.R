## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q ###################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed Q (number of groups)
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param Q number of blocks
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

    #' @description Create a new [`NB_fixed_Q`] object.
    #' @param Q required number of groups
    #' @param control structured list for specific parameters
    #' @return A new [`NB_fixed_Q`] object
    initialize = function(Y, X, Q, penalty = 0, control = NB_param()) {
      super$initialize(Y, X, Q, penalty, control)
      clustering_init <- control$clustering_init
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
    #' @param omegaQ groups inverse variance matrix
    #' @param alpha vector of groups probabilities
    #' @param tau posterior probabilities for group affectation
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB_fixed_Q`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, alpha = NA, tau = NA,
                      M = NA, S = NA, ll_list = NA) {
      super$update(B, dm1, omegaQ, ll_list)
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

    compute_loglik  = function(B, dm1, omegaQ, alpha, tau, M, S) {

      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      J <- -.5 * self$n * self$p * log(2 * pi * exp(1)) + .5 * self$n * sum(log(dm1))
      J <- J  + .5 * self$n * log_det_omegaQ
      J <- J + sum(tau %*% log(alpha))
      J <- J - sum(xlogx(tau)) + .5 * self$n * sum(log(S))

      if (self$penalty > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat and simplifies
        J <- J + self$n *self$Q / 2 - .5 * sum(diag(omegaQ %*% (crossprod(M) + self$n * diag(S, self$Q, self$Q))))
        J <- J - self$penalty * sum(abs(self$sparsity_weights * omegaQ))
      }
      J
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), omegaQ (groups precision matrix))
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
    fitted = function() self$X %*% private$B + tcrossprod(private$M, private$tau)
    ),
)


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q_diagonal ##########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' R6 class for normal-block model with fixed groups and diagonal residual covariance
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param Q number of blocks
#' @param penalty to add on blocks precision matrix for sparsity
NB_fixed_Q_diagonal <- R6::R6Class(
  classname = "NB_fixed_Q_diagonal",
  inherit = NB_fixed_Q,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    EM_initialize = function() {
      if(any(sapply(self$model_par, function(x) any(is.na(x)))) | any(sapply(self$var_par, function(x) any(is.na(x))))){
        B       <- private$XtXm1 %*% t(self$X) %*% self$Y
        R       <- self$Y - self$X %*% B
        if(is.null(self$clustering_init)){
          tau     <- as_indicator(kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster)
        }else{
          if(is.vector(self$clustering_init)){tau <- as_indicator(self$clustering_init)
          }else{ tau <- self$clustering_init}
        }
        self$clustering_init <- get_clusters(tau)
        tau     <- check_one_boundary(check_zero_boundary(tau))
        alpha   <- colMeans(tau)
        S       <- rep(0.1, self$Q)
        M       <- matrix(rep(0, self$n * self$Q), nrow = self$n)
        # dm1     <- as.vector(rep(1, self$p))
        dm1     <- as.vector(1/colMeans((self$Y - self$X %*% B)^2))
        omegaQ  <- diag(rep(1, self$Q), self$Q, self$Q)
      }else{
        B <- private$B ; dm1 <- private$dm1 ; omegaQ <- private$omegaQ
        alpha <- private$alpha
        tau <- private$tau ; M <- private$M ; S <- private$S
      }

      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, tau = tau, M = M, S = S)
    },

    EM_step = function(B, dm1, omegaQ, alpha, tau, M, S) {

      ## Auxiliary variables
      R     <- self$Y - self$X %*% B
      Gamma <- solve(omegaQ + diag(colSums(dm1 * tau), self$Q, self$Q))

      # E step
      M <- R %*% (dm1 * tau) %*% Gamma
      S <- diag(Gamma)

      if (self$Q > 1) {
        eta <- dm1 * (t(R) %*% M) - .5 * outer(dm1,  colSums(M^2) + self$n * S)
        eta <- eta + outer(rep(1, self$p), log(alpha)) - 1
        tau <- t(check_zero_boundary(check_one_boundary(apply(eta, 1, softmax))))
      }

      # M step
      MtauT <- M %*% t(tau)
      B     <- private$XtXm1 %*% t(self$X) %*% (self$Y - MtauT)
      dm1   <- 1/colMeans(R^2 - 2 * R * MtauT + (M^2 + outer(rep(1, self$n), S)) %*% t(tau))
      alpha <- colMeans(tau)
      omegaQ <- private$get_omegaQ(crossprod(M)/self$n +  diag(S, self$Q, self$Q))

      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, tau = tau, M = M, S = S)
    }
  )
)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q_spherical #########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for normal-block model with fixed groups and spherical residual covariance
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param Q number of blocks
#' @param penalty to add on blocks precision matrix for sparsity
NB_fixed_Q_spherical <- R6::R6Class(
  classname = "NB_fixed_Q_spherical",
  inherit = NB_fixed_Q,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    EM_initialize = function() {
      if(any(sapply(self$model_par, function(x) any(is.na(x)))) | any(sapply(self$var_par, function(x) any(is.na(x))))){
        B       <- private$XtXm1 %*% t(self$X) %*% self$Y
        R       <- self$Y - self$X %*% B
        if(is.null(self$clustering_init)){
          tau     <- as_indicator(kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster)
        }else{
          if(is.vector(self$clustering_init)){tau <- as_indicator(self$clustering_init)
          }else{ tau <- self$clustering_init}
        }
        self$clustering_init <- get_clusters(tau)
        tau     <- check_one_boundary(check_zero_boundary(tau))
        alpha   <- colMeans(tau)
        S       <- rep(0.1, self$Q)
        M       <- matrix(rep(0, self$n * self$Q), nrow = self$n)
        dm1     <- rep(1/mean((self$Y - self$X %*% B)^2), self$p)
        omegaQ  <- diag(rep(1, self$Q), self$Q, self$Q)
      }else{
        B <- private$B ; dm1 <- private$dm1 ; omegaQ <- private$omegaQ
        alpha <- private$alpha
        tau <- private$tau ; M <- private$M ; S <- private$S
      }

      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, tau = tau, M = M, S = S)
    },


    EM_step = function(B, dm1, omegaQ, alpha, tau, M, S) {

      ## Auxiliary variables
      R     <- self$Y - self$X %*% B
      Gamma <- solve(omegaQ + diag(colSums(dm1 * tau), self$Q, self$Q))

      # E step
      M <- R %*% (dm1 * tau) %*% Gamma
      S <- diag(Gamma)

      if (self$Q > 1) {
        eta <- dm1 * (t(R) %*% M) - .5 * outer(dm1,  colSums(M^2) + self$n * S)
        eta <- eta + outer(rep(1, self$p), log(alpha)) - 1
        tau <- t(check_zero_boundary(check_one_boundary(apply(eta, 1, softmax))))
      }

      # M step
      MtauT <- M %*% t(tau)
      B     <- private$XtXm1 %*% t(self$X) %*% (self$Y - MtauT)
      dm1_diag   <- 1/colMeans(R^2 - 2 * R * MtauT + (M^2 + outer(rep(1, self$n), S)) %*% t(tau))
      dm1   <- rep(mean(dm1_diag), self$p)
      alpha <- colMeans(tau)
      omegaQ <- private$get_omegaQ(crossprod(M)/self$n +  diag(S, self$Q, self$Q))

      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, tau = tau, M = M, S = S)
    }
  )
)

