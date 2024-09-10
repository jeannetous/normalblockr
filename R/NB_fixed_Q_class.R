## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_Q #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with fixed Q (number of groups)
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param Q number of blocks
#' @param sparsity to add on blocks precision matrix
#' @param clustering_init to propose an initial clustering
#' @export
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
    #' @return A new [`NB_fixed_Q`] object
    initialize = function(Y, X, Q, sparsity = 0, clustering_init = NULL) {
      super$initialize(Y, X, Q, sparsity)
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

    compute_loglik_old  = function(B, dm1, omegaQ, alpha, tau, M, S) {
      ## problem dimensions
      n   <- self$n; p <-  self$p; d <-  self$d; Q <-  self$Q

      R              <- t(self$Y - self$X %*% B)
      ones           <- as.vector(rep(1, n))
      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      # expectation of log(p(Y | W, C))
      elbo <- - 0.5 * n * p * log(2 * pi) + 0.5 * n * sum(log(dm1))
      elbo <- elbo - 0.5 * sum(dm1 * ((R^2 + (tau %*% t(M^2)) - 2 * R * (tau %*% t(M))) %*% ones + n * (tau %*% S)))

      # expectation of log(p(W))
      elbo <- elbo - 0.5 * n * Q * log(2 * pi) + 0.5 * n * log_det_omegaQ
      elbo <- elbo - 0.5 * sum(crossprod(ones, (M %*% omegaQ) * M))
      elbo <- elbo - 0.5 * n * S %*% diag(omegaQ)

      # expectation of log(p(C))
      elbo <- elbo + sum(crossprod(log(alpha), t(tau)))

      # Entropy term for W
      elbo <- elbo + 0.5 * n * Q * log(2 * pi * exp(1)) + .5 * n * sum(log(S))

      # Entropy term for C
      elbo <- elbo - sum(xlogx(tau))

      if (self$sparsity == 0 ) {
        elbo
      }else {
        elbo - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
    },

    compute_loglik  = function(B, dm1, omegaQ, alpha, tau, M, S) {

      log_det_omegaQ <- as.numeric(determinant(omegaQ, logarithm = TRUE)$modulus)

      J <- -.5 * self$n * self$p * log(2 * pi * exp(1)) + .5 * self$n * sum(log(dm1))
      J <- J  + .5 * self$n * log_det_omegaQ
      J <- J + sum(tau %*% log(alpha))
      J <- J - sum(xlogx(tau)) + .5 * self$n * sum(log(S))

      if (self$sparsity > 0) {
        ## when not sparse, this terms equal -n Q /2 by definition of OmegaQ_hat and simplifies
        J <- J + self$n *self$Q / 2 - .5 * sum(diag(omegaQ %*% (crossprod(M) + self$n * diag(S))))
        # J <- J - self$sparsity * sum(abs(self$sparsity_weights * omegaQ))
      }
      J
    },

    EM_initialize = function() {
      B       <- private$XtXm1 %*% t(self$X) %*% self$Y
      R       <- self$Y - self$X %*% B
      if(is.null(self$clustering_init)){
        tau     <- as_indicator(kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster)
      }else{
        if(is.vector(self$clustering_init)){tau <- as_indicator(self$clustering_init)
        }else{ tau <- self$clustering_init}
      }
      tau     <- check_one_boundary(check_zero_boundary(tau))
      alpha   <- colMeans(tau)
      S       <- rep(0.1, self$Q)
      M       <- matrix(rep(0, self$n * self$Q), nrow = self$n)
      dm1     <- as.vector(rep(1, self$p))
      omegaQ  <- diag(rep(1, self$Q))
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

      sigmaQ <- crossprod(M)/self$n +  diag(S, self$Q, self$Q)
      if (self$sparsity == 0) {
        omegaQ <- solve(sigmaQ)
      }else {
        glasso_out <- glassoFast::glassoFast(sigmaQ, rho = self$sparsity * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) stop("GLasso fails")
        omegaQ <- Matrix::symmpart(glasso_out$wi)
      }

      list(B = B, dm1 = dm1, omegaQ = omegaQ, alpha = alpha, tau = tau, M = M, S = S)
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
      ent <- 0.5 * self$n * self$Q * log(2 * pi* exp(1)) + .5 * self$n * sum(log(private$S))
      ent <- ent - sum(xlogx(private$tau))
      return(ent)
    },
    #' @field fitted Y values predicted by the model Y values predicted by the model
    fitted = function(){
      inferred_C <- t(apply(private$tau, 1, function(x) as.integer(x == max(x))))
      self$X %*% private$B + private$M  %*%  t(inferred_C)
    }
    ),
)
