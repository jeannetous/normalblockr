## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS normal_models         ########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal_fixed_sparsity model
#' @param data contains the matrix of responses (Y) and the design matrix (X).
normal_models <- R6::R6Class(
  classname = "normal_models",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of normal_data class, with responses and design matrix
    data  = NULL,

    #' @description Create a new [`normal_models`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param control structured list of more specific parameters, to generate with NB_control (useful only for NB objects)
    #' @return A new [`nb_fixed`] object
    initialize = function(data, control = NB_control()) {
      self$data <- data
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Setters    ------------------------
    #' @description
    #' Update a [`normal_models`] object
    #'
    #' All possible parameters of the child classes
    #' @param B regression matrix
    #' @param dm1 diagonal vector of inverse variance matrix (variables level)
    #' @param OmegaQ groups inverse variance matrix
    #' @param gamma  variance of  posterior distribution of W
    #' @param mu  mean for posterior distribution of W
    #' @param kappa vector of zero-inflation probabilities
    #' @param rho posterior probabilities of zero-inflation
    #' @param alpha vector of groups probabilities
    #' @param tau posterior probabilities for group affectation
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param ll_list  list of log-lik (elbo) values
    #' @return Update the current [`normal`] object
    update = function(B = NA,
                      OmegaQ = NA,
                      dm1 = NA,
                      gamma = NA,
                      mu = NA,
                      kappa = NA,
                      rho = NA,
                      alpha = NA,
                      tau = NA,
                      M = NA,
                      S = NA,
                      ll_list = NA) {
      if (!anyNA(B))       private$B       <- B
      if (!anyNA(dm1))     private$dm1     <- dm1
      if (!anyNA(OmegaQ))  private$OmegaQ  <- OmegaQ
      if (!anyNA(gamma))   private$gamma   <- gamma
      if (!anyNA(kappa))   private$kappa   <- kappa
      if (!anyNA(rho))     private$rho     <- rho
      if (!anyNA(mu))      private$mu      <- mu
      if (!anyNA(alpha))   private$alpha   <- alpha
      if (!anyNA(tau))     private$tau     <- tau
      if (!anyNA(M))       private$M       <- M
      if (!anyNA(S))       private$S       <- S
      if (!anyNA(ll_list)) private$ll_list <- ll_list
    },

    #' @description calls optimization (EM or heuristic) and updates relevant fields
    #' @param control a list for controlling the optimization proces
    #' @return optimizes the model and updates its parameters
    optimize = function(control = list(niter = 100, threshold = 1e-4)) {
      optim_out <- private$optimizer(control)
      do.call(self$update, optim_out)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods------------------
    #' @param type char for line type (see plot.default)
    #' @param log char for logarithmic axes (see plot.default)
    #' @param neg boolean plot negative log-likelihood (useful when log="y")
    #' @description plots log-likelihood values during model optimization
    plot_loglik = function(type = "b", log = "", neg = FALSE) {
      neg <- ifelse(neg, -1, 1)
      plot(seq_along(self$objective), neg * self$objective, type = type, log = log)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS -------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    B         = NA, # regression matrix
    dm1       = NA, # diagonal vector of inverse variance matrix (variables level)
    C         = NA, # the matrix of species groups
    OmegaQ    = NA, # precision matrix for clusters
    kappa     = NA, # vector of zero-inflation probabilities
    alpha     = NA, # vector of groups probabilities
    rho       = NA, # posterior probabilities of zero-inflation
    tau       = NA, # posterior probabilities for group affectation
    gamma     = NA, # variance of  posterior distribution of W
    mu        = NA, # mean for posterior distribution of W
    M         = NA, # variational mean for posterior distribution of W
    S         = NA, # variational diagonal of variances for posterior distribution of W
    optimizer = NA, # a link to the function that perform the optimization
    ll_list   = NA, # list of log-likelihoods or ELBOs

    compute_loglik  = function() {},

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MLE of MV distribution
    multivariate_normal_inference = function(){
      B     <- self$data$XtXm1 %*% self$data$XtY
      R     <- self$data$Y - self$data$X %*% B
      Sigma <- cov(R)
      list(B = B, R = R, Sigma = Sigma)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated EM inference------------------
    EM_optimize = function(control) {
      parameters <- private$EM_initialize()
      ll_list    <- do.call(private$compute_loglik, parameters)
      for (h in 2:control$niter) {
        parameters <- do.call(private$EM_step, parameters)
        ll_list    <- c(ll_list, do.call(private$compute_loglik, parameters))
        if (abs(ll_list[h] - ll_list[h - 1]) < control$threshold)
          break
      }
      c(parameters, list(ll_list = ll_list))
    },
    EM_step   = function() {},
    EM_initialize = function() {}
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples
    n = function() self$data$n,
    #' @field p number of responses per sample
    p = function() self$data$p,
    #' @field d number of variables (dimensions in X)
    d = function() self$data$d,
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), OmegaQ (groups precision matrix))
    model_par = function() list(B = private$B, dm1 = private$dm1),
    #' @field loglik (or its variational lower bound)
    loglik = function() private$ll_list[[length(private$ll_list)]],
    #' @field nb_param number of parameters in the model
    nb_param = function() as.integer(self$p * self$d),
    #' @field entropy Entropy of the variational distribution when applicable
    entropy    = function() 0,
    #' @field deviance (or its variational lower bound)
    deviance = function() - 2 * self$loglik,
    #' @field BIC (or its variational lower bound)
    BIC = function() self$deviance + log(self$n) * self$nb_param,
    #' @field AIC (or its variational lower bound)
    AIC = function() self$deviance + 2 * self$nb_param,
    #' @field ICL variational lower bound of the ICL
    ICL        = function() self$BIC + 2 * self$entropy,
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      data.frame(nb_param = self$nb_param, loglik = self$loglik,
                 deviance = self$deviance, BIC = self$BIC, AIC = self$AIC, ICL = self$ICL)
    },
    #' @field objective evolution of the objective function during (V)EM algorithm
    objective = function() private$ll_list[-1]
  )
)
