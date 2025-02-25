## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS normal_fixed_sparsity ########################
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
    #' @field inference_method which method should be used to infer parameters
    inference_method = NULL,

    #' @description Create a new [`normal_models`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param control structured list of more specific parameters, to generate with NB_control (useful only for NB objects)
    #' @return A new [`nb_fixed`] object
    initialize = function(data, control = NB_control()) {
      self$data <- data
      self$inference_method <- control$inference_method
      private$ll_list <- 0
    },

    #' @description
    #' Update a [`normal_models`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of inverse variance matrix (variables level)
    #' @param ll_list  list of log-lik (elbo) values
    #' @return Update the current [`normal`] object
    update = function(B = NA, dm1 = NA, ll_list = NA) {
      if (!anyNA(B))          private$B       <- B
      if (!anyNA(dm1))        private$dm1     <- dm1
      if (!anyNA(ll_list))    private$ll_list <- ll_list
    },

    #' @description calls EM optimization and updates relevant fields
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    #' @return optimizes the model and updates its parameters
    optimize = function(niter = 100, threshold = 1e-4) {
      if(self$inference_method == "integrated"){
        optim_out <- private$EM_optimize(niter, threshold)
      }else{
        optim_out <- private$heuristic_optimize()}
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
    B          = NA, # regression matrix
    dm1       = NA,  # diagonal vector of inverse variance matrix (variables level)
    ll_list    = NA, # list of log-likelihoods or ELBOs

    compute_loglik  = function() {},

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated EM inference------------------
    EM_optimize = function(niter, threshold) {
      parameters <- do.call(private$EM_initialize, list())
      ll_list    <- do.call(private$compute_loglik, parameters)
      for (h in 2:niter) {
        parameters <- do.call(private$EM_step, parameters)
        ll_list    <- c(ll_list, do.call(private$compute_loglik, parameters))
        if (abs(ll_list[h] - ll_list[h - 1]) < threshold)
          break
      }
      c(parameters, list(ll_list = ll_list))
    },
    EM_step = function() {},
    EM_initialize = function() {},

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference----------------------
    heuristic_optimize = function(){
      parameters <- do.call(private$get_heuristic_parameters, list())
      ll_list    <- do.call(private$compute_loglik, parameters)
      c(parameters, list(ll_list = ll_list))
    },

    get_heuristic_parameters = function(){},

    multivariate_normal_inference = function(){
      B       <- self$data$XtXm1 %*% t(self$data$X) %*% self$data$Y
      R       <- self$data$Y - self$data$X %*% B
      Sigma   <- (t(R) %*% R) / self$n
      list(B = B, R = R, Sigma = Sigma)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples
    n = function() nrow(self$data$Y),
    #' @field p number of responses per sample
    p = function() ncol(self$data$Y),
    #' @field d number of variables (dimensions in X)
    d = function() ncol(self$data$X),
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), OmegaQ (groups precision matrix))
    model_par = function() list(B = private$B, dm1 = private$dm1, OmegaQ = private$OmegaQ),
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
