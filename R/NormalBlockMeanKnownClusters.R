## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanKnownClusters #################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Mean-Block Model with Known Clustering
#'
#' R6 class for a Normal-Block-Mean model with a known clustering.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' model <- normal_block(data, blocks = ex$parameters$C, model = "mean")
#' model$plot()
#' @export
NormalBlockMeanKnownClusters <- R6::R6Class(
  classname = "NormalBlockMeanKnownClusters",
  inherit   = NormalBlockMeanBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockMeanKnownClusters`] object.
    #' @param data object of NormalMeanBlockData class, with responses and design matrix
    #' @param C clustering matrix C_jk = 1 if species j belongs to cluster k
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockMeanKnownClusters`] object
    initialize = function(data, C, sparsity = 0, control = NB_control()) {
      stopifnot("C must be a matrix" = is.matrix(C))
      stopifnot("There cannot be empty clusters" = min(colSums(C)) > 0)
      super$initialize(data, ncol(C), sparsity, control)
      private$C <- C
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    optim_initialize = function() {
      if (private$warm_started) {
        list(B = private$B, Omega = private$Omega)
      } else {
        private$get_heuristic_parameters()
      }
    },

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      Omega     <- private$omega_from_sigma(reg_res$Sigma)
      B         <- private$heuristic_cluster_B_from_variable_B(reg_res$B,
                                                               private$C)
      list(B = B, Omega = Omega)
    },


    B_estimator = function(Omega = private$Omega){
      return(self$data$XtXm1 %*% self$data$XtY %*% Omega %*% private$C %*% solve(t(private$C) %*% Omega %*% private$C))
    },

    Sigma_estimator = function(B = private$B){
      R     <- self$data$Y - self$data$X %*% B %*% t(private$C)
      Sigma <- crossprod(R) / self$n
      return(Sigma)
    },

    compute_loglik = function(B = private$B, Omega = private$Omega){
      log_det_Omega <- as.numeric(determinant(Omega, logarithm = TRUE)$modulus)
      R <- self$data$Y - self$data$X %*% B %*% t(private$C)
      l <- - self$n * self$p * log(2 * pi) / 2 + self$n * log_det_Omega / 2 - sum((R %*% Omega) * R) / 2
      ## penalized objective (see the C++ core): $loglik adds the penalty back
      if (self$sparsity > 0) l <- l - self$sparsity * sum(abs(self$sparsity_weights * Omega))
      return(as.numeric(l))
    },

    ## Moment-based fit (NB_control(heuristic = TRUE)): the heuristic
    ## initialization is the moment estimate, so there is nothing to iterate.
    ## No likelihood is computed, hence no ll_list (see NormalBlockBase's
    ## `loglik`, NA in this mode).
    heuristic_optimize = function(control) {
      private$niter <- 0
      private$get_heuristic_parameters()
    },

    ## Runs the EM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## NormalBlockMeanKnownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- NormalBlockMeanKnownClusters_fit(
        Y = self$data$Y, X = self$data$X, C = private$C,
        B0 = init$B, Omega0 = init$Omega,
        sparsity = self$sparsity, sparsity_weights = self$sparsity_weights,
        noise_covariance = private$res_covariance,
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, Omega = res$Omega, ll_list = res$objective)
    },

    ## Reference R implementation of the same recursion, kept since the C++
    ## port is being validated against it (test-cpp-normal-block-mean.R).
    EM_optimize_R = function(control){
      init_params <- private$optim_initialize()
      B           <- init_params$B
      Omega       <- init_params$Omega
      ll_prev     <- private$compute_loglik(B, Omega)
      ll_list     <- c(ll_prev)
      warm        <- NULL
      for(i in 1:control$niter){
        B     <- private$B_estimator(Omega)
        Sigma <- private$Sigma_estimator(B)
        og    <- private$omega_from_sigma_warm(Sigma, warm)
        Omega <- og$Omega; warm <- og$warm

        ll_current <- private$compute_loglik(B, Omega)
        ll_list     <- c(ll_list, ll_current)
        if(abs(ll_current - ll_prev) < control$threshold){
          break
        }else{ll_prev <- ll_current}
      }
      private$niter <- i
      list(B = B, Omega = Omega, ll_list = ll_list)
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
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function()
    {paste(private$res_covariance, "normal-block-mean model with fixed blocks")}
  )
)
