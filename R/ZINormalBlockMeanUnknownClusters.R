## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZINormalBlockMeanUnknownClusters #############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Zero-Inflated Mean-Block Model with Unknown Clustering
#'
#' R6 class for a zero-inflated Normal-Block-Mean model with a fixed number of
#' clusters (but unknown clustering), inferred by variational EM. Sigma is
#' diagonal or spherical here. See [NormalBlockMeanBase] for why a full one
#' is out of reach under a mask.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' Y <- ex$Y; Y[runif(length(Y)) < 0.2] <- 0
#' data <- NormalBlockData$new(Y, ex$X)
#' model <- normal_block(data, blocks = 3, model = "mean", zero_inflation = TRUE)
#' model$clustering
#' @export
ZINormalBlockMeanUnknownClusters <- R6::R6Class(
  classname = "ZINormalBlockMeanUnknownClusters",
  inherit   = NormalBlockMeanBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`ZINormalBlockMeanUnknownClusters`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param q number of clusters
    #' @param sparsity unused here, kept for signature symmetry (must be 0)
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`ZINormalBlockMeanUnknownClusters`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control, zero_inflation = TRUE)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    optim_initialize = function() {
      if (private$warm_started) list(B = private$B, Omega = private$Omega, tau = private$C)
      else private$get_heuristic_parameters()
    },

    get_heuristic_parameters = function() {
      fit <- private$zi_mean_inference()
      ## cluster the fitted mean trajectory rather than B itself, as in the
      ## non-ZI class -- see NormalBlockMeanUnknownClusters
      if (anyNA(private$C))
        private$C <- private$heuristic_clustering(self$data$X %*% fit$B)
      tau <- check_one_boundary(check_zero_boundary(private$C))
      tau <- tau / rowSums(tau)
      list(B     = private$heuristic_cluster_B_from_variable_B(fit$B, tau),
           Omega = fit$Omega, tau = tau)
    },

    heuristic_optimize = function(control) {
      private$niter <- 0
      init <- private$get_heuristic_parameters()
      list(B = init$B, Omega = init$Omega, C = init$tau, alpha = colMeans(init$tau))
    },

    ## Runs the VEM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## ZINormalBlockMeanUnknownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- ZINormalBlockMeanUnknownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        zeros_bar = self$data$zeros_bar, zi_cond_mean = private$ZI_cond_mean,
        B0 = init$B, Omega0 = init$Omega, tau0 = init$tau,
        noise_covariance = private$res_covariance,
        fixed_tau = isTRUE(control$fixed_tau),
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, Omega = res$Omega, C = res$C, alpha = res$alpha,
           ll_list = res$objective)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field fitted Y values predicted by the model, in Y's original units
    fitted = function() {
      private$rescale_to_original(
        (self$data$X %*% private$B %*% t(private$C)) * self$data$zeros_bar)
    },
    #' @field var_par a list with the variational parameter: tau (posterior group probabilities)
    var_par = function() list(tau = private$C),
    #' @field model_par a list with model parameters: B, Omega and kappa (zero-inflation probabilities)
    model_par = function() c(super$model_par, list(kappa = private$kappa)),
    #' @field nb_param number of parameters in the model
    nb_param = function() {
      as.integer(super$nb_param + self$q - 1 + self$p * self$data$d0) # adding alpha and kappa
    },
    #' @field entropy Entropy of the conditional distribution
    entropy = function() -sum(xlogx(private$C)),
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function()
    {paste("zero-inflated", private$res_covariance, "normal-block-mean model with unknown blocks")}
  )
)
