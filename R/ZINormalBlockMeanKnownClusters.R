## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZINormalBlockMeanKnownClusters ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Zero-Inflated Mean-Block Model with Known Clustering
#'
#' R6 class for a zero-inflated Normal-Block-Mean model with a known
#' clustering. Sigma is diagonal or spherical here. See
#' [NormalBlockMeanBase] for why a full one is out of reach under a mask.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' Y <- ex$Y; Y[runif(length(Y)) < 0.2] <- 0
#' data <- NormalBlockData$new(Y, ex$X)
#' model <- normal_block(data, blocks = ex$parameters$C, model = "mean",
#'                       zero_inflation = TRUE)
#' model$clustering
#' @export
ZINormalBlockMeanKnownClusters <- R6::R6Class(
  classname = "ZINormalBlockMeanKnownClusters",
  inherit   = NormalBlockMeanBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`ZINormalBlockMeanKnownClusters`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param C clustering matrix C_jk = 1 if species j belongs to cluster k
    #' @param sparsity unused here, kept for signature symmetry (must be 0)
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`ZINormalBlockMeanKnownClusters`] object
    initialize = function(data, C, sparsity = 0, control = NB_control()) {
      stopifnot("C must be a matrix" = is.matrix(C))
      stopifnot("There cannot be empty clusters" = min(colSums(C)) > 0)
      super$initialize(data, ncol(C), sparsity, control, zero_inflation = TRUE)
      private$C <- C
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    optim_initialize = function() {
      if (private$warm_started) list(B = private$B, Omega = private$Omega)
      else private$get_heuristic_parameters()
    },

    get_heuristic_parameters = function() {
      fit <- private$zi_mean_inference()
      list(B     = private$heuristic_cluster_B_from_variable_B(fit$B, private$C),
           Omega = fit$Omega)
    },

    heuristic_optimize = function(control) {
      private$niter <- 0
      private$get_heuristic_parameters()
    },

    ## Runs the EM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## ZINormalBlockMeanKnownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- ZINormalBlockMeanKnownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        zeros_bar = self$data$zeros_bar, zi_cond_mean = private$ZI_cond_mean,
        C = private$C, B0 = init$B, Omega0 = init$Omega,
        noise_covariance = private$res_covariance,
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, Omega = res$Omega, ll_list = res$objective)
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
    #' @field model_par a list with model parameters: B, Omega and kappa (zero-inflation probabilities)
    model_par = function() c(super$model_par, list(kappa = private$kappa)),
    #' @field nb_param number of parameters in the model
    nb_param = function() as.integer(super$nb_param + self$p * self$data$d0), # adding kappa
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function()
    {paste("zero-inflated", private$res_covariance, "normal-block-mean model with fixed blocks")}
  )
)
