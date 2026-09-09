## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZINormalBlockVarUnknownClusters ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Zero-Inflated Normal-Block Model with Unknown Clustering
#'
#' R6 class for a zero-inflated normal-block model with a fixed number of
#' clusters (but unknown clustering).
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3, kappa = rep(0.3, 20))
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' model <- normal_block(data, blocks = 3, zero_inflation = TRUE,
#'                       control = NB_control(verbose = FALSE))
#' model$clustering
#' @export
ZINormalBlockVarUnknownClusters <- R6::R6Class(
  classname = "ZINormalBlockVarUnknownClusters",
  inherit   = NormalBlockVarBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field fixed_tau whether tau should be fixed at clustering_init during optimization, useful for stability selection
    fixed_tau = NULL,

    #' @description Create a new [`ZINormalBlockVarUnknownClusters`] object.
    #' @param data object of NormalBlockVarData class, with responses and design matrix
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param q required number of groups
    #' @param control structured list of more specific parameters
    #' @return A new [`ZINormalBlockVarUnknownClusters`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control, zero_inflation = TRUE)
      self$fixed_tau  <- control$fixed_tau
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS -------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference -----------------------

    get_heuristic_parameters = function() {
      zi_diag <- private$zi_diag_normal_inference()
      if (anyNA(private$C))
        private$C <- private$heuristic_clustering(zi_diag$R)
      private$C <- clip_probabilities(private$C)
      Sigmaq <- private$heuristic_Sigmaq_from_Sigma(cov(zi_diag$R))
      Omega <- private$get_Omega(Sigmaq)
      list(B = zi_diag$B, dm1 = zi_diag$dm1, Omega = Omega,
           alpha = colMeans(private$C), kappa = private$kappa,
           C = private$C)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated inference ------------------------

    optim_initialize = function() {
      if (private$warm_started) {
        list(B = private$B, dm1 = private$dm1, Omega = private$Omega, kappa = private$kappa,
             alpha = private$alpha, C = private$C, M = private$M, S = private$S)
      } else c(private$get_heuristic_parameters(),  list(
          M = matrix(rep(0, self$n * self$q), nrow = self$n),
          S = matrix(rep(0.1, self$n * self$q), nrow = self$n)
        )
      )
    },

    ## Runs the VEM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## ZINormalBlockVarUnknownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- ZINormalBlockVarUnknownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        zeros_bar = self$data$zeros_bar, zi_cond_mean = private$ZI_cond_mean,
        B0 = init$B, dm1_0 = init$dm1, Omega0 = init$Omega,
        C0 = init$C, alpha0 = init$alpha, M0 = init$M, S0 = init$S,
        sparsity = self$sparsity, sparsity_weights = self$sparsity_weights,
        noise_covariance = private$res_covariance, fixed_tau = self$fixed_tau,
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, dm1 = res$dm1, Omega = res$Omega, alpha = res$alpha,
           C = res$C, M = res$M, S = res$S, ll_list = res$objective)
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function() super$nb_param + self$p * self$d0 + self$q - 1, # adding kappa and alpha
    #' @field var_par a list with variational parameters
    var_par  = function() {list(M = private$M, S = private$S, tau = private$C)},
    #' @field model_par a list with model parameters: B (covariates), dm1 (species variance), Omega (blocks precision matrix), kappa (zero-inflation probabilities)
    model_par  = function() {
      par       <- super$model_par
      par$kappa <- private$kappa
      par$alpha <- private$alpha
      par
    },
    #' @field entropy Entropy of the conditional distribution
    entropy    = function() {
      if (!private$approx) {
        res <- 0.5 * self$n * self$q * log(2 * pi * exp(1)) + .5 * sum(log(private$S))
        res <- res - sum(xlogx(private$C))
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model, in Y's original units
    fitted = function(){
      if (private$approx) {
        res <- self$data$X %*% private$B
      } else {
        res <- self$data$X %*% private$B + tcrossprod(private$M, private$C)
      }
      res <- res * self$data$zeros_bar
      private$rescale_to_original(res)
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function()
    {paste("zero-inflated", private$res_covariance, "normal-block-var model with", self$q, "unknown blocks")}
  )
)
