## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarUnknownClusters ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Normal-Block Model with Unknown Clustering
#'
#' R6 class for a normal-block model with a fixed number of clusters (but
#' unknown clustering).
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' model$clustering
#' @export
NormalBlockVarUnknownClusters <- R6::R6Class(
  classname = "NormalBlockVarUnknownClusters",
  inherit = NormalBlockVarBase,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field fixed_tau whether tau should be fixed at clustering_init during optimization, useful for stability selection
    fixed_tau = NULL,

    #' @description Create a new [`NormalBlockVarUnknownClusters`] object.
    #' @param data contains the matrix of responses (Y) and the design matrix (X).
    #' @param q required number of groups
    #' @param sparsity sparsity penalty to add on blocks precision matrix for sparsity
    #' @param control structured list for specific parameters
    #' @return A new [`NormalBlockVarUnknownClusters`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control)
      self$fixed_tau <- control$fixed_tau
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference -----------------------

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      if (anyNA(private$C)) # if no initial clustering provided
        private$C <- private$heuristic_clustering(reg_res$R)
      private$C <- clip_probabilities(private$C)
      Sigmaq    <- private$heuristic_Sigmaq_from_Sigma(reg_res$Sigma)
      Omega    <- private$get_Omega(Sigmaq)
      dm1       <- private$dm1_from_residuals(reg_res$R)
      list(B = reg_res$B, Omega = Omega, dm1 = dm1,
           C = private$C, alpha = colMeans(private$C))
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated inference ------------------------

    optim_initialize = function() {
      if (private$warm_started) {
        list(B = private$B, Omega = private$Omega, dm1 = private$dm1,
             C = private$C, alpha = private$alpha, M = private$M, S = private$S)
      } else c(private$get_heuristic_parameters(),  list(
            M = matrix(rep(0, self$n * self$q), nrow = self$n),
            S = rep(0.1, self$q)
          )
      )
    },

    ## Runs the VEM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## NormalBlockVarUnknownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- NormalBlockVarUnknownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        B0 = init$B, dm1_0 = init$dm1, Omega0 = init$Omega,
        C0 = init$C, alpha0 = init$alpha, M0 = init$M, S0 = init$S,
        sparsity = self$sparsity, sparsity_weights = self$sparsity_weights,
        noise_covariance = private$res_covariance, fixed_tau = self$fixed_tau,
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, Omega = res$Omega, dm1 = res$dm1, alpha = res$alpha,
           C = res$C, M = res$M, S = res$S, ll_list = res$objective)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), Omega (groups precision matrix))
    model_par  = function() {
      parameters       <- super$model_par
      parameters$alpha <- private$alpha
      parameters},
    #' @field nb_param number of parameters in the model
    nb_param = function() {as.integer(super$nb_param + self$q - 1)}, # adding alpha
    #' @field var_par a list with the matrices of the variational parameters: M (means), S (variances), tau (posterior group probabilities)
    var_par    = function() list(M = private$M,  S = private$S, tau = private$C),
    #' @field entropy Entropy of the conditional distribution
    entropy    = function() {
      if (!private$approx){
        res <- .5 * self$n * self$q * log(2 * pi * exp(1)) + .5 * self$n * sum(log(private$S))
        res <- res - sum(xlogx(private$C))
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model, in Y's original units
    fitted = function(){
      res <- if (private$approx) {
        self$data$X %*% private$B
      } else {
        self$data$X %*% private$B + tcrossprod(private$M, private$C)
      }
      private$rescale_to_original(res)
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function() {
      paste(private$res_covariance, "normal-block-var model with", self$q, "unknown blocks")
    }
  )
)
