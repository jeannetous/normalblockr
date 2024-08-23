## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @export
NB <- R6::R6Class(
  classname = "NB",
  inherit = MVEM,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field sparsity penalty on the network density
    sparsity = NULL,
    #' @field sparsity_weights distribution of sparsity on the lement of the network
    sparsity_weights = NULL,
    #' @field Q number of blocks
    Q = NULL,

    #' @description Create a new [`NB`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param Q required number of groups
    #' @param sparsity penalty on the network density
    #' @return A new [`nb_fixed`] object
    initialize = function(Y, X, Q, sparsity = 0) {
      super$initialize(Y, X)
      self$Q <- Q
      private$omegaQ <- diag(1, Q, Q)
      self$sparsity <- sparsity
      if (sparsity > 0) {
        sparsity_weights <- matrix(1, self$Q, self$Q)
        diag(sparsity_weights) <- 0
        self$sparsity_weights  <- sparsity_weights
      }
    },

    #' @description
    #' Update a [`NB`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param omegaQ groups inverse variance matrix
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, ll_list = NA) {
      super$update(B=B, dm1=dm1, ll_list=ll_list)
      if (!anyNA(omegaQ)) private$omegaQ  <- omegaQ
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    omegaQ    = NA   # groups variance matrix
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function() as.integer(super$nb_param + .5 * self$Q * (self$Q + 1)),
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), omegaQ (groups precision matrix))
    model_par  = function() list(B = private$B, dm1 = private$dm1, omegaQ = private$omegaQ),
    #' @field penalty (penalty on log-likelihood due to sparsity)
    penalty = function() - self$sparsity * sum(abs(self$sparsity_weights * private$omegaQ)),
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      res <- super$criteria
      res$Q <- self$Q
      res
    }
  )
)
