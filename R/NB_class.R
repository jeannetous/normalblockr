## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS normal #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param blocks either explicit clustering or number of groups
NB <- R6::R6Class(
  classname = "NB",
  inherit   = normal,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses
    Y  = NULL,
    #' @field X the matrix of covariates
    X = NULL,
    #' @field penalty to apply on variance matrix when calling GLASSO
    penalty = NULL,

    #' @description Create a new [`normal`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param penalty penalty on the network density
    #' @return A new [`nb_fixed`] object
    initialize = function(Y, X,  penalty = 0) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y and Xmust be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y <- Y
      self$X <- X
      self$penalty <- penalty
      private$XtXm1   <- solve(crossprod(X, X))
    },

    #' @description
    #' Update a [`normal_baseline`] object
    #' @param B regression matrix
    #' @param Sigma  p-dimensional var-covar matrix
    #' @return Update the current [`normal_baseline`] object
    update = function(B = NA) {
      if (!anyNA(B))          private$B     <- B
    },

    #' @description calls appropriate optimization and updates relevant fields
    #' @return optimizes the model and updates its parameters
    optimize = function() {
      optim_out <- private$complete_optimization()
      do.call(self$update, optim_out)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    B         = NA, # regression matrix
    XtXm1     = NA # inverse of XtX, useful for inference
  ),

  active = list(
    d = function() ncol(self$X)
  )
)
