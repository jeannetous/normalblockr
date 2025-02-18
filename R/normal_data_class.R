## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS normal_data ##################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
normal_data <- R6::R6Class(
  classname = "normal_data",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses
    Y  = NULL,
    #' @field X the matrix of covariates
    X = NULL,
    #' @field XtXm1 inverse of XtX, useful for inference
    XtXm1 = NULL,

    #' @description Create a new [`normal_data`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    initialize = function(Y, X) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y and Xmust be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y <- Y
      self$X <- X
      self$XtXm1   <- solve(crossprod(X, X))
    }
  )
)
