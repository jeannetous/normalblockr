## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS normal_data ##################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @export
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
    #' @field n sample size
    n  = NULL,
    #' @field d number of covariates
    d = NULL,
    #' @field p number of variables
    p = NULL,
    #' @field XtXm1 inverse of XtX, useful for inference
    XtXm1 = NULL,
    #' @field XtY useful for inference
    XtY = NULL,

    #' @description Create a new [`normal_data`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    initialize = function(Y, X) {
      stopifnot("Y and X must be matrices" = all(is.matrix(Y), is.matrix(X))) |> try()
      stopifnot("Y and X must have the same number of rows" = (nrow(Y) == nrow(X))) |> try()
      self$Y <- Y
      self$X <- X
      self$n <- nrow(Y)
      self$p <- ncol(Y)
      self$d <- ncol(X)
      self$XtXm1 <- solve(crossprod(X))
      self$XtY   <- crossprod(X, Y)
    }
  )
)
