## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_data ##################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param X0 zero-inflation design matrix, if applicable.
#' @export
NB_data <- R6::R6Class(
  classname = "NB_data",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses
    Y  = NULL,
    #' @field X the matrix of covariates
    X = NULL,
    #' @field X0 the matrix of zero-inflation covariates, if applicable
    X0 = NULL,
    #' @field n sample size
    n  = NULL,
    #' @field d number of covariates
    d = NULL,
    #' @field d0 number of zero-inflation covariates, if applicable
    d0 = NULL,
    #' @field p number of variables
    p = NULL,
    #' @field XtXm1 inverse of XtX, useful for inference
    XtXm1 = NULL,
    #' @field XtY useful for inference
    XtY = NULL,
    #' @field npY total number of non zeros in Y
    npY = NULL,
    #' @field nY total number of non zeros for each column/variable in Y
    nY = NULL,
    #' @field zeros where are the zero in Y
    zeros = NULL,
    #' @field zeros_bar where are the non-zeros in Y
    zeros_bar = NULL,

    #' @description Create a new [`NB_data`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param X0 zero-inflation design matrix, if applicable.
    initialize = function(Y, X, X0 = NULL) {
      stopifnot("Y and X must be matrices" = all(is.matrix(Y), is.matrix(X))) |> try()
      stopifnot("Y and X must have the same number of rows" = (nrow(Y) == nrow(X))) |> try()
      self$Y <- Y
      self$X <- X
      self$n <- nrow(Y)
      self$p <- ncol(Y)
      self$d <- ncol(X)
      if(is.null(X0)){X0 <- matrix(rep(1, self$n))}
      self$X0 <- X0
      self$d0 <- ncol(X0)
      self$XtXm1 <- solve(crossprod(X))
      self$XtY   <- crossprod(X, Y)
      self$zeros     <- 1 * (Y == 0)
      self$zeros_bar <- 1 * (Y != 0)
      self$npY <- sum(self$zeros_bar)
      self$nY  <- colSums(self$zeros_bar)
    }
  )
)
