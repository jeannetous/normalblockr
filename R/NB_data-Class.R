## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_data ##################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param X0 zero-inflation design matrix, if applicable.
#' @param formula describes the relationship between Y and X, and X0 if applicable, useful if not all of X's or X0's covariates should be used, should be formatted ~ X1 + X2... | Z1 + Z2... with the Normal formula before the | and the ZI formula after the |
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
    #' @field formula describes the relationship between Y and X, and X0 if applicable, useful if not all of X's or X0's covariates should be used, should be formatted ~ X1 + X2... | Z1 + Z2... with the Normal formula before the | and the ZI formula after the |
    formula = NULL,
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
    #' @param formula describes the relationship between Y and X, useful if not all of X's covariates should be used.
    initialize = function(Y, X, X0 = NULL, formula = NULL) {
      stopifnot("Y and X must be matrices" = all(is.matrix(Y), is.matrix(X))) |> try()
      stopifnot("Y and X must have the same number of rows" = (nrow(Y) == nrow(X))) |> try()
      self$Y <- Y
      self$n <- nrow(Y)
      self$p <- ncol(Y)
      self$formula <- formula
      fm_zi <- NULL
      if(is.null(formula)){self$X <- X
      }else{
        stopifnot("The formula should start with ~" = (formula[[1]] == "~")) |> try()
        if(formula[[2]][[1]] == "|"){
          fm <- as.formula(paste0("~", formula[[2]][2]))
          fm_zi <- as.formula(paste0("~", formula[[2]][3]))
        }else{
          fm <- formula}
        stopifnot("Covariates given in formula must be present in X" = (length(setdiff(all.vars(terms(fm)), colnames(X))) ==0))|> try()
        self$X <- matrix(model.matrix(fm, as.data.frame(X)))
      }
      if(is.null(X0)){X0 <- matrix(rep(1, self$n))}
      if(is.null(fm_zi)){
        self$X0 <- X0
      }else{
        self$X0 <- matrix(model.matrix(fm_zi, as.data.frame(X0)))
        stopifnot("Zero-inflation covariates given in the formula must be present in X0" = (length(setdiff(all.vars(terms(fm_zi)), colnames(X0))) ==0))|> try()
      }
      self$d0 <- ncol(self$X0)
      self$d <- ncol(self$X)
      self$XtXm1 <- solve(crossprod(self$X))
      self$XtY   <- crossprod(self$X, Y)
      self$zeros     <- 1 * (Y == 0)
      self$zeros_bar <- 1 * (Y != 0)
      self$npY <- sum(self$zeros_bar)
      self$nY  <- colSums(self$zeros_bar)
    }
  )
)
