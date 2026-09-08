## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockData ##################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Data Container for Normal-Block Models
#'
#' R6 class holding the responses and design matrix used to fit a
#' normal-block model.
#' @param Y the matrix of responses.
#' @param X design matrix.
#' @param X0 zero-inflation design matrix, if applicable.
#' @param formula describes the relationship between Y and X, and X0 if applicable, useful if not all of X's or X0's covariates should be used, should be formatted ~ X1 + X2... | Z1 + Z2... with the Normal formula before the | and the ZI formula after the |
#' @param zeros an optional n x p 0/1 matrix marking the structural zeros of
#' Y. By default they are read off Y itself (`Y == 0`), which is what a
#' zero-inflated model expects. Pass it explicitly when the matrix handed
#' to the model is no longer the one carrying the zeros (typically the
#' residuals of a first stage, see [normal_block_sequential()]).
#' @param scale whether to rescale each column of Y by its own standard
#' deviation before fitting (default TRUE). Columns are *not* centered: the
#' model's own intercept (the constant or group-indicator columns a user is
#' expected to include in X) already absorbs each variable's mean, so
#' centering here would be redundant -- unless X has no such column, in
#' which case the residual model is misspecified regardless of this setting.
#' Rescaling matters because the normal-block model assumes a single shared
#' covariance value within each block (`Cov(Y_j, Y_j') = Var(W_k)` for j, j'
#' in block k); on the raw scale, that assumption is swamped whenever
#' variables in the same true block have very different baseline variances
#' (e.g. species with very different total abundance). The resulting `B`,
#' `Omega` and clustering are therefore properties of the *scaled* data, not
#' directly convertible back to the original units (the per-column scaling
#' factors are not the same within a block, so there is no single way to
#' "unscale" a shared block covariance back to a p x p matrix).
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' c(n = data$n, p = data$p, d = data$d)
#' @export
NormalBlockData <- R6::R6Class(
  classname = "NormalBlockData",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses (rescaled column-wise if `scale = TRUE`)
    Y  = NULL,
    #' @field Y_scale the per-column standard deviation Y was divided by (all
    #' 1's if `scale = FALSE`)
    Y_scale = NULL,
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
    #' @field XtX useful for inference in some cases
    XtX = NULL,
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

    #' @description Create a new [`NormalBlockData`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param X0 zero-inflation design matrix, if applicable.
    #' @param formula describes the relationship between Y and X, useful if not all of X's covariates should be used.
    #' @param scale whether to rescale each column of Y by its own standard
    #' deviation (no centering). Default TRUE -- see the class-level
    #' documentation for the rationale and its limits.
    #' @param zeros an optional 0/1 matrix of structural zeros, overriding the
    #' default `Y == 0`.
    initialize = function(Y, X, X0 = NULL, formula = NULL, scale = TRUE,
                          zeros = NULL) {
      stopifnot("Y and X must be matrices" = all(is.matrix(Y), is.matrix(X)))
      stopifnot("Y and X must have the same number of rows" = (nrow(Y) == nrow(X)))
      self$Y_scale <- if (scale) pmax(apply(Y, 2, sd), .Machine$double.eps) else rep(1, ncol(Y))
      self$Y <- Y / matrix(self$Y_scale, nrow(Y), ncol(Y), byrow = TRUE)
      self$n <- nrow(Y)
      self$p <- ncol(Y)
      self$formula <- formula
      fm_zi <- NULL
      if(is.null(formula)){self$X <- X
      }else{
        stopifnot("The formula should start with ~" = (formula[[1]] == "~"))
        if(formula[[2]][[1]] == "|"){
          fm <- as.formula(paste0("~", formula[[2]][2]))
          fm_zi <- as.formula(paste0("~", formula[[2]][3]))
        }else{
          fm <- formula}
        stopifnot("Covariates given in formula must be present in X" = (length(setdiff(all.vars(terms(fm)), colnames(X))) ==0))
        self$X <- model.matrix(fm, as.data.frame(X))
      }
      if(is.null(X0)){X0 <- matrix(rep(1, self$n))}
      if(is.null(fm_zi)){
        self$X0 <- X0
      }else{
        self$X0 <- model.matrix(fm_zi, as.data.frame(X0))
        stopifnot("Zero-inflation covariates given in the formula must be present in X0" = (length(setdiff(all.vars(terms(fm_zi)), colnames(X0))) ==0))
      }
      self$d0        <- ncol(self$X0)
      self$d         <- ncol(self$X)
      self$XtX       <- crossprod(self$X)
      self$XtXm1     <- chol2inv(chol(self$XtX))
      self$XtY       <- crossprod(self$X, self$Y)
      if (is.null(zeros)) {
        self$zeros <- 1 * (self$Y == 0)
      } else {
        stopifnot("zeros must be a matrix with the same dimensions as Y" =
          is.matrix(zeros) && all(dim(zeros) == dim(self$Y)))
        self$zeros <- 1 * (zeros != 0)
      }
      self$zeros_bar <- 1 - self$zeros
      self$npY       <- sum(self$zeros_bar)
      self$nY        <- colSums(self$zeros_bar)
    },

    #' @description Ordinary-least-squares fit of Y on X, with its residuals
    #' and their covariance. Computed once and memoized: it depends only on
    #' the data, yet every model in a collection over q used to recompute it
    #' (measured at 9% of a q = 1:30 variance-block collection on `brca_rppa`).
    #' @return a list with `B` (d x p), `R` (n x p residuals) and `Sigma`
    #' (p x p residual covariance)
    ols_fit = function() {
      if (is.null(private$ols_cache)) {
        B <- self$XtXm1 %*% self$XtY
        R <- self$Y - self$X %*% B
        private$ols_cache <- list(B = B, R = R, Sigma = cov(R))
      }
      private$ols_cache
    },

    #' @description Masked counterpart of `ols_fit()`: a per-variable weighted
    #' least-squares fit of B under the zero-inflation mask, with its inverse
    #' residual variances and residuals (see `zi_weighted_fit()`). Memoized for
    #' the same reason -- every zero-inflated model in a collection over q used
    #' to redo the same IRLS.
    #' @return a list with `B` (d x p), `dm1` (p) and `R` (n x p masked residuals)
    zi_ols_fit = function() {
      if (is.null(private$zi_ols_cache)) private$zi_ols_cache <- zi_weighted_fit(self)
      private$zi_ols_cache
    },

    #' @description Zero-inflation component: `p` independent logistic
    #' regressions of each variable's zero pattern on `X0`, and the fixed
    #' contribution they make to the log-likelihood. The (V)EM never revisits
    #' these, so they are a property of the data rather than of a model --
    #' hence computed once and memoized here. Every model in a collection over
    #' q used to refit all `p` regressions (measured at 53% of a q = 2:8
    #' zero-inflated mean-block collection).
    #' @return a list with `B0` (d0 x p), `kappa` (n x p zero-inflation
    #' probabilities) and `ZI_cond_mean` (a scalar)
    zi_fit = function() {
      if (is.null(private$zi_cache)) {
        ## d0 x p, one column of coefficients per variable. Built with an
        ## explicit matrix() rather than t(sapply(...)): sapply collapses to a
        ## plain vector when d0 == 1 and returns d0 x p when d0 > 1, so the
        ## transpose that made the first case work broke the second (a
        ## non-conformable X0 %*% B0 as soon as a second zero-inflation
        ## covariate was supplied).
        no_zeros <- self$npY == self$n * self$p
        B0 <- if (!no_zeros) {
          coefs <- lapply(1:self$p, function(j) {
            df <- data.frame("zeros" = self$zeros[, j], self$X0)
            glm(zeros ~ 0 + ., family = binomial(link = "logit"), data = df)$coefficients
          })
          matrix(unlist(coefs), nrow = self$d0, ncol = self$p)
        } else {
          matrix(rep(-Inf, self$p * self$d0), nrow = self$d0)
        }
        ## With no zeros to model, kappa is 0 by construction. Taking it from
        ## the -Inf sentinel instead only works when d0 == 1: a second
        ## covariate contributes x * -Inf, which is +Inf wherever x < 0, and
        ## the linear predictor collapses to NaN.
        kappa <- if (no_zeros) matrix(0, self$n, self$p)
                 else apply(self$X0 %*% B0, MARGIN = c(1, 2), FUN = sigmoid)
        private$zi_cache <- list(
          B0 = B0, kappa = kappa,
          ZI_cond_mean = sum(xlogy(self$zeros, kappa)) +
                         sum(xlogy(self$zeros_bar, 1 - kappa)))
      }
      private$zi_cache
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## memoized by ols_fit()/zi_fit(); both are pure functions of the fields
    ## set in initialize(), none of which change afterwards
    ols_cache     = NULL,
    zi_ols_cache  = NULL,
    zi_cache      = NULL
  )
)
