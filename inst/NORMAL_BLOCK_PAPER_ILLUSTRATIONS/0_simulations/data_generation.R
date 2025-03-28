library(MASS)

#' Generate data under the model
#' @param parameters: required named list of parameters: B (regression matrix),
#' C (clustering), Sigma (cluster covariance matrix), D (diagonal variance matrix
#' for the variables), kappa (zero-inflation probability for each variable)
#' @param X: covariates matrix
normal_block_data <- function(param, X){
  n <- nrow(X) ; p <- nrow(param$C) ; Q <- ncol(param$C)
  W <- mvrnorm(n, mu=matrix(rep(0, Q), Q, 1), Sigma=param$Sigma)
  epsilon <- mvrnorm(n, mu=matrix(rep(0, p), p, 1), Sigma=param$D)
  Y = X  %*% param$B + tcrossprod(W, param$C) + epsilon
  Y[matrix(data = rbinom(n * p, size = 1, prob = rep(param$kappa, each = n)),
           nrow = n,ncol = p) == 1] <- 0
  Y
}

