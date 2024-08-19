###############################################################################
###############################################################################
n     <- 200
p     <- 50
d     <- 2
minB  <- 2
maxB  <- 4
minX  <- rep(0,  d)
maxX  <- rep(10, d)
minD  <- 0.2
maxD  <- 3
minKappa <- 0.1
maxKappa <- 0.6
B     <- matrix(rep(1, d * p), nrow = d)
for (dim in 1:d) B[dim, ] <- runif(p, min = minB, max = maxB)
X     <- matrix(rep(1, d * n), nrow = d)
for (dim in 1:d) X[dim, ] <- runif(n, min = minX[[dim]], max = maxX[[dim]])
D     <- diag(runif(p, min <- minD, max <- maxD))
Y     <-  t(B) %*% X + t(MASS::mvrnorm(n, mu = matrix(rep(0, p), p, 1), Sigma = D))
kappa <- runif(p, min = minKappa, max = maxKappa)
Y[matrix(rbinom(n * p, size <- 1, prob <- rep(kappa, n)), nrow = p, ncol = n) == 1] <- 0
Y    <- t(Y)
X    <- t(X)
###############################################################################
###############################################################################

test_that("normal_zi: check dimensions, optimization and field access", {
  model <- normal_zi$new(Y, X)
  model$optimize(niter = 60)
  params <- model$get_model_parameters()
  expect_equal(params$n, nrow(Y))
  expect_equal(params$p, ncol(Y))
  expect_equal(params$d, ncol(X))
})
