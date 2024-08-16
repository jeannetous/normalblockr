###############################################################################
###############################################################################
set.seed(3)
n     <- 200
p     <- 50
d     <- 2
Q     <- 4
Sigma <- matrix(0, Q, Q)
diag(Sigma) <- 1
Sigma[row(Sigma) != col(Sigma)] = 0.1
minB  <- -1
maxB  <- 1
minX  <- rep(-50, d)
maxX  <- rep(50, d)
minD  <- 0.2
maxD  <- 3
minKappa = 0.2
maxKappa = 0.6

B = matrix(rep(1, d*p), nrow=d)
for(dim in 1:d){B[dim,] = runif(p, min=0, max = 1)}
X = matrix(rep(1, d*n), nrow=d)
for(dim in 1:d){X[dim,] = runif(n, min=minX[[dim]], max = maxX[[dim]])}
C = matrix(rep(0, p*Q), nrow=p)
groups = sample(1 : Q, size = p, replace = TRUE)
for(dim in 1:p){C[dim, groups[[dim]]] = 1}
D = diag(runif(p, min = minD, max = maxD))
W = t(MASS::mvrnorm(n, mu=matrix(rep(0, Q), Q, 1), Sigma=Sigma))
epsilon = t(MASS::mvrnorm(n, mu=matrix(rep(0, p), p, 1), Sigma=D))
Y = t(B) %*% X + C %*% W + epsilon
kappa <- runif(p, min=minKappa, max=maxKappa)
Y[matrix(rbinom(n * p, size = 1, prob = rep(kappa, n)), nrow = p, ncol = n)==1] = 0
Y    <- t(Y)
X    <- t(X)
###############################################################################
###############################################################################

test_that("NB_unknown: check dimensions, optimization and field access", {
  model <- NB_unknown_zi$new(Y, X, c(3,6,4,5), niter = 60)
  model$optimize()
  best_model <- model$getBestModel("BIC")
  true_model <- model$get_model(Q)
  expect_equal(true_model$Q, Q)
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(best_model$BIC, 105035)
  expect_gt(true_model$loglik, -22167)
  model_sparse <- NB_unknown_zi$new(Y, X, c(3,6,4,5), c(0.01, 0.04, 0.02, 0.03), niter = 60)
  model_sparse$optimize()
})
