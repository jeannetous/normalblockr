###############################################################################
###############################################################################
set.seed(3)
n     <-  300
p     <-  150
d     <-  2
Q     <-  4
Sigma <- matrix(0, Q, Q)
diag(Sigma) <- 1
Sigma[row(Sigma) != col(Sigma)] = 0.1
minB  <-  1
maxB  <-  2
minX  <-  rep(1, d)
maxX  <-  rep(10, d)
minD  <-  0.2
maxD  <-  3
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

Y    <-  t(Y)
X    <-  t(X)
###############################################################################
###############################################################################

test_that("NB_fixed_blocks_zi: check dimensions, optimization and field access", {
  model <- NB_fixed_blocks_zi$new(Y, X, C, niter = 60)
  model$optimize()
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 161095)
  expect_gt(model$loglik, -79236)
  model_sparse <- NB_fixed_blocks_zi$new(Y, X, C, sparsity = 0.05, niter = 60)
  model_sparse$optimize()
})
