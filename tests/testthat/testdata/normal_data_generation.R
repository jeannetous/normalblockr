set.seed(3)
n     <- 200
p     <- 50
d     <- 2
Q     <- 4
Sigma <- matrix(0, Q, Q)
diag(Sigma) <- 1
Sigma[row(Sigma) != col(Sigma)] <- 0.1
minB  <- 2
maxB  <- 4
minX  <- rep(0, d)
maxX  <- rep(10, d)
minD  <- 0.2
maxD  <- 3

B <- matrix(rep(1, d * p), nrow = d)
for (dim in 1:d) B[dim, ] <- runif(p, min = 0, max = 1)
X <- matrix(rep(1, d * n), nrow = d)
for (dim in 1:d) X[dim, ] <- runif(n, min = minX[[dim]], max = maxX[[dim]])
C <- matrix(rep(0, p * Q), nrow = p)
groups <- sample(1 : Q, size = p, replace = TRUE)
for (dim in 1:p) C[dim, groups[[dim]]] <- 1
D <- diag(runif(p, min = minD, max = maxD))
W <- t(MASS::mvrnorm(n, mu = matrix(rep(0, Q), Q, 1), Sigma = Sigma))
epsilon <- t(MASS::mvrnorm(n, mu = matrix(rep(0, p), p, 1), Sigma = D))
Y    <- t(B) %*% X + C %*% W + epsilon
Y    <- t(Y)
X    <- t(X)
testdata <- list(Y = Y, X = X, C = C, B = B, D = D, Sigma = Sigma)
saveRDS(testdata, file = "tests/testthat/testdata/testdata_normal.RDS")
