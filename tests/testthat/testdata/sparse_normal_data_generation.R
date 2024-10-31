set.seed(3)
n     <- 300
p     <- 100
d     <- 2
Q     <- 4
Omega <- matrix(0, Q, Q); diag(Omega) <- 1
Omega[1,2] <-  -0.08 ; Omega[2, 1] <- -0.08
Omega[1,3] <-  -0.08 ; Omega[3, 1] <-  -0.08
Omega[2,4] <-  -0.08 ; Omega[4, 2] <-  -0.08
Sigma <- solve(Omega)
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
saveRDS(testdata, file = "tests/testthat/testdata/testdata_sparse_normal.RDS")
