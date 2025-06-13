# List of functions to generate the model parameters given a general parametrization

#' Generate D, diagonal variance matrix for individual variables
#' @param p number of variables
#' @param min_D minimum value taken by D's diagonal
#' @param max_D maximum value taken by D's diagonal
generate_D <- function(p, min_D = 0.5, max_D = 1.5){
  D <- matrix(rep(0, p*p), nrow = p)
  diag(D) <- runif(p, min_D, max_D)
  return(D)
}

#' Generate X, covariates matrix
#' @param n number of individuals
#' @param d number of covariates
#' @param min_X minimum value taken by X
#' @param max_X maximum value taken by X
generate_X <- function(n, d, min_X = 0, max_X = 10){
  X = matrix(rep(1, n * d), nrow=n)
  for(dim in 1:d){X[,dim] = runif(n, min=min_X[[dim]], max = max_X[[dim]])}
  return(X)
}

#' Generate B, regression matrix
#' @param p number of variables
#' @param X covariates matrix
#' @param Sigma variance-covariance matrix (useful to calibrate the SNR)
#' @param SNR Signal to Noise Ratio, ratio between the weight of Sigma and B
generate_B <- function(p, X, Sigma, SNR = 0.75){
  d <- ncol(X)
  B <- matrix(rep(1, d*p), nrow=d)
  for(dim in 1:d){B[dim,] = runif(p, min=0, max = 1)}
  correcting_factor <- SNR * var(as.vector(Sigma)) / (var(as.vector(X %*% B)))
  B <- sqrt(correcting_factor) * B
  return(B)
}

#' Function toenerate C, a clustering matrix
#' @param p number of variables
#' @param Q number of clusters
#' @param alpha probabilities to belong to each cluster (default is equiprobability)
generate_blocks <- function(p, Q, alpha = NULL){
  if(is.null(alpha)) alpha <- rep(1/Q, Q)
  C <- matrix(rep(0, p * Q), nrow = p)
  while(0 %in% colSums(C)){
    groups = sample(1 : Q, size = p, replace = TRUE)
    for(dim in 1:p){C[dim, groups[[dim]]] = 1}
  }
  return(C)
}

#' Function to generate the model's parametera B, C, D, Omega, Sigma, kappa
#' @param X covariates matrix
#' @param p number of variables
#' @param Q number of clusters
#' @param kappa list of zero-inflation probabilities for each variable
#' @param omega_structure either erdos-reyni, preferential attachment, community (SBM)
generate_parameters <- function(X, p, Q, kappa, omega_structure) {
  Omega <- generate_omega(Q, omega_structure)
  Sigma <- chol2inv(chol(Omega))
  list(
    B = generate_B(p, X, Sigma),
    C = generate_blocks(p, Q),
    D = generate_D(p),
    Omega = Omega,
    Sigma = Sigma,
    kappa = kappa
  )
}
