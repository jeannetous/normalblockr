generate_D <- function(p, min_D = 0.5, max_D = 1.5){
  runif(p, min_D, max_D)
}

generate_X <- function(n, d, min_X = 0, max_X = 10){
  matrix(runif(n*d, min=min_X, max = max_X), n, d)
}

generate_B <- function(p, X, Sigma, SNR = 0.75){
  d <- ncol(X)
  B <- matrix(runif(p*d), d, p)
  correcting_factor <- SNR * var(as.vector(Sigma)) / (var(as.vector(X %*% B)))
  B <- sqrt(correcting_factor) * B
  B
}

generate_blocks <- function(p, Q, alpha = rep(1/Q, Q)){
  C <- matrix(rep(0, p * Q), nrow = p)
  while(0 %in% colSums(C)){
    groups = sample(1 : Q, size = p, replace = TRUE)
    C[cbind(1:p, groups)] <- 1
  }
  C
}

# Erdos-renyi
erdos_reyni_graph <- function(Q, p = 0.5){
  as_adjacency_matrix(sample_gnp(Q, p))
}

# Preferential attachment
preferential_attachment_graph <- function(Q){
  as_adjacency_matrix(sample_pa(Q, directed = FALSE))
}

# Community structure
community_graph <- function(Q, prob = c(1/2,1/4,1/4), p_in = 0.5, p_out = 0.1) {
  pref_mat <- matrix(p_out, length(prob), length(prob))
  diag(pref_mat) <- p_in
  graph_mat <- as_adjacency_matrix(sample_sbm(Q,
                                              pref.matrix = pref_mat,
                                              block.sizes = c(rmultinom(1, Q, prob)) ))
  graph_mat
}

generate_precision_matrix <- function(Q, graph_structure = "erdos-renyi", v = 0.3, u = 0.1){
  cond <- FALSE
  while(!cond){

    if (is.matrix(graph_structure)) {
      stopifnot(isSymmetric.matrix(graph_structure))
      G <- 1 * (graph_structure != 0)
    } else {
      stopifnot(graph_structure %in% c("erdos-renyi", "preferential_attachment", "community"))
      G <- switch(graph_structure,
        "erdos-renyi" = erdos_reyni_graph(Q),
        "preferential_attachment" =  preferential_attachment_graph(Q),
        "community" = community_graph(Q)
      )
    }

    # Ensuring that the network is not empty
    if(max(G) == 0){
      off_diag_indices <- which(row(matrix(1:Q, Q, Q)) != col(matrix(1:Q, Q, Q)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      G[selected_index[["row"]], selected_index[["col"]]] <- 1
      G[selected_index[["col"]], selected_index[["row"]]] <- 1
    }
    omega_tilde <- G * v
    omega <- omega_tilde + diag(abs(min(eigen(omega_tilde)$values)) + u, Q, Q)

    # Ensuring that the network is not full for AUC to make sense
    if(min(omega) > 0){ # Ensuring that the network has 0s for AUC to make sense
      off_diag_indices <- which(row(matrix(1:Q, Q, Q)) != col(matrix(1:Q, Q, Q)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      omega[selected_index[["row"]], selected_index[["col"]]] <- 0
      omega[selected_index[["col"]], selected_index[["row"]]] <- 0
    }
    cond <- all(eigen(omega)$values > 0)
  }
  as.matrix(omega)
}

# Generate normal block set of model parameters
#
generate_normal_block_param <- function(X = matrix(rnorm(100*p), 100, p),
                                        p = 40,
                                        Q = 3,
                                        kappa = 0,
                                        omega_structure="erdos-renyi",
                                        alpha = rep(1/Q, Q),
                                        SNR = 0.75,
                                        range_D = c(0.5, 1.5)) {
  Omega <- generate_precision_matrix(Q, omega_structure)
  Sigma <- chol2inv(chol(Omega))
  list(
    B = generate_B(p, X, Sigma, SNR),
    C = generate_blocks(p, Q, alpha),
    D = generate_D(p, range_D[1], range_D[2]),
    Omega = Omega,
    Sigma = Sigma,
    kappa = kappa
  )
}

#' Generate Normal Block Data
#'
#' A function to draw data from the normal block model (see details). The function returns both the generated data and the corrresponding model parameters, in a list.
#'
#' @param n number of individuals. Default to 100.
#' @param p number of variables. Default to 40.
#' @param d number of covariates. Default to 1.
#' @param Q number of groups. Default to 3.
#' @param kappa vector (or scalar) of variable-wise probability of zero inflation. Default to 0.
#' @param omega_structure the structure of the graph on which the precision matrix between groups is built. Can be a symmetric matrix with Q rows/columns or a character picked in "erdos-renyi", "preferential_attachment", "community" in which case a graph is drawn with sensible generation parameters. See generate_precision_matrix for details.
#' @param range_X A 2-size vector defining the range of the uniform distribution used to draw values in X, the regressor matrix. Default is c(0, 10)
#' @param range_D A 2-size vector defining the range of the uniform distribution used to draw values in D, the diagonal matrix of variances of variables. Default is c(0.5, 1.5)
#' @param u_v two-size vector of positive numbers u and v controlling the generation of the precision matrix Omega: u is the off-diagonal elements of the precision matrix, controlling the magnitude of partial correlations with v a positive number being added to the diagonal elements of the precision matrix. The default value is c(0.3, 0.1).
#' @param alpha the Q-size vector of group proportion. Default to rep(1/Q, Q)
#' @param SNR Signal to noise ratio: magnitude of the regression parameters B will be adjusted so that tr(var(XB)) and tr(Sigma) match the desried SNR.
#'
#' @returns A named list with the following element
#' - Y a matrix of responses
#' - X a regressor/design matrix
#' - a list of model parameters, encompassing
#'    - B: matrix of regression coefficients
#'    - C: matrix of group membership
#'    - D: diagonal matrix of variance of the variables
#'    - Omega: precision matrix of the groups
#'    - Sigma: covariance matrix of the groups
#'    - kappa: vector of ZI inflation proabilities (one per variable)
#'
#' @importFrom igraph sample_pa sample_sbm sample_gnp as_adjacency_matrix
#' @importFrom MASS mvrnorm
#' @importFrom stats rbinom rmultinom rnorm runif var
#' @export
generate_normal_block_data <-
  function(n = 100,
           p = 40,
           d = 1,
           Q = 3,
           kappa = 0,
           omega_structure="erdos-renyi",
           u_v = c(0.3, 0.1),
           SNR = 0.75,
           alpha = rep(1/Q, Q),
           range_X = c(0, 10),
           range_D = c(0.5, 1.5)) {

  X <- matrix(runif(n*d, min=range_X[1], max = range_X[2]), n, d)
  param <- generate_normal_block_param(X, p, Q, kappa, omega_structure, alpha, SNR)

  W <- MASS::mvrnorm(n, mu = matrix(rep(0, Q), Q, 1), Sigma = param$Sigma)

  epsilon <- MASS::mvrnorm(n, mu = matrix(rep(0, p), p, 1), Sigma = diag(param$D))

  Y <- X  %*% param$B + tcrossprod(W, param$C) + epsilon

  if(any(param$kappa > 0)) {
    zero_location <- matrix(rbinom(n * p,  1, rep(param$kappa, each = n)), n, p)
    zero_location[, colSums(zero_location) == n] <- 0
    Y[zero_location == 1] <- 0
  }

  list(Y = Y, X = X, parameters = param)
}
