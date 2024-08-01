#'turns a list clustering of Q cluster labels for N elements into a matrix of
#'dimensions (N, Q) with a one-hot encoding of the clustering
#'
#'@param clustering a list of labels
#'
#'@examples
#'as_indicator(c(1,1,2,3,2))
as_indicator <- function(clustering) {
  Q <- max(clustering)
  N <- length(clustering)
  Z <- matrix(0, N, Q)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}

#' removes machine's 0 to elements equal to  1 in x
check_one_boundary <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x >= 1 - zero] <- 1 - zero
  x
}

#' adds machine's 0 to elements equal to 0 in x
check_zero_boundary <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x < zero]  <- zero
  x
}
