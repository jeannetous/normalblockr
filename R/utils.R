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

#' computes xlogx, setting it to 0 if x = 0
xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x*log(x))

#' computes softmax
softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

#' computes ARI between two clusterings
matchingGroupScores <- function(groups1, groups2){
  ari <- pdfCluster::adj.rand.index(groups1, groups2)
  # If ari is na, we want to see if a relabeling can change that
  if (is.na(ari)){if((length(unique(groups1)) == 1) & (length(unique(groups2)) == 1)){
    return(1)}}
  return(ari)
}
