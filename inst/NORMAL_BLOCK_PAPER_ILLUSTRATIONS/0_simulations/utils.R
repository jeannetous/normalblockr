library(sbm)
library(ClustOfVar)

#' For a list of edges, give corresponding (node1, node2) list.
#' @param x list of edges
#' @param n number of nodes involved in the edges
edge_to_node <- function(x, n = max(x)) {
  x <- x - 1 ## easier for arithmetic to number edges starting from 0
  n.node <- round((1 + sqrt(1 + 8*n)) / 2) ## n.node * (n.node -1) / 2 = n (if integer)
  j.grid <- cumsum(0:n.node)
  j <- findInterval(x, vec = j.grid)
  i <- x - j.grid[j]
  ## Renumber i and j starting from 1 to stick with R convention
  return(data.frame(node1 = i + 1, node2 = j + 1))
}

#' Rewrite a clustering as an indicator matrix
#' @param clustering given as a list of labels
as_indicator <- function(clustering) {
  Q <- max(clustering)
  N <- length(clustering)
  Z <- matrix(0, N, Q)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}


#' Rewrite a clustering as a list of labels (reverse function of as_indicator)
#' @param tau clustering / clustering probabilities given as an indicator / probabilities matrix
get_clusters <- function(tau) {
  return(apply(tau, 1, which.max))
}
