% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normal_block.R
\name{normal_block}
\alias{normal_block}
\title{Normal-block model}
\usage{
normal_block(
  data,
  blocks,
  sparsity = 0,
  zero_inflation = FALSE,
  control = NB_control()
)
}
\arguments{
\item{data}{contains the matrix of responses (Y, n x p) and the design matrix (X, n x d)."}

\item{blocks}{either a integer (number of blocks), a vector of integer (list of possible number of block)
or a p * Q matrix (for indicating block membership when its known)}

\item{sparsity}{either TRUE to run the optimization for different sparsity penalty values
OR float to run model with a single sparsity penalty value}

\item{zero_inflation}{boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence}

\item{control}{a list-like structure for detailed control on parameters should be
generated with NB_control().}
}
\value{
an R6 object with one of the NB classes
}
\description{
Fit a normal-block model with a variational or heuristic algorithm
}
\examples{
## Normal Data
ex_data <- generate_normal_block_data(n=50, p=50, d=1, Q=3)
data <- NBData$new(ex_data$Y, ex_data$X)
my_normal_block <- normal_block(data, blocks = 1:6)
\dontrun{
my_normal_block$plot(c("deviance", "BIC", "ICL"))
Y_hat <- my_normal_block$get_best_model()$fitted
plot(data$Y, Y_hat, log = "xy"); abline(0,1)
}
## Normal Data with Zero Inflation
ex_data_zi <- generate_normal_block_data(n=50, p=50, d=1, Q=3, kappa = rep(0.5,50))
zidata <- NBData$new(ex_data_zi$Y, ex_data_zi$X)
my_normal_block <- normal_block(zidata, blocks = 1:6, zero_inflation = TRUE)

}
