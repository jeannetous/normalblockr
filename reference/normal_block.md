# Normal-block model

Fit a normal-block model with a variational or heuristic algorithm

## Usage

``` r
normal_block(
  data,
  blocks,
  sparsity = 0,
  zero_inflation = FALSE,
  control = NB_control(),
  model = c("var", "mean")
)
```

## Arguments

- data:

  NormalBlockData object, contains the matrix of responses (Y, n x p)
  and the design matrix (X, n x d), must be created with
  NormalBlockData\$new.

- blocks:

  either a integer (number of blocks), a vector of integer (list of
  possible number of block) or a p \* q matrix (for indicating block
  membership when its known)

- sparsity:

  either TRUE to run the optimization for different sparsity penalty
  values OR float to run model with a single sparsity penalty value

- zero_inflation:

  boolean to indicate if Y is zero-inflated and adjust fitted model as a
  consequence

- control:

  a list-like structure for detailed control on parameters should be
  generated with NB_control().

- model:

  which model family to fit: "var" (the default) structures the
  clustering in the latent covariance (see \[NormalBlockVarBase\]),
  "mean" structures it in the mean (mu_i = C B' X_i, see
  \[NormalBlockMeanBase\]). The mean-block family covers the same
  collections as the variance-block one
  (\[NormalBlockMeanCollectionClusters\],
  \[NormalBlockMeanCollectionSparsity\] and their crossing,
  \[NormalBlockMeanCollectionClustersSparsity\]) – each q simply fit
  independently, without the variance-block family's SBM-path shortcut.
  Zero-inflation is supported (\[ZINormalBlockMeanKnownClusters\],
  \[ZINormalBlockMeanUnknownClusters\], and collections of those over a
  range of q), with a diagonal or spherical Sigma only, hence not along
  a sparsity path – see \[NormalBlockMeanBase\].

## Value

an R6 object with one of the model classes (or a collection of model
objects).

## Examples

``` r
## Normal Data
ex_data <- generate_normal_block_var_data(n=100, p=30, d=1, q=3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
my_normal_block <- normal_block(data, blocks = 1:6)
#> Fitting a diagonal normal-block-var model with unknown q 
#>   number of blocks = 1                number of blocks = 2                number of blocks = 3                number of blocks = 4                number of blocks = 5                number of blocks = 6           
#> DONE
my_normal_block$plot(c("deviance", "BIC", "ICL"))

Y_hat <- my_normal_block$get_best_model()$fitted
plot(data$Y, Y_hat, log = "xy"); abline(0,1)
#> Warning: 445 x values <= 0 omitted from logarithmic plot
#> Warning: 374 y values <= 0 omitted from logarithmic plot

## Normal Data with Zero Inflation
ex_data_zi <- generate_normal_block_var_data(n=50, p=50, d=1, q=3, kappa = rep(0.5,50))
zidata <- NormalBlockData$new(ex_data_zi$Y, ex_data_zi$X)
my_normal_block <- normal_block(zidata, blocks = 1:6, zero_inflation = TRUE)
#> Fitting a diagonal normal-block-var model with unknown q 
#>   number of blocks = 1                number of blocks = 2                number of blocks = 3                number of blocks = 4                number of blocks = 5                number of blocks = 6           
#> DONE
## Mean-Block model (clustering in the mean rather than the covariance)
ex_mean <- generate_normal_block_mean_data(n=50, p=20, d=1, q=3)
mean_data <- NormalBlockData$new(ex_mean$Y, ex_mean$X)
my_mean_block <- normal_block(mean_data, blocks = 3, model = "mean")
#> Fitting a diagonal normal-block-mean model with 3 unknown blocks 
#> 
#> DONE
```
