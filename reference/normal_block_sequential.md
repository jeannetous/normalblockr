# Cluster variables in the mean, then in the residual covariance

Fits the two model families one after the other: a mean-block model
(\[NormalBlockMeanBase\]) groups the variables by how they respond to
the covariates, then a variance-block model (\[NormalBlockVarBase\])
groups the residuals of that fit by how they co-vary. The two answer
different questions and generally return unrelated partitions, so
running both is often more informative than choosing one.

## Usage

``` r
normal_block_sequential(
  data,
  blocks_mean,
  blocks_var,
  crit = c("ICL", "BIC"),
  zero_inflation = FALSE,
  control_mean = NB_control(verbose = FALSE),
  control_var = NB_control(verbose = FALSE)
)
```

## Arguments

- data:

  a \[NormalBlockData\] object

- blocks_mean:

  number of clusters for the mean-block stage: an integer, a vector of
  integers to explore, or a p x q indicator matrix

- blocks_var:

  idem for the variance-block stage, run on the residuals

- crit:

  criterion used to pick a model when a range is explored, "ICL" (the
  default) or "BIC"

- zero_inflation:

  whether Y carries structural zeros. Both stages are then
  zero-inflated, sharing the same mask: the second one would otherwise
  re-derive it from residuals, which are never exactly zero.

- control_mean:

  control list for the mean-block stage, see \[NB_control()\]

- control_var:

  control list for the variance-block stage

## Value

an object of class \`normal_block_sequential\`, a list with the fitted
\`mean\` and \`var\` models and the residual matrix \`residuals\` handed
from one stage to the other.

## Details

The second stage uses an intercept-only design on purpose: the covariate
effects have already been removed by the first stage.

This is a heuristic two-stage estimator, not a joint model. On simulated
data carrying two genuinely distinct structures it recovers both exactly
(see \`inst/mean_block_analyses/sequential_mean_then_variance.R\`).

## Examples

``` r
ex   <- generate_normal_block_mean_data(n = 80, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
fit  <- normal_block_sequential(data, blocks_mean = 3, blocks_var = 2)
fit
#> A sequential mean-then-variance normal-block fit
#> ===========================================================================
#>   mean-block stage    : 3 clusters -- diagonal normal-block-mean model with 3 unknown blocks 
#>   variance-block stage: 2 clusters -- diagonal normal-block-var model with 2 unknown blocks 
#>   ARI between the two partitions: -0.056 
#>   (near 0 means the two structures are unrelated, as is usual)
#> ===========================================================================
#> * Useful fields
#>     $mean, $var (both ordinary fitted models), $residuals 
```
