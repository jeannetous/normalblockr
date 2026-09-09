# Summarize a Collection of Normal-Block Models

Summarizes a collection of normal-block models: model type, the full
criteria table, and the range of q/sparsity explored.

## Usage

``` r
# S3 method for class 'NormalBlockCollection'
summary(object, ...)
```

## Arguments

- object:

  An object inheriting from NormalBlockCollection.

- ...:

  not used, only here for S3 compatibility

## Value

An object of class \`summary.NormalBlockCollection\` (a list with the
collection's \`who_am_I\`, \`criteria\`, \`q_range\` and
\`sparsity_range\`), printed with a dedicated
\[print.summary.NormalBlockCollection()\] method.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
summary(models)
#> A diagonal normal-block-var model with unknown q 
#> ===========================================================================
#>     q ranging from 2 to 5 
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        44 2       1        0 -785.297 1570.595 1742.724 1601.310 1744.110    16
#>        48 3       3        0 -529.095 1058.190 1245.967 1002.366 1252.559     9
#>        53 4       6        0 -525.888 1051.775 1259.112  957.837 1275.748    92
#>        59 5      10        0 -525.877 1051.753 1282.562 1042.533 1314.751   146
```
