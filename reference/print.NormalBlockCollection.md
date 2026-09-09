# Print a Collection of Normal-Block Models

Print a short summary of a collection of normal-block models: model type
and the range of q/sparsity explored. See
\[summary.NormalBlockCollection()\] for the full criteria table.

## Usage

``` r
# S3 method for class 'NormalBlockCollection'
print(x, ...)
```

## Arguments

- x:

  An object inheriting from NormalBlockCollection.

- ...:

  not used, only here for S3 compatibility

## Value

Invisibly returns \`x\`.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
print(models)
#> A diagonal normal-block-var model with unknown q 
#> ===========================================================================
#>   4 model(s) explored
#>     q ranging from 2 to 5 
#> ===========================================================================
#> * Useful fields
#>     $models, $criteria
#> * Useful methods
#>     print(), summary(), plot(), $get_best_model()
```
