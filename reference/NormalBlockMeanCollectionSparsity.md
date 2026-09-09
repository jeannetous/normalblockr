# Collection of Mean-Block Models over a Sparsity Path

R6 class for a collection of mean-block models (\[NormalBlockMeanBase\])
with a fixed clustering (or a fixed number of blocks) and different
sparsity levels applied to the p x p precision matrix of the variables.
Mirrors \[NormalBlockVarCollectionSparsity\], minus the StARS/stability
selection path, which relies on \`fixed_tau\`, not supported by the
mean-block VEM.

## Super classes

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
[`NormalBlockCollectionSparsity`](NormalBlockCollectionSparsity.md) -\>
`NormalBlockMeanCollectionSparsity`

## Active bindings

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockMeanCollectionSparsity$new()`](#method-NormalBlockMeanCollectionSparsity-initialize)

- [`NormalBlockMeanCollectionSparsity$get_best_model()`](#method-NormalBlockMeanCollectionSparsity-get_best_model)

- [`NormalBlockMeanCollectionSparsity$clone()`](#method-NormalBlockMeanCollectionSparsity-clone)

Inherited methods

- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionSparsity$get_model()`](NormalBlockCollectionSparsity.html#method-get_model)
- [`NormalBlockCollectionSparsity$optimize()`](NormalBlockCollectionSparsity.html#method-optimize)
- [`NormalBlockCollectionSparsity$plot()`](NormalBlockCollectionSparsity.html#method-plot)

------------------------------------------------------------------------

### `NormalBlockMeanCollectionSparsity$new()`

Create a new \[\`NormalBlockMeanCollectionSparsity\`\] object.

#### Usage

    NormalBlockMeanCollectionSparsity$new(mydata, blocks, control = NB_control())

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `blocks`:

  either a clustering matrix (known, fixed clustering) or a single
  integer (number of blocks to infer)

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new \[\`NormalBlockMeanCollectionSparsity\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanCollectionSparsity$get_best_model()`

Extract best model in the collection

#### Usage

    NormalBlockMeanCollectionSparsity$get_best_model(
      crit = c("BIC", "EBIC", "ICL")
    )

#### Arguments

- `crit`:

  a character for the criterion used to perform the selection, either
  "BIC", "EBIC" or "ICL". Default is BIC

#### Returns

a \[\`NormalBlockMeanBase\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanCollectionSparsity$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanCollectionSparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_mean_data(n = 60, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
models <- normal_block(data, blocks = 3, sparsity = TRUE, model = "mean",
                       control = NB_control(n_sparsity_penalties = 5))
#> Fitting a normal-block-mean model with sparsity path 
#>   penalty = 0.3489394             penalty = 0.1103443             penalty = 0.03489394                penalty = 0.01103443                penalty = 0.003489394           
#> DONE
models$plot(c("BIC", "EBIC"))
```
