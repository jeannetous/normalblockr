# Mean-Block Model with Known Clustering

R6 class for a Normal-Block-Mean model with a known clustering.

## Super classes

[`NormalBlockBase`](NormalBlockBase.md) -\>
[`NormalBlockMeanBase`](NormalBlockMeanBase.md) -\>
`NormalBlockMeanKnownClusters`

## Active bindings

- `fitted`:

  Y values predicted by the model, in Y's original units

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockMeanKnownClusters$new()`](#method-NormalBlockMeanKnownClusters-initialize)

- [`NormalBlockMeanKnownClusters$clone()`](#method-NormalBlockMeanKnownClusters-clone)

Inherited methods

- [`NormalBlockBase$best_of_inits()`](NormalBlockBase.html#method-best_of_inits)
- [`NormalBlockBase$candidates_merge()`](NormalBlockBase.html#method-candidates_merge)
- [`NormalBlockBase$candidates_split()`](NormalBlockBase.html#method-candidates_split)
- [`NormalBlockBase$latent_network()`](NormalBlockBase.html#method-latent_network)
- [`NormalBlockBase$optimize()`](NormalBlockBase.html#method-optimize)
- [`NormalBlockBase$plot()`](NormalBlockBase.html#method-plot)
- [`NormalBlockBase$plot_loglik()`](NormalBlockBase.html#method-plot_loglik)
- [`NormalBlockBase$plot_network()`](NormalBlockBase.html#method-plot_network)
- [`NormalBlockBase$print()`](NormalBlockBase.html#method-print)
- [`NormalBlockBase$update()`](NormalBlockBase.html#method-update)
- [`NormalBlockMeanBase$merge()`](NormalBlockMeanBase.html#method-merge)
- [`NormalBlockMeanBase$predict()`](NormalBlockMeanBase.html#method-predict)
- [`NormalBlockMeanBase$split()`](NormalBlockMeanBase.html#method-split)
- [`NormalBlockMeanBase$warm_start_from()`](NormalBlockMeanBase.html#method-warm_start_from)

------------------------------------------------------------------------

### `NormalBlockMeanKnownClusters$new()`

Create a new \[\`NormalBlockMeanKnownClusters\`\] object.

#### Usage

    NormalBlockMeanKnownClusters$new(data, C, sparsity = 0, control = NB_control())

#### Arguments

- `data`:

  object of NormalBlockData class, with responses and design matrix

- `C`:

  clustering matrix C_jk = 1 if species j belongs to cluster k

- `sparsity`:

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`NormalBlockMeanKnownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanKnownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanKnownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
model <- normal_block(data, blocks = ex$parameters$C, model = "mean")
#> Fitting a diagonal normal-block-mean model with fixed blocks 
#> 
#> DONE
model$plot()
```
