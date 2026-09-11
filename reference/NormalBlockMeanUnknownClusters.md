# Mean-Block Model with Unknown Clustering

R6 class for a Normal-Block-Mean model with a fixed number of clusters
(but unknown clustering), inferred by variational EM.

## Super classes

[`NormalBlockBase`](NormalBlockBase.md) -\>
[`NormalBlockMeanBase`](NormalBlockMeanBase.md) -\>
`NormalBlockMeanUnknownClusters`

## Active bindings

- `fitted`:

  Y values predicted by the model, in Y's original units

- `var_par`:

  a list with the variational parameter: tau (posterior group
  probabilities)

- `nb_param`:

  number of parameters in the model

- `entropy`:

  Entropy of the conditional distribution The only latent variable is
  the clustering, so the entropy of the variational distribution reduces
  to -sum(tau \* log(tau)).

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockMeanUnknownClusters$new()`](#method-NormalBlockMeanUnknownClusters-initialize)

- [`NormalBlockMeanUnknownClusters$clone()`](#method-NormalBlockMeanUnknownClusters-clone)

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

### `NormalBlockMeanUnknownClusters$new()`

Create a new \[\`NormalBlockMeanUnknownClusters\`\] object.

#### Usage

    NormalBlockMeanUnknownClusters$new(
      data,
      q,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `data`:

  object of NormalBlockData class, with responses and design matrix

- `q`:

  number of clusters

- `sparsity`:

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`NormalBlockMeanUnknownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanUnknownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanUnknownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
model <- normal_block(data, blocks = 3, model = "mean")
#> Fitting a diagonal normal-block-mean model with 3 unknown blocks 
#> 
#> DONE
model$clustering
#>  [1] 3 2 2 2 2 3 1 2 2 1 2 2 2 2 1 3 1 3 1 3
```
