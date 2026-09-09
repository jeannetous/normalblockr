# Collection of Mean-Block Models over a Range of Cluster Counts

R6 class for a collection of mean-block models (\[NormalBlockMeanBase\])
with different numbers of clusters (q). Inherits its scaffolding
(\`print()\`/\`summary()\`/\`plot()\`/\`optimize()\`, the \`criteria\`
table) from \[NormalBlockCollection\]. Unlike
\[NormalBlockVarCollectionClusters\], there is no SBM-path shortcut
here: the shared clustering-heuristic registry's cov()/correlation-based
methods are ill-suited to the mean-block family (see
\[NormalBlockMeanBase\]'s own default), so each q is fit independently.

## Super classes

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
[`NormalBlockCollectionClusters`](NormalBlockCollectionClusters.md) -\>
`NormalBlockMeanCollectionClusters`

## Active bindings

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockMeanCollectionClusters$new()`](#method-NormalBlockMeanCollectionClusters-initialize)

- [`NormalBlockMeanCollectionClusters$clone()`](#method-NormalBlockMeanCollectionClusters-clone)

Inherited methods

- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionClusters$get_best_model()`](NormalBlockCollectionClusters.html#method-get_best_model)
- [`NormalBlockCollectionClusters$get_model()`](NormalBlockCollectionClusters.html#method-get_model)
- [`NormalBlockCollectionClusters$optimize()`](NormalBlockCollectionClusters.html#method-optimize)
- [`NormalBlockCollectionClusters$plot()`](NormalBlockCollectionClusters.html#method-plot)
- [`NormalBlockCollectionClusters$refine()`](NormalBlockCollectionClusters.html#method-refine)

------------------------------------------------------------------------

### `NormalBlockMeanCollectionClusters$new()`

Create a new \[\`NormalBlockMeanCollectionClusters\`\] object.

#### Usage

    NormalBlockMeanCollectionClusters$new(
      mydata,
      q_list,
      zero_inflation = FALSE,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `q_list`:

  list of q values (number of groups) in the collection

- `zero_inflation`:

  whether Y carries structural zeros; every model in the collection is
  then zero-inflated (Sigma diagonal or spherical only, see
  \[NormalBlockMeanBase\])

- `sparsity`:

  sparsity penalty on the network density

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`NormalBlockMeanCollectionClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanCollectionClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanCollectionClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
models <- normal_block(data, blocks = 2:5, model = "mean")
#> Fitting a normal-block-mean model with unknown q 
#>   number of blocks = 2                number of blocks = 3                number of blocks = 4                number of blocks = 5           
#> DONE
models$plot(c("BIC", "ICL"))

models$get_best_model()
#> A diagonal normal-block-mean model with 2 unknown blocks .
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        23 2       0        0 -963.653 1927.305 2017.282 2017.282 2017.282     3
#> ===========================================================================
#> * Useful fields
#>     $model_par, $posterior_par / $var_par, $clustering 
#>     $loglik, $BIC, $ICL, $objective, $nb_param, $criteria
#> * Useful S3 methods
#>     print(), summary(), plot(), coef(), sigma(), fitted(), predict() 
```
