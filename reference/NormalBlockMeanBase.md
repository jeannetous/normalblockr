# Base Class for Mean-Block Models

R6 abstract class for the Normal-Block models where the clustering
structures the mean (mu_i = C B' X_i).

## Super class

[`NormalBlockBase`](NormalBlockBase.md) -\> `NormalBlockMeanBase`

## Active bindings

- `model_par`:

  a list with the matrices of the model parameters: B (covariates), dm1
  (species variance), Omega (groups precision matrix)). On the internal
  fitting scale (\`self\$data\$Y\`, possibly column-rescaled by
  \`NormalBlockData(scale = TRUE)\`): use
  \`\$B_original\`/\`\$dm1_original\` for the same quantities converted
  back to Y's original units.

- `nb_param`:

  number of parameters in the model

- `sparsity_weights`:

  (weights associated to each pair of groups)

## Methods

### Public methods

- [`NormalBlockMeanBase$new()`](#method-NormalBlockMeanBase-initialize)

- [`NormalBlockMeanBase$predict()`](#method-NormalBlockMeanBase-predict)

- [`NormalBlockMeanBase$warm_start_from()`](#method-NormalBlockMeanBase-warm_start_from)

- [`NormalBlockMeanBase$split()`](#method-NormalBlockMeanBase-split)

- [`NormalBlockMeanBase$merge()`](#method-NormalBlockMeanBase-merge)

- [`NormalBlockMeanBase$clone()`](#method-NormalBlockMeanBase-clone)

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

------------------------------------------------------------------------

### `NormalBlockMeanBase$new()`

Create a new \[\`NormalBlockMeanBase\`\] object.

#### Usage

    NormalBlockMeanBase$new(
      data,
      q,
      sparsity = 0,
      control = NB_control(),
      zero_inflation = FALSE
    )

#### Arguments

- `data`:

  object of NormalBlockData class, with responses and design matrix

- `q`:

  number of block/cluster

- `sparsity`:

  sparsity penalty on the network density

- `control`:

  structured list of more specific parameters, to generate with
  NB_Mean_control

- `zero_inflation`:

  whether the concrete subclass models zero-inflation; set by the ZI
  subclasses themselves, not meant to be set by the end user.

#### Returns

A new \[\`NormalBlockMeanBase\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanBase$predict()`

Predicts observations Y for new covariates X, in Y's original units. The
mean-block mean is mu_i = C B' X_i, so the cluster-level predictor has
to be mapped back to the variables through C.

#### Usage

    NormalBlockMeanBase$predict(new_X)

#### Arguments

- `new_X`:

  new set of covariates.

#### Returns

A n\*p prediction matrix for new observations

------------------------------------------------------------------------

### `NormalBlockMeanBase$warm_start_from()`

Seed this model's starting parameters from another, already-optimized
model with the same q, instead of the heuristic clustering-derived
values set at construction time. Used by \[split()\]/\[merge()\].

#### Usage

    NormalBlockMeanBase$warm_start_from(other)

#### Arguments

- `other`:

  a \[NormalBlockMeanBase\] object, already optimized

#### Returns

Update the current object in place with \`other\`'s parameters

------------------------------------------------------------------------

### `NormalBlockMeanBase$split()`

Create a clone of the current \[\`NormalBlockMeanBase\`\] object after
splitting cluster \`index\`. Unlike the variance-block family, Omega and
the sparsity weights are p x p here and do not depend on q, so they
carry over unchanged; only C (tau) and B (one column per cluster) are
affected. Variables are split by their current noise variance (1 /
diag(Omega)) around its within-cluster median, the same criterion
\[NormalBlockVarBase\]'s \`split()\` uses via \`dm1\`, since
\`diag(Omega)\` plays the same per-variable-precision role here.

#### Usage

    NormalBlockMeanBase$split(index, in_place = FALSE)

#### Arguments

- `index`:

  index (integer) of the cluster to split

- `in_place`:

  should the split be applied to the object itself, or should a copy be
  sent? Default FALSE (send a copy)

#### Returns

A new \[\`NormalBlockMeanBase\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanBase$merge()`

Create a clone of the current \[\`NormalBlockMeanBase\`\] object after
merging clusters \`indices\`

#### Usage

    NormalBlockMeanBase$merge(indices, in_place = FALSE)

#### Arguments

- `indices`:

  indices (couple of integer) of the clusters to merge

- `in_place`:

  should the merge be applied to the object itself, or should a copy be
  sent? Default FALSE (send a copy)

#### Returns

A new \[\`NormalBlockMeanBase\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanBase$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# An internal abstract base class, never instantiated directly -- use
# NormalBlockMeanKnownClusters / NormalBlockMeanUnknownClusters.
```
