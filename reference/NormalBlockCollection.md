# Base Class for a Collection of Normal-Block Models

Shared scaffolding for the collections explored by \[get_model()\]/
\[normal_block()\]: a sweep over sparsity penalties
(\[\`NormalBlockVarCollectionSparsity\`\]), over the number of clusters
(\[\`NormalBlockVarCollectionClusters\`\]), or over both
(\[\`NormalBlockVarCollectionClustersSparsity\`\]). Concrete subclasses
set \`private\$progress_field\`/\`private\$progress_label\` in their
\`initialize()\` and provide their own \`get_best_model()\`, delegating
the (row of \`self\$criteria\` minimizing a criterion) lookup to
\`private\$best_id()\`.

## Public fields

- `models`:

  list of models (or sub-collections) explored by the collection

- `control`:

  store the list of user-defined model settings and optimization
  parameters

## Active bindings

- `criteria`:

  a data frame with the values of some criteria for the collection of
  models

- `loglik`:

  not defined for a collection: accessing it raises an informative
  error. Use \`logLik()\` for every model's log-likelihood, or
  \`\$get_best_model()\$loglik\` for a single one.

## Methods

### Public methods

- [`NormalBlockCollection$optimize()`](#method-NormalBlockCollection-optimize)

- [`NormalBlockCollection$print()`](#method-NormalBlockCollection-print)

- [`NormalBlockCollection$summary()`](#method-NormalBlockCollection-summary)

- [`NormalBlockCollection$clone()`](#method-NormalBlockCollection-clone)

------------------------------------------------------------------------

### `NormalBlockCollection$optimize()`

optimizes every model (or sub-collection) in the collection

#### Usage

    NormalBlockCollection$optimize(
      control = list(niter = 500, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter and threshold). When
  \`control\$clustering_init\` is \`"best_of_inits"\`, each leaf model
  is fit via its own \`best_of_inits()\` instead of a plain
  \`optimize()\`.

------------------------------------------------------------------------

### `NormalBlockCollection$print()`

User-friendly print method: model type and the range of q/sparsity
explored. See \`summary()\` for the full criteria table.

#### Usage

    NormalBlockCollection$print()

------------------------------------------------------------------------

### `NormalBlockCollection$summary()`

Summarize the collection: model type, full criteria table, and the range
of q/sparsity explored.

#### Usage

    NormalBlockCollection$summary()

#### Returns

An object of class \`summary.NormalBlockCollection\`, printed with a
dedicated \[print.summary.NormalBlockCollection()\] method.

------------------------------------------------------------------------

### `NormalBlockCollection$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockCollection$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# An internal abstract base class, never instantiated directly -- see
# normal_block() for how collections (NormalBlockVarCollectionClusters,
# NormalBlockVarCollectionSparsity, NormalBlockVarCollectionClustersSparsity)
# are actually created and fitted.
```
