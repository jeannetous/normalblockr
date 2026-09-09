# Package index

## User-facing functions

Main entry points for fitting a Normal-Block model and simulating data

- [`normal_block()`](normal_block.md) : Normal-block model
- [`normal_block_sequential()`](normal_block_sequential.md) : Cluster
  variables in the mean, then in the residual covariance
- [`NB_control()`](NB_control.md) : NB_control
- [`generate_normal_block_var_data()`](generate_normal_block_var_data.md)
  : Generate Normal Block Var Data
- [`generate_normal_block_mean_data()`](generate_normal_block_mean_data.md)
  : Generate Normal Block Mean Data

## Model and data classes

R6 classes returned by normal_block(). NormalBlockData wraps the
responses and design matrix; NormalBlockVar*Clusters (optionally
ZI-prefixed) are single fitted models, known- or unknown-clustering;
NormalBlockCollection* are collections of models explored over a range
of q and/or sparsity penalties.

- [`NormalBlockData`](NormalBlockData.md) : Data Container for
  Normal-Block Models
- [`NormalBlockVarKnownClusters`](NormalBlockVarKnownClusters.md) :
  Normal-Block Model with Known Clustering
- [`NormalBlockVarUnknownClusters`](NormalBlockVarUnknownClusters.md) :
  Normal-Block Model with Unknown Clustering
- [`NormalBlockMeanKnownClusters`](NormalBlockMeanKnownClusters.md) :
  Mean-Block Model with Known Clustering
- [`NormalBlockMeanUnknownClusters`](NormalBlockMeanUnknownClusters.md)
  : Mean-Block Model with Unknown Clustering
- [`ZINormalBlockVarKnownClusters`](ZINormalBlockVarKnownClusters.md) :
  Zero-Inflated Normal-Block Model with Known Clustering
- [`ZINormalBlockVarUnknownClusters`](ZINormalBlockVarUnknownClusters.md)
  : Zero-Inflated Normal-Block Model with Unknown Clustering
- [`ZINormalBlockMeanKnownClusters`](ZINormalBlockMeanKnownClusters.md)
  : Zero-Inflated Mean-Block Model with Known Clustering
- [`ZINormalBlockMeanUnknownClusters`](ZINormalBlockMeanUnknownClusters.md)
  : Zero-Inflated Mean-Block Model with Unknown Clustering
- [`NormalBlockVarCollectionClusters`](NormalBlockVarCollectionClusters.md)
  : Collection of Normal-Block Models over a Range of Cluster Counts
- [`NormalBlockVarCollectionSparsity`](NormalBlockVarCollectionSparsity.md)
  : Collection of Normal-Block Models over a Sparsity Path
- [`NormalBlockVarCollectionClustersSparsity`](NormalBlockVarCollectionClustersSparsity.md)
  : Collection of Normal-Block Models over Cluster Counts and Sparsity
  Levels
- [`NormalBlockMeanCollectionClusters`](NormalBlockMeanCollectionClusters.md)
  : Collection of Mean-Block Models over a Range of Cluster Counts
- [`NormalBlockMeanCollectionSparsity`](NormalBlockMeanCollectionSparsity.md)
  : Collection of Mean-Block Models over a Sparsity Path
- [`NormalBlockMeanCollectionClustersSparsity`](NormalBlockMeanCollectionClustersSparsity.md)
  : Collection of Mean-Block Models over Cluster Counts and Sparsity
  Levels

## Data sets

- [`brca_rppa`](brca_rppa.md) : Breast cancer proteomics data (TCGA,
  RPPA)
- [`onema`](onema.md) : French stream fish community data (ONEMA / OFB
  electrofishing surveys)
- [`university`](university.md) : University webpages text data (CMU "4
  Universities" / WebKB)

## S3 methods

Standard extractors and methods for a fitted model (any NormalBlockBase
subclass) and for a collection of models (any NormalBlockCollection
subclass)

- [`coef(`*`<NormalBlockBase>`*`)`](coef.NormalBlockBase.md) : Extract
  Model Coefficients
- [`fitted(`*`<NormalBlockBase>`*`)`](fitted.NormalBlockBase.md) :
  Extract Fitted Values
- [`predict(`*`<NormalBlockBase>`*`)`](predict.NormalBlockBase.md) :
  Predict Method for Variance-Block Models
- [`sigma(`*`<NormalBlockBase>`*`)`](sigma.NormalBlockBase.md) : Extract
  the Covariance Matrix
- [`print(`*`<NormalBlockBase>`*`)`](print.NormalBlockBase.md) : Print a
  Normal-Block Model
- [`summary(`*`<NormalBlockBase>`*`)`](summary.NormalBlockBase.md) :
  Summarize a Normal-Block Model
- [`print(`*`<summary.NormalBlockBase>`*`)`](print.summary.NormalBlockBase.md)
  : Print a Normal-Block Model Summary
- [`plot(`*`<NormalBlockBase>`*`)`](plot.NormalBlockBase.md) : Plot a
  Normal-Block Model
- [`logLik(`*`<NormalBlockBase>`*`)`](logLik.NormalBlockBase.md) :
  Extract Log-Likelihood of a Normal-Block Model
- [`BIC(`*`<NormalBlockBase>`*`)`](BIC.NormalBlockBase.md) : Bayesian
  Information Criterion for a Normal-Block Model
- [`print(`*`<NormalBlockCollection>`*`)`](print.NormalBlockCollection.md)
  : Print a Collection of Normal-Block Models
- [`summary(`*`<NormalBlockCollection>`*`)`](summary.NormalBlockCollection.md)
  : Summarize a Collection of Normal-Block Models
- [`print(`*`<summary.NormalBlockCollection>`*`)`](print.summary.NormalBlockCollection.md)
  : Print a Collection Summary
- [`logLik(`*`<NormalBlockCollection>`*`)`](logLik.NormalBlockCollection.md)
  : Extract Log-Likelihood of a Collection of Normal-Block Models
- [`BIC(`*`<NormalBlockCollection>`*`)`](BIC.NormalBlockCollection.md) :
  Bayesian Information Criterion for a Collection of Normal-Block Models
- [`print(`*`<normal_block_sequential>`*`)`](print.normal_block_sequential.md)
  : Print a Sequential Mean-then-Variance Fit

## Internal

Abstract base classes, low-level helpers and exploratory code kept for
reference; not part of the user-facing API

- [`NormalBlockBase`](NormalBlockBase.md) : Root Base Class for
  Normal-Block Models
- [`NormalBlockVarBase`](NormalBlockVarBase.md) : Base Class for
  Variance-Block Models
- [`NormalBlockMeanBase`](NormalBlockMeanBase.md) : Base Class for
  Mean-Block Models
- [`NormalBlockCollection`](NormalBlockCollection.md) : Base Class for a
  Collection of Normal-Block Models
- [`NormalBlockCollectionClusters`](NormalBlockCollectionClusters.md) :
  Base Class for a Collection of Models over a Range of Cluster Counts
- [`NormalBlockCollectionSparsity`](NormalBlockCollectionSparsity.md) :
  Base Class for a Collection of Models over a Sparsity Path
- [`NormalBlockCollectionClustersSparsity`](NormalBlockCollectionClustersSparsity.md)
  : Base Class for a Collection over Cluster Counts and Sparsity Levels
- [`get_model()`](get_model.md) : Create a Normal-Block Model Object
- [`isNB()`](isNB.md) : Check if an Object is a Normal-Block Model
