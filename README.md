# IncrementalRMA

## Overview

IncrementalRMA works by applying previous computed parameters from a full RMA process. This
function calculates the necessary parameters for future `incrementalRMA` use.

RMA consists of three steps: background correction, quantile normalization and median polish.
Background correction is done per-chip, so there is no need for pre-calculating parameters.
Quantile normalization computes values for each quantile in expression data and then sets the
intensity distribution equal to these quantiles. Median polish sweeps out row-wise and column-wise
medians from probe expression; column-wise effects are used as expression estimates.

Given the above description, quantiles (from quantile normalization) and row-wise effects
(from median polish) can be computed on a set of samples and stored. With this information,
`incrementalRMA` can then calculate gene expression from a new sample without
re-computing these parameters.

This function computes these parameters and stores them in a list, such that the list can then
be used. The elements of the list currently include:

- *probeEffects*: Row-wise effects from median polish step.
- *normalizationVector*: Quantile values from quantile normalization step.
- *referenceCELFiles*: An `AffyBatch` containing the reference cel
   files (for error calculations).}

A simple test of incrementalRMA is to parameterizeRMA, then apply incrementalRMA and compare with
using rma on the set.

## Installation
You can install the IncrementalRMA package from GitHub with:
```
# install.packages("devtools")
devtools::install_github("steveneschrich/IncrementalRMA")
```

## Usage

There are two steps associated with `IncrementalRMA`: calculating the parameters (`parameterizeRMA`) and
applying the parameters (`incrementalRMA`). As mentioned above, it is possible to verify that the software
works by computing the parameters from an AffyBatch, applying the parameters to the same AffyBatch and comparing
this to the canonical RMA approach.

```
BiocManager::install("affydata")
data(Dilution)
all.equal(exprs(incrementalRMA(Dilution, parameterizeRMA(Dilution))), exprs(affy::rma(Dilution)))
```
