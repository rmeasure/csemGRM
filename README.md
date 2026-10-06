# csemGRM

**Conditional Standard Errors of Measurement for the Graded Response Model**

## Overview

`csemGRM` implements analytical and bootstrap estimators of the conditional standard error of measurement (CSEM) for the Graded Response Model (GRM; Samejima, 1969). While marginal reliability summarizes measurement precision with a single coefficient, the CSEM varies across the latent trait continuum — showing where a scale measures precisely and where it does not.

## Installation

```r
remotes::install_github("rmeasure/csemGRM")
```

## Quick Start

```r
library(csemGRM)

# Step 1: Fit GRM and compute information functions
data(attitude_data)
items <- attitude_data[, 1:10]
fit   <- csem_estimate(items)
print(fit)

# Step 2: Compute analytical CSEM
csem  <- csem_analytical(fit)
print(csem)

# Step 3: Bootstrap CSEM with confidence intervals
boot  <- csem_bootstrap(fit, n_boot = 100)
print(boot)

# Step 4: Conditional reliability
rel   <- csem_reliability(fit, csem)
print(rel)

# Step 5: Visualize
csem_plot(csem, type = "csem")
csem_plot(csem, type = "info")
csem_plot(rel,  type = "reliability")
csem_plot(fit,  type = "items")
csem_plot(boot, type = "bootstrap")
```

## Functions

| Function | Description |
|---|---|
| `csem_estimate()` | Fit GRM and extract information functions |
| `csem_analytical()` | Compute analytical CSEM = 1/sqrt(I(theta)) |
| `csem_bootstrap()` | Bootstrap CSEM with confidence intervals |
| `csem_reliability()` | Conditional reliability profile |
| `csem_plot()` | Visualize CSEM, TIF, reliability, item info |

## References

- Samejima, F. (1969). *Psychometrika Monograph Supplement, 17*.
- Samejima, F. (1994). *Applied Psychological Measurement, 18*(3), 229-244.
- Chalmers, R. P. (2012). *Journal of Statistical Software, 48*(6), 1-29.
