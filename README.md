# QAL

**Quantile adaptive lasso for quantile treatment effect estimation**

R implementation of the method in Liu Y, Wei K, Huang C, Yu Y, Qin G.
[Quantile adaptive lasso: variable selection for quantile treatment effect estimation](https://doi.org/10.1080/00949655.2024.2424348).
*Journal of Statistical Computation and Simulation*, 2025;95(2):239–257.
Published online November 8, 2024.

QAL selects covariates separately for each outcome quantile, using quantile
regression coefficients as adaptive penalty weights. It tunes selection by
weighted absolute mean difference (wAMD), then estimates a marginal quantile
treatment effect using inverse probability weighted quantile regression.
**This implementation supports binary treatment coded 0/1.**

This is newly packaged research software, not a claim of reproducing all
published simulation tables or the CHARLS application. It is not on CRAN.

## Installation

R >= 4.1. Install the archived LQA dependency explicitly:

```r
install.packages(c("quantreg", "remotes"))
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/lqa/lqa_1.0-3.tar.gz",
  repos = NULL, type = "source"
)
# Available after the repository is published:
remotes::install_github("lyh199800/QAL")
```

For a downloaded source archive, after installing dependencies:

```r
install.packages("QAL_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Example: 25th and 75th quantiles

```r
library(QAL)
set.seed(2026)
d <- simulate_qal(n = 500, p = 10)
x <- paste0("X", 1:10)
f25 <- qal(d, x, "A", "Y", tau = .25)
f75 <- qal(d, x, "A", "Y", tau = .75)
f25
f75
f25$selected
f25$estimate
f25$quantiles
f25$balance
f25$tuning
predict(f25, c(0, 1))
```

The selected covariates may differ by tau; differences are not guaranteed for
every sample. `tau` is explicit and independent of global workspace variables.
To force known covariates into every adjustment set:

```r
f <- qal(d, x[-1], "A", "Y", tau = .25, known = "X1")
# Sensitivity version corresponding to the gamma5 source:
f5 <- qal(d, x, "A", "Y", tau = .25, gamma_convergence = 5,
          lambda_exponents = c(-5, -2, -1, -.75, -.5, -.25, .25, .49))
```

## Interpretation and diagnostics

- The effect is the difference between treated and control marginal outcome
  quantiles, not the quantile of individual treatment effects.
- Provide complete finite numeric data; encode categories as indicator columns.
- The unpenalized full quantile model requires n > covariate count + 2 and a
  full-rank design; document any screening performed before fitting.
- Inspect balance, overlap, effective sample size and tuning failures.
  Propensity clipping is optional (`ps_clip`) and disabled by default.
- No SE or confidence interval is automatically reported. Generic `rq`
  summaries ignore selection and propensity estimation uncertainty. Any
  inferential procedure must account for the complete pipeline.
- Causal interpretation requires consistency, no unmeasured confounding and positivity.

## Citation and provenance

Run `citation("QAL")`. See [DEVELOPMENT.md](DEVELOPMENT.md) for source mapping,
robustness changes and validation scope. Only method code and synthetic examples
are included; original manuscripts, CHARLS data and comparison-method code
are not redistributed.

Maintainer: Yahang Liu <yhliu22@m.fudan.edu.cn>.
License: GPL-3. Version: 0.1.0.
