# bekkpub

An R package implementing Gaussian QMLE for **full and diagonal BEKK(1,1)**,
with robust or inverse-Hessian standard errors and publication tables. Uses
only base/recommended R packages. This is an initial research implementation,
not a CRAN release.

## Install directly from GitHub

Requires **R 4.1 or newer**. Run these commands in R or RStudio:

```r
# Run once to install the GitHub installer.
install.packages("remotes")

# Install bekkpub directly from the r/ folder in this GitHub repository.
remotes::install_github(
  "gvpnm/bekk-R-and-python",
  subdir = "r",
  build_vignettes = FALSE,
  upgrade = "never"
)

library(bekkpub)
packageVersion("bekkpub")
```

No manual ZIP download or extraction is needed. If bekkpub was already loaded
during an update, restart R before loading the updated package.

## Basic use

```r
library(bekkpub)
# Use the included historical percentage log returns.
data(stock_returns)
returns <- tail(stock_returns[, c("DAX", "FTSE")], 500)
fit <- bekk_fit(returns, type = "diagonal", se = "robust", starts = 3)
bekk_table(fit)                         # estimates, SE, z, p, 95% CI, stars
bekk_diagnostics(fit, lags = 10)         # raw and Holm-adjusted p-values
bekk_report(fit, "results.html")
bekk_report(fit, "results.tex", format = "latex")
write.csv(bekk_table(fit), "coefficients.csv", row.names = FALSE)
write.csv(bekk_diagnostics(fit), "diagnostics.csv", row.names = FALSE)
```

Change `type = "diagonal"` to `type = "full"` for full BEKK.
For your own data, supply a numeric matrix of aligned returns, with one column
per series. Do not transform the bundled returns again.

LaTeX exports are table fragments requiring `\\usepackage{booktabs}`. Wide tables
may need landscape orientation or `\\resizebox`. HTML reports are standalone.

## Model and inference

The convention is

`H[t] = C C' + A e[t-1] e[t-1]' A' + B H[t-1] B'`.

C is lower triangular with positive diagonal. Full BEKK uses unrestricted A
and B; diagonal BEKK uses diagonal A and B. The Gaussian conditional likelihood
includes its normalizing constant. The initial covariance is the residual
sample second moment. The optimizer restricts the spectral radius of
`A tensor A + B tensor B` to less than 0.9999, ensuring covariance stationarity
with a small numerical margin. Parameters are reported in their natural scale,
with a delta-method transformation for the log C diagonal.

Means are sample means by default. Inference on covariance parameters is
conditional on these estimated means; mean SEs are not provided. Robust SEs
use a numerical score sandwich, assuming the conditional covariance is correctly
specified and the usual QMLE regularity/moment conditions. Hessian SEs assume
Gaussian innovations. Failed convergence or a non-positive/ill-conditioned
Hessian suppresses inference. Check convergence, sensitivity to starting
values, and near-boundary estimates before publication. Full models grow
rapidly in dimension and this pure-R implementation suits small systems.

## Diagnostics

* Per-series Ljung–Box tests for standardized residuals and their squares.
* Per-series ARCH LM tests on squared standardized residuals.
* Per-series Jarque–Bera normality tests.
* Joint multivariate portmanteau tests for residuals and squared residuals.
* Spectral-radius stationarity measure, log likelihood, AIC, and BIC.

Diagnostic p-values are asymptotic screening values with nominal degrees of
freedom, without adjustment for estimated BEKK parameters. They are not
bootstrap-calibrated specification tests. Holm-adjusted p-values cover all
tests returned in a diagnostic table. Cholesky whitening depends on series
ordering. Normality rejection alone does not establish covariance misspecification.

## Reproducible example and verification

Run the installed example without downloading the repository:

```r
source(system.file("examples", "demo.R", package = "bekkpub"))
```

It simulates data, fits a model, and writes coefficient/diagnostic CSVs plus
HTML and LaTeX reports into `bekk-results` in the current working directory.
The example uses simulated data, not empirical financial data.

```sh
R CMD build r
R CMD check --no-manual bekkpub_0.1.2.tar.gz
```

Tests cover recursion and likelihood identities, scalar stationarity,
whitening, diagnostics, escaping, simulation, and a fitted-model workflow.

## References

Engle, R. F. and Kroner, K. F. (1995). Multivariate Simultaneous Generalized
ARCH. *Econometric Theory*, 11(1), 122–150.

Hosking, J. R. M. (1980). The multivariate portmanteau statistic.
*Journal of the American Statistical Association*, 75(371), 602–608.

## Included historical stock-index data and Windows example

`data(stock_returns, package = "bekkpub")` loads 1,859 daily percentage log
returns for DAX, SMI, CAC, and FTSE, derived from R's `EuStockMarkets` dataset
(1991–1998). These are equity indices, not individual stocks. No exact dates
are inferred. Original closes and returns are also bundled as CSV files under
`inst/extdata`; the transformation is in `data-raw/build-stock-data.R`.

Run the Windows example with:

```r
source(system.file("examples", "stock_example_windows.R", package = "bekkpub"))
```

It fits diagonal BEKK to the last 500 DAX/FTSE observations and saves reports,
CSV tables, the fitted object, and a volatility plot in
`C:/Users/Narayan/Documents/Codex/bekk`. Full BEKK and own-price CSV examples
are included as optional code. The folder is created when the script runs
on your Windows computer.

## Examples for every function

Each exported function has a self-contained Examples section in its R help
page. View it using `?bekk_fit`, `?bekk_table`, `?bekk_diagnostics`,
`?bekk_report`, or `?bekk_simulate`. Run it with `example(bekk_fit)`
(substitute another function name as needed). Each fitting example may take
some time. Optional full-model fits are marked as long-running examples.
Report examples write to a temporary folder by default; the help page also
shows how to save results in the requested Windows folder.

