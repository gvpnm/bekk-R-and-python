# BEKK packages for R and Python

Full and diagonal BEKK(1,1) multivariate GARCH estimation, publication tables,
residual diagnostics, and historical stock-index example data.

- `r/`: R package bekkpub 0.1.2, documentation, tests and examples.
- `python/`: Python package bekkpub 0.1.0, documentation, tests and examples.
- `packages/`: installable R source archive and Python wheel.

## Step 1: Download

On GitHub select **Code > Download ZIP**, then extract it under
`C:/Users/Narayan/Documents/Codex/bekk`. Identify the extracted repository
folder containing `r`, `python`, and `packages`.

## Step 2: Install and run R

Requires R >=4.1. Open R/RStudio and run:

```r
# Select packages/bekkpub_0.1.2.tar.gz from the extracted folder.
install.packages(file.choose(), repos = NULL, type = "source")
```

Restart R, then run:

```r
library(bekkpub)
data(stock_returns)
returns <- tail(stock_returns[, c("DAX", "FTSE")], 500)
fit <- bekk_fit(returns, type = "diagonal")
print(fit)
bekk_table(fit)
bekk_diagnostics(fit, lags = 10)
folder <- "C:/Users/Narayan/Documents/Codex/bekk"
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
bekk_report(fit, file.path(folder, "R_report.html"))
bekk_report(fit, file.path(folder, "R_report.tex"), format = "latex")
```

Examples are in every function's help page:

```r
?bekk_fit
example(bekk_fit)
?bekk_table
?bekk_diagnostics
?bekk_report
?bekk_simulate
```

Run the complete Windows example using `source()`, not `install.packages()`:

```r
source(system.file("examples", "stock_example_windows.R", package = "bekkpub"))
```

## Step 3: Install and run Python

Requires Python >=3.10. Open Command Prompt in the extracted repository folder:

```bat
py -m pip install .\packages\bekkpub-0.1.0-py3-none-any.whl
py .\python\examples\stock_example.py --output "C:\Users\Narayan\Documents\Codex\bekk"
```

Missing NumPy/SciPy/pandas dependencies are downloaded during installation.
Alternatively use the library in a Python script:

```python
from pathlib import Path
from bekkpub import (load_stock_returns, bekk_fit, bekk_table,
                     bekk_diagnostics, bekk_report)
returns = load_stock_returns()[["DAX", "FTSE"]].tail(500)
fit = bekk_fit(returns, type="diagonal")
print(fit)
print(bekk_table(fit))
print(bekk_diagnostics(fit))
folder = Path("C:/Users/Narayan/Documents/Codex/bekk")
bekk_report(fit, folder / "Python_report.html")
bekk_report(fit, folder / "Python_report.tex", format="latex")
help(bekk_fit)
```

## Models and interpretation

Change `type="diagonal"` to `type="full"` for full BEKK. At least two series
are required; joint systems with three or more series and separate pairwise
fits are supported. Diagonal refers to A and B; H remains a full covariance
matrix. See each language's README for joint and pairwise examples.

The convention is `H[t] = CC' + A e[t-1] e[t-1]' A' + B H[t-1] B'`.
Initial covariance is the residual second moment. Covariance stationarity is
enforced with a spectral-radius bound of 0.9999. Means are estimated separately,
and covariance inference conditions on these means. Check optimizer convergence
and standard-error availability before interpreting coefficients.

Residual tests include Ljung-Box on residuals and squares, ARCH LM,
Jarque-Bera and joint portmanteau tests. Diagnostic p-values are asymptotic
screening values without fitted-parameter corrections. Holm-adjusted p-values
are supplied. See package documentation for assumptions and limitations.

Bundled data are real historical European equity-index returns (DAX, SMI, CAC,
FTSE) derived from R's EuStockMarkets dataset, 1991-1998. They are not individual
stocks or current market data. No exact calendar dates are inferred.

These are initial research implementations, not CRAN or PyPI releases.
R and Python optimization results can differ. See `r/README.md` and
`python/README.md` for detailed validation commands and examples.

## License

See `r/LICENSE` with `r/DESCRIPTION` for the R package's MIT license declaration,
and `python/LICENSE` for the Python MIT license. Dataset provenance is described
in the package documentation.
