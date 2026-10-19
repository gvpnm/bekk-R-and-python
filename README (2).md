# BEKK packages for R and Python

Full and diagonal BEKK(1,1) multivariate GARCH estimation, publication tables,
residual diagnostics, and historical stock-index example data.

- `r/`: R package bekkpub 0.1.2, documentation, tests and examples.
- `python/`: Python package bekkpub 0.1.0, documentation, tests and examples.
- `packages/`: installable R source archive and Python wheel.

## Install directly from GitHub

You do not need to download or extract the repository manually. Choose the
instructions for your language below. Internet access is required during installation.

## R: install and run

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

Fit the included stock-index sample and export results:

```r
library(bekkpub)
data(stock_returns)
returns <- tail(stock_returns[, c("DAX", "FTSE")], 500)
fit <- bekk_fit(returns, type = "diagonal")
print(fit)
bekk_table(fit)
bekk_diagnostics(fit, lags = 10)
folder <- path.expand("~/Documents/bekk")
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

## Python: install and run

Requires **Python 3.10 or newer**. Install the wheel directly from GitHub.

In a terminal or Command Prompt:

```sh
python -m pip install "https://raw.githubusercontent.com/gvpnm/bekk-R-and-python/main/packages/bekkpub-0.1.0-py3-none-any.whl"
```

On Windows, you can use `py -m pip` instead of `python -m pip`.

In a Jupyter notebook or Google Colab cell:

```python
!pip install "https://raw.githubusercontent.com/gvpnm/bekk-R-and-python/main/packages/bekkpub-0.1.0-py3-none-any.whl"
```

In Jupyter, `%pip install` with the same URL is also supported and targets the
current notebook kernel. Restart the kernel after upgrading a package that
was already imported.

NumPy, SciPy and pandas are installed automatically if needed. No manual
repository download, extraction, or Git installation is required. This package
is hosted on GitHub, not PyPI; use the full URL above rather than
`pip install bekkpub`. To reinstall this wheel, add
`--force-reinstall --no-deps` to the command.

After installation, run the following in Python or a notebook cell:

```python
from pathlib import Path
from bekkpub import (load_stock_returns, bekk_fit, bekk_table,
                     bekk_diagnostics, bekk_report)
returns = load_stock_returns()[["DAX", "FTSE"]].tail(500)
fit = bekk_fit(returns, type="diagonal")
print(fit)
print(bekk_table(fit))
print(bekk_diagnostics(fit))
folder = Path.home() / "Documents" / "bekk"
bekk_report(fit, folder / "Python_report.html")
bekk_report(fit, folder / "Python_report.tex", format="latex")
help(bekk_fit)
```

Reports are written to `Documents/bekk` in your home directory. Change
`folder` in either example to choose a different output directory.

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

