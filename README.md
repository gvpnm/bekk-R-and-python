# bekkpub for Python

A Python implementation of full and diagonal BEKK(1,1) multivariate GARCH,
with Gaussian QMLE, score-sandwich or inverse-Hessian standard errors,
publication tables, and residual diagnostics. This is an initial research
implementation; it is not published on PyPI. Python >=3.10 is required.

## Install directly from GitHub

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

Verify the installation:

```python
import bekkpub
print(bekkpub.__version__)
```

## Basic use

```python
from pathlib import Path
from bekkpub import (
    load_stock_returns, bekk_fit, bekk_table, bekk_diagnostics, bekk_report,
)

folder = Path.home() / "Documents" / "bekk"
folder.mkdir(parents=True, exist_ok=True)
returns = load_stock_returns()[["DAX", "FTSE"]].tail(500)
fit = bekk_fit(returns, type="diagonal", se="robust")
print(fit)
print(bekk_table(fit))
print(bekk_diagnostics(fit, lags=10))
bekk_report(fit, folder / "report.html")
bekk_report(fit, folder / "report.tex", format="latex")
bekk_table(fit).to_csv(folder / "coefficients.csv", index=False)
bekk_diagnostics(fit).to_csv(folder / "diagnostics.csv", index=False)
```

HTML is standalone; LaTeX contains table fragments requiring `booktabs`.
Wide tables may need landscape or resizing. Every public function has examples:

```python
help(bekk_fit)
help(bekk_table)
help(bekk_diagnostics)
help(bekk_report)
help(load_stock_returns)
```

Full BEKK: change `type="diagonal"` to `type="full"`. Full models are slower.
There is no fixed series-count limit: at least two columns are needed. Large
systems require more data and computing time. Start with 2–4 series.

## Joint and pairwise estimation

```python
from itertools import combinations
returns = load_stock_returns()[["DAX", "FTSE", "CAC"]].tail(500)
joint = bekk_fit(returns, type="diagonal")
pairs = {
    pair: bekk_fit(returns[list(pair)], type="diagonal")
    for pair in combinations(returns.columns, 2)
}
print(bekk_table(pairs[("DAX", "FTSE")]))
```

Diagonal refers to A and B, not H. A joint model estimates a full conditional
covariance matrix; separate pairwise models need not agree with that joint fit.

## Your own stock data

Use aligned adjusted closing prices. Align dates before computing returns,
sort by date, and resolve missing observations and duplicate dates explicitly.
Do not pass dates as model columns. If data are already returns, do not
transform them again.

```python
import numpy as np
import pandas as pd
prices = pd.read_csv(folder / "prices.csv", parse_dates=["Date"])
prices = prices.sort_values("Date").set_index("Date")
assert not prices.index.duplicated().any()
p = prices[["Stock1", "Stock2", "Stock3"]]
assert np.isfinite(p.to_numpy()).all() and (p > 0).all().all()
returns = 100 * np.log(p).diff().iloc[1:]
fit = bekk_fit(returns, type="diagonal")
```

## Model and limitations

`H[t] = C C' + A e[t-1] e[t-1]' A' + B H[t-1] B'`.
C is lower triangular with positive diagonal. The initial covariance is the
residual sample second moment. The spectral radius of
`A tensor A + B tensor B` is restricted below 0.9999. Standard errors use
finite-difference scores/Hessian and a delta transformation for C diagonal.
Means are estimated separately; covariance inference conditions on these means.
AIC/BIC include estimated means in the parameter count. Wald intervals may
extend outside parameter constraints.

Check `fit.converged`, `fit.message`, and `fit.persistence` before inference.
Unsuccessful optimization or a non-positive/ill-conditioned Hessian suppresses
inference. Robust SEs assume correct conditional covariance and standard QMLE
regularity/moment conditions; Gaussian Hessian SEs additionally assume Gaussian
innovations. Optimizer results can differ from the R implementation.
Multiple deterministic starts do not guarantee a global optimum. This pure
Python implementation is intended for small research systems.

Diagnostic tests: per-series Ljung–Box on residuals and squares, ARCH LM,
Jarque–Bera, and joint multivariate portmanteau tests on residuals and squares.
P-values use nominal asymptotic degrees of freedom without fitted-parameter
corrections and should be interpreted as screening results. Holm adjustment
covers all tests in a diagnostic table. Whitening is Cholesky-based and depends
on series order. Normality rejection alone does not establish covariance
misspecification. Bootstrap-calibrated tests are not included.

## Sample data

1,859 historical daily percentage log returns for DAX, SMI, CAC and FTSE,
derived as `100 * diff(log(close))` from R's `datasets::EuStockMarkets`
(1991–1998). These are real stock-index data, not individual stocks or
simulated returns. Index values identify trading observations, not exact
calendar dates. See the R dataset documentation for original provenance.

## Tests

From the source directory: `python -m unittest discover -s tests -v`.
No empirical conclusion should be inferred from a passing software test.

References: Engle & Kroner (1995), Econometric Theory 11, 122–150;
Hosking (1980), JASA 75, 602–608.

