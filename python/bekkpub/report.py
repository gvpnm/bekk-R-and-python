"""Coefficient tables, residual screening tests, and self-contained reports."""
from pathlib import Path
from importlib.resources import files
from html import escape
import numpy as np
import pandas as pd
from scipy import stats


NOTES = ("Diagnostic p-values use nominal asymptotic degrees of freedom without "
         "adjustment for estimated BEKK parameters; treat these as screening tests. "
         "Holm correction covers all tests in the diagnostic table. Cholesky whitening "
         "depends on series ordering. Normality rejection alone does not establish "
         "covariance misspecification.")


def load_stock_returns():
    """Load 1,859 historical percentage log returns for DAX, SMI, CAC and FTSE.

    Source: R datasets::EuStockMarkets (1991-1998). These are real equity-index
    observations, not individual stock prices. The index is an observation ID,
    not a calendar date. Returns are already 100 * diff(log(close)).

    Example
    -------
    >>> from bekkpub import load_stock_returns
    >>> returns = load_stock_returns()
    >>> returns.shape
    (1859, 4)
    >>> returns[["DAX", "FTSE"]].tail(500).shape
    (500, 2)
    """
    with files("bekkpub").joinpath("data/stock_returns.csv").open("rb") as f:
        return pd.read_csv(f, index_col="observation")


def bekk_table(result, level=.95):
    """Return coefficients, SEs, z, p-values, confidence intervals, and stars.

    Example
    -------
    >>> from bekkpub import load_stock_returns, bekk_fit, bekk_table
    >>> fit = bekk_fit(load_stock_returns().iloc[-250:, :2], type="diagonal")
    >>> table = bekk_table(fit, level=.95)
    >>> table.to_csv("coefficients.csv", index=False)
    """
    if not np.isscalar(level) or not np.isfinite(level) or not 0 < level < 1:
        raise ValueError("level must lie between zero and one")
    est = result.coefficients.to_numpy()
    se = np.sqrt(np.maximum(np.diag(result.covariance), 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = est / se
    p = 2 * stats.norm.sf(np.abs(z))
    q = stats.norm.ppf((1+level)/2)
    stars = np.select([p < .001, p < .01, p < .05], ["***", "**", "*"], default="")
    return pd.DataFrame(dict(term=result.coefficients.index, estimate=est, std_error=se,
        z=z, p_value=p, conf_low=est-q*se, conf_high=est+q*se, significance=stars))


def _ljung_box(v, lags):
    v = v - v.mean(); n = len(v)
    den = v @ v
    if den <= 0:
        return np.nan
    acf = np.array([v[h:] @ v[:-h] / den for h in range(1, lags+1)])
    return n*(n+2)*np.sum(acf**2/(n-np.arange(1, lags+1)))


def _port(v, lags):
    v = v - v.mean(axis=0); n, k = v.shape
    S = v.T @ v / n
    if np.linalg.matrix_rank(S) < k:
        return np.nan
    inv = np.linalg.inv(S)
    q = 0.
    for h in range(1, lags+1):
        ch = v[h:].T @ v[:-h] / n
        q += np.trace(ch @ inv @ ch.T @ inv)/(n-h)
    return n*n*q


def _holm(p):
    # Missing statistics remain missing and count toward the family size.
    out = np.full(len(p), np.nan); valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    out[order] = np.minimum(1, np.maximum.accumulate((len(p)-np.arange(len(order)))*p[order]))
    return out


def bekk_diagnostics(result, lags=10):
    """Residual Ljung-Box, squared Ljung-Box, ARCH LM, JB and joint tests.

    Returns a DataFrame with test, series, statistic, df, p_value, p_adjusted.
    See the 'notes' attribute for asymptotic-test limitations. No parameter-count
    adjustment is made to diagnostic degrees of freedom.

    Example
    -------
    >>> from bekkpub import load_stock_returns, bekk_fit, bekk_diagnostics
    >>> fit = bekk_fit(load_stock_returns().iloc[-250:, :2], type="diagonal")
    >>> diagnostics = bekk_diagnostics(fit, lags=10)
    >>> print(diagnostics)
    >>> diagnostics.to_csv("diagnostics.csv", index=False)
    """
    z = result.standardized_residuals.to_numpy(); n, k = z.shape
    if isinstance(lags, bool) or not isinstance(lags, (int, np.integer)) or not 1 <= lags < n/5:
        raise ValueError("lags must be a positive integer smaller than n/5")
    rows = []
    def add(test, series, stat, df):
        rows.append(dict(test=test, series=series, statistic=float(stat), df=df,
                         p_value=float(stats.chi2.sf(stat, df))))
    for i, label in enumerate(result.residuals.columns):
        v = z[:, i]
        add("Ljung-Box residuals", label, _ljung_box(v, lags), lags)
        add("Ljung-Box squared residuals", label, _ljung_box(v*v, lags), lags)
        squares = v*v
        y = squares[lags:]
        X = np.column_stack([np.ones(len(y))]+[squares[lags-h:n-h] for h in range(1, lags+1)])
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        sst = np.sum((y-y.mean())**2)
        r2 = max(0, 1-np.sum((y-X@beta)**2)/sst) if sst > 0 else np.nan
        lm = len(y)*r2 if np.linalg.matrix_rank(X) == lags+1 else np.nan
        add("ARCH LM", label, lm, lags)
        add("Jarque-Bera", label, stats.jarque_bera(v).statistic, 2)
    add("Multivariate portmanteau", "Joint", _port(z, lags), k*k*lags)
    add("Multivariate portmanteau squared", "Joint", _port(z*z, lags), k*k*lags)
    d = pd.DataFrame(rows)
    d["p_adjusted"] = _holm(d.p_value.to_numpy())
    d.attrs["notes"] = NOTES
    return d


def bekk_report(result, file=None, format="html", lags=10, digits=4):
    """Return a publication report and optionally write it as HTML or LaTeX.

    HTML is standalone. LaTeX returns table fragments requiring booktabs and
    may need landscape orientation. Parent output directories are created.

    Example
    -------
    >>> from pathlib import Path
    >>> from bekkpub import load_stock_returns, bekk_fit, bekk_report
    >>> fit = bekk_fit(load_stock_returns().iloc[-250:, :2], type="diagonal")
    >>> folder = Path.home() / "Documents" / "Codex" / "bekk"
    >>> html = bekk_report(fit, folder / "python_report.html")
    >>> latex = bekk_report(fit, folder / "python_report.tex", format="latex")
    """
    if format not in ("html", "latex"):
        raise ValueError("format must be html or latex")
    if isinstance(digits, bool) or not isinstance(digits, (int, np.integer)) or not 1 <= digits <= 10:
        raise ValueError("digits must be an integer from 1 to 10")
    coeff = bekk_table(result); diag = bekk_diagnostics(result, lags)
    info = pd.DataFrame(dict(metric=["Observations", "Parameters including means", "Log likelihood",
        "AIC", "BIC", "Covariance persistence", "Optimizer converged"],
        value=[result.nobs,result.npar,result.loglik,result.aic,result.bic,result.persistence,int(result.converged)]))
    notes = (f"BEKK(1,1), {result.type}; SE: {result.inference}. "
        "H[t] = C C' + A e[t-1] e[t-1]' A' + B H[t-1] B'. "
        "C is lower triangular with positive diagonal. Stationarity requires the spectral "
        "radius of A tensor A + B tensor B < 1. Initial H is the residual second moment. "
        "Means are estimated separately and covariance SEs condition on these means. "
        "Intervals are 95% asymptotic Wald intervals. Stars: * p<0.05, ** p<0.01, *** p<0.001. "
        "Rounded zero p-values are not exactly zero. " + NOTES)
    tables = [("Parameter estimates",coeff),("Model fit and stationarity",info),("Residual diagnostics",diag)]
    def fmt(v):
        return f"{v:.{digits}f}"
    if format == "html":
        body = "\n".join(f"<h2>{title}</h2>"+d.to_html(index=False, escape=True, float_format=fmt, na_rep="NA") for title,d in tables)
        output = ('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>BEKK report</title>'
            '<style>body{font-family:Georgia,serif;margin:40px;max-width:1400px}table{border-collapse:collapse;'
            'font-variant-numeric:tabular-nums}th,td{padding:8px;text-align:right}thead{border-top:2px solid;'
            'border-bottom:1px solid}tbody{border-bottom:2px solid}p{line-height:1.5}</style></head><body>'
            '<h1>BEKK estimation report</h1>'+body+'<p>'+escape(notes)+'</p></body></html>')
    else:
        # Manual escaping avoids optional pandas/Jinja dependencies and injection.
        def tex(s):
            replacements = {"\\":r"\textbackslash{}", "&":r"\&", "%":r"\%", "$":r"\$",
                "#":r"\#", "_":r"\_", "{":r"\{", "}":r"\}", "~":r"\textasciitilde{}", "^":r"\textasciicircum{}"}
            return "".join(replacements.get(c,c) for c in str(s))
        pieces = [r"% Requires \usepackage{booktabs}. Wide tables may need landscape."]
        for title, d in tables:
            lines = [r"\begin{table}[htbp]",r"\centering",r"\caption{"+tex(title)+"}",
                r"\begin{tabular}{l"+"r"*(len(d.columns)-1)+"}",r"\toprule",
                " & ".join(tex(c) for c in d.columns)+r" \\",r"\midrule"]
            for row in d.itertuples(index=False, name=None):
                values = ["NA" if pd.isna(v) else fmt(v) if isinstance(v,(float,np.floating)) else str(v) for v in row]
                lines.append(" & ".join(tex(v) for v in values)+r" \\")
            pieces.append("\n".join(lines+[r"\bottomrule",r"\end{tabular}",r"\end{table}"]))
        output = "\n\n".join(pieces)+"\n\\par "+tex(notes)
    if file is not None:
        path = Path(file); path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(output, encoding="utf-8")
    return output
