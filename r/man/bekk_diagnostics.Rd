\name{bekk_diagnostics}
\alias{bekk_diagnostics}
\title{Residual diagnostic tests}
\description{Residual diagnostic tests. See the package README for model conventions and inference limitations.}
\usage{
bekk_diagnostics(object, lags = 10L)
}
\arguments{
\item{object}{A fitted bekk_fit object.}
\item{lags}{Positive integer number of lags, smaller than one fifth of the sample size.}
}
\value{A data frame with test names, series, statistics, degrees of freedom, raw p-values and Holm-adjusted p-values. The notes attribute describes limitations.}
\details{Diagnostic p-values use nominal asymptotic degrees of freedom without adjustment for fitted BEKK parameters. Treat them as screening tests. Cholesky whitening depends on series order.}
\examples{
# Each help example is self-contained.
data(stock_returns, package = "bekkpub")
returns <- tail(stock_returns[, c("DAX", "FTSE"), drop = FALSE], 250)
fit <- bekk_fit(returns, type = "diagonal", se = "robust", starts = 1)
# Residual autocorrelation, remaining ARCH, and normality checks.
diagnostics <- bekk_diagnostics(fit, lags = 10)
print(diagnostics)

# Display tests significant after Holm adjustment.
subset(diagnostics, p_adjusted < 0.05)

# Read the interpretation limitations.
attr(diagnostics, "notes")

# Compare a shorter lag choice.
bekk_diagnostics(fit, lags = 5)
}
