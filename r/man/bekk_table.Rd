\name{bekk_table}
\alias{bekk_table}
\title{Coefficient table}
\description{Coefficient table. See the package README for model conventions and inference limitations.}
\usage{
bekk_table(object, level = .95)
}
\arguments{
\item{object}{A fitted bekk_fit object.}
\item{level}{Confidence level between zero and one.}
}
\value{A data frame of estimates, standard errors, z statistics, p-values, confidence intervals, and significance stars.}
\examples{
# Each help example is self-contained.
data(stock_returns, package = "bekkpub")
returns <- tail(stock_returns[, c("DAX", "FTSE"), drop = FALSE], 250)
fit <- bekk_fit(returns, type = "diagonal", se = "robust", starts = 1)
# 95 percent confidence intervals (default).
coefficients <- bekk_table(fit)
print(coefficients)

# Request 99 percent confidence intervals.
bekk_table(fit, level = 0.99)

# Save a table without changing your working directory.
output_file <- file.path(tempdir(), "bekk_coefficients.csv")
write.csv(coefficients, output_file, row.names = FALSE)
}
