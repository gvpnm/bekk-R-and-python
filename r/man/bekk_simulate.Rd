\name{bekk_simulate}
\alias{bekk_simulate}
\title{Simulate BEKK returns}
\description{Simulate BEKK returns. See the package README for model conventions and inference limitations.}
\usage{
bekk_simulate(n, C, A, B, burn = 500L, seed = NULL)
}
\arguments{
\item{n}{Positive integer sample size.}
\item{C}{Constant covariance factor, a square matrix.}
\item{A}{ARCH coefficient matrix.}
\item{B}{GARCH coefficient matrix.}
\item{burn}{Nonnegative integer burn-in length.}
\item{seed}{Optional integer random seed; preserves the calling random-number state.}
}
\value{A numeric matrix with n rows and one column per series.}
\examples{
# Two series with stationary BEKK dynamics.
C <- matrix(c(0.3, 0.06, 0, 0.25), nrow = 2)
A <- diag(c(0.22, 0.18))
B <- diag(c(0.88, 0.86))
x <- bekk_simulate(n = 500, C = C, A = A, B = B, seed = 42)
colnames(x) <- c("Stock1", "Stock2")
head(x)
dim(x)

# The same seed gives the same simulated returns.
identical(unname(x),
          bekk_simulate(n = 500, C = C, A = A, B = B, seed = 42))

\donttest{
fit <- bekk_fit(x, type = "diagonal")
bekk_table(fit)
bekk_diagnostics(fit)
}
}
