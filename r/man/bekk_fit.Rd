\name{bekk_fit}
\alias{bekk_fit}
\title{Fit a BEKK(1,1) model}
\description{Fit a BEKK(1,1) model. See the package README for model conventions and inference limitations.}
\usage{
bekk_fit(data, type = c("full", "diagonal"), demean = TRUE,
         se = c("robust", "hessian", "none"), starts = 3L,
         control = list(maxit = 1500, reltol = 1e-8))
}
\arguments{
\item{data}{Finite numeric matrix of aligned return observations, at least two columns and 30 rows.}
\item{type}{Full or diagonal BEKK parameterization.}
\item{demean}{Subtract column sample means.}
\item{se}{Sandwich robust, inverse Hessian, or no standard errors.}
\item{starts}{Positive integer number of deterministic optimizer starting points.}
\item{control}{Control list passed to stats::optim.}
}
\value{A bekk_fit object containing coefficients, covariance estimates, covariance paths, standardized residuals, likelihood, and convergence information.}
\details{The recursion is H[t] = C C′ + A e[t-1] e[t-1]′ A′ + B H[t-1] B′. C is lower triangular with positive diagonal. The spectral radius of the covariance recursion is constrained below 0.9999. Mean estimation is separate and standard errors condition on the estimated means.}
\examples{
# Each help example is self-contained.
data(stock_returns, package = "bekkpub")
returns <- tail(stock_returns[, c("DAX", "FTSE"), drop = FALSE], 250)
fit <- bekk_fit(returns, type = "diagonal", se = "robust", starts = 1)
# Inspect convergence, persistence, and coefficient estimates.
print(fit)
fit$convergence
fit$radius
coef(fit)

# Full BEKK estimates cross-series ARCH/GARCH coefficients.
\donttest{
full_fit <- bekk_fit(returns, type = "full", starts = 3)
bekk_table(full_fit)
}
}
