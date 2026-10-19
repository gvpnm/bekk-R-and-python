\name{stock_returns}
\alias{stock_returns}
\docType{data}
\title{Historical European Stock-Index Returns}
\description{Daily percentage log returns computed from the European stock-index closing values in \code{datasets::EuStockMarkets}. These are observed equity-index data, not individual stock prices or simulated data.}
\usage{data(stock_returns)}
\format{A numeric matrix with 1859 observations and four columns: DAX (Germany), SMI (Switzerland), CAC (France), and FTSE (United Kingdom). Values are 100 times the first difference of log closing values. Rows retain trading-observation order; exact calendar dates are not supplied.}
\source{\code{datasets::EuStockMarkets}, distributed with R. The original dataset contains daily stock-index closes from 1991 to 1998. See \code{help("EuStockMarkets", package = "datasets")} for original provenance.}
\examples{
data(stock_returns)
head(stock_returns)
colnames(stock_returns)
}
