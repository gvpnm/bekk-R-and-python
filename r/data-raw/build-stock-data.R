# Source: datasets::EuStockMarkets (distributed with R).
# Historical daily European stock-index closing values, 1991-1998.
prices <- matrix(as.numeric(datasets::EuStockMarkets),
                 ncol = ncol(datasets::EuStockMarkets),
                 dimnames = list(NULL, colnames(datasets::EuStockMarkets)))
stock_returns <- 100 * diff(log(prices))
attr(stock_returns, "source") <- "datasets::EuStockMarkets: daily European stock-index closes, 1991-1998"
attr(stock_returns, "units") <- "Percent log return; 100 * diff(log(index close))"
attr(stock_returns, "note") <- "Trading observations are ordered; no calendar dates are inferred from the ts time index."
save(stock_returns, file = "data/stock_returns.rda", compress = "xz")
utils::write.csv(data.frame(observation = seq_len(nrow(prices)), prices),
                 "inst/extdata/stock_index_prices.csv", row.names = FALSE)
utils::write.csv(data.frame(observation = 2:nrow(prices), stock_returns),
                 "inst/extdata/stock_returns.csv", row.names = FALSE)
