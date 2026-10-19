# Run after installing bekkpub 0.1.1. All files go to your requested folder.
library(bekkpub)
output_dir <- "C:/Users/Narayan/Documents/Codex/bekk"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Actual historical equity-index data included in the package.
data("stock_returns", package = "bekkpub")
cat("Available index series:", colnames(stock_returns), "\n")
print(head(stock_returns))

# Use the last 500 common trading observations for DAX and FTSE.
# Returns already are percentage log returns. Do not difference them again.
returns <- tail(stock_returns[, c("DAX", "FTSE"), drop = FALSE], 500)
fit <- bekk_fit(returns, type = "diagonal", se = "robust", starts = 3)
print(fit)
coefficients <- bekk_table(fit)
diagnostics <- bekk_diagnostics(fit, lags = 10)
print(coefficients)
print(diagnostics)

# Publication tables and machine-readable output.
bekk_report(fit, file.path(output_dir, "stock_report.html"), lags = 10)
bekk_report(fit, file.path(output_dir, "stock_report.tex"), format = "latex", lags = 10)
write.csv(coefficients, file.path(output_dir, "coefficients.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_dir, "diagnostics.csv"), row.names = FALSE)
write.csv(stock_returns, file.path(output_dir, "all_stock_returns.csv"), row.names = FALSE)
write.csv(returns, file.path(output_dir, "fitted_returns.csv"), row.names = FALSE)
saveRDS(fit, file.path(output_dir, "stock_fit.rds"))

# Plot estimated daily volatility in percentage-return units.
png(file.path(output_dir, "conditional_volatility.png"), width = 1200, height = 650)
volatility <- sapply(seq_len(nrow(returns)), function(t) sqrt(diag(fit$H[, , t])))
matplot(t(volatility), type = "l", lty = 1, col = c("navy", "firebrick"),
        xlab = "Trading observation in fitted sample", ylab = "Conditional volatility (%)")
legend("topright", legend = colnames(returns), col = c("navy", "firebrick"), lty = 1)
dev.off()
cat("Results saved in:", normalizePath(output_dir, winslash = "/"), "\n")
if (interactive()) browseURL(file.path(output_dir, "stock_report.html"))

# Optional full BEKK model; slower and estimates more parameters:
# full_fit <- bekk_fit(returns, type = "full", se = "robust", starts = 3)
# bekk_report(full_fit, file.path(output_dir, "full_stock_report.html"))

# To use your own CSV with columns Date, Stock1, Stock2 (adjusted close prices):
# prices <- read.csv(file.path(output_dir, "my_prices.csv"))
# prices <- prices[order(as.Date(prices$Date)), ]
# stopifnot(!anyDuplicated(prices$Date))
# p <- as.matrix(prices[, c("Stock1", "Stock2")])
# stopifnot(is.numeric(p), all(is.finite(p)), all(p > 0))
# own_returns <- 100 * diff(log(p))
# own_fit <- bekk_fit(own_returns, type = "diagonal")
# bekk_report(own_fit, file.path(output_dir, "my_stock_report.html"))
