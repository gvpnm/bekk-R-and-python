\name{bekk_report}
\alias{bekk_report}
\title{Export publication tables}
\description{Export publication tables. See the package README for model conventions and inference limitations.}
\usage{
bekk_report(object, file = NULL, format = c("html", "latex"),
            lags = 10L, digits = 4L)
}
\arguments{
\item{object}{A fitted bekk_fit object.}
\item{file}{Output path, or NULL to return the report without writing.}
\item{format}{Standalone HTML or LaTeX table fragments requiring booktabs.}
\item{lags}{Diagnostic lag count.}
\item{digits}{Number of decimal places, from one to ten.}
}
\value{The report as a character string, invisibly. Writes the report if file is supplied.}
\examples{
# Each help example is self-contained.
data(stock_returns, package = "bekkpub")
returns <- tail(stock_returns[, c("DAX", "FTSE"), drop = FALSE], 250)
fit <- bekk_fit(returns, type = "diagonal", se = "robust", starts = 1)
# Self-contained HTML report in a temporary folder.
html_file <- file.path(tempdir(), "bekk_report.html")
bekk_report(fit, file = html_file, format = "html", lags = 10)
file.exists(html_file)
if (interactive()) browseURL(html_file)

# LaTeX table fragments: add booktabs to your document preamble.
tex_file <- file.path(tempdir(), "bekk_report.tex")
bekk_report(fit, file = tex_file, format = "latex", digits = 4)

# Return the HTML text without writing a file.
html_text <- bekk_report(fit, format = "html")
nchar(html_text)

# Save permanent results to your Windows folder when you choose.
\dontrun{
folder <- "C:/Users/Narayan/Documents/Codex/bekk"
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
bekk_report(fit, file.path(folder, "stock_report.html"))
bekk_report(fit, file.path(folder, "stock_report.tex"), format = "latex")
}
}
