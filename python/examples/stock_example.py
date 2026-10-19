"""Run: python stock_example.py --output "C:\\Users\\Narayan\\Documents\\Codex\\bekk"

No downloads are needed for the bundled stock-index sample data.
"""
from pathlib import Path
import argparse
import numpy as np
from bekkpub import load_stock_returns, bekk_fit, bekk_table, bekk_diagnostics, bekk_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path.home()/"Documents"/"Codex"/"bekk")
    parser.add_argument("--observations", type=int, default=500)
    parser.add_argument("--starts", type=int, default=3)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = load_stock_returns()
    returns = data[["DAX", "FTSE"]].tail(args.observations)
    fit = bekk_fit(returns, type="diagonal", starts=args.starts, tol=1e-4)
    print(fit)
    print(bekk_table(fit).to_string(index=False))
    print(bekk_diagnostics(fit).to_string(index=False))
    bekk_report(fit, args.output/"python_stock_report.html")
    bekk_report(fit, args.output/"python_stock_report.tex", format="latex")
    bekk_table(fit).to_csv(args.output/"python_coefficients.csv", index=False)
    bekk_diagnostics(fit).to_csv(args.output/"python_diagnostics.csv", index=False)
    data.to_csv(args.output/"stock_returns.csv")
    returns.to_csv(args.output/"python_fitted_returns.csv")
    np.savez_compressed(args.output/"python_fit_arrays.npz", C=fit.C, A=fit.A, B=fit.B,
                        H=fit.H, standardized_residuals=fit.standardized_residuals.to_numpy())
    print("Saved results to", args.output.resolve())


if __name__ == "__main__":
    main()
