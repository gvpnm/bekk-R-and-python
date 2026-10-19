"""Examples for joint models and separate pairwise models (can take time)."""
from itertools import combinations
from bekkpub import load_stock_returns, bekk_fit, bekk_table, bekk_diagnostics

if __name__ == "__main__":
    returns = load_stock_returns()[["DAX", "FTSE", "CAC"]].tail(500)
    # A single three-variable covariance system.
    joint = bekk_fit(returns, type="diagonal", tol=1e-4)
    print(joint)
    print(bekk_table(joint))
    print(bekk_diagnostics(joint))
    # Separate bivariate systems; these estimates need not equal the joint fit.
    pairs = {}
    for pair in combinations(returns.columns, 2):
        fit = bekk_fit(returns[list(pair)], type="diagonal", tol=1e-4)
        pairs[pair] = fit
        print(pair, fit)
        print(bekk_table(fit))
