"""Gaussian QMLE with the convention H = CC' + Aee'A' + BHB'."""
from dataclasses import dataclass
import warnings
import numpy as np
import pandas as pd
from scipy import optimize, linalg


def _radius(A, B):
    return float(np.max(np.abs(np.linalg.eigvals(np.kron(A, A) + np.kron(B, B)))))


def _layout(k, kind):
    ci = np.tril_indices(k)
    ai = np.diag_indices(k) if kind == "diagonal" else np.indices((k, k)).reshape(2, -1)
    names = [f"{tag}[{i+1},{j+1}]" for tag, ix in (("C", ci), ("A", ai), ("B", ai)) for i, j in zip(*ix)]
    return ci, ai, names


def _unpack(p, k, ci, ai):
    nc, na = len(ci[0]), len(ai[0])
    C, A, B = (np.zeros((k, k)) for _ in range(3))
    C[ci] = p[:nc]
    with np.errstate(over="raise", invalid="raise"):
        C[np.diag_indices(k)] = np.exp(C.diagonal())
    A[tuple(ai)] = p[nc:nc+na]
    B[tuple(ai)] = p[nc+na:]
    return C, A, B


def _filter(p, e, S, ci, ai, keep=False):
    try:
        C, A, B = _unpack(p, e.shape[1], ci, ai)
        rho = _radius(A, B)
        if not np.isfinite(rho) or rho >= .9999:
            return None
        H = S.copy()
        Q = C @ C.T
        ll = np.empty(len(e))
        hs = np.empty((len(e), *H.shape)) if keep else None
        zs = np.empty_like(e) if keep else None
        for t, row in enumerate(e):
            if t:
                ar = A @ e[t-1]
                H = Q + np.outer(ar, ar) + B @ H @ B.T
            H = (H + H.T) / 2
            L = np.linalg.cholesky(H)
            z = linalg.solve_triangular(L, row, lower=True, check_finite=False)
            ll[t] = -.5 * (len(row)*np.log(2*np.pi) + 2*np.log(L.diagonal()).sum() + z@z)
            if keep:
                hs[t], zs[t] = H, z
        if not np.isfinite(ll).all():
            return None
        return (ll, hs, zs, C, A, B, rho) if keep else ll
    except (np.linalg.LinAlgError, FloatingPointError, ValueError, OverflowError):
        return None


def _hessian(fun, p):
    # Symmetric central differences; deliberately do not use BFGS's approximation.
    step = np.finfo(float).eps**.25 * (1 + np.abs(p))
    h = np.zeros((len(p), len(p)))
    f0 = fun(p)
    for i in range(len(p)):
        ei = np.zeros_like(p); ei[i] = step[i]
        h[i, i] = (fun(p+ei) - 2*f0 + fun(p-ei)) / step[i]**2
        for j in range(i):
            ej = np.zeros_like(p); ej[j] = step[j]
            h[i, j] = h[j, i] = (fun(p+ei+ej)-fun(p+ei-ej)-fun(p-ei+ej)+fun(p-ei-ej))/(4*step[i]*step[j])
    return h


@dataclass
class BEKKResult:
    """Fitted model. H has shape (observations, series, series).

    Example
    -------
    >>> from bekkpub import load_stock_returns, bekk_fit
    >>> fit = bekk_fit(load_stock_returns().iloc[-250:, :2], type="diagonal")
    >>> fit.converged, fit.persistence
    """
    coefficients: pd.Series
    covariance: pd.DataFrame
    C: np.ndarray
    A: np.ndarray
    B: np.ndarray
    H: np.ndarray
    residuals: pd.DataFrame
    standardized_residuals: pd.DataFrame
    mean: pd.Series
    loglik: float
    npar: int
    converged: bool
    message: str
    persistence: float
    inference: str
    type: str

    @property
    def nobs(self):
        return len(self.residuals)

    @property
    def aic(self):
        return -2*self.loglik + 2*self.npar

    @property
    def bic(self):
        return -2*self.loglik + np.log(self.nobs)*self.npar

    def __repr__(self):
        return (f"BEKK(1,1) {self.type}: n={self.nobs}, logLik={self.loglik:.3f}, "
                f"converged={self.converged}, persistence={self.persistence:.5f}, SE={self.inference}")


def bekk_fit(data, type="full", demean=True, se="robust", starts=3, maxiter=1500, tol=1e-5):
    """Fit stationary full or diagonal BEKK(1,1) using Gaussian QMLE.

    Parameters
    ----------
    data : pandas.DataFrame or numeric array, shape (n, k)
        Aligned return observations; n >= 30 and k >= 2. No missing values.
    type : {'full', 'diagonal'}
        Structure of A and B. Conditional covariance remains a full matrix.
    demean : bool
        Estimate means separately. Covariance inference conditions on these means.
    se : {'robust', 'hessian', 'none'}
        Score sandwich, inverse observed Hessian, or omit inference.
    starts : int
        Number of deterministic starting points.
    maxiter : int
        Maximum BFGS iterations per starting point.
    tol : float
        BFGS gradient tolerance; inspect convergence before interpreting inference.

    Returns
    -------
    BEKKResult
        Estimates, covariance path, residuals, and fit statistics. H is initialized
        at the residual second moment. Persistence is constrained below 0.9999.
        Failed convergence or an ill-conditioned Hessian suppresses inference.

    Example
    -------
    >>> from bekkpub import load_stock_returns, bekk_fit, bekk_table
    >>> x = load_stock_returns().loc[:, ["DAX", "FTSE"]].tail(500)
    >>> fit = bekk_fit(x, type="diagonal", starts=3)
    >>> print(fit)
    >>> print(bekk_table(fit))
    """
    if type not in ("full", "diagonal") or se not in ("robust", "hessian", "none"):
        raise ValueError("invalid type or se")
    if isinstance(starts, bool) or not isinstance(starts, (int, np.integer)) or starts < 1:
        raise ValueError("starts must be a positive integer")
    if not isinstance(maxiter, (int, np.integer)) or maxiter < 1 or not np.isfinite(tol) or tol <= 0:
        raise ValueError("maxiter and tol must be positive")
    x = np.asarray(data, dtype=float)
    if x.ndim != 2 or x.shape[1] < 2 or len(x) < 30 or not np.isfinite(x).all():
        raise ValueError("data must be finite, numeric, with >=30 rows and >=2 columns")
    n, k = x.shape
    columns = data.columns.copy() if isinstance(data, pd.DataFrame) else pd.Index([f"Series {j+1}" for j in range(k)])
    index = data.index.copy() if isinstance(data, pd.DataFrame) else pd.RangeIndex(n)
    mu = x.mean(axis=0) if demean else np.zeros(k)
    e = x - mu
    S = e.T @ e / n
    if np.linalg.eigvalsh(S).min() <= np.linalg.eigvalsh(S).max()*1e-12:
        raise ValueError("data covariance must be positive definite and well conditioned")
    ci, ai, names = _layout(k, type)
    if n <= len(names)+1:
        raise ValueError("too few observations for the parameter count")

    def objective(p):
        ll = _filter(p, e, S, ci, ai)
        return 1e20 if ll is None else -ll.sum()

    candidates = []
    for j in range(starts):
        a, b = .12+.08*(j % 4), .85-.05*(j % 4)
        C = np.linalg.cholesky(S*(1-a*a-b*b))
        C[np.diag_indices(k)] = np.log(C.diagonal())
        p = np.r_[C[ci], (np.eye(k)*a)[tuple(ai)], (np.eye(k)*b)[tuple(ai)]]
        candidates.append(optimize.minimize(objective, p, method="BFGS", jac="3-point", options={"maxiter": maxiter, "gtol": tol}))
    best = min(candidates, key=lambda r: r.fun)
    # Prefer an equally good converged start when the best run ended in precision loss.
    tied = [r for r in candidates if r.success and r.fun <= best.fun + 1e-6]
    if tied:
        best = min(tied, key=lambda r: r.fun)
    f = _filter(best.x, e, S, ci, ai, True)
    if f is None:
        raise RuntimeError("no stationary solution was found")
    ll, H, z, C, A, B, rho = f
    converged = bool(best.success)
    if not converged:
        warnings.warn(f"Optimizer did not converge: {best.message}. Inference is suppressed.", RuntimeWarning)
    V = np.full((len(names), len(names)), np.nan)
    inference = "not requested"
    if se != "none":
        inference = "unavailable: optimizer did not converge"
        if converged:
            h = _hessian(objective, best.x)
            ev = np.linalg.eigvalsh(h)
            inference = "unavailable: non-positive or ill-conditioned Hessian"
            if np.isfinite(h).all() and ev.min() > max(ev.max()*1e-10, 0):
                bread = np.linalg.inv(h)
                V = bread
                if se == "robust":
                    scores = np.empty((n, len(names)))
                    for i in range(len(names)):
                        step = 1e-5*(1+abs(best.x[i]))
                        delta = np.zeros(len(names)); delta[i] = step
                        lp = _filter(best.x+delta, e, S, ci, ai)
                        lm = _filter(best.x-delta, e, S, ci, ai)
                        scores[:, i] = np.nan if lp is None or lm is None else (lp-lm)/(2*step)
                    V = bread @ (scores.T @ scores) @ bread
                d = np.ones(len(names))
                ids = np.flatnonzero(ci[0] == ci[1]); d[ids] = C.diagonal()
                V *= np.outer(d, d)
                inference = se if np.isfinite(V).all() else "unavailable: boundary scores"
        if not np.isfinite(V).all():
            warnings.warn(inference, RuntimeWarning)
    estimates = np.r_[C[ci], A[tuple(ai)], B[tuple(ai)]]
    return BEKKResult(pd.Series(estimates, index=names, name="estimate"),
        pd.DataFrame(V, index=names, columns=names), C, A, B, H,
        pd.DataFrame(e, index=index, columns=columns), pd.DataFrame(z, index=index, columns=columns),
        pd.Series(mu, index=columns), float(ll.sum()), len(names)+(k if demean else 0),
        converged, str(best.message), rho, inference, type)


def bekk_simulate(n, C, A, B, burn=500, seed=None):
    """Simulate Gaussian BEKK returns without changing NumPy's global RNG.

    n and burn are integer counts. C, A, B must be finite square matrices of
    identical dimension >=2, with covariance persistence <1 and positive-definite
    unconditional covariance. Returns an array of shape (n, k).

    Example
    -------
    >>> import numpy as np
    >>> from bekkpub import bekk_simulate
    >>> x = bekk_simulate(500, np.diag([.3, .25]), np.eye(2)*.2,
    ...                   np.eye(2)*.85, seed=42)
    >>> x.shape
    (500, 2)
    """
    for count, minimum in ((n, 1), (burn, 0)):
        if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < minimum:
            raise ValueError("n and burn must be integer counts")
    C, A, B = (np.asarray(m, dtype=float) for m in (C, A, B))
    if C.ndim != 2 or C.shape[0] < 2 or C.shape[0] != C.shape[1] or any(m.shape != C.shape or not np.isfinite(m).all() for m in (C, A, B)):
        raise ValueError("C, A, B must be finite square matrices of identical size >=2")
    if _radius(A, B) >= 1:
        raise ValueError("model is not covariance stationary")
    k = len(C); Q = C @ C.T
    H = np.linalg.solve(np.eye(k*k)-np.kron(A,A)-np.kron(B,B), Q.ravel(order="F")).reshape((k,k), order="F")
    np.linalg.cholesky(H)
    rng = np.random.default_rng(seed)
    out = np.empty((n+burn, k))
    for t in range(n+burn):
        if t:
            ar = A @ out[t-1]
            H = Q + np.outer(ar,ar) + B@H@B.T
        out[t] = np.linalg.cholesky((H+H.T)/2) @ rng.standard_normal(k)
    return out[burn:]
