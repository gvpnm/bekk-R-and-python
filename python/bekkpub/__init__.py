"""BEKK covariance models. See each public function's docstring for examples."""
from .model import bekk_fit, bekk_simulate, BEKKResult
from .report import bekk_table, bekk_diagnostics, bekk_report, load_stock_returns

__version__ = "0.1.0"
__all__ = ["bekk_fit", "bekk_simulate", "BEKKResult", "bekk_table",
           "bekk_diagnostics", "bekk_report", "load_stock_returns"]
