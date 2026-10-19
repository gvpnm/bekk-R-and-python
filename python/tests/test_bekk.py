import tempfile
import unittest
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from bekkpub import (bekk_fit, bekk_simulate, bekk_table, bekk_diagnostics,
                     bekk_report, load_stock_returns)
from bekkpub.model import _filter, _layout, _radius, _hessian
from bekkpub.report import _holm, _ljung_box


class TestNumerics(unittest.TestCase):
    def test_filter_recursion_likelihood_whitening(self):
        rng = np.random.default_rng(3)
        # Exercise asymmetric full matrices and three variables.
        k = 3; ci, ai, _ = _layout(k, "full")
        C = np.array([[.4,0,0],[.1,.3,0],[.03,.02,.2]])
        A = np.diag([.2,.1,.15]); A[0,1] = .03
        B = np.diag([.8,.85,.7]); B[1,2] = .04
        cc = C.copy(); cc[np.diag_indices(k)] = np.log(C.diagonal())
        p = np.r_[cc[ci], A[tuple(ai)], B[tuple(ai)]]
        e = rng.normal(size=(90,k)); S=e.T@e/len(e)
        ll,H,z,*_ = _filter(p,e,S,ci,ai,True)
        expected = C@C.T + np.outer(A@e[0],A@e[0]) + B@S@B.T
        np.testing.assert_allclose(H[1],expected,atol=1e-12)
        for t in (0,1,40):
            direct= -.5*(k*np.log(2*np.pi)+np.linalg.slogdet(H[t])[1]+e[t]@np.linalg.solve(H[t],e[t]))
            self.assertAlmostEqual(ll[t],direct,places=11)
            np.testing.assert_allclose(np.linalg.cholesky(H[t])@z[t],e[t],atol=1e-12)
        self.assertAlmostEqual(_radius(np.eye(2)*.2,np.eye(2)*.8),.68)

    def test_hessian(self):
        Q=np.array([[4.,1.],[1.,3.]])
        actual=_hessian(lambda p: .5*p@Q@p, np.array([.2,-.4]))
        np.testing.assert_allclose(actual,Q,atol=1e-7)

    def test_simulation_and_validation(self):
        C=np.eye(3)*.3; A=np.eye(3)*.2; B=np.eye(3)*.85
        a=bekk_simulate(100,C,A,B,seed=4)
        np.testing.assert_array_equal(a,bekk_simulate(100,C,A,B,seed=4))
        self.assertEqual(a.shape,(100,3))
        with self.assertRaises(ValueError): bekk_simulate(10,C,A,np.eye(3))
        with self.assertRaises(ValueError): bekk_simulate(2.5,C,A,B)
        with self.assertRaises(ValueError): bekk_fit(np.ones((40,2)))
        with self.assertRaises(ValueError): bekk_fit(np.full((40,2),np.nan))

    def test_stock_data(self):
        x=load_stock_returns()
        self.assertEqual(x.shape,(1859,4))
        self.assertEqual(list(x.columns),["DAX","SMI","CAC","FTSE"])
        self.assertTrue(np.isfinite(x.to_numpy()).all())
        self.assertAlmostEqual(x.iloc[0,0],-.932655,places=5)

    def test_diagnostic_formulas(self):
        np.testing.assert_allclose(_holm(np.array([.01,.04,.03])),[.03,.06,.06])
        v=np.array([1.,2.,-1.,0.,3.,2.,-2.,1.])
        u=v-v.mean(); n=len(v)
        r=[np.dot(u[h:],u[:-h])/np.dot(u,u) for h in (1,2)]
        expected=n*(n+2)*sum(r[h-1]**2/(n-h) for h in (1,2))
        self.assertAlmostEqual(_ljung_box(v,2),expected)


class TestWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fit=bekk_fit(load_stock_returns()[["DAX","FTSE"]].tail(500),
                         type="diagonal",starts=1,tol=1e-4)

    def test_inference(self):
        f=self.fit
        self.assertTrue(f.converged,f.message)
        self.assertEqual(f.inference,"robust")
        self.assertTrue(np.isfinite(f.covariance.to_numpy()).all())
        self.assertTrue((np.diag(f.covariance)>0).all())
        self.assertLess(f.persistence,1)
        self.assertEqual(len(f.coefficients),7)
        self.assertEqual(f.npar,9)
        # Compare to the separately implemented R version's stock example.
        self.assertAlmostEqual(f.loglik,-1282.456,places=2)
        self.assertAlmostEqual(f.aic,-2*f.loglik+18)
        self.assertAlmostEqual(f.bic,-2*f.loglik+np.log(500)*9)

    def test_diagnostics_and_reports(self):
        f=self.fit; table=bekk_table(f); d=bekk_diagnostics(f)
        self.assertEqual(len(table),7); self.assertEqual(len(d),10)
        self.assertTrue((d.p_adjusted>=d.p_value).all())
        self.assertAlmostEqual(d.iloc[3].statistic,stats.jarque_bera(f.standardized_residuals.iloc[:,0]).statistic)
        with tempfile.TemporaryDirectory() as directory:
            target=Path(directory)/"nested"/"report.html"
            html=bekk_report(f,target)
            self.assertTrue(target.exists()); self.assertIn("<table",html)
            latex=bekk_report(f,format="latex")
            self.assertIn(r"\toprule",latex); self.assertIn(r"\_",latex)
            original=f.coefficients.index
            f.coefficients.index=["<script>&_",*original[1:]]
            try:
                self.assertNotIn("<script>",bekk_report(f))
                self.assertIn("&lt;script&gt;",bekk_report(f))
            finally: f.coefficients.index=original
        with self.assertRaises(ValueError): bekk_diagnostics(f,100)
        with self.assertRaises(ValueError): bekk_table(f,1)
        with self.assertRaises(ValueError): bekk_report(f,format="pdf")

    def test_full_fit(self):
        x=bekk_simulate(150,np.eye(2)*.3,np.eye(2)*.2,np.eye(2)*.8,seed=19)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            f=bekk_fit(x,type="full",starts=1,se="none",tol=1e-3,maxiter=300)
        self.assertEqual(len(f.coefficients),11)
        self.assertLess(f.persistence,1)
        self.assertTrue(np.isfinite(f.H).all())

    def test_failed_convergence_suppresses_inference(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            f=bekk_fit(load_stock_returns().iloc[-80:,:2],type="diagonal",maxiter=1,starts=1)
        self.assertFalse(f.converged)
        self.assertTrue(f.covariance.isna().all().all())


if __name__ == "__main__":
    unittest.main()
