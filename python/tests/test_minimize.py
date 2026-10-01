"""Tests of the Python bindings (``sqpopt.minimize``); skipped if the extension isn't built.

Build it first:  pixi run build-python
Run from the repository root:  pixi run python -m unittest discover -s python/tests
"""

import pathlib
import sys
import unittest

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

try:
    from sqpopt import Bounds, LinearConstraint, NonlinearConstraint, minimize
    from sqpopt._minimize import _native
    _native()
    BUILT = True
except ImportError:
    BUILT = False

# HS71: min x1 x4 (x1+x2+x3) + x3  s.t.  x1 x2 x3 x4 >= 25,  sum x^2 = 40,  1 <= x <= 5
X_HS71 = np.array([1.0, 4.74299963, 3.82114998, 1.37940829])
F_HS71 = 17.0140173


def hs71_f(x):
    return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]


def hs71_g(x):
    return np.array([x[3] * (2 * x[0] + x[1] + x[2]), x[0] * x[3], x[0] * x[3] + 1, x[0] * (x[0] + x[1] + x[2])])


def hs71_c(x):
    return np.array([np.sum(x ** 2), np.prod(x)])


def hs71_cj(x):
    return np.array([2 * x, [x[1] * x[2] * x[3], x[0] * x[2] * x[3], x[0] * x[1] * x[3], x[0] * x[1] * x[2]]])


def rosen(x):
    return (1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2


def rosen_g(x):
    return np.array([-2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0] ** 2), 200 * (x[1] - x[0] ** 2)])


def rosen_h(x):
    return np.array([[2 - 400 * (x[1] - 3 * x[0] ** 2), -400 * x[0]], [-400 * x[0], 200.0]])


class Coo:
    """a minimal sparse matrix (anything with ``tocoo()`` is accepted, e.g. ``scipy.sparse``)"""

    def __init__(self, row, col, data, shape):
        self.row, self.col, self.data, self.shape = np.array(row), np.array(col), np.array(data, float), shape

    def tocoo(self):
        return self


@unittest.skipUnless(BUILT, 'the sqpopt extension is not built (pixi run build-python)')
class TestMinimize(unittest.TestCase):

    def check_hs71(self, r, tol=1e-6):
        self.assertTrue(r.success, r.message)
        self.assertAlmostEqual(r.fun, F_HS71, places=5)
        np.testing.assert_allclose(r.x, X_HS71, atol=tol)

    def test_hs71_analytic(self):
        r = minimize(hs71_f, [1, 5, 5, 1], jac=hs71_g, bounds=[(1, 5)] * 4,
                     constraints=[NonlinearConstraint(hs71_c, [40, 25], [40, np.inf], jac=hs71_cj)])
        self.check_hs71(r)
        self.assertEqual(r.status, 0)
        self.assertEqual(len(r.v), 1)
        self.assertGreater(r.v[0][1], 0)                 # at its lower bound
        self.assertGreater(r.z[0], 0)                    # x1 at its lower bound
        self.assertEqual(r.derivative_switch_iteration, 0)

    def test_derivatives_required(self):
        con = NonlinearConstraint(hs71_c, [40, 25], [40, np.inf], jac=hs71_cj)
        for jac in (None, '2-point', '3-point'):
            with self.assertRaises(ValueError):
                minimize(hs71_f, [1, 5, 5, 1], jac=jac, constraints=con)
        for bad in (NonlinearConstraint(hs71_c, [40, 25], [40, np.inf]),
                    {'type': 'eq', 'fun': lambda x: np.sum(x ** 2) - 40}):
            with self.assertRaises(ValueError):
                minimize(hs71_f, [1, 5, 5, 1], jac=hs71_g, constraints=bad)

    def test_slsqp_style_dicts(self):
        cons = [{'type': 'eq', 'fun': lambda x: np.sum(x ** 2) - 40, 'jac': lambda x: 2 * x},
                {'type': 'ineq', 'fun': lambda x, k: np.prod(x) - k, 'jac': lambda x, k: hs71_cj(x)[1],
                 'args': (25,)}]
        r = minimize(lambda x: (hs71_f(x), hs71_g(x)), [1, 5, 5, 1], jac=True, bounds=[(1, 5)] * 4,
                     constraints=cons)
        self.check_hs71(r, 1e-5)
        self.assertEqual([c.size for c in r.constr], [1, 1])

    def test_sparse_jacobian_pattern(self):
        # two constraints, each on one variable: x1 + x2 >= ... written as a sparse (2, 3) Jacobian
        def c(x):
            return np.array([x[0] ** 2, x[2] ** 2])

        def cj(x):
            return Coo([0, 1], [0, 2], [2 * x[0], 2 * x[2]], (2, 3))
        r = minimize(lambda x: np.sum((x - 3) ** 2), [0.5, 0.5, 0.5], jac=lambda x: 2 * (x - 3),
                     constraints=NonlinearConstraint(c, -np.inf, [1, 4], jac=cj))
        self.assertTrue(r.success, r.message)
        np.testing.assert_allclose(r.x, [1, 3, 2], atol=1e-6)

    def test_linear_constraint(self):
        r = minimize(lambda x: x @ x, [3.0, 3.0], jac=lambda x: 2 * x, constraints=LinearConstraint([[1, 1]], 2))
        np.testing.assert_allclose(r.x, [1, 1], atol=1e-6)
        self.assertAlmostEqual(r.v[0][0], 2.0, places=5)
        r = minimize(lambda x: x @ x, [3.0, 3.0], jac=lambda x: 2 * x,
                     constraints=LinearConstraint(Coo([0, 0], [0, 1], [1, 1], (1, 2)), 2))
        np.testing.assert_allclose(r.x, [1, 1], atol=1e-6)

    def test_exact_hessian(self):
        con = NonlinearConstraint(lambda x: x @ x, -np.inf, 1.5, jac=lambda x: 2 * x,
                                  hess=lambda x, v: 2 * v[0] * np.eye(2))
        r = minimize(rosen, [-1.2, 1], jac=rosen_g, hess=rosen_h, constraints=con)
        rb = minimize(rosen, [-1.2, 1], jac=rosen_g, constraints=con)
        self.assertTrue(r.success, r.message)
        self.assertGreater(r.nhev, 0)
        np.testing.assert_allclose(r.x, rb.x, atol=1e-6)
        with self.assertRaises(ValueError):   # a nonlinear constraint without its Hessian
            minimize(rosen, [-1.2, 1], jac=rosen_g, hess=rosen_h,
                     constraints=NonlinearConstraint(lambda x: x @ x, -np.inf, 1.5, jac=lambda x: 2 * x))

    def test_inertia_control(self):
        # (min -x1 x2 on the unit disc: the Hessian is indefinite)
        con = NonlinearConstraint(lambda x: x @ x, -np.inf, 1.0, jac=lambda x: 2 * x,
                                  hess=lambda x, v: 2 * v[0] * np.eye(2))
        r = minimize(lambda x: -x[0] * x[1], [0.3, 0.6], jac=lambda x: -x[::-1],
                     hess=lambda x: np.array([[0.0, -1.0], [-1.0, 0.0]]), bounds=[(0, 2), (0, 2)], constraints=con,
                     options={'inertia_control': True})
        if 'MUMPS' in r.message:
            # (the extension was built without MUMPS: the option is invalid input)
            self.assertFalse(r.success)
            self.assertEqual(r.n_factorizations, 0)
        else:
            self.assertTrue(r.success, r.message)
            self.assertAlmostEqual(r.fun, -0.5, places=6)
            self.assertGreater(r.n_factorizations, 0)

    def test_callback(self):
        seen = []
        r = minimize(rosen, [-1.2, 1], jac=rosen_g, callback=lambda intermediate_result: seen.append(
            intermediate_result.nit))
        self.assertTrue(r.success)
        self.assertEqual(seen, list(range(1, r.nit + 1)))
        xs = []
        minimize(rosen, [-1.2, 1], jac=rosen_g, callback=xs.append)       # old style: callback(x)
        self.assertIsInstance(xs[0], np.ndarray)

        def stop(intermediate_result):
            if intermediate_result.nit == 3:
                raise StopIteration
        r = minimize(rosen, [-1.2, 1], jac=rosen_g, callback=stop)
        self.assertFalse(r.success)
        self.assertEqual(r.nit, 3)

    def test_options(self):
        r = minimize(rosen, [-1.2, 1], jac=rosen_g, options={'maxiter': 3})
        self.assertFalse(r.success)
        self.assertEqual(r.nit, 3)
        from sqpopt_options import schema
        values = schema.default_values()
        values['options']['max_iter'] = 4
        values['linesearch']['filter']['delta'] = 2.0
        r = minimize(rosen, [-1.2, 1], jac=rosen_g, options=values)   # the options dialog's dict
        self.assertEqual(r.nit, 4)
        con = NonlinearConstraint(lambda x: x @ x, -np.inf, 1.5, jac=lambda x: 2 * x)
        r = minimize(rosen, [-1.2, 1], jac=rosen_g, constraints=con,
                     options={'hessian_mode': 'sr1', 'linesearch.filter.delta': 1})
        self.assertTrue(r.success, r.message)
        for bad in ({'nope': 1}, {'max_iter': 1.5}, {'hessian_mode': 'foo'}, {'opt_tol': 1e-8}, {'ktol': -1.0},
                    {'max_evals': 10**10}):
            with self.assertRaises(ValueError, msg=str(bad)):
                minimize(rosen, [0, 0], jac=rosen_g, options=bad)

    def test_lambda0(self):
        con = NonlinearConstraint(hs71_c, [40, 25], [40, np.inf], jac=hs71_cj)
        r = minimize(hs71_f, [1, 5, 5, 1], jac=hs71_g, bounds=[(1, 5)] * 4, constraints=con)
        warm = minimize(hs71_f, r.x, jac=hs71_g, bounds=[(1, 5)] * 4, constraints=con, lambda0=r.v)
        cold = minimize(hs71_f, r.x, jac=hs71_g, bounds=[(1, 5)] * 4, constraints=con)
        self.check_hs71(warm)
        self.assertLessEqual(warm.nit, cold.nit)
        self.assertEqual(warm.nit, 1)                  # (already converged, with the right multipliers)
        flat = minimize(hs71_f, r.x, jac=hs71_g, bounds=[(1, 5)] * 4, constraints=con,
                        lambda0=np.concatenate(r.v))   # one array of all the rows
        self.assertEqual(flat.nit, warm.nit)
        with self.assertRaises(ValueError):
            minimize(hs71_f, r.x, jac=hs71_g, constraints=con, lambda0=[1.0])

    def test_output_file(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            path = pathlib.Path(d) / 'sqpopt.log'
            r = minimize(rosen, [-1.2, 1], jac=rosen_g, options={'disp': True}, output_file=path)
            self.assertTrue(r.success)
            self.assertIn('sqpopt: status 0', path.read_text())
            r = minimize(rosen, [-1.2, 1], jac=rosen_g, output_file=pathlib.Path(d) / 'no' / 'such' / 'dir')
            self.assertFalse(r.success)
            self.assertIn("can't be opened", r.message)

    def test_invalid_input(self):
        # (caught by the Fortran input validation)
        r = minimize(rosen, [0, 0], jac=rosen_g, bounds=[(1, 0), (None, None)])
        self.assertFalse(r.success)
        self.assertIn('bound', r.message)

    def test_exceptions_propagate(self):
        def boom(x):
            raise RuntimeError('user error')
        with self.assertRaises(RuntimeError):
            minimize(boom, [0.0, 0.0], jac=rosen_g)
        with self.assertRaises(RuntimeError):
            minimize(rosen, [0.0, 0.0], jac=lambda x: boom(x))

    def test_nested(self):
        def outer(x):
            inner = minimize(lambda y: (y[0] - x[0]) ** 2, [0.0], jac=lambda y: np.array([2 * (y[0] - x[0])]))
            return inner.fun + (x[0] - 1) ** 2
        r = minimize(outer, [3.0], jac=lambda x: np.array([2 * (x[0] - 1)]))   # (the inner minimum is 0)
        self.assertTrue(r.success, r.message)
        self.assertAlmostEqual(r.x[0], 1.0, places=5)


if __name__ == '__main__':
    unittest.main()
