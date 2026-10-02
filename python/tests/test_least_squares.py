"""Tests of ``sqpopt.least_squares``; skipped if the extension isn't built.

Build it first:  pixi run build-python
Run from the repository root:  pixi run python -m unittest discover -s python/tests
"""

import pathlib
import sys
import unittest

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

try:
    from sqpopt import LinearConstraint, NonlinearConstraint, least_squares
    from sqpopt._minimize import _native
    _native()
    BUILT = True
except ImportError:
    BUILT = False

# Bard's problem: 15 residuals in 3 variables, sum of squares 8.21487e-3 at the solution
BARD_Y = np.array([0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.10, 4.39])
BARD_U = np.arange(1, 16.0)
BARD_V = 16 - BARD_U
BARD_W = np.minimum(BARD_U, BARD_V)
BARD_X = np.array([0.0824106, 1.13304, 2.34370])


def bard(x):
    return BARD_Y - (x[0] + BARD_U / (BARD_V * x[1] + BARD_W * x[2]))


def bard_jac(x):
    d = (BARD_V * x[1] + BARD_W * x[2]) ** 2
    return np.column_stack([-np.ones(15), BARD_U * BARD_V / d, BARD_U * BARD_W / d])


# a line x0 + x1*t through four points
LINE_T = np.array([0.0, 1.0, 2.0, 3.0])
LINE_Y = np.array([0.1, 0.9, 2.1, 2.9])


def line(x, t=LINE_T, y=LINE_Y):
    return x[0] + x[1] * t - y


def line_jac(x, t=LINE_T, y=LINE_Y):
    return np.column_stack([np.ones(t.size), t])


class _Sparse:
    """a minimal sparse matrix (what the bindings need of one: ``tocoo()``)"""

    def __init__(self, a):
        a = np.asarray(a, dtype=float)
        self.row, self.col = np.nonzero(a)
        self.data = a[self.row, self.col]
        self.shape = a.shape

    def tocoo(self):
        return self


@unittest.skipUnless(BUILT, 'the sqpopt extension is not built (pixi run build-python)')
class TestLeastSquares(unittest.TestCase):

    def test_bard(self):
        for gauss_newton in (False, True):
            r = least_squares(bard, [1.0, 1.0, 1.0], bard_jac, gauss_newton=gauss_newton)
            self.assertEqual(r.status, 0, r.message)
            self.assertTrue(r.success)
            self.assertAlmostEqual(2 * r.cost, 8.21487e-3, places=7)
            np.testing.assert_allclose(r.x, BARD_X, atol=1e-4)
            np.testing.assert_allclose(r.fun, bard(r.x), atol=1e-8)
            # (the residuals at the starting point are evaluated once, so one evaluation per iteration)
            self.assertEqual(r.nfev, r.nit)
            self.assertEqual(r.nlp.x.size, 3 + 15)
        self.assertLess(r.nit, 15)      # (Gauss-Newton: 6 iterations, against 18 with L-BFGS)

    def test_sparse_jacobian(self):
        dense = least_squares(bard, [1.0, 1.0, 1.0], bard_jac)
        sparse = least_squares(bard, [1.0, 1.0, 1.0], lambda x: _Sparse(bard_jac(x)))
        given = least_squares(bard, [1.0, 1.0, 1.0], bard_jac, jac_sparsity=np.ones((15, 3)))
        np.testing.assert_allclose(sparse.x, dense.x, rtol=1e-12)
        np.testing.assert_allclose(given.x, dense.x, rtol=1e-12)

    def test_bounds_and_constraints(self):
        # on x0 + x1 = 1 the sum of squares is least at x1 = 5.8/6, so the bound x1 <= 0.5 is active
        for cons in (LinearConstraint([[1.0, 1.0]], 1.0, 1.0),
                     NonlinearConstraint(lambda x: x[0] + x[1], 1.0, 1.0, jac=lambda x: np.array([[1.0, 1.0]])),
                     {'type': 'eq', 'fun': lambda x: x[0] + x[1] - 1.0, 'jac': lambda x: np.array([[1.0, 1.0]])}):
            for gauss_newton in (False, True):
                r = least_squares(line, [0.0, 0.0], line_jac, bounds=[(None, None), (None, 0.5)], constraints=cons,
                                  gauss_newton=gauss_newton)
                self.assertEqual(r.status, 0, r.message)
                np.testing.assert_allclose(r.x, [0.5, 0.5], atol=1e-6)
                self.assertAlmostEqual(2 * r.cost, float(np.sum(line(r.x) ** 2)), places=8)
                self.assertEqual(len(r.constr), 1)
                self.assertEqual(len(r.v), 1)
                self.assertEqual(r.z.size, 2)
                self.assertLess(r.z[1], 0.0)        # (the multiplier of an active upper bound)
                self.assertEqual(r.z[0], 0.0)

    def test_args_and_max_step(self):
        t, y = np.array([0.0, 1.0, 2.0]), np.array([1.0, 3.0, 5.0])
        r = least_squares(line, [0.0, 0.0], line_jac, args=(t, y), max_step=0.5)
        self.assertEqual(r.status, 0, r.message)
        np.testing.assert_allclose(r.x, [1.0, 2.0], atol=1e-6)
        self.assertLess(r.cost, 1e-12)
        self.assertGreaterEqual(r.nit, 4)       # (the solution is 2 away, in steps of at most 0.5 per variable)

    def test_infeasible_start(self):
        # (the constraint doesn't hold at the starting point: the auxiliary variables start from zero)
        r = least_squares(line, [5.0, 5.0], line_jac, constraints=LinearConstraint([[1.0, 1.0]], 1.0, 1.0))
        self.assertEqual(r.status, 0, r.message)
        np.testing.assert_allclose(r.x, [1.0 - 5.8 / 6.0, 5.8 / 6.0], atol=1e-6)

    def test_needs_jacobian(self):
        with self.assertRaises(ValueError):
            least_squares(bard, [1.0, 1.0, 1.0], None)


if __name__ == '__main__':
    unittest.main()
