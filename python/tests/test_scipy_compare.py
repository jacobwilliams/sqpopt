"""SQPOPT's limited-memory BFGS method against scipy's (``L-BFGS-B``), on functions without constraints.

On such a function SQPOPT's QP subproblem has no constraints, and its solution is the quasi-Newton step, so
SQPOPT is then an L-BFGS method like scipy's: it should need about as many iterations. Without the
unconstrained step of the QP front end (``qp_solver%unconstrained_step``), each of those iterations runs
an active-set QP solver on a problem with every variable free, which costs much more (a time that isn't
tested here: timings aren't reproducible enough for a unit test).

The functions are two of ``test/scalable_functions.f90``: the chained Rosenbrock function, which needs
about 5 iterations per variable from its standard starting point, and the ill-conditioned quadratic
``trid``. Skipped if the extension isn't built, or scipy isn't installed.

Run from the repository root:  pixi run python -m unittest discover -s python/tests
"""

import pathlib
import sys
import unittest

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

try:
    from sqpopt import minimize
    from sqpopt._minimize import _native
    _native()
    BUILT = True
except ImportError:
    BUILT = False

try:
    from scipy.optimize import minimize as scipy_minimize
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

MEMORY = 20          # stored pairs, in both methods
MAX_ITER = 100000    # iteration limit, in both methods
RATIO = 1.5          # SQPOPT may take at most this many times scipy's iterations


def rosenbrock(x):
    """the chained Rosenbrock function: sum 100 (x[i+1] - x[i]^2)^2 + (x[i] - 1)^2"""
    return np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (x[:-1] - 1) ** 2)


def rosenbrock_g(x):
    t = x[1:] - x[:-1] ** 2
    g = np.zeros_like(x)
    g[:-1] = -400 * t * x[:-1] + 2 * (x[:-1] - 1)
    g[1:] += 200 * t
    return g


def rosenbrock_x0(n):
    x0 = np.full(n, -1.2)
    x0[-1] = 1.0
    return x0


def trid(x):
    """sum (x[i] - 1)^2 - sum x[i] x[i-1], term by term (see test/scalable_functions.f90)"""
    d = np.diff(x, prepend=0.0)
    return np.sum(x * d - 2 * x + 1)


def trid_g(x):
    g = 2 * (x - 1)
    g[1:] -= x[:-1]
    g[:-1] -= x[1:]
    return g


def trid_min(n):
    return -n * (n + 4) * (n - 1) / 6


@unittest.skipUnless(BUILT, 'the sqpopt extension is not built (pixi run build-python)')
@unittest.skipUnless(HAVE_SCIPY, 'scipy is not installed')
class TestScipyCompare(unittest.TestCase):

    def compare(self, fun, jac, x0, f_min, f_tol):
        """solve with both methods; both must reach the minimum, and SQPOPT must not need many more
        iterations than scipy, and must solve its QPs by the unconstrained step"""
        n = x0.size
        ref = scipy_minimize(fun, x0, jac=jac, method='L-BFGS-B',
                             options={'maxcor': MEMORY, 'maxiter': MAX_ITER, 'maxfun': 2 * MAX_ITER,
                                      'ftol': 1e-15, 'gtol': 1e-6})
        r = minimize(fun, x0, jac=jac, options={'lbfgs_memory': MEMORY, 'maxiter': MAX_ITER})
        off = minimize(fun, x0, jac=jac, options={'lbfgs_memory': MEMORY, 'maxiter': MAX_ITER,
                                                  'qp_solver%unconstrained_step': False})
        print(f'{fun.__name__:>10} n={n:4d}: iterations scipy {ref.nit:5d}, sqpopt {r.nit:5d} '
              f'(without the unconstrained step {off.nit:5d}); sqpopt time {r.execution_time:.3f} s '
              f'(without {off.execution_time:.3f} s)')
        self.assertLessEqual(abs(ref.fun - f_min), f_tol, 'scipy did not reach the minimum')
        self.assertTrue(r.success, r.message)
        self.assertTrue(off.success, off.message)
        self.assertLessEqual(abs(r.fun - f_min), f_tol)
        self.assertLessEqual(abs(off.fun - f_min), f_tol)
        self.assertLessEqual(r.nit, RATIO * ref.nit)
        # every QP is solved by the unconstrained step (there are no constraints, and no bounds)...
        self.assertGreater(r.n_qp_solves, 0)
        self.assertEqual(r.n_unconstrained_qp, r.n_qp_solves)
        self.assertEqual(r.n_qp_iterations, 0)
        # ...which gives the same steps as the QP solver, to roundoff, so about the same iterations
        self.assertEqual(off.n_unconstrained_qp, 0)
        self.assertLessEqual(r.nit, RATIO * off.nit)
        self.assertLessEqual(off.nit, RATIO * r.nit)

    def test_rosenbrock(self):
        for n in (50, 300):
            with self.subTest(n=n):
                self.compare(rosenbrock, rosenbrock_g, rosenbrock_x0(n), 0.0, 1e-8)

    def test_trid(self):
        for n in (50, 300):
            with self.subTest(n=n):
                self.compare(trid, trid_g, np.zeros(n), trid_min(n), 1e-6 * abs(trid_min(n)))


if __name__ == '__main__':
    unittest.main()
