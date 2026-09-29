"""Python bindings for SQPOPT, with a ``scipy.optimize.minimize``-like interface::

    from sqpopt import minimize, NonlinearConstraint

    res = minimize(fun, x0, jac=grad, bounds=[(1, 5)] * 4,
                   constraints=[NonlinearConstraint(cons, [40, 25], [40, np.inf], jac=cons_jac)])
    print(res.x, res.fun, res.message)

The native extension is built with f2py (see ``_build.py``)::

    pixi run python python/sqpopt/_build.py
"""

from ._constraints import Bounds, LinearConstraint, NonlinearConstraint
from ._minimize import INFINITY, OptimizeResult, minimize

__all__ = ['minimize', 'OptimizeResult', 'Bounds', 'NonlinearConstraint', 'LinearConstraint', 'INFINITY']
