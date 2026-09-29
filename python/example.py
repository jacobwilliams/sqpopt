"""A basic example of the SQPOPT Python bindings: Hock-Schittkowski problem 71.

    minimize    x1*x4*(x1 + x2 + x3) + x3
    subject to  x1^2 + x2^2 + x3^2 + x4^2 = 40
                x1*x2*x3*x4 >= 25
                1 <= x1, x2, x3, x4 <= 5

    solution:   x* = (1, 4.743, 3.821, 1.379),  f* = 17.014

Build the extension first, then run this script:

    pixi run build-python
    pixi run python python/example.py
"""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))   # (to import sqpopt from this folder)

from sqpopt import NonlinearConstraint, minimize  # noqa: E402


def objective(x):
    return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]


def gradient(x):
    return np.array([x[3] * (2 * x[0] + x[1] + x[2]), x[0] * x[3], x[0] * x[3] + 1, x[0] * (x[0] + x[1] + x[2])])


def constraints(x):
    return np.array([np.sum(x ** 2), np.prod(x)])


def constraints_jacobian(x):
    return np.array([2 * x,
                     [x[1] * x[2] * x[3], x[0] * x[2] * x[3], x[0] * x[1] * x[3], x[0] * x[1] * x[2]]])


x0 = [1.0, 5.0, 5.0, 1.0]
bounds = [(1, 5)] * 4

# 1. with analytic derivatives, printing each iteration:
con = NonlinearConstraint(constraints, lb=[40, 25], ub=[40, np.inf], jac=constraints_jacobian)
res = minimize(objective, x0, jac=gradient, bounds=bounds, constraints=con,
               callback=lambda intermediate_result: print(f'  iteration {intermediate_result.nit}: '
                                                          f'f = {intermediate_result.fun:.8f}'))
print('with analytic derivatives:')
print(res)

# 2. SLSQP-style constraint dicts ('ineq' means fun(x) >= 0), with some options:
cons = [{'type': 'eq', 'fun': lambda x: np.sum(x ** 2) - 40, 'jac': lambda x: 2 * x},
        {'type': 'ineq', 'fun': lambda x: np.prod(x) - 25, 'jac': lambda x: constraints_jacobian(x)[1]}]
res2 = minimize(objective, x0, jac=gradient, bounds=bounds, constraints=cons, options={'ktol': 1e-8, 'max_iter': 50})
print('\nwith constraint dicts:')
print(f'  x = {res2.x}, f = {res2.fun:.8f}, {res2.message} ({res2.nit} iterations, {res2.nfev} objective calls)')

assert res.success and res2.success
assert abs(res.fun - 17.0140173) < 1e-6 and abs(res2.fun - 17.0140173) < 1e-6
