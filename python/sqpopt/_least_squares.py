"""`least_squares`: nonlinear least-squares problems, solved with `minimize` in Schittkowski's form.

The problem ``min 1/2 sum(r(x)**2)`` (with bounds and constraints) is given to the solver as
``min 1/2 z'z  s.t.  r(x) - z = 0`` in the variables ``(x, z)``, as the Fortran ``sqpopt_nlls_module``
does (its documentation has the reasons). Everything here is built on `minimize`: the solver itself has no
least-squares mode.
"""

from __future__ import annotations

from typing import Any, Callable, Sequence

import numpy as np

from ._constraints import Bounds, LinearConstraint, NonlinearConstraint
from ._minimize import (INFINITY, OptimizeResult, _Constraint, _is_sparse, _pattern, _values, _variable_bounds,
                        minimize)


class _Coo:
    """a sparse matrix in coordinate form, with what `minimize` needs of one (``tocoo()``, ``toarray()``)"""

    def __init__(self, row, col, data, shape):
        self.row, self.col, self.data, self.shape = row, col, data, shape

    def tocoo(self):
        return self

    def toarray(self):
        a = np.zeros(self.shape)
        np.add.at(a, (self.row, self.col), self.data)
        return a


def least_squares(fun: Callable, x0, jac: Callable, bounds=None, constraints: Sequence | Any = (), args=(),
                  jac_sparsity=None, gauss_newton: bool = False, tol: float | None = None,
                  options: dict | None = None, output_file=None, max_step=None,
                  diagnostics_file=None) -> OptimizeResult:
    """Minimize ``1/2 sum(fun(x)**2)`` subject to bounds and constraints, with SQPOPT.

    Use it for data fitting, and for systems of equations that can't all hold (more equations than
    unknowns), which `minimize` handles badly as equality constraints.

    fun : callable
        ``fun(x, *args)`` returns the residuals (a 1-D array of ``l`` values).
    x0 : array_like
        The starting point (moved inside the bounds).
    jac : callable
        ``jac(x, *args)`` returns the residuals' Jacobian, dense ``(l, n)`` or sparse (any object with
        ``tocoo()``). Required: there are no finite differences.
    bounds, constraints, tol, options, output_file, diagnostics_file
        As for `minimize`. The constraints are on ``x``.
    args : tuple
        Extra arguments of `fun` and `jac`.
    jac_sparsity : optional
        The sparsity pattern of the residuals' Jacobian: an ``(l, n)`` array or sparse matrix whose nonzeros
        are the entries that may be nonzero. Without it, the pattern is that of ``jac``'s result at ``x0`` if
        that is sparse, else dense.
    gauss_newton : bool
        Use the Gauss-Newton Hessian (the second derivatives of the residuals and of the constraints are
        neglected) instead of the solver's quasi-Newton approximation. Usually takes fewer iterations,
        mostly so when the residuals are small at the solution. It is given to the solver as a dense matrix
        of order ``n + l``, so it is for small problems (the Fortran interface gives it as a sparse one).
    max_step : float or array, optional
        The largest change of each variable ``x`` in one major iteration (see `minimize`).

    Returns
    -------
    OptimizeResult
        With ``x``, ``fun`` (the residuals at ``x``), ``cost`` (``1/2 sum(fun**2)``), ``success``,
        ``status``, ``message``, ``nit``, ``nfev`` and ``njev`` (calls of `fun` and of `jac`), ``constr`` and
        ``v`` (the values and multipliers of each constraint object), ``z`` (the multipliers of the bounds),
        ``constr_violation``, ``kkt_error``, ``execution_time``, ``diagnostics``, and ``nlp`` (the result of
        `minimize` for the transformed problem, whose variables are ``x`` followed by the ``l`` auxiliary ones).
    """
    args = tuple(args) if isinstance(args, (tuple, list)) else (args,)
    if not callable(jac):
        raise ValueError('the Jacobian of the residuals is required: jac must be a callable')
    x0 = np.atleast_1d(np.asarray(x0, dtype=float)).ravel()
    n = x0.size
    x_lb, x_ub = _variable_bounds(bounds, n)
    x0 = np.clip(x0, x_lb, x_ub)
    count = {'fun': 0, 'jac': 0}
    last: list = [None, None]           # the last point where the residuals were evaluated, and their values

    def residuals(x: np.ndarray) -> np.ndarray:
        if last[0] is None or not np.array_equal(last[0], x):
            count['fun'] += 1
            r = np.atleast_1d(np.asarray(fun(x, *args), dtype=float)).ravel()
            last[0], last[1] = x.copy(), r
        return last[1]

    r0 = residuals(x0)
    ell = r0.size
    nn = n + ell

    # the pattern of the residuals' Jacobian, then `-1` for each auxiliary variable
    sparsity = jac_sparsity
    if sparsity is None:
        j0 = jac(x0, *args)
        count['jac'] += 1
        if _is_sparse(j0):
            sparsity = j0
    if sparsity is None:
        rows, cols = np.repeat(np.arange(ell, dtype=np.int64), n), np.tile(np.arange(n, dtype=np.int64), ell)
    else:
        rows, cols = _pattern(sparsity, ell, n)
    diag = np.arange(ell, dtype=np.int64)
    r_rows, r_cols = np.concatenate([rows, diag]), np.concatenate([cols, n + diag])

    def residual_jac(y: np.ndarray) -> _Coo:
        count['jac'] += 1
        vals = _values(jac(y[:n], *args), rows, cols, ell, n, 'jac')
        return _Coo(r_rows, r_cols, np.concatenate([vals, -np.ones(ell)]), (ell, nn))

    def no_curvature(y, v):
        return np.zeros((nn, nn))

    def extended(c: _Constraint) -> NonlinearConstraint:
        # (a constraint on `x`, as one on the variables of the transformed problem)
        def cjac(y):
            vals = c.jac_val if c.linear else _values(c.jac(y[:n]), c.rows, c.cols, c.k, n, 'constraint jac')
            return _Coo(c.rows, c.cols, vals, (c.k, nn))
        return NonlinearConstraint(lambda y: c.fun(y[:n]), c.lb, c.ub, jac=cjac, hess=no_curvature,
                                   jac_sparsity=_Coo(c.rows, c.cols, np.ones(c.rows.size), (c.k, nn)))

    if isinstance(constraints, (dict, LinearConstraint)) or hasattr(constraints, 'fun') or hasattr(constraints, 'A'):
        constraints = [constraints]
    cons = [_Constraint(c, x0, n) for c in constraints]
    ctol = 1.0e-8
    feasible = all(np.all(v >= c.lb - ctol) and np.all(v <= c.ub + ctol) for c in cons for v in [c.fun(x0)])
    all_cons = [NonlinearConstraint(lambda y: residuals(y[:n]) - y[n:], 0.0, 0.0, jac=residual_jac, hess=no_curvature,
                                    jac_sparsity=_Coo(r_rows, r_cols, np.ones(r_rows.size), (ell, nn)))]
    all_cons += [extended(c) for c in cons]

    # the auxiliary variables start from the residuals if the constraints hold at `x0`, else from zero
    y0 = np.concatenate([x0, r0 if feasible else np.zeros(ell)])
    if max_step is not None:
        max_step = np.concatenate([np.broadcast_to(np.asarray(max_step, dtype=float), (n,)), np.full(ell, np.inf)])

    def hessian(y):
        h = np.zeros((nn, nn))
        h[n:, n:] = np.eye(ell)
        return h

    res = minimize(lambda y: 0.5 * float(y[n:] @ y[n:]), y0, jac=lambda y: np.concatenate([np.zeros(n), y[n:]]),
                   hess=hessian if gauss_newton else None,
                   bounds=Bounds(np.concatenate([x_lb, np.full(ell, -INFINITY)]),
                                 np.concatenate([x_ub, np.full(ell, INFINITY)])),
                   constraints=all_cons, tol=tol, options=options, output_file=output_file, max_step=max_step,
                   diagnostics_file=diagnostics_file)
    x = res.x[:n]
    r = res.x[n:] + res.constr[0]       # (the auxiliary variables, plus what is left of `r - z`)
    return OptimizeResult(
        x=x, fun=r, cost=0.5 * float(r @ r), success=res.success, status=res.status, message=res.message,
        nit=res.nit, nfev=count['fun'], njev=count['jac'], constr=res.constr[1:], v=res.v[1:], z=res.z[:n],
        constr_violation=res.constr_violation, kkt_error=res.kkt_error, execution_time=res.execution_time,
        diagnostics=res.diagnostics, nlp=res)
