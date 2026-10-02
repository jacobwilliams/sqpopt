"""``minimize``: a ``scipy.optimize.minimize``-like interface to SQPOPT."""

from __future__ import annotations

import inspect
from typing import Any, Callable, Sequence

import numpy as np

from sqpopt_options import schema

from ._constraints import LinearConstraint

INFINITY = schema.SQPOPT_INFINITY   #: bounds of this magnitude or more are "no bound"


class OptimizeResult(dict):
    """the result of `minimize`: a dict whose keys are also attributes (as ``scipy.optimize.OptimizeResult``)"""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if not self:
            return self.__class__.__name__ + '()'
        width = max(map(len, self)) + 1
        return '\n'.join(f'{k.rjust(width)}: {v!r}' for k, v in self.items())

    def __dir__(self):
        return list(self.keys())


def _native():
    """the extension module, checked against the options schema it was built with"""
    try:
        from . import _sqpopt
    except ImportError as e:
        raise ImportError('the sqpopt extension is not built: run `pixi run python python/sqpopt/_build.py`') from e
    from ._build import options_signature
    _, _, n_options, signature = _sqpopt.sqpopt_py_info()
    if n_options != len(schema.OPTIONS) or signature != options_signature():
        raise ImportError('the sqpopt extension was built with different options: rebuild it '
                          '(`pixi run python python/sqpopt/_build.py`)')
    return _sqpopt


# ------------------------------------------------------------------------------------------------------
# options

_ALIASES = {'maxiter': 'options%max_iter', 'disp': 'options%print_level'}


def _option_path(key: str) -> str:
    """the Fortran reference of an option given by its reference (``'linesearch%filter%delta'`` or
    ``'linesearch.filter.delta'``), a scipy alias, or its (unambiguous) name (``'ktol'``)"""
    key = key.replace('.', '%')
    if key in schema.OPTIONS:
        return key
    if key in _ALIASES:
        return _ALIASES[key]
    matches = [ref for ref, o in schema.OPTIONS.items() if o.name == key]
    if 'options%' + key in matches:
        return 'options%' + key
    if len(matches) == 1:
        return matches[0]
    if matches:
        raise ValueError(f'ambiguous option {key!r}: use one of ' + ', '.join(matches))
    raise ValueError(f'unknown option {key!r}')


def _option_value(o: schema.Option, value: Any) -> Any:
    """`value` for option `o`, converted (choices may be given by their Fortran name, e.g. ``'exact'`` or
    ``'sqpopt_hessian_exact'``) and checked"""
    if o.kind == 'choice' and isinstance(value, str):
        for c in o.choices:
            if c.name and (value == c.name or c.name.endswith('_' + value)):
                value = c.value
                break
        else:
            raise ValueError(f'{o.fortran}: unknown value {value!r}')
    elif o.kind == 'bool':
        value = bool(value)
    elif o.kind == 'float' and isinstance(value, (int, np.integer, np.floating)) and not isinstance(value, bool):
        value = float(value)
    elif o.kind in ('int', 'choice') and isinstance(value, np.integer):
        value = int(value)
    msg = o.check(value)
    if not msg and o.kind in ('int', 'choice') and not -2**31 <= value < 2**31:
        msg = 'must fit in a 32-bit integer'
    if msg:
        raise ValueError(f'{o.fortran}: {msg}')
    return value


def _resolve_options(options: dict | None) -> dict[str, Any]:
    """the options as ``{fortran reference: value}``, from a nested dict (as the options dialog returns), a
    flat one, or a mix"""
    flat: dict[str, Any] = {}
    for key, value in (options or {}).items():
        if isinstance(value, dict):
            for sub, v in schema.flatten({key: value}, only_known=False).items():
                flat[_option_path(sub)] = v
        elif key == 'disp':
            flat['options%print_level'] = 1 if value is True else (0 if value is False else value)
        else:
            flat[_option_path(key)] = value
    return {ref: _option_value(schema.OPTIONS[ref], v) for ref, v in flat.items()}


def _choice(ref: str, name: str) -> int:
    """the value of choice `name` of option `ref`"""
    return next(c.value for c in schema.OPTIONS[ref].choices if c.name == name)


# ------------------------------------------------------------------------------------------------------
# problem data

def _bound_array(b, n: int, default: float) -> np.ndarray:
    """a bound as an array of size `n`, with infinite (or ``None``) values as +-`INFINITY`"""
    if b is None:
        b = default
    a = np.array([default if v is None else v for v in np.broadcast_to(np.asarray(b, dtype=object), (n,))],
                 dtype=float)
    return np.clip(a, -INFINITY, INFINITY)


def _variable_bounds(bounds, n: int) -> tuple[np.ndarray, np.ndarray]:
    """the variable bounds from ``None``, an object with ``lb``/``ub``, or a sequence of ``(lb, ub)`` pairs"""
    if bounds is None:
        return np.full(n, -INFINITY), np.full(n, INFINITY)
    if hasattr(bounds, 'lb') and hasattr(bounds, 'ub'):
        return _bound_array(bounds.lb, n, -INFINITY), _bound_array(bounds.ub, n, INFINITY)
    pairs = list(bounds)
    if len(pairs) != n:
        raise ValueError(f'bounds must have {n} (lb, ub) pairs')
    return (_bound_array([p[0] for p in pairs], n, -INFINITY), _bound_array([p[1] for p in pairs], n, INFINITY))


def _is_sparse(a) -> bool:
    return hasattr(a, 'tocoo') and not isinstance(a, np.ndarray)


def _pattern(a, k: int, n: int) -> tuple[np.ndarray, np.ndarray]:
    """the (row, column) indices of the nonzeros of a dense or sparse ``(k, n)`` matrix, in row-major order"""
    if _is_sparse(a):
        coo = a.tocoo()
        keys = np.unique(np.asarray(coo.row, dtype=np.int64) * n + np.asarray(coo.col, dtype=np.int64))
    else:
        d = np.asarray(a).reshape(k, n)
        keys = np.flatnonzero(d)
    return (keys // n).astype(np.int64), (keys % n).astype(np.int64)


def _values(a, rows: np.ndarray, cols: np.ndarray, k: int, n: int, what: str) -> np.ndarray:
    """the entries of a dense or sparse ``(k, n)`` matrix at the pattern `rows`/`cols` (in row-major order)"""
    if not _is_sparse(a):
        return np.asarray(a, dtype=float).reshape(k, n)[rows, cols]
    coo = a.tocoo()
    keys = np.asarray(coo.row, dtype=np.int64) * n + np.asarray(coo.col, dtype=np.int64)
    pattern = rows * n + cols
    idx = np.searchsorted(pattern, keys)
    inside = (idx < pattern.size) & (pattern[np.minimum(idx, pattern.size - 1)] == keys)
    data = np.asarray(coo.data, dtype=float)
    if np.any(data[~inside] != 0.0):
        raise ValueError(f'{what}: nonzero outside the sparsity pattern')
    vals = np.zeros(pattern.size)
    np.add.at(vals, idx[inside], data[inside])
    return vals


class _Constraint:
    """one constraint object (a block of rows ``lb <= fun(x) <= ub``)"""

    def __init__(self, spec, x0: np.ndarray, n: int):
        self.hess = None
        self.linear = False
        if isinstance(spec, dict):              # SLSQP style: {'type': 'eq'|'ineq', 'fun', 'jac', 'args'}
            kind = spec.get('type')
            if kind not in ('eq', 'ineq'):
                raise ValueError("a constraint dict's 'type' must be 'eq' or 'ineq'")
            cargs = tuple(spec.get('args', ()))
            f = spec['fun']
            self.fun = lambda x: np.atleast_1d(np.asarray(f(x, *cargs), dtype=float)).ravel()
            j = spec.get('jac')
            if not callable(j):
                raise ValueError("a constraint dict needs its Jacobian: 'jac'")
            self.jac = lambda x: j(x, *cargs)
            lb, ub = 0.0, (0.0 if kind == 'eq' else np.inf)
            sparsity = None
        elif hasattr(spec, 'A'):                # LinearConstraint
            a = spec.A
            if _is_sparse(a):
                coo = a.tocoo()
                rows, cols = np.asarray(coo.row), np.asarray(coo.col)
                data, k = np.asarray(coo.data, dtype=float), a.shape[0]
                self.fun = lambda x: np.bincount(rows, data * x[cols], minlength=k)
            else:
                a = np.atleast_2d(np.asarray(a, dtype=float))
                self.fun = lambda x: a @ x
            self.A = a
            self.jac = lambda x: self.A
            self.linear = True
            lb, ub = spec.lb, spec.ub
            sparsity = a
        else:                                   # NonlinearConstraint
            f = spec.fun
            self.fun = lambda x: np.atleast_1d(np.asarray(f(x), dtype=float)).ravel()
            j = getattr(spec, 'jac', None)
            if not callable(j):
                raise ValueError('a NonlinearConstraint needs its Jacobian: jac (a callable)')
            self.jac = j
            h = getattr(spec, 'hess', None)
            self.hess = h if callable(h) else None
            lb, ub = spec.lb, spec.ub
            sparsity = getattr(spec, 'jac_sparsity', None)
        c0 = self.fun(x0)
        self.k = c0.size
        self.lb = _bound_array(lb, self.k, -INFINITY)
        self.ub = _bound_array(ub, self.k, INFINITY)
        if sparsity is None:
            j0 = self.jac(x0)
            if _is_sparse(j0):
                sparsity = j0
        if sparsity is None:
            self.rows = np.repeat(np.arange(self.k, dtype=np.int64), n)
            self.cols = np.tile(np.arange(n, dtype=np.int64), self.k)
        else:
            self.rows, self.cols = _pattern(sparsity, self.k, n)
        self.jac_val = _values(self.A, self.rows, self.cols, self.k, n, 'LinearConstraint') if self.linear else None


# ------------------------------------------------------------------------------------------------------

class _Solve:
    """the state of one `minimize` call: the callbacks given to the Fortran solver"""

    def __init__(self, fun, x0, args, jac, hess, bounds, constraints, callback):
        self.x0 = np.atleast_1d(np.asarray(x0, dtype=float)).ravel().copy()
        self.n = n = self.x0.size
        self.args = tuple(args) if isinstance(args, (tuple, list)) else (args,)
        self.fun = fun
        if jac is True:
            self.grad_mode = 'combined'         # fun returns (f, g)
        elif callable(jac):
            self.grad_mode = 'callable'
        else:
            raise ValueError('the gradient is required: jac must be a callable, or True if fun returns it too')
        self.jac = jac
        if hess is not None and not callable(hess):
            raise ValueError('hess must be a callable (the Hessian of the objective)')
        self.hess_f = hess
        self.x_lb, self.x_ub = _variable_bounds(bounds, n)
        if isinstance(constraints, (dict, LinearConstraint)) or hasattr(constraints, 'fun') or \
                hasattr(constraints, 'A'):
            constraints = [constraints]
        self.cons = [_Constraint(c, self.x0, n) for c in constraints]
        if hess is not None:
            missing = [i for i, c in enumerate(self.cons) if not c.linear and c.hess is None]
            if missing:
                raise ValueError(f'with hess, every nonlinear constraint needs hess(x, v) (missing: {missing})')
        self.m = sum(c.k for c in self.cons)
        self.callback = callback
        self.callback_new_style = False
        if callback is not None:
            try:
                self.callback_new_style = 'intermediate_result' in inspect.signature(callback).parameters
            except (TypeError, ValueError):
                pass
        self.nfev = 0
        self.error: BaseException | None = None
        self.last_grad: tuple[np.ndarray, np.ndarray] | None = None

    # ---- evaluations

    def objective(self, x: np.ndarray) -> float:
        self.nfev += 1
        if self.grad_mode == 'combined':
            f, g = self.fun(x, *self.args)
            self.last_grad = (x.copy(), np.asarray(g, dtype=float).ravel().copy())
            return float(f)
        return float(self.fun(x, *self.args))

    def constraints(self, x: np.ndarray) -> np.ndarray:
        if not self.cons:
            return np.zeros(0)
        return np.concatenate([c.fun(x) for c in self.cons])

    def gradient(self, x: np.ndarray) -> np.ndarray:
        if self.grad_mode == 'callable':
            return np.asarray(self.jac(x, *self.args), dtype=float).ravel()
        if self.last_grad is None or not np.array_equal(self.last_grad[0], x):
            self.objective(x)
        return self.last_grad[1]

    def jacobian_values(self, x: np.ndarray) -> np.ndarray:
        vals = [c.jac_val if c.linear else _values(c.jac(x), c.rows, c.cols, c.k, self.n, 'constraint jac')
                for c in self.cons]
        return np.concatenate(vals) if vals else np.zeros(0)

    def lagrangian_hessian(self, x: np.ndarray, lam: np.ndarray) -> np.ndarray:
        """the Hessian of the Lagrangian f - lambda^T c (dense)"""
        h = self.hess_f(x, *self.args)
        h = np.array(h.toarray() if _is_sparse(h) else h, dtype=float).reshape(self.n, self.n)
        i = 0
        for c in self.cons:
            if not c.linear:
                hc = c.hess(x, lam[i:i + c.k])
                h -= np.asarray(hc.toarray() if _is_sparse(hc) else hc, dtype=float).reshape(self.n, self.n)
            i += c.k
        return h

    # ---- the callbacks given to the Fortran solver, with f2py's calling convention (see `fortran/_sqpopt.pyf`):
    # they write their outputs into the arrays they are given, and return only the status flag. They must not
    # raise (f2py would jump out of the Fortran code, skipping its cleanup): an exception is kept, re-raised by
    # `run`, and the solver is asked to stop (`status < 0`)

    def fc(self, status, x, f, c):
        try:
            xx = np.array(x)
            f[0] = self.objective(xx)
            c[:] = self.constraints(xx)
        except BaseException as e:  # noqa: BLE001
            self.error = e
            status = -1
        return (status,)

    def gjac(self, status, accuracy, x, g, jac_val):
        # (`accuracy` is ignored: the user's derivatives are exact)
        try:
            xx = np.array(x)
            g[:] = self.gradient(xx)
            jac_val[:] = self.jacobian_values(xx)
        except BaseException as e:  # noqa: BLE001
            self.error = e
            status = -1
        return (status,)

    def hess(self, status, x, lam, hess_val):
        try:
            h = self.lagrangian_hessian(np.array(x), np.array(lam))
            hess_val[:] = h[self.hess_rows, self.hess_cols]
        except BaseException as e:  # noqa: BLE001
            self.error = e
            status = -1
        return (status,)

    def report(self, stop, it, x, f, c, lam):
        try:
            xx = np.array(x)
            if self.callback_new_style:
                r = OptimizeResult(x=xx, fun=float(f), nit=int(it), constr=self.split(np.array(c)),
                                   v=self.split(np.array(lam)))
                result = self.callback(intermediate_result=r)
            else:
                result = self.callback(xx)
            if result is True:
                stop = 1
        except StopIteration:
            stop = 1
        except BaseException as e:  # noqa: BLE001
            self.error = e
            stop = 1
        return (stop,)

    def split(self, v: np.ndarray) -> list[np.ndarray]:
        """a vector of all the constraint rows, split by constraint object"""
        out, i = [], 0
        for c in self.cons:
            out.append(v[i:i + c.k].copy())
            i += c.k
        return out

    # ---- the solve

    def multipliers(self, lambda0) -> np.ndarray:
        """`lambda0` (one array per constraint object, as `OptimizeResult.v`, or one of size m) as one array"""
        if isinstance(lambda0, (list, tuple)) and len(lambda0) == len(self.cons) and \
                all(np.size(v) == c.k for v, c in zip(lambda0, self.cons)):
            flat = np.concatenate([np.atleast_1d(np.asarray(v, dtype=float)).ravel() for v in lambda0]) \
                if self.cons else np.zeros(0)
        else:
            flat = np.atleast_1d(np.asarray(lambda0, dtype=float)).ravel()
        if flat.size != self.m:
            raise ValueError(f'lambda0 must have one multiplier per constraint row ({self.m})')
        return flat

    def run(self, options: dict[str, Any], lambda0=None, output_file=None, max_step=None) -> OptimizeResult:
        ext = _native()
        n, m = self.n, self.m
        offsets = np.cumsum([0] + [c.k for c in self.cons])
        irow = np.concatenate([c.rows + o for c, o in zip(self.cons, offsets)]) if self.cons else np.zeros(0)
        icol = np.concatenate([c.cols for c in self.cons]) if self.cons else np.zeros(0)
        use_hess = self.hess_f is not None
        self.hess_rows, self.hess_cols = np.tril_indices(n) if use_hess else (np.zeros(0, int), np.zeros(0, int))
        refs = list(schema.OPTIONS)
        opt_id = np.array([refs.index(r) + 1 for r in options], dtype=np.int32)
        opt_val = np.array([float(v) for v in options.values()], dtype=float)
        c_lb = np.concatenate([c.lb for c in self.cons]) if self.cons else np.zeros(0)
        c_ub = np.concatenate([c.ub for c in self.cons]) if self.cons else np.zeros(0)
        lam0 = np.zeros(m) if lambda0 is None else self.multipliers(lambda0)
        if max_step is None:
            step = np.zeros(n)
        else:
            step = np.broadcast_to(np.asarray(max_step, dtype=float), (n,)).copy()   # (a scalar is for every variable)
            step[~np.isfinite(step) & (step > 0)] = INFINITY

        def pad(a, dtype=float):
            # (a contiguous array with at least one element: see `fortran/_sqpopt.pyf`)
            a = np.asarray(a, dtype=dtype).ravel()
            return np.ascontiguousarray(a if a.size else np.zeros(1, dtype=dtype))

        x, lam, z, c = pad(np.zeros(n)), pad(np.zeros(m)), pad(np.zeros(n)), pad(np.zeros(m))
        iinfo, rinfo, message = ext.sqpopt_py_solve(
            self.fc, self.gjac, self.hess, self.report, int(use_hess), int(self.callback is not None),
            n, m, irow.size, self.hess_rows.size, opt_id.size,
            pad(self.x0), pad(self.x_lb), pad(self.x_ub), pad(c_lb), pad(c_ub),
            pad(irow + 1, np.int32), pad(icol + 1, np.int32),
            pad(self.hess_rows + 1, np.int32), pad(self.hess_cols + 1, np.int32),
            pad(opt_id, np.int32), pad(opt_val), pad(lam0), int(lambda0 is not None),
            pad(step), int(max_step is not None),
            '' if output_file is None else str(output_file), x, lam, z, c)
        x, lam, c = x[:n], lam[:m], c[:m]
        if isinstance(message, bytes):
            message = message.decode(errors='replace')
        if self.error is not None:
            raise self.error
        istat = int(iinfo[0])
        return OptimizeResult(
            x=x, fun=float(rinfo[0]), success=istat <= 2, status=istat, message=str(message).strip(),
            nit=int(iinfo[1]), nfev=self.nfev, njev=int(iinfo[3]), nhev=int(iinfo[4]),
            constr=self.split(c), v=self.split(lam), z=z,
            constr_violation=float(rinfo[2]), kkt_error=float(rinfo[1]), stationarity_error=float(rinfo[3]),
            execution_time=float(rinfo[4]), n_qp_iterations=int(iinfo[5]),
            derivative_switch_iteration=int(iinfo[6]), n_factorizations=int(iinfo[13]),
            n_qp_solves=int(iinfo[14]), n_direct_qp=int(iinfo[15]), n_unconstrained_qp=int(iinfo[16]),
            time_factorization=float(rinfo[7]))


def minimize(fun: Callable, x0, args=(), jac=None, hess=None, bounds=None, constraints: Sequence | Any = (),
             tol: float | None = None, callback: Callable | None = None,
             options: dict | None = None, lambda0=None, output_file=None, max_step=None) -> OptimizeResult:
    """Minimize a function of several variables subject to bounds and constraints, with SQPOPT.

    The interface follows ``scipy.optimize.minimize``:

    fun : callable
        ``fun(x, *args) -> float`` (or ``-> (float, gradient)`` with ``jac=True``).
    x0 : array_like
        The starting point (moved inside the bounds).
    args : tuple
        Extra arguments of `fun`, `jac`, and `hess`.
    jac : callable or True
        The gradient ``jac(x, *args)``, or ``True`` if `fun` returns it too. Required: there are no finite
        differences (unlike scipy). The constraints' Jacobians are required too.
    hess : callable, optional
        The Hessian of the objective ``hess(x, *args)`` (dense or sparse). If given, the solver uses the exact
        Hessian of the Lagrangian (``options%hessian_mode = exact``), and each `NonlinearConstraint` needs
        its ``hess(x, v)``. Otherwise a quasi-Newton approximation is used.
    bounds : Bounds or sequence of (lb, ub) pairs, optional
        Variable bounds (``None`` or ``+-inf`` for no bound).
    constraints : dict, NonlinearConstraint, LinearConstraint, or a sequence of them
        Constraints ``lb <= c(x) <= ub``. A dict (as for SLSQP) has ``'type'`` (``'eq'``: ``fun(x) = 0``;
        ``'ineq'``: ``fun(x) >= 0``), ``'fun'``, ``'jac'``, and optionally ``'args'``. The scipy
        classes are accepted too.
    tol : float, optional
        The KKT tolerance (``options%ktol``).
    callback : callable, optional
        Called once per major iteration, as ``callback(intermediate_result)`` (an `OptimizeResult` with
        ``x``, ``fun``, ``nit``, ``constr``, ``v``) if its argument is named ``intermediate_result``, else
        ``callback(x)``. Return ``True`` or raise ``StopIteration`` to stop.
    options : dict, optional
        Any SQPOPT setting (see ``sqpopt_options.schema``, or the options dialog), by its Fortran reference
        (``'linesearch%filter%delta'``, or with dots), its name if unambiguous (``'ktol'``, ``'max_iter'``),
        or as the nested dict the options dialog returns. Also ``maxiter`` and ``disp`` (as scipy). Choice
        options may be given by their Fortran constant's name, e.g. ``{'hessian_mode': 'sr1'}``.
    lambda0 : optional
        Starting constraint multipliers (not in scipy): one array per constraint object (e.g. the ``v`` of a
        previous result, to warm-start from it with its ``x``), or one array with a multiplier per
        constraint row. By default they start at zero.
    max_step : float or array, optional
        The largest change of each variable in one major iteration (not in scipy): one value for every
        variable, or an array of ``n`` values (``inf``: no limit for that variable). The limits are bounds
        on the steps the solver computes. Use them for variables with different units or sensitivities.
    output_file : str or path, optional
        Write the printed output (see ``print_level``/``disp``) to this file (replacing it), instead of the
        process's standard output, which e.g. a Jupyter notebook doesn't show (not in scipy).

    Returns
    -------
    OptimizeResult
        With ``x``, ``fun``, ``success`` (converged: ``status`` 0, 1, or 2), ``status`` (SQPOPT's status
        code), ``message``, ``nit``, ``nfev`` (calls of `fun`),
        ``njev`` and ``nhev`` (gradient/Jacobian and Hessian evaluations requested by the solver),
        ``constr`` and ``v`` (the values and multipliers of each constraint object), ``z`` (the
        variable-bound multipliers), ``constr_violation``, ``kkt_error``, ``stationarity_error``,
        ``execution_time``, ``n_qp_iterations``, ``derivative_switch_iteration``, ``n_unconstrained_qp`` (QPs
        solved by the unconstrained step, without a QP solver: see ``qp_solver%unconstrained_step``),
        ``n_factorizations``,
        ``n_qp_solves``, ``n_direct_qp``, and ``time_factorization`` (the last four for the options that use
        sparse factorizations: ``inertia_control``, ``direct_qp``, and ``direct_least_squares``;
        ``time_factorization`` is the part of ``execution_time``, in seconds, spent in them). The multipliers are
        those of the Lagrangian ``f - v^T c - z^T x``: positive at a lower bound, negative at an upper one.
    """
    solve = _Solve(fun, x0, args, jac, hess, bounds, constraints, callback)
    opts = _resolve_options(options)
    if tol is not None:
        opts.setdefault('options%ktol', _option_value(schema.OPTIONS['options%ktol'], tol))
    if hess is not None:
        opts.setdefault('options%hessian_mode', _choice('options%hessian_mode', 'sqpopt_hessian_exact'))
    return solve.run(opts, lambda0=lambda0, output_file=output_file, max_step=max_step)
