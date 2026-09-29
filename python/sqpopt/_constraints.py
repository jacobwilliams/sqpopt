"""Bounds and constraints, with the same names and fields as ``scipy.optimize`` (whose objects are accepted too)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np


@dataclass
class Bounds:
    """variable bounds ``lb <= x <= ub`` (use ``-np.inf``/``np.inf`` for no bound), as ``scipy.optimize.Bounds``"""
    lb: Any = -np.inf
    ub: Any = np.inf


@dataclass
class NonlinearConstraint:
    """constraints ``lb <= fun(x) <= ub``, as ``scipy.optimize.NonlinearConstraint``.

    * ``fun(x)`` returns the constraint values (a scalar or a 1-D array).
    * ``jac(x)`` returns their Jacobian, dense ``(k, n)`` or sparse (any object with ``tocoo()``, e.g. a
      ``scipy.sparse`` matrix). Required (there are no finite differences).
    * ``hess(x, v)`` returns :math:`\\sum_i v_i \\nabla^2 c_i(x)` (only needed with an exact Hessian: see
      ``minimize``'s ``hess``).
    * ``jac_sparsity`` (not in scipy): the Jacobian's sparsity pattern, a ``(k, n)`` array or sparse matrix
      whose nonzeros are the entries that may be nonzero. Without it, the pattern is that of ``jac``'s
      result at ``x0`` if that is sparse, else dense.
    """
    fun: Callable
    lb: Any
    ub: Any
    jac: Callable | None = None
    hess: Callable | None = None
    jac_sparsity: Any = None


@dataclass
class LinearConstraint:
    """constraints ``lb <= A @ x <= ub``, as ``scipy.optimize.LinearConstraint`` (``A`` dense or sparse)"""
    A: Any
    lb: Any = -np.inf
    ub: Any = np.inf
