"""The SQPOPT settings: names, types, defaults, limits, and descriptions.

This module has no Qt dependency. It defines every user-settable option of
SQPOPT (not the problem definition, callbacks, or bounds), grouped by topic
for the dialog, and helpers to build, validate, and convert the settings
dict.

The settings dict mirrors the Fortran objects the options are set on, with
the Fortran names::

    {
        'options':      {'max_iter': 100, 'ktol': 1e-06, ...},
        'hessian':      {'damping': True, ...},
        'qp_solver':    {'max_step': 2.0, ...,
                         'dense_qp':  {'max_iter': 100, ...},
                         'sparse_qp': {'null_space': 1, ...}},
        'linesearch':   {'sigma': 0.1, ...,
                         'merit':  {'penalty': 1.0, ...},
                         'filter': {'gamma_theta': 1e-05, ...},
                         'funnel': {'beta': 0.9999, ...}},
        'trust_region': {'enabled': False, ...},
    }

so ``values['qp_solver']['dense_qp']['max_iter']`` is the Fortran
``qp_solver%dense_qp%max_iter``. Selector options (``hessian_mode``,
``qp_solver_mode``, ...) hold the integer value of the Fortran constant
(e.g. ``sqpopt_qp_auto`` = 0), or, with ``enum_names=True``, the constant's
name as a string.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any, Callable, Iterator

__all__ = [
    'Choice', 'Option', 'Section', 'Topic', 'TOPICS', 'OPTIONS', 'SQPOPT_INFINITY',
    'option', 'default_values', 'get_value', 'set_value', 'merge_values',
    'changed_values', 'flatten', 'unflatten', 'validate', 'to_fortran_values',
    'to_enum_names', 'python_literal', 'fortran_assignments',
]

SQPOPT_INFINITY = 1.0e20  #: ``sqpopt_infinity``: bounds at or beyond this are treated as absent


@dataclass(frozen=True)
class Choice:
    """one value of a selector option"""
    value: int                  #: the Fortran integer value
    label: str                  #: a short label for the dialog
    name: str | None = None     #: the Fortran constant, e.g. ``'sqpopt_qp_auto'`` (``None`` if there is none)


@dataclass(frozen=True)
class Option:
    """one setting"""
    path: tuple[str, ...]       #: the Fortran component path, e.g. ``('qp_solver', 'dense_qp', 'max_iter')``
    kind: str                   #: ``'int'``, ``'float'``, ``'bool'``, or ``'choice'``
    default: Any                #: the Fortran default
    doc: str                    #: description
    choices: tuple[Choice, ...] = ()
    minimum: float | None = None
    maximum: float | None = None
    min_open: bool = False      #: whether the lower limit is excluded (``> minimum``)
    max_open: bool = False      #: whether the upper limit is excluded (``< maximum``)
    special: dict[Any, str] = field(default_factory=dict)  #: values with a special meaning, e.g. ``{0: 'automatic'}``

    @property
    def name(self) -> str:
        """the Fortran name of the option (the last component of its path)"""
        return self.path[-1]

    @property
    def fortran(self) -> str:
        """the full Fortran reference, e.g. ``'qp_solver%dense_qp%max_iter'``"""
        return '%'.join(self.path)

    def check(self, value: Any) -> str | None:
        """the reason `value` is invalid for this option, or ``None`` if it is valid"""
        if self.kind == 'bool':
            return None if isinstance(value, bool) else 'must be True or False'
        if self.kind == 'choice':
            return None if value in [c.value for c in self.choices] else \
                'must be one of ' + ', '.join(str(c.value) for c in self.choices)
        if self.kind == 'int' and (isinstance(value, bool) or not isinstance(value, int)):
            return 'must be an integer'
        if self.kind == 'float':
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                return 'must be a number'
            if value != value or value in (float('inf'), float('-inf')):
                return 'must be finite'
        if self.minimum is not None:
            if (value <= self.minimum) if self.min_open else (value < self.minimum):
                return f"must be {'>' if self.min_open else '>='} {format_number(self.minimum)}"
        if self.maximum is not None:
            if (value >= self.maximum) if self.max_open else (value > self.maximum):
                return f"must be {'<' if self.max_open else '<='} {format_number(self.maximum)}"
        return None


@dataclass(frozen=True)
class Section:
    """a group of related options within a topic"""
    title: str
    options: tuple[Option, ...]
    note: str = ''
    #: returns the reason the section's options are not used with the given settings, or ``None`` if they are
    relevance: Callable[[dict], str | None] | None = None


@dataclass(frozen=True)
class Topic:
    """a page of the dialog"""
    title: str
    summary: str
    sections: tuple[Section, ...]

    def options(self) -> Iterator[Option]:
        for s in self.sections:
            yield from s.options


def format_number(x: float) -> str:
    """a compact representation of a number, e.g. ``1e-08`` -> ``'1e-8'``, ``1e4`` -> ``'1e4'``,
    ``100.0`` -> ``'100'``"""
    if isinstance(x, int) and not isinstance(x, bool):
        return str(x)
    if x != 0 and (abs(x) >= 1e4 or abs(x) < 1e-3):
        mant, exp = f'{x:.12e}'.split('e')
        mant = mant.rstrip('0').rstrip('.')
        return f'{mant}e{int(exp)}'
    return f'{x:.12g}'


# ---------------------------------------------------------------------------------------------------------
# the Fortran constants of the selector options

HESSIAN_MODES = (
    Choice(1, 'BFGS (limited-memory, Powell-damped)', 'sqpopt_hessian_bfgs'),
    Choice(2, 'SR1 (limited-memory symmetric rank-1)', 'sqpopt_hessian_sr1'),
    Choice(3, "exact (the user's sparse Hessian of the Lagrangian)", 'sqpopt_hessian_exact'),
)
QP_MODES = (
    Choice(0, 'automatic (dense if n <= auto_dense_max_n, else sparse)', 'sqpopt_qp_auto'),
    Choice(2, 'dense active-set QP', 'sqpopt_qp_dense'),
    Choice(3, 'sparse reduced-Hessian active-set QP', 'sqpopt_qp_reduced_hessian'),
)
LINESEARCH_MODES = (
    Choice(4, 'filter (Fletcher & Leyffer; Wächter & Biegler)', 'sqpopt_linesearch_filter'),
    Choice(5, 'funnel (Kiessling, Leyffer & Vanaret)', 'sqpopt_linesearch_funnel'),
    Choice(1, 'Armijo backtracking on a merit function', 'sqpopt_linesearch_armijo'),
    Choice(3, "Powell's watchdog on a merit function", 'sqpopt_linesearch_watchdog'),
    Choice(2, 'exact minimization of a merit function', 'sqpopt_linesearch_exact'),
)
MERIT_MODES = (
    Choice(1, 'ℓ1 exact penalty', 'sqpopt_merit_l1'),
    Choice(2, 'augmented Lagrangian', 'sqpopt_merit_augmented_lagrangian'),
)
PENALTY_UPDATES = (
    Choice(1, 'keep above the multipliers (never decreases)', 'sqpopt_penalty_multipliers'),
    Choice(2, "the merit function's own model-based rule", 'sqpopt_penalty_model'),
)
RESTORATION_MODES = (
    Choice(1, 'feasibility restoration phase', 'sqpopt_restoration_phase'),
    Choice(2, 'single Gauss-Newton step', 'sqpopt_restoration_gauss_newton'),
)
DERIVATIVE_ACCURACIES = (
    Choice(1, 'fast, then accurate near the solution', 'sqpopt_derivatives_fast'),
    Choice(2, 'accurate throughout', 'sqpopt_derivatives_accurate'),
)
NULL_SPACE_METHODS = (
    Choice(1, 'sparse LU basis (SQOPT-style)', 'sqpopt_null_space_lu'),
    Choice(2, 'orthogonal projections with LSQR', 'sqpopt_null_space_lsqr'),
)
PRINT_LEVELS = (
    Choice(0, 'none'),
    Choice(1, 'iteration log and summary'),
    Choice(2, 'also more columns (step, QP, multipliers, globalization)'),
    Choice(3, 'also the details of each iteration, and the solution'),
)
FUNNEL_UPDATES = (
    Choice(1, 'max(βτ, κθₖ + (1−κ)θ) if the violation decreased, else βτ'),
    Choice(2, 'κτ + (1−κ)θ'),
)

ALL_CHOICES = (HESSIAN_MODES, QP_MODES, LINESEARCH_MODES, MERIT_MODES, PENALTY_UPDATES,
               RESTORATION_MODES, NULL_SPACE_METHODS, DERIVATIVE_ACCURACIES)


def _o(path: str, kind: str, default: Any, doc: str, **kw) -> Option:
    return Option(tuple(path.split('%')), kind, default, doc, **kw)


def _positive(path: str, default: float, doc: str, **kw) -> Option:
    return _o(path, 'float', default, doc, minimum=0.0, min_open=True, **kw)


def _unit_open(path: str, default: float, doc: str) -> Option:
    return _o(path, 'float', default, doc, minimum=0.0, maximum=1.0, min_open=True, max_open=True)


# ---------------------------------------------------------------------------------------------------------
# relevance of the sections to the current settings

def _choice_name(values: dict, path: str, choices: tuple[Choice, ...]) -> str:
    v = get_value(values, tuple(path.split('%')))
    return next((c.name or str(c.value) for c in choices if c.value == v), str(v))


def _ls_mode(values: dict) -> int:
    return get_value(values, ('options', 'linesearch_mode'))


def _tr(values: dict) -> bool:
    return get_value(values, ('trust_region', 'enabled'))


def _needs_ls_mode(*modes: int, trust_region: bool | None = None):
    def rel(values: dict) -> str | None:
        if trust_region is False and _tr(values):
            return 'not used with the trust region (trust_region%enabled)'
        if _ls_mode(values) not in modes:
            return 'not used with linesearch_mode = ' + _choice_name(values, 'options%linesearch_mode',
                                                                    LINESEARCH_MODES)
        return None
    return rel


def _merit_used(values: dict) -> str | None:
    if _ls_mode(values) in (4, 5):
        return ('not used with linesearch_mode = '
                + _choice_name(values, 'options%linesearch_mode', LINESEARCH_MODES)
                + ' (no merit function)')
    return None


def _tr_used(values: dict) -> str | None:
    return None if _tr(values) else 'not used unless trust_region%enabled'


def _line_search_used(values: dict) -> str | None:
    return 'not used with the trust region (trust_region%enabled)' if _tr(values) else None


def _qp_used(*modes: int):
    def rel(values: dict) -> str | None:
        m = get_value(values, ('options', 'qp_solver_mode'))
        if m in modes or m == 0:
            return None
        return 'not used with qp_solver_mode = ' + _choice_name(values, 'options%qp_solver_mode', QP_MODES)
    return rel


def _quasi_newton_used(values: dict) -> str | None:
    if get_value(values, ('options', 'hessian_mode')) == 3:
        return 'not used with hessian_mode = sqpopt_hessian_exact'
    return None


def _exact_used(values: dict) -> str | None:
    if get_value(values, ('options', 'hessian_mode')) != 3:
        return 'only used with hessian_mode = sqpopt_hessian_exact'
    return None


def _inertia_used(values: dict) -> str | None:
    if get_value(values, ('options', 'hessian_mode')) not in (2, 3):
        return 'only used with hessian_mode = sqpopt_hessian_exact or sqpopt_hessian_sr1'
    return None


def _direct_qp_used(values: dict) -> str | None:
    return None if get_value(values, ('options', 'direct_qp')) else 'not used unless options%direct_qp'


def _phase_used(values: dict) -> str | None:
    if get_value(values, ('options', 'restoration_mode')) != 1:
        return 'only used with restoration_mode = sqpopt_restoration_phase'
    return None


# ---------------------------------------------------------------------------------------------------------
# the options, by topic

TOPICS: tuple[Topic, ...] = (

    Topic('Stopping criteria', 'When the solver stops.', (
        Section('Convergence', (
            _positive('options%ktol', 1e-6,
                      'Tolerance on the KKT optimality test: stationarity of the projected Lagrangian gradient, '
                      'plus the sign and complementarity of the constraint multipliers (scaled up only when the '
                      'average multiplier magnitude exceeds 100, as in IPOPT).'),
            _positive('options%ctol', 1e-8, 'Feasibility tolerance on the constraint violation.'),
            _positive('options%dual_inf_tol', 1.0,
                      'Tolerance on the stationarity residual of the unscaled problem (as in IPOPT). With automatic '
                      'scaling, ktol applies to the scaled problem, so an objective scaled far down (e.g. at a poor '
                      'starting point) would make ktol very loose in the original units; convergence also requires '
                      'this.'),
            _positive('options%acceptable_ktol', 1e-4,
                      'Looser "acceptable" KKT tolerance (as in IPOPT): if the KKT test with the acceptable '
                      'tolerances holds for acceptable_iter consecutive iterations, the solver stops with '
                      'sqpopt_acceptable.'),
            _positive('options%acceptable_ctol', 1e-6, 'Looser "acceptable" feasibility tolerance (see acceptable_ktol).'),
            _o('options%acceptable_iter', 'int', 15,
               'Consecutive acceptable iterations needed to stop with sqpopt_acceptable (0 disables the test).',
               minimum=0, special={0: 'disabled'}),
        )),
        Section('Stalled progress', (
            _o('options%ftol', 'float', 1e-8,
               'Once feasible, also stop (with sqpopt_stalled) if the relative change in both the objective and '
               'the variables from the previous iterate is below ftol and xtol for stall_iter consecutive '
               'iterations. Guards against looping to max_iter on marginal steps.', minimum=0.0),
            _o('options%xtol', 'float', 1e-8, 'See ftol: the relative change in the variables.', minimum=0.0),
            _o('options%stall_iter', 'int', 3,
               'Number of consecutive iterations the ftol/xtol stalled-progress test must hold before stopping.',
               minimum=1),
        )),
        Section('Limits', (
            _o('options%max_iter', 'int', 100, 'Maximum number of major SQP iterations.', minimum=0),
            _o('options%max_evals', 'int', 0,
               'Stop (sqpopt_max_evals_reached) after this many calls of the objective and constraints '
               '(0 = no limit).', minimum=0, special={0: 'no limit'}),
            _o('options%max_time', 'float', 0.0,
               'Stop (sqpopt_time_limit_reached) after this much wall-clock time, in seconds (0 = no limit).',
               minimum=0.0, special={0.0: 'no limit'}),
            _o('options%obj_lower_limit', 'float', -SQPOPT_INFINITY,
               'Stop (sqpopt_unbounded) if the objective falls below this at a feasible point. The default, '
               '-sqpopt_infinity (-1e20), never stops.'),
            _o('options%max_consecutive_failures', 'int', 5,
               'Stop (with sqpopt_line_search_failed or sqpopt_qp_solve_failed) after this many consecutive '
               'major iterations whose QP solve or line search/trust-region step failed.', minimum=1),
        )),
    )),

    Topic('Algorithms', 'The main algorithm choices.', (
        Section('Strategies', (
            _o('options%hessian_mode', 'choice', 1,
               'Hessian approximation: limited-memory BFGS (Powell-damped), limited-memory SR1, or the exact '
               'Hessian of the Lagrangian (requires the hess function and its sparsity pattern).',
               choices=HESSIAN_MODES),
            _o('options%qp_solver_mode', 'choice', 0,
               'QP subproblem algorithm. The dense solver forms O(n²) arrays and suits small-to-moderate problems; '
               'the sparse one never forms a dense array. Automatic picks the dense solver when '
               'n <= qp_solver%auto_dense_max_n.', choices=QP_MODES),
            _o('options%linesearch_mode', 'choice', 4,
               'Globalization. The filter (default) and funnel methods use no merit function; the Armijo, '
               'watchdog, and exact searches use the merit function selected by merit_mode. With the trust '
               'region enabled, this selects its acceptance test instead (filter, funnel, or a merit-function '
               'ratio test).', choices=LINESEARCH_MODES),
            _o('options%merit_mode', 'choice', 1,
               'Merit function for the merit-function line searches and the trust-region ratio test: the '
               'non-smooth ℓ1 exact penalty (as in SLSQP), or a smooth augmented Lagrangian (as in NPSOL/SNOPT).',
               choices=MERIT_MODES),
            _o('options%penalty_update', 'choice', 1,
               "How the merit function's penalty parameter is updated: kept above the multiplier estimates "
               "(never decreasing), or each merit function's own rule (Byrd-Nocedal model reduction for ℓ1; "
               "Gill-Murray-Saunders-Wright for the augmented Lagrangian, with a joint step in the multipliers "
               "and slacks).", choices=PENALTY_UPDATES),
            _o('options%restoration_mode', 'choice', 1,
               'Feasibility restoration when no acceptable step is found at an infeasible point: a restoration '
               'phase (feasibility QPs, as in filter-SQP methods and Uno), or a single Gauss-Newton step on the '
               'violation each time.', choices=RESTORATION_MODES),
        )),
        Section('Derivatives', (
            _o('options%derivative_accuracy', 'choice', 2,
               "The accuracy of the derivatives the solver asks gjac for (its accuracy argument). With fast, "
               "gjac may return cheaper, less accurate derivatives (e.g. forward instead of central differences) "
               "until the solver is near a solution, then it asks for accurate ones for the rest of the solve. "
               "It also switches when a step fails, when progress stalls, and before stopping.",
               choices=DERIVATIVE_ACCURACIES),
            _o('options%derivative_switch_tol', 'float', 1e-5,
               'With fast derivatives: the KKT and feasibility errors (of the scaled problem) below which the '
               'solver switches to accurate ones. Switching earlier costs accurate derivatives on iterations '
               'that do not need them.', minimum=0.0),
        )),
        Section('Degenerate points', (
            _o('options%elastic_multiplier_limit', 'float', 30.0,
               "When a constraint's multiplier diverges (its push |λᵢ|·‖∇cᵢ‖∞ exceeds this × max(1, ‖g‖∞) "
               'at two successive iterates, and keeps growing), a sign that the constraint qualification fails '
               'nearby (e.g. a constraint tangent to a bound), the QP is re-solved with that constraint elastic '
               '(an ℓ1 penalty with a bounded weight), so the iterates can leave instead of creeping toward the '
               'degenerate point. At most 3 times per solve. 0 disables this.',
               minimum=0.0, special={0.0: 'disabled'}),
        )),
        Section('Globalization', (
            _o('trust_region%enabled', 'bool', False,
               'Use trust-region radius management instead of the line search for every major iteration '
               '(see the Trust region page). linesearch_mode then selects the acceptance test.'),
        )),
    )),

    Topic('Scaling & output', 'Problem scaling and printed output.', (
        Section('Scaling', (
            _o('options%scaling', 'bool', True,
               'Gradient-based scaling (as in IPOPT): at the starting point, the objective and each constraint '
               'whose gradient has an element larger than scaling_max_gradient is scaled down so that its '
               'largest element equals it. ktol/ctol then apply to the scaled problem; all results are for the '
               'original one.'),
            _positive('options%scaling_max_gradient', 100.0, 'See scaling.'),
        )),
        Section('Output', (
            _o('options%print_level', 'choice', 0,
               '0 = no output. 1 = the problem and method, one line per iteration (objective, infeasibility, KKT '
               'error, step length, evaluations, and flags for the iteration\'s events), and a summary. 2 = also '
               'the step norm, QP iterations, evaluations per iteration, largest multiplier, unscaled stationarity, '
               'and the globalization\'s and Hessian\'s state. 3 = also the details of every iteration (QP solves, '
               'line-search or trust-region trials, restoration phases, Hessian resets) and the solution '
               '(variables and constraints with their bounds, multipliers, and which are active).',
               choices=PRINT_LEVELS),
            _o('options%output_unit', 'int', 6,
               'Fortran unit the output is written to. The Fortran default is output_unit from '
               'iso_fortran_env (standard output, unit 6 with gfortran).'),
        )),
    )),

    Topic('Hessian', 'The Hessian approximation (hessian_mode is on the Algorithms page).', (
        Section('Quasi-Newton', (
            _o('options%lbfgs_memory', 'int', 0,
               'Number of (s,y) vector pairs retained by the limited-memory Hessian. 0 (automatic) picks '
               'max(10, min(n, 100)) from the number of variables n, or 10 with direct_qp (whose cost grows with '
               'the square of the number of pairs); any positive value is used as given.',
               minimum=0, special={0: 'automatic'}),
            _positive('options%hessian_scale0', 1.0,
                      'The initial Hessian approximation is hessian_scale0 times the identity.'),
            _o('hessian%damping', 'bool', True,
               "Powell's damped BFGS update: when sᵀy < 0.2 sᵀBs, y is blended with Bs so the update keeps B "
               "positive definite while still using the new curvature information. If off, such updates are "
               "skipped instead."),
        ), relevance=_quasi_newton_used),
        Section('Inertia control (a build with MUMPS)', (
            _o('options%inertia_control', 'bool', False,
               'With the exact or the SR1 Hessian, which can be indefinite: find the shift δ of H + δI from the '
               'inertia of the KKT matrix of the QP\'s working set, by a sparse LDLᵀ factorization (MUMPS): the '
               'smallest shift tried that leaves no negative curvature. Without it, the exact Hessian is shifted '
               'tenfold whenever a QP finds negative curvature, and SR1 is not corrected at all. It needs a '
               'library built with MUMPS (the HAS_MUMPS preprocessor directive), and is invalid without it.'),
        ), relevance=_inertia_used),
        Section('Exact Hessian', (
            _positive('hessian%shift_min', 1e-4,
                      'The smallest nonzero shift δ of the inertia correction H + δI, relative to '
                      'max(1, max|Hᵢⱼ|).'),
            _positive('hessian%shift_max', 1e10, 'The largest shift δ, relative to max(1, max|Hᵢⱼ|).'),
        ), relevance=_exact_used),
    )),

    Topic('QP solver', 'The QP subproblem solvers (qp_solver_mode is on the Algorithms page).', (
        Section('Direct method (a build with MUMPS)', (
            _o('options%direct_qp', 'bool', False,
               'First try to solve each QP subproblem directly, by sparse factorizations of the KKT matrix of '
               'its working set, starting from the working set of the previous QP. The active-set QP solver is '
               'only run if that fails. Meant for large problems, where it can be orders of magnitude faster; '
               'with the exact Hessian, use it with inertia_control. With L-BFGS or SR1, the automatic memory '
               '(lbfgs_memory = 0) is 10 pairs, which keeps it cheap. It needs a library built with MUMPS (the HAS_MUMPS preprocessor '
               'directive), and is invalid without it.'),
            _o('options%direct_least_squares', 'bool', False,
               'Compute the Gauss-Newton restoration steps and the second-order corrections by a sparse '
               'factorization instead of the iterative LSQR. It pays on large problems that take such steps and '
               'whose constraints are coupled (a chain of 100,000 circle constraints: 88.7 s with LSQR, 1.1 s '
               'with this). It needs a library built with MUMPS, and is invalid without it.'),
        )),
        Section('Threads', (
            _o('options%factorization_threads', 'int', 1,
               'Number of OpenMP threads the sparse factorizations use (inertia_control, direct_qp, and '
               'direct_least_squares). 1 uses none; 0 leaves it to the OpenMP environment (OMP_NUM_THREADS, or '
               'every core). It needs MUMPS and its BLAS built with OpenMP (conda-forge\'s are). Threads only pay '
               'on large problems whose factors are dense enough (a 3-D grid: 2.3 times faster on 4 threads; '
               'banded problems: no gain), and cost a lot on small ones.',
               minimum=0, special={0: 'OpenMP environment'}),
        )),
        Section('Direct method settings', (
            _o('qp_solver%direct_max_changes', 'int', 10,
               'The direct method gives up, and the active-set QP solver is run, after this many changes of the '
               'working set (each is a factorization) without reaching the QP\'s solution.', minimum=0),
            _positive('qp_solver%direct_tol', 1e-8,
                      'The direct method\'s relative tolerance for a violated row or bound, and for the sign of a '
                      'multiplier.'),
        ), relevance=_direct_qp_used),
        Section('General', (
            _o('qp_solver%auto_dense_max_n', 'int', 200,
               'With qp_solver_mode = automatic, the dense QP solver is used for problems with at most this many '
               'variables, and the sparse one for larger problems.', minimum=0),
            _positive('qp_solver%max_step', 2.0,
                      'Trust-region-style cap on the QP step length ‖p‖₂, applied after every QP solve. The cap '
                      'starts at max(max_step, ‖x₀‖∞), doubles after a capped step that the line search accepts '
                      'in full, and halves back toward max_step after a shortened one.'),
        )),
        Section('Dense QP', (
            _o('qp_solver%dense_qp%max_iter', 'int', 100,
               'Minimum limit on active-set iterations per QP solve (the actual limit is '
               'max(max_iter, 10*(rows+1))).', minimum=1),
            _positive('qp_solver%dense_qp%active_tol', 1e-8,
                      'Relative tolerance for a row being at a bound, and for the sign of a multiplier.'),
            _positive('qp_solver%dense_qp%opt_tol', 1e-10,
                      'Relative tolerance on the reduced-gradient stationarity test.'),
            _positive('qp_solver%dense_qp%feas_tol', 1e-6,
                      'An elastic slack larger than feas_tol*max(1, |initial violation|) at the solution counts '
                      'as a violated linearized constraint.'),
            _positive('qp_solver%dense_qp%elastic_weight', 1e4,
                      'Initial elastic penalty weight, relative to max(1, ‖g‖∞).'),
            _positive('qp_solver%dense_qp%elastic_weight_max', 1e10,
                      'Largest elastic penalty weight tried (same scaling) before the linearization is declared '
                      'inconsistent. Must be >= elastic_weight.'),
            _o('qp_solver%dense_qp%warm_start', 'bool', True,
               "Start each QP from the previous QP's final working set (within one solve)."),
        ), relevance=_qp_used(2)),
        Section('Sparse QP', (
            _o('qp_solver%sparse_qp%null_space', 'choice', 1,
               'How the null space of the working set is handled: a sparse LU basis partition (SQOPT-style, with '
               'factor updates), or orthogonal projections with LSQR (much slower at scale; kept for '
               'comparison).', choices=NULL_SPACE_METHODS),
            _o('qp_solver%sparse_qp%max_iter', 'int', 100,
               'Minimum limit on active-set iterations per QP solve (the actual limit is '
               'max(max_iter, 10*(rows+1))).', minimum=1),
            _o('qp_solver%sparse_qp%max_pcg_iter', 'int', 0,
               'Maximum CG iterations per active-set face (<= 0 means twice the number of unknowns).',
               special={0: 'automatic'}),
            _o('qp_solver%sparse_qp%dense_max_ns', 'int', 50,
               'With the LU null space, a face with at most this many superbasics is solved exactly with a '
               'dense Cholesky factorization of the reduced Hessian instead of by preconditioned CG (0 means '
               'always CG).', minimum=0, special={0: 'always CG'}),
            _positive('qp_solver%sparse_qp%active_tol', 1e-8,
                      'Relative tolerance for a row being at a bound, and for the sign of a multiplier.'),
            _positive('qp_solver%sparse_qp%opt_tol', 1e-10,
                      'Relative tolerance on the projected-gradient stationarity test.'),
            _positive('qp_solver%sparse_qp%pcg_rtol', 1e-10,
                      'Projected CG stops once the projected residual has been reduced by this factor.'),
            _positive('qp_solver%sparse_qp%feas_tol', 1e-6, 'As for the dense solver.'),
            _positive('qp_solver%sparse_qp%elastic_weight', 1e4, 'As for the dense solver.'),
            _positive('qp_solver%sparse_qp%elastic_weight_max', 1e8,
                      "As for the dense solver (lower, since the iterative projections' accuracy is relative to "
                      "the weight). Must be >= elastic_weight."),
            _o('qp_solver%sparse_qp%warm_start', 'bool', True, 'As for the dense solver.'),
            _o('qp_solver%sparse_qp%lsqr_atol', 'float', 0.0,
               "LSQR relative error tolerance in A (0 = LSQR's own machine-precision default).", minimum=0.0,
               special={0.0: 'LSQR default'}),
            _o('qp_solver%sparse_qp%lsqr_btol', 'float', 0.0,
               "LSQR relative error tolerance in b (0 = LSQR's own default).", minimum=0.0,
               special={0.0: 'LSQR default'}),
            _o('qp_solver%sparse_qp%lsqr_conlim', 'float', 0.0,
               "LSQR upper limit on cond(Abar) (0 = LSQR's own default).", minimum=0.0,
               special={0.0: 'LSQR default'}),
            _o('qp_solver%sparse_qp%lsqr_itnlim', 'int', 0,
               'LSQR maximum iterations per solve (<= 0 means 2*(rows+columns)+10).',
               special={0: 'automatic'}),
        ), relevance=_qp_used(3)),
    )),

    Topic('Line search', 'The step-length search (linesearch_mode is on the Algorithms page).', (
        Section('General', (
            _positive('linesearch%major_step_limit', 2.0,
                      "Caps the initial trial step length, in all modes, so that no variable changes by more than "
                      "this factor relative to max(1, |xⱼ|) (SNOPT's \"Major step limit\")."),
        ), relevance=_line_search_used),
        Section('Backtracking', (
            _unit_open('linesearch%sigma', 0.1, 'Armijo sufficient-decrease parameter (armijo, watchdog).'),
            _unit_open('linesearch%backtrack', 0.5,
                       'Step-length reduction factor at each backtracking step (when interpolate is off, or the '
                       'interpolation has no minimizer).'),
            _o('linesearch%interpolate', 'bool', True,
               'Choose each backtracking step length by safeguarded quadratic interpolation, within '
               '[0.1α, 0.5α], as in NLPQLP (armijo, watchdog, filter, funnel).'),
            _o('linesearch%nonmonotone_len', 'int', 0,
               'If > 0, a failed search is retried non-monotonically, as in NLPQLP: against the worst merit '
               'value (filter: the worst violation and objective) of this many recent iterates (armijo, filter).',
               minimum=0, special={0: 'off'}),
            _positive('linesearch%alpha_min', 1e-10,
                      'Minimum step length: if no acceptable step is found before α would drop below this, the '
                      'search fails and no step is taken.'),
            _o('linesearch%max_ls_iter', 'int', 40, 'Maximum number of trial step lengths per search.', minimum=1),
            _positive('linesearch%tol', 1e-4, 'Desired tolerance on the minimizer (exact line search).'),
        ), relevance=_line_search_used),
        Section('Watchdog', (
            _o('linesearch%watchdog_relaxed_len', 'int', 2,
               'Number of relaxed steps tolerated before requiring a new best point.', minimum=0),
            _o('linesearch%watchdog_cooldown_len', 'int', 10,
               'Number of iterations relaxed acceptance is disabled for after a backtrack.', minimum=0),
        ), relevance=_needs_ls_mode(3, trust_region=False)),
    )),

    Topic('Merit function', 'The merit function of the Armijo, watchdog, and exact searches, and of the '
          "trust region's ratio test (merit_mode and penalty_update are on the Algorithms page).", (
        Section('Penalty parameter', (
            _o('linesearch%merit%penalty', 'float', 1.0,
               'Initial penalty parameter of the merit function (μ for ℓ1, ρ for the augmented Lagrangian).',
               minimum=0.0),
            _unit_open('linesearch%merit%penalty_rho',
                       0.1, 'penalty_update = model with the ℓ1 merit function: the fraction of the linearized '
                       'violation reduction the penalty must credit (Nocedal & Wright eq. 18.36).'),
        ), relevance=_merit_used),
    )),

    Topic('Filter', 'The filter acceptance test (linesearch_mode = filter; also the trust region with the filter). '
          'Defaults from Wächter & Biegler.', (
        Section('Filter', (
            _unit_open('linesearch%filter%gamma_theta', 1e-5,
                       'Margin γθ: a step must reduce the violation θ by this fraction...'),
            _positive('linesearch%filter%gamma_phi', 1e-5, '...or reduce the objective by γφ·θ.'),
            _positive('linesearch%filter%delta', 1.0,
                      'Switching condition α(−gᵀp)^sφ > δ·θ^sθ: when it holds (and θ <= θmin), the step must '
                      'satisfy an Armijo condition on the objective ("f-type" step).'),
            _o('linesearch%filter%s_theta', 'float', 1.1, 'Switching condition exponent sθ.',
               minimum=1.0, min_open=True),
            _o('linesearch%filter%s_phi', 'float', 2.3, 'Switching condition exponent sφ.',
               minimum=1.0, min_open=True),
            _o('linesearch%filter%eta_phi', 'float', 1e-4, 'Armijo constant for f-type steps.',
               minimum=0.0, maximum=0.5, min_open=True, max_open=True),
            _positive('linesearch%filter%theta_max_fact', 1e4,
                      'θmax as a multiple of max(1, θ(x₀)): no point with a larger violation is accepted. Must '
                      'be > theta_min_fact.'),
            _positive('linesearch%filter%theta_min_fact', 1e-4,
                      'θmin as a multiple of max(1, θ(x₀)): below it, f-type steps are allowed.'),
            _o('linesearch%filter%gamma_alpha', 'float', 0.05,
               'Safety factor in the minimum step length below which the search gives up (and a restoration '
               'step is taken).', minimum=0.0, maximum=1.0, min_open=True),
        ), relevance=_needs_ls_mode(4)),
    )),

    Topic('Funnel', 'The funnel acceptance test (linesearch_mode = funnel; also the trust region with the '
          'funnel). Defaults from the Uno solver.', (
        Section('Funnel', (
            _positive('linesearch%funnel%width_min', 1.0,
                      'Initial funnel width τ₀ = max(width_min, width_fact·θ(x₀)).'),
            _o('linesearch%funnel%width_fact', 'float', 1.5, 'See width_min.', minimum=1.0),
            _unit_open('linesearch%funnel%beta', 0.9999,
                       'An h-type step (the switching condition does not hold) must reach θ <= βτ.'),
            _unit_open('linesearch%funnel%kappa', 0.5,
                       'After an h-type step, the new width is (at least) the convex combination '
                       'κθₖ + (1−κ)θ (see update). Before a restoration step, τ = κτ + (1−κ)θₖ.'),
            _o('linesearch%funnel%update', 'choice', 1,
               'The funnel width update after an h-type step.', choices=FUNNEL_UPDATES),
            _positive('linesearch%funnel%delta', 0.999,
                      'Switching condition constant δ: an f-type step needs α(−gᵀp) > δ·θₖ^sθ, and must then '
                      'satisfy an Armijo condition on the objective.'),
            _o('linesearch%funnel%s_theta', 'float', 2.0, 'Switching condition exponent sθ.',
               minimum=1.0, min_open=True),
            _o('linesearch%funnel%eta', 'float', 1e-4, 'Armijo constant for f-type steps.',
               minimum=0.0, maximum=0.5, min_open=True, max_open=True),
            _o('linesearch%funnel%require_current', 'bool', False,
               'Also require each trial point to improve on the current point: θ < βθₖ or φ <= φₖ − γθ.'),
            _positive('linesearch%funnel%gamma', 1e-3, 'γ in require_current.'),
        ), relevance=_needs_ls_mode(5)),
    )),

    Topic('Trust region', 'Trust-region globalization, an alternative to the line search (enable it on the '
          'Algorithms page). Its steps are judged by the line search\'s acceptance test selected by '
          'linesearch_mode, so the Filter, Funnel, and Merit function settings apply to it too.', (
        Section('Radius', (
            _positive('trust_region%radius0', 1.0, 'Initial trust-region radius.'),
            _positive('trust_region%radius_min', 1e-8,
                      "Below this, a major iteration's retries give up (no step is taken)."),
            _positive('trust_region%radius_max', 1e3, 'Ceiling on the radius.'),
            _unit_open('trust_region%shrink_factor', 0.5, 'radius *= shrink_factor on a rejected step.'),
            _o('trust_region%expand_factor', 'float', 2.0,
               'radius *= expand_factor on an accepted step that used the full radius.', minimum=1.0),
            _o('trust_region%max_retries', 'int', 20,
               'Maximum QP re-solves (with a shrinking radius) per major iteration.', minimum=1),
        ), note='Must satisfy 0 < radius_min <= radius0 <= radius_max.', relevance=_tr_used),
        Section('Ratio test', (
            _unit_open('trust_region%eta1', 0.1,
                       'Ratio threshold to accept a step (merit-function ratio test only, i.e. linesearch_mode '
                       'armijo, watchdog, or exact).'),
            _unit_open('trust_region%eta2', 0.75,
                       'Ratio threshold to also grow the radius (merit-function ratio test only). Must be >= eta1.'),
        ), relevance=_tr_used),
    )),

    Topic('Restoration', 'The feasibility restoration phase (restoration_mode is on the Algorithms page).', (
        Section('Restoration phase', (
            _unit_open('options%restoration_exit_factor', 0.9,
                       'A restoration phase ends once the ℓ1 violation is below this factor times its value where '
                       'the phase started and the point is acceptable to the filter or funnel.'),
            _o('options%restoration_max_iter', 'int', 50, 'Maximum number of iterations of a restoration phase.',
               minimum=1),
        ), relevance=_phase_used),
    )),
)

#: every option, by its Fortran reference (e.g. ``'qp_solver%dense_qp%max_iter'``)
OPTIONS: dict[str, Option] = {o.fortran: o for t in TOPICS for o in t.options()}


def option(ref: str | tuple[str, ...]) -> Option:
    """the option with Fortran reference `ref` (``'a%b%c'`` or ``('a', 'b', 'c')``)"""
    return OPTIONS[ref if isinstance(ref, str) else '%'.join(ref)]


# ---------------------------------------------------------------------------------------------------------
# the settings dict

def get_value(values: dict, path: tuple[str, ...]) -> Any:
    d = values
    for k in path:
        d = d[k]
    return d


def set_value(values: dict, path: tuple[str, ...], value: Any) -> None:
    d = values
    for k in path[:-1]:
        d = d.setdefault(k, {})
    d[path[-1]] = value


def default_values() -> dict:
    """the settings dict with every option at its Fortran default"""
    values: dict = {}
    for o in OPTIONS.values():
        set_value(values, o.path, copy.copy(o.default))
    return values


def _to_value(o: Option, v: Any) -> Any:
    """`v` as the option's Python type (accepting a Fortran constant name for a selector option)"""
    if o.kind == 'choice' and isinstance(v, str):
        for c in o.choices:
            if c.name == v:
                return c.value
        raise ValueError(f'{o.fortran}: unknown value {v!r}')
    if o.kind == 'float' and isinstance(v, int) and not isinstance(v, bool):
        return float(v)
    return v


def merge_values(values: dict | None, base: dict | None = None) -> dict:
    """a complete settings dict: `base` (default: the defaults) updated with the options present in `values`.

    `values` may be nested (as returned by the dialog), possibly with only some options, or flat with Fortran
    references as keys (``{'qp_solver%max_step': 5.0}``). Selector options may be given by value or by the
    name of their Fortran constant. Unknown keys raise ``KeyError``."""
    out = copy.deepcopy(base) if base is not None else default_values()
    if not values:
        return out
    flat = values if all('%' in k for k in values) and not any(isinstance(v, dict) for v in values.values()) \
        else flatten(values, only_known=False)
    for ref, v in flat.items():
        if ref not in OPTIONS:
            raise KeyError(f'unknown SQPOPT option: {ref}')
        o = OPTIONS[ref]
        set_value(out, o.path, _to_value(o, v))
    return out


def flatten(values: dict, only_known: bool = True) -> dict[str, Any]:
    """the settings as a flat dict keyed by Fortran reference, e.g. ``{'qp_solver%dense_qp%max_iter': 100}``"""
    flat: dict[str, Any] = {}

    def walk(d: dict, prefix: tuple[str, ...]):
        for k, v in d.items():
            if isinstance(v, dict):
                walk(v, prefix + (k,))
            else:
                flat['%'.join(prefix + (k,))] = v
    walk(values, ())
    if only_known:
        flat = {k: v for k, v in flat.items() if k in OPTIONS}
    return flat


def unflatten(flat: dict[str, Any]) -> dict:
    """the nested settings dict from a flat one keyed by Fortran reference"""
    values: dict = {}
    for ref, v in flat.items():
        set_value(values, tuple(ref.split('%')), v)
    return values


def changed_values(values: dict) -> dict:
    """only the options of `values` that differ from their defaults (nested)"""
    return unflatten({ref: v for ref, v in flatten(values).items() if v != OPTIONS[ref].default})


def to_enum_names(values: dict) -> dict:
    """`values` with each selector option given by the name of its Fortran constant (where it has one)"""
    flat = flatten(values)
    for ref, v in flat.items():
        o = OPTIONS[ref]
        if o.kind == 'choice':
            flat[ref] = next((c.name for c in o.choices if c.value == v and c.name), v)
    return unflatten(flat)


def to_fortran_values(values: dict) -> dict:
    """`values` with selector options given as their Fortran integer values (the inverse of `to_enum_names`)"""
    return unflatten({ref: _to_value(OPTIONS[ref], v) for ref, v in flatten(values).items()})


def validate(values: dict) -> dict[str, str]:
    """the problems with the settings, as ``{fortran reference: message}`` (empty if they are valid).

    The same checks as the Fortran ``validate_options``, except that ``output_unit`` isn't checked to be open,
    and the exact-Hessian requirement of a user Hessian function can't be checked here."""
    errors: dict[str, str] = {}
    flat = flatten(values)
    for ref, o in OPTIONS.items():
        if ref not in flat:
            errors[ref] = 'missing'
            continue
        msg = o.check(flat[ref])
        if msg:
            errors[ref] = msg
    if errors:
        return errors

    def v(ref: str) -> Any:
        return flat[ref]

    def rule(ok: bool, ref: str, msg: str):
        if not ok and ref not in errors:
            errors[ref] = msg

    for qp in ('dense_qp', 'sparse_qp'):
        rule(v(f'qp_solver%{qp}%elastic_weight_max') >= v(f'qp_solver%{qp}%elastic_weight'),
             f'qp_solver%{qp}%elastic_weight_max', 'must be >= elastic_weight')
    rule(v('linesearch%filter%theta_max_fact') > v('linesearch%filter%theta_min_fact'),
         'linesearch%filter%theta_max_fact', 'must be > theta_min_fact')
    if v('trust_region%enabled'):
        rule(v('trust_region%radius_min') <= v('trust_region%radius0'), 'trust_region%radius0',
             'must be >= radius_min')
        rule(v('trust_region%radius0') <= v('trust_region%radius_max'), 'trust_region%radius_max',
             'must be >= radius0')
        rule(v('trust_region%eta1') <= v('trust_region%eta2'), 'trust_region%eta2', 'must be >= eta1')
    return errors


def python_literal(values: dict) -> str:
    """the settings dict as Python source (for copying into a script)"""
    def fmt(v: Any, indent: int) -> str:
        if isinstance(v, dict):
            pad = ' ' * (indent + 4)
            items = [f'{pad}{k!r}: {fmt(x, indent + 4)},' for k, x in v.items()]
            return '{\n' + '\n'.join(items) + '\n' + ' ' * indent + '}'
        if isinstance(v, float):
            return repr(v)
        return repr(v)
    return fmt(values, 0)


def fortran_assignments(values: dict, only_changed: bool = True) -> str:
    """the settings as Fortran assignment statements, e.g. ``qp_solver%max_step = 5.0_wp``, with selector
    options given by their Fortran constants"""
    lines = []
    for ref, v in flatten(values).items():
        o = OPTIONS[ref]
        if only_changed and v == o.default:
            continue
        if o.kind == 'bool':
            s = '.true.' if v else '.false.'
        elif o.kind == 'choice':
            s = next((c.name for c in o.choices if c.value == v and c.name), str(v))
        elif o.kind == 'int':
            s = str(v)
        elif v == -SQPOPT_INFINITY:
            s = '-sqpopt_infinity'
        elif v == SQPOPT_INFINITY:
            s = 'sqpopt_infinity'
        else:
            s = format_number(float(v))
            if 'e' not in s and '.' not in s:
                s += '.0'
            s += '_wp'
        lines.append(f'{ref} = {s}')
    return '\n'.join(lines)
