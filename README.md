
![Modern Fortran SQP OPTimizer](media/logo-small.png)

Modern Fortran **SQP** **OPT**imizer. A modular and extensible framework for solving large-scale nonlinear optimization problems. A work in progress.

### Overview

SQPOPT is a modern Fortran implementation of a Sequential Quadratic Programming (SQP) optimizer, for problems of the form:

```
minimize    f(x)
subject to  c_lb <= c(x) <= c_ub     (nonlinear, possibly two-sided, constraints)
            x_lb <= x    <= x_ub     (variable bounds)
```

where equality constraints are expressed as `c_lb(i) = c_ub(i)`. `f`, `c`,
and their first derivatives are supplied by the user as callbacks; the
constraint Jacobian is stored in sparse coordinate (COO) format and the
Hessian of the Lagrangian is approximated by a matrix-free limited-memory
quasi-Newton operator -- no dense `n x n` or `m x n` array is ever formed,
so the library scales to large, sparse problems.

#### Features include:
- SQP method: minimizes a nonlinear objective function subject to nonlinear equality and inequality constraints, and bounds.
- Modern Fortran implementation
- Modular architecture -- problem definition, options, Hessian approximation, QP subproblem solver, line search, and convergence checking are each their own module, so alternative algorithms can be developed and swapped in independently
- Open-source and actively maintained
- Sparse matrix support (COO constraint Jacobian; matrix-free limited-memory Hessian)
- Easy integration with existing Fortran projects (uses the FPM build system)
- Selectable real kinds (single, double, quadruple)

### Basic usage

```fortran
use sqpopt_module,         only: sqpopt_type
use sqpopt_problem_module, only: sqpopt_problem_type
use sqpopt_options_module, only: sqpopt_options_type

type(sqpopt_type)         :: solver
type(sqpopt_problem_type) :: problem
type(sqpopt_options_type) :: options
real(wp) :: x(n), lambda(m)
integer  :: istat

call problem%set_problem_size(n, m_eq, m_ineq)
call problem%set_bounds(x_lb, x_ub, c_lb, c_ub)
call problem%set_jacobian_sparsity(nnz, irow, icol)  ! fixed sparsity pattern
call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)   ! optionally also data=...

call solver%initialize(problem=problem, options=options)
call solver%solve(x0, istat)            ! optionally also lambda0=...
call solver%get_solution(x, lambda)     ! optionally also z (bound multipliers)
print *, solver%status_message()        ! e.g. 'converged successfully'
```

Each user function has the form (here, the objective):

```fortran
subroutine obj(x, f, status, data)
    real(wp), dimension(:), intent(in)    :: x
    real(wp),               intent(out)   :: f
    integer,                intent(inout) :: status  ! 0 on entry
    class(*), optional,     intent(inout) :: data
    f = ...
end subroutine obj
```

(`grad`, `cons`, and `jacv` are the same, with an array output `g(n)`,
`c(m)`, or `jac_val(nnz)`.) The two trailing arguments are:

- `status`: leave it `0` on success. Set it `> 0` if the function can't
  be evaluated at `x` (e.g. a domain error): the solver treats that point
  like one where the function returned NaN, and backs off from it. Set it
  `< 0` to stop the solver (`sqpopt_user_requested_stop`); no user
  function is called again after that.
- `data`: the object passed as `set_functions(..., data=my_data)` (absent
  if none was given), for passing any context to the functions without
  module variables. Access it with `select type`. It is *pointed to*, not
  copied, so the caller's object needs the `target` attribute; updates
  the functions make to it are seen by the caller. It is also passed to
  the `report` callback.

After a solve, `solver%get_results(results)` returns a
`sqpopt_results_type` (from `sqpopt_types_module`) containing:
- the status and message, and the number of iterations;
- the evaluation counts of each user function, and the run time;
- the final `x`, `f`, and `c`;
- the constraint multipliers `lambda` and the variable-bound multipliers
  `z`;
- the KKT and feasibility errors.

All of these are for the original problem, even when automatic scaling
is on. The Lagrangian is \( f - \lambda^T c - z^T x \), so a multiplier is
`>= 0` at a lower bound and `<= 0` at an upper bound.

Before iterating, `solve` validates the problem definition and options
(returning `istat=sqpopt_invalid_input`, with the reason in
`status_message()`, if anything is wrong), and moves `x0` inside the
variable bounds, so the user functions are never evaluated outside them.
If no acceptable step can be found along a search direction, no step is
taken (the point is never moved to one that makes the merit function
worse); the Hessian approximation is reset so the next iteration tries a
different direction. Every call to `solve` starts from the configuration given to `initialize`:
no state (penalty parameter, filter, Hessian, ...) carries over from a
previous solve.

#### Status codes (`istat`, from `sqpopt_types_module`)

| value | meaning |
|---|---|
| `sqpopt_success` (`0`) | the KKT conditions are satisfied to within `ktol`/`ctol` |
| `sqpopt_max_iter_reached` (`1`) | `max_iter` major iterations were performed |
| `sqpopt_infeasible` (`2`) | the constraints are violated at a point that is stationary for the constraint violation: the problem appears to be (locally) infeasible |
| `sqpopt_line_search_failed` (`3`) | `max_consecutive_failures` consecutive iterations failed to find an acceptable step |
| `sqpopt_qp_solve_failed` (`4`) | `max_consecutive_failures` consecutive QP subproblem solves failed |
| `sqpopt_user_requested_stop` (`5`) | the `report` callback, or a user function (`status < 0`), asked the solver to stop |
| `sqpopt_invalid_input` (`6`) | the problem definition or options are invalid (see `status_message()`) |
| `sqpopt_stalled` (`7`) | feasible, but the objective and variables have stopped changing (see `ftol`/`xtol`) before the KKT test was satisfied; usually an acceptable, if less precise, solution |
| `sqpopt_function_error` (`8`) | a problem function returned a non-finite value (NaN or Inf), or `status > 0`, at the current point (at a *trial* point, that just makes the line search/trust region reject the point and back off) |
| `sqpopt_max_evals_reached` (`9`) | `max_evals` objective evaluations were performed |
| `sqpopt_time_limit_reached` (`10`) | the `max_time` limit was reached |
| `sqpopt_unbounded` (`11`) | the objective fell below `obj_lower_limit` at a feasible point |
| `sqpopt_acceptable` (`12`) | the looser `acceptable_ktol`/`acceptable_ctol` tests held for `acceptable_iter` consecutive iterations (as in IPOPT), but the normal ones did not |

See [test/test_basic.f90](test/test_basic.f90), [test/test_hs71.f90](test/test_hs71.f90),
and [test/test_medium.f90](test/test_medium.f90) for complete worked examples.

### Configuration

`solver%initialize(problem=..., options=..., hessian=..., qp_solver=...,
linesearch=..., trust_region=..., report=...)` accepts one instance of
each sub-component, all optional (defaults are used for anything
omitted). Every setting of every component is checked when `solve`
starts; an invalid value gives `sqpopt_invalid_input`, and
`status_message()` says which setting it was. `sqpopt_options_type` covers the most commonly-tuned,
algorithm-*selecting* settings and is copied down into the other
components' `mode`-like fields at the start of every `solve()` call; the
other types (`sqpopt_hessian_type`, `sqpopt_qp_solver_type`,
`sqpopt_linesearch_type`, `sqpopt_trust_region_type`) expose further
algorithm-specific tuning parameters and are configured by constructing
them directly, e.g.:

```fortran
type(sqpopt_qp_solver_type)  :: qp_solver
type(sqpopt_linesearch_type) :: linesearch

qp_solver%max_step             = 5.0_wp
qp_solver%sparse_qp%lsqr_atol  = 5.0e-10_wp
linesearch%major_step_limit    = 1.0_wp

call solver%initialize(problem=problem, options=options, qp_solver=qp_solver, linesearch=linesearch)
```

#### Problem definition (`sqpopt_problem_type`)

| method | description |
|---|---|
| `set_problem_size(n, m_eq, m_ineq)` | number of variables `n`, equality constraints `m_eq`, and inequality constraints `m_ineq` (total `m = m_eq + m_ineq`) |
| `set_bounds(x_lb, x_ub, c_lb, c_ub)` | variable bounds and constraint bounds; use `c_lb(i)==c_ub(i)` for an equality constraint, and any value with magnitude `>= sqpopt_infinity` (`1e20`; `huge(1.0_wp)` is fine) for a one-sided/absent bound |
| `set_jacobian_sparsity(nnz, irow, icol)` | fixed 1-based COO sparsity pattern of the constraint Jacobian |
| `set_hessian_sparsity(nnz, irow, icol)` | fixed 1-based COO sparsity pattern of the Lagrangian Hessian (only needed for `sqpopt_hessian_exact`, which is not yet implemented) |
| `set_functions(f, g, c, jac, hess)` | attach the user-supplied callbacks: `f`/`g` evaluate the objective and its gradient; `c`/`jac` evaluate the constraints and the Jacobian's nonzero values; `hess` (optional) evaluates the Hessian's nonzero values (only used by `sqpopt_hessian_exact`) |

#### Solver options (`sqpopt_options_type`)

| option | default | description |
|---|---|---|
| `max_iter` | `100` | maximum number of major SQP iterations |
| `print_level` | `0` | `0` = no output; `1` = an iteration table (iteration, objective, infeasibility, KKT error, step length, and flags: `R` restoration step, `Q` QP failed, `F` no acceptable step) plus a final summary; `2` = also the penalty parameter, step norm, and QP iterations |
| `output_unit` | standard output | Fortran unit the output is written to |
| `hessian_mode` | `sqpopt_hessian_bfgs` | Hessian approximation strategy (see the [Hessian approximation](#hessian-approximation-sqpopt_hessian_type) table below) |
| `lbfgs_memory` | `10` | number of `(s,y)` vector pairs retained by the limited-memory Hessian |
| `qp_solver_mode` | `sqpopt_qp_auto` | QP subproblem algorithm (see the [QP subproblem solver](#qp-subproblem-solver-sqpopt_qp_solver_type) tables below) |
| `linesearch_mode` | `sqpopt_linesearch_filter` | line search strategy (see the [Line search & merit function](#line-search--merit-function-sqpopt_linesearch_type) tables below) |
| `merit_mode` | `sqpopt_merit_l1` | merit function, for the merit-function line searches (see the [Line search & merit function](#line-search--merit-function-sqpopt_linesearch_type) tables below) |
| `penalty_update` | `sqpopt_penalty_multipliers` | how the merit function's penalty parameter is updated (see the tables below) |
| `max_consecutive_failures` | `5` | stop (with `sqpopt_line_search_failed` or `sqpopt_qp_solve_failed`) after this many consecutive major iterations whose QP solve or line search/trust-region step failed |
| `ftol`, `xtol` | `1e-8` | once feasible, also stop (with `istat=sqpopt_stalled`) if the objective's and the variables' relative change from the previous iterate are both below these tolerances (a safeguard against looping to `max_iter` on marginal steps when the KKT test in `ktol` never quite converges) |
| `ctol` | `1e-8` | feasibility tolerance on the constraint violation |
| `acceptable_ktol`, `acceptable_ctol`, `acceptable_iter` | `1e-4`, `1e-6`, `15` | if the KKT test with these looser tolerances holds for `acceptable_iter` consecutive iterations, stop with `sqpopt_acceptable` (`acceptable_iter=0` disables this) |
| `max_evals` | `0` | stop after this many objective evaluations (`0` = no limit) |
| `max_time` | `0` | stop after this much wall-clock time, in seconds (`0` = no limit) |
| `obj_lower_limit` | `-sqpopt_infinity` | stop (`sqpopt_unbounded`) if the objective falls below this at a feasible point |
| `scaling`, `scaling_max_gradient` | `.true.`, `100` | gradient-based scaling (as in IPOPT): at the starting point, the objective and each constraint whose gradient has an element larger than `scaling_max_gradient` is scaled down so that its largest element equals it. `ktol`/`ctol` then apply to the scaled problem; all results are for the original one |
| `hessian_scale0` | `1` | the initial Hessian approximation is `hessian_scale0` times the identity |
| `ktol` | `1e-6` | tolerance on the KKT optimality test: stationarity of the projected Lagrangian gradient, plus the sign and complementarity of the constraint multipliers (scaled up only when the average multiplier magnitude exceeds 100, as in IPOPT) |

`hessian_mode`/`lbfgs_memory`, `qp_solver_mode`, `linesearch_mode`,
`merit_mode`, and `penalty_update` are copied onto the corresponding
sub-component (`hessian`, `qp_solver%mode`, `linesearch%mode`,
`linesearch%merit_mode`, `linesearch%penalty_update`) at the start
of every `solve()` call, so setting them via `options` and via a directly-
constructed sub-component are equivalent -- everything else in the tables
below (tolerances, step limits, etc.) must be set on the sub-component
directly, since `sqpopt_options_type` has no field for it.

#### Hessian approximation (`sqpopt_hessian_type`)

A limited-memory quasi-Newton approximation (never a dense `n x n`
matrix); `solve()` always (re)initializes it from `options%lbfgs_memory`/
`options%hessian_mode` at the start of every call. Apart from those two
`options` fields, only `damping` is configurable, on a directly-constructed
`sqpopt_hessian_type`:

| option | default | description |
|---|---|---|
| `hessian_mode` (`sqpopt_options_type`) | `sqpopt_hessian_bfgs` | `sqpopt_hessian_bfgs` (limited-memory Powell-damped BFGS), `sqpopt_hessian_sr1` (limited-memory symmetric rank-1), or `sqpopt_hessian_exact` (not yet implemented, falls back to BFGS) |
| `lbfgs_memory` (`sqpopt_options_type`) | `10` | number of `(s,y)` vector pairs retained, independent of the problem size `n` |
| `damping` | `.true.` | Powell's damped BFGS update: when \( s^Ty < 0.2\,s^TBs \), `y` is blended with `Bs` so the update keeps `B` positive definite while still using the new curvature information; if `.false.`, such updates are skipped instead |

#### QP subproblem solver (`sqpopt_qp_solver_type`)

| option | default | description |
|---|---|---|
| `mode` | `sqpopt_qp_auto` | which QP algorithm to use (overwritten from `options%qp_solver_mode` at the start of `solve()`) -- see the mode table below |
| `auto_dense_max_n` | `200` | `mode=sqpopt_qp_auto` uses `sqpopt_qp_dense` for problems with at most this many variables, else `sqpopt_qp_reduced_hessian` |
| `max_step` | `2.0` | initial trust-region-style cap on \|\|p\|\|₂, applied after every QP solve regardless of `mode`; the cap adapts like a trust radius (it doubles after a capped step that the line search accepts in full, and halves back toward `max_step` after a shortened one), so solutions far from the start are still reached in a few iterations |
| `dense_qp` | -- | the dense active-set QP solver's own options (used only when `mode==sqpopt_qp_dense`); see table below |
| `sparse_qp` | -- | the sparse active-set QP solver's own options (used only when `mode==sqpopt_qp_reduced_hessian`); see table below |

**`mode` values (`qp_solver%mode` / `options%qp_solver_mode`):**

| value | description |
|---|---|
| `sqpopt_qp_auto` | (default) `sqpopt_qp_dense` if `n <= auto_dense_max_n`, else `sqpopt_qp_reduced_hessian` |
| `sqpopt_qp_dense` | dense active-set QP (Householder QR null space, Cholesky of the reduced Hessian); forms `O(n^2)`/`O(mn)` dense arrays each call, so best for small-to-moderate problems |
| `sqpopt_qp_reduced_hessian` | sparse active-set QP (SQOPT-style basis partition of the working set with sparse `LUSOL` factors, conjugate gradients on the reduced Hessian, active bounds handled by fixing variables); never forms a dense array, so it scales to larger problems |

Both active-set solvers enforce the linearized constraints and bounds
exactly, and share the same robustness features:

- **Elastic mode.** Each linearized constraint violated at `p=0` gets a
  nonnegative slack with an \( \ell_1 \) penalty (SNOPT's elastic mode),
  which gives a feasible starting point with no separate phase-1 problem.
  If slacks remain positive at the solution, the penalty weight is raised
  (up to `elastic_weight_max`); if they are still positive then, the
  linearized constraints are inconsistent and the QP returns
  `sqpopt_infeasible`, and the major iteration takes a feasibility-
  restoration step instead.
- Rows join the working set only if linearly independent of it (the
  sparse solver picks its initial basis and working set with one
  rank-revealing `LUSOL` factorization, preferring to keep equality
  constraints, then bounds, then inequalities in the working set, so
  duplicated or dependent constraints are handled);
  directions of zero or negative curvature (an indefinite SR1 Hessian, or
  along an elastic slack) are followed to the nearest blocking constraint
  rather than producing a huge Newton step; and all tolerances are
  relative to the problem's scale.
- **Crash and warm starts.** The iterations start from the minimum-norm
  step satisfying an initial working-set guess -- the previous QP's final
  working set (`warm_start`, on by default: near a solution the active set
  settles and the QP finishes in one or two iterations), or else the
  equality constraints and fixed variables -- rather than from `p=0`, so
  only the constraints that guess leaves violated need elastic slacks.
- Both are checked against the KKT conditions on thousands of random
  convex, nonconvex, degenerate, and infeasible QPs (`test/test_qp_fuzz.f90`).

**Dense active-set QP options (`qp_solver%dense_qp`, used when `mode==sqpopt_qp_dense`):**

| option | default | description |
|---|---|---|
| `max_iter` | `100` | minimum limit on active-set iterations per QP solve (the actual limit is `max(max_iter, 10*(rows+1))`) |
| `active_tol` | `1e-8` | relative tolerance for a row being at a bound, and for the sign of a multiplier |
| `opt_tol` | `1e-10` | relative tolerance on the reduced-gradient stationarity test |
| `feas_tol` | `1e-6` | an elastic slack larger than `feas_tol*max(1,\|initial violation\|)` at the solution counts as a violated linearized constraint |
| `elastic_weight` | `1e4` | initial elastic penalty weight, relative to \( \max(1,\lVert g \rVert_\infty) \) |
| `elastic_weight_max` | `1e10` | largest elastic penalty weight tried (same scaling) before the linearization is declared inconsistent |
| `warm_start` | `.true.` | start each QP from the previous QP's final working set (within one `solve`) |

**Sparse active-set QP options (`qp_solver%sparse_qp`, used when `mode==sqpopt_qp_reduced_hessian`):**

| option | default | description |
|---|---|---|
| `null_space` | `sqpopt_null_space_lu` | how the null space of the working set is handled: `sqpopt_null_space_lu` (SQOPT-style) gives every general row a slack variable, so the working set is just the unknowns fixed at a bound, and splits the free unknowns into a nonsingular basis `B` and the superbasic rest `S`: the null space is `Z = [-B⁻¹S; I]`, every projection and multiplier is a direct solve with `B`'s sparse LU factors, and working-set changes update the factors (`lu8rpc`) instead of refactorizing; `sqpopt_null_space_lsqr` uses orthogonal projections, each an iterative `LSQR` least-squares solve, with projected CG (much slower at scale, kept for comparison) |
| `max_iter` | `100` | minimum limit on active-set iterations per QP solve (the actual limit is `max(max_iter, 10*(rows+1))`) |
| `max_pcg_iter` | `0` | maximum CG iterations per active-set face (`<=0` means twice the number of unknowns; CG is also stopped at the dimension of the face) |
| `dense_max_ns` | `50` | with `null_space=sqpopt_null_space_lu`, a face with at most this many superbasics is solved exactly with a dense Cholesky factorization of the reduced Hessian (if positive definite) instead of by (diagonally preconditioned) CG; `0` means always CG |
| `active_tol` | `1e-8` | relative tolerance for a row being at a bound, and for the sign of a multiplier |
| `opt_tol` | `1e-10` | relative tolerance on the projected-gradient stationarity test |
| `pcg_rtol` | `1e-10` | projected CG stops once the projected residual has been reduced by this factor |
| `feas_tol` | `1e-6` | as for the dense solver |
| `elastic_weight` | `1e4` | as for the dense solver |
| `elastic_weight_max` | `1e8` | as for the dense solver (lower, since the iterative projections' accuracy is relative to the weight) |
| `warm_start` | `.true.` | as for the dense solver |
| `lsqr_atol`, `lsqr_btol`, `lsqr_conlim` | `0.0` | `LSQR` relative error tolerances in `A`/`b`, and the upper limit on `cond(Abar)` (`0` means "let `LSQR` use its own machine-precision-based default", which is tighter than usually necessary), for the minimum-norm starting step and (with `null_space=sqpopt_null_space_lsqr`) every projection |
| `lsqr_itnlim` | `0` | `LSQR` maximum iterations per solve (`<=0` means `2*(rows+columns)+10`) |

#### Line search & merit function (`sqpopt_linesearch_type`)

| option | default | description |
|---|---|---|
| `mode` | `sqpopt_linesearch_filter` | line search strategy (overwritten from `options%linesearch_mode` at the start of `solve()`) -- see the mode table below |
| `merit_mode` | `sqpopt_merit_l1` | merit function (overwritten from `options%merit_mode` at the start of `solve()`) -- see the mode table below |
| `penalty` | `1.0` | current penalty parameter used in the merit function (\( \mu \) for `sqpopt_merit_l1`, \( \rho \) for `sqpopt_merit_augmented_lagrangian`) |
| `penalty_update` | `sqpopt_penalty_multipliers` | how the penalty is updated (overwritten from `options%penalty_update` at the start of `solve()`) -- see the table below |
| `penalty_rho` | `0.1` | `sqpopt_penalty_model` with `sqpopt_merit_l1`: the fraction of the linearized violation reduction the penalty must credit |
| `major_step_limit` | `2.0` | caps the *initial* trial step length (before any backtracking), used by all three modes, so that no variable changes by more than this factor relative to \( \max(1,\lvert x_j\rvert) \) (SNOPT's "Major step limit" option); guards against divergence from a QP step that is technically feasible but unreasonably large |
| `sigma` | `0.1` | Armijo sufficient-decrease parameter, \( 0<\sigma<1 \) (`sqpopt_linesearch_armijo`/`sqpopt_linesearch_watchdog` modes) |
| `backtrack` | `0.5` | step-length reduction factor at each backtracking step (when `interpolate` is off, or the interpolation has no minimizer) |
| `interpolate` | `.true.` | choose each backtracking step length by safeguarded quadratic interpolation of the merit function (or, in the filter search, of the violation if the trial made it worse, else the objective), within \( [0.1\alpha, 0.5\alpha] \), as in NLPQLP (`armijo`/`watchdog`/`filter` modes) |
| `nonmonotone_len` | `0` | if `> 0`, a failed search is retried non-monotonically, as in NLPQLP: against the worst merit value (filter search: the worst violation and objective) of the last `nonmonotone_len` iterates instead of the current one (`armijo`/`filter` modes). Helps the merit-function searches on hard problems; not the filter search |
| `alpha_min` | `1e-10` | minimum step length: if no acceptable step is found before `alpha` would drop below this, the search fails and no step is taken (`armijo`/`watchdog`/`filter` modes) |
| `max_ls_iter` | `40` | maximum number of trial step lengths per search (`armijo`/`watchdog`/`filter` modes) |
| `tol` | `1e-4` | desired tolerance on the minimizer (`sqpopt_linesearch_exact` mode) |
| `watchdog_relaxed_len` | `2` | number of relaxed steps tolerated before requiring a new best point (`sqpopt_linesearch_watchdog` mode) |
| `watchdog_cooldown_len` | `10` | number of iterations relaxed acceptance is disabled for after a backtrack (`sqpopt_linesearch_watchdog` mode) |
| `filter_gamma_theta`, `filter_gamma_phi` | `1e-5`, `1e-5` | filter margins \( \gamma_\theta,\gamma_\varphi \): a step must reduce the violation \( \theta \) by the fraction \( \gamma_\theta \), or the objective by \( \gamma_\varphi\theta \) (`filter` mode, and the trust region with it) |
| `filter_delta`, `filter_s_theta`, `filter_s_phi` | `1`, `1.1`, `2.3` | switching condition \( \alpha(-g^Tp)^{s_\varphi} > \delta\theta^{s_\theta} \): when it holds (and \( \theta\le\theta_{min} \)) the step must satisfy an Armijo condition on the objective |
| `filter_eta_phi` | `1e-4` | Armijo constant for those ("f-type") steps |
| `filter_theta_max_fact`, `filter_theta_min_fact` | `1e4`, `1e-4` | \( \theta_{max} \), \( \theta_{min} \) as multiples of \( \max(1,\theta(x_0)) \): no point with violation above \( \theta_{max} \) is accepted |
| `filter_gamma_alpha` | `0.05` | safety factor in the minimum step length below which the search gives up (and a restoration step is taken) |

**`mode` values (`linesearch%mode` / `options%linesearch_mode`):**

| value | description |
|---|---|
| `sqpopt_linesearch_armijo` | standard backtracking line search with an Armijo-type sufficient-decrease test on the merit function (as used by default in `slsqp`) |
| `sqpopt_linesearch_exact` | (approximately) minimizes the merit function along the search direction using the derivative-free `fmin` routine |
| `sqpopt_linesearch_watchdog` | Powell's watchdog technique: tracks the best point found so far and, for a short window after a genuine improvement, relaxes the sufficient-decrease test (accepting the full step outright) rather than stalling near a curved/simultaneously-active constraint boundary (the Maratos effect); backtracks to the best point and disables relaxed acceptance for `watchdog_cooldown_len` iterations if the window is used up without a new best point |
| `sqpopt_linesearch_filter` | (default) filter line search (Fletcher & Leyffer's filter, with the globally convergent line-search rules of Wächter & Biegler, as in IPOPT): no merit function or penalty parameter; a trial point is judged by its (violation, objective) pair, which must not be dominated by the filter; a switching condition decides whether the step must reduce the objective (Armijo) or may trade it for feasibility, and only the latter steps enlarge the filter; if no acceptable step is found, a feasibility-restoration step is taken. Ignores `merit_mode`/`penalty` |

**`merit_mode` values (`linesearch%merit_mode` / `options%merit_mode`):**

| value | description |
|---|---|
| `sqpopt_merit_l1` | (default) non-smooth \( \ell_1 \) exact penalty function (as in `slsqp`) |
| `sqpopt_merit_augmented_lagrangian` | smooth augmented Lagrangian merit function (Gill, Murray, Saunders & Wright; the merit function used in NPSOL and, in spirit, SNOPT) -- twice continuously differentiable, which helps avoid the Maratos effect |

**`penalty_update` values (`linesearch%penalty_update` / `options%penalty_update`):**

| value | description |
|---|---|
| `sqpopt_penalty_multipliers` | (default) keep the penalty above the multiplier estimates, \( \mu \ge \lVert\lambda\rVert_\infty + 1 \) (Han/Powell, as in `slsqp`); it never decreases, so it can grow very large at degenerate points |
| `sqpopt_penalty_model` | each merit function's principled rule. `sqpopt_merit_l1`: Byrd-Nocedal model reduction -- the penalty only increases as needed for the step's predicted merit reduction to credit a fraction `penalty_rho` of its linearized violation reduction. `sqpopt_merit_augmented_lagrangian`: Gill-Murray-Saunders-Wright -- the line search moves the multipliers (toward the QP's) and the slacks (toward the QP's linearized constraint values) along with the variables, the penalty is set just large enough for the merit's slope to be at most \( -\tfrac12 p^THp \), and it can also decrease (a limited number of times) |

On the Hock-Schittkowski test set (`test/test_hs_suite.f90`, 305 problems):

| line search / merit / penalty | solved | local | failed | `f` evaluations (solved) |
|---|--:|--:|--:|--:|
| **filter** (default, with interpolation) | **274** | 31 | **0** | **10,369** |
| filter, without interpolation | 273 | 32 | 0 | 10,593 |
| Armijo / \( \ell_1 \) / multipliers | 269 | 32 | 4 | 28,555 |
| Armijo / \( \ell_1 \) / model (Byrd-Nocedal) | 271 | 32 | 2 | 31,154 |
| Armijo / augmented Lagrangian / multipliers | 270 | 33 | 2 | 10,960 |
| Armijo / augmented Lagrangian / multipliers, interpolation + non-monotone (10) | 273 | 31 | 1 | 11,442 |
| Armijo / augmented Lagrangian / model (GMSW) | 266 | 34 | 5 | 14,198 |
| watchdog / \( \ell_1 \) / model | 271 | 32 | 2 | 34,286 |

(The merit-function rows are without interpolation, except where noted. With interpolation, Armijo / \( \ell_1 \) / multipliers needs 18,317 `f` evaluations instead of 28,555, but fails on one more problem.)

**Second-order correction.** In every line search, and in the trust
region, when the first (full) trial step
is rejected without reducing the constraint violation -- the signature of
the Maratos effect near curved constraints -- a second-order-corrected
step (a minimum-norm correction for the true constraint values at the
trial point, using the current Jacobian, kept within the bounds) is tried
before backtracking. See `sqpopt_soc_module`.

#### Trust-region globalization (`sqpopt_trust_region_type`)

An opt-in *alternative* to the line search above (disabled by default):
instead of solving the QP once and searching for a step length `alpha`
along the resulting `p`, the QP is re-solved as needed with a shrinking
trust-region radius (enforced by temporarily tightening the variable
bounds passed to whichever `qp_solver_mode` is selected -- no QP solver
changes needed) until a step is accepted or the retries are exhausted.
See `plan/TRUST_REGION_PLAN.md` for the full design.

| option | default | description |
|---|---|---|
| `enabled` | `.false.` | if `.true.`, use trust-region radius management instead of `linesearch%search` for every major iteration |
| `radius0` | `1.0` | initial trust-region radius |
| `radius_min` | `1e-8` | below this, a major iteration's retries give up and accept the last (smallest-radius) trial point anyway |
| `radius_max` | `1e3` | ceiling on the radius |
| `eta1` | `0.1` | ratio threshold to accept a step (merit-ratio acceptance only -- see below) |
| `eta2` | `0.75` | ratio threshold to also grow the radius (merit-ratio acceptance only) |
| `shrink_factor` | `0.5` | `radius *= shrink_factor` on a rejected step |
| `expand_factor` | `2.0` | `radius *= expand_factor` on an accepted step that used the full radius |
| `max_retries` | `20` | maximum QP re-solves (with a shrinking radius) per major iteration |

When enabled, `linesearch%mode` is reinterpreted as *which acceptance
test* to use, not which line search to run (there is no `alpha` to
search): `sqpopt_linesearch_filter` uses the same filter acceptance test
as the filter line search (with the quadratic model's predicted decrease
in the switching condition) -- a trust-region filter-SQP method in the
spirit of Fletcher & Leyffer -- while `armijo`/`exact`/`watchdog` all collapse to
the same classical trust-region-SQP ratio test on the merit function
selected by `merit_mode`.

See [plan/PLAN.md](plan/PLAN.md) for the full architecture write-up, algorithm
details, and backlog of future work.

### Developing

Use the `pixi` environment and the Fortran Package Manager (FPM):

```
pixi shell
fpm build --profile release
fpm test --profile release
```

A scalable benchmark (function evaluations and run time on a nonlinear
optimal-control problem and a constrained chained-Rosenbrock problem, at
sizes that exercise both QP solvers) is in `example/benchmark.f90`:

```
fpm run --example benchmark --profile release
```

### Dependencies of this package

This package depends on the following external libraries (which will be automatically fetched and built by FPM):

* [LSQR](https://github.com/jacobwilliams/LSQR) -- iterative solver for sparse linear systems and least-squares problems
* [lusol](https://github.com/jacobwilliams/lusol) -- sparse LU factorization library (the sparse QP's basis factors and updates, and its rank-revealing basis choice)
* [fmin](https://github.com/jacobwilliams/fmin.git) -- derivative-free minimization routine used for exact line search


### Other Fortran SQP Solvers

 * [SNOPT](https://ccom.ucsd.edu/~optimizers/solvers/snopt/) -- Large-scale SQP solver developed by Philip Gill, Walter Murray, and Michael Saunders. A commercial product.
 * [VF13AD](https://www.hsl.rl.ac.uk/archive/) -- Classic SQP method from the HSL Archive.
 * [SLSQP](https://github.com/jacobwilliams/slsqp) -- Originally by Dieter Kraft, one of the optimization methods in SciPy.
 * [PSQP](https://github.com/jacobwilliams/psqp) -- Another SQP code, originally by Ladislav Luksan.

### References

 * Gill, P. E., Murray, W., Saunders, M. A. (2002). SNOPT: An SQP Algorithm for Large-Scale Constrained Optimization, SIAM Journal on Optimization, 12(4), 979-1006.
 * Kraft, D. (1988). A software package for sequential quadratic programming. Forschungsbericht Deutsche Forschungs- und Versuchsanstalt für Luft- und Raumfahrt.
 * Nocedal, J., & Wright, S. J. (2006). Numerical Optimization. Springer.
 * Fletcher, R., & Leyffer, S. (2002). Nonlinear programming without a penalty function. Mathematical Programming, 91(2), 239-269.