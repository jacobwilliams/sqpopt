
![Modern Fortran SQP OPTimizer](media/logo.png)
Modern Fortran SQP OPTimizer

### Goals

A modern Fortran implementation of a Sequential Quadratic Programming (SQP)
optimizer, for problems of the form:

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

Features include:
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
call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)

call solver%initialize(problem=problem, options=options)
call solver%solve(x0, istat)
call solver%get_solution(x, lambda)
```

See [test/test_basic.f90](test/test_basic.f90), [test/test_hs71.f90](test/test_hs71.f90),
and [test/test_medium.f90](test/test_medium.f90) for complete worked examples.

### Configuration

`solver%initialize(problem=..., options=..., hessian=..., qp_solver=...,
linesearch=..., report=...)` accepts one instance of each sub-component,
all optional (defaults are used for anything omitted). `sqpopt_options_type`
covers the most commonly-tuned, algorithm-*selecting* settings and is
copied down into the other components' `mode`-like fields at the start of
every `solve()` call; the other types (`sqpopt_hessian_type`,
`sqpopt_qp_solver_type`, `sqpopt_linesearch_type`) expose further
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
| `set_bounds(x_lb, x_ub, c_lb, c_ub)` | variable bounds and constraint bounds; use `c_lb(i)==c_ub(i)` for an equality constraint, and a large sentinel value (e.g. `1e20`) for a one-sided/absent bound |
| `set_jacobian_sparsity(nnz, irow, icol)` | fixed 1-based COO sparsity pattern of the constraint Jacobian |
| `set_hessian_sparsity(nnz, irow, icol)` | fixed 1-based COO sparsity pattern of the Lagrangian Hessian (only needed for `sqpopt_hessian_exact`, which is not yet implemented) |
| `set_functions(f, g, c, jac, hess)` | attach the user-supplied callbacks: `f`/`g` evaluate the objective and its gradient; `c`/`jac` evaluate the constraints and the Jacobian's nonzero values; `hess` (optional) evaluates the Hessian's nonzero values (only used by `sqpopt_hessian_exact`) |

#### Solver options (`sqpopt_options_type`)

| option | default | description |
|---|---|---|
| `max_iter` | `100` | maximum number of major SQP iterations |
| `print_level` | `0` | amount of diagnostic printing (currently unused) |
| `hessian_mode` | `sqpopt_hessian_bfgs` | Hessian approximation strategy (see the [Hessian approximation](#hessian-approximation-sqpopt_hessian_type) table below) |
| `lbfgs_memory` | `10` | number of `(s,y)` vector pairs retained by the limited-memory Hessian |
| `linear_solver_mode` | `sqpopt_linsolve_lusol` | reserved for future use -- not currently wired to any solver dispatch |
| `qp_solver_mode` | `sqpopt_qp_composite` | QP subproblem algorithm (see the [QP subproblem solver](#qp-subproblem-solver-sqpopt_qp_solver_type) tables below) |
| `linesearch_mode` | `sqpopt_linesearch_armijo` | line search strategy (see the [Line search & merit function](#line-search--merit-function-sqpopt_linesearch_type) tables below) |
| `merit_mode` | `sqpopt_merit_l1` | merit function (see the [Line search & merit function](#line-search--merit-function-sqpopt_linesearch_type) tables below) |
| `ftol`, `xtol`, `ctol` | `1e-8` | convergence tolerances on the relative change in the objective/variables, and on constraint feasibility |
| `ktol` | `1e-6` | tolerance on the KKT optimality (projected-gradient) test |

`hessian_mode`/`lbfgs_memory`, `qp_solver_mode`, `linesearch_mode`, and
`merit_mode` are copied onto the corresponding sub-component (`hessian`,
`qp_solver%mode`, `linesearch%mode`, `linesearch%merit_mode`) at the start
of every `solve()` call, so setting them via `options` and via a directly-
constructed sub-component are equivalent -- everything else in the tables
below (tolerances, step limits, etc.) must be set on the sub-component
directly, since `sqpopt_options_type` has no field for it.

#### Hessian approximation (`sqpopt_hessian_type`)

A limited-memory quasi-Newton approximation (never a dense `n x n`
matrix); `solve()` always (re)initializes it from `options%lbfgs_memory`/
`options%hessian_mode` at the start of every call, so the only way to
configure it is via those two `options` fields:

| option | default | description |
|---|---|---|
| `hessian_mode` (`sqpopt_options_type`) | `sqpopt_hessian_bfgs` | `sqpopt_hessian_bfgs` (limited-memory damped BFGS), `sqpopt_hessian_sr1` (limited-memory symmetric rank-1), or `sqpopt_hessian_exact` (not yet implemented, falls back to BFGS) |
| `lbfgs_memory` (`sqpopt_options_type`) | `10` | number of `(s,y)` vector pairs retained, independent of the problem size `n` |

#### QP subproblem solver (`sqpopt_qp_solver_type`)

| option | default | description |
|---|---|---|
| `mode` | `sqpopt_qp_composite` | which QP algorithm to use (overwritten from `options%qp_solver_mode` at the start of `solve()`) -- see the mode table below |
| `max_step` | `2.0` | trust-region-style cap on \|\|p\|\|₂ applied after every QP solve, regardless of `mode` |
| `active_tol` | `1e-6` | tolerance used by the composite step (`mode=sqpopt_qp_composite`) to decide whether an inequality constraint is part of the active set |
| `bound_enforcement` | `sqpopt_bounds_scalar` | how the composite step (`mode=sqpopt_qp_composite` only) corrects a bound violation in its computed step -- see the mode table below; not used by `sqpopt_qp_dense`/`sqpopt_qp_reduced_hessian`, which enforce bounds exactly as part of the QP solve itself |
| `max_iter`, `linear_solver_mode` | `0`, `sqpopt_linsolve_lusol` | reserved for future use -- not currently wired to any solver dispatch |
| `dense_qp` | -- | the dense active-set QP solver's own options (used only when `mode==sqpopt_qp_dense`); see table below |
| `sparse_qp` | -- | the sparse active-set QP solver's own options (used only when `mode==sqpopt_qp_reduced_hessian`); see table below |

**`mode` values (`qp_solver%mode` / `options%qp_solver_mode`):**

| value | description |
|---|---|
| `sqpopt_qp_composite` | (default) matrix-free composite-step heuristic (multiplier estimate + normal step + tangential step, all via `LSQR`); does not enforce the linearized general-constraint bounds exactly, relying on the outer major iterations to converge to feasibility |
| `sqpopt_qp_dense` | dense active-set QP (Householder QR null-space + modified Cholesky reduced-Hessian solve); enforces bounds/constraints exactly; forms `O(n^2)`/`O(mn)` dense arrays each call, so best for small-to-moderate problems |
| `sqpopt_qp_reduced_hessian` | sparse/matrix-free active-set QP (projected conjugate gradients, `LSQR`-based null-space projections); enforces bounds/constraints exactly; needs more major iterations than `sqpopt_qp_dense` (an `LSQR`-iterative-tolerance cost) but never forms a dense array, so it scales to larger problems |

**`bound_enforcement` values (`qp_solver%bound_enforcement`, `sqpopt_qp_composite` mode only):**

| value | description |
|---|---|
| `sqpopt_bounds_scalar` | (default) clip only the violating components of `x+p` to their bound; the other components of `p` are left unchanged |
| `sqpopt_bounds_vector` | rescale the *entire* step `p` by the same factor so that `x+p` just touches the first bound it would otherwise violate, preserving `p`'s direction exactly |

**Dense active-set QP options (`qp_solver%dense_qp`, used when `mode==sqpopt_qp_dense`):**

| option | default | description |
|---|---|---|
| `max_iter` | `100` | maximum number of active-set changes (add/drop a row) allowed per QP solve |
| `active_tol` | `1e-8` | tolerance used to detect an (in)active/equality row |
| `opt_tol` | `1e-8` | tolerance on the reduced-gradient stationarity test |

**Sparse (projected-CG) active-set QP options (`qp_solver%sparse_qp`, used when `mode==sqpopt_qp_reduced_hessian`):**

| option | default | description |
|---|---|---|
| `max_iter` | `100` | maximum number of active-set changes (add/drop a row) allowed per QP solve |
| `max_pcg_iter` | `0` | maximum projected-CG iterations per active-set face (`<=0` means "use `n`") |
| `active_tol` | `1e-8` | tolerance used to detect an (in)active/equality row |
| `opt_tol` | `1e-8` | tolerance on the projected-residual stationarity test |
| `lsqr_atol`, `lsqr_btol`, `lsqr_conlim` | `0.0` | `LSQR` relative error tolerances in `A`/`b`, and the upper limit on `cond(Abar)` (`0` means "let `LSQR` use its own machine-precision-based default", which is tighter than usually necessary); loosening these is the main lever for trading QP-solve accuracy for speed in this mode |
| `lsqr_itnlim` | `100` | `LSQR` maximum iterations per solve |

#### Line search & merit function (`sqpopt_linesearch_type`)

| option | default | description |
|---|---|---|
| `mode` | `sqpopt_linesearch_armijo` | line search strategy (overwritten from `options%linesearch_mode` at the start of `solve()`) -- see the mode table below |
| `merit_mode` | `sqpopt_merit_l1` | merit function (overwritten from `options%merit_mode` at the start of `solve()`) -- see the mode table below |
| `penalty` | `1.0` | current penalty parameter used in the merit function (\( \mu \) for `sqpopt_merit_l1`, \( \rho \) for `sqpopt_merit_augmented_lagrangian`) |
| `major_step_limit` | `2.0` | caps the *initial* trial step length (before any backtracking), used by all three modes, so that no variable changes by more than this factor relative to \( \max(1,\lvert x_j\rvert) \) (SNOPT's "Major step limit" option); guards against divergence from a QP step that is technically feasible but unreasonably large |
| `sigma` | `0.1` | Armijo sufficient-decrease parameter, \( 0<\sigma<1 \) (`sqpopt_linesearch_armijo`/`sqpopt_linesearch_watchdog` modes) |
| `backtrack` | `0.5` | step-length reduction factor at each backtracking step (`sqpopt_linesearch_armijo`/`sqpopt_linesearch_watchdog` modes) |
| `alpha_min` | `0.1` | minimum step length; backtracking never goes below this, and accepts it outright if the floor is reached without satisfying the sufficient-decrease test (`sqpopt_linesearch_armijo`/`sqpopt_linesearch_watchdog` modes) |
| `max_ls_iter` | `20` | maximum number of backtracking steps (`sqpopt_linesearch_armijo`/`sqpopt_linesearch_watchdog` modes) |
| `tol` | `1e-4` | desired tolerance on the minimizer (`sqpopt_linesearch_exact` mode) |
| `watchdog_relaxed_len` | `2` | number of relaxed steps tolerated before requiring a new best point (`sqpopt_linesearch_watchdog` mode) |
| `watchdog_cooldown_len` | `10` | number of iterations relaxed acceptance is disabled for after a backtrack (`sqpopt_linesearch_watchdog` mode) |

**`mode` values (`linesearch%mode` / `options%linesearch_mode`):**

| value | description |
|---|---|
| `sqpopt_linesearch_armijo` | (default) standard backtracking line search with an Armijo-type sufficient-decrease test on the merit function (as used by default in `slsqp`) |
| `sqpopt_linesearch_exact` | (approximately) minimizes the merit function along the search direction using the derivative-free `fmin` routine |
| `sqpopt_linesearch_watchdog` | Powell's watchdog technique: tracks the best point found so far and, for a short window after a genuine improvement, relaxes the sufficient-decrease test (accepting the full step outright) rather than stalling near a curved/simultaneously-active constraint boundary (the Maratos effect); backtracks to the best point and disables relaxed acceptance for `watchdog_cooldown_len` iterations if the window is used up without a new best point |

**`merit_mode` values (`linesearch%merit_mode` / `options%merit_mode`):**

| value | description |
|---|---|
| `sqpopt_merit_l1` | (default) non-smooth \( \ell_1 \) exact penalty function (as in `slsqp`) |
| `sqpopt_merit_augmented_lagrangian` | smooth augmented Lagrangian merit function (Gill, Murray, Saunders & Wright; the merit function used in NPSOL and, in spirit, SNOPT) -- twice continuously differentiable, which avoids the Maratos effect without needing a second-order correction |

See [PLAN.md](PLAN.md) for the full architecture write-up, algorithm
details, and backlog of future work.

### to Build

Use the `pixi` environment and FPM:

```
pixi shell
fpm build --profile release
fpm test --profile release
```


