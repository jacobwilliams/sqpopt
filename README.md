
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

### Available options (`sqpopt_options_type`)

| option | default | description |
|---|---|---|
| `max_iter` | `100` | maximum number of major SQP iterations |
| `hessian_mode` | `sqpopt_hessian_bfgs` | Hessian approximation: `sqpopt_hessian_bfgs` (limited-memory BFGS), `sqpopt_hessian_sr1` (limited-memory SR1), or `sqpopt_hessian_exact` (not yet implemented, falls back to BFGS) |
| `lbfgs_memory` | `10` | number of `(s,y)` vector pairs retained by the limited-memory Hessian |
| `linear_solver_mode` | `sqpopt_linsolve_lusol` | sparse linear solver used by the QP subproblem: `sqpopt_linsolve_lusol`, `sqpopt_linsolve_lsqr`, or `sqpopt_linsolve_lsmr` |
| `qp_solver_mode` | `sqpopt_qp_composite` | QP subproblem algorithm: `sqpopt_qp_composite` (matrix-free composite-step heuristic, does not enforce linearized bounds exactly), `sqpopt_qp_dense` (dense active-set QP, enforces bounds exactly, best for small/moderate problems), or `sqpopt_qp_reduced_hessian` (sparse/matrix-free projected-CG active-set QP, enforces bounds exactly, scales to larger problems) |
| `linesearch_mode` | `sqpopt_linesearch_armijo` | line search strategy: `sqpopt_linesearch_armijo` (backtracking, as in `slsqp`), `sqpopt_linesearch_exact` (derivative-free 1-D minimization via `fmin`), or `sqpopt_linesearch_watchdog` (Powell's VF13 watchdog technique: relaxed step acceptance with a best-point safety net) |
| `merit_mode` | `sqpopt_merit_l1` | merit function: `sqpopt_merit_l1` (non-smooth exact penalty) or `sqpopt_merit_augmented_lagrangian` (smooth NPSOL/SNOPT-style augmented Lagrangian) |
| `ftol`, `xtol`, `ctol` | `1e-8` | convergence tolerances on the objective, variables, and constraint feasibility |
| `ktol` | `1e-6` | tolerance on the KKT optimality (projected-gradient) test |
| `print_level` | `0` | amount of diagnostic printing (currently unused) |

The QP subproblem solver (`sqpopt_qp_solver_type`) and line search
(`sqpopt_linesearch_type`) each expose further tuning parameters (e.g.
`qp_solver%max_step`, a trust-region-style cap on the step norm;
`qp_solver%active_tol`, the active-set filter tolerance;
`qp_solver%bound_enforcement` (`sqpopt_qp_composite` mode only), how the
composite step's bound violations are corrected -- `sqpopt_bounds_scalar`
(default, clip only the violating components) or `sqpopt_bounds_vector`
(rescale the whole step uniformly, preserving its direction);
`linesearch%alpha_min`/`sigma`/`backtrack`, the Armijo parameters) that can
be set by constructing them directly and passing them to
`solver%initialize(problem=..., options=..., qp_solver=..., linesearch=...)`.

See [PLAN.md](PLAN.md) for the full architecture write-up, algorithm
details, and backlog of future work.

### to Build

Use the `pixi` environment and FPM:

```
pixi shell
fpm build --profile release
fpm test --profile release
```


