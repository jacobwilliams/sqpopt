
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
call problem%set_functions(fc=fc, gjac=gjac)         ! optionally also data=...

call solver%initialize(problem=problem, options=options)
call solver%solve(x0, istat)            ! optionally also lambda0=...
call solver%get_solution(x, lambda)     ! optionally also z (bound multipliers)
print *, solver%status_message()        ! e.g. 'converged successfully'
```

The problem is defined by two user functions: `fc` evaluates the
objective and the constraints together, and `gjac` evaluates the
objective gradient and the nonzero values of the constraint Jacobian
together (in the order of the `irow`/`icol` sparsity pattern):

```fortran
subroutine fc(x, f, c, status, data)
    real(wp), dimension(:), intent(in)    :: x
    real(wp),               intent(out)   :: f
    real(wp), dimension(:), intent(out)   :: c        ! dimension(m)
    integer,                intent(inout) :: status   ! 0 on entry
    class(*), optional,     intent(inout) :: data
    f = ...
    c = ...
end subroutine fc

subroutine gjac(x, g, jac_val, status, data)
    real(wp), dimension(:), intent(in)    :: x
    real(wp), dimension(:), intent(out)   :: g        ! dimension(n)
    real(wp), dimension(:), intent(out)   :: jac_val  ! dimension(nnz)
    integer,                intent(inout) :: status   ! 0 on entry
    class(*), optional,     intent(inout) :: data
    g = ...
    jac_val = ...
end subroutine gjac
```

The two trailing arguments are:

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
- the number of calls of `fc` and `gjac` (`n_eval_fc`, `n_eval_gjac`), and the run time;
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
| `sqpopt_max_evals_reached` (`9`) | `max_evals` calls of `fc` were performed |
| `sqpopt_time_limit_reached` (`10`) | the `max_time` limit was reached |
| `sqpopt_unbounded` (`11`) | the objective fell below `obj_lower_limit` at a feasible point |
| `sqpopt_acceptable` (`12`) | the looser `acceptable_ktol`/`acceptable_ctol` tests held for `acceptable_iter` consecutive iterations (as in IPOPT), but the normal ones did not |

See [test/test_basic.f90](test/test_basic.f90), [test/test_hs71.f90](test/test_hs71.f90),
and [test/test_medium.f90](test/test_medium.f90) for complete worked examples.

### Configuration

`solver%initialize(problem=..., options=..., hessian=..., qp_solver=...,
linesearch=..., trust_region=..., report=...)` accepts one instance of
each sub-component, all optional (defaults are used for anything
omitted). `sqpopt_options_type` covers the most commonly tuned settings,
including the algorithm selectors; the other types expose further
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

The full reference for every option (problem definition, solver options,
Hessian approximation, QP subproblem solvers, line search and merit
functions, and trust region), along with benchmark results, is in the
[user guide](web/index.html).

See [plan/PLAN.md](plan/PLAN.md) for the full architecture write-up, algorithm
details, and backlog of future work.

### Developing

Use the `pixi` environment and the Fortran Package Manager (FPM):

```
pixi shell
fpm build --profile release
fpm test --profile release
```

The user guide is in [web/](web/). To generate the FORD API documentation
that it links to:

```
ford ford.md --output_dir web/api
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

### Other Optimization Libraries

 * [IPOPT](https://github.com/coin-or/Ipopt) -- Interior Point OPTimizer for large-scale nonlinear optimization.
 * [Uno](https://github.com/cvanaret/Uno) -- Uno (Unifying Nonlinear Optimization) is a C++ framework for solving nonlinearly constrained optimization problems

### References

 * Gill, P. E., Murray, W., Saunders, M. A. (2002). SNOPT: An SQP Algorithm for Large-Scale Constrained Optimization, SIAM Journal on Optimization, 12(4), 979-1006.
 * Kraft, D. (1988). A software package for sequential quadratic programming. Forschungsbericht Deutsche Forschungs- und Versuchsanstalt für Luft- und Raumfahrt.
 * Nocedal, J., & Wright, S. J. (2006). Numerical Optimization. Springer.
 * Fletcher, R., & Leyffer, S. (2002). Nonlinear programming without a penalty function. Mathematical Programming, 91(2), 239-269.
 * D. Kiessling, S. Leyffer, C. Vanaret, A Unified Funnel Restoration SQP Algorithm, Mathematical Programming, Volume 217, pages 323-367 (2026)
 * Gill, P E; Murray, W; Saunders, M A; Wright, M H, Some Theoretical Properties of an Augmented Lagrangian Merit Function, SOL-86-6, 1 April 1986
 * P. E. Gill and E. Wong, User's Guide for SQOPT Version 7.7: Software for Large-Scale Linear and Quadratic Programming, Mar 2021
 * R. Fletcher and S. Leyffer, User manual for filterSQP, University of Dundee, April 1998
 * C. Vanaret1, S. Leyffer, Implementing a unified solver for nonlinearly constrained optimization, Mathematical Programming Computation, 10 June 2026