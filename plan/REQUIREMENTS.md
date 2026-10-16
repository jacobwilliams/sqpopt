# SQPOPT Software Requirements Specification

| | |
|---|---|
| **Product** | SQPOPT, a sequential quadratic programming solver for sparse, nonlinearly constrained optimization |
| **Version** | 0.1.0 (`fpm.toml`) |
| **Baseline** | commit `ac7eb7c5` (2026-10-15) |
| **Status** | Derived from the implementation: each requirement states what the code does at the baseline |

## 1. Introduction

### 1.1 Purpose

This document states the requirements that SQPOPT satisfies, for use in its development process: as the
reference for reviews, for the verification of changes, and for tracing each behavior to the code that
implements it and the test that verifies it. The requirements were derived from the source code, its
documentation, and its tests at the baseline above. A change that alters a behavior stated here also
changes this document.

### 1.2 Scope

The requirements cover:

* the Fortran library (`src/`): its programming interface, its algorithm's externally visible behavior,
  its outputs, and the constraints on its implementation;
* the optional least-squares interface (`sqpopt_nlls_module`);
* the Python bindings and options dialog (`python/`);
* the build configurations and the verification assets (`test/`, `tools/`, the CI workflow).

They do not cover the user guide's content (`web/`), the benchmarking programs (`example/`,
`tools/cutest_benchmark.py`), or the design documents (`plan/`), except where a requirement refers to them.

### 1.3 Conventions

* **"shall"** marks a requirement. Each requirement has a unique identifier `AREA-NNN`, which is not reused
  if the requirement is deleted.
* **Verification method** (column **V**): **T** test (an automated test in `test/` or `python/tests/`),
  **I** inspection (review of the code or of a tool's output), **A** analysis, **D** demonstration (running
  a program or tool and observing the result).
* **Trace** names the main implementing source (module or file) and, for **T**, the verifying test.
* Code entities are written as in the source: `options%ktol` is the component `ktol` of an
  `sqpopt_options_type`; status codes are referred to by their named constants, never by value.

### 1.4 Definitions

| Term | Meaning |
|---|---|
| NLP | The nonlinear program \( \min f(x) \) s.t. \( c_l \le c(x) \le c_u \), \( x_l \le x \le x_u \), with \( x \in \mathbb{R}^n \), \( c: \mathbb{R}^n \to \mathbb{R}^m \) |
| QP subproblem | The quadratic program solved at each major iteration for the search direction |
| Major iteration | One SQP iteration: evaluation, convergence tests, Hessian update, QP solve, and step |
| KKT error | The stationarity and complementarity residual measured by `check_convergence` |
| COO | Coordinate sparse storage: 1-based (row, column) index pairs with a value per nonzero |
| Infinite bound | A bound of magnitude at least `sqpopt_infinity` (\( 10^{20} \)) |
| User function | A procedure supplied by the caller: `fc`, `gjac`, `hess`, or `report` |
| HS suite | The 305 Hock–Schittkowski problems of `test/hs_problems.f90`, run by `test/test_hs_suite.f90` |
| Default build | The library built by `fpm` with no preprocessor directives |
| MUMPS build / LAPACK build | The library built with `HAS_MUMPS` / `HAS_LAPACK` defined |

### 1.5 References

* `README.md`; the user guide, `web/index.html`; `CLAUDE.md` (development rules).
* `plan/ROADMAP.md` (open issues and ideas, e.g. F34).
* The FORD API documentation generated from the source (`ford.md`).

## 2. Overall description

### 2.1 Product perspective

SQPOPT is a library. A calling program defines an NLP through an `sqpopt_problem_type`, configures the solver
through option and component types, and calls `sqpopt_type%solve`. The library has no executable of its own.
It is layered as follows: the solver (`sqpopt_module`) runs major iterations (`sqpopt_iterate_module`), each
of which uses separately replaceable components (Hessian, QP solvers, line search, merit function, filter,
funnel, trust region, restoration, second-order correction, convergence test, sparse factorizations,
diagnostics, log). The least-squares interface and the Python bindings are layers on top of the solver that
the solver does not depend on.

### 2.2 Users

* Fortran programmers who embed the solver in an application.
* Python users, through a `scipy.optimize`-like interface.
* Developers of SQPOPT, through the tests, tools, and CI.

### 2.3 Operating environment

* A Fortran 2008 compiler (with the C preprocessor for the `.F90` files); developed with gfortran 15 on macOS
  (arm64) and Linux (x86-64); the pixi environment also declares Windows (x86-64).
* Build system: fpm 0.13, inside the pixi environment (`pixi.toml`, `pixi.lock`).
* Python 3.14 with NumPy and SciPy for the bindings; Qt (through qtpy and PySide6) for the options dialog.

### 2.4 Dependencies

| Dependency | Use | Required |
|---|---|---|
| LSQR | Minimum-norm least-squares solves (restoration, SOC, sparse QP null space) | Yes (fpm) |
| LUSOL | Sparse LU (sparse QP basis, independent columns) | Yes (fpm) |
| fmin | The exact line search | Yes (fpm) |
| qdldl-fortran | The default sparse \( LDL^T \) factorization | Yes (fpm) |
| daqp-fortran | The DAQP QP mode | Yes (fpm) |
| slsqp | The SLSQP comparison test | Tests only (fpm dev-dependency) |
| MUMPS (sequential, double) | Optional sparse solver | No (`HAS_MUMPS`) |
| LAPACK and BLAS | Optional dense factorization, eigensolver, products | No (`HAS_LAPACK`) |

### 2.5 Assumptions

* The user functions are deterministic, and their derivatives are consistent with their values (the solver
  has no finite differences of its own).
* The sparsity patterns of the Jacobian and of the Hessian do not change during a solve.

## 3. Functional requirements

### 3.1 Problem definition (PRB)

| ID | Requirement | V | Trace |
|---|---|---|---|
| PRB-001 | The library shall accept NLPs with \( n > 0 \) variables and \( m \ge 0 \) general constraints \( c_l \le c(x) \le c_u \), and bounds \( x_l \le x \le x_u \) on the variables. | T | `sqpopt_problem_module`; `test_basic`, `test_hs_suite` |
| PRB-002 | The problem size shall be set by `set_problem_size(n, m)`, which shall allocate the bound arrays and remove any step limits previously set. | T | `set_problem_size`; `test_input_validation` |
| PRB-003 | The bounds shall be set by `set_bounds(x_lb, x_ub, c_lb, c_ub)`. A constraint with `c_lb(i) == c_ub(i)` shall be an equality, any other an inequality; equalities and inequalities may be in any order. | T | `set_bounds`; `test_basic`, `test_hs71` |
| PRB-004 | A bound of magnitude at least `sqpopt_infinity` shall be treated as absent, and shall be clamped to \( \pm \)`sqpopt_infinity` so that `-huge(1.0_wp)` and `huge(1.0_wp)` can be used without floating-point overflow. | T | `validate`; `test_infinite_bounds` |
| PRB-005 | The constraint Jacobian's sparsity pattern shall be given once, as 1-based COO index pairs, by `set_jacobian_sparsity(nnz, irow, icol)`; the user's `gjac` shall return only the values of those nonzeros, in that order. | T | `set_jacobian_sparsity`; `test_large_sparse` |
| PRB-006 | The Hessian of the Lagrangian's sparsity pattern shall be given by `set_hessian_sparsity(nnz, irow, icol)`. Each off-diagonal element shall be given once, in either triangle, and stand for both \( (i,j) \) and \( (j,i) \); repeated entries shall be summed. | I | `set_hessian_sparsity`, `hessian_dense` |
| PRB-007 | The user shall be able to limit the change of each variable in one major iteration by `set_max_step(max_step)`; the limits shall be enforced as bounds on every step the solver computes (QP, second-order correction, restoration). A limit of at least `sqpopt_infinity` shall mean no limit. | T | `set_max_step`, `step_bounds`; `test_max_step` |
| PRB-008 | No dense \( m \times n \) Jacobian or \( n \times n \) Hessian shall be required from, or formed for, the user: derivatives shall be exchanged as COO values. | I | `sqpopt_problem_module` |

### 3.2 User functions (USR)

| ID | Requirement | V | Trace |
|---|---|---|---|
| USR-001 | The problem shall be evaluated by two user procedures attached by `set_functions(fc, gjac, hess, data)`: `fc(x, f, c, status, data)` for the objective and constraints, and `gjac(x, g, jac_val, accuracy, status, data)` for the gradient and the Jacobian's nonzeros. | T | `sqpopt_fc_func`, `sqpopt_gjac_func`; `test_basic` |
| USR-002 | An optional `hess(x, lambda, hess_val, status, data)` shall supply the nonzeros of the Hessian of the Lagrangian \( \nabla^2 f - \sum_i \lambda_i \nabla^2 c_i \), for `options%hessian_mode = sqpopt_hessian_exact`. | T | `sqpopt_hessian_func`; `test_inertia`, `test_hs_suite --hessian=exact` |
| USR-003 | An optional user data object (`class(*)`, with the `target` attribute) given to `set_functions` shall be passed by reference to every user function and to the `report` callback; changes the functions make to it shall be visible to the caller. | T | `set_functions`; `test_callbacks` |
| USR-004 | Every user function shall receive `status = 0`. If it returns `status > 0` (cannot evaluate at `x`), its outputs shall be treated as NaN. | T | `check_status`, `raw_fc`; `test_callbacks` |
| USR-005 | If `fc` cannot be evaluated, or returns a non-finite value, at a trial point, the solver shall reject the point and try a shorter or different step. | T | line search, trust region; `test_nonfinite`, `test_callbacks` |
| USR-006 | If any user function returns `status < 0`, the solver shall stop with `sqpopt_user_requested_stop`, and shall not call any user function again during that solve. | T | `sqpopt_problem_module`, `sqpopt_iterate`; `test_callbacks` |
| USR-007 | The solver shall pass `gjac` an `accuracy` argument, `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`, as set by `options%derivative_accuracy` (REQ ALG-012). | T | `set_derivative_accuracy`; `test_derivative_accuracy` |
| USR-008 | An optional `report(iter, x, f, c, lambda, user_stop, data)` callback, given to `initialize`, shall be called once per major iteration, before any step is taken, including at the starting point and at the final point, with the values of the original (unscaled) problem. If it sets `user_stop = .true.`, the solver shall stop with `sqpopt_user_requested_stop`. | T | `sqpopt_iterate`; `test_callbacks` |
| USR-009 | The solver shall never evaluate the user functions at a point outside the variable bounds: the starting point shall be projected onto them, and every step shall satisfy them. | T | `sqpopt_solve`; `test_input_validation` |
| USR-010 | The solver shall not call `fc` again at a point at which one of its last four evaluations was made, nor `gjac` at the point of its last evaluation (bitwise-equal points). The reported evaluation counts shall equal the number of calls of the user functions. | I, T | `raw_fc`, `raw_gjac`; `test_results` (the counts against the user functions' own counters) |

### 3.3 Solver interface and life cycle (API)

| ID | Requirement | V | Trace |
|---|---|---|---|
| API-001 | The solver shall be the type `sqpopt_type`, with the public procedures `initialize`, `solve`, `get_solution`, `get_results`, `status_message`, and `destroy`; its components shall be private. | I | `sqpopt_module` |
| API-002 | `initialize(problem, options, hessian, qp_solver, linesearch, trust_region, report)` shall take every argument as optional, and use a default-initialized object for each one that is absent. | T | `sqpopt_initialize`; `test_basic` |
| API-003 | `solve(x0, istat, lambda0)` shall solve the NLP from the starting point `x0`, optionally with starting multipliers `lambda0` (sign convention of RES-006), and return a status code in `istat`. | T | `sqpopt_solve`; `test_basic`, `test_results` |
| API-004 | Every `solve` shall start from the configuration given to `initialize`: no state (scaling, evaluation caches, penalty parameter, filter, funnel, watchdog, trust-region radius, Hessian approximation, QP warm start, restoration) shall carry over from a previous solve. Two solves from the same starting point shall give identical results. | T | `sqpopt_solve`; `test_resolve`, `test_nlls` |
| API-005 | `get_solution(x, lambda, z)` shall return the final point, the constraint multipliers, and optionally the variable-bound multipliers of the last solve, for the original problem. | T | `sqpopt_get_solution`; `test_results` |
| API-006 | `get_results(results)` shall return the complete `sqpopt_results_type` of the last solve (section 3.12). | T | `sqpopt_get_results`; `test_results` |
| API-007 | `status_message()` shall return the description of the last solve's status, including, for `sqpopt_invalid_input`, what was invalid. | T | `sqpopt_get_status_message`; `test_input_validation` |
| API-008 | `destroy()` shall release the solver's allocated memory and return it to the default configuration. | I | `sqpopt_destroy` |
| API-009 | Separate `sqpopt_type` objects shall be independent: each shall keep its own configuration, state, and log settings, with no module-level mutable state. | I | `sqpopt_module`, `sqpopt_log_module` |
| API-010 | The status codes, the results type, the callback interface `sqpopt_report_func`, and `sqpopt_infinity` shall be public in `sqpopt_types_module`, with `sqpopt_status_message(istat)` giving each code's description. | T | `sqpopt_types_module`; `test_results` |

### 3.4 Input validation (VAL)

| ID | Requirement | V | Trace |
|---|---|---|---|
| VAL-001 | `solve` shall check the problem definition and every option before evaluating any user function. On the first invalid input found it shall return `sqpopt_invalid_input`, with a message naming the input, without calling a user function. | T | `validate`, `validate_options`; `test_input_validation` |
| VAL-002 | The problem shall be invalid if: \( n \le 0 \); \( m < 0 \); the bounds are not set or have the wrong sizes; a bound is NaN; a lower bound exceeds its upper bound; a lower bound is \( +\infty \) or an upper bound \( -\infty \); a step limit has the wrong size or is not positive; the Jacobian pattern is missing (with \( m > 0 \) or nonzeros), has the wrong size, or has an index outside \( 1..m \) / \( 1..n \); the Hessian pattern is missing (with nonzeros), has the wrong size, or has an index outside \( 1..n \); `fc` or `gjac` is not set. | T | `validate`; `test_input_validation` |
| VAL-003 | `x0` shall be invalid if its size is not `n` or it is not finite; `lambda0` (if present) shall be invalid if its size is not `m` or it is not finite. Finiteness of `x0` shall be checked before the projection onto the bounds. | T | `sqpopt_solve`; `test_input_validation` |
| VAL-004 | The options shall be invalid if outside the ranges of Table A.2, or if an enumerated option is not one of its named constants. | T | `validate_options`; `test_input_validation` |
| VAL-005 | `options%hessian_mode = sqpopt_hessian_exact` shall be invalid input unless the `hess` function is set. | T | `validate_options`; `test_input_validation` |
| VAL-006 | `options%linear_solver = sqpopt_linear_solver_mumps` shall be invalid input in a library built without MUMPS, and `sqpopt_linear_solver_lapack` (and `hessian%eigen_solver = sqpopt_eigen_lapack`) in a library built without LAPACK. | T | `validate_options`; `test_input_validation` |
| VAL-007 | When a factorization option is in use (`inertia_control` with a non-BFGS Hessian, `direct_qp`, or `direct_least_squares` with \( m > 0 \)), a linear solver whose order limit is below \( n + m \) shall be invalid input (the dense solvers' limit is 2,000). | T | `validate_options`; `test_input_validation` |
| VAL-008 | With `print_level > 0`, an `output_unit` that is not open shall be invalid input; with `diagnostic_level >= 2`, a `diagnostics_unit` other than `-1` that is not open shall be invalid input. | T | `validate_options`; `test_input_validation` |
| VAL-009 | The component settings shall be validated as in Table A.3 (line search, filter, funnel, trust region when enabled, QP solvers, Hessian). | T | `validate_options`; `test_input_validation` |
| VAL-010 | On invalid input, the results shall have zero evaluation counts, zero `f`, and zero-filled `c`, `lambda`, and `z` of the problem's sizes, and nothing shall be printed. | I | `finish` |

### 3.5 Major iteration and algorithm (ALG)

| ID | Requirement | V | Trace |
|---|---|---|---|
| ALG-001 | Each major iteration shall evaluate \( f, \nabla f, c, J \) at the current point, test for convergence and the other stopping conditions, update the Hessian approximation, solve the QP subproblem for a search direction, and take a step by the selected globalization. | I | `sqpopt_iterate` |
| ALG-002 | The QP subproblem shall be \( \min_p \tfrac12 p^T H p + g^T p \) s.t. \( c_l - c \le J p \le c_u - c \) and the variable bounds (tightened by the step limits), with its multipliers used as the new multiplier estimates. | I | `sqpopt_qp_solver_module` |
| ALG-003 | The QP step shall be capped at \( \lVert p \rVert_2 \le \) `qp_solver%max_step` \( \times \) `step_scale`, where `step_scale` starts at \( \max(1, \lVert x_0 \rVert_\infty / \) `max_step` \( ) \) and adapts during the solve. | I | `sqpopt_solve`, `sqpopt_qp_solver_module` |
| ALG-004 | If a QP solve reaches its iteration limit, the iteration shall use the QP's last step. If no acceptable step is found, the point shall be left unchanged and the Hessian approximation reset (or the exact Hessian's shift increased). After `options%max_consecutive_failures` consecutive iterations with either failure, the solve shall stop with that iteration's status (`sqpopt_qp_solve_failed` or `sqpopt_line_search_failed`). | I | `sqpopt_solve`, `sqpopt_iterate` |
| ALG-005 | If a constraint's multiplier diverges (its \( |\lambda_i| \lVert \nabla c_i \rVert_\infty \) exceeds `options%elastic_multiplier_limit` \( \times \max(1, \lVert g \rVert_\infty) \) and keeps growing), the QP shall be re-solved with that constraint elastic. `elastic_multiplier_limit = 0` shall disable this. | T | `sqpopt_iterate`; `test_degenerate` |
| ALG-006 | Before reporting `sqpopt_infeasible` at a point that is stationary for the constraint violation, the solver shall look for a second-order decrease of the violation (an escape step) and continue from it if found, at most 3 times per solve. | T | `escape_step`; `test_escape` |
| ALG-007 | With a feasible point whose objective is below `options%obj_lower_limit` (default \( -\infty \): off), the solve shall stop with `sqpopt_unbounded`. | T | `sqpopt_iterate`; `test_termination` |
| ALG-008 | If \( f \), \( \nabla f \), \( c \), \( J \), or the exact Hessian is not finite at the current (accepted or starting) point, the solve shall stop with `sqpopt_function_error`. | T | `sqpopt_iterate`; `test_nonfinite` |
| ALG-009 | With the exact Hessian, after a step whose Hessian shift dominated the step, the multipliers shall be replaced by a least-squares estimate before the next Hessian evaluation. | T | `estimate_multipliers`; `test_multipliers` |
| ALG-010 | With a limited-memory BFGS Hessian and the QP front end's `unconstrained_step` enabled (default), the solver shall first try the unconstrained minimizer \( -H^{-1} g \) of the QP's objective, and accept it without running a QP solver if it satisfies the bounds and the linearized constraints. With SR1, the same shall apply when the matrix is positive definite, with the step limited to the step cap. | T | `sqpopt_qp_solver_module`; `test_unconstrained_step` |
| ALG-011 | The QP solver shall be selected by `options%qp_solver_mode` (section 3.7). | T | `sqpopt_qp_solver_module`; `test_hs_suite --qp=...` |
| ALG-012 | With `options%derivative_accuracy = sqpopt_derivatives_fast`, the solver shall request fast derivatives until the KKT and feasibility errors are below `derivative_switch_tol`, or a stopping condition would hold, or a step fails, or progress stalls; it shall then request accurate derivatives for the rest of the solve, re-evaluate the derivatives at the current point, and repeat the convergence tests with them. | T | `sqpopt_iterate`; `test_derivative_accuracy` |
| ALG-013 | No dense \( n \times n \) or \( m \times n \) array shall be formed in a major iteration, except by the dense QP solver, the DAQP mode, and the dense linear solvers (section 4.2). | I | `sqpopt_iterate_module`, `tools/stack_check.sh` |

### 3.6 Convergence and termination (CNV)

| ID | Requirement | V | Trace |
|---|---|---|---|
| CNV-001 | The solve shall stop with `sqpopt_success` at a point where, for the scaled problem: the projected stationarity residual of the Lagrangian \( g - J^T \lambda \) and the complementarity residual are at most `ktol` \( \times s_d \); the constraint and bound violations are at most `ctol`; and the unscaled stationarity residual is at most `dual_inf_tol`. | T | `check_convergence`; `test_convergence`, `test_basic` |
| CNV-002 | The stationarity residual shall be projected onto the variable bounds: at a bound, only a component pointing out of the feasible region shall count; a fixed variable shall have no stationarity requirement. | T | `projected_inf_norm`; `test_convergence` |
| CNV-003 | For an inequality constraint, a multiplier of the wrong sign for the bound the constraint is at (\( \lambda_i > 0 \) only at the lower bound, \( \lambda_i < 0 \) only at the upper bound) shall count as a violation of \( \lambda_i \min(1, \text{distance to that bound}) \); equality multipliers shall be free. | T | `check_convergence`; `test_convergence`, `test_kkt_sign` |
| CNV-004 | The scale \( s_d = \max(100, \lVert \lambda \rVert_1 / m)/100 \) shall loosen the stationarity and complementarity tolerances only when the multipliers are large on average. | I | `check_convergence` |
| CNV-005 | The solve shall stop with `sqpopt_acceptable` when, for `acceptable_iter` (default 15) consecutive iterations, the KKT and feasibility errors are below `acceptable_ktol` and `acceptable_ctol` and the objective changed by at most `acceptable_obj_change_tol` relative to \( \max(1, |f|) \). `acceptable_iter = 0` shall disable this. | T | `sqpopt_iterate`; `test_termination` |
| CNV-006 | The solve shall stop with `sqpopt_stalled` when, at a feasible point, for `stall_iter` (default 3) consecutive iterations, the relative change of the objective is at most `ftol` and the relative change of the variables, \( \lVert x - x_{prev} \rVert_2 / \max(1, \lVert x \rVert_2) \), is at most `xtol`. This test shall be skipped on the first iteration and right after a failed step. | T | `check_convergence`, `sqpopt_iterate`; `test_hs_suite` |
| CNV-007 | The solve shall stop with `sqpopt_infeasible` at an infeasible point (violation above `ctol`) that has not reduced the violation by 1% since the previous iterate and is stationary for the squared constraint violation subject to the bounds (projected \( J^T r_c \) at most `ktol` \( \lVert r_c \rVert_2 \)), after the escape steps of ALG-006. The test shall not be made at the starting point. | T | `check_convergence`; `test_infeasible` |
| CNV-008 | The solve shall stop with `sqpopt_max_iter_reached` after `max_iter` (default 100) major iterations; with `sqpopt_max_evals_reached` once the calls of `fc` reach `max_evals` (\( > 0 \)); and with `sqpopt_time_limit_reached` once the wall-clock time reaches `max_time` (\( > 0 \) seconds). | T | `sqpopt_solve`; `test_termination` |
| CNV-009 | The tolerances `ktol` and `ctol` shall apply to the scaled problem (section 3.10). | T | `check_convergence`; `test_results` |

### 3.7 QP subproblem solvers (QP)

| ID | Requirement | V | Trace |
|---|---|---|---|
| QP-001 | `sqpopt_qp_auto` (default) shall select the dense QP solver for \( n \le \) `qp_solver%auto_dense_max_n` (default 200), and the sparse QP solver otherwise. | T | `sqpopt_qp_solver_module`; `test_large_sparse` |
| QP-002 | `sqpopt_qp_dense` shall be a primal active-set method with a null-space step (Householder QR and Cholesky), treating the general constraints and bounds as two-sided rows, adding only linearly independent rows to the working set. | T | `sqpopt_qp_dense_module`; `test_qp_dense`, `test_qp_fuzz` |
| QP-003 | The dense QP solver shall handle inconsistent linearized constraints by an elastic mode: an \( \ell_1 \)-penalized slack for each violated row, with the weight increased up to `elastic_weight_max`; if slacks remain positive, it shall return `sqpopt_infeasible` with the elastic solution. | T | `sqpopt_qp_dense_module`; `test_qp_dense`, `test_infeasible` |
| QP-004 | The dense QP solver shall follow a direction of nonpositive curvature of the reduced Hessian to the nearest blocking row instead of taking a Newton step, and report negative curvature in the variables. | I | `sqpopt_qp_dense_module` |
| QP-005 | The dense QP solver shall warm-start from the previous solve's working set when the problem size is unchanged (`warm_start`, default on), else start from the equality constraints and fixed variables. | T | `sqpopt_qp_dense_module`; `test_qp_fuzz` |
| QP-006 | If the dense QP solver cannot allocate its matrices, it shall return `sqpopt_out_of_memory` with \( p = 0 \), \( \lambda = 0 \), and the solve shall stop with `sqpopt_out_of_memory`. | T | `sqpopt_qp_dense_module`; `test_out_of_memory`, `test_restoration_out_of_memory` |
| QP-007 | `sqpopt_qp_reduced_hessian` shall be a sparse active-set method using matrix-free Hessian products and projected conjugate gradients on the reduced Hessian, with a null-space basis from a sparse LU partition (`sqpopt_null_space_lu`, default) or from LSQR projections (`sqpopt_null_space_lsqr`). | T | `sqpopt_qp_reduced_hessian_module`; `test_qp_reduced_hessian`, `test_qp_fuzz`, `test_sparse_options` |
| QP-008 | `sqpopt_qp_daqp` shall solve convex QPs with the DAQP dual active-set solver, and solve with the dense QP solver every QP that DAQP does not solve (nonconvex, inconsistent, or a forced elastic re-solve), counting these in `n_daqp_fallbacks`. | T | `sqpopt_qp_daqp_module`; `test_qp_daqp` |
| QP-009 | With `options%direct_qp`, each QP shall first be tried by sparse factorizations of the KKT matrix of its working set, starting from the previous QP's working set; the active-set solver shall run only if that does not give the solution within `qp_solver%direct_max_changes` (default 10) changes of the working set. | T | `sqpopt_qp_direct_module`; `test_direct` |

### 3.8 Hessian of the Lagrangian (HES)

| ID | Requirement | V | Trace |
|---|---|---|---|
| HES-001 | `options%hessian_mode` shall select a limited-memory BFGS approximation (`sqpopt_hessian_bfgs`, default), a limited-memory SR1 approximation (`sqpopt_hessian_sr1`), or the user's exact sparse Hessian (`sqpopt_hessian_exact`). | T | `sqpopt_hessian_module`; `test_hessian_consistency`, `test_hs_suite --hessian=...` |
| HES-002 | The limited-memory approximations shall store only the last \( k \) pairs \( (s, y) \), where \( k \) is `options%lbfgs_memory`, or, if that is 0 (default), \( \max(10, \min(n, 100)) \), or 10 with `direct_qp`. Hessian-vector products shall cost \( O(nk) \) without forming an \( n \times n \) matrix. | T | `lbfgs_memory`, `hessian_vector_product`; `test_hessian_consistency` |
| HES-003 | The BFGS update shall use Powell's damping when \( s^T y \) is small or negative (`hessian%damping`, default on), and skip ill-conditioned updates. The SR1 update shall be skipped when the standard SR1 safeguard fails. | T | `hessian_update_bfgs`, `hessian_update_sr1`; `test_hessian_consistency` |
| HES-004 | The initial approximation shall be `options%hessian_scale0` \( \times I \), replaced by the usual quasi-Newton scaling after the first update. | I | `sqpopt_hessian_module` |
| HES-005 | With the exact Hessian, a shift \( \delta I \) shall correct indefiniteness: increased tenfold (from `shift_min` times the Hessian's size) where a quasi-Newton approximation would be reset, divided by 3 after each good step, and set to 0 below its minimum. | T | `hessian_reset`; `test_hs_suite --hessian=exact` |
| HES-006 | With `options%inertia_control` and the exact or SR1 Hessian, the shift shall be found from the inertia of the KKT matrix of the QP's working set, by a sparse factorization. It shall not be used with BFGS. | T | `sqpopt_inertia_module`; `test_inertia` |
| HES-007 | With `hessian%convexify` and SR1, an indefinite matrix shall be shifted by its smallest eigenvalue before each QP, from a spectral decomposition of its compact representation that forms no \( n \times n \) or \( n \times k \) array. | T | `sqpopt_spectral_module`; `test_sr1_convexify`, `test_spectrum` |
| HES-008 | The small eigenproblems shall be solved by `hessian%eigen_solver`: QL (any precision), Jacobi, or LAPACK's `DSYEV` (LAPACK build only); `sqpopt_eigen_auto` (default) shall select LAPACK where available, else QL. | T | `sqpopt_eigen_module`; `test_spectrum` |

### 3.9 Globalization (GLB)

| ID | Requirement | V | Trace |
|---|---|---|---|
| GLB-001 | `options%linesearch_mode` shall select one of: a backtracking Armijo search on a merit function (`armijo`); an approximate exact minimization of the merit function by `fmin` (`exact`); Powell's watchdog technique (`watchdog`); the filter method of Fletcher and Leyffer with the rules of Wächter and Biegler (`filter`, default); or the funnel method of Kiessling, Leyffer and Vanaret (`funnel`). | T | `sqpopt_linesearch_module`; `test_acceptance`, `test_hs_suite --linesearch=...` |
| GLB-002 | The filter shall judge a trial point by its \( \ell_1 \) constraint violation and objective, require it not to be dominated by the filter, apply a switching condition between an Armijo (f-type) test and a violation-reduction test, and add points to the filter only on non-f-type steps. | T | `sqpopt_filter_module`; `test_acceptance` |
| GLB-003 | The funnel shall bound the \( \ell_1 \) violation by a width that shrinks on h-type steps, with the same switching between f-type and h-type steps. | T | `sqpopt_funnel_module`; `test_acceptance` |
| GLB-004 | `options%merit_mode` shall select the \( \ell_1 \) exact penalty (default) or a smooth augmented Lagrangian merit function, for the merit-based searches and the trust region's ratio test; `options%penalty_update` shall select the penalty rule (above the multipliers, never decreasing, by default; or each merit's model-based rule). | T | `sqpopt_merit_module`; `test_acceptance`, `test_hs_suite --merit=...` |
| GLB-005 | When a full step is rejected, a second-order correction shall be tried: a minimum-norm correction for the constraints active in the linearization or violated at the trial point, reusing the Jacobian at the current point, projected onto the bounds, and discarded if longer than the step. | T | `sqpopt_soc_module`; `test_maratos` |
| GLB-006 | With `trust_region%enabled` (default off), the step shall be found by re-solving the QP within an \( \ell_\infty \) trust region intersected with the bounds, accepted by the filter or funnel test (with those line-search modes) or by a ratio of actual to predicted merit decrease, with the radius shrunk on rejection (at most `max_retries` per iteration) and grown on a full accepted step. | T | `sqpopt_trust_region_module`; `test_hs_suite --trust-region` |
| GLB-007 | When a line search fails, it shall be retried non-monotonically if `linesearch%nonmonotone_len > 0`. | T | `sqpopt_linesearch_module`; `test_hs_suite --nonmonotone=N` |
| GLB-008 | A trial point with a non-finite \( f \) or \( c \), or one that could not be evaluated, shall never be accepted. | T | `sqpopt_linesearch_module`, `sqpopt_trust_region_module`; `test_nonfinite` |

### 3.10 Feasibility restoration (RST)

| ID | Requirement | V | Trace |
|---|---|---|---|
| RST-001 | When the QP's linearized constraints are inconsistent, the solver shall take a Gauss-Newton step on the constraint violation subject to the bounds, falling back to the QP's elastic step. | T | `restoration_step`; `test_overdetermined`, `test_infeasible` |
| RST-002 | With `options%restoration_mode = sqpopt_restoration_phase` (default), when the filter or funnel search or the trust region finds no acceptable step at an infeasible point, the solver shall enter a restoration phase that minimizes the violation by feasibility QPs with a proximal term, ending when the violation is below `restoration_exit_factor` (default 0.9) times its starting value and the point is acceptable to the filter or funnel, or the point is feasible, or after `restoration_max_iter` (default 50) iterations, or when a phase step fails. The starting point of the phase shall be added to the filter (or the funnel tightened toward it). | T | `sqpopt_restoration_module`; `test_hs_suite --restoration=...` |
| RST-003 | With `sqpopt_restoration_gauss_newton`, each restoration shall be a single Gauss-Newton step. | T | `sqpopt_restoration_module`; `test_hs_suite --restoration=gauss-newton` |

### 3.11 Scaling (SCL)

| ID | Requirement | V | Trace |
|---|---|---|---|
| SCL-001 | With `options%scaling` (default on), the objective and each constraint whose gradient at the starting point has an element larger than `scaling_max_gradient` (default 100) shall be scaled down so that its largest element is that value, but by no factor smaller than `scaling_min_value` (default \( 10^{-8} \)); the finite constraint bounds shall be scaled to match. | T | `compute_scaling`; `test_scaling`, `test_results` |
| SCL-002 | Every reported value (objective, constraints, multipliers, stationarity error, values passed to `report`, `results`) shall be of the original, unscaled problem, except `kkt_error`, which shall be of the scaled problem. | T | `finish`; `test_results` |
| SCL-003 | Starting multipliers `lambda0` shall be converted to the scaled problem. | T | `sqpopt_solve`; `test_results` (warm start) |

### 3.12 Sparse factorizations (FAC)

| ID | Requirement | V | Trace |
|---|---|---|---|
| FAC-001 | The factorization-based options (`inertia_control`, `direct_qp`, `direct_least_squares`) shall use the linear solver of `options%linear_solver`: QDLDL (`sqpopt_linear_solver_qdldl`, default, always available), MUMPS (MUMPS build), SQPOPT's dense \( LDL^T \) (`sqpopt_linear_solver_dense`), or LAPACK's `DSYTRF` (`sqpopt_linear_solver_lapack`, LAPACK build). | T | `sqpopt_symmetric_solver_module`; `test_kkt`, `test_direct`, `test_inertia` |
| FAC-002 | Each backend shall report the inertia it can determine; QDLDL, which does not pivot, shall report a zero pivot that is not a zero eigenvalue as negative curvature. | T | `sqpopt_qdldl_ldl_module`, `sqpopt_kkt_module`; `test_kkt` |
| FAC-003 | With `options%direct_least_squares`, the Gauss-Newton restoration steps and the second-order corrections shall be computed by a sparse factorization instead of LSQR. | T | `sqpopt_least_squares_module`; `test_direct` |
| FAC-004 | `options%factorization_threads` shall set the number of OpenMP threads of MUMPS's factorizations (1 by default; 0: the OpenMP environment's setting). | I | `sqpopt_mumps_ldl_module` |
| FAC-005 | If a sparse solver cannot be started, the solve shall continue without the option that needs it. | I | `sqpopt_solve` |

### 3.13 Results (RES)

| ID | Requirement | V | Trace |
|---|---|---|---|
| RES-001 | Each solve shall produce an `sqpopt_results_type`, reset at the start of the solve, with the status code and its message. | T | `sqpopt_types_module`, `finish`; `test_results` |
| RES-002 | The results shall include the final point `x`, the constraint values `c`, the objective `f`, the constraint multipliers `lambda`, and the variable-bound multipliers `z`. | T | `finish`; `test_results` |
| RES-003 | The results shall include the KKT error (scaled problem), the feasibility error (largest violation of a constraint or bound, original problem), and the unscaled stationarity error. | T | `finish`; `test_results` |
| RES-004 | The results shall include the counts of major iterations, of calls of `fc`, `gjac`, and `hess`, of QP iterations, QP solves, direct QP solves, unconstrained-step QPs, DAQP fallbacks, sparse factorizations, second-order corrections, Hessian resets, restoration steps and phases, elastic re-solves, and escape steps, and the iteration of the switch to accurate derivatives (0 if none). | T | `finish`, `count_events`; `test_results`, `test_termination` |
| RES-005 | The results shall include the wall-clock time of the solve, and the parts of it spent in the user functions, the QP solver, and the factorizations. | T | `finish`; `test_termination` |
| RES-006 | The multipliers shall be those of the Lagrangian \( f - \lambda^T c - z^T x \): non-negative at a lower bound, non-positive at an upper bound, and \( z_j = 0 \) for a variable not at a bound. | T | `sqpopt_results_type`; `test_results` |
| RES-007 | With `options%diagnostic_level >= 1`, the results shall include the diagnosis (section 3.15) in `results%diagnosis`. | T | `sqpopt_diagnostics_module`; `test_diagnostics` |

### 3.14 Status codes (STS)

| ID | Requirement | V | Trace |
|---|---|---|---|
| STS-001 | The solver shall return exactly the status codes of Table A.1, grouped by value: 0–2 success, 10–19 limits and user stop, 20–29 failures. | I | `sqpopt_types_module` |
| STS-002 | The status shall be returned in both `istat` and `results%istat`, and `results%message` shall be `sqpopt_status_message(istat)`, followed for invalid input by a description of what is invalid. | T | `finish`; `test_input_validation`, `test_termination` |
| STS-003 | `sqpopt_status_message` shall return "unknown status code" for a value that is not a status code. | I | `sqpopt_status_message` |

### 3.15 Printed output (OUT)

| ID | Requirement | V | Trace |
|---|---|---|---|
| OUT-001 | With `options%print_level = 0` (default), the solver shall write nothing to `output_unit`, at any diagnostic level (the iteration history of DIA-003 is not printed output). | T | `sqpopt_solve`; `test_termination`, `test_diagnostics` |
| OUT-002 | All printed output shall go to `options%output_unit` (default standard output). | T | `sqpopt_solve`; `test_termination` |
| OUT-003 | With `print_level >= 1`, the solver shall print a header (problem size, method, scaling, step limits, tolerances), one line per major iteration (objective, infeasibility, KKT error, step length, `fc` calls, and flags), and a summary (status, objective, errors, multipliers, active set, iterations, evaluations, events, times) with a legend of the columns and flags. | T | `print_header`, `print_iteration`, `print_summary`; `test_termination` |
| OUT-004 | With `print_level >= 2`, each iteration line shall also have the step norm, QP iterations, evaluations in the iteration, the largest multiplier, the unscaled stationarity error, the globalization's state, and the Hessian's measure. | I | `print_iteration` |
| OUT-005 | With `print_level >= 3`, the solver shall also print the details of each iteration (QP solves, line-search and trust-region trials, restoration phases, Hessian resets) after that iteration's line, the smallest constraint scale factors, and the solution's variables and constraints with their bounds, multipliers, and activity. | T | `sqpopt_log_module`, `print_details`, `print_solution`; `test_termination` |
| OUT-006 | The iteration flags shall be: R restoration step, P restoration phase, S second-order correction, H Hessian reset, E elastic QP re-solve, X escape step, N non-monotone step, W watchdog relaxed step, D switched to accurate derivatives, Q QP failed, F no acceptable step. | I | `sqpopt_iteration_flags`, `print_legend` |
| OUT-007 | No failure of a write shall stop the solver or change its results: every write of the printed output shall use `iostat=`. | I | `sqpopt_module`, `sqpopt_log_module`, `sqpopt_diagnostics_module` |
| OUT-008 | Printing at any level shall not change the iterates or the results. | T | `test_termination` (level 3 against level 0) |

### 3.16 Diagnostics (DIA)

| ID | Requirement | V | Trace |
|---|---|---|---|
| DIA-001 | `options%diagnostic_level` 0 (default) shall compute no diagnostics, at no cost. | T | `sqpopt_diagnostics_module`; `test_diagnostics` |
| DIA-002 | Level 1 shall diagnose the final point: the (at most 5) constraints with the largest violations and multipliers, the variables with the largest stationarity residuals, the numbers of active, weakly active, wrong-sign, and linearly dependent active constraints, the elastic constraints of the last QP, and a verdict for the final status. | T | `diagnostics_finish`; `test_diagnostics` |
| DIA-003 | Level 2 shall also report on the starting point (gradient size ratios, constant and single-variable constraints, dependent equalities) and on the iterations (convergence rate, iterations needed at that rate, active-set flips, constraints whose changes disagree with their derivatives, the slowest iteration); with `diagnostics_unit`, it shall write one comma-separated line per iteration, with a header line, whatever the `print_level`. | T | `diagnostics_start`, `diagnostics_record`; `test_diagnostics` |
| DIA-004 | Level 3 shall also evaluate `fc` at two points next to the starting point, and report whether either value is not finite. | T | `sqpopt_diagnostics_module`; `test_diagnostics` |
| DIA-005 | Levels 0 to 2 shall call the user functions exactly as often as level 0, and no level shall change the iterates. | T | `test_diagnostics`; `test_hs_suite --diagnostics=2` (same counts as the baseline) |
| DIA-006 | The diagnostics shall be returned in `sqpopt_diagnosis_type` as values and as the texts `report` and `problem_report`, which are printed only with `print_level >= 1`. | T | `sqpopt_types_module`; `test_diagnostics` |

### 3.17 Least-squares interface (NLS)

| ID | Requirement | V | Trace |
|---|---|---|---|
| NLS-001 | `sqpopt_nlls_type` shall solve \( \min \tfrac12 \sum_i r_i(x)^2 \) s.t. \( c_l \le c(x) \le c_u \), \( x_l \le x \le x_u \), by solving with `sqpopt_type` the equivalent problem in \( (x, z) \) with the equality constraints \( r(x) - z = 0 \) (Schittkowski's form). | T | `sqpopt_nlls_module`; `test_nlls` |
| NLS-002 | Its `initialize` shall take the sizes, a residual function `residuals(x, r, c, status, data)`, a Jacobian function `jacobian(x, rjac_val, cjac_val, accuracy, status, data)`, the residual Jacobian's COO pattern, optional bounds and constraint pattern, optional solver components, and optional user data. | T | `nlls_initialize`; `test_nlls` |
| NLS-003 | With `options%hessian_mode = sqpopt_hessian_exact`, the Gauss-Newton Hessian (identity on \( z \), zero on \( x \)) shall be used. | T | `nlls_hess`; `test_nlls` |
| NLS-004 | `solve` shall start \( z \) from the residuals at the (projected) starting point if the constraints hold there, else from zero, and reuse that evaluation for the solver's first one. | I | `nlls_solve`, `nlls_fc` |
| NLS-005 | `get_solution` shall return the variables and optionally the residuals, the sum of squares, and the constraint multipliers; `get_results` shall return the results of the transformed problem. | T | `nlls_get_solution`; `test_nlls` |
| NLS-006 | Invalid input (sizes, indices out of range, no `initialize`) shall be reported by `solve` as `sqpopt_invalid_input`. | T | `nlls_solve`; `test_nlls` |
| NLS-007 | The solver shall not depend on the least-squares interface: no module other than `sqpopt_nlls_module` shall refer to it. | I | `src/` |

### 3.18 Python bindings (PY)

| ID | Requirement | V | Trace |
|---|---|---|---|
| PY-001 | The package `sqpopt` shall provide `minimize(fun, x0, args, jac, hess, bounds, constraints, tol, callback, options, lambda0, output_file, max_step, diagnostics_file)` following the interface of `scipy.optimize.minimize`, returning an `OptimizeResult`. | T | `python/sqpopt/_minimize.py`; `test_minimize`, `test_scipy_compare` |
| PY-002 | `minimize` shall require the objective's gradient and the constraints' Jacobians (no finite differences), accept dense or sparse Jacobians and Hessians, and accept constraints as SLSQP-style dicts, `NonlinearConstraint`, or `LinearConstraint` (SciPy's or its own). | T | `_minimize.py`, `_constraints.py`; `test_minimize` |
| PY-003 | `minimize` shall use the exact Hessian of the Lagrangian when `hess` is given (each nonlinear constraint then needing its `hess(x, v)`), and a quasi-Newton approximation otherwise. | T | `_minimize.py`; `test_minimize` |
| PY-004 | Every Fortran option shall be settable through `options`, by its Fortran reference, its unambiguous name, or a nested dict; choice options by the name of their constant; and `maxiter` and `disp` as in SciPy. | T | `_minimize.py`, `schema.py`; `test_minimize`, `test_schema` |
| PY-005 | The result's `success` shall be true exactly when `status <= 2`; `status` and `message` shall be the solver's. | T | `_minimize.py`; `test_minimize` |
| PY-006 | The callback shall be called once per major iteration, as `callback(intermediate_result)` if its argument is so named, else `callback(x)`; returning `True` or raising `StopIteration` shall stop the solve. | T | `_minimize.py`; `test_minimize` |
| PY-007 | `least_squares(fun, x0, jac, bounds, constraints, args, jac_sparsity, gauss_newton, tol, options, output_file, max_step, diagnostics_file)` shall minimize \( \tfrac12 \sum \) `fun(x)**2` on top of `minimize`. | T | `_least_squares.py`; `test_least_squares` |
| PY-008 | `minimize` shall refuse to run with an extension built from a different options schema than the installed one. | I | `_native`, `_build.options_signature` |
| PY-009 | The schema (`sqpopt_options.schema`) shall describe every Fortran option with its default and limits, and shall be importable without Qt; `python/tests/test_schema.py` shall fail if a Fortran field and the schema disagree. | T | `schema.py`; `test_schema` |
| PY-010 | The options dialog (`SqpoptOptionsDialog`, `edit_options`) shall edit every option of the schema and return the settings as a nested dict with the Fortran names. | T | `dialog.py`; `test_dialog` |

## 4. Non-functional requirements

### 4.1 Memory and stack (MEM)

| ID | Requirement | V | Trace |
|---|---|---|---|
| MEM-001 | The library shall have no automatic arrays: every local array whose size depends on the problem shall be `allocatable` and allocated at run time. | D | `tools/stack_check.sh` (fails if one is found) |
| MEM-002 | The library shall not use the `matmul` intrinsic; dense products shall be computed by `dense_matvec` and `dense_matmul` into arrays allocated by the caller. | D | `sqpopt_dense_product_module`; `tools/stack_check.sh` (fails on a use) |
| MEM-003 | Outside the dense QP solver and its linear algebra, no array temporary whose size grows with the problem shall be created (array expressions as actual arguments, constructors, or array-valued function results of that size shall be assigned to work arrays first). | I | `tools/stack_check.sh` (lists the temporaries for review) |
| MEM-004 | Work arrays of a routine called in an inner loop shall be kept between calls (in a type or by the host procedure), so that the call allocates nothing. | I | `src/` |
| MEM-005 | No dense \( n \times n \) or \( m \times n \) array shall be formed, except in the dense QP solver, the DAQP mode, and the dense linear solvers (`sqpopt_linear_solver_dense`, `sqpopt_linear_solver_lapack`, limited to an order of 2,000). | I | `src/` |

### 4.2 State and reentrancy (STA)

| ID | Requirement | V | Trace |
|---|---|---|---|
| STA-001 | No local variable shall be initialized in its declaration (which would give it the `save` attribute). | I | `src/` (review rule, `CLAUDE.md`) |
| STA-002 | The objects that hold a sparse solver (`sqpopt_kkt_type`, `sqpopt_least_squares_type`) shall live for one solve and be freed at its end. | I | `sqpopt_solve`, `finish` |

### 4.3 Portability and precision (PRT)

| ID | Requirement | V | Trace |
|---|---|---|---|
| PRT-001 | The working precision shall be double (`real64`) by default, and selectable as single (`REAL32`) or quadruple (`REAL128`) precision by a preprocessor directive. | I | `sqpopt_kinds.F90` |
| PRT-002 | The MUMPS and LAPACK builds shall be double precision only; defining `REAL32` or `REAL128` with `HAS_MUMPS` or `HAS_LAPACK` shall stop the compilation with an error. | I | `#error` in `sqpopt_mumps_ldl_module.F90`, `sqpopt_dense_ldl_module.F90`, `sqpopt_eigen_module.F90`, `sqpopt_dense_product_module.F90` |
| PRT-003 | Every program unit shall have `implicit none`; reals shall use the kind `wp` from `sqpopt_kinds`. | D | `pixi run fortitude check` |

### 4.4 Build and dependencies (BLD)

| ID | Requirement | V | Trace |
|---|---|---|---|
| BLD-001 | The default build shall need nothing but fpm and its fetched dependencies (section 2.4). | D | CI `Build` job (`pixi run fpm test`) |
| BLD-002 | MUMPS shall be referenced only inside `#ifdef HAS_MUMPS`, in `sqpopt_mumps_ldl_module.F90`. `sqpopt_has_mumps` shall say whether the library was built with it. | I | `sqpopt_mumps_ldl_module.F90` |
| BLD-003 | LAPACK and BLAS shall be referenced only inside `#ifdef HAS_LAPACK`, in `sqpopt_dense_ldl_module.F90`, `sqpopt_eigen_module.F90`, and `sqpopt_dense_product_module.F90`. `sqpopt_has_lapack` shall say whether the library was built with them. | I | those modules |
| BLD-004 | In the LAPACK build, `dense_matvec` and `dense_matmul` shall call BLAS's `DGEMV` and `DGEMM` for non-empty products of contiguous arrays, and their own loops otherwise. | T | `sqpopt_dense_product_module`; `test_dense_product` (run with `pixi run test-lapack`) |
| BLD-005 | Every test shall pass in the default build, the MUMPS build, and the LAPACK build. | D | `pixi run fpm test`, `pixi run test-mumps`, `pixi run test-lapack` |
| BLD-006 | `fpm install` shall install the static library and module files, from which the Python bindings' extension is built (`pixi run build-python`). | D | `fpm.toml`, `python/sqpopt/_build.py` |

### 4.5 Robustness (ROB)

| ID | Requirement | V | Trace |
|---|---|---|---|
| ROB-001 | The solver shall not crash on invalid input (it shall report `sqpopt_invalid_input`), on non-finite function values (it shall reject the point or report `sqpopt_function_error`), or on a failure to allocate the dense QP's matrices (it shall report `sqpopt_out_of_memory`). | T | `test_input_validation`, `test_nonfinite`, `test_out_of_memory` |
| ROB-002 | With floating-point trapping enabled for overflow, invalid operations, and division by zero, infinite bounds given as \( \pm \)`huge` shall cause no exception. | T | `test_infinite_bounds` |
| ROB-003 | The QP solvers shall return a correct solution, or a failure status, on randomized small QPs. | T | `test_qp_fuzz` |

### 4.6 Performance and regression (PRF)

| ID | Requirement | V | Trace |
|---|---|---|---|
| PRF-001 | With the default options, every HS-suite problem that is not in the `known_unsolved` list of `test_hs_suite.f90` (25 problems at the baseline) shall be solved: feasible to \( 10^{-6} \), with an objective within \( 10^{-4} \) (relative) of the validated optimum or better. | T | `test_hs_suite` (fails on a regression) |
| PRF-002 | The HS suite's reference solutions shall be valid: at each problem's recorded optimum, the objective shall match and the constraints shall hold. | T | `test_hs_solutions` |
| PRF-003 | The solver shall solve large sparse problems (with the sparse QP solver) with the default options and with each combination of the options that matter for the sparse QP. | T | `test_large_sparse`, `test_sparse_options`, `test_scalable` |
| PRF-004 | A change that can affect convergence shall be compared with the HS baseline, in a release build, before it is accepted; new failures shall be investigated, not added to the baseline. | I | development process (`CLAUDE.md`) |

### 4.7 Code standards and documentation (STD)

| ID | Requirement | V | Trace |
|---|---|---|---|
| STD-001 | The Fortran sources of `src/`, `example/`, and the Python bindings shall pass `fortitude check` with the repository's `fortitude.toml` (4-space indentation, lines of at most 132 columns, single-quoted strings). | D | CI `Lint` job |
| STD-002 | Every module shall have a `!>` header; every procedure, type, and type component a `!!` docstring; and every dummy argument (except the passed object) a trailing `!!` comment. | I | FORD documentation build (CI) |
| STD-003 | Status codes shall be referred to by their named constants, never by their numeric values. | I | `src/`, `test/` |
| STD-004 | The API documentation shall be generated by FORD from the source, and published with the user guide and the coverage report by CI on the `master` branch. | D | CI `Build` job |

## Appendix A. Reference tables

### A.1 Status codes

| Constant | Value | Class | Meaning |
|---|--:|---|---|
| `sqpopt_success` | 0 | success | converged successfully |
| `sqpopt_acceptable` | 1 | success | converged to the acceptable level (CNV-005) |
| `sqpopt_stalled` | 2 | success | feasible, but no further progress (CNV-006) |
| `sqpopt_max_iter_reached` | 10 | limit | `max_iter` major iterations |
| `sqpopt_user_requested_stop` | 11 | limit | `report` or a user function asked to stop |
| `sqpopt_max_evals_reached` | 12 | limit | `max_evals` calls of `fc` |
| `sqpopt_time_limit_reached` | 13 | limit | `max_time` seconds |
| `sqpopt_infeasible` | 21 | failure | the problem appears to be (locally) infeasible (CNV-007) |
| `sqpopt_line_search_failed` | 22 | failure | no acceptable step, `max_consecutive_failures` times |
| `sqpopt_qp_solve_failed` | 23 | failure | the QP solver failed, `max_consecutive_failures` times |
| `sqpopt_invalid_input` | 24 | failure | invalid problem definition or options (section 3.4) |
| `sqpopt_function_error` | 25 | failure | non-finite or failed evaluation at the current point (ALG-008) |
| `sqpopt_unbounded` | 26 | failure | the objective fell below `obj_lower_limit` at a feasible point |
| `sqpopt_out_of_memory` | 27 | failure | an array could not be allocated (QP-006) |

### A.2 Options (`sqpopt_options_type`): defaults and valid values

| Option | Default | Valid values |
|---|---|---|
| `max_iter` | 100 | \( \ge 0 \) |
| `print_level` | 0 | 0, 1, 2, 3 (OUT) |
| `output_unit` | standard output | an open unit, if `print_level > 0` |
| `diagnostic_level` | 0 | 0, 1, 2, 3 |
| `diagnostics_unit` | \(-1\) (none) | \(-1\), or an open unit, if `diagnostic_level >= 2` |
| `hessian_mode` | `sqpopt_hessian_bfgs` | `sqpopt_hessian_bfgs`, `_sr1`, `_exact` (the last needs `hess`) |
| `inertia_control` | `.false.` | |
| `direct_qp` | `.false.` | |
| `direct_least_squares` | `.false.` | |
| `factorization_threads` | 1 | \( \ge 0 \) |
| `linear_solver` | `sqpopt_linear_solver_qdldl` | `_qdldl`, `_dense`, `_mumps` (MUMPS build), `_lapack` (LAPACK build); order limit, VAL-007 |
| `lbfgs_memory` | 0 (automatic) | \( \ge 0 \) |
| `qp_solver_mode` | `sqpopt_qp_auto` | `sqpopt_qp_auto`, `_dense`, `_reduced_hessian`, `_daqp` |
| `linesearch_mode` | `sqpopt_linesearch_filter` | `_armijo`, `_exact`, `_watchdog`, `_filter`, `_funnel` |
| `merit_mode` | `sqpopt_merit_l1` | `sqpopt_merit_l1`, `_augmented_lagrangian` |
| `penalty_update` | `sqpopt_penalty_multipliers` | `sqpopt_penalty_multipliers`, `_model` |
| `restoration_mode` | `sqpopt_restoration_phase` | `sqpopt_restoration_phase`, `_gauss_newton` |
| `restoration_exit_factor` | 0.9 | \( (0, 1) \) |
| `restoration_max_iter` | 50 | \( \ge 1 \) |
| `derivative_accuracy` | `sqpopt_derivatives_accurate` | `sqpopt_derivatives_fast`, `_accurate` |
| `derivative_switch_tol` | \( 10^{-5} \) | \( \ge 0 \) |
| `elastic_multiplier_limit` | 30 | \( \ge 0 \) (0: off) |
| `max_consecutive_failures` | 5 | \( \ge 1 \) |
| `ftol`, `xtol` | \( 10^{-8} \) | \( \ge 0 \) |
| `ctol` | \( 10^{-8} \) | \( > 0 \) |
| `ktol` | \( 10^{-6} \) | \( > 0 \) |
| `dual_inf_tol` | 1 | \( > 0 \) |
| `acceptable_ktol` | \( 10^{-4} \) | \( > 0 \) |
| `acceptable_ctol` | \( 10^{-6} \) | \( > 0 \) |
| `acceptable_iter` | 15 | \( \ge 0 \) (0: off) |
| `acceptable_obj_change_tol` | \( 10^{-3} \) | \( \ge 0 \) |
| `stall_iter` | 3 | \( \ge 1 \) |
| `max_evals` | 0 (no limit) | \( \ge 0 \) |
| `max_time` | 0 (no limit) | \( \ge 0 \) seconds |
| `obj_lower_limit` | \( -10^{20} \) (no limit) | any |
| `scaling` | `.true.` | |
| `scaling_max_gradient` | 100 | \( > 0 \) |
| `scaling_min_value` | \( 10^{-8} \) | \( [0, 1] \) |
| `hessian_scale0` | 1 | \( > 0 \) |

### A.3 Component settings: valid values

| Component | Constraint checked by `validate_options` |
|---|---|
| `linesearch` | `sigma`, `backtrack` in \( (0,1) \); `alpha_min` (default \( 10^{-10} \)), `tol`, `major_step_limit` (default 2) \( > 0 \); `max_ls_iter` (default 40) \( \ge 1 \); `merit%penalty`, `watchdog_relaxed_len`, `watchdog_cooldown_len`, `nonmonotone_len` \( \ge 0 \); `merit%penalty_rho` in \( (0,1) \) |
| `linesearch%filter` | `gamma_theta` in \( (0,1) \); `gamma_phi`, `delta` \( > 0 \); `s_theta`, `s_phi` \( > 1 \); `eta_phi` in \( (0, 0.5) \); `0 < theta_min_fact < theta_max_fact`; `gamma_alpha` in \( (0, 1] \) |
| `linesearch%funnel` | `width_min` \( > 0 \); `width_fact` \( \ge 1 \); `beta`, `kappa` in \( (0,1) \); `delta` \( > 0 \); `s_theta` \( > 1 \); `eta` in \( (0, 0.5) \); `gamma` \( > 0 \); `update` 1 or 2 |
| `trust_region` (if `enabled`) | `0 < radius_min <= radius0 <= radius_max` (defaults \( 10^{-8} \), 1, \( 10^3 \)); `0 < eta1 <= eta2 < 1`; `shrink_factor` in \( (0,1) \); `expand_factor` \( \ge 1 \); `max_retries` (default 20) \( \ge 1 \) |
| `qp_solver` | `max_step` (default 2) \( > 0 \); `auto_dense_max_n` (default 200) \( \ge 0 \); `direct_max_changes` (default 10) \( \ge 0 \); `direct_tol` (default \( 10^{-8} \)) \( > 0 \) |
| `qp_solver%dense_qp` | `active_tol`, `opt_tol`, `feas_tol`, `elastic_weight` \( > 0 \); `elastic_weight_max` \( \ge \) `elastic_weight`; `max_iter` \( \ge 1 \) |
| `qp_solver%sparse_qp` | `active_tol`, `opt_tol`, `feas_tol`, `pcg_rtol`, `elastic_weight` \( > 0 \); `elastic_weight_max` \( \ge \) `elastic_weight`; `lsqr_atol`, `lsqr_btol`, `lsqr_conlim` \( \ge 0 \); `max_iter` \( \ge 1 \); `dense_max_ns` \( \ge 0 \); `null_space` `sqpopt_null_space_lu` or `_lsqr` |
| `qp_solver%daqp_qp` | `primal_tol`, `dual_tol` \( > 0 \); `max_iter` \( \ge 1 \) |
| `hessian` | `eigen_solver` `sqpopt_eigen_auto`, `_jacobi`, `_ql`, or `_lapack` (LAPACK build) |

## Appendix B. Known limitations at the baseline

These are behaviors of the baseline that a future requirement may change (see `plan/ROADMAP.md`):

* The stalled-progress test (CNV-006) measures the change of the variables relative to \( \lVert x \rVert_2 \),
  so on a badly scaled problem it can stop the solve while the iteration is recovering from a quasi-Newton
  pair of large curvature, far from a solution (ROADMAP F34, case 1; `test_nlls` sets `ftol = 0` for that
  problem).
* A crawl at a constant step length just above `xtol` is not detected as a stall (ROADMAP F34, case 2).
* 25 HS problems end at a local solution with the default options (PRF-001).
