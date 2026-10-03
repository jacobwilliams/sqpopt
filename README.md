![Modern Fortran SQP OPTimizer](media/logo-small.png)

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![Build Status](https://github.com/jacobwilliams/sqpopt/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/sqpopt/actions)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/sqpopt)](https://github.com/jacobwilliams/sqpopt/commits/master)
[![Docs](https://img.shields.io/badge/docs-user%20guide-blue)](https://jacobwilliams.github.io/sqpopt/)


Modern Fortran **SQP** **OPT**imizer: a modular sequential quadratic
programming solver for large, sparse nonlinear optimization problems. A work
in progress.

**Using SQPOPT?** See the **[User Guide](https://jacobwilliams.github.io/sqpopt/)**.
It covers installation as an FPM dependency, the API, every option, the
status codes, worked examples, and benchmark results. The
[API documentation](https://jacobwilliams.github.io/sqpopt/api/),
[test coverage](https://jacobwilliams.github.io/sqpopt/coverage/), and
[Hock-Schittkowski results](https://jacobwilliams.github.io/sqpopt/performance.html#hs-results)
are published there too.

This README is for working on SQPOPT itself.

## Development environment

Everything (gfortran, fpm, lcov, FORD, Python with Qt, ...) comes from the
[pixi](https://pixi.sh) environment in `pixi.toml`, which CI uses too. Run
commands with `pixi run ...`, or open a shell in the environment with
`pixi shell`. (On macOS, use the pixi gfortran: a system gfortran may fail
to link.)

## Building and testing

```sh
pixi run fpm build                     # debug build of the library
pixi run fpm build --profile release
pixi run fpm test                      # all the tests
pixi run fpm test test_hs71            # one test
```

The default precision is double. For single or quadruple precision, define
`REAL32` or `REAL128` (see `src/sqpopt_kinds.F90`):

```sh
pixi run fpm test --flag "-DREAL128"
```

### Sparse factorizations: QDLDL, and MUMPS (optional)

Three options, meant for large problems, use a sparse LDLᵀ factorization,
by the solver of `options%linear_solver`: [QDLDL](https://github.com/jacobwilliams/qdldl-fortran)
(the default: no pivoting, one thread, very little overhead; an fpm
dependency, so always available) or [MUMPS](https://mumps-solver.org)
(with pivoting and threads; optional, and chosen with
`options%linear_solver = sqpopt_linear_solver_mumps`):

| option | what it does | module |
|---|---|---|
| `options%inertia_control` | finds the shift of an indefinite Hessian (exact or SR1) from the inertia of the KKT matrix | `sqpopt_inertia_module` |
| `options%direct_qp` | solves the QP subproblems directly, by factoring the KKT matrix of the working set | `sqpopt_qp_direct_module` |
| `options%direct_least_squares` | computes restoration steps and second-order corrections directly instead of with `LSQR` | `sqpopt_least_squares_module` |

All three are built on `sqpopt_kkt_module` (the KKT matrix of a working
set) and `sqpopt_symmetric_solver_module`, which holds both solvers, and is
the only source file that refers to MUMPS, and only inside
`#ifdef HAS_MUMPS`. So the default build still needs nothing but fpm, and
has all three options, with QDLDL. QDLDL is exact for the least-squares
systems and the quasi-Newton Hessians, and much faster than MUMPS on banded
and chained problems. Choose MUMPS for the exact Hessian with inertia
control, and for problems coupled in two or three dimensions (see the
guide's "Sparse solver" section, which compares them;
`example/sparse_solvers.f90` times the two solvers on grid matrices).
`options%factorization_threads` sets the number of
OpenMP threads MUMPS uses (1 by default; conda-forge's `mumps-seq` is built
with OpenMP). With `HAS_MUMPS`, the library must be compiled in
double precision (the default: `REAL32` and `REAL128` are a compile error
with it). The pixi environment has the sequential MUMPS library
(conda-forge's `mumps-seq`), and tasks that build with it:

```sh
pixi run build-mumps                   # fpm build, with MUMPS
pixi run test-mumps                    # fpm test, with MUMPS
pixi run test-mumps test_kkt           # (other arguments go to fpm)
pixi run test-mumps test_hs_suite --profile release -- results.md --hessian=exact --inertia --direct
pixi run run-mumps --example benchmark_large --profile release -- --scale=10
```

They run fpm with
`--flag "-DHAS_MUMPS -I$CONDA_PREFIX/include" --link-flag "-ldmumps_seq"`:
the include path is for MUMPS's `dmumps_struc.h`. A program that uses a
library built this way must link with `-ldmumps_seq` too. Without
`HAS_MUMPS`, `sqpopt_has_mumps` is false and only
`linear_solver = sqpopt_linear_solver_mumps` is rejected as invalid input.
Changes to this code must be tested in both builds: `test_kkt` checks each
solver of the build, and `test_inertia`, `test_direct`, and `test_qp_fuzz`
the features with each solver of the build (`test_qp_fuzz` with the default,
QDLDL), and, without MUMPS, the rejection of MUMPS. The HS suite's,
`test_scalable`'s, and `benchmark_large`'s `--linear-solver=qdldl|mumps`
compare the two.

### Tests

The tests are in `test/`:

- **Unit tests** of the components on their own: `test_qp_dense`,
  `test_qp_reduced_hessian`, and `test_qp_fuzz` (both QP solvers on thousands
  of random convex, nonconvex, degenerate, and infeasible QPs, checked against
  the KKT conditions), `test_hessian_consistency`, `test_acceptance` (merit
  function, filter, funnel), `test_convergence`,
  `test_independent_columns`, and, for the sparse factorizations (with
  each solver of the build), `test_kkt` (the sparse solver, the KKT matrix against a dense
  reference, and the least-squares solver) and `test_inertia` (the shift's
  search). `test_qp_fuzz` gives its QPs to the direct method too.
- **Solver tests** on small problems with known solutions (`test_basic`,
  `test_hs71`, `test_medium`, `test_maratos`, `test_degenerate` (a
  constraint tangent to a bound), `test_direct` (the options that factor
  matrices, with every Hessian mode), and `test_multipliers` (the
  least-squares multiplier estimate, and a hanging chain that the exact
  Hessian only solves with it), larger sparse ones
  (`test_large_sparse`), a two-variable problem whose path is drawn in the
  guide (`test_rosenbrock_disk`), 17 bound-constrained functions of any size
  (`test_scalable`), tests of single features (`test_unconstrained_step`,
  the QP front end's shortcut for a QP with nothing active;
  `test_max_step`, the limits on each variable's step; `test_scaling`, the
  safeguards for a badly scaled starting point), and regression tests of the interface and edge cases
  (`test_callbacks`, `test_input_validation`, `test_infeasible`,
  `test_nonfinite`, `test_resolve`, `test_results`, `test_termination`, ...).
- **`test_hs_suite`**: the 305 Hock-Schittkowski problems
  (`test/schittkowski_problems.f90`). It is also the regression baseline:
  the test fails if a problem not listed in `known_unsolved` isn't solved,
  and reports problems in that list that now are. Its command-line options
  select other configurations (these skip the regression check):

  ```sh
  pixi run fpm test test_hs_suite --profile release -- results.md --linesearch=funnel
  pixi run fpm test test_hs_suite -- /dev/null --problem=71 --print=2   # one problem, with its iterations
  ```

  The options are `--linesearch=`, `--merit=`, `--penalty=`, `--hessian=`,
  `--inertia`, `--direct`, `--direct-ls`, `--qp=`, `--restoration=`,
  `--trust-region`, `--no-interpolate`,
  `--nonmonotone=N`, `--lbfgs-memory=N`, `--problem=N`, `--print=L`, and
  `--web-data=FILE` (see the header of `test/test_hs_suite.f90`).
- **`test_hs_solutions`** checks the collection's reference solutions: at
  each problem's recorded optimal point, the objective must equal the
  recorded optimal value and the constraints must hold, up to the rounding
  of the recorded digits. It guards the problem definitions, and lists the
  collection's own inconsistencies (`known_inconsistent`, with the reason
  for each).
- **`test_hs_slsqp`** runs the same problems with
  [SLSQP](https://github.com/jacobwilliams/slsqp) (a dev-dependency), for
  the comparison on the results page.

### Benchmark

`example/benchmark.f90` measures function evaluations and run time on a
discretized optimal-control problem and a constrained chained-Rosenbrock
problem, at sizes that use each QP solver:

```sh
pixi run fpm run --example benchmark --profile release
```

`example/benchmark_large.f90` solves larger problems with analytic second
derivatives (those two, a nonconvex chain of double wells, and a chain of
circle constraints that needs second-order corrections) with each
Hessian mode and with the options that use sparse factorizations (with
`--linear-solver=` choosing their solver). It reports where the time goes (the QP solver, the
factorizations) and how many QPs were solved directly. `--scale=S`
multiplies the sizes (`S = 100` gives a million variables; see the header of
the file for the other options):

```sh
pixi run fpm run --example benchmark_large --profile release
pixi run run-mumps --example benchmark_large --profile release -- --scale=10 --no-active-set
```

## Tools

| | |
|---|---|
| `coverage.sh [--mumps]` | runs the tests with `--coverage` and makes an lcov HTML report in `coverage/html` (dark-mode aware). With `--mumps` (what CI uses), the build with MUMPS is tested, so the report covers the code that uses it |
| `pixi run fortitude check` | lints the Fortran sources with [Fortitude](https://fortitude.readthedocs.io) (configured in `fortitude.toml`; also run by the VS Code extension) |
| `tools/hs_performance_table.sh --mumps` | runs the HS suite in every configuration of the Performance page's table (`web/performance.html`) and prints the table rows. It also regenerates the data of that page's results of every problem (`web/js/hs_results_data.js`, `web/js/hs_slsqp_data.js`). Run it whenever a change affects the HS results. (`--mumps` builds with MUMPS, for the rows of the options that need it; without it those rows are left out.) |
| `tools/rosenbrock_disk_figure.py` | redraws the figure of the Performance page's first example, from the iterates that `test_rosenbrock_disk` writes with `--path=FILE` (the script's header has both commands) |
| `tools/hs_compare.sh N [options]` | runs HS problem `N` with SQPOPT and SLSQP, printing both solvers' iterations, for investigating a difference |
| `python/` | Python bindings with a `scipy.optimize.minimize`-like interface, and a Qt options dialog for SQPOPT (see [python/README.md](python/README.md)) |

All are run from the repository root, e.g. `pixi run tools/hs_compare.sh 220`.

## Documentation and website

The website is `web/`: the user guide (`web/index.html`), the Performance
page (`web/performance.html`: the Hock-Schittkowski results, the large
problems, and the scalable test functions), the settings page
(`web/choosing_settings.html`), and their CSS and JavaScript. Update the
guide along with any change to the API, options, or behavior. When the HS
results change, regenerate the Performance page's table with
`tools/hs_performance_table.sh`.

The API documentation is generated from the source comments with
[FORD](https://github.com/Fortran-FOSS-Programmers/ford):

```sh
pixi run ford ford.md --output_dir web/api
```

The design documents and the roadmap are in [plan/](plan/): the
architecture ([PLAN.md](plan/PLAN.md)), the backlog
([ROADMAP.md](plan/ROADMAP.md)), and the latest code review
([CODE_REVIEW.md](plan/CODE_REVIEW.md)).

## Source layout

| module | |
|---|---|
| `sqpopt_module` | the solver object: `initialize`, `solve`, results, and the validation of the inputs |
| `sqpopt_iterate_module` | the major iterations |
| `sqpopt_problem_module` | the problem definition, user-function interface, evaluation caching, and scaling |
| `sqpopt_options_module` | the solver options |
| `sqpopt_types_module` | status codes, the sparse matrix and results types, and small utilities |
| `sqpopt_hessian_module` | the limited-memory BFGS/SR1 approximations, and the exact Hessian |
| `sqpopt_symmetric_solver_module` | the sparse symmetric indefinite solver: chooses a backend (`options%linear_solver`), and refines the solves |
| `sqpopt_sparse_ldl_module` | the abstract interface of a sparse LDLᵀ backend |
| `sqpopt_qdldl_ldl_module` | the QDLDL backend (the default; no pivoting) |
| `sqpopt_mumps_ldl_module` | the MUMPS backend (only in a build with `HAS_MUMPS`) |
| `sqpopt_kkt_module` | the KKT matrix of a QP working set, factored with that solver (its inertia, and solves) |
| `sqpopt_inertia_module` | inertia control of the exact and SR1 Hessians |
| `sqpopt_qp_direct_module` | the direct QP method (a primal-dual active-set method on the KKT matrix) |
| `sqpopt_least_squares_module` | direct minimum-norm solves, for the restoration steps and second-order corrections |
| `sqpopt_qp_solver_module` | the QP subproblem front end, which chooses a solver (and tries the unconstrained quasi-Newton step first) |
| `sqpopt_qp_dense_module` | the dense active-set QP |
| `sqpopt_qp_reduced_hessian_module` | the sparse active-set QP (LUSOL basis, reduced-Hessian CG) |
| `sqpopt_linesearch_module` | the line searches (Armijo, exact, watchdog, filter, funnel) |
| `sqpopt_merit_module`, `sqpopt_filter_module`, `sqpopt_funnel_module` | the acceptance tests the line searches and the trust region use |
| `sqpopt_trust_region_module` | the trust-region globalization |
| `sqpopt_restoration_module` | feasibility restoration |
| `sqpopt_nlls_module` | an optional interface for least-squares problems (`sqpopt_nlls_type`), on top of the solver: it gives the solver the problem in Schittkowski's form, `min 1/2 z'z` subject to `r(x) - z = 0` |
| `sqpopt_soc_module` | second-order corrections |
| `sqpopt_convergence_module` | the KKT convergence test |
| `sqpopt_log_module` | the detailed iteration log, and the formatting of the printed output |
| `sqpopt_linalg_module`, `sqpopt_dense_linalg_module` | sparse and dense linear algebra |
| `sqpopt_kinds` | the real kind (precision) |

## Continuous integration

On every push, [CI](.github/workflows/CI.yml) builds the pixi environment
from the locked `pixi.lock`, runs the tests of the default build, then
those of the build with MUMPS with coverage, and builds the FORD
documentation. On `master`, it deploys `web/` (with the coverage
report and API docs) to GitHub Pages.

## Dependencies

Fetched and built by fpm:

- [LSQR](https://github.com/jacobwilliams/LSQR): iterative sparse least-squares solver
- [lusol](https://github.com/jacobwilliams/lusol): sparse LU factorization (the sparse QP's basis factors and updates)
- [fmin](https://github.com/jacobwilliams/fmin): derivative-free 1-D minimization (the exact line search)
- [qdldl-fortran](https://github.com/jacobwilliams/qdldl-fortran): sparse LDLᵀ factorization without pivoting (the default sparse solver without MUMPS)
- [slsqp](https://github.com/jacobwilliams/slsqp) (tests only): the SLSQP comparison

Optional, from the pixi environment (see "Sparse factorizations"):

- [MUMPS](https://mumps-solver.org): sparse symmetric indefinite factorization (inertia control, the direct QP method, and direct least-squares solves)

## License

[MIT](LICENSE)
