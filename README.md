![Modern Fortran SQP OPTimizer](media/logo-small.png)

Modern Fortran **SQP** **OPT**imizer: a modular sequential quadratic
programming solver for large, sparse nonlinear optimization problems. A work
in progress.

**Using SQPOPT?** See the **[User Guide](https://jacobwilliams.github.io/sqpopt/)**.
It covers installation as an FPM dependency, the API, every option, the
status codes, worked examples, and benchmark results. The
[API documentation](https://jacobwilliams.github.io/sqpopt/api/),
[test coverage](https://jacobwilliams.github.io/sqpopt/coverage/), and
[Hock–Schittkowski results](https://jacobwilliams.github.io/sqpopt/hs_results.html)
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

### Tests

The tests are in `test/`:

- **Unit tests** of the components on their own: `test_qp_dense`,
  `test_qp_reduced_hessian`, and `test_qp_fuzz` (both QP solvers on thousands
  of random convex, nonconvex, degenerate, and infeasible QPs, checked against
  the KKT conditions), `test_hessian_consistency`, `test_acceptance` (merit
  function, filter, funnel), `test_convergence`, and
  `test_independent_columns`.
- **Solver tests** on small problems with known solutions (`test_basic`,
  `test_hs71`, `test_medium`, `test_maratos`), larger sparse ones
  (`test_large_sparse`), and regression tests of the interface and edge cases
  (`test_callbacks`, `test_input_validation`, `test_infeasible`,
  `test_nonfinite`, `test_resolve`, `test_results`, `test_termination`, ...).
- **`test_hs_suite`**: the 305 Hock–Schittkowski problems
  (`test/schittkowski_problems.f90`). It is also the regression baseline:
  the test fails if a problem not listed in `known_unsolved` isn't solved,
  and reports problems in that list that now are. Its command-line options
  select other configurations (these skip the regression check):

  ```sh
  pixi run fpm test test_hs_suite --profile release -- results.md --linesearch=funnel
  pixi run fpm test test_hs_suite -- /dev/null --problem=71 --print=2   # one problem, with its iterations
  ```

  The options are `--linesearch=`, `--merit=`, `--penalty=`, `--hessian=`,
  `--qp=`, `--restoration=`, `--trust-region`, `--no-interpolate`,
  `--nonmonotone=N`, `--lbfgs-memory=N`, `--problem=N`, `--print=L`, and
  `--web-data=FILE` (see the header of `test/test_hs_suite.f90`).
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

## Tools

| | |
|---|---|
| `coverage.sh` | runs the tests with `--coverage` and makes an lcov HTML report in `coverage/html` (dark-mode aware) |
| `tools/hs_performance_table.sh` | runs the HS suite in every configuration of the guide's Performance table and prints the table rows. It also regenerates the data of the results page (`web/js/hs_results_data.js`, `web/js/hs_slsqp_data.js`). Run it whenever a change affects the HS results. |
| `tools/hs_compare.sh N [options]` | runs HS problem `N` with SQPOPT and SLSQP, printing both solvers' iterations, for investigating a difference |
| `python/` | a Qt options dialog for SQPOPT, for use in other programs (see [python/README.md](python/README.md)) |

All are run from the repository root, e.g. `pixi run tools/hs_compare.sh 220`.

## Documentation and website

The website is `web/`: the user guide (`web/index.html`), the
Hock–Schittkowski results page (`web/hs_results.html`), and their CSS and
JavaScript. Update the guide along with any change to the API, options, or
behavior. When the HS results change, regenerate its Performance table with
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
| `sqpopt_qp_solver_module` | the QP subproblem front end, which chooses a solver |
| `sqpopt_qp_dense_module` | the dense active-set QP |
| `sqpopt_qp_reduced_hessian_module` | the sparse active-set QP (LUSOL basis, reduced-Hessian CG) |
| `sqpopt_linesearch_module` | the line searches (Armijo, exact, watchdog, filter, funnel) |
| `sqpopt_merit_module`, `sqpopt_filter_module`, `sqpopt_funnel_module` | the acceptance tests the line searches and the trust region use |
| `sqpopt_trust_region_module` | the trust-region globalization |
| `sqpopt_restoration_module` | feasibility restoration |
| `sqpopt_soc_module` | second-order corrections |
| `sqpopt_convergence_module` | the KKT convergence test |
| `sqpopt_linalg_module`, `sqpopt_dense_linalg_module` | sparse and dense linear algebra |
| `sqpopt_kinds` | the real kind (precision) |

## Continuous integration

On every push, [CI](.github/workflows/CI.yml) builds the pixi environment
from the locked `pixi.lock`, runs the tests with coverage, and builds the
FORD documentation. On `master`, it deploys `web/` (with the coverage
report and API docs) to GitHub Pages.

## Dependencies

Fetched and built by fpm:

- [LSQR](https://github.com/jacobwilliams/LSQR): iterative sparse least-squares solver
- [lusol](https://github.com/jacobwilliams/lusol): sparse LU factorization (the sparse QP's basis factors and updates)
- [fmin](https://github.com/jacobwilliams/fmin): derivative-free 1-D minimization (the exact line search)
- [slsqp](https://github.com/jacobwilliams/slsqp) (tests only): the SLSQP comparison

## License

[MIT](LICENSE)
