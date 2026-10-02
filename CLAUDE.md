# CLAUDE.md

How this library is developed, and what has to be kept in sync when anything changes.

## The project

`sqpopt` is a modern Fortran SQP (sequential quadratic programming) solver for sparse, nonlinearly constrained problems. It is built with [fpm](https://fpm.fortran-lang.org), inside a [pixi](https://pixi.sh) environment.

| Path | What it holds |
|---|---|
| `src/` | The library. `sqpopt_module.F90` has the solver (`initialize`/`solve`), input validation, and all the printed output. `sqpopt_iterate_module.f90` has one major iteration. There is one module per component: problem, options, Hessian, QP solvers, line search, merit, filter, funnel, trust region, restoration, SOC, convergence, and log. The optional sparse factorizations are layered: `sqpopt_symmetric_solver_module.F90` is the only code that uses MUMPS, and only inside `#ifdef HAS_MUMPS`; `sqpopt_kkt_module` is the KKT matrix of a QP working set on top of it; and `sqpopt_inertia_module` (the Hessian's shift), `sqpopt_qp_direct_module` (direct QP steps), and `sqpopt_least_squares_module` (direct least-squares solves) use that. `sqpopt_types_module.f90` has the status codes and the results type. |
| `test/` | Unit and regression tests (fpm auto-tests: every `test/*.f90` program is a test). It also has the Hock–Schittkowski (HS) harnesses: `test_hs_suite.f90` (305 problems, the main regression and benchmark test), `test_hs_slsqp.f90` (the SLSQP comparison), and `test_hs_solutions.f90`; and `test_scalable.f90`, which solves the bound-constrained functions of any size in `scalable_functions.f90` (its header has the command-line options for larger sizes). |
| `example/` | `hs71.f90` (mirrored in the guide's worked example), `benchmark.f90`, and `benchmark_large.f90` (large problems, for the timing tables of the guide's Performance section and of `plan/MUMPS_PLAN.md`). |
| `web/` | The user guide (`index.html`), the interactive HS results page (`hs_results.html`, with data in `web/js/*_data.js`), the settings study (`choosing_settings.html`), the results of `test_scalable` (`scalable_results.html`: its tables are written by hand from the test's output, so update them when those results change), and CSS/JS. CI deploys it to GitHub Pages with the FORD API docs (`web/api`) and coverage (`web/coverage`). |
| `python/` | The Python bindings (`sqpopt/`: a scipy-like `minimize`, whose extension f2py builds from the Fortran shim and signature file in `sqpopt/fortran/`), a Qt options dialog (`sqpopt_options/schema.py` describes every option), and their tests. |
| `tools/` | `hs_performance_table.sh` regenerates the guide's Performance table and the results-page data. `hs_compare.sh` runs one HS problem with both SQPOPT and SLSQP. |
| `plan/` | Design documents. `plan/ROADMAP.md` is the current roadmap and backlog. |

## Ground rules

- **Never `git commit`** (or push). The user commits changes themselves.
- **Build and test through pixi.** The system gfortran fails to link on the dev Mac.
  ```bash
  pixi run fpm build
  pixi run fpm test                                   # everything
  pixi run fpm test test_hs_suite --profile release   # the HS regression test
  pixi run fortitude check                            # lint (src/, example/, and the Python bindings' Fortran)
  pixi run test-mumps                                 # everything, in the build with MUMPS (HAS_MUMPS)
  pixi run test-mumps test_hs_suite --profile release -- /dev/null --hessian=exact --inertia --direct
  pixi run run-mumps --example benchmark_large --profile release -- --scale=10 --no-active-set
  pixi run build-python                               # rebuild the Python bindings' extension
  pixi run test-python                                # the Python tests (bindings and options dialog)
  ```
  fpm sometimes runs a stale build after edits (the old output appears, or a test crashes with a segmentation fault after a derived type changed). If results look unchanged when they shouldn't, delete `build/gfortran_*` and rebuild. Always do that before a measurement that will be written down.
- **Never initialize a local variable in its declaration** (`integer :: n = 0`). In Fortran that gives it the implicit `save` attribute, so it keeps its value between calls. Declare it, then assign it in the executable code. (Default values on derived-type *components* are fine.)
- Every program unit has `implicit none`. Reals use `wp` (from `sqpopt_kinds`), and literals are written `1.0_wp`.
- Style: 4-space indentation, lines up to 132 columns, single-quoted strings (see `fortitude.toml`). Match the density and tone of the surrounding comments.
- No dense `n x n` or `m x n` arrays in the library, except in the dense QP solver.
- Every `solve` must start from the configuration given to `initialize`. No state may carry over between solves (`test_resolve` checks this). New component state must be reset at the start of `solve`.
- Printed output must never stop the solver. Every `write` in the logging and printing code uses `iostat=`, and writes `****` if it fails.
- Status codes are always referred to by their named constants (`sqpopt_success`, …), never by their numeric values.
- **MUMPS is optional.** The default build must need nothing but fpm, so MUMPS is only referenced inside `#ifdef HAS_MUMPS` in `sqpopt_symmetric_solver_module.F90`. Other code tests `sqpopt_has_mumps` or an object's `enabled` flag (`kkt%enabled`, `inertia%enabled`, `least_squares%enabled`), and every test must pass in both builds. A change to the factorization-based options (`inertia_control`, `direct_qp`, `direct_least_squares`), or to code they share with the matrix-free paths (the QP front end `solve_qp_subproblem`, the QP re-solve loop in `sqpopt_iterate`, `hessian%shift`, the Hessian's compact form), is tested with `pixi run fpm test` and `pixi run test-mumps`. The MUMPS build is double precision only.
- The objects that hold a sparse solver (`sqpopt_kkt_type`, `sqpopt_least_squares_type`) have pointers inside: they are never copied, live for one `solve` (as locals of `sqpopt_solve`), and are freed in `finish`.

## Documentation conventions (FORD)

- Every module has a `!>` header describing what it does. Every procedure and type has a `!!` docstring.
- **Every dummy argument gets a trailing `!!` comment**, including arguments of callbacks and test/example routines. The `me` (passed-object) argument doesn't need one.
- Type components get `!!` comments that give their meaning and units or allowed values.
- `[[name]]` makes a FORD link, and it must name a real entity. For anything else (e.g. external modules), use plain `` `code` ``.
- Avoid comment lines that start with `word:` (e.g. `stopped:`). FORD parses them as metadata.
- Keep docstrings true when behavior changes. Stale docstrings are treated as bugs.

## What to update with every change

Before calling a change done, go through the items that apply.

### Always
1. **Docstrings and argument comments** in every touched routine, and in the module header if the module's role changed.
2. **Tests.** Add or extend a test in `test/` for new behavior (a new `test/test_<feature>.f90` is picked up automatically). Run the full `pixi run fpm test`, and report failures faithfully.
3. **Lint** with `pixi run fortitude check`.
4. **User guide** (`web/index.html`). Update any text, table, or code snippet that describes what changed. See the specific cases below.

### New or changed option (a field of `sqpopt_options_type` or of a component type)
- Add the field, with its default and a `!!` doc, to the type.
- Add a check to `validate_options` in `sqpopt_module.F90`, or to the component's own validation. Add a case to `test_input_validation` if the check is new.
- Add it to `print_header` if it affects the method or tolerances shown in the log.
- Add it to the Python schema (`python/sqpopt_options/schema.py`): an `_o(...)` entry in the right section, with the same default and limits. For an enum-like option, add a `Choice` tuple of the integer values, and add that tuple to `ALL_CHOICES`. `python/tests/test_schema.py` fails if a Fortran field and the schema disagree.
- Add a row to the guide's options table for that type, and add or update the prose section that explains the feature.
- If a component is added or moved, update the guide's "How configuration works" section and diagram, and the Python layout.

- Rebuild the Python bindings (`pixi run build-python`). The option setter is generated from the schema, and `minimize` refuses to run with a stale extension.

### New results field (`sqpopt_results_type`)
- Add it to the type in `sqpopt_types_module.f90`, with a doc. Results are reset with `sqpopt_results_type()` at the start of `solve`.
- Set it (in `finish`, or `count_events`, in `sqpopt_module.F90`).
- Add it to the printed summary if it is useful there.
- Add a row to the guide's Results table.
- If it's useful from Python, return it through the bindings: add it to `iinfo`/`rinfo` in `python/sqpopt/fortran/sqpopt_python.f90` (and their sizes `sqpopt_python_n_iinfo`/`sqpopt_python_n_rinfo` there, and in `python/sqpopt/fortran/_sqpopt.pyf`), and to the `OptimizeResult` in `python/sqpopt/_minimize.py`.

### New or renamed status code
- Add the constant to `sqpopt_types_module.f90`, keeping the numeric grouping: `0–2` success, `1x` limits and user stop, `2x` failures. Add a case to `sqpopt_status_message`.
- Make the code reach the caller. A code that ends the solve is returned by `sqpopt_iterate` with `done = .true.` (its docstring lists the codes it can return, and so does the comment at the `finish` call in `solve`'s loop). Any other non-success code from an iteration is only counted toward `max_consecutive_failures`.
- If a component returns it (a QP solver, the line search, …), follow it through every caller. A QP status goes through `solve_qp_subproblem` to `sqpopt_iterate`, and also to the trust region and the restoration phase, which solve QPs themselves. Add the code to the docstring of each routine that can return it, and to `qp_status_text` in `sqpopt_log_module.f90` if a QP solve can return it.
- Check every `select case` and comparison on status codes in `src/` and in the tests (`grep sqpopt_<an existing code>` finds them), including the flags in `print_iteration`.
- Update the guide's Status codes table (the row's class is `ok`, `warn`, or `bad`) and its grouping sentence, and mention the code where the guide describes the feature that returns it.
- Add a test that produces the code, and checks both `istat` and `results%istat`.
- Python: `minimize` passes the code and its message through (`status`, `message`), and `success` is `status <= 2`, so nothing changes unless the code is a new kind of success. Rebuild the bindings and run `pixi run test-python`. If a component type gained a state field to carry the code, list it as internal state in `python/tests/test_schema.py` (the test fails on any field that is neither an option nor listed there).
- The README doesn't list the codes. The HS results data (`web/js/hs_results_data.js`) stores status messages, so regenerate it if a message that appears there is reworded.

### New iteration-log flag or detail line
- In `sqpopt_module.F90`, set the flag in `print_iteration` and describe it in `print_legend`. Add it to the guide's `print_level` description.
- Detail lines (`print_level >= 3`) go through the log (`lg%put(sqpopt_log_detail, ...)`). `test_termination` checks that the level-3 log doesn't change results.

### Change to a user-facing interface (callbacks, `set_*` routines, public constants)
- Update **every** implementation in `test/`, `example/`, and the HS harnesses.
- Update the guide's code snippets: Quick start, User functions, and Scaling (`fc_y`/`gjac_y`). Keep the worked example identical to `example/hs71.f90`.
- Update the README if it shows the interface.
- Update the Python bindings' Fortran shim (`python/sqpopt/fortran/`) if the change affects it, then rebuild and run `pixi run test-python`.
- Note in the summary to the user that the change breaks existing user code.

### Change that can affect convergence (algorithm, defaults, tolerances)
Run the HS suite in release mode, and compare it with the baseline recorded in the `known_unsolved` comment of `test/test_hs_suite.f90`. That baseline is currently 280 solved, 25 local, 0 failed, and 9,173 `fc` calls.
- If the change affects the exact Hessian, also run `--hessian=exact` (270 solved, 32 local, 3 failed, 11,267 `fc`). If it affects the factorization-based options, run, in the build with MUMPS: `--hessian=exact --inertia` (274, 29, 2, and 9,472), `--hessian=exact --inertia --direct` (274, 29, 2, and 9,402), `--hessian=sr1 --inertia` (274, 27, 4, and 10,454), and `--direct` (279, 26, 0, and 9,544; its automatic L-BFGS memory is 10 pairs); and `benchmark_large` for the timings. None of these is regression-tested, so compare them by hand.
- If problems newly fail, it is a regression. Investigate it, don't just update the baseline.
- If results change, do all of the following:
  - Update the `known_unsolved` list and the counts and date in its comment.
  - Regenerate the Performance table with `pixi run tools/hs_performance_table.sh --mumps` and paste its rows into the guide's Performance section (`--mumps` adds the rows of the options that need MUMPS). The script also regenerates `web/js/hs_results_data.js` and `web/js/hs_slsqp_data.js`.
  - Update any numbers quoted elsewhere in the guide or README.
- If results *don't* change, revert the regenerated files whose diffs are only timestamps: `test/hs_suite_results.md` and `web/js/*_data.js`. A debug-profile run also rewrites `test/hs_suite_results.md` with slightly different counts, so the committed report must come from `--profile release`.
- When comparing alternatives (an option's value, a new rule), use the harness's command-line options (see `test_hs_suite.f90`'s header). Report the numbers, and give the reasons for the chosen default in its docstring or in the guide.

### Planning
- Ideas that are discussed but not implemented go in `plan/ROADMAP.md` as a roadmap entry. Mark implemented items there as done.
- Don't start implementing a roadmap item unless asked.
