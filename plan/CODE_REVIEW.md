# sqpopt code review (2026-09-27)

A review of the library as of commit `d050014`, covering correctness,
robustness, performance, maintainability, and ease of use.

**Method.** I read every library module (`src/`, 8,600 lines), with the
closest reading of the public API (`sqpopt_module`), the problem and
evaluation layer, the iteration, convergence, line-search, Hessian,
restoration, trust-region, and second-order-correction modules, and the
dense QP; the sparse QP (1,700 lines) was reviewed at the level of its
interfaces and key paths. Evidence from running the code: the full test
suite (55 tests, all passing), the Hock–Schittkowski and SLSQP harnesses,
the benchmark, the existing coverage data, and a set of edge-case probes
written for this review (a one-variable problem, a NaN starting point, fixed
variables, `lbfgs_memory = 1`, `max_iter = 0`, a constraint with no Jacobian
entries, a free constraint row, duplicate Jacobian entries, an exact-Hessian
linear problem, and a warm start at the solution).

## Summary

The library is in good shape. The code is careful, heavily documented,
free of leftover `TODO`s, `stop` statements, module-level `save` state, and
stray output (so separate solver objects are thread-safe, as long as the
user's functions are), and it validates its inputs thoroughly. Tests cover
93.6% of the library's lines. Every edge-case probe but one behaved
correctly.

The findings, by priority:

| # | Finding | Kind | Priority |
|---|---|---|---|
| B1 | A NaN in `x0` (or `lambda0`) is silently accepted | bug | high (easy fix) |
| U1 | No derivative checker for users | ease of use | high |
| U2 | The `report` callback must be passed as a procedure *pointer* | ease of use | high (easy fix) |
| P1 | The automatic L-BFGS memory makes some per-iteration costs grow with the square of the memory | performance | medium |
| R1 | Variables are not scaled; the absolute step cap slows problems with large variables | robustness | medium |
| R2 | A restoration phase whose steps keep failing never gives the optimality QP another chance | robustness | medium |
| U3 | Solver settings are spread over five objects | ease of use | medium |
| M1 | `sqpopt_iterate` takes 25 arguments of loose state | maintainability | medium |
| M2 | The line-search module is a 1,700-line grab-bag of strategies | maintainability | medium |
| T1 | Unused, unpinned dev-dependencies; a dead test file | tooling | low (easy fix) |

The details, and smaller items, follow.

## 1. Bugs and correctness

**B1. A NaN in `x0` is silently accepted** (high; easy fix). *(Done
2026-09-27: a non-finite `x0` or `lambda0` is now rejected with
`sqpopt_invalid_input`, "x0 must be finite" / "lambda0 must be finite";
tests in `test_input_validation`.)*
[sqpopt_module.F90:197](../src/sqpopt_module.F90#L197) projects the
starting point onto the bounds with `min(max(x0, x_lb), x_ub)`. With
gfortran, `max(NaN, x_lb)` returns `x_lb`, so a NaN component silently
becomes a bound, and the solver reports success from a point the user never
gave. (Probe: `x0 = [NaN]` with bounds `[-10, 10]` "converged" from
`x = -10`.) The same applies to a non-finite `lambda0`.
*Fix:* in the input checks (lines 177–188), reject a non-finite `x0` or
`lambda0` with `sqpopt_invalid_input` ("x0 must be finite").

**B2. `m_eq`/`m_ineq` are never checked against the bounds** (low).
*(Done 2026-09-27: `set_problem_size(n, m)` now takes the total number of
constraints; the `m_eq`/`m_ineq` components are gone, and every caller,
the README, and the guide are updated. This is an API change.)*
[`set_problem_size`](../src/sqpopt_problem_module.f90#L189) takes the
numbers of equality and inequality constraints, but only their sum is used:
whether a row is an equality is decided by `c_lb == c_ub`. A caller who
passes `m_eq = 2, m_ineq = 1` with bounds describing one equality and two
inequalities gets no warning, and the numbers in the problem object are then
wrong. *Fix:* either check them in `validate`, or (simpler for users) take
a single `m` and deprecate the split.

**B3. Smaller correctness points** (low).
*(Done 2026-09-27:
- LSQR's `istop` is now checked in the SOC and the Gauss-Newton
  restoration step. An iteration-limit exit counts as a failure, and the
  limit is `2*(rows+columns)+10`, as in the QP, instead of LSQR's default
  of 100.
- An SOC correction larger than the step is rejected.
- The non-monotone interpolation now fits through the reference value.
- `sqpopt_error` is removed.
- The SOC bound handling is unchanged: holding the variables at bounds
  fixed, in either of two variants, sent TP13 to 793 `fc` calls or TP116
  to a local solution. A note in the module explains this.
- `hessian_inverse_vector_product` is kept.
- HS default (release): 275/30/0, `fc` 9,099 → 9,059.)*
- The second-order correction
  ([sqpopt_soc_module.f90:114](../src/sqpopt_soc_module.f90#L114)) solves
  for the correction ignoring the variable bounds that are active, then
  clips the result, which can undo part of the correction; and it ignores
  LSQR's stopping reason (`istop`), as does the Gauss-Newton restoration
  step ([sqpopt_restoration_module.f90:142](../src/sqpopt_restoration_module.f90#L142)).
  Treat an LSQR iteration-limit exit as a failed correction.
- The non-monotone retry's step-length interpolation
  ([sqpopt_linesearch_module.f90:790](../src/sqpopt_linesearch_module.f90#L790))
  fits the quadratic through φ(0) rather than the relaxed reference value it
  tests against, so it can shorten steps the retry would accept.
- `sqpopt_error` (status −1) is defined and documented but never returned.
  Remove it, or use it.
- `hessian_inverse_vector_product` (and its CG solve) is never called by
  the solver; it is only unit-tested. Keep it only if it will be used (e.g.
  as a preconditioner).

## 2. Robustness

**R1. Variables are not scaled, and the QP step cap is absolute** (medium).
The objective and constraints are scaled (gradient-based, as in IPOPT), but
the variables are not, and `qp_solver%max_step` (2, doubling after capped
full steps) is an absolute length. Problems whose variables are large or
start far from the solution spend many iterations just growing the cap: on
TP220 (start 25,000 from the solution) the first 14 iterations are this
ramp. *Suggestions:* start the cap relative to `max(1, ‖x₀‖∞)`, or let
users supply variable scale factors (IPOPT's `x_scaling`), or both. (A
relative cap was tried once and hurt TP26/27/375, so this needs
benchmarking.)

**R2. A restoration phase whose steps keep failing stays in the phase**
(medium). In [`restoration_phase_iteration`](../src/sqpopt_iterate_module.f90#L570),
the phase's exit test (`restoration%done`, including its iteration limit)
runs only after a *successful* step. If the phase step and its Gauss-Newton
fallback both fail, the phase stays active, every later iteration retries
it from the same point, and the solve stops after
`max_consecutive_failures` with `sqpopt_line_search_failed`, without ever
trying the optimality QP again. *Fix:* on a failed phase step, count it
toward the phase's limit and end the phase after a failure (or after a
second one), so the next iteration re-solves the optimality QP.

**R3. The scale factors have no floor** (low).
[`compute_scaling`](../src/sqpopt_problem_module.f90#L469) divides by the
largest gradient at `x₀` without a lower bound, so a huge gradient at a poor
starting point (e.g. `1e12`) gives a scale of `1e-10`, and the tolerances,
which apply to the scaled problem, become meaninglessly loose in the
original units. IPOPT floors the scale at `1e-8` (`nlp_scaling_min_value`).

**R4. `max_evals` and `max_time` are only checked between iterations**
(low). A single iteration can make up to `max_ls_iter` (40) trial
evaluations, plus second-order corrections and restoration steps, so the
limits can be overshot noticeably. Checking inside the evaluation layer
(`raw_fc`) would make them hard limits; at least document the overshoot.

**R5. SR1 is much weaker than BFGS** (low; noted earlier). On the HS suite
`hessian_mode = sqpopt_hessian_sr1` solves 239 problems (BFGS: 275). The
negative-curvature shift added for the exact Hessian could be applied to
SR1 too, since SR1 matrices are also indefinite.

## 3. Performance

**P1. Some per-iteration costs grow with the square of the L-BFGS memory**
(medium). The new automatic memory (up to 100 pairs for large `n`) made
two costs 10–100× larger:
- [`factor_middle_matrix`](../src/sqpopt_hessian_module.f90#L405) recomputes
  all the `SᵀS` and `SᵀY` inner products (`O(nk²)`) whenever a pair is
  added, instead of updating them in `O(nk)`.
- [`hessian_diagonal`](../src/sqpopt_hessian_module.f90#L348) (the sparse
  QP's preconditioner, recomputed for each QP) costs `O(n(2k)²)`.

On the benchmark's chained Rosenbrock problem (`n = 2000`) the run time rose
from about 0.3 s to 0.5–0.7 s. *Fix:* update the inner products
incrementally (the standard compact L-BFGS bookkeeping), and form the
diagonal from `Ψ M⁻¹` once per Hessian change rather than per QP.

**P2. The dense QP is `O(n³)` per active-set iteration and dense in memory**
(low). [sqpopt_qp_dense_module.f90:183](../src/sqpopt_qp_dense_module.f90#L183)
builds an `(m+n+n_v) × (n+n_v)` row matrix and a null-space basis at every
iteration. The automatic mode only uses it for `n ≤ 200`, but a user who
forces `sqpopt_qp_dense` on a large problem gets no warning. *Fix:* warn,
or refuse, above a size limit.

**P3. Many per-call automatic arrays** (low). Evaluation and iteration
routines declare work arrays sized `n`, `m`, or `n + m` on every call
(e.g. [`raw_fc`](../src/sqpopt_problem_module.f90#L575)); for very large
problems these can exhaust the stack, depending on compiler flags. Moving
the larger ones into the solver object would make memory use predictable.

## 4. Maintainability

**M1. `sqpopt_iterate` passes loose state through 25 arguments**
([sqpopt_iterate_module.f90:82](../src/sqpopt_iterate_module.f90#L82)).
`x_prev`, `gl_prev`, `f_prev`, `viol_prev`, `jac`, `n_acceptable`,
`n_stalled`, `n_escape`, and `restoration` are all iteration state; each new
feature (the escape counter, the restoration phase) added another argument.
Bundle them in one `sqpopt_iteration_state` type owned by the solver.

**M2. The line-search type is a grab-bag**
([sqpopt_linesearch_module.f90:164](../src/sqpopt_linesearch_module.f90#L164),
1,700 lines). It holds the merit functions, both penalty rules, the
augmented Lagrangian's joint-step state, the filter, the funnel, the
watchdog, the non-monotone queue, and five search strategies, and the trust
region reaches into it for the filter and funnel. Splitting the acceptance
tests (merit, filter, funnel) from the step-length search, as Uno does
(see `plan/UNO_COMPARISON.md`), would make each piece testable on its own and
new strategies cheaper to add.

**M3. Two null-space methods in the sparse QP.** The LU basis method and the
LSQR method (kept as a fallback and for comparison) roughly double the
sparse QP's code (1,700 lines). If the LSQR path is only a fallback,
consider whether it earns its maintenance cost.

**M4. Implicit typing is enabled for the whole project**
([fpm.toml](../fpm.toml)) because about 40 routines in the legacy test file
`test/schittkowski_problems.f90` have no `IMPLICIT` statement. Every library
file declares `implicit none`, but a new file that forgets it would compile
silently. Adding `IMPLICIT DOUBLE PRECISION (A-H,O-Z)` to those routines would
let the project turn implicit typing off.

## 5. Ease of use

**U1. No derivative checker** (high). Wrong derivatives are the most common
cause of failure in practice (the HS collection itself shipped 18 problems
with wrong analytic derivatives), and sqpopt gives users no way to check
theirs. The test harness already has a robust checker
([test/hs_derivatives.f90](../test/hs_derivatives.f90): central differences
at three step sizes, one-sided at bounds). Promote it to a public routine,
e.g. `solver%check_derivatives(x, report)`, which compares `gjac` (and
`hess`) with finite differences and reports the worst entries. (F8 on the
roadmap was dropped; this is the smaller, most valuable part of it.)

**U2. The `report` callback must be a procedure pointer** (high; easy fix).
[`initialize`](../src/sqpopt_module.F90#L104) declares `report` as an
optional procedure *pointer* without `intent(in)`, so a user can't write
`report=my_report`; they must declare a pointer variable and point it at
their routine first (as `test/test_callbacks.f90` does). Declaring it
`procedure(sqpopt_report_func), optional :: report` (or adding
`intent(in)`) accepts a plain procedure name.

**U3. Settings are spread over five objects** (medium). Users configure
`options`, plus fields on the `linesearch`, `qp_solver`, `trust_region`, and
`hessian` objects, each passed separately to `initialize`; some options on
those objects are overwritten from `options` (e.g. `linesearch%mode`,
`qp_solver%mode`), which is surprising. The tested combinations in the
guide's Performance table could be offered as named presets (see
`plan/UNO_COMPARISON.md`), e.g. `options%preset = 'robust'`, and the most
common component settings could move into `options`.

**U4. Smaller ease-of-use points** (low).
- The Jacobian must be given in sparse (COO) form even for small dense
  problems. A helper such as `problem%set_dense_jacobian()` (pattern: every
  row and column, row-major values) would remove boilerplate for the common
  small case.
- There is no finite-difference fallback for problems without derivatives
  (the HS harness implements one; F8 dropped it from the library).
- There is no Hessian-vector-product callback for the exact-Hessian mode,
  which would suit problems whose Hessian is dense or expensive to form.
- `results%kkt_error` is for the scaled problem while every other result is
  for the original one; the guide says so, but a user comparing it with
  `ktol` in original units will be misled. Consider reporting both.
- The `hessian` argument of `initialize` is re-initialized from `options`
  on every solve, so only three of its fields (`damping`, `shift_min`,
  `shift_max`) matter; that is documented, but moving those three into
  `options` would remove the argument.
- There is no C API or Python wrapper yet (F12).

## 6. Testing and tooling

**T1. Dependencies** (low; easy fix). In [fpm.toml](../fpm.toml):
- `nlesolver-fortran` and `psqp` are dev-dependencies that nothing uses,
  and `NumDiff` is used only by `test/test_dg.f90`, which is almost entirely
  commented out. Remove them (and `test_dg.f90`, or revive it).
- None of the dev-dependencies is pinned to a tag, so a change upstream
  (e.g. in `slsqp`, which the SLSQP comparison now runs in every `fpm test`)
  can break the tests or change the results page without any change here.
  Pin them like the runtime dependencies.

**T2. Coverage gaps.** The coverage data (which predates the exact-Hessian
work) shows `sqpopt_iterate_module` at 83.8% of lines, the lowest; the
exact-Hessian paths, the restoration phase's failure branches, and the
trust region's restoration entry are the likely gaps. Re-run `coverage.sh`
and add targeted tests for what it shows.

**T3. The stale-build problem.** Twice this session, fpm did not recompile
modules that depended on a changed derived type (the solver then crashed
with a malloc error, or failed to build with "Mismatch in components"),
and only moving the build directory aside fixed it. It seems to be an fpm
dependency-tracking issue with `.mod` files. Worth a note in the guide's
"Developing" section ("if you see a malloc error or a component mismatch
after changing a type, rebuild from scratch"), and possibly an upstream
report.

**T4. Generated files in the repository.** `web/js/hs_results_data.js`,
`web/js/hs_slsqp_data.js`, and `test/hs_suite_results.md` are regenerated
by the harnesses and committed. Either generate them in CI (as the coverage
report is) or check in CI that they are current, so the results page can't
drift from the code.

## 7. Suggested order

1. The quick fixes: B1 (non-finite `x0`), U2 (`report` argument), T1
   (dependencies).
2. U1 (a public derivative checker): the highest-value addition for users.
3. P1 (incremental L-BFGS bookkeeping), since the new default memory made it
   matter for large problems.
4. R2 and R1, each benchmarked on the HS suite.
5. The refactorings M1 and M2, before the next large feature (e.g. F11 or
   the CUTEst harness), since both would otherwise add to them.
