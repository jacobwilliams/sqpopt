# sqpopt Roadmap (v2): review findings and plan

Written 2026-09-25 after a full review of all 16 modules in `src/`
(~4,000 lines) and the test suite. Every item tagged **[verified]** was
confirmed with a probe program built against the library under
`-fcheck=all` (the probes are described in §6 and should become
regression tests). Untagged items come from reading the code.
[PLAN.md](PLAN.md) remains the historical design log. This document
supersedes its §5 backlog.

## 1. Where things stand

**Strengths.** The architecture is good. It has clean module boundaries
(problem, options, Hessian, QP, line search, trust region, convergence),
pluggable modes for each component, a sparse COO Jacobian, and a
matrix-free Hessian. It builds cleanly, and all current tests pass under
`-fcheck=all`. The dense and reduced-Hessian QP modes already solve HS71
well.

**Main weaknesses.**
1. There are several correctness bugs in the quasi-Newton and
   convergence code. Two of them silently degrade every run of the
   real-QP modes, and one makes the default mode report success at
   non-optimal points.
2. The default QP mode (`sqpopt_qp_composite`) is a heuristic that is not
   a real QP solve. It limit-cycles on HS71 and, as it turns out, only
   "works" there because of bug B2 below.
3. There are no failure paths. QP failures, line-search failures, and
   infeasibility never reach the caller. Every non-converged run ends as
   `sqpopt_max_iter_reached`.
4. Robustness plumbing is missing: infinite bounds, input validation,
   evaluation errors, scaling, and resetting state between `solve()` calls.
5. There is no CI and no benchmark suite. Some tests accept a 0.5
   solution error on runs that never converge.

**Measured impact of fixing just B1 and B2** (HS71, objective
evaluations `i_obj` and max x-error):

| mode | today | after B1 fix | after B1 + B2 fix |
|---|---|---|---|
| dense QP | 158 evals, 1.2e-6 | 38, 3e-7 | **22, 2e-7** |
| reduced-Hessian QP | 26, 3.4e-5 | 26, 3.4e-5 | **26, 2e-7** |
| filter + dense QP | 3203, 1.9e-4 | 166, 2e-7 | **22, 2e-7** |
| composite (default) | 21960, 0.35 (max_iter) | unchanged | **diverges, 4.6** |

## Phase 0 status: done (2026-09-25)

B1–B7, B9, B10, and B16 are fixed. The short-term part of F1 is also
done: the new `sqpopt_qp_auto` default picks the dense QP for `n ≤ 200`
and the reduced-Hessian QP otherwise, and composite is demoted to
legacy. All of the §6 regression tests exist except the random-QP fuzz
test, which moves to Phase 1 with B11. The whole suite passes under
`-fcheck=all -ffpe-trap=invalid,zero,overflow`.

What changed, and where:
- **B1, B4:** [sqpopt_hessian_module.f90](../src/sqpopt_hessian_module.f90).
  The BFGS `L`/`Lᵀ` blocks are fixed. SR1 now uses the compact L-SR1 form
  (the cached `w`/`denom` arrays are removed).
- **B2, B6, B7:** [sqpopt_iterate_module.f90](../src/sqpopt_iterate_module.f90)
  and [sqpopt_module.F90](../src/sqpopt_module.F90).
  - `y` now uses consistent multipliers.
  - `sqpopt_iterate` returns `done` + `istat`.
  - `solve` stops after `options%max_consecutive_failures` (default 5)
    failed iterations in a row.
  - The stalled-progress exit returns the new code `sqpopt_stalled` (7).
- **B3:** [sqpopt_convergence_module.f90](../src/sqpopt_convergence_module.f90).
  Adds the multiplier sign and complementarity test, with IPOPT-style
  multiplier scaling of the tolerances (`s_max=100`). An earlier
  `ktol·max(1,‖λ‖∞)` scaling let a diverging λ mask a non-KKT point, so
  it was not used. Also adds a local-infeasibility test: the projected
  `Jᵀr_c` must be ≤ `ktol·‖r_c‖`, which is relative to `‖r_c‖` so that
  nearly feasible points don't trigger it.
- **B5, B16:** `sqpopt_infinity` (1e20) and `sqpopt_problem_type%validate`.
  Bounds are clamped, and the problem and options are checked, returning
  `sqpopt_invalid_input` (6) with a message. `x0` is projected onto the
  bounds. Added `sqpopt_status_message` and `solver%status_message()`.
- **B9:** `initialize` keeps pristine copies of the line-search and
  trust-region objects, and `solve` restores them.
- **B10:** both active-set QPs keep the multipliers paired with the
  working set they were computed for. They also report
  `sqpopt_infeasible` when the final step violates a linearized row
  (`feas_tol`, default 1e-6).
- **Infeasible linearizations:** these now trigger a Gauss–Newton
  feasibility-restoration step
  ([sqpopt_restoration_module.f90](../src/sqpopt_restoration_module.f90)).
  This is a lightweight precursor to F2.
- **Composite QP:** now releases active rows whose multiplier has the
  wrong sign.

Results: every HS71 variant reaches `sqpopt_success` in 22–31 objective
evaluations with 2e-7 error. The exception is the tuned-LSQR variant,
which ends as `sqpopt_stalled` at 5e-6 error. `test_medium` went from 83
to 65 evaluations.

**New findings during Phase 0** (added to later phases):
- **Degenerate attraction (Phase 1, with F2/B11).** On
  `min x₁+x₂` s.t. `x₁²+x₂²=1` in the box `[-1,1]²`, starting from a
  box corner, the dense and reduced-Hessian modes converge to a point
  where the circle is tangent to a bound: (0,−1) or (−1,0). LICQ fails
  there and λ diverges, so they now honestly report `sqpopt_stalled`
  rather than a false success, but they still don't reach the optimum
  (−0.707, −0.707). An elastic-mode QP or multiplier safeguarding should
  fix this. The case is worth adding to the benchmark suite.
- ~~**Composite can't detect infeasibility.**~~ Moot: the composite mode
  was removed (2026-09-26).
- **`REAL32` builds fail inside the `LSMR` dependency** (`Real constant
  overflows its kind` in `lsmrModule.f90`), which nothing in sqpopt
  calls. Add this to the dependency-trim decision (§8.4) and to the CI
  precision matrix.

## Phase 1 status: done (2026-09-25)

All seven items are done. Items 5–7 are described after items 1–4
below.

### Items 1–4

Done: F3 (damping), B8 (no forced acceptance), F5/B13 (SOC done right),
and non-finite handling. The trust region's merit-model mismatch (part
of B15) is also fixed.

What changed:
- **Damped BFGS (F3).** `sqpopt_hessian_type%damping` (default on)
  applies Powell's damping, per Nocedal & Wright Proc. 18.2. A
  negative-curvature pair is now used, not skipped, and `B` stays
  positive definite. Tested in `test_hessian_consistency`.
- **No forced acceptance (B8).**
  - A failed line search or trust-region step takes **no step**
    (`x_new=x`). This applies to Armijo, exact, the watchdog's
    no-window case, the filter, and the trust region.
  - The iteration then resets the Hessian, so the next direction
    differs, and skips the stalled-progress test once, since "no move"
    is not convergence.
  - New defaults: `alpha_min` 0.1 → 1e-10 and `max_ls_iter` 20 → 40.
  - The Armijo test uses `min(dphi0,0)`, so a non-descent direction
    can't loosen it.
  - Merit comparisons allow a roundoff slack of `10·ε·max(1,|φ₀|)`
    (IPOPT's `Compare_le`). Without it, near a solution the line search
    backtracked about 29 times to α≈1e-9 on every iteration.
- **SOC (F5/B13).** [sqpopt_soc_module.f90](../src/sqpopt_soc_module.f90)
  is rewritten as `soc_step`. It corrects only rows that are active in
  the linearization or violated at the trial point, and projects the
  corrected step onto the variable bounds.
  - It is tried only after the first trial step is rejected **and** did
    not reduce the constraint violation.
  - The line searches receive it as an optional `soc` procedure (Armijo,
    watchdog, filter). The trust region calls it directly, within its
    box.
  - It no longer runs on every iteration, which saves 2 `f` and 2 `c`
    evaluations per iteration.
- **Steps stay within the bounds.** The QP dispatcher now clips every
  QP step onto the bounds it was given. The active-set QPs had been
  returning steps up to about 1e-11 outside them. `test_hs71` now
  asserts that no function is ever evaluated outside the bounds.
- **Non-finite values.** Trial points with a non-finite `f` or `c` are
  rejected everywhere: all line searches, the trust region, and the
  restoration step. Checks use `ieee_is_finite`, so there are no NaN
  comparisons, and the suite passes under `-ffpe-trap=invalid`. A
  non-finite `f`, `g`, `c`, or `J` at the current point stops the run
  with the new status **`sqpopt_function_error` (8)**. Tested in
  `test_nonfinite`.
- **Trust-region `pred` (part of B15).** It is now
  `φ(x) − φ(model)`, with the model's `f−q` and `c+Jp`, for either
  merit function. Previously it was ℓ1-based even with the AL merit,
  which made trust region + AL fail on the Maratos test.
- **`major_step_limit = huge` overflowed** in `initial_step_length`
  even though the docs recommend that value to disable the limit.
  Fixed.
- **New test `test_maratos`** (Nocedal & Wright Ex. 15.4). Every line
  search × merit combination, with and without the trust region,
  converges in 8–9 evaluations.

Results: HS71 needs 14–17 evaluations (22 after Phase 0),
`test_medium` 41 (65), and the constrained Rosenbrock problem in
`test_resolve` 70–87 (73–195).

**Finding carried into item 5.** The QP solvers' `opt_tol` was
absolute, so near a solution the reduced-Hessian QP's multipliers were
too imprecise for the KKT test. Item 5 made the tolerances relative, and
`test_hs71` again requires `sqpopt_success` from every variant.

### Items 5–7: QP robustness, elastic mode, Wächter–Biegler filter

**Random-QP fuzz test (item 5).** `test/test_qp_fuzz.f90` solves 400
seeded random QPs with each solver and checks the full KKT conditions:
primal feasibility, stationarity with correctly signed bound
multipliers, and multiplier sign and complementarity. The QPs cover
convex, nonconvex (indefinite SR1), degenerate (dependent, duplicated,
or bound-parallel rows; fixed variables; more equalities than
variables), and infeasible cases. The previous dense solver failed
**41/400**, with 16 wrong answers (not stationary, or wrong-sign
multipliers). The previous reduced-Hessian solver failed 44/400. Both
new solvers pass 400/400. They also passed 3,000/3,000 in a one-off
stress run with n ≤ 12 and m ≤ 15.

**Dense QP, rewritten** ([sqpopt_qp_dense_module.f90](../src/sqpopt_qp_dense_module.f90)):
- **Elastic mode (F2).** Only the rows violated at `p=0` get a slack
  with an ℓ1 penalty, so the problem stays about `n` wide. This gives a
  feasible start with no phase 1. The weight ρ starts at
  `1e4·max(1,‖g‖∞)` and is raised ×100 while slacks stay positive, up
  to `1e10` (same scaling), after which the QP reports
  `sqpopt_infeasible`.
- **Independent working set (B11).** The initial set is built with
  Gram–Schmidt independence checks. A row that blocks a step along the
  null space is independent by construction.
- **Nonpositive curvature (B12).** The new `dense_cholesky_curvature`
  either factors `ZᵀHZ` or returns a direction of nonpositive
  curvature, which is followed downhill to the nearest blocking row.
  The old `1e-10` pivot floor is gone.
- **Relative tolerances.** Ratio-test rates are measured relative to
  `‖a‖·‖d‖`, and optimality relative to `1+‖Hu+g‖∞`. Infinite bounds
  are skipped.
- **Multipliers.** The multiplier sign tolerance is scaled while
  excluding the elastic slacks' own bounds; including them had inflated
  the tolerance to about 1e-3. Multipliers get one step of iterative
  refinement, since the normal equations square the conditioning.

**Reduced-Hessian QP, rewritten** ([sqpopt_qp_reduced_hessian_module.f90](../src/sqpopt_qp_reduced_hessian_module.f90)):
- It uses the same elastic formulation, with `elastic_weight_max` 1e8.
- **Rows are stored in CSR form (E3).** Row products cost only that
  row's nonzeros.
- **Active bounds fix variables** instead of being LSQR rows, as SQOPT
  does. The fixed coordinates of every projection are then exactly
  zero, and bound multipliers are read from the residual. This was the
  key fix: projecting gradients dominated by ρ-sized slack components
  left roundoff in zero-curvature directions, which CG amplified into
  wild steps.
- **Proximal curvature on the slacks.** `δ = max(1,‖g‖∞)` keeps CG
  steps along slacks bounded. It doesn't change feasible solutions.
- **Projected CG hardening:**
  - it uses `gp·gp` instead of `r·gp`;
  - it is capped at the null-space dimension;
  - its step is re-projected at the end;
  - a curvature direction is re-projected before it is followed;
  - the step accumulated before truncation is now ratio-tested (it
    wasn't before).
- **Tolerances.** They are relative to the gradient on the free
  unknowns. `pcg_rtol` was added, and `lsqr_itnlim` now defaults to
  automatic.

**Exact line search.** It now also tries SOC when its minimizer falls
short of the full step. Without it, the exact search converged only
linearly (Maratos) once the always-on SOC was removed in item 3. The
constrained Rosenbrock problem went from `max_iter` at 2,011
evaluations to success in 207.

**Wächter–Biegler filter (F6, B14)**
([sqpopt_linesearch_module.f90](../src/sqpopt_linesearch_module.f90)):
- This replaces the Fletcher–Leyffer envelope with the switching
  condition, Armijo on f-type steps, filter margins `γθ`/`γφ`,
  augmentation only on non-f-type steps, `θ_max`/`θ_min`, and W&B's
  `α_min`. The switching powers are computed in log space, so they
  can't overflow.
- When the search fails, the iteration adds the current point to the
  filter and, if infeasible, takes a restoration step.
- The trust region's filter mode uses the same test, with `−q` in place
  of `gᵀp`.
- Old options were removed: `filter_beta`, `filter_alpha1/2`,
  `filter_ubd`, `filter_tt`, `filter_feas_tol`, and
  `filter_penalty_estimate`. This is an API change, but pre-1.0.
- Results: the constrained Rosenbrock problem went from 74 to 35
  evaluations. With the filter, even the composite mode now detects
  infeasibility, via restoration; `test_infeasible` runs both line
  searches.

**Open issue: attraction to points where constraint qualifications
fail.** This carries over from Phase 0 and stays open. On
`min x₁+x₂` s.t. `x₁²+x₂²=1` in the box `[-1,1]²`, from a box corner,
the active-set modes approach (0,−1), where the circle is tangent to a
bound. They now report `sqpopt_success` there, whereas after Phase 0
they reported `stalled`. At that point the constraint is violated by only
ε² < `ctol`, and λ ≈ 8e3 with λx₁ → ½, so it is an approximate KKT
(Fritz–John) point with a diverging multiplier, not a minimizer. The
composite mode reaches the true solution. The fix is an SQP-level
elastic strategy, like SNOPT's: keep the elastic weight bounded, and
enter elastic mode when the multipliers exceed it, instead of raising
ρ inside the QP. Also consider a separate status, or a warning, when the
multipliers are huge at "convergence". Add this problem to the
benchmark suite.

## Phase 3 status: done (2026-09-25)

**Decisions.**
- **Callbacks (§8.2):** keep procedure pointers and add `status` and
  `class(*)` data arguments. The API is intentionally broken for this.
- **F8 (the derivative checker and the finite-difference fallback):**
  both dropped at the user's request.

What changed:
- **Callback API.** Every user function (`f`, `g`, `c`, `jac`, `hess`)
  now takes `(..., status, data)`, and `report` takes `data`:
  - `status` is `0` on entry;
  - `> 0` means the function can't be evaluated at this `x`, so the
    point is treated like a NaN and rejected;
  - `< 0` requests a stop (`sqpopt_user_requested_stop`), after which no
    user function is called again (uncached evaluations return NaN);
  - `data` is `class(*)`, optional, and is the object given to
    `set_functions(..., data=)`. It is held by pointer, so the caller's
    object needs `target`, and updates to it are visible to the caller.
- **Evaluation layer.** `sqpopt_problem_type`'s `f`/`c`/`g`/`jac`
  methods now:
  - handle `status` and `data`;
  - count calls (`n_eval_f`/`g`/`c`/`jac`);
  - cache results (4 entries for `f` and `c`, 1 for `g` and `jac`);
  - apply scaling.

  The line search now uses internal `(x, f)`/`(x, c)` interfaces.
- **Results.** `sqpopt_results_type` (in `sqpopt_types_module`) and
  `solver%get_results` report the status and message, the iterations,
  the evaluation counts, the time, `x`, `f`, `c`, `lambda`, the bound
  multipliers `z` (new), and the KKT and feasibility errors. Everything
  is unscaled. `get_solution(x, lambda, z)` now takes an optional `z`.
- **Log.** `print_level` 1 prints an iteration table (objective,
  infeasibility, KKT error, α, and flags `R`/`Q`/`F`) plus a summary; 2
  adds the penalty, step norm, and QP iterations. Output goes to
  `output_unit`. The QP solvers now report `n_iter`.
- **Termination (F10).** New options `max_evals`, `max_time`,
  `obj_lower_limit`, and `acceptable_ktol`/`acceptable_ctol`/
  `acceptable_iter` (IPOPT-style). They come with new statuses 9–12:
  `max_evals_reached`, `time_limit_reached`, `unbounded`, `acceptable`.
- **Scaling (§4).** `options%scaling` is on by default, with
  `scaling_max_gradient` = 100. It uses IPOPT's gradient-based rule at
  `x0`: `s_f`, and each `s_c,i`, is `min(1, 100/‖∇‖∞)`, and the finite
  constraint bounds are scaled to match. The evaluations it needs are
  cached, so they aren't repeated. `lambda`, `z`, `f`, and `c` are
  unscaled on output and in the `report` callback. Tolerances apply to
  the scaled problem.
- **Warm start (F9).** `solve(x0, istat, lambda0=...)` and
  `options%hessian_scale0` (B₀ = scale·I; `hessian%initialize` has a new
  `scale0` argument). QP warm-start state still does not carry over
  between `solve` calls, which keeps B9's determinism.
- **Validation.** Every component's settings are checked at `solve`:
  line search, filter, trust region, the QP solver and both QP backends,
  and the new options (including that `output_unit` is open when
  printing). The problem is restored from a pristine copy (`problem0`)
  on each `solve`, since scaling changes the working copy's bounds.
- **Tests.**
  - `test_callbacks`: user data reaching every function and `report`,
    `status > 0` avoidance, and `status < 0` stopping with no further
    calls. This caught the line search calling `f` again after a stop
    request.
  - `test_results`: λ and `z` against hand-computed values, for problem
    scales 1 and 1e4 with scaling on and off; evaluation counts against
    the callbacks' own counters; a `lambda0` warm start converging in 1
    iteration.
  - `test_termination`: all four new criteria, plus logging to a scratch
    unit.
  - `test_input_validation` gained 6 cases for component settings and
    `lambda0` size.

  All 45 test programs pass under `-fcheck=all -ffpe-trap=...`.
- **Not done:** consolidating the knobs spread across the components
  into `sqpopt_options_type`. They are validated and documented, but
  still live on their components.

## Phase 2 status: done (2026-09-25)

A new scalable benchmark, `example/benchmark.f90`
(`fpm run --example benchmark --profile release`), measures evaluations
and time on two problems:
- a nonlinear optimal-control problem with `n = 2N+1` and `m = N+1`
  equalities, plus bounded controls;
- a chained Rosenbrock problem with `n/2` coupled circle constraints.

Each is run at a size where `sqpopt_qp_auto` picks the dense QP and at
one where it picks the reduced-Hessian QP.

| problem | before Phase 2 | after |
|---|---|---|
| control N=50 (n=101, dense QP) | 0.86 s, 50 f-evals | **0.038 s** (23×), 25 f-evals |
| control N=150 (n=301, RH QP) | 9.66 s, 63 f-evals, `stalled` | **0.51 s** (19×), 23 f-evals, `success` |
| rosenbrock n=40 (dense QP) | 0.018 s, 104 f-evals | 0.016 s, 54 f-evals |
| rosenbrock n=300 (RH QP) | 0.066 s, 111 f-evals | 0.041 s, 64 f-evals |

Test-suite evaluation counts also dropped. HS71 went from 14–18
evaluations to 8–10, `test_medium` from 41 to 29, and the Maratos problem
from 8 to 5.

What changed:
- **E1: evaluation cache.** `sqpopt_problem_type` keeps the last 4
  evaluations of `f` and of `c` (the `f`/`c` methods), matched on the
  exact bits of `x`. The major iteration, all line searches, the trust
  region, and the restoration step go through it. The accepted point is
  one of the line search's trial points, so this removes one `f` and one
  `c` evaluation per iteration. It also provides the evaluation counters
  `n_eval_f`/`n_eval_c` (for Phase 3's statistics). The cache is reset at
  every `solve`.
- **E2: cached middle matrix.** The L-BFGS/L-SR1 compact-representation
  middle matrix is LU-factored once per pair or `gamma` change, not per
  product, so `hv_product` costs O(nk). The singularity test is now
  relative. This turned out not to be the bottleneck (about 4% here), but
  it removes an O(k²n+k³) cost from every CG iteration.
- **E4: crash and warm starts.** This was the big win. Profiling showed
  that each reduced-Hessian QP took about 150 active-set iterations,
  because every dynamics equality started with its own elastic slack.
  Both QPs now start from the minimum-norm step satisfying an initial
  working-set guess:
  - the guess is the previous QP's final working set (a warm start, on
    by default, with `warm_start`), or else the equalities and fixed
    variables (a crash start);
  - violated bounds are added to the guess over up to 4 rounds, then the
    step is clipped;
  - elastic slacks are created only for rows still violated after that.

  In the reduced-Hessian QP, the initial working set's per-row LSQR
  independence check is skipped for rows of the previous working set,
  when that solve had no elastic slacks. Rows that had slacks were only
  independent in the extended space; trusting those broke the fuzz test,
  which is how this was found. The QP solver's state is reset on each
  `solve` (`qp_solver0`, like the line search's). `test_qp_fuzz` now
  solves every QP twice with the same object, so it also checks warm
  starts. It passes 400/400, and 3,000/3,000 in the n ≤ 12, m ≤ 15
  stress run.
- **E5: storage.** The L-BFGS pairs are in a circular buffer, so no
  columns are shifted. `test_hessian_consistency` checks a wrapped BFGS
  buffer. The Jacobian's sparsity structure is set once per `solve`, and
  only its values are refreshed each iteration.
- **E3** was already done in Phase 1 (the reduced-Hessian rows are in
  CSR form).
- **E6 (dense QR updates) is deferred.** At `n ≤ 200` the dense QP isn't
  a bottleneck; the dense control problem takes 38 ms in total.

**Where the time went after Phase 2** (control N=150; see F1 in "Phase 4
status" for since): about 90% was in LSQR,
through the reduced-Hessian QP's null-space projections, both in
projected CG and in the remaining independence checks. Options for
later:
- a looser default LSQR tolerance, tied to the QP's own relative
  tolerance;
- reusing LSQR work across CG iterations;
- the sparse direct KKT factorization of F1, which would replace most of
  these iterative projections.

## Phase 4 status: in progress

**F13 done (2026-09-26).** The reduced-Hessian QP's initial working set
is now picked by one rank-revealing `lu1fac` (threshold complete
pivoting, `keepLU=1` for the relative singularity test) on all the
candidate rows at once, via the new `independent_columns` in
`sqpopt_linalg_module`:
- Each candidate row is normalized and scaled by its priority (equality
  rows 1, variable bounds 1e-2, inequality rows 1e-4). Complete pivoting
  prefers large elements, so this steers which row of a dependent group
  is kept; LUSOL also keeps exact unit columns (equality bounds) first.
  A row is dependent if its diagonal of `U` is ≤ 1e-8 times the largest
  element of its column of `U` (as the old LSQR test's 1e-8).
- This replaces the per-row LSQR independence check and the Phase 2
  "trust the previous working set" shortcut (and `warm_independent`),
  so dependent and duplicated constraints are handled the same way on
  warm and cold starts. The old LSQR loop is kept only as a fallback if
  the factorization fails.
- Results: `test_qp_fuzz` passes 5,000/5,000 for both solvers (400 in the
  regular suite). The benchmark's iterations are unchanged and its time
  is within noise (the shortcut had already removed most of the checks'
  cost). The HS suite with the sparse QP forced has identical outcomes on
  all 305 problems (245 solved; evaluation counts changed on 4). New
  `test_independent_columns` unit test.
- Not done: flagging redundant equality constraints at the start of a
  solve with the same factorization (the QP already copes with them).
- Since F1, this path (`initial_working_set`) is only used by the LSQR
  null-space method; the basis method picks its initial basis and
  working set the same way (one rank-revealing LU, with priority
  weights), see below.

**F1 done (2026-09-26): SQOPT-style basis method in the sparse QP.** New
`null_space` option on `sqpopt_reduced_hessian_qp_type`:
`sqpopt_null_space_lu` (the default) or `sqpopt_null_space_lsqr` (the old
method, kept for comparison and as a fallback).
- **Slack formulation.** Every general row gets a slack (`J p + E e − s
  = 0`, with the elastic slacks `e`), so all constraints are bounds and the
  working set is just the fixed unknowns. The free unknowns are split into
  `m` basic ones (a nonsingular `B`, factorized by `lu1fac` with rook
  pivoting) and the superbasics `S`; `Z = [−B⁻¹S; I]`.
- **Direct solves.** Products with `Z`/`Zᵀ` are one `lu6sol` each; CG runs
  on `ZᵀHZ` in the superbasics (steps stay on the constraints by
  construction, so no reprojection); multipliers are `Bᵀy = g_B` exactly.
- **Updates, not refactorizations.** Fixing or freeing a superbasic, or
  freeing a fixed unknown, leaves `B` unchanged. Fixing a basic unknown
  swaps in the superbasic with the largest pivot in its row of `B⁻¹S`,
  via `lu8rpc` (Bartels–Golub); `B` is refactorized every 100 updates or
  if an update fails. The basics are recomputed from the constraints
  after each change, so there is no drift.
- **Initial basis and working set** from one rank-revealing `lu1fac`
  (TCP) on the **row-normalized** `[J E −I]`, with free row slacks as
  exact unit columns (LUSOL takes them first) and the candidates for the
  working set scaled down by priority (inequality slacks 1e-2, bounds
  1e-4, equalities 1e-6). Candidates not picked are fixed; those picked
  stay free (dropped from the working set only if dependent).
  Normalizing *columns* instead was tried and is bad: e.g. the control
  problem's control columns (one entry `−h`) became unit columns and
  were made basic, giving `B⁻¹S ~ 1/h` and ~20× more CG iterations.
- **Trusting CG on a face.** After an unblocked CG step the face is taken
  as solved (as the LSQR loop does); re-testing the recomputed reduced
  gradient can stall, because of cancellation with the large elastic
  multipliers (seen on rosenbrock N=6000: QPs hitting the iteration
  limit).
- New `sqpopt_lu_type` in `sqpopt_linalg_module` (`factorize`, `solve`,
  `replace_column`).
- **Results.** `test_qp_fuzz` now runs three solvers (dense, LU, LSQR):
  5,000/5,000 each, with the LU method held to the dense solver's 1e-6
  KKT tolerance (LSQR 1e-5). HS suite with the sparse QP forced: 243
  solved (245 with the LSQR method); TP13 moves from 7.9e-5 to 1.2e-4
  relative error (degenerate optimum, borderline), and TP372 now hits
  `max_iter` with every step capped at `max_step`, as the dense QP also
  does (it is in `known_unsolved`). Default results unchanged (246).
  Timings (release):

  | problem | LSQR | LU (row basis, interim) | LU (slack basis) |
  |---|--:|--:|--:|
  | control N=150 | 0.48 s | 0.02 s | 0.02 s |
  | control N=500 | 27.7 s | 0.38 s | 0.31 s |
  | control N=1500 | > 3 min (stopped) | 5.2 s | 4.1 s |
  | rosenbrock N=2000 | 5.7 s | 3.7 s | 2.6 s |
  | rosenbrock N=6000 | – | 78 s | 58 s |

- **Where the time goes now:** at large `m`, the one-change-per-iteration
  active-set method itself. E.g. rosenbrock N=6000 starts with 3,000
  violated constraints, and each QP takes ~9,000 active-set iterations
  (each row blocks and is released one at a time), each with a few CG
  iterations (Hessian products, about 70% of the time). SQOPT's CG
  option would take the same steps. Options: a better initial working
  set on infeasible starts, a preconditioner for `ZᵀHZ`, or the dense
  reduced-Hessian (SQOPT "Cholesky") option when `nS` is small.
- **Follow-ups (2026-09-26)**, each measured on the benchmark (release):
  1. *Basic starting step* (SQOPT-style: the guessed working set fixed at
     its bounds, superbasics zero, `B v_B = −N v_N`; one LU solve) instead
     of the minimum-norm LSQR one, which was 45% of control N=1500's time.
  2. *Free all wrongly-signed fixed unknowns at once* at a face optimum
     (the worst one only in the second half of the iteration limit, as a
     safeguard against cycling), instead of one per active-set iteration.
  3. *Diagonal (Jacobi) preconditioner* for CG on `ZᵀHZ`: the diagonal
     of the L-BFGS matrix (new `sqpopt_hessian_type%diagonal`, O(n k²))
     on the superbasics. Halves CG iterations per face on control.
  4. *Gradient-projection step* when CG is blocked by a superbasic: step
     as far as the basics allow, clip and fix every superbasic taken past
     a bound, recompute the basics; accepted only if the basics stay
     feasible and the QP objective decreases, with backtracking (halving,
     up to 8 times). Fixes e.g. the elastic slacks reaching zero together
     (rosenbrock: ~1,000 faces per QP → ~3).
  5. *Dense reduced Hessian* (SQOPT's "Cholesky" option, without updates)
     for faces with `nS ≤ dense_max_ns` (default 50), CG otherwise. No
     speed gain (forming `ZᵀHZ` costs about as much as CG), but with the
     sparse QP forced on the HS suite it solves 245 problems instead of
     243.

  Tried and rejected: making free *elastic* slacks unit columns in the
  basis choice (so the variables would be superbasic, and projectable):
  rosenbrock N=2000 went from 0.25 s back to 1.3 s.

  | problem | before (F1) | 1 | +2 | +3 | +4 (all) |
  |---|--:|--:|--:|--:|--:|
  | control N=500 | 0.31 s | 0.17 s | 0.12 s | 0.12 s | 0.06 s |
  | control N=1500 | 4.1 s | 2.0 s | 1.55 s | 0.96 s | 0.39 s |
  | rosenbrock N=2000 | 2.6 s | 2.6 s | 1.3 s | 1.3 s | 0.26 s |
  | rosenbrock N=6000 | 58 s | – | 27.5 s | – | 2.4 s |

  What remains on rosenbrock N=6000 is mostly its first QP (infeasible
  start): ~6,000 faces, each a *basic* variable reaching a bound (a basis
  swap each), at ~1.5 CG iterations per face. The dense option could be
  made much cheaper by updating its factors between faces, as SQOPT does.

## 2. Bugs: correctness (fix first)

| # | Issue | Where | Evidence |
|---|---|---|---|
| B1 | **Compact L-BFGS forward product has `L` and `Lᵀ` swapped** in the BNS middle matrix. `B·v` fails the secant condition for the newest pair (error 0.72) and disagrees with the two-loop inverse (`‖B·H⁻¹v − v‖ = 0.71`). Every consumer of `hv_product` sees the wrong Hessian: the dense QP (densifies `H`), the reduced-Hessian QP's PCG, the filter's `q`, and the trust region's `pred`. The fix (upper-right block `mid(p,k+q)=s_pᵀy_q` for `p>q`, lower-left its transpose) brings the secant error to 0 and consistency to 7e-16. | [sqpopt_hessian_module.f90:234-238](../src/sqpopt_hessian_module.f90#L234-L238) | **[verified]** |
| B2 | **Quasi-Newton `y` mixes two multiplier estimates.** `gl_prev = g_k − J_kᵀλ_k` but `gl = g_{k+1} − J_{k+1}ᵀλ_{k+1}`. The correct value is `y = ∇L(x_{k+1},λ_{k+1}) − ∇L(x_k,λ_{k+1})`. Fix: at the end of the iteration, store `gl_prev = g − J_kᵀ·new_lambda`. | [sqpopt_iterate_module.f90:118-124, 213-214](../src/sqpopt_iterate_module.f90#L118-L124) | **[verified]** (table above) |
| B3 | **The KKT test ignores the sign of the multipliers and complementarity for general constraints**, so it reports success at non-KKT points. `min −x₁−x₂` s.t. `0≤xᵢ≤10` (as constraints), started at `x=0`: the composite mode returns `istat=0` at `x=(0,0)`, `λ=(−1,−1)`. The true optimum is `(10,10)`. The test needs: `λᵢ ≥ −tol` at a lower bound, `λᵢ ≤ tol` at an upper bound, `|λᵢ| ≤ tol` when strictly inactive, and free for equalities. It should also use scaled tolerances (SNOPT-style, relative to `max(1,‖λ‖)`). | [sqpopt_convergence_module.f90:62-85](../src/sqpopt_convergence_module.f90#L62-L85) | **[verified]** |
| B4 | **The L-SR1 product is inconsistent.** The cached `wᵢ = yᵢ − B sᵢ` are computed against the *old* `γ`, and they go stale when `γ` changes or the oldest pair is dropped. Newest-pair secant error is 0.32. Fix: use the compact L-SR1 form (`B = θI + Ψ M⁻¹ Ψᵀ` with `Ψ = Y − θS`), rebuilt from `S`,`Y` on update, or freeze `θ` for SR1. | [sqpopt_hessian_module.f90:144-150, 166-172](../src/sqpopt_hessian_module.f90#L144-L150) | **[verified]** |
| B5 | **Infinite bounds overflow.** `x_lb=-huge` makes `x_ub-x_lb = +Inf`, so `print_iterations` crashes under `-ffpe-trap=overflow`. The same pattern appears in the convergence test, both QP solvers, and the composite active set. There needs to be a documented `sqpopt_infinity` (e.g. `1e20`), and `\|b\| ≥ infinity` should be treated as "no bound" everywhere. | [sqpopt_convergence_module.f90:70](../src/sqpopt_convergence_module.f90#L70), [sqpopt_qp_dense_module.f90:132-133](../src/sqpopt_qp_dense_module.f90#L132-L133), ... | **[verified]** |
| B6 | **Failure statuses never reach the caller.** `solve` only checks `converged` and user stop. The first `qp_istat` is never acted on, and line-search failures are ignored. `sqpopt_infeasible`, `sqpopt_qp_solve_failed`, and `sqpopt_line_search_failed` are never returned. An infeasible problem runs all 100 iterations and returns `istat=1` in all three QP modes. (This contradicts the "DONE" note in `notes.txt`.) | [sqpopt_module.F90:138-152](../src/sqpopt_module.F90#L138-L152), [sqpopt_iterate_module.f90:148](../src/sqpopt_iterate_module.f90#L148) | **[verified]** |
| B7 | **The stalled-progress test can report success after a failed step.** If the line search or trust region fails and returns `x_new≈x`, then `rel_f, rel_x ≤ tol` holds with feasibility and the run returns `sqpopt_success` while KKT is not satisfied. This needs a distinct status (e.g. `sqpopt_stalled`) that is not reported as success. | [sqpopt_convergence_module.f90:87-96](../src/sqpopt_convergence_module.f90#L87-L96) | |
| B8 | **Non-descent steps are forcibly accepted.** Armijo and filter accept `alpha_min=0.1` unconditionally, and the trust region accepts a rejected step. This breaks global convergence and can increase the merit function. Also, when `max_ls_iter` runs out first, the reported `alpha=alpha_min` does not match the returned `x_trial`. | [sqpopt_linesearch_module.f90:399-403, 650-656](../src/sqpopt_linesearch_module.f90#L399-L403), [sqpopt_trust_region_module.f90:223-229](../src/sqpopt_trust_region_module.f90#L223-L229) | |
| B9 | **Solver state leaks across `solve()` calls.** `linesearch%penalty`, the filter (`filter_ready`), the watchdog state, and `trust_region%ready/radius` are never reset. Only the Hessian and the mode fields are re-initialized. | [sqpopt_module.F90:132-136](../src/sqpopt_module.F90#L132-L136) | |
| B10 | **The QP failure path returns mismatched multipliers.** If the active-set loop exhausts its iterations, `coeff` comes from an older working set (or has size 0) while `orig_idx` is current. The result is an out-of-bounds read or wrong `λ`, and the caller uses `p` anyway. | [sqpopt_qp_dense_module.f90:290-295](../src/sqpopt_qp_dense_module.f90#L290-L295), [sqpopt_qp_reduced_hessian_module.f90:286-291](../src/sqpopt_qp_reduced_hessian_module.f90#L286-L291) | |
| B11 | **The active-set QPs are not robust to degenerate or infeasible subproblems.** Phase 1 is two projections, not a real feasibility phase. An inconsistent linearization is never detected. `dense_null_space` assumes the working set has full row rank (wrong `n_z` when rows are dependent or `n_active>n`). The Gram-matrix Cholesky is silently perturbed. There is no anti-cycling. | [sqpopt_qp_dense_module.f90:141-160](../src/sqpopt_qp_dense_module.f90#L141-L160), [sqpopt_dense_linalg_module.f90:36-96](../src/sqpopt_dense_linalg_module.f90#L36-L96) | |
| B12 | **An indefinite reduced Hessian gives huge steps.** Modified Cholesky floors non-positive pivots at `1e-10`, so a negative-curvature direction is scaled by ~1e10 and only `max_step` catches it. This matters for SR1 and a future exact Hessian. | [sqpopt_dense_linalg_module.f90:122](../src/sqpopt_dense_linalg_module.f90#L122) | |
| B13 | **SOC is always run and not bound-safe.** It runs every iteration (2 extra `f` and 2 extra `c` evaluations plus an LSQR solve) even when the full step would be accepted. `x+p_soc` is not projected onto the variable bounds, it corrects inactive inequality rows too, and the line search then backtracks along `p_soc` rather than the SOC arc. The standard approach: try SOC only after the full step is rejected, use active rows only, and project onto the bounds. | [sqpopt_soc_module.f90](../src/sqpopt_soc_module.f90), [sqpopt_iterate_module.f90:185](../src/sqpopt_iterate_module.f90#L185) | |
| B14 | **The filter deviates from the theory.** Every accepted *trial* point is added (instead of `x_k`, and only on h-type iterations). There is no switching condition or Armijo test on f-type iterations, and the stored `q` is from the previous step. This explains the filter's high evaluation counts. Replace it with the Wächter–Biegler rules (see F6). | [sqpopt_linesearch_module.f90:626-656](../src/sqpopt_linesearch_module.f90#L626-L656) | |
| B15 | Smaller issues. The descent safeguard re-solves the QP without re-updating the penalty and returns early without updating `x_prev`. After a watchdog backtrack to `x_opt`, the multipliers still come from the current point. The trust region's `pred` is ℓ1-based even with the AL merit (the ratio is inconsistent) and it runs SOC on every retry. The penalty only ever increases. | [sqpopt_iterate_module.f90:167-176](../src/sqpopt_iterate_module.f90#L167-L176), [sqpopt_trust_region_module.f90:146-171](../src/sqpopt_trust_region_module.f90#L146-L171) | |
| B16 | **No input validation.** Nothing checks `n≤0`, `lb>ub`, Jacobian indices out of range, `size(x0)/=n`, missing callbacks (a null procedure pointer segfaults), or `lbfgs_memory≤0` (`push_pair` writes column 0). `x0` is not projected onto its bounds, so the functions can be evaluated outside them. | [sqpopt_problem_module.f90](../src/sqpopt_problem_module.f90), [sqpopt_hessian_module.f90:166-177](../src/sqpopt_hessian_module.f90#L166-L177) | |
| B17 | Documentation and behavior disagree. `sqpopt_hessian_exact` silently falls back to BFGS (`eval_hess` is never called). The README says "damped BFGS", but only a skip rule exists. `m_eq`/`m_ineq` are unused (only `m` is), and the "equalities first" ordering is neither required nor used. | various | |

## 3. Efficiency

| # | Issue | Fix |
|---|---|---|
| E1 | `f` and `c` are re-evaluated at the start of each major iteration at the point the line search just evaluated. SOC adds 2+2 more evaluations per iteration. | Carry `f`/`c` for the accepted point out of the line search. Make SOC conditional (B13). |
| E2 | `hv_product` rebuilds and factors the `2k×2k` middle matrix on **every call**, which costs `O(k²n+k³)`. PCG and the dense QP's densify step (`n` calls) make it hot. | Cache `SᵀS`, `L`, `D` and the factorization, and update them only in `update_*`. |
| E3 | `sparse_dot_row` scans all `nnz` for each row, so the reduced-Hessian ratio test costs `O((m+n)·nnz)` per active-set step. This defeats the large-scale mode. | Build a CSR row-pointer copy of `J` once per major iteration. |
| E4 | Both active-set QPs cold-start from an empty working set every major iteration. | Warm-start from the previous major iteration's working set. This is the single biggest QP speedup near convergence (SNOPT/SLSQP do it). |
| E5 | The Jacobian `irow`/`icol` are copied and reallocated every iteration, and `push_pair` shifts every column (`O(nk)`). | Keep the structure in solver state. Use a circular buffer index. |
| E6 | The dense QP recomputes a full QR for each working-set change (`O(n³)`). | Acceptable for small `n`. Use QR updates later if dense-mode sizes grow. |

## 4. Architecture and API gaps (production-readiness)

- **Results and diagnostics.** Return iteration count, function and
  derivative evaluation counts, final objective, `c(x)`, KKT and
  feasibility residuals, and **multipliers for variable bounds** (not
  returned today). Add a `status_message(istat)` function.
- **Status codes.** Add `stalled`, `infeasible` (locally infeasible, with
  the point that minimizes infeasibility), `unbounded` (objective below
  a limit), `eval_error`, and `invalid_input`. Make `max_iter` the only
  "ran out of budget" code.
- **Callback design.** The user cannot signal an evaluation failure
  (NaN or a domain error), and there is no way to pass user data except
  module globals or internal procedures. Options:
  (a) add an optional `status` argument and pass `class(*)` user data;
  (b) switch to an **abstract `sqpopt_problem_class`** with deferred
  type-bound `eval_*` methods. Option (b) is the idiomatic modern-Fortran
  choice and also enables combined `f+c` / `g+J` evaluation
  (see `notes.txt`). *Decision needed: this breaks the API.*
- **Non-finite handling.** If a trial point yields NaN/Inf, backtrack
  (line search) or shrink (trust region) instead of propagating it.
- **Options.** Tuning knobs are spread across five types, and some are
  silently overwritten from `options` at `solve()`. Consolidate the
  commonly used ones, validate every option, and document precedence.
- **Output.** Add a configurable output unit and a proper iteration table
  with a header (iter, f, ‖c‖, KKT, α, penalty, #active, QP iterations,
  flags), plus a final summary. Today it prints two lines per iteration,
  both labelled `f`.
- **Scaling.** Nothing is scaled and all tolerances are absolute. Add
  automatic gradient-based scaling of the objective and constraints
  (IPOPT-style) and optional user variable scaling.
- **Dependencies.** *Partly done (2026-09-26):* the unused
  `solve_sparse_linear_system` wrapper and the `sqpopt_linsolve_*`
  constants are removed, and `lbfgsb` is dropped. `lusol` and `LSMR` are
  still dependencies but unused until F13/F14 put them to work (or they
  are dropped).

## 5. Features toward state of the art

- **F1: a real QP as the default, and a real large-scale QP.** *(Done
  2026-09-26; see "Phase 4 status".)* The
  composite step was removed (2026-09-26; it couldn't survive the B2
  fix, see §1). Short term (done): make `sqpopt_qp_dense` the default for
  small `n` and `sqpopt_qp_reduced_hessian` the default for large `n`
  (auto-select by `n`/`nnz`). Longer term: build a **sparse KKT
  active-set QP**. PLAN.md §3's "fundamental mismatch" is not
  fundamental: the compact L-BFGS form is `B = θI − W M Wᵀ` (sparse plus
  rank-`2k`), so the KKT system can be solved directly by a sparse
  factorization of the `θI` + `J_A` augmented system (`lusol`, which is
  already a dependency) plus a `2k×2k` Woodbury correction, or by
  bordering. This gives direct rather than LSQR-iterative accuracy at
  scale.

  **The recommended approach: SQOPT-style basis partitioning with
  `lusol`** (reviewed 2026-09-26).
  - Split the active general rows, restricted to the free variables,
    into a square nonsingular basis `B` plus the rest `S`. The null-space
    basis is then `Z = [−B⁻¹S; I]`, so every projection is one exact
    solve with `B`'s LU factors, instead of an LSQR solve to a
    tolerance.
  - When a variable enters or leaves the basis, update the factors with
    `lu8rpc` (column replacement, a Bartels–Golub-type update) instead
    of refactoring.
  - Handle the reduced Hessian `ZᵀHZ` densely when the null space is
    small (SQOPT's "Cholesky" option) and with CG otherwise (its "CG"
    option).

  This targets the measured hotspot: after Phase 2, about 90% of the
  large control benchmark's time is LSQR inside the reduced-Hessian QP.
  It also removes the iterative-accuracy problems patched in Phases 1–2
  (penalty-scaled roundoff, CG losing null-space membership, LSQR's
  early stop on dependent rows). It needs only `lu1fac`, `lu6sol`, and
  `lu8rpc`, which the dependency exports.

  LUSOL does not export row/column additions, so an augmented-system
  approach that adds constraint rows would have to refactor on every
  working-set change. That favors the basis partition over the
  augmented system above. LUSOL also gives no inertia, so it can't
  provide F7's inertia control (§8.3).
- **F2: elastic mode and infeasibility detection.** Use SNOPT-style ℓ1
  elastic QPs when the linearization is inconsistent, plus a feasibility
  phase. This produces a real `sqpopt_infeasible` status (B6, B11).
- **F3: Powell-damped BFGS.** It is cheap once B1 is fixed (one
  `hv_product` gives `sᵀBs`). Keep the skip rule as a fallback. This is
  needed on non-convex problems, where skipped updates leave `H` stale.
- **F4: a principled merit and penalty.** Implement the full
  Gill–Murray–Saunders–Wright augmented Lagrangian: joint `(x, λ, s)`
  step, the `ρ̂` threshold from Lemma 4.3, and allowing ρ to decrease.
  Real-QP multipliers now make this possible (PLAN.md §6.1's caveat is
  resolved). Alternatively, use the Byrd–Nocedal model-reduction penalty
  update for ℓ1.
- **F5: SOC done right** (B13): only after a rejected full step, active
  rows only, bound-projected, and integrated into both the merit and the
  filter tests.
- **F6: the Wächter–Biegler filter line search**: switching condition,
  Armijo on f-type steps, filter augmentation only on h-type steps,
  `θ_min/θ_max` margins, and a feasibility restoration phase (shared with
  F2).
- **F7: exact Hessian mode.** A user sparse Hessian or a Hessian-vector
  callback. With the matrix-free option it works immediately in the
  reduced-Hessian QP's PCG (which already truncates on negative
  curvature). A sparse-factorization path needs inertia control, which
  means an LDLᵀ solver; that is an optional external dependency and a
  decision point.
- **F8: derivative checking and finite differences.** Add a derivative
  verifier (SNOPT "Verify level") and a finite-difference fallback with
  automatic sparsity detection. `NumDiff` is already a dev-dependency;
  consider promoting it to a runtime dependency.
- **F9: warm start.** Accept user-supplied `λ₀` (and bound multipliers)
  and an initial Hessian scaling. Hot-start the QP working set across
  major iterations (E4) and across repeated `solve()` calls.
- **F10: termination options.** Scaled KKT tolerances, IPOPT-style
  "acceptable level" termination, a maximum number of function
  evaluations, maximum wall time, and an objective lower limit
  (unboundedness).
- **F11: linear constraints.** Flag rows as linear so the solver keeps
  them satisfied once feasible, skips SOC on them, and never
  re-linearizes them.
- **F12: interoperability.** A `bind(c)` C API, then a thin Python
  wrapper. This is how SLSQP-style solvers get adopted.

- **F13: LUSOL rank detection for the working set.** *(Done
  2026-09-26; see "Phase 4 status".)* This is a small,
  standalone first use of `lusol`, and a stepping stone to F1:
  - one `lu1fac` with threshold rook or complete pivoting (TRP/TCP,
    which reveal rank) on the candidate rows picks a linearly
    independent subset directly. That replaces the reduced-Hessian QP's
    one-LSQR-solve-per-candidate independence check in
    `initial_working_set`, which was about 50% of QP time before the
    Phase 2 "trust the previous working set" shortcut;
  - it handles duplicated and dependent constraints properly. Those
    broke both that shortcut and LSQR (its early stop on rank-deficient
    systems);
  - the same factorization can flag redundant equality constraints at
    the start of a solve.

  It needs the low-level `lu1fac`/`lu6sol` interface: `lusol_ez%solve`
  refactors on every call, so it only suits one-off solves.
- **F14: try LSMR for the reduced-Hessian QP's projections.** This is an
  optional, benchmarked experiment. LSMR solves the same least-squares
  problems as LSQR with the same COO interface (`lsmr_ez`), so it is
  close to a drop-in swap.
  - Its monotonically decreasing quantity is the normal-equation
    residual. For a projection that is `J·(projected vector)`, exactly
    the quantity that must be ≈ 0, so LSMR may stop earlier at the same
    projection accuracy.
  - It also offers local reorthogonalization (`localSize`) for
    ill-conditioned working sets.

  It is still iterative, with the same rank-deficiency guard, so it
  softens LSQR's problems rather than removing them (F1 does that). It
  is also what breaks the `REAL32` build (Phase 0 finding): adopting it
  means fixing that upstream or dropping the `REAL32` option. If the
  benchmark gain is small, drop LSMR, which also resolves §8.4 for it.

## 6. Testing and infrastructure

- **Promote the review probes to regression tests** (they were throwaway
  programs built against the library):
  - `test_hessian_consistency`: secant condition for the newest pair and
    `‖B·H⁻¹v−v‖` for BFGS and SR1 (catches B1, B4).
  - `test_kkt_sign`: the `min −x₁−x₂`, `0≤xᵢ≤10` problem from `x=0`,
    where every mode must reach `(10,10)` (catches B3).
  - `test_infeasible`: `x₁∈[0,1]` and `x₁∈[2,3]`, which must return
    `sqpopt_infeasible` (B6, F2).
  - `test_infinite_bounds`: `±huge` and `±1e20` bounds under
    `-ffpe-trap=invalid,zero,overflow` (B5).
  - `test_resolve`: two `solve()` calls on one object must give
    identical results and evaluation counts in every line-search and
    trust-region mode (B9).
  - A random-QP fuzz test that checks the QP solvers against a dense
    reference on primal feasibility, dual feasibility, complementarity,
    and multiplier signs, including degenerate and dependent-row cases
    (B10, B11).
- **Tighten existing tests.** `test_hs71` accepts `x_error<0.5` with
  `istat=max_iter` for the composite, AL, and watchdog variants. After
  F1 every mode should be required to reach `sqpopt_success`.
- **A benchmark suite.** Port roughly 40–60 Hock–Schittkowski problems
  plus a few scalable sparse problems (a discretized optimal-control
  problem, chained Rosenbrock with constraints). Record success rate and
  evaluation counts per mode as a checked-in baseline, compare against
  `slsqp`/`psqp` (already dev-dependencies), and finish `test_dg.f90`
  (Duran–Grossmann from the FilterSQP manual).
- **CI** (none today): GitHub Actions running gfortran (two versions)
  and ifx, with debug `-fcheck=all -ffpe-trap=invalid,zero,overflow`
  and release builds, plus the `REAL32` and `REAL128` precision
  variants. Run the benchmark on a schedule.
- **Documentation.** Once the fixes land, update the README (defaults,
  "damped BFGS" claim, status codes) and PLAN.md §3 (where the
  "mismatch" claim and the lessons learned under B2 are now outdated).

## 7. Proposed order of work

**Phase 0: correctness (small diffs, very large payoff).** B1, then B2
together with switching the default QP (F1, short-term part), then B3,
B4, B5, B6/B7 (status plumbing), B9, B10, B16. Land each one with its
§6 regression test. Expected result: every mode converges on HS71 in
about 20–30 evaluations, and failures are reported honestly.

**Phase 1: robustness.** F3 (damping), B8 plus F5 (no forced acceptance,
proper SOC), F2 (elastic mode, infeasibility), B11/B12 (degenerate and
indefinite QP handling), non-finite evaluation handling, and F6.

**Phase 2: efficiency.** E1–E5, then measure on the benchmark suite.

**Phase 3: API and usability.** Results and diagnostics object, bound
multipliers, iteration log and output unit, the callback redesign (§4
decision), option consolidation and validation, scaling, F8, F9, F10.
*(Done; F8 dropped. See "Phase 3 status".)*

**Phase 4: large scale and advanced.** F13 first (LUSOL rank detection:
small, and it exercises the low-level LUSOL interface), then the sparse
basis-factorization active-set QP (F1, long-term part), plus F4, F7,
F11, and F12. Try F14 (LSMR) as a benchmarked experiment along the way.

**Phase 5 (runs alongside every phase):** CI, the benchmark suite, and
documentation.

## 8. Open decisions

1. **Default QP solver.** *Decided:* auto-select dense for small `n` and
   reduced-Hessian for large `n` (`sqpopt_qp_auto`); the composite mode
   was removed (2026-09-26).
2. **Callback API.** An abstract problem class (breaks the API, cleaner)
   or procedure pointers plus a `class(*)` context argument (additive).
3. **External sparse LDLᵀ** (e.g. MUMPS, as an optional dependency)
   for exact-Hessian inertia control. Alternatively, stay matrix-free
   with PCG only.
4. **Dependency trim.** Keep `lusol` for F1's sparse KKT solve, or drop
   `lusol`/`LSMR`/`lbfgsb`. *Recommendation (2026-09-26):*
   - keep `lusol`, for F13 *(done: now used)* and then F1;
   - keep `LSMR` only if F14 shows a real benchmark gain, otherwise drop
     it (which also fixes the `REAL32` build);
   - drop `lbfgsb`, which is unused and has no identified role *(done)*;
   - either way, remove the unused `solve_sparse_linear_system` wrapper
     *(done)*.
