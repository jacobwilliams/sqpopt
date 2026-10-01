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
- ~~**`REAL32` builds fail inside the `LSMR` dependency**~~ Fixed: `LSMR`
  was dropped (2026-09-26, F14), and the `REAL32` build now succeeds.
  Still to add to the CI precision matrix.

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
fail.** *(Resolved 2026-09-27: `options%elastic_multiplier_limit`
re-solves the QP with a constraint elastic, at a bounded weight, when its
multiplier is large at two successive iterates and still growing, at most
3 times per solve. All 18 cases of `test/test_degenerate.f90` (three line
searches × two QP solvers × three corners) now reach the minimum. Before
this, the sparse QP reported a false success at the tangent points
with λ ≈ 1e6. HS default: 278/27/0 unchanged, `fc` 8,920 → 8,915
(TP221 42 → 8, TP13 68 → 42). The trust region doesn't use this yet.)* This carries over from Phase 0 and stays open. On
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

**HS suite bugs fixed (2026-09-26).** The Schittkowski/HS benchmark
(`test/test_hs_suite.f90`, 305 problems) went from 246 solved / 27 local /
32 failed to **269 / 32 / 4** with the default options, with fewer
evaluations on the problems solved both before and after (13,581 →
12,990 `f`):
- **ℓ1 merit slope overstated for non-QP steps** (the main one). The
  directional derivative was `gᵀp − μ‖viol‖₁`, exact only for a step that
  satisfies the linearized constraints. For a step shortened by the
  `max_step` cap it overstates the decrease, so no step length passes the
  Armijo test: every line search failed from the start point (TP59, 74,
  75, 83, 87, 109, 116, 236–239, 373, 392: all 176 `f` / 1 `g`). Now the
  exact one-sided derivative of `‖viol(c + αJp)‖₁` at `α = 0⁺`.
- **Absolute step cap.** `max_step = 2` made problems whose solution is
  far away (TP74/75 at `x ≈ 1000`) walk there in steps of 2 (700 `f`).
  The cap now adapts like a trust radius (`qp_solver%step_scale`:
  doubled after a capped step accepted in full, halved back toward 1
  after a shortened one): TP74 now takes 31 `f` (NLPQLP 10). A cap
  relative to `‖x‖` (SNOPT-style) was tried first, but the larger early
  steps raised the ℓ1 penalty and made TP26/27/375 crawl.
- **False infeasibility at a stationary point of the violation** (TP316–
  321, TP61). At a start with `J = 0` (a maximum of the violation), the
  linearization is inconsistent and the Gauss-Newton restoration step is
  zero, so the infeasibility test fired at iteration 2. Now (a) if the
  restoration step can't decrease the violation, the QP's elastic step
  is tried as the restoration direction (any decrease accepted), and (b)
  the infeasibility test also requires that the violation has stopped
  decreasing (by < 1% since the previous iterate, `viol_prev`), so a
  point that is stationary for the violation but not a minimum of it
  (TP61, just after a restoration step) isn't a stopping point.
- **Crawling line searches**: after 3 consecutive steps with `α < 0.01`
  the quasi-Newton Hessian is reset (as after a failed step). This fixed
  TP375 (which the adaptive cap had sent into a crawl) and TP61.
- **Harness**: a non-finite analytic derivative where the finite
  difference is finite now counts as a mismatch (TP25's gradient takes a
  fractional power of a negative number, a bug in the original code).
  TP25 then stops at its start, a flat plateau (every term underflows),
  as NLPQLP does.

**Remaining**: TP116, 332, 335, 355 hit `max_iter` crawling along the
constraints with a large ℓ1 penalty (8e7 on TP116), and 32 problems end
at other local solutions. **The filter line search solves all four**
(and TP61 in 12 `f` instead of ~7,400): with `linesearch_mode =
sqpopt_linesearch_filter` the suite solves **273 / 32 local / 0 failed**
with 10,593 `f` (vs 28,555 with ℓ1), and the scalable benchmark is as
fast or faster. The augmented-Lagrangian merit solves 270 (2 new
failures). So the ℓ1 penalty's monotone growth is now the main weakness:
either make the filter the default, or do F4 (a penalty that can
decrease).

**F4 done, and the filter is now the default (2026-09-26).** New option
`penalty_update`: `sqpopt_penalty_multipliers` (the old rule, still the
default for the merit-function line searches) or `sqpopt_penalty_model`,
each merit's principled rule (`update_penalty_parameter` in the line
search module):
- **ℓ1: Byrd–Nocedal model reduction** (N&W eq. 18.36): the penalty only
  increases to `(gᵀp + ½max(pᵀHp,0)) / ((1−ρ)Δv)`, with `Δv` the
  linearized violation reduction, instead of tracking `‖λ‖∞`.
- **Augmented Lagrangian: Gill–Murray–Saunders–Wright.** A joint step: the
  multipliers move from `λ₀` toward the QP's, and the slacks from the
  merit's minimizer `s₀` toward the QP's linearized constraint values
  `clip(c+Jp)`, so that `r = c−s` changes by `d = Jp − q` (`= −r₀` for a
  consistent QP). `ρ̂ = (A + ½pᵀHp)/(−B)` with `A = gᵀp − ξᵀr₀ − λ₀ᵀd`,
  `B = r₀ᵀd`; `ρ` increases to `max(ρ̂, 2ρ)` when below it, and decreases
  to `max(ρ̂, ρ_f, √(ρ·max(ρ̂, ρ_f)))` when above 4× that, with the floor
  `ρ_f` doubling after each decrease. The new multipliers are those
  along the step. (A first version kept the slacks' closed-form
  minimizer at each trial point: with a small `ρ` that `r` measures
  `λ/ρ`, not the violation, `rᵀJp` could be positive, and no penalty
  gave descent — 6 line-search failures on HS.)

HS results for every combination (the README has the table): the filter
solves 273 / 32 local / 0 failed with 10,593 `f`; the best merit-function
variants 270–271. GMSW solves 266 (265–267 with the decreases disabled
or with the QP's multipliers after the step), so the classical rule is
not worse in practice here. **The filter line search is now the default**
(`options%linesearch_mode`), which also avoids the merit functions'
penalty growth; the scalable benchmark is as fast or faster.

**NLPQLP line-search options (2026-09-26)**, from the NLPQLP 4.2 user's
guide (`references/NLPQLP.pdf`), on the line search type:
- `interpolate` (**default on**): each backtracking step length is the
  minimizer of the quadratic through `φ(0)`, `φ'(0)`, `φ(α)`, safeguarded
  to `[0.1α, 0.5α]` (NLPQLP's Algorithm 2.1), instead of `α/2`. In the
  filter search it interpolates the violation if the trial made it worse
  (slope `−θ₀`), else the objective. HS suite, filter: 273 → **274**
  solved (TP259), `f` 10,593 → 10,369, no regressions (TP332 250 → 123,
  TP355 501 → 391). With the ℓ1 merit: `f` 28,555 → 18,317 (TP61 7,368 →
  466), but TP268 and TP375 then fail. Benchmark neutral or better
  (rosenbrock N=2000 80 → 65 `f`).
- `nonmonotone_len` (default 0 = off): a failed search is retried against
  the worst merit value (filter: the worst `θ` and `f`, within `θ_max`) of
  the last `nonmonotone_len` iterates (NLPQLP's MAXNM, used only in the
  error situation). No gain for the filter search (10: same results, 10%
  more `f`; 40: TP214 fails); with the merit searches it helps: ℓ1 270
  solved (TP332), augmented Lagrangian with interpolation 270 → 273.
  NLPQLP's case for it is noisy functions (its Table 2), which the HS
  suite doesn't test yet.

**API change: combined user functions (2026-09-26).** The four user
callbacks `f`, `g`, `c`, `jac` are replaced by two:
`fc(x, f, c, status, data)` (objective and constraints) and
`gjac(x, g, jac_val, status, data)` (gradient and Jacobian values), since
the pairs usually share intermediate results; `set_functions(fc, gjac,
hess, data)`. The evaluation layer keeps one cache per routine (so asking
for `c` where `f` was just evaluated is free), and the results report
`n_eval_fc`/`n_eval_gjac` instead of `n_eval_f`/`_g`/`_c`/`_jac`;
`max_evals` counts `fc` calls. HS suite unchanged (276 solved); `fc`
counts are ~0.1% above the old `f` counts (where the solver needed only
`c`, e.g. in restoration).

**Second-order escape before declaring infeasibility (2026-09-26).** TP88
failed with the sparse QP (`--qp=sparse`): its functions are even in
`x₂`, an exact QP step landed on `x₂ = 0`, and no first-order step leaves
that plane, on which the problem is infeasible. The dense QP escaped only
through roundoff. New `escape_step` (restoration module): before stopping
with `sqpopt_infeasible`, probe each variable whose Jacobian column is
negligible in the violated rows (up to 10) by ±0.1%, 1%, 10%, and continue
from the first point with a lower violation (at most 3 escapes per solve).
HS suite: `--qp=sparse` 275/29/1 → 276/29/0, `--qp=sparse-lsqr`
273/31/1 → 274/31/0; the default is unchanged. The harness has new
`--problem=N` and `--print=L` options for debugging a single problem.

**Funnel method (2026-09-26).** New line search mode
`sqpopt_linesearch_funnel` (Kiessling, Leyffer & Vanaret, as implemented
in Uno; see `plan/UNO_COMPARISON.md`): the filter is replaced by one
number, the funnel width `τ`. A trial point must satisfy `θ ≤ τ`; if the
switching condition `α(−gᵀp) > δθᵏ^s` holds it needs an Armijo decrease in
`f` (f-type), otherwise `θ ≤ βτ` (h-type, which shrinks `τ`). A failed
search shrinks `τ` toward `θₖ` and takes the usual restoration step. The
trust region supports it too (acceptance with the model decrease `q`).
Options `funnel_*` on the line search type, with Uno's defaults. HS suite:
line search 276/29/0 with 10,709 `fc` (filter: 276/29/0, 9,793); trust
region 253/41/11 (filter: 251/36/18). Uno's other width update rule
(`funnel_update=2`) fails 2 problems; `funnel_require_current`,
`s_θ=1.1`, and `κ=0.9` change only the `fc` count (±2%). The harness has
new `--linesearch=funnel` and `--trust-region` options.

**F2 feasibility restoration phase (2026-09-26).** New option
`restoration_mode`: `sqpopt_restoration_phase` (default) or
`sqpopt_restoration_gauss_newton` (the old single step). Modeled on Uno's
`FeasibilityRestoration` (see `plan/UNO_COMPARISON.md`):
- Entered when the filter or funnel line search, or the trust region
  (which had no restoration before), finds no acceptable step at an
  infeasible point. The point is added to the filter (or the funnel is
  tightened toward it).
- Each phase iteration solves a feasibility QP with the regular QP solvers
  (a copy, so the optimality QP's warm start is kept): a proximal objective
  `ζ(x−x_ref)ᵀp + ½ζ‖p‖²` toward the phase's starting point, with the
  linearized constraints enforced or, if inconsistent, their ℓ1 violation
  minimized by the elastic mode. Then an Armijo backtracking search on the
  ℓ1 violation against the linearized decrease. Falls back to the
  Gauss-Newton step if that fails.
- Exit: `θ ≤ restoration_exit_factor·θ_ref` (0.9) and acceptable to the
  filter or funnel, or feasible, or `restoration_max_iter` (50)
  iterations. Multipliers are kept.
- An inconsistent QP still takes the Gauss-Newton step (then the QP's
  elastic step). Using the phase there too loses TP61 in every
  line-search configuration: from `x₂=x₃=0` the pure feasibility phase
  stays in that plane (the constraints' gradients have no `x₂`, `x₃`
  components there) and ends at a worse local solution, while the elastic
  step, which includes the objective's gradient, leaves it.
- HS suite: every line-search configuration is unchanged (their searches
  practically never fail at an infeasible point first; filter without
  interpolation: 10,037 → 10,061 `fc`, same outcomes). Trust region with
  the filter: 251/36/18 → 252/40/13; with the funnel unchanged
  (253/41/11). `test_infeasible` now also covers the funnel search and the
  trust region (which, without the phase, stopped with
  `sqpopt_line_search_failed`). The Performance table has trust-region
  rows, and the harness a `--restoration=` option.
- Not done: a quasi-Newton model of the constraints' curvature in the
  feasibility QP (it uses `ζI`); the trust region's own radius inside the
  phase (the phase uses the line search).

**F7 exact Hessian (2026-09-26).** `options%hessian_mode =
sqpopt_hessian_exact` now uses the user's sparse Hessian of the
Lagrangian (`hess` in `set_functions`, pattern from
`set_hessian_sparsity`, each off-diagonal element given once), evaluated
once per major iteration at the current multipliers (scaled: the user sees
the original problem's multipliers). It is matrix-free in the solver:
`sqpopt_hessian_type` has an exact mode whose products and diagonal come
from the nonzeros, so both QPs and the trust region use it unchanged.
- **Indefinite Hessians** (a primal inertia correction without a
  factorization): the QPs now report negative curvature in the variables
  (`qp_solver%negative_curvature`; the dense QP's reduced-Hessian Cholesky,
  the sparse QP's CG). In exact mode, when the QP finds negative curvature,
  fails, or gives a non-descent step, the Hessian is shifted by `δI`
  (`δ` ×10, from `shift_min`·max|H|) and the QP re-solved (up to 15
  times); `δ` also grows where a quasi-Newton Hessian would be reset (a
  failed step, a run of short steps), and is divided by 3 after each good
  step. Without the negative-curvature test (shifting only on failures)
  the HS suite gave 267/29/9; without the retry loop, 244/28/33.
- **HS suite** (new harness option `--hessian=exact`, with Hessians by
  central differences of the analytic gradients; the 16 FD-gradient
  problems keep BFGS): 268/33/4 (BFGS: 276/29/0). On the 258 problems both
  solve: median iterations −29% (TP302 491 → 10, TP301 225 → 8, TP116
  166 → 25), but `fc` 14,315 vs 8,023, from a few problems where Newton
  steps crawl along a curved valley (TP210 20 → 420 iterations, TP281,
  TP380 to `max_iter`). Failures: TP99 (converges with violation 8e-5,
  likely the FD Hessian's accuracy at `f ≈ −8e8`), TP103, TP111 (line
  search), TP238 (`max_iter`). Trust region + exact: 254/37/14; funnel +
  exact: 267/34/4.
- **Large problems** (`test_large_sparse`, analytic Hessians): control
  N=500: 47 → 8 `fc`, and full convergence instead of the acceptable
  level; chained Rosenbrock n=2000: 65 → 14 `fc`, to a different local
  solution (1978.26 vs 1974.67; the problem is nonconvex).
- Also: `results%n_eval_hess`; input validation (exact mode without
  `hess`, Hessian indices out of range); `test_hs71` runs with the exact
  Hessian (6 `fc` vs 7–9) with both QPs.
- Side finding, not investigated: SR1 (`--hessian=sr1`) is much weaker than
  BFGS on the HS suite (239/35/31); the same negative-curvature shift could
  apply to it.
- Not done: a factorization-based path with true inertia control (needs an
  LDLᵀ solver, §8.3); a Hessian-vector-product callback (for problems
  whose Hessian is dense or expensive).

**Trust-region step cap bug (2026-09-27).** Found with the new
`tools/hs_compare.sh` (one problem, SQPOPT and SLSQP side by side) on
TP220 (starts 25000 from its solution): in trust-region mode every step was
cut to the QP solver's line-search step cap `max_step*step_scale` (2),
because `step_scale` is only adapted in the line-search branch; the steps
then never used the full radius, so the radius couldn't grow either, and
TP220 crawled 2 per iteration to `max_iter`. The trust region now raises
`step_scale` so the cap can't bind inside its box. Trust region with the
filter: 252/40/13 → 261/37/7 (`fc` 11,742 → 9,898); with the funnel:
253/41/11 → 264/39/2; with the Armijo ratio test: 225/39/41 → 231/39/35;
with the exact Hessian: 254/37/14 → 260/35/10. The line search is
unaffected. (TP220's line-search cost, 89 `fc` vs SLSQP's 19, is mostly
structural: at its solution the constraint's gradient is opposite to the
active bound's, and the linearization then only allows `x₁−1` to shrink by
⅔ per iteration; SLSQP escapes that by jumping to `x₁ = 1` early.)

**Automatic L-BFGS memory (2026-09-27).** Found with `tools/hs_compare.sh`
on TP300–302 (unconstrained quadratics, `n` = 20, 50, 100): SLSQP (dense
BFGS) takes `n+1` iterations, NLPQLP `2n+1` evaluations, and SQPOPT, with
10 `(s,y)` pairs, about `5n` (74, 225, 491): the limited memory loses
BFGS's finite termination on quadratics. `options%lbfgs_memory` now
defaults to `0`, automatic: `max(10, min(n, 100))`. A fixed larger memory
is not better: more pairs than variables keep stale curvature (TP355, `n`
= 4: 391 → 1,129 `fc` with 100 pairs; TP332, `n` = 2, fails with 20).
HS suite: 276/29/0 → 275/30/0, `fc` 9,793 → 9,099 (−7%), `gjac` 8,247 →
7,560 (−8%); TP302 501 → 150, TP301 233 → 83, TP300 79 → 48; worst TP116
+34, TP380 +31. TP391 (finite-difference derivatives, `f* = 0`) now stops
at `f = 2.8e-3`, just outside `rel_tol`, and is in `known_unsolved`. The
harness has a new `--lbfgs-memory=N` option. (The same investigation found
that TP294–299's NLPQLP counts are anomalous: 72 and 114 evaluations for
`n` = 6, 10, then 31–33 for every `n` from 16 to 100, while SQPOPT's and
SLSQP's grow linearly with `n`, as expected on the chained Rosenbrock
function.)

**Escape step in `real128` (2026-09-30).** With `-DREAL128`, TP88 stopped
as infeasible next to its symmetry plane `x₂ = 0`: without `real64`'s
round-off the iterates stay near the plane (`x₂` ≈ 1e-23, growing to 1e-8
by the time the violation is stationary), and `escape_step` didn't probe
`x₂`, because its Jacobian column (~1e-8) was far above the relative
threshold `sqrt(epsilon)` (1e-17 in `real128`). A column is now also
negligible below `ktol`, the tolerance of the stationarity test that
triggers the escape (new test `test_escape`). `real64` HS suite unchanged
(280/25/0, 8,963 `fc`; `--qp=sparse` and `--qp=sparse-lsqr` +6 `fc`, same
outcomes). `real128` (debug): 278/26/1 → 279/26/0. Also fixed
in the HS harness: the COMMON blocks shared with the `DOUBLE PRECISION`
problem code are declared with that kind in every build, and
finite-difference steps are based on the accuracy of the function values
(`hs_epsilon`). Still open in `real128`: TP61 ends at another local
solution, so `test_hs_suite`'s regression check (a `real64` baseline)
fails there. TP299 takes ~90 s (the dense QP in software quad
arithmetic), which looks like a hang.

**Dense Hessian for the dense QP (2026-09-30).** Found while looking at
why TP299 (`n` = 100, so 100 pairs) takes ~90 s in `real128`: it is the
slowest problem in `real64` too (2.1 of the suite's 3 s), with the time in
forming the dense Hessian from `n` Hessian-vector products, each with a
solve with the `2k x 2k` middle matrix. New `hessian%dense`: for BFGS it
applies the stored updates to `θI` directly (`2n²` per pair instead of
`4n² + 4nk`), for the exact Hessian it copies the nonzeros. TP299: 2.1 →
0.8 s (`real64`), 86 → 40 s (`real128`); whole suite 2.95 → 1.22 s
(`-O2`). The matrix is the same up to round-off, which is enough to move
the evaluation counts: default 280/25/0 unchanged, `fc` 8,963 → 9,173
(18 problems change, TP332 alone +163); the Performance table is
regenerated (the Armijo and trust-region rows move by a problem or two
either way). The automatic L-BFGS memory stays: with 10, 20, 50 pairs the
suite needs 9,638, 9,521, 9,846 `fc` (and 50 solves one fewer), against
9,024 with `min(n,100)` (same build), and only 0.5–0.9 s instead of 1.2 s.

Not done:

- What is left of TP299's time is mostly the LU of the `2k x 2k` middle
  matrix (once per iteration, for the products in the damping and the
  descent test). Eliminating its `-D` block leaves a `k x k` positive
  definite matrix (`θSᵀS + L D⁻¹ Lᵀ`, as in L-BFGS-B), about 8 times
  cheaper to factor.
- `lu_factor` takes the middle matrix to be singular when a pivot is below
  `1e-14` times its largest element, and the products then fall back to
  `θI`, dropping all the curvature. This happens with pairs stored on
  several HS problems (TP54, 87, 109, 220, 322, 333, 373, 376), and may be
  bad scaling between the blocks rather than singularity (not checked): without the fallback in `hessian%dense`, TP87 is
  solved (281/24/0). The block elimination above would avoid it.

**Out-of-memory status (2026-09-30).** New status code
`sqpopt_out_of_memory` (27): the dense QP solver checks the allocation of
its three large matrices (the dense Jacobian, the constraint rows, and the
Hessian), and the solve stops with that code, at the current point,
instead of aborting the program. The QP solver type keeps a flag so that a
failure in any QP solve of an iteration (also the trust region's and the
restoration phase's) ends the solve. New test `test_out_of_memory` (a
dense QP with `n` = 6,000,000; it uses about 1 GB itself). The solver's
other allocations are still unchecked (see "Allocation failures" in §4).

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
  constants are removed, and `lbfgsb` is dropped. `lusol` is now used
  (F13, F1), and `LSMR` was dropped (F14).
- **Allocation failures.** Only the dense QP solver's three large
  matrices are checked so far (`sqpopt_out_of_memory`, see "Out-of-memory
  status" above). Go through every allocation whose size grows with the
  problem (`n`, `m`, the Jacobian and Hessian nonzeros, the L-BFGS
  memory) and make it end the solve with that status instead of aborting
  the program. That includes explicit `allocate` statements, automatic
  (re)allocation on assignment, and automatic arrays (which can't be
  checked, so the large ones would have to become allocatable). The
  candidates are the rest of the dense QP (its working-set and
  null-space matrices), the sparse QP and its `LUSOL` factors, the
  L-BFGS storage, the work vectors of `solve` and `sqpopt_iterate`, the
  problem type's bound, scaling, and cache arrays, and the Python shim.
  Each component needs a way to report the failure to its caller, as
  the QP solver type's `out_of_memory` flag does.

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
- **F2: elastic mode and infeasibility detection.** *(Done 2026-09-26:
  elastic QPs in Phase 1, the restoration phase in Phase 4.)* Use SNOPT-style ℓ1
  elastic QPs when the linearization is inconsistent, plus a feasibility
  phase. This produces a real `sqpopt_infeasible` status (B6, B11).
- **F3: Powell-damped BFGS.** It is cheap once B1 is fixed (one
  `hv_product` gives `sᵀBs`). Keep the skip rule as a fallback. This is
  needed on non-convex problems, where skipped updates leave `H` stale.
- **F4: a principled merit and penalty.** *(Done 2026-09-26; see "Phase 4
  status".)* Implement the full
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
- **F7: exact Hessian mode.** *(Done 2026-09-26, the matrix-free part; see
  "Phase 4 status".)* A user sparse Hessian or a Hessian-vector
  callback. With the matrix-free option it works immediately in the
  reduced-Hessian QP's PCG (which already truncates on negative
  curvature). A sparse-factorization path needs inertia control, which
  means an LDLᵀ solver; that is an optional external dependency and a
  decision point.
- **F8: derivative checking and finite differences.** Add a derivative
  verifier (SNOPT "Verify level") and a finite-difference fallback with
  automatic sparsity detection. `NumDiff` could be added as a runtime
  dependency for this (it was a dev-dependency until 2026-09-27, when it
  was removed as unused).
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
  *(Python part done 2026-09-28, without a C API: `python/sqpopt`, a
  `scipy.optimize.minimize`-like `minimize` built with f2py; see
  `python/README.md`. There are no finite differences, so the gradient and
  the constraint Jacobians are required.)*

  **Follow-up: make the sparsity interface scipy-like** (not started).
  Today the Jacobian pattern comes from a non-scipy
  `NonlinearConstraint.jac_sparsity` field, or from the structure of
  `jac(x0)` if it is sparse, or it is dense. The Hessian pattern is always
  the dense lower triangle: sparse `hess` results are densified to `n × n`
  on every call. In scipy, sparsity is expressed by the *return types*: `jac` and
  `hess` may return dense arrays, sparse matrices, or (for `hess`) a
  `LinearOperator`, and `hess` may instead be a `HessianUpdateStrategy`
  (`BFGS()`, `SR1()`). The only explicit patterns scipy has are the
  finite-difference ones (`NonlinearConstraint.finite_diff_jac_sparsity`,
  and `jac_sparsity` in `least_squares`). The plan:
  - **Patterns from the return types, for both derivatives.** Take the
    Jacobian pattern of each constraint from its `jac(x0)`, and the Hessian
    pattern (lower triangle) from `hess(x0)` and each constraint's
    `hess(x0, v)`, if sparse; dense otherwise. Assemble the Lagrangian
    Hessian's values sparsely (no `n × n` array), and pass only the
    pattern's entries.
  - **Structural zeros.** A pattern read at `x0` misses entries that happen
    to be zero there (e.g. `2*x[j]` at `x[j] = 0`, a common starting point).
    Options: evaluate the structure at `x0` and at a perturbed point inside
    the bounds and take the union; and/or document that sparse results must
    keep a fixed structure (explicit zeros allowed), as `scipy.sparse`
    matrices built with a fixed `indices`/`indptr` do. A value outside the
    pattern should keep raising a clear `ValueError` naming the constraint
    and entry.
  - **Explicit patterns, scipy names.** For structures that can't be
    inferred, accept scipy's name `finite_diff_jac_sparsity` on
    `NonlinearConstraint` (it describes the same thing, even without finite
    differences), and a `hess_sparsity` for `minimize` (not in scipy). Drop
    or deprecate the current `jac_sparsity` field.
  - **`hess` as a strategy.** Map `hess=BFGS()` / `SR1()` (scipy's
    `HessianUpdateStrategy` instances, or our own equivalents) to
    `options%hessian_mode`, so scipy code ports unchanged. A
    `LinearOperator` result needs a Hessian-vector-product mode in the
    library first (the F7 "not done" item), so it would be rejected with a
    clear error until then.
  - **Tests with real scipy objects.** scipy isn't in the pixi environment,
    so the scipy classes and `scipy.sparse` results are only duck-typed
    today (the tests use a minimal `tocoo()` stand-in). Add scipy as a test
    dependency and run the tests with `scipy.optimize.Bounds`,
    `NonlinearConstraint`, `LinearConstraint`, and `scipy.sparse` Jacobians
    and Hessians (CSR, CSC, and COO, with duplicate entries).
  - **A sparse benchmark from Python.** Port `test_large_sparse`'s control
    problem to Python (sparse `jac` and `hess`), to check that the
    bindings stay O(nnz) per evaluation, and to compare against the Fortran
    run.

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
- **F14: try LSMR for the reduced-Hessian QP's projections.** *(Dropped
  2026-09-26: F1 replaced the LSQR projections with direct LU solves, so
  LSQR is only used by the fallback `sqpopt_null_space_lsqr` method, and
  `LSMR` was removed as a dependency, which also fixed the `REAL32`
  build.)* This was an
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

- **F15: a persistent SNOPT-style elastic phase.** *(Possible future
  update, not started.)* SQPOPT uses SNOPT's elastic idea in two places
  today, but only one QP solve at a time:
  - an **inconsistent QP**: the QP's elastic slacks start at weight
    `elastic_weight` × max(1, ‖g‖∞) (1e4; SNOPT's "Elastic weight" is
    1e6) and are raised ×100 (SNOPT: ×10) up to `elastic_weight_max`.
    If slacks remain, the QP reports `sqpopt_infeasible`, and the solver
    enters the **restoration phase** (Uno / filter-SQP style) instead of
    staying elastic;
  - **diverging multipliers**: `options%elastic_multiplier_limit`
    (2026-09-27) re-solves one QP with the offending constraints elastic
    at a fixed weight, at most 3 times per solve (see the resolved "open
    issue" under Phase 1).

  SNOPT instead enters an elastic *mode* that persists across major
  iterations: it solves the elastic problem
  \( \min f + \gamma \sum_i (v_i + w_i) \) s.t.
  \( c_l \le c(x) - v + w \le c_u \), \( v, w \ge 0 \), with
  \( \gamma \) = elastic weight × ‖g‖, leaving it once the constraints are
  satisfied, and raising \( \gamma \) (×10) while the elastic solution is
  infeasible for the original problem. A converged elastic problem with
  positive violation is SNOPT's infeasibility certificate.

  What it would take:
  - state in the iterate loop (elastic on/off, \( \gamma \)), with the
    elastic slacks carried through the QP (the forced elastic mode added
    for `elastic_multiplier_limit` already lets the QP make chosen rows
    elastic at a fixed weight), and the acceptance tests applied to the
    elastic objective (the merit function adds \( \gamma \) × violation;
    the filter and funnel would need the elastic objective as φ);
  - entry rules (an inconsistent QP; diverging multipliers, replacing
    the one-shot re-solve), exit rules (feasible again), and the
    weight-increase rule;
  - the infeasibility test and status (converged elastic problem with
    positive violation → `sqpopt_infeasible`).

  Open questions: whether it should replace the restoration phase for
  inconsistent QPs or complement it (e.g. elastic first, restoration if
  the elastic phase stalls); its interaction with the trust region; and
  the default weight. **Benchmark before adopting**: the restoration phase
  does well on the HS suite (0 failures with the default options), so
  compare solved/local/failed and `fc` on every configuration of the
  guide's Performance table, plus `test_infeasible` and `test_degenerate`.

## 6. Testing and infrastructure

- **CUTEst test problems, as a pure-Fortran harness** *(planned,
  deferred 2026-09-27)*. A test like `test_hs_suite` on problems from the
  CUTEst collection ([ralna/CUTEst](https://github.com/ralna/CUTEst): the
  library, LGPL v3, built with Meson; the problems are SIF files decoded by
  SIFDecode), without depending on CUTEst or SIFDecode at build time:
  translate the SIF files ourselves into self-contained Fortran.
  - **Why it's feasible.** SIF describes each problem in group partially
    separable form (objective and constraints as sums of groups, each a
    function of a linear term plus weighted element functions), and the
    element and group functions are written in a Fortran-like expression
    syntax *with their first and second derivatives supplied*, so the
    translated problems get exact gradients (and exact Hessians, for
    `sqpopt_hessian_exact`). Precedent: S2MPJ (Gratton & Toint, 2024)
    translates the CUTEst SIF files into self-contained Python, Julia, and
    MATLAB files plus a small runtime that evaluates the group structure;
    we would do the same with Fortran as the output.
  - **Pieces.** (1) A translator (a Python script in `tools/`, run once,
    not part of the build) that interprets SIF's parameter language (loops,
    indexed names, parameter arithmetic, size parameters) and writes one
    Fortran module per problem: the data tables (variables, bounds, start,
    groups and their linear terms, element uses, constraint types) and the
    element/group functions. (2) A small Fortran runtime assembling `f`,
    `c`, the gradient, and the sparse Jacobian (whose pattern comes from
    the structure, so the sparse QP is exercised naturally) from those
    tables. (3) A harness like `test_hs_suite` (sqpopt, optionally SLSQP;
    Markdown report and web data, e.g. a CUTEst tab on
    `web/hs_results.html`).
  - **Hard or uncertain parts.** SIF is a large language (S2MPJ needed a
    substantial translator to cover all of CUTEst), so start with a subset
    and grow it file by file. Most CUTEst problems have no validated
    optimum (some SIF files note one in a comment), so classify results by
    convergence status and by comparison with the other solvers rather than
    against known optima. Many problems have size parameters (up to 10⁵
    variables and more), so sizes must be chosen per problem. The license
    of the SIF problem files must be checked before committing generated
    code.
  - **Plan.** Start with 20–30 small, well-known CUTEst problems not
    already in HS (e.g. BT, ROSENBR, HIMMELxx, and SIF versions of HS
    problems, which can be cross-checked against `test_hs_suite`), and grow
    the translator and runtime until they all pass the derivative check
    (as in `hs_derivatives.f90`) and match S2MPJ's output or published
    values; then widen the set. Rough effort: a few sessions for a working
    translator on common problems, more to reach 100–300 validated
    problems.

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
  `slsqp` (a dev-dependency; `psqp` was removed on 2026-09-27), and add
  the Duran–Grossmann problem from the FilterSQP manual (the unfinished
  `test_dg.f90` was removed on 2026-09-27).
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
F11, and F12. *(F14, LSMR, dropped.)*

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
   with PCG only. See [INERTIA_CONTROL.md](INERTIA_CONTROL.md).
4. **Dependency trim.** Keep `lusol` for F1's sparse KKT solve, or drop
   `lusol`/`LSMR`/`lbfgsb`. *Recommendation (2026-09-26):*
   - keep `lusol`, for F13 *(done: now used)* and then F1;
   - keep `LSMR` only if F14 shows a real benchmark gain, otherwise drop
     it (which also fixes the `REAL32` build) *(done: dropped)*;
   - drop `lbfgsb`, which is unused and has no identified role *(done)*;
   - either way, remove the unused `solve_sparse_linear_system` wrapper
     *(done)*.

## 9. Cleanups

Code-health items from a review of the library (2026-09-28): duplicated
code, and places that make changes harder than they need to be. None
changes behavior. Items C1 and C5–C7 are low-risk and mechanical; C2 and
C3 are the biggest wins for future changes; C4 closes a gap in the
consistency checks. Each should leave the HS results unchanged (see
CLAUDE.md), which is the check that a cleanup is behavior-neutral.

- **C1: merge the filter and funnel line searches.** `filter_line_search`
  and `funnel_line_search` (`sqpopt_linesearch_module.f90`) are ~80-line
  near-copies. They differ only in the acceptance call, the minimum step
  (`alpha_lim` vs `alpha_min`), what is recorded after an accepted step, and
  the non-monotone retry, which only the filter has. Merge them into one
  routine with mode-specific hooks, dispatched like
  `globalization_acceptable`. (Decide whether the funnel should get the
  non-monotone retry too, and measure it.)
- **C2: a common parent type for the two QP solvers.**
  `sqpopt_dense_qp_type` and `sqpopt_reduced_hessian_qp_type` share 14
  fields:
  - the settings `max_iter`, `active_tol`, `opt_tol`, `feas_tol`,
    `elastic_weight`, `elastic_weight_max`, and `warm_start`;
  - the forced elastic mode's `force_sign` and `force_weight`;
  - the outputs `n_iter`, `n_working`, `n_slacks`, `negative_curvature`,
    and `warm_status`.

  The `forced()` function is identical in both, and the ×100
  elastic-weight escalation appears three times. The dispatcher
  (`sqpopt_qp_solver_module.f90`, `solve`) also copies the forcing inputs
  in, and the four outputs out, once per solver. Hold these once in a
  parent type, and read the outputs through it. The two solvers' defaults
  differ (`elastic_weight_max` is 1e10 in the dense QP and 1e8 in the
  sparse one), so the parent must allow per-solver defaults. The Python
  schema and `test_schema.py`'s `TYPES` table follow the new layout.
- **C3: a helper module for the tests' user functions.** 15 test programs
  each define the same `fc_obj_cons`/`gjac_grad_jacv` wrappers around their
  own `obj`/`grad`/`cons`/`jacv`. Adding `gjac`'s `accuracy` argument meant
  editing 19 files. Add a helper module whose `fc` and `gjac` call four
  simple procedures, passed through the library's `data` argument, so a
  test supplies only those. The next change to the user-function interface
  then touches one file.
- **C4: check the option limits against each other.** The limits are in
  three places:
  - Fortran's `validate_options` (`sqpopt_module.F90`, about 31 checks);
  - the Python schema's `minimum`/`maximum`;
  - the guide's option tables.

  `test_schema.py` checks the names, defaults, and types, but not the limits.
  Add a Python test that, for each option with a limit, solves a trivial
  problem through the bindings with a value just inside and just outside
  it, and expects success and `sqpopt_invalid_input`. It would also check
  the guide's table defaults against the schema, by parsing
  `web/index.html`.
- **C5: small helpers for repeated computations:**
  - `max_violation(c, c_lb, c_ub)` next to `l1_violation`, replacing the
    three inline max-norm violations (in `sqpopt_convergence_module.f90`,
    and twice in `sqpopt_iterate_module.f90`);
  - `sqpopt_merit_module.f90` re-implements `l1_violation` inline in the
    ℓ1 merit value: use `l1_violation`;
  - `lagrangian_gradient(jac, g, lambda)` replacing the six
    `sparse_matvec_transpose` + `g - jtlam` blocks (in
    `sqpopt_iterate_module.f90`, `sqpopt_convergence_module.f90`, and
    `sqpopt_merit_module.f90`).
- **C6: one evaluation check in the iteration.** The
  `if (problem%stop_requested) ... if (.not. sqpopt_all_finite(...)) ...`
  block appears four times in `sqpopt_iterate_module.f90`. Replace it with
  one internal `evaluation_ok()` that sets `done` and `istat`.
- **C7: name the sparse QP's hard-coded tolerances.** In
  `sqpopt_qp_reduced_hessian_module.f90`:
  - the LU pivot tolerance `1.0e-12` is repeated in two `factorize` calls;
  - the curvature test `1.0e-10` is repeated in both CG routines;
  - the working-set priority weights are named parameters in
    `initial_working_set` (`w_bound`, `w_ineq`), but literals (`1.0e-6`,
    `1.0e-4`, `1.0e-2`, in a different order) in `choose_basis`.

  Make them named module parameters, with comments on what they are for.
  Check whether the two priority orderings are meant to differ.
- **C8: move the printed output out of `sqpopt_solve`.** About 500 of
  `sqpopt_module.F90`'s ~1,100 lines are printing routines contained in
  `sqpopt_solve`. That includes `print_header`, `print_iteration`,
  `print_summary`, `print_solution`, the legend, and their formatting
  helpers. Move them to a `sqpopt_report_module`, which leaves `solve` short
  enough to read. This is medium effort: they use host association (the
  solver, the iteration info, the results), so they need an explicit context
  argument.
- **C9: one way for the step routines to evaluate the problem.** The line
  search takes `eval_f`/`eval_c` procedure arguments, with its own three
  abstract interfaces (`sqpopt_ls_objective_func`,
  `sqpopt_ls_constraint_func`, `sqpopt_soc_func`), and the iteration passes
  wrappers of the problem's cached evaluations. The trust region and the
  restoration take the problem object directly. Both styles make separate
  `f` and `c` calls, though the user API (and the cache) is one combined
  `fc`. Pick one style, and evaluate `f` and `c` together. The procedure
  arguments are only worth keeping if standalone testing of the line
  search (as in `test_acceptance`) needs them.
- **C10: decide whether to retire the LSQR null-space method.** The sparse
  QP carries two complete null-space methods. `sqpopt_null_space_lu` is the
  default; `sqpopt_null_space_lsqr` is kept "for comparison and as a
  fallback", and has its own CG, projection, and working-set code
  (`projected_cg`, `project_null`, `initial_working_set`,
  `add_independent_rows_lsqr`). Retiring it would remove a few hundred lines
  and the `LSQR` dependency. Before deciding, measure how often the
  fallback actually fires (on the HS suite with `--qp=sparse`, and on
  `test_large_sparse`), and whether `test_qp_fuzz` needs it as a reference.
- **C11: the Python bindings' result layout in one place.** The
  `iinfo`/`rinfo` layout is spread across three places:
  - the array constructors in `python/sqpopt/fortran/sqpopt_python.f90`;
  - the fixed sizes (13, 7) in `_sqpopt.pyf`;
  - the positional indexing in `python/sqpopt/_minimize.py`.

  Keep one list of field names in Python, checked against
  `sqpopt_py_info()`, so that adding a result is a one-line change on the
  Python side. Also, `_build.py` hard-codes `gfortran` and `ar`: take them
  from the environment (`FC`, `AR`), with those as the defaults.
- **C12: a version number in the code.** The version (0.1.0) is only in
  `fpm.toml` and `pixi.toml`. Add a `sqpopt_version` constant to the
  library, print it in the log header, and expose it as
  `sqpopt.__version__` in Python. A test would check that it matches
  `fpm.toml`.
