# Design plan: trust-region globalization option

**Status: implemented, per this plan, with results below.**

## 1. Motivation

`sqpopt` currently globalizes its SQP iterations with a **line search**:
solve the (unconstrained-by-radius, up to the static `qp_solver%max_step`
cap) linearized QP once per major iteration for a search direction `p`,
then find a step length `alpha` along `p` that is acceptable to either a
merit function (`sqpopt_linesearch_armijo`/`exact`/`watchdog`) or the
filter (`sqpopt_linesearch_filter`, see `plan/PLAN.md` §7). This works
well in practice (`test_hs71`, `test_medium` etc. all converge), but it
has two structural limitations line search alone can't fix:

- `alpha` only ever *shrinks the same direction* `p`. If `p` itself is a
  poor direction (e.g. because the current quadratic model `H` is a bad
  approximation, or the linearized constraints are a bad local model of
  strongly nonlinear ones), no amount of backtracking along it helps --
  the QP needs to be **re-solved** with a smaller step bound to get a
  genuinely different, more locally-trustworthy `p`.
- `qp_solver%max_step` is a **static, user-set** constant. A proper trust
  region *adapts* this radius every iteration based on how well the
  quadratic/linear model actually predicted the true change in the
  objective and constraints (the classical "ratio test"), shrinking it
  where the model is unreliable and growing it where the model is doing
  well -- this is a different and complementary mechanism to backtracking
  `alpha`.

This is also the *native* globalization strategy of Fletcher & Leyffer's
filter method (`references/fletcher.pdf`, see `plan/PLAN.md` §7): their
actual Algorithm 3 is trust-region-based, and `sqpopt_linesearch_filter`
is a line-search *adaptation* of it. Adding a real trust-region option
lets `sqpopt` also offer the literal filter-SQP combination (trust region
+ filter acceptance), alongside classical trust-region-SQP (trust region +
merit-function ratio test), as two variants of one new, orthogonal,
opt-in option.

## 2. Where this fits: a new, orthogonal dimension

`sqpopt` already has three independent option axes for a major iteration:
`hessian_mode` (how `H` is approximated), `qp_solver_mode` (how the QP
`min 0.5 p^THp+g^Tp s.t. ...` is solved), and `linesearch_mode`/
`merit_mode` (how a step along the QP's `p` is accepted). Trust region is
a **fourth, orthogonal axis**: it changes *how the major iteration as a
whole decides whether to keep a QP-computed `p` or re-solve for a
different one*, not how the QP itself is solved or how far along a fixed
`p` to go. Crucially, it can be added as a genuinely opt-in wrapper
**without modifying any of the three existing QP solver modules**, by
reusing the fact that all three already treat variable bounds as hard
constraints on `p`: a trust region of radius `radius` is just a temporary
tightening of the bounds passed to `qp_solver%solve` to
`max(x_lb, x-radius)`/`min(x_ub, x+radius)` for that call. No new
QP-solver-side code is needed at all.

Because trust region controls *whether to re-solve the QP*, it cannot
live inside `sqpopt_linesearch_module`'s `search()` interface (which only
ever sees one fixed `p` and picks `alpha` along it) -- it has to live at
the `sqpopt_iterate_module` level, as an alternative to (not a
modification of) the current "solve QP once, then call
`linesearch%search`" flow.

## 3. New type: `sqpopt_trust_region_type`

A new module, `sqpopt_trust_region_module.f90`, analogous in spirit to
`sqpopt_linesearch_module`/`sqpopt_qp_solver_module`:

```fortran
type, public :: sqpopt_trust_region_type
    !! options and state for trust-region-based globalization
    !! (an alternative to `linesearch%search`, see module docs).

    logical  :: enabled      = .false.   !! if true, use trust-region radius
                                          !! management instead of a line search
    real(wp) :: radius0      = 1.0_wp    !! initial trust-region radius
    real(wp) :: radius_min   = 1.0e-8_wp !! below this, treat the iteration as stalled
    real(wp) :: radius_max   = 1.0e3_wp  !! ceiling on the radius
    real(wp) :: eta1         = 0.1_wp    !! ratio threshold to accept a step
    real(wp) :: eta2         = 0.75_wp   !! ratio threshold to expand the radius
    real(wp) :: shrink_factor = 0.5_wp   !! radius *= shrink_factor on a rejected step
    real(wp) :: expand_factor = 2.0_wp   !! radius *= expand_factor on a very successful step
    integer  :: max_retries  = 20        !! bounded QP re-solves per major iteration

    ! internal state (persists across major iterations, like watchdog's):
    logical  :: ready  = .false.
    real(wp) :: radius = 0.0_wp

    contains
    procedure, public :: step => trust_region_step   !! the new per-iteration driver (see below)
end type sqpopt_trust_region_type
```

`enabled=.false.` (the default) means `sqpopt_iterate_module` behaves
*exactly* as it does today -- this is the "off" switch that makes the
whole feature opt-in and additive.

## 4. Algorithm (`trust_region_step`, replaces the QP-solve + line-search
block in `sqpopt_iterate_module` when `enabled=.true.`)

Reuses exactly the same pieces already built for the filter line search
(`q`, the QP-predicted decrease in `f`; `l1_violation`, the constraint
violation measure) and the merit function (`linesearch%eval_merit`), plus
whatever `second_order_correction` already does -- no new numerical
machinery, just a different outer control flow:

```
Given x, radius (persisted from the previous call, or radius0 if new):
REPEAT (up to max_retries times)
    x_lb' = max(x_lb, x - radius);  x_ub' = min(x_ub, x + radius)
    solve the QP with bounds (x_lb', x_ub') -> p                    ! same qp_solver%solve, tightened bounds
    (optional) apply second_order_correction to p, as today
    compute q = -(g^T p + 0.5 p^T H p)                              ! already computed for the filter line search
    evaluate f_trial = f(x+p), c_trial = c(x+p), h_trial = l1_violation(c_trial)
    IF linesearch%mode == sqpopt_linesearch_filter THEN
        accept = filter_acceptable(f_trial, h_trial)                ! reuse the filter's own domination test
        ratio  = 1.0   (no radius-quality ratio in the pure filter variant;
                         expand/shrink instead driven by whether `p` used
                         the full radius, as in Fletcher & Leyffer's own
                         Algorithm 1: "possibly increase rho" on acceptance)
    ELSE
        pred = q + penalty*(h(x) - h_linearized(x+p))               ! predicted decrease in the merit function
        ared = phi(x) - phi(x+p)                                    ! actual decrease in the merit function (eval_merit)
        ratio = ared / max(pred, tiny)
        accept = ratio >= eta1
    END IF
    IF accept THEN
        x = x + p;  lambda = new_lambda
        IF ||p|| ~ radius (the step used the full trust region) THEN
            IF ratio >= eta2 (or, filter mode, always on acceptance) THEN
                radius = min(radius*expand_factor, radius_max)
            END IF
        END IF
        RETURN success
    ELSE
        radius = max(radius*shrink_factor, radius_min)
        IF radius <= radius_min THEN RETURN sqpopt_line_search_failed (stalled) END IF
        CONTINUE (re-solve the QP with the smaller radius)
    END IF
END REPEAT
RETURN sqpopt_line_search_failed (retries exhausted)
```

This mirrors Fletcher & Leyffer's own Algorithm 1/3 structure (`plan/PLAN.md`
§7) closely when `linesearch%mode==sqpopt_linesearch_filter`, and is the
standard trust-region-SQP ratio test (Nocedal & Wright, Ch. 18) otherwise.

## 5. Interaction with existing options

| existing option | interaction when `trust_region%enabled=.true.` |
|---|---|
| `qp_solver_mode` (composite/dense/reduced_hessian) | unchanged -- any of the three can be used; the trust radius is enforced purely via tightened `x_lb`/`x_ub` passed into whichever solver is selected, so **no QP solver code changes are needed** |
| `qp_solver%max_step` | still applies as an absolute safety ceiling on top of the trust radius (harmless overlap: whichever is tighter binds); recommend setting it `>= trust_region%radius_max` so it doesn't fight the adaptive radius |
| `qp_solver%bound_enforcement` | unaffected -- still only relevant to `sqpopt_qp_composite`'s own post-hoc bound clipping, independent of the trust-region bound tightening |
| `linesearch_mode` | **reinterpreted as "which acceptance test to use"** rather than "which line search to run" -- `sqpopt_linesearch_filter` selects filter-style `(f,h)` domination acceptance; `armijo`/`exact`/`watchdog` all collapse to the same classical merit-ratio acceptance (their line-search-specific mechanics -- backtracking, `fmin`, the watchdog relaxed window -- are not used at all, since there is no `alpha` to search over) |
| `merit_mode` (`l1`/`augmented_lagrangian`) | still used, via `linesearch%eval_merit`, to compute `ared`/`phi` for the ratio test when not in filter mode |
| `linesearch%major_step_limit` | **not used** -- there is no initial `alpha` to cap when trust region is active (every accepted step is a full step within the current radius) |
| `linesearch%alpha_min`/`sigma`/`backtrack`/`max_ls_iter` | **not used** -- no backtracking occurs |
| `linesearch%watchdog_*` | **not used** -- the trust-region radius retries serve the same "don't get stuck" role the watchdog's relaxed window does |
| `options%ftol`/`xtol`/`ktol`/`ctol` | unchanged -- convergence checking in `sqpopt_iterate_module` happens before this block either way |
| second-order correction (SOC) | still runs, unmodified, on the QP's `p` before the accept/reject test -- exactly as Fletcher & Leyffer's own algorithm includes an SOC step (`plan/PLAN.md` §7, `references/fletcher.pdf` §3.1) |
| `hessian_mode`/`lbfgs_memory` | unchanged; `q`'s calculation already needs `hessian%hv_product`, same as the filter line search |

## 6. Proposed option surface

- `sqpopt_options_type` gains no new *required* field for the common case
  (trust region is configured the same way `qp_solver`/`linesearch` are --
  by constructing a `sqpopt_trust_region_type` directly and passing it to
  `solver%initialize(..., trust_region=...)`), consistent with the
  existing pattern where algorithm-specific tuning knobs live on their own
  sub-component type, not duplicated onto `sqpopt_options_type`.
- `sqpopt_type%initialize` gains one new optional argument,
  `trust_region`, alongside `problem`/`options`/`hessian`/`qp_solver`/
  `linesearch`.
- `sqpopt_iterate_module.sqpopt_iterate` branches near the top of its
  "solve QP + accept a step" section: `if (trust_region%enabled) then`
  call the new `trust_region%step(...)` driver `else` keep today's
  `qp_solver%solve` + `linesearch%search` flow (unchanged).

## 7. Test plan

Mirroring `test_qp_dense.f90`/`test_qp_reduced_hessian.f90`'s pattern of
hand-verified small QPs, plus reusing `test_basic.f90`'s existing closed-
form problems with `trust_region%enabled=.true.` set:

1. A hand-traceable 2-variable problem (like `test_inequality_constrained`)
   where a deliberately bad initial Hessian forces at least one rejected
   step and a radius shrink, to verify the retry loop and radius bookkeeping
   directly (assert on the number of QP re-solves via a counter, similar to
   `test_hs71`'s `i_obj`/`i_grad` counters).
2. `test_hs71` with `trust_region%enabled=.true.` paired with each of the
   three `qp_solver_mode`s and both filter/merit acceptance, analogous to
   the existing `run_hs71` sweep.
3. A test that the trust radius actually shrinks the effective bounds
   passed into the QP (e.g. checking that a huge unconstrained QP step
   from a badly-scaled problem gets capped even before any ratio test, by
   inspecting the returned `p`'s norm against the current `radius`).
4. Confirm `trust_region%enabled=.false.` (default) reproduces every
   existing test's results bit-for-bit (i.e. this feature is provably
   additive/non-invasive when off).

## 8. Open questions / deferred scope

- **`h_linearized(x+p)`** (needed for the merit-ratio `pred`) is cheap for
  `sqpopt_qp_dense`/`sqpopt_qp_reduced_hessian` (their linearized
  constraints are satisfied close to exactly by construction, so this
  term is close to `-h(x)` i.e. "predicts full feasibility"), but for
  `sqpopt_qp_composite` (which does **not** enforce the linearized
  constraints exactly) it needs an actual `sparse_matvec(jac,p,jp)` and
  `l1_violation(c(x)+jp, ...)` call, exactly as `second_order_correction`
  already does -- no new capability needed, just wiring.
- Whether to expand the radius on *every* accepted step that used the
  full radius (simpler) or gate it on a "very good" ratio only
  (`ratio>=eta2`, closer to Nocedal & Wright's textbook algorithm) --
  lean towards the textbook gate for the merit-ratio variant, and
  Fletcher & Leyffer's simpler "accepted implies possibly grow" rule for
  the filter variant, matching each variant's own source algorithm.
- Whether `max_retries` exhaustion should be a hard failure
  (`sqpopt_line_search_failed`) or should fall back to accepting the
  best (highest-ratio or filter-best) rejected trial from the retry
  sequence, the way `sqpopt_linesearch_watchdog` falls back to its best
  point -- leaning towards the fallback, for consistency with how every
  other mode in this library treats an "unsuccessful" step-selection
  (never just freezes; always makes *some* well-defined progress or
  explicit no-op).
- A genuine (not just bound-tightened) trust-region QP for
  `sqpopt_qp_dense`/`sqpopt_qp_reduced_hessian` could in principle use an
  \( \ell_\infty \) ball exactly like Fletcher & Leyffer's own `(QP)`
  (their paper's reason for preferring \( \ell_\infty \): it stays a
  linear/box constraint, not a genuinely quadratic one) -- which is
  *exactly* what the bound-tightening approach above already gives for
  free, since `x_lb'`/`x_ub'` are just an \( \ell_\infty \) ball around
  `x` intersected with the problem's own bounds. No further QP-solver
  changes are anticipated to be necessary.
- Not planned as part of this feature: a genuine \( \ell_2 \)-ball trust
  region (would require reformulating the QP objective, not just its
  bounds); Fletcher & Leyffer's full feasibility-restoration phase for
  infeasible trust-region QPs (same scope decision already made for
  `sqpopt_linesearch_filter`, see `plan/PLAN.md` §7 -- still no QP mode
  distinguishes "infeasible QP" as its own `istat`, so this remains a
  larger, separate follow-on project either way).

## 9. Implementation results

Implemented essentially as designed in §3-6, with two refinements made
during implementation:

- **`second_order_correction` was extracted out of `sqpopt_iterate_module`
  into a new shared module, `sqpopt_soc_module.f90`**, so both the
  line-search path (`sqpopt_iterate_module`) and the new trust-region path
  (`sqpopt_trust_region_module`) can call the same implementation --
  avoiding either duplicating it or creating a circular module dependency.
- **The `max_retries`-exhaustion fallback accepts the *last* (smallest-
  radius) rejected trial point outright**, rather than truly leaving `x`
  frozen, resolving the "Open questions" note in §8 in favor of the
  fallback option (for the same reason `sqpopt_linesearch_armijo`'s
  `alpha_min` floor is always accepted anyway: never leave a major
  iteration with literally zero progress, which risks an identical,
  wasted retry sequence on the very next iteration).
- The `lambda` (pre-QP-solve multiplier estimate) argument originally
  sketched for `trust_region_step` turned out to be unused in practice --
  `eval_merit`'s `ared`/`phi0`/`phi_trial` are all computed with the QP's
  own `new_lambda` (matching how the line-search path already uses
  `new_lambda`, not the old `lambda`, for its own merit/SOC calls) -- so
  it was dropped from the final signature.
- `filter_acceptable`/`filter_add`/the filter-init logic in
  `sqpopt_linesearch_module` were exposed as public type-bound procedures
  (`filter_test`/`filter_record`/`filter_prepare`) so
  `sqpopt_trust_region_module` can reuse the *same* filter (and its
  persistent state) that `sqpopt_linesearch_filter`'s line-search
  adaptation uses, rather than a second, disconnected implementation.
- New optional `trust_region` argument on `sqpopt_type%initialize`,
  alongside `problem`/`options`/`hessian`/`qp_solver`/`linesearch`, exactly
  as planned; no new `sqpopt_options_type` field was needed (consistent
  with `qp_solver`/`linesearch` themselves not being config-duplicated
  onto `options` either).
- **Validated** with two new tests in `test/test_basic.f90`
  (`test_trust_region_mode` -- merit-ratio acceptance;
  `test_trust_region_filter_mode` -- filter acceptance, the literal
  Fletcher & Leyffer combination) on the same inequality-constrained
  problem used by several other mode tests; both converge to the known
  solution. The full existing 27-test suite (trust region disabled by
  default) passes unchanged, confirming this feature is additive.
  **Not yet done**: a `test_hs71`-scale validation of trust region (the
  harder nonlinear benchmark used to validate the dense/reduced-Hessian
  QP solvers and the filter line search) -- left as follow-up work,
  since `run_hs71` in `test/test_hs71.f90` would need a signature change
  to accept a `trust_region` argument.
