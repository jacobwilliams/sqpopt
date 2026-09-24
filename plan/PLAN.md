# sqpopt Implementation Plan

This document lays out the plan for turning the current `sqpopt` skeleton into a
working, sparse-native Sequential Quadratic Programming (SQP) solver for the
general problem:

$$ \min_x f(x) \quad \text{s.t.} \quad c_l \le c(x) \le c_u, \quad x_l \le x \le x_u $$

where equality constraints are `c_l = c_u` and inequality constraints are
`c_l \ne c_u` (using `\pm\infty` for one-sided bounds).

## 1. Current state (v1 implemented)

**v1 is implemented and passing tests** (see [test/test_basic.f90](test/test_basic.f90),
run with `fpm test`): equality constraints, inequality constraints, and
variable bounds all solve correctly, with both the BFGS and SR1
limited-memory Hessian modes.

| Module | Responsibility | Status |
|---|---|---|
| `sqpopt_kinds` | selectable real working precision (32/64/128-bit) | done |
| `sqpopt_types_module` | status codes, `sqpopt_sparse_matrix` (1-based COO) | done |
| `sqpopt_problem_module` | NLP definition: sizes, bounds, Jacobian/Hessian sparsity patterns, user callbacks | done |
| `sqpopt_options_module` | tolerances, iteration limits, Hessian/linear-solver mode selectors | done |
| `sqpopt_hessian_module` | limited-memory (L-BFGS/L-SR1) Lagrangian Hessian, matrix-free | done |
| `sqpopt_linalg_module` | sparse matvec + dispatch to `lusol`/`LSQR`/`LSMR` | done |
| `sqpopt_qp_solver_module` | QP direction finder: v1 composite-step (default), dense active-set (§6.4), or sparse reduced-Hessian active-set (§6.2), all opt-in except the default | done |
| `sqpopt_qp_dense_module` | opt-in dense active-set QP solver (see §6.4/`DENSE_QP_PLAN.md`) | done |
| `sqpopt_qp_reduced_hessian_module` | opt-in sparse projected-CG active-set QP solver (see §6.2/`REDUCED_HESSIAN_QP_PLAN.md`) | done |
| `sqpopt_dense_linalg_module` | dense QR/modified-Cholesky helpers, used only by `sqpopt_qp_dense_module` | done |
| `sqpopt_linesearch_module` | \( \ell_1 \)/augmented-Lagrangian merit functions + Armijo (default)/exact/watchdog line search | done |
| `sqpopt_convergence_module` | projected-gradient KKT test + feasibility | done |
| `sqpopt_iterate_module` | orchestrates one major SQP iteration | done |
| `sqpopt_module` | public `sqpopt_type` facade (`initialize`/`set_problem`/`set_options`/`solve`/`get_solution`/`destroy`) | done |

Key design decision, upheld throughout **by default**: **no dense `n×n` or
`m×n` arrays, ever**. The Jacobian is sparse COO triplets with a fixed
sparsity pattern; the Hessian approximation is a matrix-free limited-
memory operator (only `O(n * lbfgs_memory)` storage, `lbfgs_memory` a
small constant). The one exception is the *explicitly opt-in*
`sqpopt_qp_dense` mode (§6.4), which forms dense arrays on purpose for
users who know their problem is small enough that the tight convergence
it buys is worth the `O(n^2)`/`O(mn)` memory -- never the default. For
larger sparse problems where that tradeoff isn't worth it,
`sqpopt_qp_reduced_hessian` (§6.2) gets the same exact-QP-solve benefit
while staying fully sparse/matrix-free (at the cost of needing more major
SQP iterations, since its `LSQR`-based null-space projections are
iterative rather than one-shot direct factorizations).

Deliberately **not yet implemented** (see §5 for the "optional/advanced"
backlog): `sqpopt_hessian_exact` mode (falls back to BFGS) and full
Powell damping for the BFGS update. The dense (§6.4) and sparse (§6.2)
rigorous active-set QP solvers, and the `ftol`/`xtol` no-progress
stopping tests and `options%print_level` diagnostic printing, are all
now implemented.

## 2. Dependency inventory & reuse strategy

We now have six reference dependencies fetched under `build/dependencies/`:
`slsqp`, `psqp`, `NumDiff`, `nlesolver-fortran`, `lbfgsb`, plus the sparse
linear algebra trio `LSQR`, `LSMR`, `lusol`, and `fmin`.

### Directly reusable (linked into `sqpopt` itself, `[dependencies]`) -- actually used in v1

- **`LSQR`** (`lsqr_module`, `lsqr_solver_ez`) — used directly inside
  [sqpopt_qp_solver_module](src/sqpopt_qp_solver_module.f90) for: (a) the
  Lagrange multiplier least-squares estimate, (b) the minimum-norm "normal
  step" that reduces linearized constraint violation, and (c) the
  null-space projection of the quasi-Newton "tangential step" (see §3).
- **`fmin`** (`fmin_module.fmin`) — a derivative-free 1-D minimizer (golden
  section / successive parabolic interpolation), used directly in
  [sqpopt_linesearch_module](src/sqpopt_linesearch_module.f90) to minimize
  the \\( \\ell_1 \\) merit function along the search direction, instead of a
  hand-written backtracking search.
- **`lusol`** (`lusol_ez_module.solve`) / **`LSMR`** (`lsmrModule.lsmr_ez`) —
  wired up in [sqpopt_linalg_module](src/sqpopt_linalg_module.f90) as
  `solve_sparse_linear_system`, a standalone general-purpose sparse
  square-system solver dispatching on the `sqpopt_linsolve_*` constants,
  for future use by a more rigorous QP solver (not currently called by
  any of the three QP solver modes, which each use `LSQR`/dense Cholesky
  directly; the `linear_solver_mode` option that used to sit on top of
  it was removed since nothing dispatched through it -- see §5).

`lusol`/`LSQR`/`LSMR` all take/accept 1-based COO `irow`/`icol`/`val`
triplets, which is exactly the convention used by `sqpopt_sparse_matrix` —
no conversion layer needed.

### `lbfgsb`: investigated, NOT reusable as originally hoped

**Correction to an earlier plan revision**: `lbfgsb_module` only exports
`setulb` (its monolithic bound-constrained reverse-communication driver);
`bmv` (the compact-BFGS matrix-vector product) and every other internal
routine are `private` module procedures with no `public ::` statement
making them accessible. There is therefore no way to reuse just its
two-loop-recursion/compact-BFGS machinery without copying private source
code. **v1 implements its own limited-memory BFGS (two-loop recursion for
the inverse, and the Byrd-Nocedal-Schnabel compact representation for the
forward product) and limited-memory SR1 (sequential rank-1 application)
directly in [sqpopt_hessian_module](src/sqpopt_hessian_module.f90)** — this
was unavoidable, not a case of reinventing something already reusable.
`lbfgsb` remains a dependency (added by the user) but is currently unused;
it could still be useful later for bound-only subproblems.

### Algorithmic inspiration only (NOT linked as dependencies — both are fully dense)

- **`slsqp`** (Kraft's SLSQP, BSD): dense damped-BFGS update, QP subproblem
  solved as bounded-variable least-squares (BVLS). Its \\( \\ell_1 \\)
  merit-function idea carried over into `sqpopt_linesearch_module` (adapted
  to the two-sided `c_l <= c(x) <= c_u` form); the damped-BFGS update was
  simplified in v1 to a standard curvature-condition "skip rule" instead
  (see §3), since implementing the compact-representation forward product
  for damping added complexity without a demonstrated need yet.
- **`psqp`** (Lukšan's PSQP, LGPL): dense active-set QP + Hoshino update.
  Not used in v1 (see §3 for why); remains the reference for a future
  proper active-set QP implementation.
- **`nlesolver-fortran`**: not an SQP algorithm; its sparse-mode plumbing
  (`irow`/`icol` sparsity pattern set once, dispatch to `lusol`/`LSQR`/
  `LSMR`) is the template already followed by `sqpopt_problem_module` and
  `sqpopt_linalg_module`. Kept as a dev-dependency for reference only.
- **`NumDiff`**: finite-difference derivatives + automatic sparsity-pattern
  detection (`dsm`). Kept as a dev-dependency; not yet used in the test
  suite (see §5) but a natural way to sanity-check hand-written analytic
  derivatives in future tests.

## 3. Algorithm actually implemented in v1

The originally-planned "primal-dual active-set QP with a sparse KKT solve
via `lusol`" turned out to have a fundamental mismatch with the matrix-free
Hessian design: an active-set method needs to *factorize* systems involving
`H`, but `H` is a limited-memory operator with no explicit sparse nonzero
structure to factorize (it's effectively dense, just never formed). This
was discovered empirically -- see below -- while building and testing v1,
and the design was corrected accordingly.

**v1 uses a simplified composite-step method** (implemented in
[sqpopt_qp_solver_module](src/sqpopt_qp_solver_module.f90)), entirely
matrix-free/sparse:

1. **Multiplier estimate**: least-squares solve of \\( J^T \\lambda \\approx g \\) via `LSQR`.
2. **Normal step** \\( p_n \\): minimum-norm `LSQR` solution of \\( J p_n = \\text{viol} \\),
   where `viol` is the linearized constraint violation (drives feasibility).
3. **Tangential step**: the unconstrained quasi-Newton direction
   \\( v = H^{-1} g \\) (matrix-free two-loop recursion), **projected onto
   the null space of \\( J \\)** via another `LSQR` solve (\\( p_t = v - J^T z \\),
   `z` minimizing \\( \\lVert J^T z - v \\rVert_2 \\)) so it does not reintroduce
   infeasibility. **This null-space projection turned out to be essential**:
   an early version without it (using the raw `v` as the tangential step)
   failed a basic equality-constrained test (converged to a nearby but
   wrong point) because the unprojected quasi-Newton direction has a
   component that fights the normal step. Confirmed fixed by testing against
   3 small NLPs with known closed-form solutions (equality, inequality, and
   bounds-only cases) plus a repeat with the SR1 Hessian mode -- all 4 now
   converge to the exact solution.
4. `p = p_n - p_t` is clipped component-wise so `x+p` respects bounds.

This is a legitimate simplification of a full active-set/composite-step
SQP method (similar in spirit to Byrd-Omojokun-style normal/tangential
decompositions), but it does **not** enforce general inequality-constraint
*bounds* exactly within the QP subproblem itself -- it relies on the outer
major SQP iterations (re-linearization + line search) to converge to
feasibility, which the test suite confirms works for these simple cases.
A rigorous active-set or interior-point QP solver remains a natural v2
enhancement (see §5).

**Convergence test**: also had to be corrected during testing -- a naive
`||g - J^T*lambda|| <= ktol` stationarity test never succeeds when a bound
is active (since `g` is nonzero there by design). Fixed with the standard
*projected-gradient* KKT test (as used by L-BFGS-B's `pgtol`): a nonzero
reduced-gradient component only counts as a violation if it points into
the feasible region at an active bound.

**BFGS update**: uses a simple curvature-condition skip rule
(`s^T y > eps`) rather than full Powell damping (which would need the
forward compact-BFGS product just to decide whether to damp) -- adequate
for the convex test problems tried so far; full damping is a possible
future refinement if non-convex problems prove to need it.

**Line search**: two modes are available, `sqpopt_linesearch_armijo`
(**default**) -- standard backtracking with an Armijo sufficient-decrease
test on the merit function, as used by default in `slsqp` -- and
`sqpopt_linesearch_exact` -- (approximate) exact 1-D minimization via
`fmin`. An exact line search is usually overkill (many extra function
evaluations for marginal benefit), so Armijo is the recommended default.

Switching the default to Armijo exposed several latent robustness bugs in
the v1 composite step that the (more forgiving) exact search had been
masking -- all fixed and covered by the test suite:

1. **Tiny-step Hessian corruption**: when a line search step is very small,
   the resulting `(s,y)` pair carries almost no reliable curvature
   information and can produce a huge/ill-conditioned `rho`, corrupting
   the Hessian approximation. Fixed by skipping the BFGS/SR1 update
   outright when `norm2(s) <= 1e-10`.
2. **Step overshoot**: the composite step can occasionally be much larger
   than reasonable. Fixed with a simple trust-region-style cap
   (`qp_solver%max_step`, default `10.0`): `p` is rescaled if
   `norm2(p)` exceeds it.
3. **Penalty parameter too small**: the \( \ell_1 \) exact penalty theory
   (Han/Powell) requires the penalty `mu` to exceed `||lambda||_inf` for
   the merit function's minimizer to coincide with the true constrained
   optimum; a too-small fixed penalty can make an infeasible point look
   better than the true solution. Fixed by adaptively updating
   `linesearch%penalty = max(penalty, max(|lambda|)+1)` every iteration
   (mirrors `slsqp`'s own multiplier-based penalty update).
4. **Non-descent step**: unlike a true QP solution, the v1 composite step
   has no guarantee of being a descent direction for the merit function.
   Fixed with an `slsqp`-style safeguard: if the (approximate) directional
   derivative `dot(g,p) - mu*viol(x) >= 0`, reset the Hessian to the
   identity and recompute `p` once before the line search.
5. **False convergence at inactive constraints**: with no complementarity
   enforcement, the least-squares multiplier estimate could "explain away"
   the gradient using an *inactive* constraint's direction whenever the
   gradient happened to be parallel to it, satisfying the KKT test at a
   non-optimal point. Fixed by adding a simple **active-set filter**
   (`qp_solver%active_tol`, default `1e-6`): only constraint rows at/beyond
   their bound are used for the multiplier estimate and null-space
   projection (the normal step still uses the full Jacobian, which is
   harmless since inactive rows have zero linearized violation).
   (Tried widening `active_tol` to `1.0` to help the HS71 nonlinear test
   below converge faster -- this *broke* `test_inequality_constrained`,
   which false-converged early because points well short of the boundary
   got flagged "active". Reverted to `1e-6`: the right tolerance is
   inherently problem-scale-dependent, and a small/exact one is safer.)
6. **Line search cascade**: `sqpopt_linesearch_armijo` now falls back to
   `sqpopt_linesearch_exact` if no backtracking step satisfies the Armijo
   test, rather than accepting a useless near-zero step.
7. **Second-order correction (SOC)**: for strongly nonlinear constraints,
   the linearized constraint prediction at `x+p` can differ enough from
   the true (nonlinear) value that a genuinely good step gets rejected by
   the merit function (the classic "Maratos effect"). `sqpopt_iterate_module`
   now computes an extra small correction step from the true constraint
   residual at `x+p` (another `LSQR` solve reusing the same Jacobian) and
   uses the corrected step if it improves the merit function value.

## 4. Test suite

[test/test_basic.f90](test/test_basic.f90) (`fpm test`) currently covers,
each against a known closed-form optimum:

1. `test_equality_constrained` — linear equality constraint.
2. `test_inequality_constrained` — linear inequality constraint (active at the solution).
3. `test_bounds_only` — a variable bound only, no general constraints (`m=0` edge case).
4. `test_sr1_hessian_mode` — repeats (1) with `hessian_mode = sqpopt_hessian_sr1`.
5. `test_exact_linesearch_mode` — repeats (2) with `linesearch_mode = sqpopt_linesearch_exact`.

All five pass with `istat = sqpopt_success` (clean convergence, not just
hitting `max_iter`).

[test/test_hs71.f90](test/test_hs71.f90) adds a genuinely **nonlinear** test:
Hock-Schittkowski problem 71 (nonlinear objective, one nonlinear equality and
one nonlinear inequality constraint, both simultaneously active at the
solution, plus variable bounds). **This originally exposed a real v1
limitation**: with the default `sqpopt_qp_composite` QP solver, the solver
settles into a small, stable oscillation (limit cycle) near the true
solution rather than converging tightly to it -- confirmed by running up to
3000 iterations with no further improvement, so it is not merely slow
convergence. `istat` ends as `sqpopt_max_iter_reached`, not `sqpopt_success`,
for that mode (and for the augmented Lagrangian merit / watchdog line search
variants, §6.1/§6.3, which only partially help), so those cases use a
generous `0.5` tolerance rather than requiring `sqpopt_success`. **Both
real QP solves fix this properly**: with `qp_solver_mode=sqpopt_qp_dense`
(§6.4) or `sqpopt_qp_reduced_hessian` (§6.2), `test_hs71` now reaches
`sqpopt_success` and requires it (tight `1e-4` tolerance) -- confirming
the long-standing hypothesis that a real QP solve, not another
merit-function/line-search patch, was the actual fix needed.

## 5. Backlog ("optional/advanced" work, deferred from v1)

- ~~A rigorous **active-set or interior-point QP solver** that enforces
  linearized general-constraint bounds exactly (the current composite-step
  heuristic relies on outer-iteration convergence instead).~~ **Implemented
  twice**: `sqpopt_qp_dense` (see §6.4/[DENSE_QP_PLAN.md](DENSE_QP_PLAN.md),
  a dense active-set QP) and `sqpopt_qp_reduced_hessian` (see
  §6.2/[REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md), a sparse/
  matrix-free projected-CG active-set QP) -- both get `test_hs71` to
  `sqpopt_success`. Use `sqpopt_qp_dense` for small-to-moderate problems
  (fewer major iterations needed, simpler direct linear algebra);
  `sqpopt_qp_reduced_hessian` for problems too large for `sqpopt_qp_dense`'s
  `O(n^2)`/`O(mn)` dense arrays (more major iterations needed, due to
  `LSQR`'s iterative rather than direct null-space projections).
- ~~A **smooth augmented Lagrangian merit function** as an alternative to the
  current \( \ell_1 \) merit function, to avoid the Maratos effect without
  needing the ad hoc second-order-correction patch.~~ **Implemented** as
  `sqpopt_merit_augmented_lagrangian` (see §6.1) -- but it did not clear up
  `test_hs71`'s limit cycle on its own (see the "Status" note in §6.1), so
  the SOC patch is still needed and the top-priority backlog item above
  (a real QP solve, §6.2) remains the higher-leverage fix for that case.
- A **watchdog line search** (`sqpopt_linesearch_watchdog`, see §6.3), from
  Powell's VF13 / the Chamberlain-Lemarechal-Pedersen-Powell watchdog
  technique: a best-point-so-far safety net plus a periodic relaxed
  acceptance criterion, targeting the same `test_hs71` limit cycle via a
  different mechanism than §6.1/SOC. **Implemented** -- gives the smallest
  `test_hs71` error of the three line-search/merit options tried so far,
  but (like §6.1) does not fully resolve the limit cycle on its own; see
  §6.3 for the full result, including an unrelated control-flow bug fix
  in `sqpopt_iterate_module` found and fixed along the way.
- A **dense QP solver option** (`sqpopt_qp_dense`, see §6.4), prompted by
  `slsqp` itself solving `test_hs71` to machine precision in 6 iterations
  using a dense BFGS + dense active-set QP. Design-in-progress, see
  [DENSE_QP_PLAN.md](DENSE_QP_PLAN.md) -- recommended as the faster,
  lower-risk path to validate that a real QP solve is what's needed,
  *before* investing in §6.2's sparse version.
- **`sqpopt_hessian_exact`** mode (user-supplied sparse Hessian of the
  Lagrangian) — currently falls back to BFGS in `sqpopt_iterate_module`.
- **Full Powell damping** for the BFGS update (currently a simpler
  curvature-condition skip rule).
- Whether to expose `NumDiff`'s automatic sparsity-pattern detection as a
  convenience path in `sqpopt_problem_module`, or keep sparsity patterns
  strictly user-supplied (current behavior).
- Default `lbfgs_memory` (currently 10), `qp_solver%max_step`
  (currently a fixed `2.0`, not problem-scale-aware), and
  `qp_solver%active_tol` (currently a fixed `1e-6`, also not problem-scale-
  aware -- see the note in §3) are none of them currently exposed on
  `sqpopt_options_type` either; worth revisiting together once there's a
  more rigorous QP solver to tune.
- A **trust-region globalization option**, as an alternative to the
  current line-search-only approach: adaptively re-solves the QP with a
  radius-tightened box around `x` (ratio-test/filter-based accept/reject
  of the radius itself) instead of only backtracking `alpha` along one
  fixed QP step. Design-in-progress, see
  [TRUST_REGION_PLAN.md](TRUST_REGION_PLAN.md) -- notably, this is also
  what unlocks the *literal* Fletcher & Leyffer filter-SQP combination
  (trust region + filter acceptance, §7), of which
  `sqpopt_linesearch_filter` is currently only a line-search adaptation.

## 6. SNOPT-family design ideas (`references/merit.pdf`, `sqdoc7.pdf`, `sndoc7.pdf`)

Three reference documents on the Stanford SOL group's SNOPT/SQOPT/NPSOL
family were added to `references/` (read via `pdftotext -layout`, not code
dependencies -- SNOPT/SQOPT are proprietary and only their PDF manuals /
papers are available here, not usable/linkable Fortran source):

- `merit.pdf` -- Gill, Murray, Saunders & Wright, *"Some Theoretical
  Properties of an Augmented Lagrangian Merit Function"* (SOL 86-6R), the
  paper describing the smooth merit function used in NPSOL/NPSQP (and, in
  spirit, SNOPT).
- `sqdoc7.pdf` -- the SQOPT 7 User's Guide (the large-scale sparse
  active-set QP solver used as SNOPT's QP engine).
- `sndoc7.pdf` -- the SNOPT User's Guide (large-scale SQP for nonlinear
  problems, built on SQOPT + the NPSOL-style merit function/line search).

Both ideas below are natural **user-selectable options** (following the
same `options%..._mode` pattern already used for `hessian_mode`/
`linesearch_mode`), not replacements for the v1
defaults -- consistent with "most common algorithms first, optional ones
later".

### 6.1 Smooth augmented Lagrangian merit function (from `merit.pdf`)

NPSOL's `NPSQP` algorithm treats the Lagrange multiplier estimate \(
\lambda \) as an extra variable (not just a by-product of the QP), and adds
non-negative slacks \( s \) for the inequality constraints so everything
can be included smoothly in the linesearch:

$$ L(x,\lambda,s,\rho) = f(x) - \lambda^T\!\left(c(x)-s\right) + \tfrac{1}{2}\rho \lVert c(x)-s \rVert_2^2, \quad s \ge 0 $$

Key mechanics (Section 2-3 of the paper):

- The search direction for \( \lambda \) is \( \xi = \mu - \lambda \) (where
  \( \mu \) is the QP's own multiplier), so a full step (\( \alpha=1 \))
  sets \( \lambda \to \mu \) exactly.
- The slacks are set in closed form each iteration:
  \( s_i = \max(0, c_i - \lambda_i/\rho) \) (or \( \max(0,c_i) \) if \( \rho=0 \)).
- The penalty parameter \( \rho \) is increased **only when needed** to
  guarantee \( \phi'(0,\rho) \le -\tfrac{1}{2}p^THp \) (a closed-form
  threshold \( \hat\rho \), Lemma 4.3) -- this is a more principled version
  of our current "`penalty = max(penalty, max|lambda|+1)`" heuristic
  (Han/Powell rule for the non-smooth \( \ell_1 \) merit function).
- Because \( L \) is **twice continuously differentiable** (unlike the
  \( \ell_1 \) merit function's kinks at the constraint boundaries), the
  line search can use a Wolfe-type test with both a value condition (3.1a)
  and a *derivative* condition (3.1b) -- and the paper proves this avoids
  the Maratos effect and allows full steps near the solution, which is
  **exactly the failure mode our second-order-correction patch was added to
  work around** for `test_hs71`. A correctly implemented augmented
  Lagrangian merit function would likely make that patch unnecessary.

**Proposed v2 design**: add `sqpopt_merit_l1` (current, default) and
`sqpopt_merit_augmented_lagrangian` (new) as a `merit_mode` option on
`sqpopt_linesearch_type`/`sqpopt_options_type`. The augmented Lagrangian
mode would need: (a) slack variables `s` (size `m`, only for inequality
rows), (b) the multiplier `lambda` promoted to real linesearch state
(already partially true -- we carry `lambda` across iterations, just not
as a linesearch variable with its own search direction `xi`), and (c) the
closed-form `rho` update above in place of the current heuristic.
**Caveat**: the paper's global convergence theory assumes the QP (1.3) is
solved to KKT optimality, giving a well-defined multiplier `mu` satisfying
(1.4); our v1 composite-step QP only produces an *approximate* least-squares
`lambda`, not a true QP multiplier -- so this pairs best with §6.2 (a real
QP solve), not as a drop-in replacement for the current merit function alone.

**Status: implemented** as `sqpopt_merit_augmented_lagrangian` in
[sqpopt_linesearch_module](src/sqpopt_linesearch_module.f90), with two
deliberate simplifications relative to the full NPSQP theory above (both
consistent with the caveat just noted, and clearly commented in the code):

- The slack `s` generalizes \( s_i=\max(0,c_i-\lambda_i/\rho) \) to
  `sqpopt`'s two-sided bounds: \( s = \text{clip}(c-\lambda/\rho,\, c_l,\,
  c_u) \) (equality rows are handled automatically, since `s` is then
  clipped to the single value `c_lb=c_ub` regardless of `lambda`/`rho`).
- `lambda`/`s` are recomputed fresh at each trial point during the line
  search (using the *current* `lambda` estimate, held fixed) rather than
  being advanced along a joint `(x,lambda,s)` step with its own search
  direction `xi`/`q` -- much simpler, at the cost of not exactly matching
  the paper's line-search theory.
- The penalty parameter reuses the same `penalty = max(penalty,
  max|lambda|+1)` heuristic as `sqpopt_merit_l1`, rather than the paper's
  closed-form `rho-hat` threshold (Lemma 4.3) -- implementing that exactly
  requires the QP's own multiplier `mu` (distinct from `lambda`) and the
  `xi`/`q` machinery above, which pairs better with §6.2 (a real QP solve).

`test/test_basic.f90`'s `test_augmented_lagrangian_merit` confirms the new
mode converges correctly on the existing inequality-constrained test.
**However, trying it on `test_hs71` did *not* clearly fix the limit-cycle**
(similar residual error to the default `sqpopt_merit_l1` + SOC combination)
-- consistent with the caveat above: without a real QP multiplier and the
full joint-step line search, smoothness alone doesn't fully deliver the
paper's Maratos-avoidance guarantee. §6.2 (a real QP solve) remains the
higher-leverage next step for `test_hs71` specifically.

### 6.2 Reduced-Hessian active-set QP (from `sqdoc7.pdf`, SQOPT)

SQOPT solves (convex) QPs of the same two-sided-bounds form we already use
(`l <= (x, Ax) <= u`, i.e. exactly `x_lb<=x<=x_ub` and `c_lb<=c(x)<=c_ub`
after linearizing `c`), via a **reduced-Hessian active-set / reduced-gradient
method**. Two things make it a strong architectural fit for `sqpopt`:

- **The Hessian is never formed**: SQOPT requires only a user subroutine
  `qpHx(x) -> Hx` that returns the Hessian-vector product -- exactly our
  `sqpopt_hessian_type%hv_product`. A reduced-Hessian method only ever
  needs `Hv` products restricted to the *current null-space basis* `Z`
  (i.e. `Z^T H Z v`), which composes naturally with our existing
  `hv_product`.
- **The constraint matrix is sparse** (their CSC `Acol`/`indA`/`locA`
  triplet, vs. our COO `sqpopt_sparse_matrix` -- trivial to convert
  between, or just build both representations from the same problem data).
- It uses the **same slack-variable device** as SQOPT's `sndoc7.pdf`
  algorithm and our own bound handling: general constraints become
  equalities `Ax - s = 0` with bounds moved onto `s`, unifying variable and
  constraint bounds into one bounded-variable framework -- this is
  basically the qpOASES-style "generalized bounds" idea considered (and set
  aside for complexity) back in the original architecture discussion.
- SQOPT's **elastic mode** (relaxing bounds with a penalty when a linearized
  QP subproblem would otherwise be infeasible) is directly relevant to a
  known soft spot in our own design: our v1 QP heuristic side-steps this by
  never enforcing bounds exactly inside the QP at all, relying on the outer
  iteration instead. A proper elastic-mode active-set QP would let us
  enforce bounds properly *and* always have a feasible subproblem.

**Proposed v2 design**: a new `sqpopt_qp_solver_mode` option,
`sqpopt_qp_composite` (current v1 heuristic, default) vs.
`sqpopt_qp_reduced_hessian` (new): a genuine active-set QP working on the
same `x_lb<=x+p<=x_ub`, `c_lb<=c+Jp<=c_ub` subproblem, using `hv_product`
for all Hessian-vector products (never forming `H`), `sparse_matvec`/
`sparse_matvec_transpose` for `J`, and `lusol` (already a dependency, unused
by the v1 QP path) to factorize the small working-set systems that arise as
constraints/bounds enter and leave the active set -- this is the concrete
version of the "rigorous active-set QP solver" already flagged as the
top-priority backlog item in §5, and (combined with §6.1) is expected to fix
the `test_hs71` limit-cycle behavior properly, since a real QP solve
guarantees its solution is a descent direction for a correctly-parameterized
merit function (the theoretical property Lemma 4.1(a) relies on, and which
our heuristic composite step cannot guarantee -- see the safeguards in §3).

**Status: implemented.** A detailed, staged implementation plan was
written up in [REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md)
and then implemented as designed:
- New [sqpopt_qp_reduced_hessian_module](src/sqpopt_qp_reduced_hessian_module.f90)
  (`sqpopt_reduced_hessian_qp_type`) -- a **projected-conjugate-gradient**
  active-set QP (Gould, Hribar & Nocedal 1998). `lusol` ends up unused
  here: null-space projections are obtained by re-solving a small
  least-squares problem with `LSQR` each time (`project_null`, exactly
  the technique v1's own composite step already uses for its tangential
  step), not by factorizing a maintained basis with `lusol` as the
  original proposal above suggested -- see REDUCED_HESSIAN_QP_PLAN.md §3
  for why the maintained-basis approach was judged too large/risky to
  build from scratch. Works directly in `p`-space (general constraints
  and bounds unified as `m+n` two-sided rows), the same simplification
  made in [DENSE_QP_PLAN.md](DENSE_QP_PLAN.md)'s dense sibling (§6.4)
  rather than the `w=(p,s)` slack padding originally proposed above.
- Validated in isolation first
  ([test/test_qp_reduced_hessian.f90](test/test_qp_reduced_hessian.f90)),
  the exact same 4 hand-verified QPs used for the dense solver's tests --
  all passed exactly on the first attempt (matching `p` and Lagrange
  multiplier signs), a good sign the active-set control logic (shared in
  spirit with the dense solver, see the comparison table in
  DENSE_QP_PLAN.md §6) transferred correctly to the sparse/LSQR backend.
- Wired in via the same `mode` field on `sqpopt_qp_solver_type`
  (`sqpopt_qp_reduced_hessian=3`, alongside `sqpopt_qp_composite=1` and
  `sqpopt_qp_dense=2`) plus `options%qp_solver_mode`.
- **Result on `test_hs71`**: also reaches `sqpopt_success`, converging to
  within `7e-7` -- matching the dense solver's result while staying fully
  sparse/matrix-free. It needs noticeably more major SQP iterations to
  get there, though (`test_hs71`'s shared `max_iter` was raised to `3000`
  to accommodate it; the dense solver succeeds well within `300`) -- the
  cost of `LSQR`'s iterative tolerances vs. the dense solver's one-shot
  direct factorizations, precisely the tradeoff the plan called out in
  advance (§3/§9 there). For problems too large for `sqpopt_qp_dense`'s
  `O(n^2)`/`O(mn)` dense arrays, this is the mode that scales.
- **`LSQR` tuning matters a lot, and is exposed as user-settable fields**
  (`lsqr_atol`/`lsqr_btol`/`lsqr_conlim`/`lsqr_itnlim` on
  `sqpopt_reduced_hessian_qp_type`) rather than hard-coded: by default,
  every `LSQR` solve in this module passes `atol=btol=0`, which `LSQR`
  itself treats as "use machine precision" (see `lsqr.f90`) -- tighter
  than necessary and part of why so many major iterations were needed
  above. A sweep on `test_hs71` (using the function-call counters
  `i_obj`/`i_grad`/`i_cons`/`i_jac` added to `test_hs71.f90` as the metric
  to reduce, not just major-iteration count) found `atol=btol=5e-10` cuts
  `i_obj` from `8484` to `1237` (~85% fewer function calls) while still
  reaching `sqpopt_success` to the same accuracy. **This is not a safe new
  default, though**: the sweep is sharply non-monotonic near the
  boundary -- `1e-10` helps a little (`6573`), `5e-10` helps a lot
  (`1237`), but `8e-10` already breaks convergence entirely, and `1e-9`/
  `2e-9` are worse still (one run even produced `NaN`). The "sweet spot"
  is real but narrow and almost certainly problem-specific (a different
  active-set path gets taken at different precisions), so it's kept as a
  documented, user-tunable knob and a worked example in
  `test/test_hs71.f90` (`'rh, tuned LSQR (atol=btol=5e-10)'`), not a
  library-wide default.

### 6.4 Dense QP solver option (from comparing against `slsqp` directly)

Prompted by adding `test/slsqp_test_71.f90` alongside `test_hs71.f90`
(same underlying HS71 problem, confirmed by direct comparison): `slsqp`
converges to machine precision in 6 iterations, where `sqpopt`'s v1
heuristic only manages a loose limit cycle. `slsqp` does this with an
entirely **dense** BFGS Hessian + dense active-set-style QP, which
confirmed (by reading `slsqp_core.f90` directly) that its own
constraint-Jacobian and Hessian-Cholesky-factor arguments are genuinely
dense (`dimension(la,n+1)`/packed dense Cholesky factor) -- there is no
sparse entry point to reuse, consistent with `PLAN.md` §2's note that
`slsqp` is "algorithmic inspiration only, not linked as a dependency."

**Status: implemented.** A staged plan for an *opt-in* `sqpopt_qp_dense`
mode -- forms dense `J`/`H` from the existing sparse/matrix-free
representations each iteration, and reuses the **same** active-set
control logic as §6.2/[REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md)
§6, just with a dense QR/modified-Cholesky backend instead of sparse
`LSQR`/projected-CG -- was written up in
[DENSE_QP_PLAN.md](DENSE_QP_PLAN.md) and then implemented as designed:
- New [sqpopt_dense_linalg_module](src/sqpopt_dense_linalg_module.f90)
  (`dense_null_space` via Householder QR, `dense_modified_cholesky`,
  `dense_solve_cholesky`) and
  [sqpopt_qp_dense_module](src/sqpopt_qp_dense_module.f90) (the dense
  active-set QP itself, working directly in `p`-space rather than the
  `w=(p,s)` padding described in the plan -- general constraints and
  variable bounds are just treated uniformly as `m+n` two-sided rows on
  `p`, which is mathematically equivalent and simpler to implement).
- Validated in isolation first ([test/test_qp_dense.f90](test/test_qp_dense.f90),
  4 hand-verified small QPs including a two-simultaneously-active-
  constraints case), all matching hand-derived KKT solutions exactly
  (including Lagrange multiplier signs) on the first attempt.
- Wired in via a `mode` field on `sqpopt_qp_solver_type`
  (`sqpopt_qp_composite`/`sqpopt_qp_dense`) plus
  `options%qp_solver_mode`, exactly as planned -- no changes needed to
  `sqpopt_iterate_module` or `sqpopt_module`'s public API.
- **Result on `test_hs71`**: with `qp_solver_mode=sqpopt_qp_dense`, the
  solver now reaches **`istat=sqpopt_success`**, converging to within
  `2e-7` of the known solution -- essentially matching `slsqp`'s own
  6-iteration machine-precision result on the same problem, and fully
  resolving the limit cycle that the v1 composite step, the augmented
  Lagrangian merit function (§6.1), and the watchdog line search (§6.3)
  could each only partially work around. `sqpopt`'s test suite now
  enforces this tightly (`test_hs71`'s `dense QP` case requires
  `sqpopt_success` and `1e-4` accuracy, unlike the other three modes'
  loose `0.5` tolerance).
- One unrelated, pre-existing latent bug was found and fixed along the
  way (via `-fcheck=all`, after an initial segfault): `exact_line_search`'s
  nested `merit_along_direction` function declared its trial constraint
  array `c_trial` with `dimension(size(x))` (`n`) instead of `dimension
  (size(c_lb))` (`m`) -- harmless when `n==m`, but a real out-of-bounds
  array access whenever `n/=m` (as in `test_hs71`, `n=4`,`m=2`). Always
  worth re-running the full suite with `-fcheck=all` after any change
  that touches shape-sensitive array code.

### 6.3 Watchdog line search (from `references/vf13`, Powell's VF13 / HSL archive)

`references/vf13` contains the HSL archive package spec (`vf13_Fortran.pdf`)
plus its Fortran 77 source (`vf13d.f`/`vf13s.f`, `ddeps.f`/`sdeps.f`) for
**VF13**, M.J.D. Powell's variable-metric SQP method -- a direct ancestor of
the `slsqp` family already in `[dev-dependencies]`. Two ideas stand out:

- **`VE17AD`: a dense Goldfarb-Idnani dual active-set QP solver.** This is
  what VF13 calls each iteration to solve the linearized QP subproblem
  (BFGS Hessian `B`, linearized constraints). It's a smaller, simpler
  cousin of SQOPT's reduced-Hessian method (§6.2) -- useful as a reference
  for active-set bookkeeping (how constraints/bounds enter and leave the
  working set, the KKT system update on each change) -- but its workspace
  is dense (`B` is `n x n`, `CN` is `(n+1) x m`, cost `~5n^2/2`), which
  conflicts with `sqpopt`'s "never form a dense `n x n`" design constraint.
  Like `slsqp`/`psqp`, this is algorithmic inspiration only, not portable
  as-is; §6.2's sparse/matrix-free design remains the right target.
- **The "watchdog technique"** (Chamberlain, Lemarechal, Pedersen & Powell,
  *Math. Prog. Study 16* (1982)) -- a line-search relaxation mechanism
  aimed at **exactly** `test_hs71`'s failure mode: when constraint
  boundaries are curved and several are simultaneously active, a strict
  merit-function line search can force the iterates to hug the boundary
  and zigzag (the Maratos effect) instead of taking the good, nearly-full
  quasi-Newton step that would actually make progress. VF13BD's mechanism
  (traced through `vf13d.f` lines ~230-450):
  - Tracks the **best point seen so far** (`XOPT`/`FOPT`/`WOPT`, using an
    \( \ell_1 \)-style merit `W = F + sum_k VMU(k)*violation(c_k)`, with
    **per-constraint** penalty weights `VMU(k)` that are only ever
    increased, never decreased -- more granular than our single scalar
    `penalty`).
  - On most iterations it uses the standard sufficient-decrease test (like
    our Armijo/`sqpopt_linesearch_armijo`). But for up to `NWDOGT=2`
    consecutive iterations after the merit value has been "good enough"
    once (`ISWDOG=0`, the *relaxed* criterion), it accepts steps using a
    much weaker test based on the actual Lagrangian value rather than the
    merit function -- allowing the merit function to **temporarily get
    worse** so the iteration isn't trapped hugging a curved boundary.
  - If, after those relaxed iterations, the merit value still hasn't beaten
    `WOPT`, it **backtracks** all the way to `XOPT` (the safety net) and
    forces the strict criterion for the next `NWDXXX=10` iterations (a
    "cooldown" before trying relaxed acceptance again).
  - This is a different, well-established mechanism for the same problem
    our second-order-correction patch (§3) and the augmented Lagrangian
    merit function (§6.1) both target, and it doesn't require a true QP
    multiplier to work (unlike the AL merit function's full theoretical
    guarantee) -- it only needs a merit value and a "best point so far",
    both of which the v1 composite-step method already has.

**Proposed v3 design**: a new `sqpopt_linesearch_watchdog` mode, additive
to (not a replacement for) the existing `armijo`/`exact` modes and
orthogonal to `merit_mode`:
  - New persistent state carried across major iterations (alongside the
    existing `x_prev`/`gl_prev` quasi-Newton state in
    `sqpopt_iterate_module`): best-point-so-far `x_opt`/`f_opt`/`w_opt`, a
    relaxed-mode countdown, and a cooldown counter.
  - Each iteration: accept the full QP step immediately if it improves on
    `w_opt`; otherwise fall back to the current strict Armijo test *unless*
    the relaxed window is open, in which case accept a much weaker test for
    up to 2 iterations before forcing a backtrack to `x_opt` and starting a
    10-iteration cooldown.
  - Worth trying on `test_hs71` specifically, since it targets that exact
    limit-cycle failure mode with a different (and historically effective)
    mechanism than what's already been tried in §6.1.

**Status: implemented**, as `sqpopt_linesearch_watchdog` in
[sqpopt_linesearch_module](src/sqpopt_linesearch_module.f90) (a simplified
version of the design above; no per-constraint `VMU` weights, just the
single scalar `penalty` shared with the other modes), together with two
other changes made along the way:

- **A real (and unrelated) bug fix in `sqpopt_iterate_module`**: the major
  iteration used to do `if (istat /= sqpopt_success) return` right after
  the line search call, which *skipped* the `x = x + alpha*p` update
  whenever `armijo_line_search` reported `sqpopt_line_search_failed` --
  even though that mode's own doc comment says it "accepts `alpha_min`
  anyway". Since the next major iteration would then recompute the
  *identical* `p` from the *identical*, unchanged `x`, a single failed
  line search could silently freeze the whole rest of a run. Fixed by
  always applying the line search's returned point. **Empirically, on
  `test_hs71`, this fix alone made the (already loose) `sqpopt_merit_l1`/
  `sqpopt_merit_augmented_lagrangian` results *worse*** (x-error grew from
  ~0.15-0.24 to ~0.42-0.45) -- the old "freeze" bug had apparently been
  acting as an accidental stabilizer, freezing at a reasonably good point
  instead of continuing to wander. This is a good illustration of why
  `test_hs71` matters: it is a much more sensitive probe of small control-
  flow changes than `test_basic`'s simpler problems (which were and remain
  unaffected either way).
- **The `line_search`/`armijo_line_search`/`exact_line_search` interface
  grew an `x_new` output** (the actual next point), instead of only
  `alpha` (with the caller always computing `x + alpha*p` itself) --
  needed so that `sqpopt_linesearch_watchdog`'s backtrack can return an
  *arbitrary earlier point* (`x_opt`, from a previous major iteration, not
  reachable as `x + alpha*p` along the *current* search direction `p`).
- The watchdog's "reward" for opening a fresh relaxed-acceptance window
  was tuned to only trigger after a **genuinely good, (nearly) full**
  standard-test-passing step (`alpha >= 0.99`) -- mirroring VF13's
  `ISWDOG` logic, which is gated on the predicted merit reduction, not
  just any improvement. An earlier, more permissive version (reward on
  *any* improvement, however tiny) made `test_hs71` noticeably worse (the
  relaxed window was open almost constantly, so it behaved close to
  "always accept the full step", which just wandered more).

**Result on `test_hs71`** (all three modes still end at
`istat=sqpopt_max_iter_reached`, none reach `sqpopt_success`): the tuned
watchdog gives the *smallest* max-component error of the three options
tried so far (`~0.32`, vs `~0.42` for `sqpopt_merit_l1` and `~0.45` for
`sqpopt_merit_augmented_lagrangian`, all measured after the bug fix above)
-- a genuine, if partial, improvement, and the best evidence yet that a
better line-search/step-acceptance strategy alone can move the needle on
this benchmark. It still doesn't fully resolve the limit cycle, though --
consistent with the recurring conclusion across §6.1 and this section that
a real QP solve (§6.2) is needed for `sqpopt_qp_solver_module`'s composite
step to have the guaranteed-descent property these line-search techniques
assume.

## 7. Filter method (Fletcher & Leyffer, from `references/fletcher.pdf`)

`references/fletcher.pdf` is Fletcher & Leyffer, *"Nonlinear programming
without a penalty function"* (Math. Program. 91 (2002)) -- the paper that
introduced the **filter** concept as a globalization strategy for SQP that
dispenses with a merit function/penalty parameter entirely. Instead of
combining the objective `f` and constraint violation `h(c(x))` into a
single scalar via a penalty parameter (with the well-known difficulties of
choosing that parameter -- too small and the penalty function may not have
a local minimum at the solution, too large and it damps out the objective
near a curved constraint boundary), a trial point is accepted if its
`(f,h)` pair is not *dominated* by any previously-accepted iterate's
`(f,h)` pair: i.e. it is better than every prior accepted point in *at
least one* of the two objectives, akin to a Pareto/multi-objective
acceptance test (the "filter"). A sufficient-reduction envelope (their
eqs. 3-4, using each filter entry's own QP-predicted objective decrease
`q` and an estimated penalty-parameter scale `mu`) excludes points that
are only trivially better than an existing entry, preventing cycling.

The paper's full algorithm (their Algorithm 3) is **trust-region**-based:
an \( \ell_\infty \)-norm-bounded QP subproblem, with a rejected step
handled by shrinking the trust-region radius, and an infeasible QP
triggering a whole separate **feasibility restoration phase** (a nested
SQP-like iteration minimizing `h(c(x))` with its own "phase I filter"),
plus North-West/South-East filter corner rules (which the paper itself
notes, §3.6, can later be dispensed with). `sqpopt` is a **line-search**
SQP, not trust-region, so a literal port isn't a good fit -- the natural
line-search analogue of "shrink the trust-region radius and retry" is
"backtrack `alpha` and retry" (this is essentially what IPOPT's filter
line-search method, Wächter & Biegler 2005, does relative to the original
trust-region filter-SQP).

**Status: implemented**, as `sqpopt_linesearch_filter` in
[sqpopt_linesearch_module](src/sqpopt_linesearch_module.f90) -- a scoped
port that keeps the core filter concept (domination test + eqs. 3-4
sufficient-reduction envelope + the §3.2 upper bound on `h`) but, matching
the paper's own permission to simplify (§3.6), omits the NW/SE corner
rules, and omits the restoration phase entirely (infeasible-QP handling is
left to the existing outer SQP machinery/`istat` codes, same as every
other line-search mode -- none of `sqpopt`'s three QP solver modes
currently distinguish "infeasible QP" from other failure modes as a
trigger for a dedicated restoration phase; this would be a substantial
follow-on project, not a line-search-module change). Design notes:

- The filter stores four parallel arrays (`filter_f`/`filter_h`/
  `filter_q`/`filter_mu`), one entry appended per accepted iterate,
  dominated entries removed on each acceptance -- directly mirroring the
  paper's "later on four" scalars-per-entry design (§2).
- `q` (the QP's predicted decrease in `f`, \( q=-(g^Tp+\tfrac12 p^THp)
  \)) needs the Hessian, which `sqpopt_linesearch_module` doesn't have
  access to -- computed in `sqpopt_iterate_module` (which does) via
  `hessian%hv_product(p,hp)` and threaded through as a new argument to
  `linesearch%search`/`line_search`. `mu` (the per-entry penalty-parameter
  estimate, §3.5: least power of ten larger than \( \lVert\lambda
  \rVert_\infty\), clipped to \( [10^{-6},10^6] \)) needs no new plumbing
  since `lambda` was already passed to `line_search`.
- **Found via testing, an important edge case not obvious from the paper
  (which assumes `m>0` throughout)**: if `h` stays at/near zero for both
  the current and every trial point (an unconstrained problem, or simply
  an already-fully-feasible run of iterates), the domination test alone
  is *vacuous* -- eq. (3) (`h_trial <= beta*h_l`) is trivially satisfied
  whenever `h_trial` and `h_l` are both ~0, regardless of `f`, so *every*
  point would be "accepted" with no globalization on `f` at all. Fixed
  with an explicit fallback (`filter_feas_tol`, default `1e-8`): when both
  the current and trial point are below this violation threshold, plain
  monotonic descent in `f` is also required. Verified with a dedicated
  test (`test_filter_linesearch_mode_equality` in `test/test_basic.f90`,
  a feasible-starting-point equality-constrained problem where `h`stays
  ~0 throughout) that this fallback is actually exercised and needed.
- SOC (`second_order_correction` in `sqpopt_iterate_module`) and the
  "is `p` a descent direction" safeguard both needed **zero** changes: they
  already operate purely on `linesearch%eval_merit`/`merit_mode`, which is
  completely independent of `linesearch%mode` -- the filter mode just
  never uses the resulting merit *value* as its own acceptance test, but
  is happy to let those two mechanisms keep using it internally as a
  ranking/safeguard heuristic on `p` before the filter-based line search
  ever runs.
- **Result on `test_hs71`** (paired with `sqpopt_qp_dense`, the same real
  active-set QP the other two `sqpopt_success`-reaching modes use):
  reaches `sqpopt_success`, x-error `~1.9e-4` -- looser than the Armijo/
  watchdog+dense-QP combination's `~1e-6`-`3e-5` (needs a relaxed `5e-4`
  tolerance in `test_hs71.f90`'s check, vs `1e-4` for the others) and
  needs substantially more function evaluations (`i_obj~3200` vs
  `~26-160` for the other dense/reduced-Hessian runs) -- consistent with
  the filter method's more permissive, non-monotone acceptance test
  producing a noisier endgame without the paper's own SOC-integrated-
  into-filter refinement (§3.1, which this port omits, reusing the
  existing merit-based SOC instead -- see above). Still a genuine, useful,
  independent globalization strategy, just not (yet) tuned to match the
  other modes' precision on this particular hard benchmark.

