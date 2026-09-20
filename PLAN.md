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
| `sqpopt_qp_solver_module` | v1 composite-step direction finder (see §3) | done (simplified) |
| `sqpopt_linesearch_module` | \( \ell_1 \) merit function + Armijo (default) / `fmin`-exact line search | done |
| `sqpopt_convergence_module` | projected-gradient KKT test + feasibility | done |
| `sqpopt_iterate_module` | orchestrates one major SQP iteration | done |
| `sqpopt_module` | public `sqpopt_type` facade (`initialize`/`set_problem`/`set_options`/`solve`/`get_solution`/`destroy`) | done |

Key design decision, upheld throughout: **no dense `n×n` or `m×n` arrays,
ever**. The Jacobian is sparse COO triplets with a fixed sparsity pattern;
the Hessian approximation is a matrix-free limited-memory operator (only
`O(n * lbfgs_memory)` storage, `lbfgs_memory` a small constant).

Deliberately **not yet implemented** (see §5 for the "optional/advanced"
backlog): `sqpopt_hessian_exact` mode (falls back to BFGS), a rigorous
active-set/interior-point QP solver (v1 uses a simplified composite-step
heuristic instead, see §3), `ftol`/`xtol` no-progress stopping tests, and
diagnostic printing (`options%print_level` is currently unused).

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
  wired up in [sqpopt_linalg_module](src/sqpopt_linalg_module.f90)
  (`solve_sparse_linear_system`, dispatched by `linear_solver_mode`) as
  general-purpose sparse square-system solvers, for future use by a more
  rigorous QP solver (not yet called by the v1 composite-step QP solver,
  which only needs `LSQR`'s rectangular least-squares capability).

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

## 5. Backlog ("optional/advanced" work, deferred from v1)

- A rigorous **active-set or interior-point QP solver** that enforces
  linearized general-constraint bounds exactly (the current composite-step
  heuristic relies on outer-iteration convergence instead).
- **`sqpopt_hessian_exact`** mode (user-supplied sparse Hessian of the
  Lagrangian) — currently falls back to BFGS in `sqpopt_iterate_module`.
- **Full Powell damping** for the BFGS update (currently a simpler
  curvature-condition skip rule).
- **`ftol`/`xtol`** no-progress stopping tests (only the KKT/feasibility
  test is implemented).
- **Diagnostic printing** (`options%print_level` is defined but unused).
- Nonlinear test problems with curvature (the current test suite is all
  convex quadratic objectives with linear constraints); a good next step is
  porting an actual Hock-Schittkowski problem (e.g. HS71, already used by
  `slsqp`'s and `psqp`'s own test suites) to exercise nonlinear Jacobians
  and multiple BFGS/SR1 updates per solve.
- Whether to expose `NumDiff`'s automatic sparsity-pattern detection as a
  convenience path in `sqpopt_problem_module`, or keep sparsity patterns
  strictly user-supplied (current behavior).
- Default `lbfgs_memory` (currently 10) and default `linear_solver_mode`
  (currently `lusol`, though not yet exercised by the v1 QP path) —
  reasonable defaults, worth revisiting once benchmarked on larger problems.

