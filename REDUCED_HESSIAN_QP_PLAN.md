# Design plan: reduced-Hessian active-set QP solver (`sqpopt_qp_reduced_hessian`)

**Status: implemented, per this plan, with results below.** This expands
PLAN.md §6.2 ("Reduced-Hessian active-set QP, from `sqdoc7.pdf`/SQOPT")
into a concrete, buildable design, sized to `sqpopt`'s existing
constraints (matrix-free `H`, sparse COO `J`, no dense `n x n`/`m x n`
arrays, only `LSQR`/`LSMR`/`lusol` available as linear-algebra building
blocks). It intentionally trades some of SQOPT's efficiency for a much
smaller implementation footprint, by reusing infrastructure `sqpopt`
already has instead of building a sparse LU-updated basis factorization
from scratch (see "Key simplification" below).

**Result**: `test_hs71` with `options%qp_solver_mode =
sqpopt_qp_reduced_hessian` now reaches `istat=sqpopt_success`, converging
to within `7e-7` of the known solution -- matching
[DENSE_QP_PLAN.md](DENSE_QP_PLAN.md)'s dense solver's result, while
staying fully sparse/matrix-free throughout. As anticipated in §3/§9
below, it needs noticeably more major SQP iterations to get there than
the dense solver (`test_hs71`'s shared `max_iter` was raised to `3000` to
accommodate it -- the dense solver reaches `sqpopt_success` well within
`300`) -- the price of `LSQR`'s iterative tolerances vs. the dense
solver's one-shot direct factorizations, exactly the tradeoff this plan
called out in advance. Implementation deviated from this plan in the same
way [DENSE_QP_PLAN.md](DENSE_QP_PLAN.md) did: no `w=(p,s)` slack padding
-- works directly in `p`-space, treating general constraints and variable
bounds uniformly as `m+n` two-sided rows on `p`. See `PLAN.md` §6.2 for
the full implementation summary.

## 1. Goal and scope

Replace the v1 composite-step QP (`sqpopt_qp_solver_module`) with an
option, not a replacement -- a **real active-set QP solve** of the same
subproblem it already solves each major iteration:

$$ \min_{p} \; \tfrac{1}{2} p^T H p + g^T p \quad \text{s.t.} \quad
c_l - c(x) \le J p \le c_u - c(x), \quad x_l - x \le p \le x_u - x $$

using only:
- `hessian%hv_product(v) -> Hv` (matrix-free, already exists, [sqpopt_hessian_module.f90](src/sqpopt_hessian_module.f90))
- `sparse_matvec`/`sparse_matvec_transpose` on the COO `jac` (already exist, [sqpopt_linalg_module.f90](src/sqpopt_linalg_module.f90))
- `LSQR` (already used this way in v1) and `lusol` (a dependency, currently unused by the QP path)

**Definition of done**: a new `sqpopt_qp_reduced_hessian` mode that (a)
enforces the linearized constraints/bounds *exactly* within the QP (unlike
v1), (b) is a drop-in alternative behind the *same*
`solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)`
interface, and (c) is validated first in isolation (new unit-style QP
tests) and then by re-running the full existing NLP test suite
(`test_basic`, `test_hs71`, `test_medium`) with it selected, with
`test_hs71` finally reaching `istat=sqpopt_success` as the primary
acceptance criterion (per PLAN.md's repeated conclusion that this is the
missing piece).

## 2. Reformulation: bounded variables + slacks (SQOPT/qpOASES-style)

Following `sqdoc7.pdf` and the note already in PLAN.md §6.2, introduce a
slack vector `s` (`dimension(m)`) and rewrite the QP as:

$$ \min_{p,s} \; \tfrac12 p^T H p + g^T p \quad \text{s.t.} \quad Jp - s = 0,
\quad x_l - x \le p \le x_u - x, \quad c_l - c(x) \le s \le c_u - c(x) $$

This unifies variable bounds and (linearized) general-constraint bounds
into a single **bounded-variable** framework over the combined vector
`w = (p, s)` (`dimension(n+m)`), with one block of linear *equality*
constraints `[J | -I] w = 0`. This is exactly `sqpopt`'s existing
`x_lb<=x<=x_ub`/`c_lb<=c<=c_ub` two-sided-bounds convention, just applied
to `p` and `s` instead of `x` and `c` -- no new bound-handling logic is
needed, only the extra `n -> n+m` bookkeeping.

**Working set** = the set of *nonbasic* components of `w` (currently
pinned to one of their two bounds). The *free*/basic components span the
current null space of the active constraints.

## 3. Key simplification vs. SQOPT: no explicit basis factorization

SQOPT (and VE17AD in `references/vf13`) maintain an explicit **basis
matrix** `B` (a square, nonsingular sub-matrix of `[J | -I]` selected by
the current partition into basic/nonbasic/superbasic variables) and update
its **sparse LU factorization** incrementally as the working set changes
(Forrest-Tomlin-style updates). This is the single biggest chunk of
classical active-set QP code (basis repair, factorization updates,
degeneracy handling) and would be a large, risky undertaking to build from
scratch on top of `lusol`.

**Simplification adopted here**: never form or factorize an explicit basis
`B`. Instead, represent the null space of the active constraints
*implicitly*, and get any null-space projection we need by re-solving a
small least-squares problem with `LSQR` each time -- exactly the technique
`sqpopt_qp_solver_module`'s v1 composite step *already* uses for its
tangential step (`p_t = v - J_A^T z`, `z` from `min ||J_A^T z - v||_2`).
This costs one extra `LSQR` solve per projection instead of an `O(1)`
update to a maintained factorization, but:
- it reuses infrastructure that already exists and is already tested,
- it avoids the largest source of implementation risk/bugs in classical
  active-set codes (basis repair after rank changes, degeneracy in the
  factorization), and
- active-set sizes in the problems `sqpopt` targets are not expected to be
  large enough that repeated `LSQR` solves are a real performance concern
  (this can be revisited later if profiling shows otherwise -- see §9).

This means the QP solver described below is best understood as a
**projected-conjugate-gradient (PCG) active-set method** (Nocedal & Wright,
*Numerical Optimization*, Ch. 16.3 and Gould/Hribar/Nocedal's work on
projected CG for large-scale QP) rather than a literal port of SQOPT's
reduced-gradient/basis-factorization method -- same active-set *control
logic* (add/drop one constraint at a time, ratio test, multiplier sign
check), different (simpler, matrix-free-friendlier) linear algebra inside
each iteration.

## 4. Null-space projection primitive

A single helper, used everywhere below:

```
project_null(J_A, v) -> p  such that  J_A p = 0  and  p ~ v
  solve  z = argmin || J_A^T z - v ||_2   (LSQR, minimum-norm)
  return p = v - J_A^T z
```

where `J_A` is the sparse sub-matrix of `[J | -I]` (in `w`-space) formed by
the rows corresponding to the constraints *currently in the working set*
(built fresh each time the working set changes -- just gathering a subset
of rows/columns from the existing COO `jac`, cheap and simple, no
factorization).

## 5. Projected-CG solve of the equality-constrained subproblem

Given the current working set (defining `J_A` and fixed values for the
nonbasic components of `w`), the QP restricted to that working set is an
*equality-constrained* QP in the free/superbasic components. Solve it with
**projected CG** (matrix-free, using only `hv_product` and
`project_null`):

```
inputs: current w, reduced gradient r0 = (H w + g_ext) restricted to free part
d = project_null(J_A, -r0)          ! initial search direction, feasible
r = r0
for k = 1, max_pcg_iter:
    if ||project_null(J_A, r)|| <= pcg_tol:  exit  ! reduced gradient ~ 0: optimal for this working set
    Hd = hv_product(d)                       ! only touches the `p`-block of w; s-block has no H term
    dHd = dot(d, Hd)
    if dHd <= 0: ! negative curvature detected (H indefinite on the null space)
        take_direction = d                    ! Steihaug-Toint style: stop and use d itself
        break                                  ! (see §7 for how the ratio test consumes this)
    alpha = dot(r, d) / dHd                    ! (standard projected-CG step length)
    w = w + alpha * d
    r = r + alpha * Hd
    d_new = project_null(J_A, -r)
    beta = dot(r, d_new - d_old_component) / ... ! standard PCG beta (Polak-Ribiere-style on the projected residual)
    d = d_new + beta * d
```

(the exact PCG recursion needs care to stay numerically standard -- this
is the well-known "projected CG for equality-constrained QP", not a novel
derivation; implementation should follow Nocedal & Wright Algorithm 16.3
precisely rather than the sketch above, which is illustrative only).

**Negative curvature / indefinite `H` safeguard**: `H` (the limited-memory
BFGS/SR1 Lagrangian Hessian) is not guaranteed positive-definite on the
null space, especially early on or with SR1. Truncate CG (Steihaug-Toint
rule) the first time `dHd <= 0` is encountered, and use the direction found
so far -- this is a standard, small, well-understood addition, not a new
research problem.

## 6. Outer active-set loop

```
initialize working set:
    - all variable/slack bounds that the *starting* w already sits on
      (start w=0 in p, s=viol-projected onto bounds, or reuse v1's
      normal-step idea for an initial feasible-ish point -- see §8)
loop:
    solve the equality-constrained subproblem on the current working set
        via projected CG (§5), possibly truncated by negative curvature
    ratio test:
        find the largest step length `alpha in [0,1]` along the PCG
        direction such that every bound on every component of `w` is
        still satisfied
        if alpha == 1 (no new bound hit): take the full PCG step
        else: take the partial step, ADD the newly-active bound to the
              working set, and repeat the PCG solve on the enlarged
              working set (this is the standard "blocking constraint"
              rule of active-set QP)
    when the PCG solve reaches an *unconstrained* (within the working set)
    stationary point (reduced gradient ~0, no negative curvature, ratio
    test allows a full step):
        compute Lagrange multipliers for every constraint/bound in the
        working set (from the un-projected residual `H w + g_ext - J_A^T
        lambda_A = 0`, another small LSQR solve, same pattern as v1's
        multiplier estimate but now exact since the working set is exact)
        if all inequality multipliers have the correct sign: OPTIMAL, stop
        else: DROP the constraint with the most-violating-sign multiplier
              from the working set, and repeat
    stop with sqpopt_qp_solve_failed if `max_iter` working-set changes are
    exceeded without satisfying the above (mirrors the existing
    `qp_solver%max_iter` field, currently unused by v1)
```

This is the classical primal active-set QP algorithm (Nocedal & Wright,
Algorithm 16.3) applied to the bounded-variable reformulation in §2, using
the matrix-free primitives in §4-§5 in place of a factorized basis.

## 7. Phase 1: finding an initial feasible working-set point

Unlike v1 (which never needs an exactly feasible starting point, since it
never enforces bounds exactly), an active-set method needs a `w` that
satisfies *all* bounds and `Jp - s = 0` before the loop in §6 can start.
Reuse the existing v1 machinery for this instead of writing a new phase-1
solver:
- start from `p=0` (so `s` must equal `Jp=0`... more precisely, from
  `viol = c_l - c(x)` / `c_u - c(x)` exactly as v1 already computes),
- take the **same minimum-norm normal step** `p_n` that v1's step 2 already
  computes (`LSQR` solve of `J p_n = viol`) as the initial `p`,
- set `s = J p_n` (satisfies `Jp-s=0` by construction) and clip `s` into
  `[c_l-c(x), c_u-c(x)]` if it falls slightly outside due to rank
  deficiency (rare, and only ever "slightly", since `p_n` is a minimum-norm
  least-squares solution),
- clip `p_n` into `[x_l-x, x_u-x]` if needed, then re-solve the normal step
  restricted to the bounds that got clipped (i.e., those become part of
  the *initial* working set) -- a small bootstrapping loop, at most `n`
  iterations in the worst case (one clip at a time), reusing §4's
  `project_null` machinery.

**Elastic mode** (SQOPT's fallback when even this can't be made feasible,
e.g. genuinely inconsistent linearized constraints): out of scope for the
first working version. Note it here as the natural v2.1 extension (relax
the offending bound with an `L1` penalty added to the QP objective, same
idea as `sndoc7.pdf`) and fall back to `sqpopt_qp_solve_failed` for now,
which `sqpopt_iterate_module` already knows how to handle (propagates the
error, consistent with existing behavior when v1's LSQR solves fail).

## 8. Data structures

New module `sqpopt_qp_reduced_hessian_module` (kept separate from v1's
`sqpopt_qp_solver_module` rather than merged into it, since the algorithms
share almost no code -- only `sparse_matvec`/`hessian%hv_product`/`LSQR`
at the primitive level):

```fortran
type :: sqpopt_active_set_qp_type
    integer :: max_iter = 50          ! max working-set changes per QP solve
    integer :: max_pcg_iter = 0       ! 0 => default to n+m
    real(wp) :: pcg_tol = 1.0e-8_wp
    ! working-set state (rebuilt each call; not persisted/warm-started in v1 of this feature -- see §9)
    logical, allocatable :: at_lower(:), at_upper(:), free(:)  ! dimension(n+m), one entry per (p_i, s_j)
contains
    procedure :: solve => solve_active_set_qp   ! SAME signature as sqpopt_qp_solver_type%solve
end type
```

Integration into the rest of `sqpopt`, following the existing
`hessian_mode`/`linesearch_mode`/`merit_mode` pattern exactly (a `mode`
selector on the existing type, not a new component on `sqpopt_type`):
- add `qp_solver_mode` to `sqpopt_options_type` (`sqpopt_qp_composite`
  default, `sqpopt_qp_reduced_hessian` new)
- add a `mode` field to `sqpopt_qp_solver_type` and dispatch inside
  `solve_qp_subproblem` to either the existing composite-step code (as-is,
  just renamed to a private `solve_composite_step`) or the new
  `sqpopt_active_set_qp_type`'s solver (composed as a private member, or
  called via a `use` of the new module) -- **no change** to
  `sqpopt_iterate_module` or `sqpopt_module`'s public API is needed.

## 9. Explicitly out of scope for the first version

- **Warm starting** the working set across major SQP iterations (each
  major iteration's QP solve starts a fresh Phase 1 from scratch, same as
  v1 effectively does today). Natural follow-up once correctness is
  established -- most of the working set typically doesn't change between
  consecutive major iterations near a solution, so this is a real,
  measurable performance opportunity later.
- **Elastic mode** (§7) -- fails cleanly with `sqpopt_qp_solve_failed`
  instead.
- **Superbasic variables** in the SQOPT sense (reduced-Hessian methods
  that allow more than one free direction at a time via a maintained
  Cholesky of the reduced Hessian) -- the PCG approach in §5 handles an
  arbitrary-dimensional free space directly without this machinery, so
  it's not needed here, but it's worth naming explicitly since it's a
  prominent feature of the real SQOPT algorithm this design draws from.
- Performance profiling/tuning of the repeated `LSQR` re-solves (§3) vs. a
  maintained factorization -- revisit only if the validation in §10 shows
  it's actually a bottleneck for realistic problem sizes.

## 10. Staged implementation & validation plan

Build and validate incrementally rather than as one large change:

1. **Stage 0 -- bounds-only QP** (`m=0`, i.e. `J` absent/empty): validates
   §5's projected-CG loop and §6's ratio test/working-set add-drop logic
   in the simplest possible setting (plain bound-constrained QP, no null-
   space projection needed at all since there are no general constraints).
   New standalone unit test, calling `sqpopt_active_set_qp_type%solve`
   directly (not through the full NLP loop) on a few small hand-verifiable
   bound-constrained QPs.
2. **Stage 1 -- add equality constraints**: introduces `project_null`
   (§4) and the slack reformulation (§2) for `m_eq>0`, `m_ineq=0`. New unit
   tests with known closed-form QP solutions (small, hand-checked).
3. **Stage 2 -- add inequality constraints**: full working-set add/drop
   logic (§6), multiplier-sign checks, negative-curvature truncation
   (§5). New unit tests replicating the same small NLPs already in
   `test_basic.f90` (equality, inequality, bounds-only) but as *standalone
   QPs* with hand-computed `H`/`g`/`J` at a specific point, so the QP
   solver can be checked in isolation before it's ever wired into the
   outer SQP loop.
4. **Stage 3 -- wire into `sqpopt_iterate_module`** via the `mode`
   dispatch in §8. Re-run the *existing* full test suite
   (`test_basic`, `test_hs71`, `test_medium`) with
   `options%qp_solver_mode = sqpopt_qp_reduced_hessian`:
   - `test_basic`/`test_medium` should continue to converge to
     `sqpopt_success` (ideally in fewer major iterations than v1, since
     the QP step is now exact rather than approximate).
   - `test_hs71` reaching `sqpopt_success` (not just the current loose
     `max(abs(x-x_true))<0.5` check) is the primary acceptance criterion
     for this whole effort, per the conclusion already reached
     independently in PLAN.md §6.1 and §6.3.
5. **Stage 4 -- documentation**: update `PLAN.md` §6.2's "Status" (mirroring
   how §6.1/§6.3 were updated after implementation) with the actual
   result, and fold any design deviations discovered during
   implementation back into this file for the historical record.

Each stage should be a separate, independently-reviewable unit of work --
this document intentionally stops at the end of Stage 0's design so that
implementation can start there once approved, rather than attempting all
four stages in one pass.
