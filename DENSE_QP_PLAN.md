# Design plan: dense QP solver option (`sqpopt_qp_dense`)

**Status: design only -- no code written yet.** Prompted by comparing
`test_hs71` against the newly-added `test/slsqp_test_71.f90`: `slsqp`
converges to machine precision in 6 iterations on the same problem where
`sqpopt`'s v1 composite-step QP only manages a loose limit cycle. `slsqp`
gets that tight convergence using an entirely **dense** BFGS Hessian +
dense active-set-style QP subproblem (see the investigation in §1 below).
This plan adds an *opt-in* `sqpopt_qp_dense` mode that forms dense
matrices from `sqpopt`'s existing sparse/matrix-free representations each
iteration and solves the QP subproblem with dense linear algebra -- for
users who know their problem is small enough that "dense" is a fine
trade for the very tight, very reliable convergence it buys.

This is the second of two backlog QP options (see
[REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md) for the sparse
one) and, as §9 argues, is very likely the one to **build and validate
first** -- it's simpler and lower-risk than the sparse plan, since dense
linear algebra (one-shot QR/Cholesky factorizations) sidesteps the
iterative-tolerance subtleties of the sparse plan's projected-CG backend.

## 1. What's actually reusable from `slsqp`/`bvls` -- an honest assessment

Before designing anything new, it's worth being precise about what the
`slsqp` dependency actually offers, since it's tempting to assume its
proven, tested dense QP code can just be called directly.

- **`slsqp_core.f90`'s `lsq`/`lsei`/`lsi`/`ldp`/`nnls`/`hfti`/`h12`/`g1`**
  -- the actual dense-QP-subproblem-forming machinery (Cholesky-factor the
  BFGS matrix, transform to a least-squares problem, eliminate equality
  constraints, reduce general inequalities to a "least distance
  programming" problem, solve via NNLS) -- are **all private** module
  procedures inside `slsqp_core`, with only the top-level `slsqp`
  subroutine itself `public`. This is the exact same situation already
  documented for `lbfgsb` in `PLAN.md` §2 (`bmv` private, only `setulb`
  public): there is no way to reuse just the QP-forming piece without
  copying private source code.
- **`bvls_module`'s `bvls`/`bvls_wrapper`** (Lawson-Hanson Bounded-Variable
  Least Squares) *is* public, small, and self-contained (only depends on
  `slsqp_kinds`/`slsqp_support` for `wp`/a few constants) -- a genuinely
  reusable building block, **but it only solves a bound-constrained
  least-squares problem** (`min ||Ax-b||_2` s.t. simple box bounds
  `lo<=x<=hi`). It was tempting to use it as the dense QP engine after
  eliminating equality constraints via a change of basis, but **that
  doesn't work cleanly**: once you eliminate the equality constraints
  `Jp-s=0` via a change of variables `w = w_0 + Z u` (§3), the *original*
  simple bounds on `w` become **general linear inequalities in `u`**
  (`w_lb <= w_0+Zu <= w_ub`), not simple box bounds on `u` -- exactly the
  reduction that `slsqp`'s private `lsi`/`ldp` chain exists to handle, and
  precisely why `bvls` alone can't substitute for it.

**Conclusion**: don't try to reuse `slsqp`/`bvls` internals for this.
`slsqp` stays a dev-dependency (useful for reference/comparison tests like
`test/slsqp_test_71.f90`, not linked into the library). This plan instead
writes a small amount of new, self-contained dense linear algebra (§4) and
reuses the **same active-set control logic** already designed for the
sparse QP option (§6), just with a dense backend.

## 2. Goal and scope

Same QP subproblem as always (see `sqpopt_qp_solver_module`'s docs and
[REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md) §1):

$$ \min_{p} \; \tfrac12 p^T H p + g^T p \quad \text{s.t.} \quad
c_l - c(x) \le J p \le c_u - c(x), \quad x_l - x \le p \le x_u - x $$

`sqpopt_qp_dense` is a **third** `sqpopt_qp_solver_type` mode (alongside
`sqpopt_qp_composite`, the current v1 default, and `sqpopt_qp_reduced_hessian`,
the planned sparse option), selected explicitly by the user, that:
- forms the dense Jacobian `J` (`m x n`) and dense Hessian `H` (`n x n`)
  from the existing sparse COO `jac` and matrix-free `hessian%hv_product`
  each time they're needed (§5),
- solves the QP subproblem with dense one-shot QR/Cholesky factorizations
  instead of the v1 heuristic or the sparse plan's iterative projected CG,
- is a drop-in alternative behind the same
  `solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)`
  interface, so no changes are needed to `sqpopt_iterate_module` or
  `sqpopt_module` beyond the same `mode`-dispatch pattern already planned
  for the sparse option.

**Explicitly opt-in, not a new default**: forming `H` and `J` densely is
`O(n^2)`/`O(m n)` memory and the QR/Cholesky factorizations are `O(n^3)`-ish
per working-set change -- fine for the small-to-moderate `n` this targets
(roughly the same range where `slsqp` itself is the right tool), wrong for
the large sparse problems `sqpopt`'s default composite-step/sparse modes
are designed for. `sqpopt` should never form a dense array *unless the
user explicitly asked for this mode*.

## 3. Reformulation (shared with the sparse plan)

Reuse [REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md) §2's
bounded-variable-plus-slack reformulation verbatim: `w=(p,s)`
(`dimension(n+m)`), one block of linear equality constraints `[J|-I]w=0`,
simple bounds on every component of `w`. This is what unifies variable
bounds and linearized general-constraint bounds into one bounded-variable
framework, and is *not* specific to sparse vs. dense -- only what happens
next (computing a null-space basis and solving the reduced problem)
differs.

## 4. New dense linear-algebra primitives needed

None of `sqpopt`'s current dependencies (`LSQR`, `LSMR`, `lusol`, `fmin`)
provide dense QR or Cholesky factorizations (they're all sparse/iterative
by design). A new small, self-contained module,
`sqpopt_dense_linalg_module`, is needed:

- **`dense_qr(A) -> Q, R`** (or just enough of a Householder QR to get an
  orthonormal basis for `null(A)`): given the (`m_A x (n+m)`) active
  constraint matrix `J_A` (a sub-selection of rows of `[J|-I]`), compute an
  orthonormal `Z` (`dimension(n+m, n+m-m_A)`) spanning its null space. This
  is the dense analogue of the sparse plan's repeated `LSQR`
  `project_null` calls (§4 there) -- classic, well-documented, textbook
  Householder QR (e.g. Golub & Van Loan, or the same technique underlying
  `slsqp_core`'s own private `h12`/`g1` Householder helpers, which we
  cannot reuse directly per §1 but can freely reimplement from the
  well-known algorithm).
- **`dense_cholesky_modified(A) -> L, ok`**: a standard **modified
  Cholesky** (Gill-Murray-Wright-style: add a small multiple of `I`, or
  perturb a pivot, whenever it would otherwise go non-positive) --
  factorizes the reduced Hessian `Z^T H Z` (`dimension(n+m-m_A,
  n+m-m_A)`), which is not guaranteed positive-definite (the limited-memory
  BFGS/SR1 `H` isn't guaranteed PD on an arbitrary subspace, same caveat as
  the sparse plan's negative-curvature handling in its projected-CG loop
  -- §5 there). Small, self-contained, no external dependency needed.

Both are small (a few hundred lines total), standard, easily unit-testable
in isolation (feed in a known matrix, check `Q^T Q = I`/`A=LL^T`) before
ever being wired into the QP solver -- much lower implementation risk than
the active-set control logic itself.

## 5. Densifying `H` and `J` from the existing representations

- **`J` (dense, `m x n`)**: scatter the existing sparse COO `jac%val` into
  a dense array -- trivial, one loop over `jac%nnz`.
- **`H` (dense, `n x n`)**: apply `hessian%hv_product` (already exists,
  matrix-free, unchanged) to each of the `n` unit basis vectors `e_i` and
  collect the results as columns. This **reuses the existing
  `sqpopt_hessian_module` unmodified** -- `hessian_mode` (BFGS/SR1) keeps
  working exactly as it does today, and no second, separate dense-BFGS
  update path needs to be written and validated. Cost: `n` extra
  `hv_product` calls whenever `H` needs to be (re-)densified -- acceptable
  given this mode already targets small `n` (the same range where `n`
  calls to a cheap matrix-free operator is a non-issue).

Both are only ever formed for this QP solve when `sqpopt_qp_dense` is
selected; the rest of `sqpopt` (Hessian approximation, line search,
merit function, convergence test, reporting) is completely unaware of
the difference and continues to operate on the sparse/matrix-free
representations as always.

## 6. Active-set algorithm: same control logic as the sparse plan

This is the key simplification: **the outer active-set loop doesn't need
to be designed twice.** [REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md)
§6 already specifies it in a backend-agnostic way (initialize a working
set, solve the equality-constrained subproblem on it, ratio-test to find a
blocking bound, add it and repeat, or check multiplier signs and drop a
constraint) -- the *only* thing that differs here is what "solve the
equality-constrained subproblem on the current working set" means:

| step | sparse plan (`sqpopt_qp_reduced_hessian`) | dense plan (`sqpopt_qp_dense`) |
|---|---|---|
| null-space basis of `J_A` | implicit; get a projection via one `LSQR` solve each time (§4 there) | explicit orthonormal `Z`; one dense Householder QR per working-set change (§4 here) |
| solve reduced problem | projected CG (matrix-free `hv_product`), truncated on negative curvature (§5 there) | direct: form `Z^T H Z` (dense, small), modified Cholesky, one triangular solve (§4 here) |
| cost per working-set change | several `LSQR` solves + a CG loop | one QR + one Cholesky (both `O((n+m-m_A)^3)`-ish, but `n+m` here is small by construction) |
| robustness | iterative; needs a tolerance (`pcg_tol`) and truncation rule | direct; no iterative tolerance, generally the more numerically robust choice for small dense systems (part of why `slsqp` itself converges so tightly) |

Phase 1 (finding an initial feasible working-set point) and the
elastic-mode/out-of-scope notes in
[REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md) §7/§9 carry over
unchanged -- reuse the same minimum-norm-normal-step bootstrapping idea
(here, computed with a dense least-squares solve instead of `LSQR`, e.g.
via the same `dense_qr` primitive).

## 7. Data structures & integration

New module `sqpopt_qp_dense_module`:

```fortran
type :: sqpopt_dense_qp_type
    integer :: max_iter = 50   ! max working-set changes per QP solve (mirrors the sparse plan)
contains
    procedure :: solve => solve_dense_qp   ! SAME signature as sqpopt_qp_solver_type%solve
end type
```

plus the new `sqpopt_dense_linalg_module` (§4), used only by this module
(and its own unit tests, §8) -- no other part of `sqpopt` needs dense
linear algebra.

Dispatch: the **same** `mode` field on `sqpopt_qp_solver_type` already
planned for the sparse option
([REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md) §8), now a
3-way choice --

```fortran
integer, parameter :: sqpopt_qp_composite       = 1  ! v1 heuristic (current default)
integer, parameter :: sqpopt_qp_reduced_hessian = 2  ! sparse active-set (planned)
integer, parameter :: sqpopt_qp_dense           = 3  ! dense active-set (this plan)
```

-- dispatched inside `solve_qp_subproblem`, so `sqpopt_iterate_module` and
`sqpopt_module`'s public API need no changes at all beyond what the sparse
plan already requires.

No `fpm.toml` changes are needed: per §1, this plan does not deepen the
dependency on `slsqp`/`bvls` beyond its current dev-dependency role.

## 8. Staged implementation & validation plan

1. **Stage 0 -- dense linear algebra unit tests**: `dense_qr` and
   `dense_cholesky_modified` validated in isolation against known small
   matrices (`Q^T Q = I`, `A = L L^T`, a deliberately indefinite matrix to
   exercise the modified-Cholesky safeguard) before anything else is
   built on top of them.
2. **Stage 1 -- bounds-only dense QP**: same idea as
   [REDUCED_HESSIAN_QP_PLAN.md](REDUCED_HESSIAN_QP_PLAN.md) Stage 0 (no
   general constraints, `J` absent), validating the active-set add/drop
   loop with the dense Cholesky backend in the simplest setting.
3. **Stage 2 -- full QP with equality + inequality constraints**: new
   standalone unit tests (small, hand-verifiable QPs, called directly via
   `sqpopt_dense_qp_type%solve`, not through the full NLP loop yet).
4. **Stage 3 -- wire into `sqpopt_iterate_module`** via the `mode`
   dispatch (§7). Re-run the existing full test suite with
   `options%qp_solver_mode = sqpopt_qp_dense`:
   `test_basic`/`test_medium` should keep converging to `sqpopt_success`,
   and **`test_hs71` reaching `sqpopt_success`** is again the primary
   acceptance criterion -- directly testable against the now-concrete
   comparison point of `slsqp_test_71` converging in 6 iterations.
5. **Stage 4 -- documentation**: update `PLAN.md` with the result, same as
   was done for §6.1/§6.3 after their implementations.

## 9. Recommendation: build this one before the sparse plan

Both QP plans target the exact same problem (`test_hs71`'s limit cycle)
via the exact same active-set control logic (§6); they differ only in the
linear-algebra backend. Given that:

- the dense backend (one-shot QR + modified Cholesky) has no iterative
  tolerances or truncation rules to get right, unlike the sparse plan's
  projected-CG loop (fewer numerical judgment calls, easier to get
  correct the first time),
- `test_hs71` (`n=4`, `m=2`) is comfortably within the size range where
  "dense" is a completely reasonable default anyway, and
- `slsqp`'s own dense approach is *proof* (not just theory) that a dense
  active-set-style QP fixes exactly this failure mode on exactly this
  problem,

**this plan is the faster, lower-risk way to validate the "a real QP
solve fixes `test_hs71`" hypothesis** that's been the recurring conclusion
throughout `PLAN.md` §6.1/§6.3. Once validated here, the same active-set
control logic (§6) carries over directly to the sparse plan, which mainly
needs its *backend* (projected CG replacing dense QR/Cholesky) validated
separately -- de-risking that larger effort too.
