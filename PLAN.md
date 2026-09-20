# sqpopt Implementation Plan

This document lays out the plan for turning the current `sqpopt` skeleton into a
working, sparse-native Sequential Quadratic Programming (SQP) solver for the
general problem:

$$ \min_x f(x) \quad \text{s.t.} \quad c_l \le c(x) \le c_u, \quad x_l \le x \le x_u $$

where equality constraints are `c_l = c_u` and inequality constraints are
`c_l \ne c_u` (using `\pm\infty` for one-sided bounds).

## 1. Current state

The skeleton (all bodies are `! TODO: implement`) already defines the module
boundaries and public interfaces:

| Module | Responsibility |
|---|---|
| `sqpopt_kinds` | selectable real working precision (32/64/128-bit) |
| `sqpopt_types_module` | status codes, `sqpopt_sparse_matrix` (1-based COO) |
| `sqpopt_problem_module` | NLP definition: sizes, bounds, Jacobian/Hessian sparsity patterns, user callbacks |
| `sqpopt_options_module` | tolerances, iteration limits, Hessian/linear-solver mode selectors |
| `sqpopt_hessian_module` | limited-memory (L-BFGS/L-SR1) Lagrangian Hessian, matrix-free |
| `sqpopt_qp_solver_module` | sparse QP subproblem solver |
| `sqpopt_linesearch_module` | merit function + line search |
| `sqpopt_convergence_module` | KKT/feasibility stopping tests |
| `sqpopt_linalg_module` | sparse matvec + dispatch to `lusol`/`LSQR`/`LSMR` |
| `sqpopt_iterate_module` | orchestrates one major SQP iteration |
| `sqpopt_module` | public `sqpopt_type` facade |

Key design decision already locked in: **no dense `n×n` or `m×n` arrays,
ever**. The Jacobian/Hessian are sparse COO triplets with a fixed sparsity
pattern; the Hessian approximation is a matrix-free limited-memory operator.

## 2. Dependency inventory & reuse strategy

We now have six reference dependencies fetched under `build/dependencies/`:
`slsqp`, `psqp`, `NumDiff`, `nlesolver-fortran`, `lbfgsb`, plus the sparse
linear algebra trio `LSQR`, `LSMR`, `lusol`.

### Directly reusable (linked into `sqpopt` itself, `[dependencies]`)

- **`lusol`** — sparse LU factorization (`lu1fac`/`lu6sol`, or the
  `lusol_ez_module.solve` convenience wrapper). Used as the default direct
  solver for the KKT/Newton systems in the QP subproblem.
- **`LSQR`** / **`LSMR`** — matrix-free iterative sparse least-squares
  solvers. Useful as an alternative to `lusol` for very large/ill-conditioned
  KKT systems, and directly usable for the initial least-squares estimate of
  Lagrange multipliers (a common SQP initialization trick, see below).
- **`lbfgsb`** (BSD, Jacob Williams' modernization of Nocedal/Morales'
  L-BFGS-B) — this is a genuinely good match for `sqpopt_hessian_module`.
  Its limited-memory BFGS machinery is *already matrix-free*: history is kept
  in `Ws`/`Wy` (`dimension(n,m)`, `m` a small constant) plus `Sy`/`Ss`/`Wt`
  (`dimension(m,m)`), and the compact-formula matrix-vector product is
  implemented in the standalone `bmv` subroutine (`lbfgsb_module`). Rather
  than write our own two-loop recursion from scratch, we will reuse `bmv`
  (and the S/Y/`Sy`/`Wt` bookkeeping it expects) directly inside
  `sqpopt_hessian_module` for the `sqpopt_hessian_bfgs` mode. Note `lbfgsb`
  itself only solves *box-constrained* problems (no general linear/nonlinear
  constraints) — we are not depending on its outer driver (`setulb`/`mainlb`,
  Cauchy point, subspace minimization), only its inner compact-BFGS
  matrix-vector product machinery.

`lusol`/`LSQR`/`LSMR` already take/accept 1-based COO `irow`/`icol`/`val`
triplets, which is exactly the convention used by `sqpopt_sparse_matrix` — no
conversion layer needed.

### Algorithmic inspiration only (NOT linked as dependencies — both are fully dense)

- **`slsqp`** (Kraft's SLSQP, BSD): dense damped-BFGS update (with the
  Powell/L1-test damping condition), and QP subproblem solved as a
  **bounded-variable least-squares (BVLS)** problem. We will port the *damped
  BFGS update formula* and the *L1 merit-function / multiplier-averaging
  logic* (`slsqpb`) algorithmically into `sqpopt_hessian_module` and
  `sqpopt_linesearch_module`, adapted to (a) limited-memory storage and (b)
  the general two-sided `c_l <= c(x) <= c_u` constraint form instead of
  SLSQP's `c_j = 0` / `c_j >= 0` split.
- **`psqp`** (Lukšan's PSQP, LGPL): dense variable-metric update (BFGS or
  Hoshino) plus a **dual range-space active-set QP** solver and an
  **extended line search without derivatives**. This is the best available
  reference for the *active-set logic* (constraint addition/deletion, working
  set management) that our sparse QP solver needs to replicate — we will
  reimplement the active-set bookkeeping using sparse factorizations
  (`lusol`) in place of PSQP's dense LDLᵀ range-space matrix.
- **`nlesolver-fortran`**: not an SQP algorithm (it's a Newton-type nonlinear
  *equation* solver), but its sparse-mode plumbing (`irow`/`icol` sparsity
  pattern set once, sparse Broyden updates, dispatch to `lusol`/`LSQR`/`LSMR`)
  is the direct template already followed by `sqpopt_problem_module` and
  `sqpopt_linalg_module`. Kept as a dev-dependency for reference/testing only.
- **`NumDiff`**: finite-difference Jacobian/gradient computation with
  automatic **sparsity pattern detection** (`dsm` column-partitioning
  algorithm). Kept as a dev-dependency; used in tests to (a) verify
  user-supplied analytic derivatives against finite differences, and (b) as
  an optional convenience path for users who don't want to hand-derive a
  sparsity pattern (future `sqpopt` feature, not required for v1).

### Nothing else has a ready-made sparse SQP implementation

No existing Fortran library in the dependency graph implements a sparse SQP
method directly — this is genuinely new work. The plan below is a from-scratch
implementation guided by the algorithms above.

## 3. Algorithm design decisions

1. **Hessian approximation**: limited-memory BFGS by default
   (`sqpopt_hessian_bfgs`), with limited-memory SR1 as an alternative
   (`sqpopt_hessian_sr1`) and user-supplied exact sparse Hessian
   (`sqpopt_hessian_exact`) for problems where it's cheap/available. BFGS
   update uses SLSQP-style damping to guarantee a positive-definite update
   even with non-convex problems.
2. **QP subproblem**: primal-dual active-set method. Working-set
   add/drop logic follows PSQP's dual range-space approach conceptually, but
   the linear algebra (solving the reduced KKT system for the current working
   set) is done sparsely via `lusol` (default) or iteratively via `LSQR`/
   `LSMR` (for large/ill-conditioned problems), selected by
   `options%linear_solver_mode`.
3. **Merit function / line search**: \( \ell_1 \) exact penalty merit
   function (as in SLSQP/PSQP) with dynamically updated penalty parameter,
   backtracking (Armijo) line search by default. A filter-based alternative
   is a possible v2 enhancement but not required for v1.
4. **Multiplier initialization**: use a sparse least-squares solve
   (`LSQR`/`LSMR` on \( J^T \lambda \approx -g \)) to get a good initial
   Lagrange multiplier estimate at iteration 0, rather than starting from
   zero — mirrors a common trick in dense SQP codes, adapted to sparse.
5. **Convergence**: KKT stationarity (\( \| \nabla_x \mathcal{L} \| \le
   ktol \)), constraint feasibility (\( \le ctol \)), and no-progress checks
   on \( f \) and \( x \) (`ftol`, `xtol`), matching `sqpopt_options_type`.

## 4. Implementation roadmap

Bottom-up order, so each layer can be unit-tested before the next depends on it.

1. **`sqpopt_types_module`** — finalize `sqpopt_sparse_matrix` helpers if
   needed (e.g. a constructor/validation routine). Mostly done.
2. **`sqpopt_linalg_module`** — implement `sparse_matvec`,
   `sparse_matvec_transpose` (trivial COO loops), and
   `solve_sparse_linear_system` dispatch to `lusol_ez_module`, `lsqr_module`,
   `LSMRmodule`. Unit test against small hand-built sparse systems with known
   solutions.
3. **`sqpopt_problem_module`** — implement `set_problem_size`,
   `set_jacobian_sparsity`, `set_hessian_sparsity`, `set_functions`
   (allocation + pointer assignment only, no numerics).
4. **`sqpopt_hessian_module`** — for the BFGS mode, reuse `lbfgsb`'s compact
   two-loop recursion (`bmv`) and its `Ws`/`Wy`/`Sy`/`Ss`/`Wt` storage
   convention instead of writing the recursion from scratch; wrap it behind
   `hv_product`/`inverse_vector_product` so callers never see the `lbfgsb`
   internals. Implement the circular buffer for `(s,y)` pairs and
   `update_bfgs` (with SLSQP-style damping) to keep `Sy`/`Ss`/`Wt` up to
   date. `update_sr1` (limited-memory SR1) has no `lbfgsb` equivalent and
   must be written from scratch. Unit test: verify quadratic convergence on
   a small unconstrained quadratic (no QP solver needed yet — just Hessian ≈
   true Hessian for a quadratic function).
5. **`sqpopt_qp_solver_module`** — implement the sparse active-set QP solver.
   This is the largest piece of new work. Suggested sub-steps:
   - a. Equality-constrained QP only (no inequalities/bounds) — solve the
     sparse KKT system directly via `lusol`. Validates the linear algebra path.
   - b. Add simple bound handling (active-set on `x_l <= x <= x_u`).
   - c. Add general inequality constraints (`c_l <= J p <= c_u`), full
     active-set add/drop iteration.
   - Test against small QPs with known solutions, and cross-check against
     dense reference solves (e.g. a temporary dense LAPACK path used only in
     tests, never in the library).
6. **`sqpopt_linesearch_module`** — implement the \( \ell_1 \) merit function
   and backtracking line search with SLSQP-style multiplier averaging for the
   penalty parameter.
7. **`sqpopt_convergence_module`** — implement KKT residual and feasibility
   checks.
8. **`sqpopt_iterate_module`** — wire steps 3–7 together into one major
   iteration.
9. **`sqpopt_module`** — implement `initialize`, `set_problem`, `set_options`,
   `solve` (drives `sqpopt_iterate` in a loop with `check_convergence`),
   `get_solution`, `destroy`.
10. **Test suite** (`test/`) — port the classic SLSQP/PSQP test problems
    (Hock-Schittkowski-style problems already present in
    `build/dependencies/slsqp/test/` and `build/dependencies/psqp/test/`) as
    `sqpopt` test cases with known optimal `x*`/`f*`, expressed via sparse
    Jacobian/Hessian patterns. Use `NumDiff` in tests to sanity-check any
    analytic derivatives written for the test problems.
11. **Documentation/examples** — FORD docs (a `ford.md` + `pixi add ford` is
    already set up), plus 2-3 worked examples in `example/` showing equality
    constraints, inequality constraints, and variable bounds together.

## 5. Open decisions to revisit with the user

- Default `lbfgs_memory` (currently 10) and default `linear_solver_mode`
  (currently `lusol`) — reasonable defaults, but worth revisiting once
  benchmarked.
- Whether an interior-point QP alternative to active-set is worth adding
  later for problems with many inequality constraints (active-set QP can be
  slow when the working set changes a lot).
- Whether to expose `NumDiff`'s automatic sparsity-pattern detection as a
  convenience path in `sqpopt_problem_module`, or keep sparsity patterns
  strictly user-supplied for v1 (simpler, faster to implement).
