# sqpopt: basis repair and basis-selection settings for the sparse QP's basis method

Changes to the reduced-Hessian (SQOPT-style) sparse QP solver. Its basis method (`null_space = sqpopt_null_space_lu`) picks a basis with a rank-revealing LUSOL factorization. The changes are in `src/sqpopt_qp_reduced_hessian_module.f90` and `src/sqpopt_linalg_module.f90`, on top of `d8ff0e0`.

The new behavior is opt-in. With the defaults, every QP is solved as before. The one exception is the LUSOL workspace retry count (item 3), which only matters where the old code gave up.

## Motivation

These came from the eclipse optimization of a 20-year NRHO multiple-shooting trajectory in halo. The problem has 31,367 variables and 26,888 constraints (26,880 equalities), with 376,376 Jacobian nonzeros.

With the defaults, the basis method failed on every QP, and every QP fell back on the LSQR method, which takes hours at this size:

1. **The basis came out 2–5 columns short.** `choose_basis` reported `count(chosen) < m`, but the matrix isn't rank deficient. LUSOL found full rank, and a check outside sqpopt confirmed it. The missing columns are pivot columns that `lu1fac` flags with `w(j) <= 0`, because their diagonal of `U` is small relative to the rest of their column of `U` (the `rel_tol = 1e-8` test). `independent_columns` drops them, so their pivot rows have no basic column.

   MINOS and SNOPT handle this case with a *basis repair*: each such row's slack replaces the rejected column. The basis is then square and nonsingular, and the active-set iterations proceed from it.
2. **With equilibrated (better-conditioned) problem data, `lu1fac` returned `inform = 7` three times in a row.** That code means insufficient storage. The workspace is doubled on each retry, and three retries weren't enough for this basis.

With the repair and more retries, the basis method solves each of these QPs in 2–64 iterations. On a 100-revolution version of the problem, sqpopt converges in 2.8 s with the default tolerances.

## Changes

### 1. New `sqpopt_reduced_hessian_qp_type` components

All of these are used only when `null_space = sqpopt_null_space_lu`.

| Component | Default | Meaning |
|---|---|---|
| `basis_repair` | `.false.` | If the basis picked has fewer than `m` columns, complete it with the slacks of the rows that need them, instead of giving up on the basis method. |
| `basis_rel_tol` | `1.0e-8` | The relative singularity tolerance of the factorization that picks the basis. This value was previously hard-coded. |
| `basis_eq_slack_weight` | `1.0e-6` | The column weight of the equality rows' slacks (the unknowns with equal bounds) when the basis is picked, relative to the free columns' `1`. This value was previously hard-coded. |

Both call sites of `choose_basis` (in the basic-solution and crash code, and in `basis_active_set`) pass these components.

### 2. `choose_basis`: optional `rel_tol`, `eq_slack_weight` and `repair` arguments

- The optional arguments default to the old hard-coded values (`1e-8`, `1e-6`, no repair).
- With `repair`, `independent_columns` returns the rows that need a slack (`repair_rows`). For each such row, the slack column with its single entry in that row is added to `chosen`. A basis that is still short is reported as before (`istat = -1`).
- The routine's docstring describes the repair.

### 3. `independent_columns`: optional `repair_rows`, and more workspace retries

- **`repair_rows`** has dimension `nrows`. It is `.true.` for every row that isn't the pivot row of an accepted column: LUSOL's k-th pivot is in row `p(k)` and column `q(k)`, for `k = 1..luparm(16)`, the rank, and a column is accepted if `w(q(k)) > 0`. So it covers both the pivot rows of the rejected columns and the rows without any pivot. It is `.false.` everywhere on a factorization error.
- **The workspace retry loop** goes from 3 attempts to 8, with the workspace doubled each time `lu1fac` returns `inform = 7`. This one change is not opt-in. Before, the method failed in this case; now it gets more memory. Successful factorizations are unaffected.

## Remaining for the merge request

Per the checklist in `CLAUDE.md`, these haven't been done yet:

- [ ] **Validation:** checks for the new components in the component's validation (`basis_rel_tol > 0`, `basis_eq_slack_weight > 0`), with a case in `test_input_validation`.
- [ ] **Python schema:** `_o(...)` entries in `python/sqpopt_options/schema.py` for the three components, then `pixi run build-python` and `pixi run test-python`.
- [ ] **User guide:** rows in the sparse-QP options table, and a paragraph on the basis repair (`web/index.html`).
- [ ] **Requirements:** rows in Table A.3, and a requirement for the repair, with its test (`plan/REQUIREMENTS.md`).
- [ ] **Test:** a `test/test_basis_repair.f90`. It would use a QP whose basis LUSOL rejects a small-pivot column of, and check that the QP fails over to LSQR without `basis_repair` and is solved by the basis method with it. It would also cover a case that needs more than 3 workspace retries, if a small one can be built.
- [ ] **Lint:** `pixi run fortitude check`.
- [ ] **HS suite:** `pixi run fpm test test_hs_suite --profile release`. The defaults are unchanged, so the baseline (280 / 25 / 0, 9,129 `fc`) should be too.
