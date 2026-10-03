# Plan: one function vector, one derivative matrix, constant elements

Status: plan only (2026-10-02). Nothing here is implemented.

## Goal

1. **One function vector.** The user evaluates `F(x) = [f(x); c(x)]` in one
   routine, and its derivatives `G(x) = F'(x)` in another, as one sparse
   matrix whose first row is `∇f`. This is SNOPT's `snOptA` layout with its
   `ObjRow = 1`.
2. **Constant elements.** Optionally, the user lists elements of `G` that
   are constant, as values with their row and column indices. The user's
   derivative routine then doesn't compute them.
3. **Use the structure.** The solver derives what the structure implies
   (linear constraints, a linear objective, linear variables, the pattern
   of the Hessian of the Lagrangian) and uses it where that is measured to
   help.

Parts 1 and 2 change no results: the HS suite must give the same counts
after them. Part 3 is a set of separate, measured changes.

## What exists now

- The user gives `fc(x, f, c, status, data)` and
  `gjac(x, g, jac_val, accuracy, status, data)`, optionally
  `hess(x, lambda, hess_val, status, data)`, through `set_functions`, and
  the constraint Jacobian's pattern through `set_jacobian_sparsity`. The
  gradient is dense: `g` has `n` values.
- Everything else in the solver gets `f`, `c`, `g`, and the Jacobian's
  values through the evaluation layer in `sqpopt_problem_module`
  (`problem%f`, `%c`, `%g`, `%jac`, with caches and scaling). Nothing
  outside that module calls the user's functions. This is what keeps parts
  1 and 2 contained.
- The interface is implemented 72 times (`set_functions` calls in 38 files:
  tests, examples, the HS harnesses, `sqpopt_nlls_module`, the Python shim),
  with about 100 derivative routines.
- Linearity is known nowhere in the Fortran solver. Python's
  `LinearConstraint` is evaluated as `A x` with a fixed Jacobian, but the
  solver isn't told. Roadmap F11 (flag linear rows) is not started.

## Decisions to make first

Each has a recommendation.

**D1. Which row is the objective.** *Recommended: row 1.* `F(1) = f`,
`F(1+i) = c_i`, and the pattern's rows run from `1` to `m+1`. That is what
"the first row" means in Fortran, and it is SNOPT's convention. The
alternative is rows `0:m`, which keeps a constraint's row equal to its
index. It is tempting, but a user who declares the dummy argument as
`dimension(:)` instead of `dimension(0:)` would silently get the objective
in `F(1)` and every constraint shifted by one. With row 1, a declaration
can't go wrong. The cost is that constraint `i` is row `i+1` in the
pattern, while `c_lb(i)`, `c_ub(i)`, the multipliers `lambda(i)`, and
`results%c(i)` stay indexed by constraint. The guide must say this
plainly, and `validate` should reject rows above `m+1`.

**D2. What a constant element means.** *Recommended: only its derivative
is constant, and `F` still includes its term.* The user's `F(x)` is the
whole function, as now, and the constant elements are simply left out of
what the derivative routine returns. SNOPT's `snOptA` does it the other
way: `F(x) = f(x) + A x`, where the user's routine computes only the
nonlinear part `f(x)` and the solver adds `A x`. That lets the solver
evaluate a fully linear row without calling the user. But it changes what
`F` means, and a user who leaves a linear term in `f(x)` as well gets it
counted twice, silently. The gain is small, since `F` is evaluated
anyway. The recommended meaning is also exactly "constant elements of the
Jacobian", as asked.

**D3. Whether to keep the old interface.** *Recommended: a clean break.*
The guide already calls the API a work in progress, and keeping
`fc`/`gjac` beside the new routines doubles the evaluation layer and every
guide snippet. During the migration (stage 2), the old path stays only so
that the tests pass at every step, and it is removed at the end of that
stage, before anything is released.

**D4. Names.** *Recommended:*
- the callbacks `sqpopt_fun_func(x, F, status, data)` and
  `sqpopt_deriv_func(x, G_val, accuracy, status, data)`;
- `problem%set_derivative_sparsity(nnz, irow, icol)`, which replaces
  `set_jacobian_sparsity`;
- `problem%set_constant_derivatives(nnz, val, irow, icol)`;
- `problem%set_functions(fun, deriv, hess, data)`.

`set_jacobian_sparsity` should be removed rather than given the new
meaning: old code that passes a constraints-only pattern would otherwise
compile and silently get a zero gradient.

**D5. Overlap.** *Recommended: as in SNOPT, an element is either in the
varying pattern or a constant, never both.* `validate` reports an overlap,
and duplicates within either list, as `sqpopt_invalid_input`.

**D6. The Hessian callback.** *Recommended: unchanged.* It stays the
Hessian of the Lagrangian with `lambda(m)` per constraint. Changing it to
SNOPT's `(σ, λ)` form gains nothing here.

## Stage 1: the new interface, behind the existing evaluation layer

Only `sqpopt_problem_module` changes. The solver still sees a dense `g`
and the constraint Jacobian in its own COO pattern.

- The new callbacks, setters, and type components: the varying pattern
  (`G_irow`, `G_icol`, rows `1..m+1`), and the constant pattern and values
  (empty in this stage, except for validation).
- An internal map, built once in `validate`, from the user's pattern to
  what the solver uses. The entries of row 1 go to positions of the dense
  `g`, those of rows `2..m+1` to the internal Jacobian `(i-1, j)`, whose
  pattern is the constraint rows of the varying set (plus, in stage 3, the
  constants). The internal Jacobian keeps the user's order of entries, so
  the results stay the same.
- `raw_fc` becomes `raw_fun`: the user writes `F` into a work vector of
  length `m+1` kept in the type (one vector, so nothing new on the stack),
  and that is split into the caches `cache_f` and `cache_c`.
- `raw_gjac` becomes `raw_deriv`: the user writes the varying values into
  a work vector, which is scattered into `cache_g` (zero outside row 1's
  pattern) and `cache_jac`.
- `validate`: rows in `1..m+1`, columns in `1..n`, no duplicates, and the
  checks of D1 and D5. New cases in `test_input_validation`.
- The old `fc`/`gjac` path stays during stage 2 (D3).

**Done when:** a new test (`test_structure.f90`) solves HS71 and a sparse
problem through the new interface with the same iterates as through the
old one, and `pixi run fpm test` passes.

## Stage 2: migrate every caller, then remove the old interface

This is mechanical but large. It is mostly the same edit 72 times, which a
script can draft: merge `f` and `c` into `F`, put `g` in front of
`jac_val`, and add row 1 to each pattern.

- **Tests** (31 files), **examples** (`hs71.f90`, `benchmark.f90`,
  `benchmark_large.f90`, `settings_study.f90`, `nlls_bard.f90`), and the
  **HS harnesses** (only their wrappers `fc_obj_cons`/`gjac_grad_jacv`
  change, not the problem modules).
- **`sqpopt_nlls_module`:** its own `fun`/`deriv` for the transformed
  problem. The user's interface for residuals is unchanged.
- **Python shim** (`python/sqpopt/fortran/sqpopt_python.f90`,
  `_sqpopt.pyf`, `sqpopt_python_f2py.f90`) and `_minimize.py`'s callbacks.
  The Python API (`minimize`, `least_squares`) is unchanged: its users
  never see `F` or `G`. Rebuild with `pixi run build-python`.
- **Docs:** the guide's Quick start, User functions, the worked example
  (kept identical to `example/hs71.f90`), Scaling (`fc_y`/`gjac_y`), and
  Least-squares sections, the README if it shows the interface, and
  CLAUDE.md. The guide gets a short migration note ("what changed and how
  to convert"), since this breaks user code.
- Then remove `sqpopt_fc_func`, `sqpopt_gjac_func`, `set_jacobian_sparsity`,
  and the old evaluation path.

**Done when:** all tests pass in both builds (`fpm test`, `test-mumps`),
the HS suite in release gives exactly the current baseline (281 / 24 / 0,
9,121 `fc`, and the other five configurations in CLAUDE.md), and
`test-python` passes. Any difference in the counts is a bug in the
migration.

## Stage 3: constant elements

- `set_constant_derivatives(nnz, val, irow, icol)`, rows `1..m+1`
  (row 1: linear terms of the objective). The values are fixed for the
  solve, given in the user's units, and scaled with their row like any
  other entry.
- The internal Jacobian's pattern becomes the varying entries plus the
  constants. The constants' values are written into `cache_jac` (and
  `cache_g`) once, when the caches are set up, and the user's routine
  fills only the varying positions. With no constants, nothing changes.
- **Checking them.** A wrongly declared constant corrupts the solve
  silently, so the diagnostics check them:
  - at level 2, the passive derivative check (which compares the change
    in `F` along each step with the derivatives at both ends) includes the
    constant entries, at no extra cost;
  - at level 3, a check at the starting point by differences of `F` in
    the constant entries' columns. This costs function calls, so it
    belongs at that level.
- **Python:** each `LinearConstraint`'s elements become constants
  automatically, since their values are known exactly. Its rows then need
  no Jacobian call in the shim.
- **`sqpopt_nlls_module`:** the `-1` entries of the auxiliary variables
  are constants.
- **Tests:** a problem with a linear objective term and linear rows, solved
  with and without the constants declared, giving the same results; a
  wrong constant caught by the level-2 diagnostics; and validation of
  overlaps.

**Done when:** as for stage 2 (no change to any result), plus the new
tests.

## Stage 4: derive the structure

Computed once in `validate`, from the two patterns and their constants:

- **linear rows:** every element of the row is constant (row 1: a linear
  objective);
- **linear variables:** every element of the column, row 1 included, is
  constant;
- **the structural Hessian pattern:** the union, over the rows that are
  not linear, of the cliques of each row's varying columns (a constant
  element `(i, j)` means row `i` has no second derivatives involving `x_j`).

These go into the printed header ("3 of 10 constraints linear, 12 of 40
variables linear") and into `results`, and the diagnostics report a
variable that appears nowhere (an all-zero column of `G`). Nothing in the
algorithm changes yet. There is a test of the derived quantities on small
patterns.

## Stage 5: use the structure

Each item is a separate change: measured, defaulted or made opt-in on the
measurements, and documented with its numbers, following CLAUDE.md's rules
for changes that can affect convergence. They are listed in a sensible
order of work.

**5a. No second-order correction on linear rows** (small). Their
linearization is exact, so `sqpopt_soc_module` drops them from the rows it
corrects. Saves least-squares work. Results should be unchanged up to
roundoff.

**5b. A starting point that satisfies the linear constraints** (medium).
SNOPT first finds a point that satisfies the linear constraints and the
bounds, and then every iterate keeps them satisfied, because their
linearization is exact.
- Here: before the first iteration, move `x0` to the nearest such point
  (a QP with `H = I` in the existing solvers).
- If none exists, the problem is infeasible, which can be reported at once
  (`sqpopt_infeasible`).
- Restoration steps act on the whole violation, so the plan must check
  that they don't break linear feasibility (or restrict them to steps
  that keep it).
- Measure on CUTEst: PyCUTEst flags linear constraints
  (`is_linear_cons`), and the benchmark script can pass them as constants.

**5c. No elastic relaxation of linear rows** (small to medium). Once the
linear constraints hold (5b), an inconsistent QP comes from the nonlinear
rows only. Relaxing only those keeps the QP's step consistent with the
linear constraints. Depends on 5b.

**5d. Quasi-Newton curvature only on the nonlinear variables** (large).
- **Why:** today L-BFGS starts from `γ I` on every variable, so a linear
  variable gets curvature it doesn't have, which restrains its steps until
  the updates learn better.
- **What changes:** an initial matrix `diag(γ on nonlinear, δ on linear)`,
  which needs the compact L-BFGS form to take a diagonal `B0` (`γ` is used
  in 15 places in `sqpopt_hessian_module`).
- **Why it is large:** with `δ = 0` the reduced Hessian can be singular,
  which the QP solvers don't handle (roadmap B12: the modified Cholesky
  then gives huge steps). So either a small `δ`, or making the QPs handle
  semidefinite reduced Hessians, which SNOPT's QP does.
- **Risk:** the most likely of these items to change results. Measure on
  HS and CUTEst.

**5e. A finite-difference Hessian of the Lagrangian with colouring**
(medium).
- **What:** with the structural pattern of stage 4, `∇²L` can be built from
  differences of the derivative routine, with the columns grouped by a
  greedy colouring of the pattern (Curtis, Powell, and Reid). It costs one
  call per colour, roughly the bandwidth plus one, instead of the `2n`
  calls of the HS harness's `hess_fd`.
- **Why:** the exact-Hessian, inertia-control, and direct-QP options would
  then work for users who have only first derivatives.
- **How:** a new option (e.g. `hessian_mode = sqpopt_hessian_differences`)
  with its schema entry, validation, and guide row. The HS harness's
  `hess_fd` can be replaced by it.
- **Expectation:** on HS the exact Hessian is behind L-BFGS (270 against
  281 solved), so the payoff is on large problems where L-BFGS struggles.
  Measure with `benchmark_large` and the hanging chain.

**Not in this plan** (roadmap): partitioned quasi-Newton updates per
element (LANCELOT), which the row structure makes possible but which is
research-grade; and SNOPT's additive `f(x) + A x` form (D2).

## Testing and measurement

- **Stages 1–3:** no result may change. The HS suite (release, all six
  configurations in CLAUDE.md), `fpm test`, `test-mumps`, `test-python`,
  `fortitude check`, and `tools/stack_check.sh` must all pass at the end of
  each stage.
- **Stage 5:** each item is compared on:
  - the HS suite. The harness can declare linear rows by checking the
    Jacobian at a few random points, which is good enough for a test
    harness and documented as such;
  - CUTEst (`pixi run cutest run`, with `is_linear_cons` passed as
    constants);
  - `benchmark_large`, for 5d and 5e.
- **Python:** the schema test covers any new options (5d's `δ`, 5e's mode).

## Effort and order

| stage | effort | changes results | depends on |
|---|---|---|---|
| 1. interface behind the evaluation layer | 1 session | no | D1–D6 |
| 2. migrate the callers, remove the old interface | 1–2 sessions (mostly mechanical) | no | 1 |
| 3. constant elements | 1 session | no | 2 |
| 4. derived structure | under 1 session | no | 3 |
| 5a. no SOC on linear rows | small | barely | 4 |
| 5b. linear-feasible start | 1 session | yes | 4 |
| 5c. no elastic relaxation of linear rows | small | yes | 5b |
| 5d. curvature only on nonlinear variables | 2+ sessions | yes | 4 (and B12) |
| 5e. coloured finite-difference Hessian | 1–2 sessions | new option | 4 |

Stages 1–4 are worth doing together: they are the interface change and its
groundwork. Stage 5 should be decided item by item on its measurements.
5b and 5e are the most likely to pay off. 5d is the most likely to cost
robustness.

## Risks

- **Breaking change.** Every user's code changes once. The guide gets a
  migration note, and the change is announced in the summary of the
  release.
- **The row offset (D1).** Constraint `i` is row `i+1` in the pattern but
  `i` everywhere else. Validation catches rows out of range, but not a
  pattern shifted by one inside the range. A guide example, and a test
  that gives a deliberately shifted pattern and checks the derivative
  diagnostics flag it, mitigate this.
- **Wrong constants.** These are caught only by the diagnostics (level 2
  passively, level 3 actively), which are opt-in.
- **Ground rules.** Every new work vector of size `n`, `m`, or `nnz` is
  allocatable and kept in the type (`tools/stack_check.sh`), and no state
  carries over between solves (`test_resolve`).
