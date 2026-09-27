# Comparison with Uno

[Uno](https://github.com/cvanaret/Uno) (Charlie Vanaret and Sven Leyffer,
MIT license) is a modular C++ solver for nonlinearly constrained
optimization. This note compares its design with sqpopt's and lists ideas
worth borrowing.

Sources: the Uno source code (commit `822a55e`, September 2026), its
documentation (<https://unosolver.readthedocs.io>), and the paper
*Implementing a unified solver for nonlinearly constrained optimization*
(Vanaret and Leyffer, [arXiv:2406.13454](https://arxiv.org/abs/2406.13454)).
Uno was not built or benchmarked for this comparison.

## Summary

Uno is a broader framework than sqpopt: one set of building blocks covers
SQP, interior-point and SLP (sequential linear programming) methods, using
external subproblem and linear solvers. sqpopt is a self-contained
quasi-Newton SQP in pure Fortran. The ideas most worth borrowing are Uno's
full feasibility restoration phase, its funnel method, and named presets.

## How Uno is organized

Uno splits every method into interchangeable "ingredients". Each is an
abstract C++ class, chosen by a string option at run time:

| Ingredient | Uno's choices | sqpopt's equivalent |
|---|---|---|
| What to do when the QP's linearized constraints are inconsistent | none, or a feasibility restoration phase | elastic slacks inside the QP, plus a Gauss-Newton restoration step and the second-order escape probe |
| How inequalities are handled | QP active set, or interior point (log barrier) | QP active set only |
| Hessian | exact, L-BFGS, L-SR1, identity, zero | L-BFGS, L-SR1 (exact is roadmap item F7) |
| Inertia correction | none, primal, primal-dual (IPOPT-style) | none (quasi-Newton, plus the QP's negative-curvature handling) |
| Step acceptance | merit function, Fletcher filter, Wächter filter, funnel | filter, Armijo, watchdog, exact; ℓ1 or augmented Lagrangian merit; two penalty rules; non-monotone retry |
| Line search or trust region | both | both |
| Subproblem and linear solvers | BQPD, HiGHS, MA27/57/86, MUMPS, SSIDS (external) | its own dense QP and sparse QP (LUSOL basis method, LSQR) |

**Presets** bundle whole configurations:

- `filtersqp`: trust region, Fletcher filter, exact Hessian;
- `ipopt`: line search, Wächter filter, barrier method, primal-dual inertia
  correction;
- `funnelsqp`: trust region, funnel method, exact Hessian;
- `filterslp`: `filtersqp` with a zero Hessian;
- `auto`: picks `ipopt` when n+m ≥ 2000 or the Jacobian and Hessian have
  at least 50,000 nonzeros, else `filtersqp`.

Uno has C, Python, Julia, AMPL and Fortran interfaces.

## Comparison

- **Similar philosophy.** Both treat line search vs trust region, filter
  vs merit, and the quasi-Newton flavor as selectable. sqpopt has more
  merit-function options (augmented Lagrangian, the GMSW and Byrd–Nocedal
  penalty updates, the non-monotone retry).
- **Where Uno is more factored.**
  - Its restoration is a full strategy, not a helper.
  - Its step acceptance works from predicted-reduction models shared by
    the filter, funnel and merit tests.
  - In sqpopt, the line-search type bundles the filter, merit, penalty and
    watchdog state together, and restoration is a single function.
- **Where sqpopt is different on purpose.**
  - It is pure Fortran with no licensed dependencies; Uno's strongest
    configurations need BQPD or HSL.
  - Its large-scale path is a matrix-free sparse QP with L-BFGS, where Uno
    uses sparse LDLᵀ factorizations.
  - That is sqpopt's niche: an SLSQP-like drop-in that scales.
- **Scope.** Uno's interior-point path and exact Hessians with inertia
  correction are what make it strong on large CUTEst problems. sqpopt has
  neither.

## Ideas worth borrowing

1. **A real feasibility restoration phase** (roadmap items F2 and F6).
   This is the biggest win.
   - When the QP is infeasible or the line search fails, Uno switches to
     minimizing the ℓ1 violation plus a proximal term around the
     pre-restoration point. That subproblem uses the same QP machinery
     and its own filter.
   - It blocks returning to the pre-restoration point (the optimality
     filter or funnel is told to avoid it).
   - It switches back only when the violation has dropped by a factor
     (0.9 by default) and, optionally, the linearized constraints are
     satisfied.
   - It recomputes the bound multipliers on the way out.
   - sqpopt's Gauss-Newton step plus escape probe is a lightweight
     stand-in. Uno's design is a clean template for replacing it, and it
     would make "infeasible" verdicts more reliable.
2. **The funnel method** (Kiessling, Leyffer and Vanaret). It replaces the
   filter's list of points with one number, a shrinking bound on the
   violation. Steps that reduce the objective get an Armijo test; steps
   that reduce the violation shrink the bound. It is about 150 lines in
   Uno and fits sqpopt as another line-search mode.
3. **Finer termination statuses.** Uno separates KKT points, Fritz–John
   points (stationary where the constraint gradients are degenerate),
   infeasible stationary points, and feasible or infeasible "small step"
   points. sqpopt folds Fritz–John points into success or stalled.
   Detecting them (the objective multiplier going to 0) would give users a
   better diagnosis.
4. **Protecting the step-acceptance test from roundoff.** Uno's `ipopt`
   preset adds about 10·ε·|f| to the actual reduction (an IPOPT trick).
   Near convergence this stops the filter or merit test from rejecting
   good steps because of cancellation. It is a few lines.
5. **Presets.** sqpopt's configurations are spread across several objects
   (`options%linesearch_mode`, `merit_mode`, `penalty_update`,
   `linesearch%interpolate`, `qp_solver%sparse_qp%null_space`, …). A
   preset option (for example `default`, `slsqp-like`, `robust`, `large`)
   would make the tested combinations from the Performance table
   one-liners.
6. **Filter details.** A cap on the filter's size, and filter resets after
   a number of rejections (with a limit on how many). Both are cheap.

Longer term, only if the scope grows:

- Exact Hessians with inertia correction (F7), and an interior-point path
  for large problems with many inequalities. Uno shows how to structure
  both, but they need a sparse LDLᵀ solver, i.e. an external dependency.
- Uno's Woodbury solver (L-BFGS as diagonal plus low-rank, used with a
  sparse factorization) is the matching trick if sqpopt ever adds a direct
  sparse KKT solver.

Not worth copying: the full eight-ingredient C++ class hierarchy. The
Fortran module split is adequate; the one refactor worth doing is making
restoration its own strategy, as part of idea 1.

## Suggested order

Idea 1, then ideas 2–4, measuring each on the HS suite as with the earlier
options.

## Status

- **Idea 1 (feasibility restoration phase):** done; see "F2 feasibility
  restoration phase" in `plan/ROADMAP.md`'s Phase 4 status.
- **Idea 2 (funnel method):** see the "Funnel method" section of
  `plan/ROADMAP.md`'s Phase 4 status.
