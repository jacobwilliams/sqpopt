# Comparison with Yukon (GMAT)

Yukon is the SQP optimizer of NASA's General Mission Analysis Tool
([GMAT](https://sourceforge.net/projects/gmat/), Apache 2.0 license), in its
`YukonOptimizerPlugin`. It was written by Steven Hughes (NASA GSFC) in 2017
and converted to C++ by Joshua Raymond (Thinking Systems). This note compares
it with sqpopt and lists the ideas worth exploring.

Sources: the plugin's source code (the copy in GMAT R2026a: `Yukon.cpp`,
`MinQP.cpp`, `NLPFunctionGenerator.cpp`, and the GMAT glue `Yukonad.cpp`,
about 9,700 lines). Yukon was not built or run for this comparison, so
nothing here is a measurement.

## Summary

Yukon is a small dense SQP for the problems GMAT poses: a handful to a few
dozen variables, functions that are whole mission-sequence runs, and
derivatives by finite differences. sqpopt covers far more ground (sparse
problems, limited-memory and exact Hessians, filter and funnel methods, a
trust region, restoration, sparse factorizations). What Yukon has that
sqpopt lacks comes from the application, not the algorithm: it is driven by
its host one function evaluation at a time, it limits the step of each
variable separately, and it asks nothing of the user but function values.
Those three are the ideas worth exploring.

## How Yukon works

- **Reverse communication.** The optimizer is a state machine
  (`CheckStatus` says what it needs, `RespondToData` takes the values and
  advances: "ReadyForLineSearch", "StepTaken", "LineSearchConverged", ...).
  It never calls a user function. GMAT needs that: a function evaluation is
  a run of the mission sequence, which GMAT's own solver loop controls.
- **QP subproblem.** Its own dense active-set solver, `MinQP` (Nocedal and
  Wright's algorithm 16.1): a phase I that minimizes the largest constraint
  violation, a hot start from the previous active set, and a test for
  linearly dependent constraints. Dependent constraints are removed from the
  problem for the rest of the solve.
- **Hessian.** A dense BFGS matrix, with Powell's damping (thresholds 0.1
  and 0.9 instead of 0.2 and 0.8) or, the default, a self-scaled BFGS
  update: when the new curvature `s'y` is below `s'Bs`, the matrix is
  scaled by their ratio before the update. No limited memory, no exact
  Hessian.
- **Globalization.** A line search on the l1 merit function with one
  penalty parameter per constraint, raised from the QP's multipliers when
  the predicted reduction is too small. Rejected steps are shortened by an
  interpolation, at most 20 times; then the Hessian is reset to the
  identity, and a second failure ends the solve.
- **Non-monotone steps.** After a step with a sufficient decrease, the next
  steps are accepted without the test. If three in a row fail to improve on
  the best merit value, the solver goes back to the best point, with the
  Hessian saved there, and requires monotone steps for 10 iterations. This
  is a form of the watchdog technique.
- **Elastic mode.** When the QP fails, the whole problem is reformulated
  with elastic variables added to the decision vector (as SNOPT does), for
  the rest of the solve. Their weight is multiplied by 10 while they stay
  positive, up to `MaximumElasticWeight`; with the weight at its maximum
  and elastic variables still positive, the problem is reported infeasible.
- **Bounds and step limits.** Variable bounds are passed to the QP as
  general constraints. The whole step is then scaled back so that no
  variable moves by more than its own maximum step (GMAT's `MaxStep` of
  each `Vary` command) or beyond its bounds.
- **Convergence.** Feasible, and either the gradient of the Lagrangian is
  below a tolerance or the relative change of the objective is.

## Comparison

| | Yukon | sqpopt |
|---|---|---|
| problem size | small, dense | small dense to a million variables, sparse |
| interface | reverse communication | callbacks |
| derivatives | finite differences, by the host | supplied by the user |
| QP | dense active set, phase I | dense and sparse active set, elastic; direct (MUMPS); unconstrained step |
| Hessian | dense BFGS (damped or self-scaled) | L-BFGS, L-SR1, exact with a shift or inertia control |
| globalization | l1 merit line search, non-monotone | filter (default), funnel, merit (l1, augmented Lagrangian), watchdog, trust region |
| penalty | one per constraint | one for all constraints |
| second-order correction | no | yes |
| infeasible subproblem | elastic variables in the problem | elastic QP, then a restoration phase |
| step limit | per variable | one cap on the step's length, adapted |
| scaling | none in the optimizer | gradient-based, automatic |

sqpopt is ahead on every algorithmic row. The rows where Yukon does
something sqpopt doesn't are the interface, the derivatives, the penalty,
and the step limit.

## Ideas worth exploring

1. **A reverse-communication interface.** The most valuable, and the most
   work. A host like GMAT, or any code where the "function" is a simulation
   that the host must drive (or run in another process, or in parallel),
   can't hand sqpopt a callback. `sqpopt_iterate` already does one major
   iteration per call, but the line search, the second-order correction,
   and the restoration call the user functions from inside it. A
   reverse-communication layer would have to suspend at each of those
   calls: either those loops are rewritten as a state machine, or the
   solver runs in its own thread and hands over at each evaluation. It
   needs a design of its own before any code.

2. **A maximum step for each variable** (`MaxStep`). sqpopt has one cap on
   the step's 2-norm (`qp_solver%max_step`, adapted like a trust radius).
   In trajectory problems the variables have different units and
   sensitivities (a burn component in km/s, an epoch in days), and the
   user often knows how far each can move before the linearization is
   useless. A vector of limits is easy to add where it belongs: as bounds
   on the QP's step (`x_lb - x <= p <= x_ub - x` intersected with
   `+-max_var_step`), so the QP's solution respects them, instead of
   scaling the whole step back afterwards as Yukon does. The trust region
   already intersects the bounds with a box in this way.

3. **Finite-difference derivatives in the library.** sqpopt requires
   `gjac`. The HS harness has its own differencing for the 16 problems
   without analytic derivatives, and the Python bindings require `jac`
   (unlike scipy). A built-in option (forward or central, with a relative
   and an absolute perturbation for each variable, respecting the bounds)
   would let a user supply `fc` alone, as Yukon's users do, and would use
   the existing `derivative_accuracy` switch (forward differences far from
   the solution, central near it). With the sparsity pattern known, the
   columns can be grouped so that a sparse Jacobian needs only a few
   evaluations.

4. **One penalty parameter per constraint**, for the merit-function line
   searches. With constraints of very different scales, a single penalty
   is set by the worst one and over-penalizes the rest. The filter, which
   is the default, has no penalty, so this only matters to the merit
   modes, and the gradient-based scaling already evens the constraints out
   at the starting point. Low priority.

5. **Keep the Hessian of the best point in the watchdog.** sqpopt's
   watchdog stores the best point (`watchdog_x_opt`); Yukon also stores
   the Hessian approximation there and restores both. Returning to an old
   point with a quasi-Newton matrix built from later, rejected steps is a
   mismatch. For L-BFGS that means saving the pairs (or their count at the
   best point). A small change to try on the watchdog rows of the HS
   table.

6. **Self-scaled BFGS.** sqpopt's L-BFGS already rescales its initial
   matrix from the latest pair at every update, which does the same job
   for a limited-memory method. Nothing to add, unless a dense BFGS mode
   is ever wanted for very small problems.

Not worth copying: removing linearly dependent constraints for the rest
of the solve (dependence is a property of the current point; sqpopt's QP
solvers handle it at each iteration), the phase I on the largest violation
(the elastic l1 form does the same job and gives multipliers), and elastic
variables in the problem itself (the restoration phase covers that case
and can report local infeasibility).

## Suggested order

2 (small, self-contained, useful for the trajectory problems sqpopt is
meant for), then 3, then 5 as an experiment. 1 needs a design of its own
first. 4 only if the merit modes get more use.

## Status

Written 2026-10-02. Idea 2 is implemented (2026-10-02): `problem%set_max_step`, a limit on the change of
each variable per major iteration, enforced as bounds on every step (the QP's, the second-order
correction, the restoration steps, and the trust region's box); `max_step=` in the Python bindings;
tests `test_max_step` and `test_minimize.py`. The others are not.
