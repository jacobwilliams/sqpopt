# Comparison with OPTGRA

OPTGRA is the European Space Agency's optimizer for trajectory design,
written by Johannes Schoenmaekers (2008, Fortran 77), and distributed with
[pyoptgra](https://github.com/esa/pyoptgra). The copy examined is
<https://github.com/jacobwilliams/optgra>, a modern Fortran refactoring
(one module, `src/optgra.F90`, 3,300 lines). This note compares it with
sqpopt and lists the ideas worth exploring.

Sources: that module's code and comments, and its README. OPTGRA was not
built or run for this comparison.

## Summary

OPTGRA is a gradient-projection method "for near-linear optimization
problems with many constraints": it first corrects the constraints, then
moves along the boundary of the feasible region to improve the objective,
correcting violations as they appear. It is a small dense code with no
second-order model beyond conjugate-gradient directions, and its
documentation says it does less well on very nonlinear problems. As an
algorithm it has nothing sqpopt lacks. As a tool built for trajectory
designers it has several conveniences that sqpopt doesn't: scale factors
for each variable and constraint, a sensitivity analysis of the solution,
finite differences built in, and names in its output.

## How OPTGRA works

- **Two phases.** A correction phase (`ogcorr`) brings the iterate to the
  linearized constraints, with an active set of constraints that grows and
  shrinks (`ogincl`, `ogexcl`). An optimization phase (`ogopti`) then moves
  in the null space of the active constraints.
- **Search direction.** Steepest descent, conjugate gradients, or a
  spectral conjugate-gradient method (`optmet`), on the reduced gradient. A
  second-order estimate along the direction comes from a finite difference
  with the perturbation `varsnd`. No Hessian approximation is stored.
- **Step limit.** A maximum distance per iteration in scaled variables
  (`varmax`).
- **Scaling by the user.** Each variable has a scale factor (`varsca`) and
  each constraint a convergence threshold (`consca`, `delcon`). The user's
  functions work in physical units; the solver divides by the scales.
- **Constraint priorities** (`conpri`): during the first correction phase,
  iteration `k` only considers the constraints of priority up to `k`, so
  the important ones are satisfied first.
- **Derived data.** A "constraint" of type -2 is evaluated and reported
  but not enforced: a quantity the designer wants to see at the solution.
- **Derivatives.** Supplied by the user, or by forward or central
  differences with a perturbation for each variable (`varder`, `varper`).
- **Sensitivity analysis** (`ogsens`): after a solve, the sensitivities of
  the constraints and the objective to the active constraints and to
  designated parameters, and of the variables to both (variables of type 1
  are parameters, held fixed in the optimization). A "sensitivity
  optimization" mode (`senopt`) then re-optimizes from the stored state
  after a change. The state can be saved and restored (`oggsst`, `ogssst`).
- **Output.** Names for the variables and constraints, a log, a table of
  the iterations on its own unit, and printing every N-th iteration.
- **Maximization.** The objective's type says whether to minimize or
  maximize.

## Comparison

| | OPTGRA | sqpopt |
|---|---|---|
| method | gradient projection, CG directions | SQP |
| second-order information | a finite difference along the direction | L-BFGS, L-SR1, or the exact Hessian |
| problem size | small, dense | small dense to a million variables, sparse |
| globalization | correction, then a bounded step | filter or funnel line search, merit functions, trust region |
| scaling | by the user: each variable and constraint | automatic, of the objective and constraints; the variables are the user's job |
| derivatives | user, or forward or central differences | user |
| sensitivities | of the solution, to constraints and parameters | the multipliers only |
| output | named variables and constraints | numbered |

## Ideas worth exploring

1. **Scale factors for the variables.** sqpopt scales the objective and
   the constraints, but not the variables: the guide tells the user to
   solve in `y = x/d` and shows the wrappers to write (`fc_y`, `gjac_y`).
   Yet much of the solver depends on the variables' scale: the initial
   L-BFGS matrix, the step cap, the trust region, the restoration's
   proximal term. A `problem%set_variable_scaling(d)` would do the
   wrappers' work inside the evaluation layer (which already applies the
   objective and constraint scale factors): the solver works in `y`, and
   the user's functions, bounds, starting point, step limits, and results
   stay in `x`. This was already noted in the roadmap (§4, "optional user
   variable scaling"). For trajectory problems, where the variables have
   very different units, it is the idea here with the most effect.

2. **Sensitivity of the solution.** sqpopt returns the multipliers, which
   are the sensitivities of the optimal objective to the constraint
   bounds. OPTGRA also gives the sensitivities of the *variables*, and to
   chosen parameters. At a solution these come from one solve with the KKT
   matrix of the active set for each parameter, and sqpopt's
   `sqpopt_kkt_module` can already factor that matrix (in a build with
   MUMPS, with the exact Hessian or the quasi-Newton one). A routine called
   after `solve` could return `dx/db` for the active constraints' bounds,
   and `dx/dp` for parameters given the derivatives of the gradient and
   the constraints with respect to them. It is what a mission designer
   uses to see what a constraint costs, or to update a solution after a
   small change without re-optimizing.

3. **Finite-difference derivatives in the library**, with a perturbation
   for each variable. The same idea as from Yukon (see
   [YUKON_COMPARISON.md](YUKON_COMPARISON.md), idea 3); OPTGRA is a second
   trajectory tool that has it.

4. **Names for the variables and constraints.** Optional names, used in
   the detailed log, in the final table of the variables, and in messages
   ("constraint 17" becomes "periapsis altitude"). Small, and it makes the
   output of a real problem readable.

5. **Smaller conveniences.** A flag to maximize the objective instead of
   negating it by hand; printing every N-th iteration; an iteration table
   on its own unit (the `report` callback already lets a user write one).

Not worth copying: the method itself; constraint priorities (the
restoration phase and the elastic QP treat the constraints together, and a
priority order would need a different restoration); and "derived data"
constraints, which sqpopt already expresses as a constraint with infinite
bounds (though it would be worth checking that such free rows cost
nothing in the QP).

## Suggested order

1, then 4 (small). 2 and 3 are larger and independent. 5 as wanted.

## Status

Written 2026-10-02. Nothing here is implemented.
