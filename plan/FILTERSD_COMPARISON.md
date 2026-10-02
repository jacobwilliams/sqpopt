# Comparison with filterSD

filterSD is Roger Fletcher's last nonlinear programming code (University of
Dundee, 2010-2011; Fortran 77, Eclipse Public License). The copy examined
is <https://github.com/jacobwilliams/FilterSD>, a start on a modern Fortran
version of it. This note compares it with sqpopt and lists the ideas worth
exploring.

Sources: the source (`filterSD.f90`, `glcpd.f90`, `l1sold.f90`, `checkd.f90`,
`schurQR.f90`; about 14,000 lines) and Fletcher's notes `docs/filterSD.pdf`
and `docs/glcpd.pdf`. filterSD was not built or run for this comparison.

## Summary

filterSD is not an SQP method. Its aims, in Fletcher's words, are "to avoid
the use of second derivatives, and to avoid storing an approximate reduced
Hessian matrix". Each major iteration minimizes the nonlinear objective
itself, subject to the linearized constraints and a trust region
(Robinson's method), and the minimization uses only gradients. So it has no
quasi-Newton matrix to offer, and its subproblem is a different one. What
it does offer sqpopt is at the edges: a derivative checker, a report of
which constraints are infeasible, and some trust-region rules.

## How filterSD works

- **Subproblem.** A linearly constrained problem (LCP): minimize `f(x)`
  subject to the linearized constraints and a box of radius `rho`. It is
  solved by `glcpd`, which calls the user's gradient many times within one
  major iteration.
- **No Hessian.** `glcpd` minimizes in the null space of the active
  constraints by a limited-memory steepest-descent method: a few (6 or 7)
  recent reduced gradients give Ritz values, estimates of the reduced
  Hessian's eigenvalues, whose inverses are the next step lengths. The Ritz
  values are passed from one subproblem to the next, and can be passed from
  one run to the next.
- **Active set.** A recursive active-set method with Wolfe's method for
  degeneracy and steepest-edge coefficients, on sparse LU factors updated
  by a Schur complement (`schurQR`) or by Fletcher-Matthews updates.
- **Globalization.** A filter with a trust region. An LP is solved first:
  it tells whether the linearized constraints are compatible with the trust
  region, and gives multiplier estimates. A rejected step is followed by a
  projection step (toward the constraints) before the radius is reduced.
  On a rejection the radius becomes a fraction (0.1 to 0.5, by the ratio of
  the violations) of the *step's* length; after an accepted step that
  reached the boundary it doubles.
- **Feasibility restoration.** If the LP is infeasible, an l1 feasibility
  problem is solved (phase 1), with the constraints split into those that
  are violated and relaxed, and the others. A post-processor (`l1sold`)
  finds the best l1 solution, so that the code converges fast to a point of
  local infeasibility, and reports which side of which constraint is
  infeasible there (`cstype`).
- **Safeguards.** A bound `ubd` on the constraint violation of any iterate,
  a lower bound `fmin` to detect an unbounded problem, and a limit `mxgr`
  on the gradient calls of one subproblem.
- **Derivative checker.** `checkd` checks the user's first derivatives
  before a solve: for each variable, the difference quotient over a
  user-given interval must lie between the derivatives at the two ends of
  the interval. No tolerance is needed.

Fletcher's own assessment (section 8 of his notes): when the null space is
not small, the lack of reduced Hessian information costs many more gradient
evaluations than codes with an approximate or exact reduced Hessian.

## Comparison

| | filterSD | sqpopt |
|---|---|---|
| method | Robinson (nonlinear objective, linearized constraints) | SQP (quadratic model) |
| second-order information | none: Ritz values of the reduced Hessian | L-BFGS, L-SR1, or the exact Hessian |
| globalization | filter with a trust region | filter or funnel line search (default), merit functions, trust region |
| gradient calls per major iteration | many | one |
| feasibility | l1 phase 1, with the infeasible constraints reported | elastic QP, restoration phase; `sqpopt_infeasible` |
| derivative checking | `checkd` | none in the library |
| safeguards | `ubd`, `fmin`, `mxgr` | `theta_max` of the filter, `obj_lower_limit`, QP iteration limits |

## Ideas worth exploring

1. **A derivative checker in the library.** The most useful thing here.
   Wrong derivatives are the most common cause of a solver "failing", and
   sqpopt has nothing for the user (its tests have their own checks). A
   routine such as `problem%check_derivatives(x, ...)` would compare
   `gjac` (and `hess`) with differences of `fc` (and of `gjac`), using the
   sparsity patterns, and report the entries that disagree. Fletcher's test
   is attractive because it needs no tolerance: the difference quotient
   over an interval must lie between the derivative values at its ends
   (true for a function that is monotone in its derivative there, so it can
   give a false alarm on an interval with an inflection, which the report
   should say). The Python bindings should have it too.

2. **Say which constraints are infeasible.** When sqpopt stops with
   `sqpopt_infeasible`, the user gets the point and the constraint values,
   and must work out which constraints are violated. filterSD returns that
   directly (`cstype`). The printed summary could list the violated
   constraints (index, side, amount), and the results could carry a count.
   Small.

3. **Start a solve from the previous solve's quasi-Newton approximation.**
   filterSD accepts the Ritz values of a previous run. The analogue for
   sqpopt is to keep the L-BFGS pairs: when a sequence of closely related
   problems is solved (a trajectory re-optimized with new data), each solve
   now starts from the identity. It must be an explicit request, since
   every `solve` starts from the configuration given to `initialize` by
   design, and the pairs are only valid if the problem changed little.

4. **The trust region's radius after a rejection.** sqpopt halves the
   radius; filterSD sets it to a fraction of the rejected step's length,
   which also cuts a step that was shorter than the radius. *Tried
   2026-10-02* in its simplest form (`radius = shrink_factor * min(radius,
   |p|_inf)`), on the HS suite: with the filter 261/40/4 with 8,235 `fc`
   (now 263/39/3, 9,251); with the funnel 265/38/2 with 11,443 (now
   263/40/2, 10,676); with the augmented Lagrangian merit function
   252/40/13 with 13,626 (now 256/39/10, 17,633). Mixed, so not adopted.
   filterSD's fraction also depends on how much the violation grew, and it
   takes a projection step first; neither was tried. The trust region is
   sqpopt's weakest globalization on this suite, so this may be worth a
   more careful look.

Not worth copying: the Hessian-free subproblem (Fletcher's own notes say
what it costs, and sqpopt's L-BFGS with the unconstrained step already
handles problems without second derivatives), the LP before each
subproblem (the elastic QP answers the same question), and the recursive
degeneracy scheme (no cycling has been seen in sqpopt's QP solvers, which
have their own safeguard).

## Suggested order

1, then 2. 3 if a use for it comes up. 4 only as part of a wider look at
the trust region.

## Status

Written 2026-10-02. Nothing here is implemented.
