# Comparison with OpenSQP

OpenSQP is an SQP method written in Python by Anugrah Jo Joshy and John
T. Hwang (UC San Diego), distributed in their
[modOpt](https://github.com/LSDOlab/modopt) library. This note compares it
with sqpopt and lists the ideas worth exploring.

Source: the paper, "OpenSQP: A Reconfigurable Open-Source SQP Algorithm in
Python for Nonlinear Optimization" (arXiv:2512.05392v1, December 2025;
`references/2512.05392v1.pdf`). The code was not read or run, so what is
said about OpenSQP here is what the paper says. The statements about sqpopt
marked "checked" were run for this note.

## Summary

OpenSQP is a textbook dense line-search SQP, built from parts that can be
swapped: a damped BFGS matrix, a dense QP solver, and SNOPT's augmented
Lagrangian merit function. It is meant for problems of up to about 100
variables, and its authors list as future work what sqpopt already has (a
limited-memory Hessian, exact Hessians, sparse QP solvers). So there is no
algorithm in it for sqpopt to adopt. What is worth taking is how it was
tested: 575 CUTEst problems through PyCUTEst, against SNOPT, IPOPT, SLSQP,
and scipy's `trust-constr`, reported as success rates and performance
profiles. sqpopt has only the 305 Hock-Schittkowski problems and its own
large examples, and has never been compared with SNOPT or IPOPT.

## How OpenSQP works

- **Hessian.** A dense BFGS matrix with Powell's damping, started from the
  identity, and reset to the identity if it becomes indefinite by
  round-off.
- **QP.** The dense dual active-set method of Goldfarb and Idnani
  (`quadprog`). A dual method needs no feasibility phase, which the paper
  argues is what costs most in the first iterations, but it can't be
  warm-started from the previous working set.
- **Inconsistent QP.** Powell's device, as in SLSQP: one extra variable
  `eta` in `[0, 1]` that scales the violated constraints toward
  consistency, with the penalty `gamma*eta**2/2` (`gamma` from `1e6`, times
  10 after 25 inconsistent iterations in a row, up to `1e12`). That QP is
  solved by HiGHS. Equality constraints are first split into two
  inequalities, which the paper says is what lets it solve some problems
  with more equality constraints than variables.
- **Globalization.** A line search on the augmented Lagrangian of Gill,
  Murray, and Saunders, in the variables, the multipliers, and slacks,
  with one penalty parameter per constraint (the smallest vector that
  gives the required descent, and SNOPT's rule that lets it decrease
  slowly). The step satisfies the strong Wolfe conditions (MINPACK's
  search). If that search fails, a backtracking search on the `l1` merit
  function is tried instead.
- **Starting point.** The user's point is projected onto the bounds. If
  the functions can't be evaluated there, the user's point itself is used.
- **Evaluation failures.** A trial point where a function fails is backed
  away from, along the step.
- **Convergence.** SNOPT's test: the tolerances are multiplied by
  `1 + |x|_inf` (feasibility) and `1 + |lambda|_inf` (optimality).

## Comparison

| | OpenSQP | sqpopt |
|---|---|---|
| language | Python (NumPy; QP solvers in C and C++) | Fortran, with Python bindings |
| problem size | up to about 100 variables in the paper | small dense to a million variables, sparse |
| Hessian | dense damped BFGS | L-BFGS (damped), L-SR1, exact with a shift or inertia control |
| QP | dense dual active set; HiGHS when inconsistent | dense and sparse primal active set, elastic; direct (MUMPS); unconstrained step |
| inconsistent QP | one variable for all the constraints (Powell) | an elastic slack per violated constraint, then a restoration phase |
| globalization | augmented Lagrangian line search (strong Wolfe), `l1` as a fallback | filter (default), funnel, `l1` and augmented Lagrangian merit, watchdog, trust region |
| penalty | one per constraint | one for all constraints |
| second-order correction | no (the smooth merit function doesn't need one) | yes |
| scaling | none | gradient-based, automatic |
| convergence test | tolerances relative to `|x|` and `|lambda|` (SNOPT) | IPOPT's scaling by the mean multiplier, on the scaled problem |
| tested on | 575 CUTEst problems, against four solvers | 305 HS problems against SLSQP and NLPQLP's published counts; large examples |

## What the paper measured

575 CUTEst problems with at most 100 variables and 100 constraints, at
most 250 iterations, full-memory quasi-Newton Hessians in every solver
(no second derivatives), and the tolerances of Gill, Saunders, and Wong's
benchmark (optimality `1.22e-4`, feasibility `2e-6`).

| solver | solved |
|---|---|
| SNOPT | 479 (83.3%) |
| OpenSQP | 479 (83.3%) |
| IPOPT (limited-memory Hessian) | 475 (82.6%) |
| SLSQP | 463 (80.5%) |
| scipy `trust-constr` | 378 (65.7%) |

SLSQP needed the fewest function evaluations on the problems it solved.
43 of the problems have more equality constraints than variables: IPOPT,
SLSQP, and `trust-constr` solved none of them, and SNOPT and OpenSQP about
a fifth.

## Ideas worth exploring

1. **A CUTEst benchmark through PyCUTEst and the Python bindings, with
   performance profiles.** The most valuable, and now cheap, because sqpopt
   has a `minimize` for Python. The roadmap's CUTEst item plans a
   pure-Fortran harness, translating the SIF files (several sessions of
   work, and deferred). PyCUTEst gives the same problems with their
   derivatives to any Python solver, so a script could run sqpopt, SLSQP,
   `trust-constr`, and IPOPT (cyipopt) on the paper's set, with its limits
   and tolerances, and its Table 2 would then be a published reference for
   SNOPT, which we can't run. The outputs would be the success counts, a
   performance profile (Dolan and Moré: the fraction of problems solved
   within a factor of the best solver's time) and a data profile (the same
   for evaluations), for the Performance page. It would also show how
   sqpopt does on kinds of problems the HS set has few of: the 43
   overdetermined ones, unbounded ones, and ones whose functions are
   undefined in places. It needs CUTEst and its problem files installed
   (PyCUTEst builds each problem with a Fortran compiler), so it is a tool
   for the developer, not part of `fpm test`. The pure-Fortran harness
   stays the way to a regression test that needs nothing installed.

2. **If the functions can't be evaluated at the projected starting
   point.** sqpopt projects `x0` onto the bounds and stops with
   `sqpopt_function_error` if a function fails there (checked: a start
   outside its bounds whose projection is on the edge of a logarithm's
   domain ends with status 25). OpenSQP goes back to the user's point.
   sqpopt can't do that, because it promises never to evaluate outside the
   bounds. What it could do is move the projected point a little inside
   the bounds that it is on (as interior-point codes push the starting
   point off the bounds) and try again. Small, and only for this case; the
   diagnostics' probe (level 3) already reports a start at the edge of the
   functions' domain.

3. **One penalty parameter per constraint**, for the augmented Lagrangian
   merit function. The same idea as in the Yukon comparison: with
   constraints of very different sizes, one penalty is set by the worst
   and over-penalizes the rest. OpenSQP's rule is SNOPT's: the vector of
   smallest norm that gives the descent the line search needs, and a
   damped decrease. sqpopt's filter, the default, has no penalty, and the
   gradient-based scaling evens the constraints out at the start, so this
   only matters to the merit modes. Low priority, as before; this is a
   second solver that does it.

4. **Another globalization as a fallback.** When its line search fails,
   OpenSQP tries a different one (a simpler merit function and test)
   before giving up. sqpopt resets the Hessian and, with the filter, takes
   a restoration step, and stops after `max_consecutive_failures`
   iterations. Trying the `l1` merit search on a step the filter rejected
   would be easy to add. But there is nothing on the HS set to measure it
   with: none of the 305 problems ends with a line-search failure in the
   default configuration. Worth a look only if the CUTEst run (idea 1)
   shows such failures.

## Not worth copying

- **The dual active-set QP.** Its advantage is having no feasibility
  phase when there is no warm start. sqpopt's QP solvers start from a
  crash basis or the previous working set and use elastic slacks instead
  of a separate phase, and on problems this small the QP's time doesn't
  matter (the whole HS suite takes 0.7 s).
- **The strong Wolfe line search.** Its use is to keep the BFGS update's
  curvature positive without damping, on unconstrained problems. sqpopt's
  L-BFGS already takes about as many iterations as scipy's L-BFGS-B,
  which uses the same MINPACK search (273 against 269 on the chained
  Rosenbrock function with 50 variables, 1,524 against 1,502 with 300:
  `python/tests/test_scipy_compare.py`).
- **Powell's single variable for an inconsistent QP.** It relaxes every
  violated constraint by the same fraction. An elastic slack per
  constraint, which sqpopt has, relaxes only the ones that have to be.
- **Splitting equalities into two inequalities.** The paper credits this
  for its results on overdetermined problems. sqpopt's elastic slacks
  already relax an equality on the side it is violated, and it solves
  small overdetermined problems that are consistent (checked: three
  nonlinear equalities in two variables, from three starting points, and
  four linear ones in three variables). Whether it does as well on
  CUTEst's 43 is part of idea 1.
- **SNOPT's relative tolerances.** A feasibility tolerance that grows
  with `|x|` is looser for problems with large variables. sqpopt's own
  additions to its test went the other way (`dual_inf_tol`, and the
  objective-change test of the acceptable stop, make it stricter in the
  original units).
- **Swappable components.** sqpopt's components are already separate
  types set at `initialize`. Letting user code replace one would mean
  making them abstract types, which is a redesign with no request behind
  it.

## Suggested order

1, as a script under `tools/` or `python/`, first on a few dozen problems
to see that PyCUTEst and the bindings work together. 2 if a CUTEst problem
shows the need. 3 and 4 only with evidence from 1.

## Status

Written 2026-10-02. Idea 1 is started (2026-10-02): `tools/cutest_benchmark.py`
runs sqpopt, SLSQP, and `trust-constr` on 664 CUTEst problems through
PyCUTEst; sqpopt reports success on 525 of them, SLSQP on 493, and
`trust-constr` on 412. The numbers, and what they say to look at next, are
in the roadmap ("CUTEst through PyCUTEst and the Python bindings"). IPOPT
and the performance profiles are not done. The other ideas are not
implemented.
