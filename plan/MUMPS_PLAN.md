# More uses of MUMPS: a plan for large problems

*Status: implemented 2026-10-01, except the thread and ordering study (G3) and the two proposals that were to wait for evidence (M4, M5). §8 says what was built, what was measured, and where the implementation differs from this plan. The plan itself (§1–§7) is kept as written 2026-10-01, after the inertia control of [INERTIA_CONTROL.md](INERTIA_CONTROL.md) §9 was built.*

MUMPS is in the library to help with **large, sparse problems**. The inertia control uses it only to count negative pivots: it factors the KKT matrix, reads one number, and discards the factors. This document proposes what else to do with a sparse symmetric factorization, in the order that looks most useful, and how to tell whether each step paid off.

Everything here stays optional, behind `HAS_MUMPS`, and the default build must keep needing only fpm.

## 1. Where the time goes on large problems

With the exact Hessian, almost all of a large solve is the QP subproblem. One run each, release build, `test_large_sparse`'s problems at larger sizes (2026-10-01):

| problem | n | inertia control | iterations | time (s) | of which QP (s) | factorizations |
|---|--:|---|--:|--:|--:|--:|
| control | 10,001 | off | 9 | 0.444 | 0.439 | – |
| | | on | 9 | 0.441 | 0.386 | 14 |
| control | 100,001 | off | 11 | 26.5 | 26.4 | – |
| | | on | 10 | 22.3 | 21.7 | 14 |
| rosenbrock | 20,000 | off | 16 | 31.2 | 31.1 | – |
| | | on | 17 | 21.1 | 21.0 | 20 |

- The QP takes 97% of the time or more.
- A factorization is cheap next to it: about 40 ms at order 150,000 (`control`, n = 100,001), and about 8 ms at order 30,000 (`rosenbrock`, n = 20,000).
- QP time grows much faster than the problem: 60 times more for 10 times the variables on `control`.

ROADMAP.md (F1, "Where the time goes now") says what the QP time is: active-set iterations that change one constraint at a time, each with conjugate-gradient (CG) iterations on the reduced Hessian ZᵀHZ. On `rosenbrock` N = 6000, the Hessian products inside CG were about 70% of the time.

So a factorization can only make a large solve faster by replacing work inside the QP. That is the main proposal (M1). The others are about robustness, or about bringing M1 to the other Hessian modes.

## 2. Groundwork (needed by every proposal)

**G1. A general sparse symmetric solver type.** `sqpopt_inertia_type` mixes two things: the MUMPS instance, and the search for the shift. Split out a `sqpopt_kkt_type` (name open) that owns the MUMPS instance and offers:

- `initialize` (the sparsity pattern), `factor` (values in; negative and null pivot counts out), `solve` (one or several right-hand sides), `destroy`;
- the fixed-pattern convention the inertia control already uses: the matrix always has order `n+m`, a variable at a bound gets the identity's row and column, a constraint outside the working set minus the identity's. The pattern is then analysed once per solve.

The inertia control becomes a user of this type. MUMPS stays referenced in one file only.

**G2. Solves, not only factorizations.** Today nothing calls MUMPS's solve phase. Add it, with one or two steps of iterative refinement done in working precision (`wp`) against the true residual. MUMPS is double precision only, so in a `REAL128` build the refinement is what recovers the accuracy; without it, direct steps would be limited to double precision.

**G3. Threads and ordering.** MUMPS is pinned to one thread (`ICNTL(16)=1`), because its threads made the 305 small HS problems take 24.5 s instead of 1.4 s. That choice has never been tested on a large problem. Measure threads (1, 2, 4, 8) and the fill-reducing ordering (AMD, METIS, SCOTCH: `ICNTL(7)`) on the benchmarks of §5, then either pick the thread count from the matrix order or make it an option.

**G4. Time and failure accounting.**
- Add `results%time_factorization` (factor and solve time), next to `time_functions` and `time_qp`, and print it in the summary. Without it the cost of every proposal below is hidden in "other".
- A MUMPS allocation failure should end the solve with `sqpopt_out_of_memory`, as the dense QP's does. Today a failed factorization silently turns inertia control off for the rest of the solve.

**G5. Benchmarks that exercise it.** `test_large_sparse`'s problems need at most one shift, so they show overhead, not benefit. Add to `example/benchmark.f90` (or a new example):
- the existing problems at n = 10⁴ and 10⁵;
- a large, strongly nonconvex problem, where the Hessian of the Lagrangian is indefinite at most iterates (for instance a chained problem with concave terms and bounds);
- a problem with many degrees of freedom (`elliptic` of `test_sparse_options`, scaled up), where CG has the most work.

## 3. Proposals

### M1. Solve the QP's face directly, with the factorization we already have *(exact Hessian; the main one)*

**What.** Before each QP, the inertia control factors

    K = [ H+δI  Jₐᵀ ]
        [ Jₐ    0   ]

for the working set the QP will start from. One solve with those factors gives the minimizer of the QP on that face, and its multipliers. This is INERTIA_CONTROL.md §5.2 and the "sparse KKT active-set QP" of ROADMAP.md F1.

**How.**
1. Right-hand side: `−g` for the free variables, corrected for the variables held at their bounds (`−H p_fixed`, `−Jₐ p_fixed`); the distance to its bound for each active row and each fixed variable; zero for the rows outside the working set.
2. Solve (G2). The result is the step `p` and the multipliers of the working set's rows; the bound multipliers follow from `g + Hp − Jᵀλ`.
3. **Accept** the step as the QP's solution if it satisfies every constraint and bound outside the working set, and every multiplier has the right sign (both to the QP's own tolerances). The QP solver is then not called at all.
4. **Otherwise** pass `p` to the existing QP solver as its starting step (a new optional argument of `solve_qp_subproblem`; the solvers already build a starting step from a guessed working set). Nothing is lost: the QP does what it does today, from a better point.

**Stage 2 (only if stage 1 pays off).** When the test in step 3 fails by a few rows, change the working set (add the most violated row, or drop the row with the worst multiplier), refactor, and solve again, up to a small number of times, before handing over. That is a primal active-set QP on the KKT matrix. It is only sensible while a refactorization is much cheaper than a QP solve, which §1 says it is.

**Expected benefit.** Near a solution the working set settles, so most major iterations would cost one factorization and one solve (tens of milliseconds at n = 10⁵) instead of a QP solve (seconds). The direct step is also exact, where CG stops at a tolerance.

**Limits and risks.**
- MUMPS can't update its factors when the working set changes: every change is a refactorization. On an infeasible start where thousands of constraints change status (ROADMAP.md: the first QP of `rosenbrock` N = 6000), stage 1 does nothing and stage 2 must give up quickly.
- When the linearized constraints are inconsistent the QP needs its elastic slacks. The direct step should only be tried when the previous QP used none.
- A rank-deficient working set makes K singular. MUMPS's null-pivot detection reports it; then skip the direct step.
- It changes the QP's solution only within the tolerances, but that can still change iteration counts. The HS suite (`--hessian=exact --inertia`) must be compared.

**Option.** `options%direct_qp_step` (logical, default off at first; requires `inertia_control`). Name open.

**Measure.** QP time and total time on G5's benchmarks, with and without; the share of major iterations whose QP was skipped (a new count in the results); HS suite counts unchanged or better.

### M2. Direct least-squares solves in the Jacobian *(any Hessian mode)*

**What.** Three places still use `LSQR`, an iterative solver, for minimum-norm or least-squares problems in J:

| where | problem | today's failure path |
|---|---|---|
| `restoration_step` ([sqpopt_restoration_module.f90](../src/sqpopt_restoration_module.f90)) | Gauss-Newton step: minimum-norm `J p = −r` | `LSQR` at its iteration limit: no restoration step |
| `soc_step` ([sqpopt_soc_module.f90](../src/sqpopt_soc_module.f90)) | second-order correction: minimum-norm on the active rows | `LSQR` at its iteration limit: no correction |
| `sqpopt_null_space_lsqr` (the sparse QP's fallback) | projections onto the null space | slow (2–200× the LU method) |

Each is one solve with

    [ I   Jₐᵀ ] [ p ]   [ 0 ]
    [ Jₐ  −εI ] [ y ] = [ r ]

where a small ε > 0 makes the matrix nonsingular when the rows are dependent (the solution is then the regularized least-squares one).

**Expected benefit.** Robustness more than speed: a direct solve doesn't hit an iteration limit, and its accuracy doesn't depend on the conditioning of J. I have not measured how much time these solves take on large problems, and expect it to be small next to the QP. Measure that first (G4); if it is negligible and the `LSQR` failure paths are never taken on the benchmarks, this proposal can wait.

The first two are worth doing; replacing the `LSQR` null-space method is not, since the LU method already replaced it.

**A by-product.** The same system gives least-squares multiplier estimates at the starting point (`min ‖g − Jᵀλ‖`), which the solver doesn't compute today (it starts from `λ = 0` unless `lambda0` is given).

**Option.** These don't depend on the Hessian mode, so they need their own switch, e.g. `options%linear_solver = sqpopt_linear_solver_iterative | sqpopt_linear_solver_mumps`. Name open.

**Measure.** `test_infeasible`, `test_maratos`, `test_degenerate`, and the HS suite in every configuration of the Performance table: outcomes unchanged or better, and the counts of restoration steps and corrections that failed for lack of a step.

### M3. Quasi-Newton Hessians through a low-rank correction *(L-BFGS and SR1)*

**What.** L-BFGS and SR1 are not sparse matrices, so M1 and the inertia control can't use them as they are. Both have the compact form `B = θI − W M⁻¹ Wᵀ`, with `W` of `2k` (BFGS) or `k` (SR1) columns. So the KKT matrix is a sparse matrix K₀ (with `θI` in place of H) minus a low-rank term, and:

- **solves** use the Sherman–Morrison–Woodbury formula: `2k+1` solves with K₀'s factors and one small dense system;
- **the inertia** follows from the inertia of K₀ and of that small dense matrix (the Haynsworth inertia additivity formula; to be checked carefully for the signs before coding).

This is INERTIA_CONTROL.md §5.3 and the long-term item of ROADMAP.md F1. IPOPT does the same for its L-BFGS option.

**Expected benefit.**
- **L-BFGS** (the default; most users have no second derivatives): it needs no inertia control, being positive definite, but it would get M1's direct face solve.
- **SR1:** an inertia correction, which it lacks entirely. It is much weaker than BFGS on the HS suite (244 solved, 30 failed), plausibly because its indefinite updates are never corrected.

**Limits and risks.**
- Cost: with the automatic memory of up to 100 pairs, that is up to 200 extra solves per factorization for L-BFGS. MUMPS solves several right-hand sides at once, but this still needs measuring; a short memory may be the practical choice with it.
- It is the most new machinery of the proposals, and only worth it if M1 pays off for the exact Hessian first.

**Measure.** As M1, with `--hessian=bfgs` and `--hessian=sr1`; SR1's HS counts against 244/31/30.

### M4. A constraint preconditioner for CG on the reduced Hessian

**What.** ROADMAP.md F1 lists "a preconditioner for ZᵀHZ" as open (CG has only a diagonal one). A factorization of the KKT matrix with a simple approximation of H (its diagonal, or `θI`) is the standard *constraint preconditioner* for projected CG.

**When.** Only if M1's refactorizations turn out too expensive on problems whose working set keeps changing: the preconditioner's matrix depends on the working set too, but it can be kept while the working set changes a little. M1 mostly subsumes it, so decide after M1.

### M5. Direct solves in the restoration phase

**What.** The feasibility QP of the restoration phase has the Hessian `ζI`, so its KKT matrix is the one of M2 scaled. The phase's QPs could use M1's direct step with it.

**When.** Small, and only matters on problems that spend many iterations in restoration. Do it only if a benchmark shows restoration-phase QPs taking real time.

## 4. Not proposed

- **The sparse QP's LUSOL basis.** It relies on column-replacement updates (`lu8rpc`), which MUMPS doesn't have.
- **Rank detection and the initial working set** (F13). One LUSOL factorization already does it.
- **The dense QP, and small problems generally.** On the HS suite a factorization is no faster than what is there, and MUMPS's fixed costs are comparable to the whole solve.
- **A second factorization library** (HSL MA57/MA97, PARDISO). G1's type is the place where one could be added later, but nothing here needs it.

## 5. Order of work

1. **G1–G5** (groundwork and benchmarks). Without the benchmarks and the time accounting, none of the rest can be judged.
2. **M1, stage 1.** The smallest change with the largest expected effect: it reuses a factorization that is already computed, and falls back on the current QP.
   - *Decision point:* if the QP is skipped in fewer than about half of the major iterations on the large benchmarks, or the total time doesn't drop, stop and write down why.
3. **M1, stage 2**, if stage 1's skipped share is high but the first iterations still dominate.
4. **M2** (restoration and correction steps), if G4's accounting or the failure counts justify it.
5. **M3**, if M1 paid off: first L-BFGS solves, then SR1 with its inertia.
6. **M4, M5** only on evidence.

## 6. Success measures

- **Large problems** (G5, n = 10⁴–10⁵): total time and QP time, against the table in §1. M1's target is that the QP no longer dominates once the working set has settled.
- **HS suite:** the default (L-BFGS) results must not change while the new options are off: 280 solved, 25 local, 0 failed, 9,173 `fc`. With them on, compare against 274/29/2 and 9,581 `fc` (exact Hessian with inertia control).
- **Both builds** pass every test, and each new option is rejected as invalid input without MUMPS, as `inertia_control` is.
- **Memory:** the factors' size (MUMPS reports it) on the largest benchmark, since fill is the known risk of a direct method (INERTIA_CONTROL.md §6).

## 7. Open questions

- **Option names and grouping.** One switch per feature (`direct_qp_step`, …), or one `linear_solver` option that turns on everything MUMPS can do? The second is simpler for users; the first is needed while the features are being evaluated.
- **Threads:** a fixed rule from the matrix order, or a user option (G3)?
- **How many refactorizations** M1's stage 2 may spend before handing over to the QP.
- **Regularization ε** in M2: fixed, or raised only when MUMPS reports null pivots?
- **`REAL32` and `REAL128` builds:** MUMPS has a single-precision library, but only the double one is used. Is refinement in `wp` (G2) enough for `REAL128`, or should direct steps be turned off there?

## 8. What was implemented and measured (2026-10-01)

### What was built

| plan item | status | where |
|---|---|---|
| G1 solver type | done | `sqpopt_symmetric_solver_module.F90` (the only file that refers to MUMPS) and `sqpopt_kkt_module` (the KKT matrix of a working set); `sqpopt_inertia_module` is now a user of them |
| G2 solves with refinement | done | `symmetric_solver_solve` (up to 2 refinement steps). `REAL128` was dropped instead: a build with MUMPS is double precision only (a compile error otherwise) |
| G3 threads and ordering | threads done (2026-10-02); the ordering **not done** | `options%factorization_threads` (default 1). Threads gave nothing on the banded benchmarks and on a 2-D grid, and 2.3× on a 3-D grid of order 216,000 with 4 threads. MUMPS's default ordering |
| G4 time and failure accounting | done | `results%time_factorization`, `n_qp_solves`, `n_direct_qp`; an out-of-memory factorization ends the solve with `sqpopt_out_of_memory` |
| G5 benchmarks | done | `example/benchmark_large.f90`: `control`, `rosenbrock`, and the nonconvex `wells` |
| M1 direct QP | done, both stages | `sqpopt_qp_direct_module`, called from `solve_qp_subproblem`; `options%direct_qp` |
| M2 direct least squares | done | `sqpopt_least_squares_module`; `options%direct_least_squares` |
| M3 quasi-Newton Hessians | done, both parts | the low-rank form in `sqpopt_kkt_module`; `direct_qp` with L-BFGS and SR1, and `inertia_control` with SR1 |
| M4 CG preconditioner | not done | no evidence for it: with the direct QP, the QP no longer dominates |
| M5 direct restoration-phase QPs | not done | the benchmarks never enter a restoration phase |

### Where it differs from the plan

- **M1 has no "starting step" for the active-set QP.** The plan's step 4 was to pass a rejected direct step to the active-set solver. Instead, stage 2 was built at once, as a primal-dual active-set method (add every violated row and bound, drop every wrongly-signed one, refactor), and if that gives up, the active-set solver starts as it always did.
- **A nonconvex face doesn't end the direct method.** With inertia control it raises the shift and continues. Without this, one fallback QP on `wells` at n = 100,000 took 517 s.
- **A singular face is regularized** (`-εI` in the rows' block) instead of ending the method. Adding every violated bound at once often leaves a row with no free variable; on `wells` at n = 100,000 that happened with 10,719 null pivots.
- **A stationarity check** on every direct solve (the accuracy of the solve), which the plan didn't have.
- **`direct_qp` doesn't require `inertia_control`.** It checks the face's convexity itself, and gives up on a nonconvex face without it.
- **The inertia of T in M3** is counted from a Householder tridiagonalization, and T is solved with by LU. A Jacobi eigenvalue method was tried first: the HS suite took 66 s instead of 6.7 s with L-BFGS and the direct QP.

### Results

Large problems (release, one run each with nothing else running, re-measured 2026-10-02; total time in seconds).

| problem | n | exact Hessian | exact, inertia | exact, inertia, direct | QPs solved directly |
|---|--:|--:|--:|--:|--:|
| control | 10,001 | 0.53 | 0.52 | 0.13 | 8 of 8 |
| rosenbrock | 5,000 | 0.91 | 0.93 | 0.04 | 15 of 15 |
| wells | 10,000 | 0.12 | 0.30 | 0.17 | 19 of 19 |
| circles | 10,000 | 0.95 | 0.98 | 0.07 | 6 of 6 |
| control | 100,001 | 26.8 | 22.1 | 1.3 | 10 of 10 |
| rosenbrock | 50,000 | 92.6 | 100.1 | 0.48 | 15 of 15 |
| wells | 100,000 | 1.1 | not finished in 10 min | 1.9 | 14 of 14 |
| circles | 100,000 | – | – | 1.1 | 8 of 8 |
| control | 1,000,001 | – | – | 13.2 | 11 of 11 |
| rosenbrock | 500,000 | – | – | 5.8 | 17 of 17 |
| wells | 1,000,000 | – | – | 18.8 | 20 of 20 |
| circles | 1,000,000 | – | – | 10.6 | 9 of 9 |

- **M1's decision point is passed:** the QP was skipped in every major iteration, and the time dropped by a factor of 2 to 200.
- **L-BFGS with the direct QP** (M3): with 10 pairs, 5.2 / 0.95 / 0.64 / 0.12 s at the smallest sizes (active-set QP, same memory: 17.8 / 1.5 / 1.8 / 0.95 s), and 122 / 11 / 12 / 2.2 s at the middle sizes. With 100 pairs: about 18 / 3.1 / 1.2 s on the first three at the smallest sizes (active-set: 27.8 / 1.9 / 1.9 s), so it only pays with a short memory. The cost is two solves per pair for each factorization, plus work of order n × pairs². The active-set reference for L-BFGS on `control` at n = 40,001 was stopped after 30 minutes without finishing.
- **HS suite** (solved / local / failed, `fc`): the defaults are unchanged (280/25/0, 9,173). Exact, inertia, direct: 273/30/2, 9,515 (without direct: 274/29/2, 9,581). L-BFGS, direct: 279/26/0, 8,920. SR1 with inertia control: 275/27/3, 11,714 (without: 244/31/30, 47,454); with the direct QP as well: 267/33/5. Trust region: SR1 with inertia 258/39/8 (without: 265/27/13); exact, inertia, direct 262/36/7.
- **M2 showed no benefit at first** (see "Settled since" for the problem that does show one). HS defaults with `direct_least_squares`: 279/25/1, 9,498. The new failure is TP106, where the direct correction is more accurate than `LSQR`'s (violation 7e-2 instead of 141 after it) and the iterates then crawl to the iteration limit. The large benchmarks never take a restoration step or a correction. It stays off by default.
- **Cost on small problems:** the 305 HS problems take 1.5 s with the exact Hessian and inertia control, 7.5 s with L-BFGS and the direct QP, and 13 s with SR1 and inertia control (1.0 s for the defaults).

### Settled since (2026-10-02)

- **CI.** The "Run tests with MUMPS" step passes on Linux (GitHub Actions, commit b7689a3), so the build flags and the library name work there too.
- **The ordering (G3).** MUMPS's automatic choice (`ICNTL(7)=7`) is kept, with no option: it was the best or within noise everywhere. Solver alone on a 3-D grid of order 125,000: automatic 4.2 s, METIS 4.2, PORD 5.1, SCOTCH 5.4, AMF 5.9, AMD 9.5, QAMD 9.6. `control` at n = 100,001 (whole solve, direct QP): automatic 1.38 s, AMF 1.35, QAMD 1.36, AMD 1.42, METIS 1.66, SCOTCH 1.73, PORD 1.78. `wells` at n = 100,000: 1.82 to 2.07 s.
- **Inertia control without the direct QP on `wells`.** It is the active-set QP, not the factorizations: with the smallest shift that makes the face convex, the QPs take more active-set iterations (160 instead of 32 at n = 20,000), and the solve takes 4.7 s instead of 2.6 s there, 69 s instead of 10.6 s at n = 40,000. The active-set solver's time on this problem is erratic anyway (the plain exact Hessian takes 10.6 s at n = 40,000 and 1.2 s at n = 100,000). With the direct QP it is 0.34 s and 0.55 s. So the advice stands: on large problems, use `inertia_control` with `direct_qp`.
- **Tests of the direct method's special paths** (`test_direct`): a nonconvex face without and with inertia control, a singular face that is regularized, the limit on the changes, and an infeasible QP.

- **The automatic L-BFGS memory with `direct_qp` is 10 pairs** (`lbfgs_memory = 0`; an explicit value is used as given). HS suite with `--direct`: 279/26/0 with 9,555 `fc` in 1.8 s (with 100 pairs: 8,920 `fc` in 7.5 s). SR1 with inertia control and the direct QP: 270/31/4 (267/33/5 with 100 pairs). SR1 with inertia control alone keeps the usual automatic memory (10 pairs gave 277/26/2 with 14,124 `fc`, 100 pairs 275/27/3 with 11,714: no clear winner).
- **Coverage.** CI now runs the default build's tests, and then the MUMPS build's with coverage (`coverage.sh --mumps`), so the report covers the code that uses MUMPS: 90.5% of the lines overall, and 90% or more of each new module.

- **A large problem for M2** (`circles` in `benchmark_large`, and in `test_direct` at n = 2000): a chain of `n-1` coupled circle constraints with the Maratos example's objective, whose steps need second-order corrections. Each constraint shares a variable with the next, so the corrections' least-squares problems are ill-conditioned and `LSQR` needs many iterations. With the direct QP, the solve takes 0.59 s with `LSQR` and 0.07 s with `direct_least_squares` at n = 10,000, and 88.7 s and 1.1 s at n = 100,000. A second problem, `hyperbolas` (only with `--problem=hyperbolas`), takes Gauss-Newton restoration steps: 9 of them cost 2.3 s with `LSQR` and 0.02 s directly at n = 10,000, but its total time (over 40 s) goes to the active-set QPs that find the linearization inconsistent, which the direct method can't replace.
- **Two bugs found with these problems, and fixed:**
  - *Workspace.* On a hanging-chain problem, whose Hessian of the Lagrangian is very indefinite, MUMPS needed many times its estimated workspace (delayed pivots). The factorization was retried with twice the allowance only 6 times, and the solve then ended as out of memory, on a matrix of order 602. It is now retried up to 20 times, as IPOPT does (`test_direct` has the regression test).
  - *Regularization of a singular face.* Near the solution of `circles` at n = 100,000, MUMPS reported the KKT matrix as singular (it isn't: half of the Hessian's diagonal is zero there). The regularized solve left the rows 7e-8 short, so the direct method gave up on a QP it had solved, and the active-set solver took minutes. A face whose regularized step only fails to satisfy its rows is now solved again with a much smaller regularization (1e-11 instead of 1e-8, relative). Using the smaller one from the start broke `wells`, where the rows really are dependent and the multipliers then blow up. In `test_qp_fuzz` the direct method now solves 276 of the 300 convex QPs instead of 244. HS with exact, inertia, direct: 273/30/2 with 9,515 `fc` (was 9,500).
- **The guide's reference timings** were re-measured with nothing else running (the table above).

### Still open

- The hanging chain itself: with the exact Hessian, with or without inertia control, the solver stops as "stalled" at a point that is not the solution (objective −65.2 for 200 links; L-BFGS finds −91.1). Not investigated. It is why that problem is only a regression test and not a benchmark.
- An inconsistent linearization is only found by solving an elastic QP with the active-set solver, which dominates `hyperbolas`. A cheaper test (or a direct elastic QP) would be needed for large problems that start far from feasible.
- The direct method isn't used for the feasibility QPs of a restoration phase (M5) or for an elastic re-solve, and the trust region makes no inertia test after its QPs.
- The Windows build with MUMPS has never been tried (the pixi tasks assume the library is called `dmumps_seq`).
