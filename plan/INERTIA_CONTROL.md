# Inertia control and a sparse LDLᵀ solver (MUMPS)

*Status: not started. This is open decision §8.3 of [ROADMAP.md](ROADMAP.md). Written 2026-09-28.*

This document is for future reference. It records what inertia control is, how sqpopt copes without it today, what an external sparse LDLᵀ solver such as MUMPS would add (for the exact Hessian and beyond), and what it would cost.

## 1. The problem: indefinite Hessians

The QP subproblem of each major iteration is

    min ½ pᵀHp + gᵀp   subject to the linearized constraints and bounds.

A step p is only a minimizer if H is positive definite on the null space of the active constraints: ZᵀHZ ≻ 0, where the columns of Z span that null space. Otherwise the QP is unbounded, or the solver stops at a saddle point.

| Hessian mode | Can H be indefinite? |
|---|---|
| L-BFGS (default) | No. Kept positive definite by construction (Powell damping). |
| SR1 | Yes. |
| Exact (`sqpopt_hessian_exact`, F7) | Yes. It is the user's Hessian of the Lagrangian. |

### Inertia

The standard test is the **inertia** of the KKT matrix, the count of its positive, negative and zero eigenvalues:

    K = [ H  Jᵀ ]    (n variables, m active constraints, J of full rank)
        [ J  0  ]

By Sylvester's law of inertia, In(K) = In(ZᵀHZ) + (m, m, 0). So ZᵀHZ ≻ 0 exactly when **In(K) = (n, m, 0)**.

A symmetric indefinite factorization K = LDLᵀ (Bunch–Kaufman pivoting, with 1×1 and 2×2 blocks in D) gives the inertia almost for free, from the signs of D's eigenvalues.

### Inertia control (as in IPOPT)

1. Factor K = LDLᵀ.
2. If K has more than m negative eigenvalues, replace H by H + δI and refactor. Increase δ until the inertia is right, and reuse the last successful δ as the next starting guess.
3. Solve with that factorization.

Each correction costs one refactorization, and the shift found is close to the smallest one that works.

## 2. What sqpopt does now: matrix-free, with indirect detection

sqpopt never forms or factors K. The QPs work with H only through products (and its diagonal), so negative curvature is found indirectly:

- the dense QP's reduced-Hessian Cholesky fails;
- CG on ZᵀHZ in the sparse QP meets a direction with pᵀHp ≤ 0;
- together these are reported as `qp_solver%negative_curvature`;
- or the QP fails, or returns a step that is not a descent direction for the merit function.

In exact mode, when any of these happens, H is shifted by δI and **the whole QP is re-solved**. δ starts at `shift_min`·max|H| and grows ×10 per retry, up to 15 retries. δ also grows where a quasi-Newton Hessian would be reset, and is divided by 3 after each good step. See ROADMAP.md, "F7 exact Hessian".

Weaknesses of this approach:

- **Incomplete detection.** CG only sees the curvature along the directions it explores. The dense QP only sees it in the null space of the current working set. Some negative curvature is found late or missed.
- **Expensive retries.** Each correction is a full QP solve, not one refactorization.
- **A coarse shift.** The ×10 steps overshoot the shift that's needed, which damps the Newton steps.
- **Measured cost.** On the HS suite, exact mode (with Hessians by finite differences) gives 268 solved, 33 local, 4 failed, against 276/29/0 for BFGS. It needs far fewer iterations on many problems, but more `fc` calls overall. Part of the gap is this correction scheme; the Hessians' finite-difference accuracy is another part.
- **SR1 has no correction at all.** It is much weaker than BFGS on the HS suite (239/35/31), plausibly because of indefinite updates.

## 3. Why the current dependencies can't provide it

| Solver | What it does | Inertia? |
|---|---|---|
| LUSOL | Sparse unsymmetric LU. Used by the sparse QP's basis method (F1, F13). | No |
| LSQR | Iterative least squares. The sparse QP's fallback null-space method. | No |
| CG | Iterative, on ZᵀHZ. | No (it only detects negative curvature along its own directions) |

The inertia needs a **symmetric indefinite** factorization, which none of them provides.

## 4. MUMPS

[MUMPS](https://mumps-solver.org) is a mature, open-source (CeCILL-C) sparse direct solver based on the multifrontal method.

- It has a symmetric indefinite LDLᵀ mode (`SYM=2`), and reports the **number of negative pivots** (`INFOG(12)`). That number is the inertia check.
- It handles large sparse problems well: fill-reducing orderings (AMD, METIS, SCOTCH), optional MPI/OpenMP, and reuse of the analysis phase when only the values change (which is our case, since the sparsity patterns are fixed).
- It has a native Fortran interface (`dmumps_struc`, with `JOB` = analyse / factorize / solve).
- It is widely used for exactly this purpose. IPOPT supports it as a linear solver for its inertia correction.

## 5. What it could be used for

### 5.1 Inertia control (the main reason)

- **Exact Hessian.** H is sparse, so it goes straight into K. This path would factor K once, check the inertia, shift minimally, refactor, and then take the step or feed the QP. It would replace the "shift ×10 and re-solve the QP" loop.
- **SR1.** The same correction would apply, but the limited-memory form needs the low-rank handling in §5.3.
- **L-BFGS** doesn't need inertia control.

### 5.2 Other sparse solves (in any Hessian mode)

These don't need the inertia, and LUSOL already covers several of them. MUMPS would compete with LUSOL here rather than add something new:

- **The sparse QP.** Since F1, products with Z use a LUSOL basis factorization. What remains iterative is CG on ZᵀHZ in the superbasic variables. With an explicit H, a direct factorization of the KKT system (or of ZᵀHZ, when the null space is small) could replace CG. That removes CG's accuracy issues and gives the inertia at the same time.
- **Least-squares systems in J.** Multiplier estimates, second-order corrections, and restoration/Gauss-Newton steps are minimum-norm or least-squares problems. They can be written as sparse symmetric systems such as `[I Jᵀ; J 0]`.

### 5.3 The catch with limited-memory Hessians

L-BFGS and SR1 are not sparse matrices. They are operators of the form σI + UCUᵀ, with a low-rank U. To use a factorization with them:

1. factor the sparse part, `[σI Jᵀ; J 0]`;
2. handle the low-rank part with the Sherman–Morrison–Woodbury formula, which is a small dense system in the rank.

IPOPT does this for its L-BFGS option. It works, but it is extra machinery. The inertia of the full matrix also has to be tracked through the low-rank correction.

## 6. Costs and risks

- **The build.** MUMPS needs BLAS/LAPACK and an ordering library (METIS or SCOTCH, or its built-in AMD). It also needs either MPI or its sequential stub library (`libseq`). None of this is fpm-native. The "build with fpm and nothing else" simplicity would be lost for anyone who enables it.
- **It must stay optional.** A preprocessor flag in the `.F90` files, with the matrix-free path always available, would work. Tests and CI would have to cover both builds; pixi (conda-forge has `mumps-seq`) could provide the optional one.
- **More code paths.** There would be a factorization-based QP, or KKT path, alongside the existing dense and matrix-free ones.
- **Dense rows or columns.** A few dense rows or columns in J or H cause a lot of fill in K.
- **The design principle.** "No dense `n×n`/`m×n` arrays" still holds, but the factors' fill can be large.

## 7. Alternatives

1. **Stay matrix-free** (the status quo), and improve it:
   - smarter shifts (IPOPT-style: reuse the last successful δ, grow ×8 at first and ×100 when there's no prior shift, and decrease gradually);
   - Lanczos-based curvature estimates on ZᵀHZ;
   - the same negative-curvature shift for SR1 (already noted in ROADMAP.md).
2. **HSL MA57 / MA27 / MA97.** These are the classic LDLᵀ solvers with inertia, and IPOPT's defaults. They are free for academic use, but their licences make redistribution harder.
3. **A Hessian-vector-product callback** (ROADMAP.md, F7 "not done"). It is useful when H is dense or expensive to form, but it is incompatible with a factorization. It favors staying matrix-free.

## 8. Suggested path, if pursued

1. First, cheap improvements to the matrix-free path (§7.1). Measure the exact and SR1 modes with `--hessian=exact|sr1` on the HS suite, to see how much of the gap to BFGS they close.
2. If a gap remains, prototype MUMPS as an optional dependency for **exact-Hessian inertia control only**: factor K with the current working set, check the inertia, shift, and pass the shifted H to the existing QPs. This isolates the benefit of accurate shifts from any change to the QP algorithms.
3. Only if that pays off, consider a factorization-based sparse QP (§5.2) and the low-rank extension for SR1 (§5.3).

Success measures:
- The HS suite with `--hessian=exact`: the solved/local/failed counts and `fc` against BFGS (currently 268/33/4 vs 276/29/0).
- `test_large_sparse` with analytic Hessians: the time and `fc`.
- The default (BFGS) results must be unchanged.
