# Comparison with HiOp

[HiOp](https://github.com/llnl/hiop) (Cosmin Petra, Nai-Yuan Chiang, Jingyi
Wang and others, LLNL, BSD-3 license) is a C++ solver for large nonlinear
programs on HPC machines: MPI-distributed memory and GPUs (CUDA, HIP,
through RAJA and Umpire). This note compares its design with sqpopt's and
lists ideas worth borrowing.

Sources: the HiOp source code (commit `fb37b2c`, October 2026), its user
manual (`doc/hiop_usermanual.pdf`), and the papers it cites: Petra, *A
memory-distributed quasi-Newton solver for nonlinear programming problems
with a small number of general constraints* (JPDC 133, 2019); Petra et al.,
*On the implementation of a quasi-Newton interior-point method for
PDE-constrained optimization using finite element discretizations* (OMS 38,
2023); Chiang and Zavala, *An inertia-free filter line-search algorithm for
large-scale nonlinear programming* (COAP 64, 2016); Wang and Petra, *A
sequential quadratic programming algorithm for nonsmooth problems with
upper-C² objective* (SIOPT 33, 2023). HiOp was not built or benchmarked for
this comparison.

## Summary

HiOp is an interior-point method, not an SQP method: its main algorithm is
the Wächter–Biegler filter line-search IPM (as in IPOPT), and its design is
driven by where the linear algebra runs (several MPI ranks, a GPU), not by
small-to-medium problems on one core. Most of its machinery (condensed and
GPU-friendly KKT forms, the barrier parameter's rules, bound relaxation) has
no counterpart in an active-set SQP and isn't worth porting. Its separate
PriDec solver (a decomposition method for two-stage stochastic problems) and
`hiopbbpy` (Bayesian optimization in Python) are different products.

What is worth borrowing is small and is all about the factorization-based
options: the **inertia-free curvature test** (for QDLDL, whose inertia
sqpopt can't trust with the exact Hessian), a **"safe mode" that switches to
a more robust linear solver** after a failure, and the
**least-squares multiplier estimate**. Two cheaper experiments are the
**L-BFGS initial scaling** and, further off, **weighted inner products**
for discretized function spaces.

## How HiOp is organized

HiOp has three problem interfaces, each with its own KKT linear algebra:

| Interface | For | Hessian | KKT solve |
|---|---|---|---|
| `hiopInterfaceDenseConstraints` | many variables, few (dense) constraints, distributed over MPI ranks | L-BFGS only (compact form, diagonal plus low rank) | Sherman–Morrison–Woodbury: only `m x m` and `2l x 2l` dense systems |
| `hiopInterfaceSparse` | general sparse problems | exact | sparse LDLᵀ (MA57, PARDISO, STRUMPACK, cuSOLVER/ReSolve, Ginkgo), or condensed or normal-equation forms solved by Cholesky on a GPU (HyKKT) |
| `hiopInterfaceMDS` | mixed dense-sparse: a sparse block plus a few dense "coupling" variables | exact | Schur complement onto the dense block, dense LDLᵀ (LAPACK or MAGMA on a GPU) |

Around them:

- **Step acceptance:** the Wächter–Biegler filter (switching condition,
  Armijo, second-order correction, `theta_min`/`theta_max` from the initial
  violation), and a feasibility restoration phase that solves a separate
  NLP (`hiopFRProb*`) with HiOp itself.
- **Factorization acceptance** (`fact_acceptor`): IPOPT-style inertia
  correction (`hiopPDPerturbation`: primal shift δw, dual regularization
  δc), or the inertia-free curvature test of Chiang and Zavala (see idea 1).
- **"Safe mode":** the fast condensed KKT form (Cholesky, no inertia) is
  used until it fails; HiOp then switches for good to the stabler
  `XDYcYd` form with inertia correction, and changes the barrier rule.
- **Problem transforms** (`hiopNlpTransforms`): fixed variables removed or
  their bounds relaxed; gradient-based scaling of the objective and
  constraints (as in IPOPT; sqpopt's `options%scaling` is the same idea);
  bound relaxation and an "elastic mode" that tightens it as the barrier
  parameter falls.
- **Weighted spaces** (`VectorSpace`): for a problem discretized by finite
  elements, norms, inner products, the L-BFGS update and the termination
  test use a lumped mass (or mass-plus-stiffness) matrix instead of the
  Euclidean ones, so the iteration counts don't grow as the mesh is refined.
- **Other:** a finite-difference derivative checker; checkpoints of the
  whole solver state to a file every N iterations (through Axom Sidre), to
  restart a long run; warm start; a C and a Fortran interface.

## Comparison

- **Method.** HiOp's barrier method is at its best on very large problems
  with many inequalities, where one factorization per iteration (of a
  matrix whose pattern never changes) beats an active-set QP. sqpopt's
  active-set SQP with quasi-Newton Hessians needs fewer function
  evaluations on small and medium problems, and warm-starts well; it
  matches the use case of expensive functions (trajectories, simulations)
  that HiOp's HPC focus doesn't target.
- **Quasi-Newton.** Both keep L-BFGS in the compact form of Byrd, Nocedal
  and Schnabel and never form an `n x n` matrix. HiOp skips an update with
  `sᵀy ≤ √ε ‖s‖‖y‖`; sqpopt damps it (Powell) by default. HiOp has no SR1.
  HiOp's quasi-Newton mode exists only for the dense-constraints interface;
  sqpopt's works with any sparsity.
- **Linear algebra.** sqpopt's layering (an abstract LDLᵀ backend, QDLDL by
  default, MUMPS and LAPACK optional) is HiOp's on a smaller scale. HiOp
  has many more backends, but each is an external library; sqpopt's
  default needs none.
- **Robustness features both have:** a filter line search with SOC, a
  feasibility restoration phase, gradient-based scaling, acceptable-level
  termination, a maximum step per variable (HiOp's `moving_lim_abs`/`rel`;
  sqpopt's `problem%set_max_step`), and fixed variables (sqpopt's QPs fix
  them directly, so neither of HiOp's transforms is needed).
- **What sqpopt has and HiOp doesn't:** SR1, trust region and funnel
  options, the augmented Lagrangian merit, elastic QPs, the least-squares
  interface, and a dependency-free default build.

## Ideas worth borrowing

1. **The inertia-free curvature test** (`fact_acceptor = inertia_free`,
   `neg_curv_test_fact`, default 1e-11). Instead of requiring the right
   inertia, accept a factorization if the step it gives has enough
   curvature, `dᵀ(H + δI)d ≥ κ dᵀd` (Chiang and Zavala give the exact form
   with the constraints' part), and raise the shift only if it fails. This
   is aimed at solvers that can't report inertia, which is QDLDL's case in
   sqpopt: with the exact Hessian a zero pivot that isn't a zero eigenvalue
   is counted as negative curvature, and `kkt_negatives` has to guess. The
   gap it could close is measurable: `--hessian=exact --inertia` is 273
   solved, 3 failed and 13,654 `fc` with QDLDL, against 274, 2 and 9,400
   with MUMPS. The test would apply where the factorization's step is used
   (`inertia_correct`, the direct QP step); the active-set QPs keep their
   own negative-curvature handling. Caveat: the test accepts a direction,
   not a convex QP, so it fits the direct step better than the shift that
   the iterative QP solvers then use. Try it as a third option of the
   inertia control (`inertia`, `inertia_free`), with QDLDL first.

2. **A safe mode for the linear solver.** HiOp runs the fast, less stable
   KKT form until it fails, then switches to the robust one for the rest
   of the solve. In sqpopt that would be: if QDLDL's factorization is
   singular, or the inertia correction reaches `shift_max`, refactor with
   a pivoting backend (MUMPS if built, else the dense LDLᵀ when the order
   is at most 2,000) and keep it for the rest of the solve. This belongs
   in `sqpopt_symmetric_solver_module` (nothing above it changes), behind
   an option (e.g. `linear_solver_fallback`), and is measured with the
   same `--hessian=exact --inertia` runs as idea 1. Today a failed
   factorization falls back on the matrix-free methods instead.

3. **Least-squares multiplier estimates** (`duals_init = lsq`,
   `duals_update_type = lsq`, capped by `duals_lsq_ini_max`, default 1e3). The
   multipliers that minimize `‖∇f + Jᵀλ‖`, discarded if their norm is too
   large. For sqpopt this matters only where there is no QP's λ to use:
   the exact Hessian's first iteration (which uses `lambda0`, zero by
   default), and the exit from the restoration phase. The solve is the
   transpose of what `sqpopt_least_squares_module` already does (the same
   augmented matrix, with a right-hand side of `[-∇f; 0]`), or LSQR
   without it. Check with `--hessian=exact` and `--hessian=sr1 --inertia`.

4. **The L-BFGS initial scaling** (`sigma_update_strategy`). HiOp offers
   five choices of `B₀ = σI`: `sᵀy/sᵀs` (its default), `yᵀy/sᵀy`
   (sqpopt's, the Shanno–Phua choice), `‖y‖/‖s‖`, their average, and a
   constant, all clamped to `[1e-8, 1e8]`. A one-line change in
   `hessian_update_bfgs` to try on the HS suite; adopt only if it
   improves the baseline, and record the numbers either way.

5. **Weighted inner products** (`WeightedSpaceType`, `VectorSpace`). For
   problems discretized on a non-uniform mesh (collocation with segments of
   different lengths, for example), using the lumped mass matrix `M` in the
   L-BFGS update (`sᵀMs`, `yᵀM⁻¹y`), in `B₀ = σM`, and in the stationarity
   norm makes the method independent of the mesh. A diagonal `M` is a
   generalization of the variable scaling idea from the OPTGRA comparison
   (`problem%set_variable_scaling`), and should be designed with it, not
   separately. Larger, and only worth it if a user has such problems.

6. **Checkpoints** (`checkpoint_save`, `checkpoint_save_every_N_iter`,
   `checkpoint_load_on_start`). Saving x, λ, the L-BFGS pairs, the penalty
   and the filter, to restart a long run after a crash. Low priority: it
   needs a file format, and it collides with the rule that no state
   carries over between solves unless the user passes it in explicitly
   (as the filterSD idea of starting from the previous solve's L-BFGS
   pairs would). If done, make it an explicit argument of `solve`, not
   hidden state.

Not worth borrowing:

- The barrier method's machinery: the condensed, normal-equation and
  HyKKT KKT forms, the rules for μ, bound relaxation and the elastic mode
  that tightens it, the fraction-to-the-boundary rule.
- GPU and MPI execution (RAJA, Umpire), and the MDS interface (a Schur
  complement onto a dense block: sqpopt's sparse backends handle such a
  block directly).
- PriDec and `hiopbbpy`: separate algorithms for other problem classes.
- The derivative checker: F8, dropped at the user's request (sqpopt's
  level-2 diagnostics already check the derivatives without evaluations).
- Randomized regularization (`regularization_method`, off by default in
  HiOp).

## Suggested order

Ideas 1 and 2 first, together: they address the same measured gap
(QDLDL against MUMPS with the exact Hessian) and are tested with the same
runs. Then 4 (an afternoon), then 3. Ideas 5 and 6 only if a user asks.

## Status

Ideas 1 and 2 done (2026-11-08); see the HS results below. Ideas 3 to 6
not started. Listed in `plan/ROADMAP.md` §5 as "Ideas from HiOp".

Measured on the HS suite (release build, solved / local / failed / `fc`):

| configuration | QDLDL alone | QDLDL with the fallback | curvature test, where the inertia is unknown | curvature test, always |
|---|---|---|---|---|
| exact, inertia control | 273 / 28 / 4 / 13,563 | 273 / 30 / 2 / 9,431 | 272 / 30 / 3 / 9,508 | 277 / 27 / 1 / 12,599 |
| exact, inertia control, direct QP | 275 / 28 / 2 / 13,812 | 277 / 27 / 1 / 9,467 | 274 / 30 / 1 / 9,444 | 271 / 32 / 2 / 14,073 |
| SR1, inertia control | 275 / 29 / 1 / 11,546 | unchanged | unchanged | 268 / 31 / 6 / 15,130 |

- **Idea 2, the fallback**, is `options%linear_solver_fallback`, on by
  default: it gives the dense solver's (and MUMPS's) results, and changes
  nothing where QDLDL can tell the inertia (the quasi-Newton Hessians). In a
  build without MUMPS it is a dense solver, so it is limited to an order of
  500.
- **Idea 1, the curvature test**, is `options%inertia_test`, off by
  default (`sqpopt_inertia_test_count`): always testing the curvature
  solved the most problems with the exact Hessian (the same with every
  linear solver, since it doesn't use the inertia), but with 34% more
  evaluations, and it was worse with the direct QP method and with SR1.
  Testing it only where the inertia is unknown does nothing once the
  fallback is on.
