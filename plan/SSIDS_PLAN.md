# Plan: a pure-Fortran SSIDS package (fpm)

A plan for a standalone fpm package, `ssids-fortran`: the sparse symmetric
indefinite direct solver SSIDS, extracted from
[SPRAL](https://github.com/ralna/spral) (STFC; BSD 3-clause), with its C++
CPU code translated to Fortran, its GPU code removed, and its external
dependencies (METIS, BLAS, LAPACK, hwloc) made optional, so that it builds
with fpm alone. Like [DAQP_PLAN.md](DAQP_PLAN.md), it is written to be read
on its own: section 10 says what sqpopt needs from it (roadmap item F25),
but the package must not depend on sqpopt.

The facts below about SPRAL's sources (files, line counts, dependencies,
options) are from a clone of SPRAL's `master` at commit `3d96a48`
(2026-08-14); the latest release tag is `v2025.09.18`. Pin a release tag
(stage 0) and correct this plan against it.

---

## 1. What SSIDS is, and why extract it

- **What it does.** Solves \( Ax = b \) for a sparse symmetric matrix,
  positive definite or indefinite, by a multifrontal \( LDL^T \)
  factorization: an analyse phase (ordering, elimination tree, supernode
  amalgamation, from the pattern only), a factorize phase (new values, same
  pattern), and a solve phase (any number of right-hand sides). For an
  indefinite matrix it uses 1×1 and 2×2 pivots with threshold pivoting,
  and reports the **inertia** (the number of negative pivots, and the rank),
  which is what an optimization solver's inertia control needs.
- **How.** Its CPU factorization uses *a posteriori threshold pivoting*
  (APP): a block of a front is factored without pivoting, then checked, and
  the columns that fail the threshold test are delayed. This keeps the work
  in level-3 BLAS and parallel OpenMP tasks, where traditional threshold
  pivoting (also present, as `ldlt_tpp`, for small fronts) is sequential
  per column. (J. D. Hogg, E. Ovtchinnikov, J. A. Scott, *A sparse symmetric
  indefinite direct solver for GPU architectures*, ACM TOMS 42(1), 2016;
  I. S. Duff, J. D. Hogg, F. Lopez, *A new sparse LDLᵀ solver using a
  posteriori threshold pivoting*, SIAM J. Sci. Comput. 42(2), 2020: check
  the references.)
- **Why a pure-Fortran package.** It is one of IPOPT's linear solvers, the
  best open-source (BSD) sparse indefinite solver with pivoting and inertia
  besides MUMPS, and mostly Fortran already. But SPRAL is built with Meson,
  needs METIS, BLAS, and LAPACK installed, and its CPU kernels are C++
  (with AVX intrinsics and OpenMP tasks). An fpm package that needs nothing
  else would give Fortran projects (sqpopt first) a pivoting \( LDL^T \)
  that is always available, as qdldl-fortran is for the non-pivoting case.
- **What it isn't.** A GPU solver (the GPU code is dropped), an
  unsymmetric solver, or as fast as SPRAL on large problems without an
  optimized BLAS and METIS (both stay available as options; section 7).

---

## 2. Inventory of SPRAL's sources (what SSIDS uses)

Line counts include comments.

| Part | Files | Lines | Language | In the package |
|---|---|---|---|---|
| SSIDS driver and data | `src/ssids/{ssids,anal,fkeep,akeep,datatypes,inform,contrib,contrib_free,subtree,profile_iface}.f90/.F90` | 4,300 | Fortran | **keep** (drop profiling) |
| CPU numeric subtree | `src/ssids/cpu/{cpu_iface,subtree}.f90` | 600 | Fortran | keep (adapt to the translated kernels) |
| CPU numeric subtree | `src/ssids/cpu/*.cxx, *.hxx` (NumericSubtree, SymbolicSubtree, SmallLeaf*, ThreadStats, Workspace, factor.hxx) | 2,200 | C++ | **translate** |
| CPU allocators | `BuddyAllocator.hxx`, `AppendAlloc.hxx`, `BlockPool.hxx`, `SimpleAlignedAlloc.hxx` | 700 | C++ | **replace** (Fortran allocatable pools) |
| CPU kernels | `src/ssids/cpu/kernels/` (ldlt_app 2,591, assemble 445, block_ldlt 415, ldlt_tpp 317, SimdVec 273, cholesky 215, verify 183, common 143, calc_ld 120, ldlt_nopiv 114, wrappers 92, ...) | 5,200 | C++ | **translate** (the core of the work) |
| Analyse support | `core_analyse.f90` (1,100), `matrix_util.f90` (3,304; only a part is used) | 4,400 | Fortran | keep (trim `matrix_util`) |
| Scaling and matching | `scaling.f90` (1,698: Hungarian/MC64-like, auction, MC77-like equilibration), `match_order.f90` (632; needs METIS) | 2,300 | Fortran | keep |
| Orderings | `metis5_wrapper.F90` (466), `metis4_wrapper.F90` (246) | 700 | Fortran + C METIS | **optional** (section 7) |
| BLAS/LAPACK | `blas_iface.f90`, `lapack_iface.f90` (interfaces); calls of `dgemm`, `dtrsm`, `dsyrk`, `dtrsv`, `dgemv`, `dpotrf`, `dsytrf` | 200 | Fortran | interfaces kept; implementations vendored, external optional |
| Machine topology | `hw_topology/` (hwloc wrapper, a guessing fallback) | 340 | Fortran + C++ | **replace** (OpenMP only) |
| OpenMP helpers | `omp.cxx/.hxx`, `compat.cxx` | 130 | C++ | replace |
| GPU | `src/ssids/gpu/`, `src/cuda/` | 14,400 | CUDA, C++, Fortran | **drop** |
| Debug dumps | `rutherford_boeing.f90` (1,022; `options%rb_dump`) | 1,000 | Fortran | optional (tests only) |
| Tests | `tests/ssids/ssids.f90`, `tests/ssids/kernels/*.cxx` | 4,700 | Fortran, C++ | port (translate the kernel tests) |

So: about **11,000 lines of Fortran to keep** (lightly modernized), about
**8,000 lines of C++ to translate** (5,200 of them kernels), and about
**15,000 lines to drop** (GPU, profiling).

**Features of the C++ that need a Fortran equivalent:**

- templates, instantiated for `double` (and some on a `bool` for the
  pivoting variant): write the double-precision versions (and real kinds
  by `wp`, section 3);
- OpenMP tasks: 25 `task`, 7 `taskgroup`, `taskwait`, `taskyield`, 13
  `cancel` / 14 `cancellation point`, 46 `atomic`, and `depend` clauses (38
  in `ldlt_app`, 9 in `cholesky`, 4 in `NumericSubtree`). Fortran OpenMP
  has all of these (4.0 for `depend` and `cancel`, 4.5 for `taskloop` if
  wanted). SPRAL requires `OMP_CANCELLATION=TRUE` in the environment (else
  error `SSIDS_ERROR_OMP_CANCELLATION`, -53);
- AVX/AVX2 intrinsics (`SimdVec.hxx`, used in `block_ldlt`, `calc_ld`,
  `common`), with a scalar fallback already in the code: use plain loops
  with `!$omp simd`;
- aligned allocation and custom allocators (buddy allocator for the
  factors, append-only allocator, block pool): Fortran allocatables with a
  simple pool per subtree, unless measurements say otherwise;
- C++ exceptions (32 `throw`/`catch`, mostly allocation failures): status
  flags (`stat=` on every allocation);
- classes with member functions (NumericSubtree, SymbolicSubtree,
  Workspace): derived types with type-bound procedures; the Fortran side
  now holds them as `type(c_ptr)` (30 `bind(C)` interfaces), which go away.

---

## 3. Package layout

```
ssids-fortran/
  fpm.toml
  LICENSE            BSD 3-clause, with STFC's copyright notice (section 9)
  NOTICE             credits SPRAL and the papers
  README.md
  ford.md
  upstream/spral/    git submodule: SPRAL, pinned to a release tag
  src/
    ssids_kinds.F90          real and integer kinds (preprocessor-selected)
    ssids_blas.F90           the BLAS routines used, in Fortran (reference
                             implementations), or external BLAS (-DSSIDS_EXTERNAL_BLAS)
    ssids_lapack.F90         dpotrf and dsytrf (and their callees), likewise
    ssids_ordering.F90       AMD (default), user-supplied, METIS (-DHAS_METIS)
    ssids_scaling.f90        from scaling.f90 (Hungarian, auction, equilibration)
    ssids_match_order.F90    from match_order.f90 (only with METIS)
    ssids_core_analyse.f90   from core_analyse.f90
    ssids_matrix_util.f90    the used part of matrix_util.f90
    ssids_datatypes.f90 ... ssids_anal.F90, ssids_fkeep.F90, ...  the SSIDS modules
    ssids_cpu_kernels_*.f90  the translated kernels (ldlt_app, ldlt_tpp, ldlt_nopiv,
                             cholesky, block_ldlt, calc_ld, assemble, wrappers)
    ssids_cpu_subtree*.f90   the translated numeric and symbolic subtrees
    ssids_module.f90         the public API (SPRAL's, plus an object interface, section 6)
  test/
    test_kernels_*.f90       ports of tests/ssids/kernels/*.cxx
    test_ssids.f90           port of tests/ssids/ssids.f90
    test_inertia.f90         matrices of known inertia, singular ones included
    test_random.f90          random KKT-like matrices, residuals and inertia
    test_kinds.f90           run in each real kind by CI
  compare/                   NOT part of the package build (section 5)
    fpm.toml                 depends on the package by path, links upstream SPRAL
    compare.f90
  tools/
    build_upstream.sh        builds upstream/spral with Meson (METIS, OpenBLAS from conda-forge)
  example/
    example_kkt.f90          a KKT matrix: analyse once, factor twice, inertia, solve
```

Conventions (as in qdldl-fortran and daqp-fortran):

- `implicit none` everywhere; reals `real(wp)`, literals `1.0_wp`; integer
  kinds for the pattern (`int32`, and `int64` for pointers, as SPRAL's
  `ptr64` interfaces).
- Kind selection: `-DREAL32`, `-DREAL64` (default), `-DREAL128`. SPRAL is
  double only: the kernels' tolerances (pivot threshold `u`, `small`) must
  be reviewed for the other kinds.
- **No local variable initialized in its declaration.**
- **Nothing sized by the problem on the stack:** fronts and workspaces
  allocatable, allocated with `stat=`, reported as out of memory.
- FORD comments: `!>` module headers, `!!` for every procedure, type,
  component, and dummy argument; mark ported files as changed (section 9).
- No C, no C++, no `bind(C)` in the default build.

---

## 4. Stages

### Stage 0: upstream as a submodule, built separately

- `git submodule add https://github.com/ralna/spral upstream/spral`, pinned
  to a release tag (`v2025.09.18` or later).
- `tools/build_upstream.sh`: Meson build of SPRAL (CPU only,
  `-Dgpu=false`), with METIS, OpenBLAS, and hwloc from conda-forge (a pixi
  environment), into `build/upstream/`. Record compilers, flags, versions.
- Run upstream's tests (`OMP_CANCELLATION=TRUE`); they must pass before
  anything is compared against them.
- Correct section 2 against the pinned tag.

**Done when:** the upstream library and its tests build and pass from a
clean clone with one command.

### Stage 1: the Fortran parts under fpm, with the C++ kernels as they are

A working fpm package first, before any translation:

- Copy the Fortran sources of section 2 (without GPU, profiling, hwloc) into
  `src/`, and the C++ CPU sources too (fpm builds C and C++ sources).
- BLAS/LAPACK external for now (`--link-flag "-lopenblas"`), METIS external
  (`HAS_METIS`), hwloc removed (one NUMA region; `ignore_numa` is already
  the default).
- Port `tests/ssids/ssids.f90` to an fpm test.

**Done when:** `fpm test` passes, with results identical to upstream's (the
same code, built differently).

### Stage 2: translate the kernels, one at a time

In dependency order, each with its own test comparing the Fortran kernel
with the C++ one (still in the tree) on random fronts, bit for bit where
the arithmetic is the same (fused multiply-adds disabled in both), else to
round-off:

1. `wrappers` (the BLAS/LAPACK calls), `common`, `SimdVec` (as plain loops);
2. `cholesky` (with its task `depend` clauses);
3. `ldlt_nopiv`, `ldlt_tpp`;
4. `calc_ld`, `block_ldlt`;
5. `ldlt_app` (2,591 lines: the a posteriori pivoting, its block column
   structure, the tasks with dependencies, the cancellation on failure, and
   the backup and restore of failed blocks): the largest single piece;
6. `assemble` (the extend-add of child contributions).

Translate the C++ kernel tests (`tests/ssids/kernels/*.cxx`) alongside.

**Done when:** each kernel passes its comparison, and `fpm test` passes
with the Fortran kernels in place of the C++ ones.

### Stage 3: translate the subtrees, remove the C++

- `SymbolicSubtree`, `NumericSubtree`, `SmallLeafSymbolicSubtree`,
  `SmallLeafNumericSubtree`, `Workspace`, `ThreadStats`, `factor.hxx`:
  derived types and procedures; `cpu_iface.f90` and `subtree.f90` call them
  directly instead of through `type(c_ptr)`.
- Allocators: replace the buddy allocator, append allocator, and block
  pool by Fortran allocatable storage (one contiguous array per subtree for
  the factors, as the buddy allocator gives, if the measurements need it).
- `omp.cxx` and `compat.cxx`: Fortran OpenMP calls (`omp_lib`), under `!$`
  sentinels so that the package builds and runs **without OpenMP** too.
- Decide how to handle cancellation: keep `!$omp cancel` (and document
  `OMP_CANCELLATION=TRUE`), or replace it by a shared failure flag checked
  by the tasks, so that the environment variable is no longer needed.
- Remove the C++ sources and the `bind(C)` interfaces.

**Done when:** the package has no C or C++, and the comparison with stage
1's build (section 5) passes on every set.

### Stage 4: dependencies made optional

- **Ordering.** Default: AMD in Fortran (from qdldl-fortran's
  `qdldl_amd.f90`, itself a port of SuiteSparse's AMD, BSD 3-clause), and
  a user-supplied order (`options%ordering = 0`, as now). Optional: METIS
  through its C interface under `-DHAS_METIS` (SPRAL's wrapper, as it is).
  Note that `match_order` (ordering 2, matching-based) needs METIS. A
  Fortran port of METIS's nested dissection is a separate project, only if
  the measurements of section 5 show AMD's fill is a problem.
- **BLAS/LAPACK.** Vendor the reference implementations (Netlib, modified
  BSD) of the seven routines and their callees (`dsytrf` pulls in
  `dlasyf`, `dsytf2`, `idamax`, `dswap`, `dscal`, `dger`, `dsyr`, `ilaenv`,
  `lsame`, `xerbla`), modernized, in `ssids_blas`/`ssids_lapack`, with
  `-DSSIDS_EXTERNAL_BLAS` to call an optimized library instead. Consider
  a dependency on fortran-lang's stdlib, whose `stdlib_linalg_blas` and
  `stdlib_linalg_lapack` provide the same in pure Fortran with the same
  external switch; weigh its size against vendoring seven routines.
- **hwloc:** gone (stage 1).

**Done when:** `fpm build` with no flags and nothing installed gives a
working solver, and each option (`HAS_METIS`, external BLAS, OpenMP) builds
and passes the tests.

### Stage 5: the object interface, docs, kinds, release

- The interface of section 6, the README (what it is, the options, the
  inertia, `OMP_CANCELLATION` if kept, how to get speed: an optimized BLAS
  and METIS), FORD docs.
- The kinds (section 3), with the tolerances reviewed for single and
  quadruple precision.
- CI: gfortran (and ifx if convenient), with and without OpenMP, with and
  without external BLAS and METIS; a separate job that builds upstream and
  runs the comparison on a small set.
- Release `v0.1.0`.

---

## 5. Comparing the port with upstream

- `compare/` is a separate fpm project, depending on the package by path,
  linked with the upstream library of stage 0 (`bind(C)`-free: SPRAL's own
  Fortran API, from its `.mod` files, under different module names, or a
  thin wrapper program built against upstream and run as a separate
  process, which avoids module-name clashes).
- **Matrices:** random symmetric indefinite matrices with known inertia
  (built as \( Q D Q^T \), sparse), singular ones (rank deficiency, zero
  rows), KKT matrices \( \begin{bmatrix} H & J^T \\ J & 0 \end{bmatrix} \)
  (from sqpopt's `sparse_solvers` example and its HS and CUTEst runs), the
  2-D and 3-D grid matrices of `sparse_solvers`, and small matrices of the
  SuiteSparse collection.
- **Same ordering for the comparison:** give both the same user-supplied
  order (`ordering = 0`) so that differences come from the factorization
  only; then compare AMD (the port's default) with METIS (upstream's) for
  fill (`inform%num_factor`) and time.
- **Pass when:** the same inertia (`num_neg`, `matrix_rank`) and the same
  number of delayed pivots and 2×2 pivots; the scaled residual
  \( \lVert Ax - b \rVert / (\lVert A \rVert \lVert x \rVert + \lVert b \rVert) \)
  of each within a small multiple of the other's; bit-for-bit factors where
  the arithmetic is the same.
- **Speed:** analyse, factor, and solve times separately, at 1 and 4
  threads, with an optimized BLAS in both. Target: the port's factorization
  within 1.5× of upstream's with the same BLAS and ordering. Also record the
  vendored reference BLAS's times, which will be much slower on large
  fronts: that is the cost of "nothing installed", and the README must
  say so.
- `compare/RESULTS.md`: the tables, the environment, the upstream tag.

---

## 6. Public API (target)

SPRAL's procedural API, unchanged for its users (`ssids_analyse`,
`ssids_analyse_coord`, `ssids_factor`, `ssids_solve`, `ssids_free`,
`ssids_enquire_posdef`, `ssids_enquire_indef`, `ssids_alter`, with
`ssids_akeep`, `ssids_fkeep`, `ssids_options`, `ssids_inform`), under the
package's module name, plus an object interface:

```fortran
type(ssids_type) :: solver
call solver%analyse(n, ptr, row, istat [, order])       ! the pattern (CSC, lower triangle)
call solver%factor(val, istat, posdef=.false.)          ! new values, same pattern
call solver%inertia(n_negative, n_null)                 ! from inform%num_neg and matrix_rank
call solver%solve(x, istat)                             ! x: b on entry, the solution on exit (one or more columns)
call solver%destroy()
```

- Options as components with SPRAL's defaults (`u` (pivot threshold),
  `small`, `nemin`, `scaling`, `ordering`, `action` on singularity, ...).
- Status codes as named constants (SPRAL's flags, plus out of memory).
- Copyable after stage 3 (allocatable components only, no C pointers).

---

## 7. Pitfalls to watch

- **`ldlt_app`'s task graph.** Its correctness depends on the `depend`
  clauses and on cancellation when a block fails; a translation error
  shows up as a rare race, not a wrong answer every time. Test with many
  threads, many seeds, and repeated runs; compare with one thread.
- **Pivoting decisions are sensitive to rounding.** A different order of
  operations (e.g. fused multiply-adds) can delay a different column, which
  changes the factors but not the inertia or the accuracy. Compare inertia
  and residuals, and bit for bit only with matched arithmetic.
- **Performance without an optimized BLAS and METIS.** SSIDS's speed is in
  level-3 BLAS on the fronts and in METIS's orderings. With the vendored
  reference BLAS and AMD it will be correct but slower on large problems;
  that is a documented trade-off, with both options a flag away.
- **OpenMP in Fortran.** Task `depend` clauses on array sections, `cancel`
  semantics, and `atomic` updates differ in detail from C++; check against
  the OpenMP specification, and test with gfortran and ifx.
- **Allocators.** The buddy allocator keeps a subtree's factors contiguous
  and avoids many small allocations; a naive translation (one allocatable
  per node) may cost time or memory. Measure before optimizing.
- **`matrix_util` is large** (3,304 lines) and only partly used: extract
  what SSIDS calls, rather than carrying all of it.
- **64-bit pointers.** SPRAL has `ptr32` and `ptr64` versions of its
  interfaces (generic interfaces over both); keep both.
- **Licences of the vendored code.** Netlib BLAS/LAPACK (modified BSD) and
  SuiteSparse AMD (BSD 3-clause) are compatible; METIS (Apache-2.0) stays
  an external, optional library.

---

## 8. Validation matrix (what "done" means)

| Check | Problems | Pass when |
|---|---|---|
| Upstream tests | `tests/ssids` | pass (ported) |
| Kernels vs C++ | random fronts | bit for bit with matched arithmetic, else round-off |
| Port vs upstream | section 5 sets | same inertia, delays, 2×2 pivots; residuals within a small multiple |
| Inertia | known-inertia and singular matrices | exact `num_neg` and rank |
| Threads | 1 to 8 threads, repeated runs | identical inertia, no failures or hangs |
| Without OpenMP | all tests | pass |
| Options | `HAS_METIS`, external BLAS | build and pass |
| Kinds | REAL32, REAL64, REAL128 | tests pass (tolerances scaled) |
| Memory | huge `n` | out-of-memory status, no crash |
| Stack | the whole library | no automatic arrays or problem-sized temporaries |
| Speed | grid and KKT matrices | within 1.5× of upstream with the same BLAS and ordering |

---

## 9. Licensing

- SPRAL is BSD 3-clause, © STFC. The package is a derivative work: release
  it under BSD 3-clause too, keeping STFC's notice in `LICENSE` beside the
  port's, credit SPRAL and the papers in `NOTICE` and the README, and mark
  each ported file as changed.
- Vendored Netlib BLAS/LAPACK routines (modified BSD) and SuiteSparse AMD
  (through qdldl-fortran; BSD 3-clause) keep their notices.
- METIS (Apache-2.0) and external BLAS are linked, not included.
- BSD 3-clause is compatible with sqpopt's MIT licence. (Not legal advice.)

---

## 10. The first user: sqpopt (roadmap F25)

sqpopt's sparse factorizations go through one abstract backend
(`sqpopt_sparse_ldl_type`, in `src/sqpopt_sparse_ldl_module.f90`), with five
procedures: `start` (analyse a pattern once), `refactor` (factor new values,
and report what inertia it can tell: `n_negative`, `n_null`, and
`inertia_known`), `back_solve`, `multiply` (the product with the last
values, for iterative refinement), and `free`.
QDLDL (no pivoting, always available), MUMPS (optional), and the dense
solvers (opt-in, small problems) implement it. SSIDS would be a fourth
backend:

| sqpopt needs | SSIDS gives |
|---|---|
| analyse once, factor many times (`start`, `refactor`) | `ssids_analyse` once, `ssids_factor` per KKT matrix |
| the product with the matrix (`multiply`) | not in SSIDS: the backend keeps the values and multiplies itself, as the QDLDL backend does |
| the inertia, also of singular matrices | `inform%num_neg`, `inform%matrix_rank` (null pivots = `n - rank`), with `action = .true.` to continue on singularity |
| pivoting (for exact Hessians with zeros on the diagonal, where QDLDL over-shifts) | threshold 1×1/2×2 pivoting with a posteriori checks |
| any real kind, no external libraries | the package's kinds, AMD, vendored BLAS |
| threads | OpenMP (`options%factorization_threads` maps to the thread count) |
| out of memory reported | the package's status |

Steps (a roadmap item in sqpopt, not part of this package):

1. Add the package as an fpm dependency (pinned tag).
2. `sqpopt_ssids_ldl_module.f90`: a backend extending
   `sqpopt_sparse_ldl_type`; a new `options%linear_solver =
   sqpopt_linear_solver_ssids`; a case in `symmetric_solver_initialize`.
   Since the package is pure Fortran, it needs no preprocessor guard.
3. Run `test_kkt`, `test_direct`, `test_inertia` with it, and the HS
   protocol of CLAUDE.md with `--linear-solver=ssids` (`--hessian=exact
   --inertia`, `--hessian=exact --inertia --direct`, `--hessian=sr1
   --inertia`, `--direct`, `--direct-ls`), compared with MUMPS and the dense
   solver; and `benchmark_large` and the `sparse_solvers` example (grid and
   KKT matrices) for the times.
4. Decide from those whether SSIDS should be the recommended solver for the
   exact Hessian with inertia control (in place of "MUMPS, or the dense
   solver for small problems"), and update the guide's Linear solver
   section and the Choosing settings page.

---

## References

- SPRAL: <https://github.com/ralna/spral>; documentation of SSIDS:
  <https://ralna.github.io/spral/> (check the URL).
- J. D. Hogg, E. Ovtchinnikov, J. A. Scott, *A sparse symmetric indefinite
  direct solver for GPU architectures*, ACM Transactions on Mathematical
  Software 42(1), 2016 (check).
- I. S. Duff, J. D. Hogg, F. Lopez, *A new sparse LDLᵀ solver using a
  posteriori threshold pivoting*, SIAM Journal on Scientific Computing
  42(2), 2020 (check).
- P. R. Amestoy, T. A. Davis, I. S. Duff, *An approximate minimum degree
  ordering algorithm*, SIAM J. Matrix Anal. Appl. 17(4), 1996 (AMD).
