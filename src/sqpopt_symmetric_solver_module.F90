!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  A sparse direct solver for symmetric indefinite linear systems
!  \( A x = b \), by an \( LDL^T \) factorization, with one of two sparse
!  solvers (`options%linear_solver`):
!
!  * [MUMPS](https://mumps-solver.org): multifrontal, with pivoting, and
!    threads. This is the library's only interface to MUMPS, and the only
!    source file that depends on the `HAS_MUMPS` preprocessor directive:
!    without it, `sqpopt_has_mumps` is false and MUMPS can't be chosen. With
!    it, the library must be linked with the sequential, double precision
!    MUMPS library, and compiled in double precision (the default kind, see
!    [[sqpopt_kinds]]).
!  * [QDLDL](https://github.com/jacobwilliams/qdldl-fortran) (a Fortran
!    port of OSQP's solver, an fpm dependency, so always available, in any
!    real kind): an up-looking factorization in a fill-reducing order (AMD),
!    *without pivoting*, single-threaded, with very little overhead per
!    call. It is exact for quasi-definite matrices, and for the KKT matrices
!    \( \begin{bmatrix} H & J^T \\ J & 0 \end{bmatrix} \) with \( H \)
!    positive definite once each constraint's row comes after its variables
!    (the `signs` of [[symmetric_solver_initialize]], see
!    [[delay_negative_rows]]): the least-squares systems, and the matrices
!    of the quasi-Newton Hessians. For an indefinite \( H \) its inertia is
!    still exact as long as no pivot is zero, but a zero pivot that isn't a
!    zero eigenvalue (e.g. from a Hessian with zeros on its diagonal, which
!    would need a 2x2 pivot) leaves the inertia unknown: it is then reported
!    as a negative eigenvalue (see [[symmetric_solver_factor]]), so that the
!    inertia control shifts the Hessian, by more than MUMPS would need.
!
!  QDLDL is the default, in every build. Measured (release; see the user
!  guide's "Sparse solver" section): on the banded and chained problems of
!  `example/benchmark_large.f90` (5,000 to 1,000,000 variables) both take
!  the same iterations, and QDLDL's factorizations take 7 to 30 times less
!  time; on the Hock-Schittkowski problems the quasi-Newton Hessians give
!  the same results with both (QDLDL's faster), but the exact Hessian with
!  inertia control needs about 45% more evaluations with QDLDL. On grids
!  (`example/sparse_solvers.f90`), MUMPS refactors 2 to 10 times faster in
!  2-D and 9 to 34 times in 3-D, where the factors have large dense blocks
!  (and its threads help), though its first factorization, with its
!  analysis, can cost much more (72 s against 4.5 s for a 2-D KKT matrix of
!  order 735,000). So MUMPS is the better choice for the exact Hessian with
!  inertia control, and for problems coupled in two or three dimensions that
!  refactor many times.
!
!  How it is used:
!
!  1. [[symmetric_solver_initialize]] with the sparsity pattern (the entries
!     of one triangle, in coordinate form; entries given more than once are
!     added together). QDLDL analyses (orders) it here, MUMPS in the first
!     factorization; either way only once, so a matrix whose values change
!     but whose pattern doesn't is cheap to refactor.
!  2. [[symmetric_solver_factor]] with the values of those entries. The
!     factorization gives the *inertia* of the matrix: `n_negative`
!     negative and `n_null` zero eigenvalues.
!  3. [[symmetric_solver_solve]], any number of times, for right-hand sides.
!     Each solution is improved by iterative refinement.
!  4. [[symmetric_solver_destroy]] to free it. The type holds pointers (for
!     MUMPS), so it must not be copied.
!
!  **Threads.** MUMPS and the BLAS under it can use OpenMP threads, if they
!  were built with OpenMP (conda-forge's `mumps-seq` is). The number of
!  threads is the `threads` argument of [[symmetric_solver_initialize]]:
!  `1` by default, a larger number for that many threads, or `0` to leave it
!  to the OpenMP environment (`OMP_NUM_THREADS`, or every core). The default
!  is one thread because on small matrices the threads cost far more time
!  than they save (the 305 Hock-Schittkowski problems took 24.5 s with 10
!  threads, and 1.4 s with one). Whether more threads help a large matrix
!  depends on its factors: the threads work inside the dense blocks of the
!  factorization (its frontal matrices), so they need those to be large.
!  Measured (analysis, two factorizations, and a solve; 1, 2, 4, and 8
!  threads): a 3-D grid matrix of order 216,000 took 11.5, 7.0, 4.9, and
!  4.8 s; a 2-D grid matrix of order 490,000 took 2.0 s with any number;
!  and the banded KKT matrices of `example/benchmark_large.f90` gained
!  nothing either. MUMPS sets the number of threads when it is called and
!  restores it when it returns, so the rest of the program is not affected.
!  QDLDL is single-threaded.
!
!  A singular matrix is not an error: its null pivots are counted
!  (`n_null`), and a solve then returns one of the solutions of a consistent
!  system.

    module sqpopt_symmetric_solver_module

    use, intrinsic :: iso_fortran_env, only: int64
    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_all_finite
    use qdldl_module, only: qdldl_type, qdldl_wp, qdldl_success, qdldl_error_out_of_memory, qdldl_order_amd, &
                            qdldl_error_zero_pivot

    implicit none

    private

#ifdef HAS_MUMPS
#if defined(REAL32) || defined(REAL128)
#error "HAS_MUMPS needs the default (double precision) real kind: don't define REAL32 or REAL128 with it"
#endif
    include 'dmumps_struc.h'

    interface
        subroutine dmumps(id)
            !! the MUMPS driver (double precision): initializes, analyses, factors,
            !! solves, or frees, as `id%job` says
            import :: dmumps_struc
            implicit none
            type(dmumps_struc), intent(inout) :: id !! the MUMPS instance (its inputs and outputs)
        end subroutine dmumps
    end interface

    logical, parameter, public :: sqpopt_has_mumps = .true.  !! whether the library was built with MUMPS (the
                                                             !! `HAS_MUMPS` preprocessor directive), which
                                                             !! `sqpopt_linear_solver_mumps` requires
#else
    logical, parameter, public :: sqpopt_has_mumps = .false. !! whether the library was built with MUMPS (the
                                                             !! `HAS_MUMPS` preprocessor directive), which
                                                             !! `sqpopt_linear_solver_mumps` requires
#endif

    ! the sparse solvers (`options%linear_solver`):
    integer, parameter, public :: sqpopt_linear_solver_mumps = 1 !! MUMPS (needs a library built with it)
    integer, parameter, public :: sqpopt_linear_solver_qdldl = 2 !! (default) QDLDL (always available; see the
                                                                 !! module documentation)

    integer, parameter :: kind_check = 1/merge(1, 0, qdldl_wp == wp) !! (a compile-time error if QDLDL was built
                                                                     !! with another real kind)
    real(wp), parameter :: qdldl_zero_tol = 1.0e-5_wp*epsilon(1.0_wp) !! QDLDL: a pivot at most this times the
                                                                      !! largest element is a null pivot

    public :: sqpopt_linear_solver_available, sqpopt_linear_solver_name

    type, public :: sqpopt_symmetric_solver_type
        !! a sparse symmetric indefinite solver for matrices with one sparsity
        !! pattern (see the module documentation). It holds pointers, so it
        !! must not be copied, and [[symmetric_solver_destroy]] must be called
        !! when it is no longer needed.

        private

        logical,  public :: ready = .false.         !! whether it is initialized
        logical,  public :: factored = .false.      !! whether it holds a factorization (to solve with)
        integer,  public :: n_negative = 0          !! number of negative eigenvalues of the matrix last factored
        integer,  public :: n_null = 0              !! number of zero eigenvalues of the matrix last factored
        integer,  public :: n_factor = 0            !! number of factorizations so far
        integer,  public :: n_solve = 0             !! number of solves so far (including those of the refinement)
        real(wp), public :: time = 0.0_wp           !! wall-clock time spent in the solver so far (seconds)
        logical,  public :: out_of_memory = .false. !! whether an analysis, factorization, or solve ran out of
                                                    !! memory (it stays set)
        integer :: n = 0                            !! order of the matrix
        logical :: analysed = .false.               !! whether the sparsity pattern has been analysed
        integer, public :: backend = 0              !! the solver in use, once initialized (`0`: none)
                                                    !! (`sqpopt_linear_solver_mumps` or `sqpopt_linear_solver_qdldl`)
#ifdef HAS_MUMPS
        type(dmumps_struc), pointer :: id => null() !! the MUMPS instance (with the pattern and the values)
#endif
        type(qdldl_type) :: ldl                     !! the QDLDL solver (with the pattern, the values, and the factors)

        contains

        procedure, public :: initialize => symmetric_solver_initialize
        procedure, public :: factor     => symmetric_solver_factor
        procedure, public :: solve      => symmetric_solver_solve
        procedure, public :: multiply   => symmetric_solver_multiply
        procedure, public :: destroy    => symmetric_solver_destroy

    end type sqpopt_symmetric_solver_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  whether the sparse solver `solver` (a `sqpopt_linear_solver_*` value)
!  can be used in this build: QDLDL always, MUMPS only in a build with it.

    pure logical function sqpopt_linear_solver_available(solver) result(available)

    integer, intent(in) :: solver !! the solver (`sqpopt_linear_solver_*`)

    select case (solver)
    case (sqpopt_linear_solver_qdldl)
        available = .true.
    case (sqpopt_linear_solver_mumps)
        available = sqpopt_has_mumps
    case default
        available = .false.
    end select

    end function sqpopt_linear_solver_available
!*******************************************************************************

!*******************************************************************************
!>
!  the name of the sparse solver `solver` (a `sqpopt_linear_solver_*`
!  value).

    pure function sqpopt_linear_solver_name(solver) result(name)

    integer, intent(in) :: solver !! the solver (`sqpopt_linear_solver_*`)
    character(len=:), allocatable :: name

    select case (solver)
    case (sqpopt_linear_solver_mumps)
        name = 'MUMPS'
    case default
        name = 'QDLDL'
    end select

    end function sqpopt_linear_solver_name
!*******************************************************************************

!*******************************************************************************
!>
!  start the solver for symmetric matrices of order `n` with the sparsity
!  pattern `irow`/`icol`: each entry stands for itself and (off the
!  diagonal) its mirror image, so give each off-diagonal element in one
!  triangle only; entries given more than once are added together.
!  `threads` is the number of OpenMP threads to factor and solve with (see
!  the module documentation; default 1; MUMPS only). `solver` is the
!  sparse solver (a `sqpopt_linear_solver_*` value; default QDLDL).
!  `signs` (QDLDL only) is the expected sign of each row's pivot: `-1` for
!  the rows of a block that should be negative definite (the constraints of
!  a KKT matrix), else `+1` or `0`; each row with `-1` is ordered after all
!  of its neighbours (see the module documentation). `ok` is false if the
!  solver isn't available in this build, or couldn't be started.

    subroutine symmetric_solver_initialize(me, n, irow, icol, ok, threads, solver, signs)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me
    integer,               intent(in)  :: n    !! order of the matrix
    integer, dimension(:), intent(in)  :: irow !! row indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)  :: icol !! column indices of the entries `dimension(nnz)`
    logical,               intent(out) :: ok   !! whether the solver is ready
    integer, optional,     intent(in)  :: threads !! number of OpenMP threads (default `1`; `0`: as the OpenMP
                                                  !! environment says)
    integer, optional,     intent(in)  :: solver  !! the sparse solver (`sqpopt_linear_solver_*`, default QDLDL)
    integer, dimension(:), optional, intent(in) :: signs !! QDLDL: the expected sign of each row's pivot
                                                         !! `dimension(n)` (see above)

#ifdef HAS_MUMPS
    integer :: nnz, alloc_stat
#endif
    integer :: which

    call me%destroy()
    ok = .false.
    which = sqpopt_linear_solver_qdldl
    if (present(solver)) which = solver
    if (.not. sqpopt_linear_solver_available(which)) return

    if (which == sqpopt_linear_solver_qdldl) then
        call qdldl_start()
        return
    end if

#ifdef HAS_MUMPS
    nnz = size(irow)
    allocate(me%id, stat=alloc_stat)
    if (alloc_stat /= 0) then
        nullify(me%id)
        me%out_of_memory = .true.
        return
    end if
    me%id%comm = 0   ! (the sequential library's MPI is a stub: no communicator is needed)
    me%id%sym  = 2   ! general symmetric
    me%id%par  = 1
    me%id%job  = -1
    call dmumps(me%id)
    if (me%id%infog(1) < 0) then
        deallocate(me%id)
        return
    end if

    me%id%icntl(1:3) = -1   ! no printed output
    me%id%icntl(4)   = 0
    me%id%icntl(13)  = 1    ! (needed for the number of negative pivots to be exact)
    me%id%icntl(16)  = 1    ! the number of OpenMP threads (see the module documentation)
    if (present(threads)) me%id%icntl(16) = max(threads, 0)
    me%id%icntl(24)  = 1    ! detect null pivots, rather than fail on a singular matrix

    me%id%n   = n
    me%id%nnz = int(nnz, kind(me%id%nnz))
    allocate(me%id%irn(nnz), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%id%jcn(nnz), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%id%a(nnz), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%id%rhs(n), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call me%destroy()
        me%out_of_memory = .true.
        return
    end if
    me%id%irn = irow
    me%id%jcn = icol

    me%n       = n
    me%backend = sqpopt_linear_solver_mumps
    me%ready   = .true.
    ok = .true.
#endif

    contains

        subroutine qdldl_start()
        !! analyse the pattern with QDLDL: AMD, then (with `signs`) each row
        !! with the sign `-1` moved after all of its neighbours
        integer :: istat
        integer, dimension(:), allocatable :: perm
        if (present(threads)) continue   ! (QDLDL is single-threaded)
        call me%ldl%analyze(n, irow, icol, istat, ordering=qdldl_order_amd)
        if (istat == qdldl_success .and. present(signs)) then
            if (size(signs) == n) then
                if (any(signs == -1)) then
                    call me%ldl%get_permutation(perm)
                    call delay_negative_rows(n, irow, icol, signs, perm)
                    call me%ldl%analyze(n, irow, icol, istat, perm=perm)
                end if
            end if
        end if
        if (istat /= qdldl_success) then
            if (istat == qdldl_error_out_of_memory) me%out_of_memory = .true.
            call me%ldl%destroy()
            return
        end if
        ! (see [[symmetric_solver_factor]] for the null pivots; a pivot of the
        ! wrong sign is kept, so the inertia is that of the matrix)
        me%ldl%reg_eps        = 0.0_wp
        me%ldl%zero_pivot_tol = qdldl_zero_tol
        me%ldl%max_refine     = 0   ! (the refinement is this module's, for both solvers)
        me%n       = n
        me%backend = sqpopt_linear_solver_qdldl
        me%ready   = .true.
        ok = .true.
        end subroutine qdldl_start

    end subroutine symmetric_solver_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  change the fill-reducing order `perm` (row `perm(k)` of the matrix is
!  the `k`-th to be eliminated) so that each row `i` with `signs(i) == -1`
!  comes after all of its neighbours (the rows `j` with an entry `(i,j)`),
!  and otherwise as early as it can. In a KKT matrix
!  \( \begin{bmatrix} H & J^T \\ J & 0 \end{bmatrix} \) with \( H \) positive
!  definite, the pivot of a constraint's row is then an element of
!  \( -J H^{-1} J^T \) (less the earlier constraints' parts), which is zero
!  only if the rows of \( J \) are dependent. Without the delay, a
!  constraint's row eliminated before its variables has a zero pivot even
!  when the matrix is nonsingular, which a solver without pivoting (QDLDL)
!  can't get around.

    subroutine delay_negative_rows(n, irow, icol, signs, perm)

    integer,               intent(in)    :: n     !! order of the matrix
    integer, dimension(:), intent(in)    :: irow  !! row indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)    :: icol  !! column indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)    :: signs !! the expected sign of each row's pivot `dimension(n)`
    integer, dimension(:), intent(inout) :: perm  !! the order, changed in place `dimension(n)`

    integer, dimension(:), allocatable :: ptr, adj, waiting, new_perm, fill
    logical, dimension(:), allocatable :: done
    integer :: k, i, j, v, c, p, nnew

    ! the neighbours of each row (both directions), without the diagonal:
    allocate(ptr(n+1), fill(n), waiting(n), done(n), new_perm(n))
    ptr = 0
    do k = 1, size(irow)
        i = irow(k)
        j = icol(k)
        if (i == j) cycle
        ptr(i) = ptr(i) + 1
        ptr(j) = ptr(j) + 1
    end do
    p = 1
    do i = 1, n
        fill(i) = p
        p = p + ptr(i)
        ptr(i) = fill(i)
    end do
    ptr(n+1) = p
    allocate(adj(max(p-1, 0)))
    do k = 1, size(irow)
        i = irow(k)
        j = icol(k)
        if (i == j) cycle
        adj(fill(i)) = j
        fill(i) = fill(i) + 1
        adj(fill(j)) = i
        fill(j) = fill(j) + 1
    end do

    ! how many neighbours of each delayed row are still to come (duplicates
    ! counted, consistently, in both places):
    waiting = 0
    do i = 1, n
        if (signs(i) /= -1) cycle
        do p = ptr(i), ptr(i+1) - 1
            if (signs(adj(p)) /= -1) waiting(i) = waiting(i) + 1
        end do
    end do

    done = .false.
    nnew = 0
    do k = 1, n
        v = perm(k)
        if (done(v)) cycle
        if (signs(v) == -1) then
            if (waiting(v) == 0) call emit(v)   ! (else it is emitted after its last neighbour)
        else
            call emit(v)
            do p = ptr(v), ptr(v+1) - 1
                c = adj(p)
                if (signs(c) /= -1) cycle
                waiting(c) = waiting(c) - 1
                if (waiting(c) == 0 .and. .not. done(c)) call emit(c)
            end do
        end if
    end do
    ! (anything left: rows whose neighbours were all delayed rows themselves)
    do k = 1, n
        if (.not. done(perm(k))) call emit(perm(k))
    end do
    perm = new_perm

    contains

        subroutine emit(row)
        !! append `row` to the new order
        integer, intent(in) :: row !! the row
        nnew = nnew + 1
        new_perm(nnew) = row
        done(row) = .true.
        end subroutine emit

    end subroutine delay_negative_rows
!*******************************************************************************

!*******************************************************************************
!>
!  factor the matrix with the values `val` (in the order of the pattern
!  given to [[symmetric_solver_initialize]]). On success, `n_negative` and
!  `n_null` are its numbers of negative and zero eigenvalues, and
!  [[symmetric_solver_solve]] can be called. With MUMPS, the pattern is
!  analysed in the first call, with that call's values; if MUMPS runs out
!  of its estimated workspace, the factorization is repeated with twice
!  the allowance, up to 20 times (as IPOPT does), and the larger allowance
!  is kept for the following factorizations. With QDLDL (analysed by
!  [[symmetric_solver_initialize]]), a pivot of at most `qdldl_zero_tol`
!  times the largest element counts as a null pivot.

    subroutine symmetric_solver_factor(me, val, ok)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: val !! values of the entries `dimension(nnz)`
    logical,                intent(out) :: ok  !! whether the factorization succeeded

#ifdef HAS_MUMPS
    integer :: attempt
    integer, parameter :: max_attempts = 20 !! (each doubles the workspace increase; a very indefinite matrix,
                                            !! whose pivots are delayed, can need many times the estimate)
#endif
    integer(int64) :: t0, t1, rate
    integer :: istat
    integer :: null_row !! QDLDL: the row of the first null pivot (0: none)
    real(wp) :: amax

    ok = .false.
    me%factored = .false.
    if (.not. me%ready) return

    call system_clock(t0, rate)
    if (me%backend == sqpopt_linear_solver_qdldl) then
        ! first without replacing null pivots: the first one stops the
        ! factorization, and tells its row
        me%ldl%regularize = .false.
        call me%ldl%factor(val, istat)
        me%n_factor = me%n_factor + 1
        null_row = 0
        if (istat == qdldl_error_zero_pivot) then
            null_row = me%ldl%zero_pivot_column
            ! again, with each null pivot replaced by a large one (which makes its
            ! row and column of the factors negligible, so that a solve returns
            ! one of the solutions of a consistent system, as with MUMPS):
            amax = 0.0_wp
            if (size(val) > 0) amax = maxval(abs(val))
            me%ldl%regularize = .true.
            me%ldl%reg_delta  = max(amax, 1.0_wp)/sqrt(epsilon(1.0_wp))
            call me%ldl%factor(val, istat)   ! (counted as the same factorization)
        end if
        call system_clock(t1)
        me%time = me%time + real(t1 - t0, wp)/real(rate, wp)
        if (istat /= qdldl_success) return
        me%analysed   = .true.
        me%factored   = .true.
        me%n_negative = me%ldl%n_negative
        me%n_null     = me%ldl%n_zero
        if (null_row > 0) then
            if (.not. faithful()) then
                ! a null pivot doesn't always mean that the matrix is singular: this
                ! order may meet a zero on the way (e.g. a Hessian with a zero
                ! diagonal, which needs a 2x2 pivot). Then the factors, with that
                ! pivot replaced, are not those of the matrix, and without pivoting
                ! its inertia is unknown: it is reported as a negative eigenvalue
                ! (the inertia control then shifts the Hessian, and the direct QP
                ! sees a nonconvex face)
                me%n_negative = me%n_negative + me%n_null
                me%n_null     = 0
            end if
        end if
        ok = .true.
        return
    end if

#ifdef HAS_MUMPS
    me%id%a = val
    do attempt = 1, max_attempts
        me%id%job = merge(2, 4, me%analysed)   ! (4: analyse, then factor)
        call dmumps(me%id)
        if (me%id%infog(1) /= -8 .and. me%id%infog(1) /= -9) exit
        me%id%icntl(14) = 2*max(me%id%icntl(14), 20)   ! more workspace (a percentage of the estimate)
    end do
    me%n_factor = me%n_factor + 1
    call system_clock(t1)
    me%time = me%time + real(t1 - t0, wp)/real(rate, wp)
    if (me%id%infog(1) < 0) then
        ! (-13: an allocation failed; -5, -7: an allocation of the analysis failed;
        ! -8, -9: still not enough workspace)
        if (any(me%id%infog(1) == [-13, -5, -7, -8, -9])) me%out_of_memory = .true.
        return
    end if
    me%analysed   = .true.
    me%factored   = .true.
    me%n_negative = me%id%infog(12)
    me%n_null     = me%id%infog(28)
    ok = .true.
#endif

    contains

        logical function faithful()
        !! QDLDL, after null pivots were replaced: whether the factors are still
        !! those of the matrix (as they are if the matrix is singular, since the
        !! rest of a null pivot's row is then zero too): for a fixed `x` and
        !! `b = A x`, the solve with the factors satisfies `A x' = b` (`x'` may
        !! differ from `x` along the null space)
        real(wp), parameter :: rtol = 1.0e-6_wp !! relative tolerance on the residual
        real(wp), dimension(:), allocatable :: x, b, ax
        integer :: i, solve_istat
        allocate(x(me%n), b(me%n), ax(me%n))
        do i = 1, me%n
            x(i) = 1.0_wp + real(mod(i, 7), wp)/7.0_wp
        end do
        call me%ldl%multiply(x, b)
        x = b
        call me%ldl%solve(x, solve_istat, refine=0)
        faithful = solve_istat == qdldl_success
        if (.not. faithful) return
        call me%ldl%multiply(x, ax)
        faithful = maxval(abs(ax - b)) <= rtol*(amax*maxval(abs(x)) + maxval(abs(b)))
        end function faithful

    end subroutine symmetric_solver_factor
!*******************************************************************************

!*******************************************************************************
!>
!  solve \( A x = b \) with the factorization of [[symmetric_solver_factor]]:
!  `b` is overwritten by `x`. The solution is improved by iterative
!  refinement (up to `max_refine` steps, each one more solve, while the
!  residual \( b - Ax \) is above roundoff and each step at least halves
!  it), unless `refine` is false. `ok` is false if there is no
!  factorization, the solver reports that the solve failed (`b` is then
!  unchanged, and `out_of_memory` is set if an allocation failed), or the
!  solution isn't finite. If a solve of the refinement fails, the
!  refinement stops, and the solution found so far is returned.

    subroutine symmetric_solver_solve(me, b, ok, refine)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: b      !! the right-hand side, overwritten by the solution `dimension(n)`
    logical,                intent(out)   :: ok     !! whether a solution was found
    logical, optional,      intent(in)    :: refine !! whether to refine the solution (default `.true.`)

    integer, parameter :: max_refine = 2 !! maximum number of refinement steps
    real(wp), dimension(:), allocatable :: rhs, r, ax, x_new
    real(wp) :: rnorm, rnorm_new, tol
    integer :: step
    logical :: do_refine, solved
    integer(int64) :: t0, t1, rate

    ok = .false.
    if (.not. (me%ready .and. me%factored)) return

    call system_clock(t0, rate)
    do_refine = .true.
    if (present(refine)) do_refine = refine

    rhs = b
    call back_solve(b, solved)
    ok = solved
    if (.not. solved) then
        b = rhs
    else if (do_refine) then
        allocate(r(me%n), ax(me%n))
        tol = 10.0_wp*epsilon(1.0_wp)*max(maxval(abs(rhs)), tiny(1.0_wp))
        call me%multiply(b, ax)
        r = rhs - ax
        rnorm = maxval(abs(r))
        do step = 1, max_refine
            if (.not. rnorm > tol) exit
            call back_solve(r, solved)
            if (.not. solved) exit
            x_new = b + r
            call me%multiply(x_new, ax)
            r = rhs - ax
            rnorm_new = maxval(abs(r))
            if (.not. rnorm_new <= 0.5_wp*rnorm) exit
            b = x_new
            rnorm = rnorm_new
        end do
    end if
    if (ok) ok = sqpopt_all_finite(b)
    call system_clock(t1)
    me%time = me%time + real(t1 - t0, wp)/real(rate, wp)

    contains

        subroutine back_solve(v, success)
        !! one solve with the factors: `v` is overwritten by \( A^{-1} v \)
        real(wp), dimension(:), intent(inout) :: v       !! the right-hand side, overwritten by the solution
        logical,                intent(out)   :: success !! whether the solver solved it (if not, `v` is unchanged)
        integer :: istat
        success = .false.
        me%n_solve = me%n_solve + 1
        if (me%backend == sqpopt_linear_solver_qdldl) then
            call me%ldl%solve(v, istat, refine=0)
            success = istat == qdldl_success
            return
        end if
#ifdef HAS_MUMPS
        me%id%rhs = v
        me%id%job = 3
        call dmumps(me%id)
        success = me%id%infog(1) >= 0
        if (success) then
            v = me%id%rhs
        else if (me%id%infog(1) == -13) then
            me%out_of_memory = .true.   ! (an allocation failed)
        end if
#endif
        end subroutine back_solve

    end subroutine symmetric_solver_solve
!*******************************************************************************

!*******************************************************************************
!>
!  the product \( y = A x \) with the matrix last given to
!  [[symmetric_solver_factor]] (used for the residuals of the refinement,
!  and for testing).

    subroutine symmetric_solver_multiply(me, x, y)

    class(sqpopt_symmetric_solver_type), intent(in) :: me
    real(wp), dimension(:), intent(in)  :: x !! the vector `dimension(n)`
    real(wp), dimension(:), intent(out) :: y !! the product `dimension(n)`

#ifdef HAS_MUMPS
    integer :: k, i, j
#endif
    integer :: istat

    y = 0.0_wp
    if (.not. me%ready) return
    if (me%backend == sqpopt_linear_solver_qdldl) then
        call me%ldl%multiply(x, y, istat)
        if (istat /= qdldl_success) y = 0.0_wp
        return
    end if
#ifdef HAS_MUMPS
    do k = 1, size(me%id%a)
        i = me%id%irn(k)
        j = me%id%jcn(k)
        y(i) = y(i) + me%id%a(k)*x(j)
        if (i /= j) y(j) = y(j) + me%id%a(k)*x(i)
    end do
#endif

    end subroutine symmetric_solver_multiply
!*******************************************************************************

!*******************************************************************************
!>
!  free the solver, and return it to its initial state. It does nothing
!  to a solver that was never initialized.

    subroutine symmetric_solver_destroy(me)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me

#ifdef HAS_MUMPS
    if (associated(me%id)) then
        if (associated(me%id%irn)) deallocate(me%id%irn)
        if (associated(me%id%jcn)) deallocate(me%id%jcn)
        if (associated(me%id%a))   deallocate(me%id%a)
        if (associated(me%id%rhs)) deallocate(me%id%rhs)
        me%id%job = -2
        call dmumps(me%id)
        deallocate(me%id)
    end if
#endif
    call me%ldl%destroy()
    me%ready      = .false.
    me%factored   = .false.
    me%analysed   = .false.
    me%backend    = 0
    me%n          = 0
    me%n_negative = 0
    me%n_null     = 0
    me%n_factor   = 0
    me%n_solve    = 0
    me%time       = 0.0_wp
    me%out_of_memory = .false.

    end subroutine symmetric_solver_destroy
!*******************************************************************************

    end module sqpopt_symmetric_solver_module
!*******************************************************************************
