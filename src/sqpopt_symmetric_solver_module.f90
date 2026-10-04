!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  A direct solver for symmetric indefinite linear systems
!  \( A x = b \), with one of four backends
!  (`options%linear_solver`), each in a module of its own, behind the
!  interface of [[sqpopt_sparse_ldl_module]]:
!
!  * `sqpopt_linear_solver_qdldl` (the default, in every build): QDLDL
!    ([[sqpopt_qdldl_ldl_module]]), without pivoting, single-threaded, with
!    very little overhead, always available;
!  * `sqpopt_linear_solver_mumps`: MUMPS ([[sqpopt_mumps_ldl_module]]), with
!    pivoting and threads, in a library built with it (`sqpopt_has_mumps`);
!  * `sqpopt_linear_solver_dense` (opt-in): the matrix as a dense array
!    ([[sqpopt_dense_ldl_module]]), with pivoting and the exact inertia, for
!    small matrices only (\( O(n^3) \) per factorization; it refuses an
!    order above `dense_max_order`);
!  * `sqpopt_linear_solver_lapack` (opt-in): the same dense matrix, factored
!    by LAPACK's `DSYTRF` ([[sqpopt_lapack_ldl_type]]), in a library built
!    with LAPACK (`sqpopt_has_lapack`).
!
!  This module chooses the backend, and does what is the same for both:
!  the iterative refinement of the solves, the counts, and the timing.
!
!  **Which one.** Measured (release; see the user guide's "Sparse solver"
!  section): on the banded and chained problems of
!  `example/benchmark_large.f90` (5,000 to 1,000,000 variables) both take
!  the same iterations, and QDLDL's factorizations take 7 to 30 times less
!  time; on the Hock-Schittkowski problems the quasi-Newton Hessians give
!  the same results with both (QDLDL's faster), but the exact Hessian with
!  inertia control needs about 45% more evaluations with QDLDL, which can't
!  tell the inertia of a Hessian with zeros on its diagonal (see
!  `inertia_known`). On grids (`example/sparse_solvers.f90`), MUMPS
!  refactors 2 to 10 times faster in 2-D and 9 to 34 times in 3-D, where
!  the factors have large dense blocks (and its threads help), though its
!  first factorization, with its analysis, can cost much more (72 s against
!  4.5 s for a 2-D KKT matrix of order 735,000). So MUMPS is the better
!  choice for the exact Hessian with inertia control, and for problems
!  coupled in two or three dimensions that refactor many times. For a small
!  problem with the exact Hessian and inertia control, the dense backend
!  gives MUMPS's results (on the Hock-Schittkowski problems, 9,452
!  evaluations against MUMPS's 9,400 and QDLDL's 13,654) at QDLDL's speed.
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
!     negative and `n_null` zero eigenvalues, if `inertia_known` (else those
!     are the counts of the factors, which are not the matrix's).
!  3. [[symmetric_solver_solve]], any number of times, for right-hand sides.
!     Each solution is improved by iterative refinement.
!  4. [[symmetric_solver_destroy]] to free it. A backend may hold pointers
!     (MUMPS), so it must not be copied.
!
!  A singular matrix is not an error: its null pivots are counted
!  (`n_null`), and a solve then returns one of the solutions of a consistent
!  system.

    module sqpopt_symmetric_solver_module

    use, intrinsic :: iso_fortran_env, only: int64
    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_all_finite
    use sqpopt_sparse_ldl_module, only: sqpopt_sparse_ldl_type
    use sqpopt_qdldl_ldl_module, only: sqpopt_qdldl_ldl_type
    use sqpopt_mumps_ldl_module, only: sqpopt_mumps_ldl_type, sqpopt_has_mumps
    use sqpopt_dense_ldl_module, only: sqpopt_dense_ldl_type, sqpopt_lapack_ldl_type, sqpopt_has_lapack

    implicit none

    private

    public :: sqpopt_has_mumps   ! (from [[sqpopt_mumps_ldl_module]])
    public :: sqpopt_has_lapack  ! (from [[sqpopt_dense_ldl_module]])

    ! the sparse solvers (`options%linear_solver`):
    integer, parameter, public :: sqpopt_linear_solver_mumps = 1 !! MUMPS (needs a library built with it)
    integer, parameter, public :: sqpopt_linear_solver_qdldl = 2 !! (default) QDLDL (always available; see the
                                                                 !! module documentation)
    integer, parameter, public :: sqpopt_linear_solver_dense = 3 !! dense, for small matrices (always available;
                                                                 !! see [[sqpopt_dense_ldl_module]])
    integer, parameter, public :: sqpopt_linear_solver_lapack = 4 !! dense, factored by LAPACK's `DSYTRF`, for small
                                                                  !! matrices (needs a library built with LAPACK)

    public :: sqpopt_linear_solver_available, sqpopt_linear_solver_name

    type, public :: sqpopt_symmetric_solver_type
        !! a sparse symmetric indefinite solver for matrices with one sparsity
        !! pattern (see the module documentation). Its backend may hold
        !! pointers, so it must not be copied, and [[symmetric_solver_destroy]]
        !! must be called when it is no longer needed.

        private

        logical,  public :: ready = .false.         !! whether it is initialized
        logical,  public :: factored = .false.      !! whether it holds a factorization (to solve with)
        integer,  public :: n_negative = 0          !! of the last factorization: the number of negative eigenvalues
        integer,  public :: n_null = 0              !! of the last factorization: the number of zero eigenvalues
        logical,  public :: inertia_known = .true.  !! of the last factorization: whether `n_negative` and `n_null`
                                                    !! are the matrix's (a backend without pivoting may not be able
                                                    !! to tell: see [[sqpopt_qdldl_ldl_module]])
        integer,  public :: n_factor = 0            !! number of factorizations so far
        integer,  public :: n_solve = 0             !! number of solves so far (including those of the refinement)
        real(wp), public :: time = 0.0_wp           !! wall-clock time spent in the solver so far (seconds)
        logical,  public :: out_of_memory = .false. !! whether an analysis, factorization, or solve ran out of
                                                    !! memory (it stays set)
        integer,  public :: backend = 0             !! the solver in use, once initialized (`0`: none)
                                                    !! (`sqpopt_linear_solver_mumps` or `sqpopt_linear_solver_qdldl`)
        integer :: n = 0                            !! order of the matrix
        class(sqpopt_sparse_ldl_type), allocatable :: ldl !! the backend

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
    case (sqpopt_linear_solver_qdldl, sqpopt_linear_solver_dense)
        available = .true.
    case (sqpopt_linear_solver_mumps)
        available = sqpopt_has_mumps
    case (sqpopt_linear_solver_lapack)
        available = sqpopt_has_lapack
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
    case (sqpopt_linear_solver_dense)
        name = 'dense'
    case (sqpopt_linear_solver_lapack)
        name = 'LAPACK'
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
!  `solver` is the backend (a `sqpopt_linear_solver_*` value; default
!  QDLDL). `threads` is the number of OpenMP threads to factor and solve
!  with (MUMPS only; see [[sqpopt_mumps_ldl_module]]; default 1). `signs`
!  is the expected sign of each row's pivot (QDLDL only): `-1` for the rows
!  of a block that should be negative definite (the constraints of a KKT
!  matrix), else `+1` or `0` (see [[delay_negative_rows]]). `ok` is false
!  if the solver isn't available in this build, or couldn't be started.

    subroutine symmetric_solver_initialize(me, n, irow, icol, ok, threads, solver, signs)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me
    integer,               intent(in)  :: n    !! order of the matrix
    integer, dimension(:), intent(in)  :: irow !! row indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)  :: icol !! column indices of the entries `dimension(nnz)`
    logical,               intent(out) :: ok   !! whether the solver is ready
    integer, optional,     intent(in)  :: threads !! number of OpenMP threads (default `1`; `0`: as the OpenMP
                                                  !! environment says)
    integer, optional,     intent(in)  :: solver  !! the sparse solver (`sqpopt_linear_solver_*`, default QDLDL)
    integer, dimension(:), optional, intent(in) :: signs !! the expected sign of each row's pivot `dimension(n)`
                                                         !! (see above)

    integer :: which, alloc_stat

    call me%destroy()
    ok = .false.
    which = sqpopt_linear_solver_qdldl
    if (present(solver)) which = solver
    if (.not. sqpopt_linear_solver_available(which)) return

    select case (which)
    case (sqpopt_linear_solver_mumps)
        allocate(sqpopt_mumps_ldl_type :: me%ldl, stat=alloc_stat)
    case (sqpopt_linear_solver_dense)
        allocate(sqpopt_dense_ldl_type :: me%ldl, stat=alloc_stat)
    case (sqpopt_linear_solver_lapack)
        allocate(sqpopt_lapack_ldl_type :: me%ldl, stat=alloc_stat)
    case default
        allocate(sqpopt_qdldl_ldl_type :: me%ldl, stat=alloc_stat)
    end select
    if (alloc_stat /= 0) then
        me%out_of_memory = .true.
        return
    end if

    call me%ldl%start(n, irow, icol, ok, threads=threads, signs=signs)
    if (.not. ok) then
        me%out_of_memory = me%ldl%out_of_memory
        call me%ldl%free()
        deallocate(me%ldl)
        return
    end if
    me%n       = n
    me%backend = which
    me%ready   = .true.

    end subroutine symmetric_solver_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  factor the matrix with the values `val` (in the order of the pattern
!  given to [[symmetric_solver_initialize]]). On success, `n_negative`,
!  `n_null`, and `inertia_known` describe its inertia (see the module
!  documentation), and [[symmetric_solver_solve]] can be called.

    subroutine symmetric_solver_factor(me, val, ok)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: val !! values of the entries `dimension(nnz)`
    logical,                intent(out) :: ok  !! whether the factorization succeeded

    integer(int64) :: t0, t1, rate

    ok = .false.
    me%factored = .false.
    if (.not. me%ready) return

    call system_clock(t0, rate)
    call me%ldl%refactor(val, ok)
    me%n_factor = me%n_factor + 1
    call system_clock(t1)
    me%time = me%time + real(t1 - t0, wp)/real(rate, wp)
    if (me%ldl%out_of_memory) me%out_of_memory = .true.
    if (.not. ok) return
    me%factored      = .true.
    me%n_negative    = me%ldl%n_negative
    me%n_null        = me%ldl%n_null
    me%inertia_known = me%ldl%inertia_known

    end subroutine symmetric_solver_factor
!*******************************************************************************

!*******************************************************************************
!>
!  solve \( A x = b \) with the factorization of [[symmetric_solver_factor]]:
!  `b` is overwritten by `x`. The solution is improved by iterative
!  refinement (up to `max_refine` steps, each one more solve, while the
!  residual \( b - Ax \) is above roundoff and each step at least halves
!  it), unless `refine` is false. `ok` is false if there is no
!  factorization, the backend reports that the solve failed (`b` is then
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
        logical,                intent(out)   :: success !! whether it was solved (if not, `v` is unchanged)
        me%n_solve = me%n_solve + 1
        call me%ldl%back_solve(v, success)
        if (me%ldl%out_of_memory) me%out_of_memory = .true.
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

    y = 0.0_wp
    if (.not. me%ready) return
    call me%ldl%multiply(x, y)

    end subroutine symmetric_solver_multiply
!*******************************************************************************

!*******************************************************************************
!>
!  free the solver, and return it to its initial state. It does nothing
!  to a solver that was never initialized.

    subroutine symmetric_solver_destroy(me)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me

    if (allocated(me%ldl)) then
        call me%ldl%free()
        deallocate(me%ldl)
    end if
    me%ready         = .false.
    me%factored      = .false.
    me%backend       = 0
    me%n             = 0
    me%n_negative    = 0
    me%n_null        = 0
    me%inertia_known = .true.
    me%n_factor      = 0
    me%n_solve       = 0
    me%time          = 0.0_wp
    me%out_of_memory = .false.

    end subroutine symmetric_solver_destroy
!*******************************************************************************

    end module sqpopt_symmetric_solver_module
!*******************************************************************************
