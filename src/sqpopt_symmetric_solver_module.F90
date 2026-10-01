!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  A sparse direct solver for symmetric indefinite linear systems
!  \( A x = b \), by an \( LDL^T \) factorization. This is the library's
!  only interface to [MUMPS](https://mumps-solver.org), and the only source
!  file that depends on the `HAS_MUMPS` preprocessor directive: without it,
!  `sqpopt_has_mumps` is false, [[symmetric_solver_initialize]] always
!  fails, and everything built on this solver (see [[sqpopt_kkt_module]])
!  is unavailable. With it, the library must be linked with the
!  sequential, double precision MUMPS library, and compiled in double
!  precision (the default kind, see [[sqpopt_kinds]]).
!
!  How it is used:
!
!  1. [[symmetric_solver_initialize]] with the sparsity pattern (the entries
!     of one triangle, in coordinate form; entries given more than once are
!     added together).
!  2. [[symmetric_solver_factor]] with the values of those entries. The
!     pattern is analysed (ordered) in the first call only, so a matrix
!     whose values change but whose pattern doesn't is cheap to refactor.
!     The factorization gives the *inertia* of the matrix: `n_negative`
!     negative and `n_null` zero eigenvalues.
!  3. [[symmetric_solver_solve]], any number of times, for right-hand sides.
!     Each solution is improved by iterative refinement.
!  4. [[symmetric_solver_destroy]] to free it. The type holds pointers, so
!     it must not be copied.
!
!  MUMPS is run on one OpenMP thread, whatever `OMP_NUM_THREADS` is: on
!  small matrices its threads cost far more time than they save (the 305
!  Hock-Schittkowski problems took 24.5 s with 10 threads, and 1.4 s with
!  one). A singular matrix is not an error: its null pivots are counted
!  (`n_null`), and a solve then returns one of the solutions of a consistent
!  system.

    module sqpopt_symmetric_solver_module

    use, intrinsic :: iso_fortran_env, only: int64
    use sqpopt_kinds, only: wp => sqpopt_module_wp

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
                                                             !! the options that factor a matrix require
#else
    logical, parameter, public :: sqpopt_has_mumps = .false. !! whether the library was built with MUMPS (the
                                                             !! `HAS_MUMPS` preprocessor directive), which
                                                             !! the options that factor a matrix require
#endif

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
        logical,  public :: out_of_memory = .false. !! whether an analysis or factorization ran out of memory (it
                                                    !! stays set)
        integer :: n = 0                            !! order of the matrix
        logical :: analysed = .false.               !! whether the sparsity pattern has been analysed
#ifdef HAS_MUMPS
        type(dmumps_struc), pointer :: id => null() !! the MUMPS instance (with the pattern and the values)
#endif

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
!  start the solver for symmetric matrices of order `n` with the sparsity
!  pattern `irow`/`icol`: each entry stands for itself and (off the
!  diagonal) its mirror image, so give each off-diagonal element in one
!  triangle only; entries given more than once are added together. `ok`
!  is false if the library was built without MUMPS, or MUMPS couldn't be
!  started.

    subroutine symmetric_solver_initialize(me, n, irow, icol, ok)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me
    integer,               intent(in)  :: n    !! order of the matrix
    integer, dimension(:), intent(in)  :: irow !! row indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)  :: icol !! column indices of the entries `dimension(nnz)`
    logical,               intent(out) :: ok   !! whether the solver is ready

#ifdef HAS_MUMPS
    integer :: nnz, alloc_stat
#endif

    call me%destroy()
    ok = .false.

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
    me%id%icntl(16)  = 1    ! one OpenMP thread (see the module documentation)
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

    me%n     = n
    me%ready = .true.
    ok = .true.
#endif

    end subroutine symmetric_solver_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  factor the matrix with the values `val` (in the order of the pattern
!  given to [[symmetric_solver_initialize]]). On success, `n_negative` and
!  `n_null` are its numbers of negative and zero eigenvalues, and
!  [[symmetric_solver_solve]] can be called. The pattern is analysed in the
!  first call, with that call's values. If MUMPS runs out of its estimated
!  workspace, the factorization is repeated with more.

    subroutine symmetric_solver_factor(me, val, ok)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: val !! values of the entries `dimension(nnz)`
    logical,                intent(out) :: ok  !! whether the factorization succeeded

#ifdef HAS_MUMPS
    integer :: attempt
    integer, parameter :: max_attempts = 6 !! (each doubles the workspace increase)
    integer(int64) :: t0, t1, rate
#endif

    ok = .false.
    me%factored = .false.
    if (.not. me%ready) return

#ifdef HAS_MUMPS
    call system_clock(t0, rate)
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
        ! (-13: an allocation failed; -8, -9: still not enough workspace)
        if (any(me%id%infog(1) == [-13, -8, -9])) me%out_of_memory = .true.
        return
    end if
    me%analysed   = .true.
    me%factored   = .true.
    me%n_negative = me%id%infog(12)
    me%n_null     = me%id%infog(28)
    ok = .true.
#endif

    end subroutine symmetric_solver_factor
!*******************************************************************************

!*******************************************************************************
!>
!  solve \( A x = b \) with the factorization of [[symmetric_solver_factor]]:
!  `b` is overwritten by `x`. The solution is improved by iterative
!  refinement (up to `max_refine` steps, each one more solve, while the
!  residual \( b - Ax \) is above roundoff and each step at least halves
!  it), unless `refine` is false. `ok` is false if there is no
!  factorization, or the solution isn't finite.

    subroutine symmetric_solver_solve(me, b, ok, refine)

    class(sqpopt_symmetric_solver_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: b      !! the right-hand side, overwritten by the solution `dimension(n)`
    logical,                intent(out)   :: ok     !! whether a solution was found
    logical, optional,      intent(in)    :: refine !! whether to refine the solution (default `.true.`)

#ifdef HAS_MUMPS
    integer, parameter :: max_refine = 2 !! maximum number of refinement steps
    real(wp), dimension(:), allocatable :: rhs, r, ax, x_new
    real(wp) :: rnorm, rnorm_new, tol
    integer :: step
    logical :: do_refine
    integer(int64) :: t0, t1, rate
#endif

    ok = .false.
    if (.not. (me%ready .and. me%factored)) return

#ifdef HAS_MUMPS
    call system_clock(t0, rate)
    do_refine = .true.
    if (present(refine)) do_refine = refine

    rhs = b
    call back_solve(b)
    if (do_refine) then
        allocate(r(me%n), ax(me%n))
        tol = 10.0_wp*epsilon(1.0_wp)*max(maxval(abs(rhs)), tiny(1.0_wp))
        call me%multiply(b, ax)
        r = rhs - ax
        rnorm = maxval(abs(r))
        do step = 1, max_refine
            if (.not. rnorm > tol) exit
            call back_solve(r)
            x_new = b + r
            call me%multiply(x_new, ax)
            r = rhs - ax
            rnorm_new = maxval(abs(r))
            if (.not. rnorm_new <= 0.5_wp*rnorm) exit
            b = x_new
            rnorm = rnorm_new
        end do
    end if
    ok = all(abs(b) <= huge(1.0_wp))
    call system_clock(t1)
    me%time = me%time + real(t1 - t0, wp)/real(rate, wp)

    contains

        subroutine back_solve(v)
        !! one solve with the factors: `v` is overwritten by \( A^{-1} v \)
        real(wp), dimension(:), intent(inout) :: v !! the right-hand side, overwritten by the solution
        me%id%rhs = v
        me%id%job = 3
        call dmumps(me%id)
        v = me%id%rhs
        me%n_solve = me%n_solve + 1
        end subroutine back_solve
#endif

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

    y = 0.0_wp
#ifdef HAS_MUMPS
    if (.not. me%ready) return
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
    me%ready      = .false.
    me%factored   = .false.
    me%analysed   = .false.
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
