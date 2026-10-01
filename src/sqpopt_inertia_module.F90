!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Inertia control for the exact Hessian (`options%inertia_control`), with
!  the sparse symmetric indefinite \( LDL^T \) factorization of
!  [MUMPS](https://mumps-solver.org). It is only available in a build with
!  the `HAS_MUMPS` preprocessor directive (see `sqpopt_has_mumps`), which
!  must be linked with the sequential, double precision MUMPS library.
!
!  The QP subproblem's step is only a minimizer if the Hessian \( H \) is
!  positive definite on the null space of the working set's constraints,
!  \( Z^THZ \succ 0 \). That is a property of the *inertia* (the numbers of
!  positive, negative, and zero eigenvalues) of the KKT matrix
!
!  $$ K = \begin{bmatrix} H + \delta I & J_a^T \\ J_a & 0 \end{bmatrix} $$
!
!  where \( J_a \) holds the working set's rows of the Jacobian, in the
!  variables that are not at a bound of the working set. By Sylvester's law
!  of inertia, \( K \) has as many negative eigenvalues as \( J_a \) has
!  (independent) rows, plus one for each negative eigenvalue of
!  \( Z^T(H+\delta I)Z \). The factorization gives that number (its negative
!  pivots), so the smallest shift \( \delta \) that leaves no negative
!  curvature is found by refactoring ([[inertia_correct]]), much as in
!  IPOPT: first the current shift, then `shift_min` times the size of the
!  Hessian (or a third of the last shift that was needed, if that is
!  larger), then 8 times as much at each attempt, up to `shift_max`. (IPOPT
!  multiplies by 100 until a shift has been needed once. On the
!  Hock-Schittkowski problems that was no better: 273 solved and 3 failed,
!  against 274 and 2, with 1% more function evaluations. Factors of 2, 4,
!  and 16 gave 271 or 272 solved.)
!
!  So that the sparsity pattern is the same for every working set (it is
!  analysed only once per solve), \( K \) always has order `n+m`: a variable
!  at a bound of the working set has its row and column replaced by those
!  of the identity, and a constraint outside it by those of minus the
!  identity. Then no negative curvature is left when the factorization has
!  at most `m` negative pivots.
!
!  This only decides the shift. The QP subproblems are still solved by the
!  matrix-free solvers of [[sqpopt_qp_solver_module]], with the shifted
!  Hessian (see [[sqpopt_iterate]]).
!
!  MUMPS is run on one OpenMP thread, whatever `OMP_NUM_THREADS` is: on
!  small matrices its threads cost far more time than they save (the 305
!  Hock-Schittkowski problems took 24.5 s with 10 threads, and 1.4 s with
!  one).

    module sqpopt_inertia_module

    use, intrinsic :: iso_fortran_env, only: real64
    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix
    use sqpopt_hessian_module, only: sqpopt_hessian_type

    implicit none

    private

#ifdef HAS_MUMPS
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
                                                             !! `options%inertia_control` requires
#else
    logical, parameter, public :: sqpopt_has_mumps = .false. !! whether the library was built with MUMPS (the
                                                             !! `HAS_MUMPS` preprocessor directive), which
                                                             !! `options%inertia_control` requires
#endif

    type, public :: sqpopt_inertia_type
        !! the inertia control of one solve: the MUMPS instance, and the state
        !! of the shift's search (see the module documentation). It holds
        !! pointers, so it must not be copied, and [[inertia_destroy]] must be
        !! called when the solve ends.

        private

        logical, public  :: enabled = .false.     !! whether inertia control is in use (set by [[inertia_initialize]],
                                                  !! and unset if a factorization fails)
        integer, public  :: n_factor = 0          !! number of factorizations so far
        real(wp), public :: shift_last = 0.0_wp   !! the last shift that a correction ended with (`0` if none has
                                                  !! been needed yet)
        integer  :: n = 0                         !! number of variables
        integer  :: m = 0                         !! number of constraints
        logical  :: analysed = .false.            !! whether the sparsity pattern has been analysed
        logical  :: checked = .false.             !! whether `in_set`, `shift_checked`, and `ok_checked` hold the
                                                  !! outcome of a correction with the current matrices
        logical, dimension(:), allocatable :: in_set !! the working set of that correction `dimension(m+n)`
        real(wp) :: shift_checked = 0.0_wp        !! the shift it ended with
        logical  :: ok_checked = .false.          !! whether it left no negative curvature
#ifdef HAS_MUMPS
        type(dmumps_struc), pointer :: id => null() !! the MUMPS instance
#endif

        contains

        procedure, public :: initialize => inertia_initialize
        procedure, public :: new_matrices => inertia_new_matrices
        procedure, public :: correct => inertia_correct
        procedure, public :: raise => inertia_raise
        procedure, public :: destroy => inertia_destroy

    end type sqpopt_inertia_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  start the inertia control of a solve, for a problem with the sparsity
!  patterns of `hessian` (in exact mode) and `jac_irow`/`jac_icol`. This
!  creates the MUMPS instance, and sets the pattern of the KKT matrix (see
!  the module documentation): its diagonal, then the Hessian's elements,
!  then the Jacobian's, all in the lower triangle. `ok` is false, and
!  inertia control is left disabled, if the library was built without MUMPS
!  or MUMPS couldn't be started.

    subroutine inertia_initialize(me, hessian, m, jac_irow, jac_icol, ok)

    class(sqpopt_inertia_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(in)    :: hessian  !! the exact Hessian (for its sparsity pattern)
    integer,                    intent(in)    :: m        !! number of constraints
    integer, dimension(:),      intent(in)    :: jac_irow !! row indices of the Jacobian's nonzeros
    integer, dimension(:),      intent(in)    :: jac_icol !! column indices of the Jacobian's nonzeros
    logical,                    intent(out)   :: ok       !! whether inertia control is now enabled

#ifdef HAS_MUMPS
    integer :: n, nh, nj, k, alloc_stat
#endif

    call me%destroy()
    ok = .false.

#ifdef HAS_MUMPS
    n  = hessian%n
    nh = size(hessian%h_irow)
    nj = size(jac_irow)
    me%n = n
    me%m = m

    allocate(me%id, stat=alloc_stat)
    if (alloc_stat /= 0) then
        nullify(me%id)
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
    me%id%icntl(16)  = 1    ! one OpenMP thread
    me%id%icntl(16)  = 1    ! one OpenMP thread (see the module documentation)
    me%id%icntl(24)  = 1    ! detect null pivots, rather than fail on a singular matrix

    me%id%n   = n + m
    me%id%nnz = int(n + m + nh + nj, kind(me%id%nnz))
    allocate(me%id%irn(n+m+nh+nj), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%id%jcn(n+m+nh+nj), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%id%a(n+m+nh+nj), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%in_set(m+n), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call me%destroy()
        return
    end if
    do k = 1, n + m
        me%id%irn(k) = k
        me%id%jcn(k) = k
    end do
    do k = 1, nh
        me%id%irn(n+m+k) = max(hessian%h_irow(k), hessian%h_icol(k))
        me%id%jcn(n+m+k) = min(hessian%h_irow(k), hessian%h_icol(k))
    end do
    do k = 1, nj
        me%id%irn(n+m+nh+k) = n + jac_irow(k)
        me%id%jcn(n+m+nh+k) = jac_icol(k)
    end do

    me%enabled = .true.
    ok = .true.
#endif

    end subroutine inertia_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  note that the Hessian's or the Jacobian's values have changed (a new
!  major iteration), so the outcome of the last correction no longer holds.

    subroutine inertia_new_matrices(me)

    class(sqpopt_inertia_type), intent(inout) :: me

    me%checked = .false.

    end subroutine inertia_new_matrices
!*******************************************************************************

!*******************************************************************************
!>
!  increase the Hessian's shift (`hessian%shift`), if necessary, until
!  \( H + \delta I \) has no negative curvature on the null space of the
!  working set `status` (see the module documentation for the shifts that
!  are tried). The shift is never decreased: the search starts from its
!  current value.
!
!  `ok` is true if no negative curvature is left. It is false if inertia
!  control isn't enabled, if a factorization failed (inertia control is then
!  disabled for the rest of the solve), or if the largest shift,
!  `shift_max`, isn't enough. A call with the working set and the shift that
!  the last call ended with (and the same matrices, see
!  [[inertia_new_matrices]]) returns its outcome without factoring again.

    subroutine inertia_correct(me, hessian, jac, status, changed, ok, n_negative)

    class(sqpopt_inertia_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! the exact Hessian (its `shift` is updated)
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! the constraint Jacobian `dimension(m,n)`
    integer, dimension(:),      intent(in)    :: status  !! the working set: nonzero for each general row, then each
                                                         !! variable bound, that is in it `dimension(m+n)`
    logical,                    intent(out)   :: changed !! whether the shift was increased
    logical,                    intent(out)   :: ok      !! whether no negative curvature is left (see above)
    integer, optional,          intent(out)   :: n_negative !! the number of negative eigenvalues of
                                                            !! \( Z^THZ \) found with the shift on entry
                                                            !! (`-1` if nothing was factored)

    real(wp) :: shift_max
    integer  :: n_neg
    logical  :: first, factored

    changed = .false.
    ok      = .false.
    if (present(n_negative)) n_negative = -1
    if (.not. me%enabled) return

    if (me%checked) then
        if (hessian%shift == me%shift_checked .and. all((status /= 0) .eqv. me%in_set)) then
            ok = me%ok_checked
            return
        end if
    end if

    shift_max = hessian%shift_max*hessian%magnitude()
    first = .true.
    do
        call factor(me, hessian, jac, status, n_neg, factored)
        if (.not. factored) then
            ! (fall back on the matrix-free tests from here on)
            me%enabled = .false.
            return
        end if
        if (first .and. present(n_negative)) n_negative = max(n_neg - me%m, 0)
        if (n_neg <= me%m) then
            ok = .true.
            exit
        end if
        if (hessian%shift >= shift_max) exit
        call me%raise(hessian)
        changed = .true.
        first   = .false.
    end do
    if (ok .and. changed) me%shift_last = hessian%shift

    me%checked       = .true.
    me%in_set        = status /= 0
    me%shift_checked = hessian%shift
    me%ok_checked    = ok

    end subroutine inertia_correct
!*******************************************************************************

!*******************************************************************************
!>
!  increase the Hessian's shift to the next one to try (see the module
!  documentation): 8 times the current one, and at least `shift_min` times
!  the size of the Hessian, or a third of the last shift that was needed;
!  at most `shift_max` times the size of the Hessian.

    subroutine inertia_raise(me, hessian)

    class(sqpopt_inertia_type), intent(in)    :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! the exact Hessian (its `shift` is updated)

    real(wp) :: shift

    shift = hessian%shift_min*hessian%magnitude()
    if (me%shift_last > 0.0_wp) shift = max(shift, me%shift_last/3.0_wp)
    shift = max(shift, 8.0_wp*hessian%shift)
    hessian%shift = min(shift, hessian%shift_max*hessian%magnitude())

    end subroutine inertia_raise
!*******************************************************************************

!*******************************************************************************
!>
!  factor the KKT matrix of the working set `status`, with the Hessian's
!  current values and shift (see the module documentation), and return its
!  number of negative pivots. The sparsity pattern is analysed in the first
!  call, with that call's values. If MUMPS runs out of its estimated
!  workspace, the factorization is repeated with more.

    subroutine factor(me, hessian, jac, status, n_neg, ok)

    class(sqpopt_inertia_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(in)    :: hessian !! the exact Hessian
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! the constraint Jacobian `dimension(m,n)`
    integer, dimension(:),      intent(in)    :: status  !! the working set (see [[inertia_correct]]) `dimension(m+n)`
    integer,                    intent(out)   :: n_neg   !! the number of negative pivots
    logical,                    intent(out)   :: ok      !! whether the factorization succeeded

#ifdef HAS_MUMPS
    integer :: n, m, nh, k, i, j, attempt
    integer, parameter :: max_attempts = 6 !! (each doubles the workspace increase)
#endif

    n_neg = 0
    ok    = .false.

#ifdef HAS_MUMPS
    n  = me%n
    m  = me%m
    nh = size(hessian%h_val)
    associate (a => me%id%a)
        ! the diagonal: the shift, or 1 for a variable at a bound of the working
        ! set; and 0, or -1 for a constraint outside it:
        do j = 1, n
            a(j) = merge(1.0_real64, real(hessian%shift, real64), status(m+j) /= 0)
        end do
        do i = 1, m
            a(n+i) = merge(0.0_real64, -1.0_real64, status(i) /= 0)
        end do
        do k = 1, nh
            i = hessian%h_irow(k)
            j = hessian%h_icol(k)
            a(n+m+k) = merge(0.0_real64, real(hessian%h_val(k), real64), status(m+i) /= 0 .or. status(m+j) /= 0)
        end do
        do k = 1, jac%nnz
            a(n+m+nh+k) = merge(real(jac%val(k), real64), 0.0_real64, &
                                status(jac%irow(k)) /= 0 .and. status(m+jac%icol(k)) == 0)
        end do
    end associate

    do attempt = 1, max_attempts
        me%id%job = merge(2, 4, me%analysed)   ! (4: analyse, then factor)
        call dmumps(me%id)
        if (me%id%infog(1) /= -8 .and. me%id%infog(1) /= -9) exit
        me%id%icntl(14) = 2*max(me%id%icntl(14), 20)   ! more workspace (a percentage of the estimate)
    end do
    me%n_factor = me%n_factor + 1
    if (me%id%infog(1) < 0) return
    me%analysed = .true.
    n_neg = me%id%infog(12)
    ok    = .true.
#endif

    end subroutine factor
!*******************************************************************************

!*******************************************************************************
!>
!  end the inertia control of a solve: free the MUMPS instance and the
!  arrays, and return to the initial (disabled) state. It does nothing if
!  inertia control was never started.

    subroutine inertia_destroy(me)

    class(sqpopt_inertia_type), intent(inout) :: me

#ifdef HAS_MUMPS
    if (associated(me%id)) then
        if (associated(me%id%irn)) deallocate(me%id%irn)
        if (associated(me%id%jcn)) deallocate(me%id%jcn)
        if (associated(me%id%a))   deallocate(me%id%a)
        me%id%job = -2
        call dmumps(me%id)
        deallocate(me%id)
    end if
#endif
    if (allocated(me%in_set)) deallocate(me%in_set)
    me%enabled    = .false.
    me%n_factor   = 0
    me%shift_last = 0.0_wp
    me%n          = 0
    me%m          = 0
    me%analysed   = .false.
    me%checked    = .false.
    me%shift_checked = 0.0_wp
    me%ok_checked    = .false.

    end subroutine inertia_destroy
!*******************************************************************************

    end module sqpopt_inertia_module
!*******************************************************************************
