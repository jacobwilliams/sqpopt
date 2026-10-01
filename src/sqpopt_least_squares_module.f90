!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Direct minimum-norm solves with rows of the constraint Jacobian
!  (`options%direct_least_squares`), by a sparse factorization (see
!  [[sqpopt_kkt_module]]), so it is only available in a build with MUMPS.
!
!  The Gauss-Newton restoration step ([[restoration_step]]) and the
!  second-order correction ([[soc_step]]) both need the minimum-norm
!  solution `d` of
!
!  $$ J_S\, d = r $$
!
!  for a set `S` of the Jacobian's rows (in the least-squares sense, if the
!  rows are dependent or the system is inconsistent). Without this module
!  they use the iterative solver `LSQR`, and take no step if it stops at its
!  iteration limit. Here the solution comes from one factorization and one
!  solve of
!
!  $$ \begin{bmatrix} I & J_S^T \\ J_S & -\epsilon I \end{bmatrix}
!     \begin{bmatrix} d \\ y \end{bmatrix} =
!     \begin{bmatrix} 0 \\ r \end{bmatrix} $$
!
!  whose solution is \( d = J_S^T (J_S J_S^T + \epsilon I)^{-1} r \). The
!  small regularization \( \epsilon \) (relative to the square of the
!  Jacobian's largest element) keeps the matrix nonsingular when the rows
!  are dependent; as \( \epsilon \to 0 \), `d` tends to the minimum-norm
!  least-squares solution.
!
!  The same matrix, with another right-hand side, gives the least-squares
!  *multiplier estimate* \( \lambda = \arg\min \lVert g - J_S^T \lambda
!  \rVert \) of [[multiplier_estimate]], which [[sqpopt_iterate]] uses
!  when the QP's multipliers can't be trusted. That routine works in every
!  build: it uses `LSQR` unless this module's direct solver is enabled.

    module sqpopt_least_squares_module

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_sparse_matrix, sqpopt_all_finite
    use sqpopt_kkt_module,   only: sqpopt_kkt_type
    use lsqr_module,         only: lsqr_solver_ez

    implicit none

    private

    public :: multiplier_estimate

    real(wp), parameter :: reg = 1.0e-8_wp !! the regularization \( \epsilon \), relative to the square of the
                                           !! Jacobian's largest element (it limits the relative accuracy of
                                           !! the solution of a consistent system to about this)

    type, public :: sqpopt_least_squares_type
        !! the direct least-squares solver of one solve (see the module
        !! documentation). It holds a sparse solver, so it must not be copied,
        !! and [[least_squares_destroy]] must be called when the solve ends.

        private

        logical, public :: enabled = .false.  !! whether it can be used (set by [[least_squares_initialize]], and
                                              !! unset if a factorization fails)
        type(sqpopt_kkt_type), public :: kkt  !! the matrix above, as the KKT matrix of the Hessian `I` (see
                                              !! [[sqpopt_kkt_module]]; public for its solver's counts and time)
        integer :: n = 0  !! number of variables
        integer :: m = 0  !! number of constraints

        contains

        procedure, public :: initialize   => least_squares_initialize
        procedure, public :: new_matrices => least_squares_new_matrices
        procedure, public :: min_norm     => least_squares_min_norm
        procedure, public :: multipliers  => least_squares_multipliers
        procedure, public :: destroy      => least_squares_destroy

    end type sqpopt_least_squares_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  start the solver for a problem with `n` variables, `m` constraints, and
!  the Jacobian sparsity pattern `jac_irow`/`jac_icol`. `ok` is false, and
!  the solver is left disabled, if the library was built without MUMPS or
!  the sparse solver couldn't be started.

    subroutine least_squares_initialize(me, n, m, jac_irow, jac_icol, ok, threads)

    class(sqpopt_least_squares_type), intent(inout) :: me
    integer,               intent(in)  :: n        !! number of variables
    integer,               intent(in)  :: m        !! number of constraints
    integer, dimension(:), intent(in)  :: jac_irow !! row indices of the Jacobian's nonzeros
    integer, dimension(:), intent(in)  :: jac_icol !! column indices of the Jacobian's nonzeros
    logical,               intent(out) :: ok       !! whether the solver is now enabled
    integer, optional,     intent(in)  :: threads  !! number of OpenMP threads of the sparse solver (see
                                                   !! [[symmetric_solver_initialize]])

    me%n = n
    me%m = m
    call me%kkt%initialize(n, m, jac_irow, jac_icol, ok, threads=threads)
    me%enabled = ok

    end subroutine least_squares_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  note that the Jacobian's values have changed (a new major iteration).

    subroutine least_squares_new_matrices(me)

    class(sqpopt_least_squares_type), intent(inout) :: me

    call me%kkt%new_matrices()

    end subroutine least_squares_new_matrices
!*******************************************************************************

!*******************************************************************************
!>
!  the (regularized) minimum-norm solution `d` of \( J_S d = r \) (see the
!  module documentation), for the rows `S` of `jac` selected by `rows`.
!  `ok` is false if the solver isn't enabled, or the factorization or the
!  solve failed (the solver is then disabled for the rest of the solve):
!  the caller then falls back on its iterative method.

    subroutine least_squares_min_norm(me, jac, rows, r, d, ok)

    class(sqpopt_least_squares_type), intent(inout) :: me
    type(sqpopt_sparse_matrix), intent(in)  :: jac  !! the constraint Jacobian `dimension(m,n)`
    logical,  dimension(:),     intent(in)  :: rows !! whether each row is in `S` `dimension(m)`
    real(wp), dimension(:),     intent(in)  :: r    !! the right-hand side (used for the rows in `S`) `dimension(m)`
    real(wp), dimension(:),     intent(out) :: d    !! the solution `dimension(n)`
    logical,                    intent(out) :: ok   !! whether `d` was computed

    integer,  dimension(me%m + me%n) :: status
    real(wp), dimension(me%n + me%m) :: v
    real(wp) :: jmax

    d  = 0.0_wp
    ok = .false.
    if (.not. me%enabled) return

    ! (the working set is the rows of `S`, and no variable is fixed)
    status = 0
    where (rows) status(1:me%m) = 1
    jmax = 1.0_wp
    if (jac%nnz > 0) jmax = max(1.0_wp, maxval(abs(jac%val(1:jac%nnz))))

    call me%kkt%factor_identity(1.0_wp, jac, status, reg*jmax**2, ok)
    if (ok) then
        v(1:me%n) = 0.0_wp
        v(me%n+1:) = merge(r, 0.0_wp, rows)
        call me%kkt%solve(v, ok)
        if (ok) ok = sqpopt_all_finite(v(1:me%n))
    end if
    if (.not. ok) then
        me%enabled = .false.
        return
    end if
    d = v(1:me%n)

    end subroutine least_squares_min_norm
!*******************************************************************************

!*******************************************************************************
!>
!  the (regularized) least-squares multipliers of the rows `S` of `jac`
!  selected by `rows`, in the variables selected by `free`:
!  \( \lambda = (J_S J_S^T + \epsilon I)^{-1} J_S\, g \), from the system
!  of the module documentation with the right-hand side \( (g, 0) \).
!  `lambda` is only set for the rows in `S`. `ok` is false if the solver
!  isn't enabled, or the factorization or the solve failed (the solver is
!  then disabled for the rest of the solve).

    subroutine least_squares_multipliers(me, jac, rows, free, g, lambda, ok)

    class(sqpopt_least_squares_type), intent(inout) :: me
    type(sqpopt_sparse_matrix), intent(in)    :: jac    !! the constraint Jacobian `dimension(m,n)`
    logical,  dimension(:),     intent(in)    :: rows   !! whether each row is in `S` `dimension(m)`
    logical,  dimension(:),     intent(in)    :: free   !! whether each variable is free (not at a bound) `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g      !! the objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(inout) :: lambda !! the multipliers `dimension(m)` (set for the rows in `S`)
    logical,                    intent(out)   :: ok     !! whether they were computed

    integer,  dimension(me%m + me%n) :: status
    real(wp), dimension(me%n + me%m) :: v
    real(wp) :: jmax

    ok = .false.
    if (.not. me%enabled) return

    ! (the working set is the rows of `S`, and the variables that are not free)
    status = 0
    where (rows) status(1:me%m) = 1
    where (.not. free) status(me%m+1:) = 1
    jmax = 1.0_wp
    if (jac%nnz > 0) jmax = max(1.0_wp, maxval(abs(jac%val(1:jac%nnz))))

    call me%kkt%factor_identity(1.0_wp, jac, status, reg*jmax**2, ok)
    if (ok) then
        v(1:me%n)  = merge(g, 0.0_wp, free)
        v(me%n+1:) = 0.0_wp
        call me%kkt%solve(v, ok)
        if (ok) ok = sqpopt_all_finite(v(me%n+1:))
    end if
    if (.not. ok) then
        me%enabled = .false.
        return
    end if
    where (rows) lambda = v(me%n+1:)

    end subroutine least_squares_multipliers
!*******************************************************************************

!*******************************************************************************
!>
!  a first-order estimate of the constraint multipliers at a point: the
!  least-squares solution of the stationarity condition
!
!  $$ \lambda = \arg\min \lVert (g - J_S^T \lambda)_{\text{free}} \rVert_2 $$
!
!  over the rows `S` of `jac` selected by `rows`, in the variables selected
!  by `free` (those not at a bound, whose bound multipliers are not zero).
!  Unlike the multipliers of a QP subproblem, it doesn't depend on the QP's
!  Hessian. `lambda` is only set for the rows in `S`, and only if `ok`.
!
!  It is computed by the direct solver `least_squares`, if that is given
!  and enabled (see `options%direct_least_squares`), and otherwise by the
!  iterative solver `LSQR`; `ok` is false if that stops at its iteration
!  limit.

    subroutine multiplier_estimate(jac, rows, free, g, lambda, ok, least_squares)

    type(sqpopt_sparse_matrix), intent(in)    :: jac    !! the constraint Jacobian `dimension(m,n)`
    logical,  dimension(:),     intent(in)    :: rows   !! whether each row is in `S` `dimension(m)`
    logical,  dimension(:),     intent(in)    :: free   !! whether each variable is free (not at a bound) `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g      !! the objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(inout) :: lambda !! the multipliers `dimension(m)` (set for the rows in `S`)
    logical,                    intent(out)   :: ok     !! whether they were computed
    type(sqpopt_least_squares_type), optional, intent(inout) :: least_squares !! the direct least-squares solver

    integer, parameter :: lsqr_itnlim_stop = 5 !! `LSQR`'s `istop` for "iteration limit reached"
    type(lsqr_solver_ez) :: lsqr
    integer,  dimension(:), allocatable :: irow, icol, col_of_row
    real(wp), dimension(:), allocatable :: val, lam
    integer :: m, n, m_s, nnz_s, i, k, istop

    m = size(rows)
    n = size(free)
    ok = .false.
    m_s = count(rows)
    if (m_s == 0) return

    if (present(least_squares)) then
        call least_squares%multipliers(jac, rows, free, g, lambda, ok)
        if (ok) return
    end if

    ! the transposed sub-Jacobian (free variables x rows of S), for LSQR:
    allocate(col_of_row(m))
    col_of_row = 0
    k = 0
    do i = 1, m
        if (rows(i)) then
            k = k + 1
            col_of_row(i) = k
        end if
    end do
    allocate(irow(jac%nnz), icol(jac%nnz), val(jac%nnz), lam(m_s))
    nnz_s = 0
    do k = 1, jac%nnz
        if (rows(jac%irow(k)) .and. free(jac%icol(k))) then
            nnz_s = nnz_s + 1
            irow(nnz_s) = jac%icol(k)
            icol(nnz_s) = col_of_row(jac%irow(k))
            val(nnz_s)  = jac%val(k)
        end if
    end do

    call lsqr%initialize(n, m_s, val(1:nnz_s), irow(1:nnz_s), icol(1:nnz_s), itnlim=4*(m_s+n)+10)
    call lsqr%solve(merge(g, 0.0_wp, free), 0.0_wp, lam, istop)
    if (istop == lsqr_itnlim_stop .or. .not. sqpopt_all_finite(lam)) return
    do i = 1, m
        if (rows(i)) lambda(i) = lam(col_of_row(i))
    end do
    ok = .true.

    end subroutine multiplier_estimate
!*******************************************************************************

!*******************************************************************************
!>
!  free the solver, and return to the initial (disabled) state. It does
!  nothing if the solver was never started.

    subroutine least_squares_destroy(me)

    class(sqpopt_least_squares_type), intent(inout) :: me

    call me%kkt%destroy()
    me%enabled = .false.
    me%n = 0
    me%m = 0

    end subroutine least_squares_destroy
!*******************************************************************************

    end module sqpopt_least_squares_module
!*******************************************************************************
