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

    module sqpopt_least_squares_module

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_sparse_matrix, sqpopt_all_finite
    use sqpopt_kkt_module,   only: sqpopt_kkt_type

    implicit none

    private

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

    subroutine least_squares_initialize(me, n, m, jac_irow, jac_icol, ok)

    class(sqpopt_least_squares_type), intent(inout) :: me
    integer,               intent(in)  :: n        !! number of variables
    integer,               intent(in)  :: m        !! number of constraints
    integer, dimension(:), intent(in)  :: jac_irow !! row indices of the Jacobian's nonzeros
    integer, dimension(:), intent(in)  :: jac_icol !! column indices of the Jacobian's nonzeros
    logical,               intent(out) :: ok       !! whether the solver is now enabled

    me%n = n
    me%m = m
    call me%kkt%initialize(n, m, jac_irow, jac_icol, ok)
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
