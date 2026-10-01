!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The KKT matrix of a QP subproblem's working set, factored by the sparse
!  direct solver of [[sqpopt_symmetric_solver_module]] (so it is only
!  available in a build with MUMPS):
!
!  $$ K = \begin{bmatrix} H & J_a^T \\ J_a & -\epsilon I \end{bmatrix} $$
!
!  where \( H \) is the Hessian in the *free* variables (those not at a
!  bound of the working set), \( J_a \) holds the working set's rows of the
!  Jacobian in those variables, and \( \epsilon \ge 0 \) is an optional
!  regularization (zero, unless a singular matrix is being regularized).
!
!  It answers two questions about the QP restricted to that working set
!  (a *face* of the QP):
!
!  * **Is it convex there?** By Sylvester's law of inertia, \( K \) has one
!    negative eigenvalue for each (independent) row of \( J_a \), plus one
!    for each negative eigenvalue of the reduced Hessian \( Z^THZ \). So
!    `n_negative`, the number of negative eigenvalues beyond the number of
!    rows, counts the directions of negative curvature on the face. This is
!    what the inertia control uses (see [[sqpopt_inertia_module]]).
!  * **Where is its minimizer?** One solve with the factors gives the step
!    to the minimizer of the face, and its multipliers (see
!    [[sqpopt_qp_direct_module]]).
!
!  **A fixed pattern.** So that the sparsity pattern is the same for every
!  working set (it is analysed only once per solve), \( K \) always has
!  order `n+m`: a variable at a bound of the working set has its row and
!  column replaced by those of the identity, and a constraint outside it by
!  those of minus the identity. The unknowns of a solve are then all `n`
!  step components followed by all `m` multipliers (see [[kkt_solve]]).
!
!  **One factorization at a time.** The last factorization is kept, with
!  the working set it was made for. Asking for the same one again (while
!  the matrices haven't changed, see [[kkt_new_matrices]]) costs nothing,
!  so the inertia control and the direct QP step share their
!  factorizations.
!
!  The Hessian is the solver's ([[kkt_factor]]), or a multiple of the
!  identity ([[kkt_factor_identity]], for the least-squares problems of
!  [[sqpopt_least_squares_module]]).
!
!  **Quasi-Newton Hessians.** The user's exact Hessian is sparse, and goes
!  into \( K \) as it is. A limited-memory quasi-Newton Hessian is not
!  sparse, but it is a multiple of the identity plus a matrix of low rank
!  `r` (see [[sqpopt_hessian_module]]),
!  \( B = \theta I + \sigma\, U M^{-1} U^T \). So
!  \( K = K_0 + \sigma\, E M^{-1} E^T \), where \( K_0 \) is the sparse
!  matrix with \( \theta I \) for the Hessian, and \( E \) is \( U \)
!  in the free variables (zero elsewhere). Only \( K_0 \) is factored, and:
!
!  * a solve uses the Sherman-Morrison-Woodbury formula,
!    \( K^{-1} = K_0^{-1} - K_0^{-1} E\, T^{-1} E^T K_0^{-1} \), with the
!    small dense matrix \( T = \sigma M + E^T K_0^{-1} E \): two solves
!    with \( K_0 \) each, after `r` solves to form \( T \) when the
!    matrix is factored (as IPOPT does for its limited-memory option);
!  * the inertia follows from the two ways of eliminating a block of
!    \( \begin{bmatrix} K_0 & E \\ E^T & -\sigma M \end{bmatrix} \)
!    (the Haynsworth inertia additivity formula): the number of negative
!    eigenvalues of \( K \) is that of \( K_0 \), plus the number of
!    positive eigenvalues of \( T \), minus that of \( \sigma M \).

    module sqpopt_kkt_module

    use sqpopt_kinds,                   only: wp => sqpopt_module_wp
    use sqpopt_types_module,            only: sqpopt_sparse_matrix
    use sqpopt_hessian_module,          only: sqpopt_hessian_type
    use sqpopt_symmetric_solver_module, only: sqpopt_symmetric_solver_type
    use sqpopt_dense_linalg_module,     only: dense_lu_factor, dense_lu_solve, dense_symmetric_inertia

    implicit none

    private

    type, public :: sqpopt_kkt_type
        !! the KKT matrix of a working set, and its factorization (see the
        !! module documentation). It holds the sparse solver, so it must not
        !! be copied, and [[kkt_destroy]] must be called when the solve ends.

        private

        logical, public :: enabled = .false.  !! whether it can be used (set by [[kkt_initialize]], and unset if
                                              !! a factorization fails)
        integer, public :: n_negative = 0     !! of the last factorization: the number of directions of negative
                                              !! curvature on the working set's face (the negative eigenvalues of
                                              !! \( K \) beyond `m`)
        logical, public :: singular = .false. !! of the last factorization: whether \( K \) is singular (the
                                              !! working set's rows are dependent, or the face has a direction
                                              !! of zero curvature)
        type(sqpopt_symmetric_solver_type), public :: solver !! the sparse solver (with its counts and its time)

        integer :: n = 0   !! number of variables
        integer :: m = 0   !! number of constraints
        integer, dimension(:), allocatable :: h_irow !! row indices of the Hessian's pattern
        integer, dimension(:), allocatable :: h_icol !! column indices of the Hessian's pattern
        real(wp), dimension(:), allocatable :: val   !! the values of \( K \): its diagonal, then the Hessian's
                                                     !! elements, then the Jacobian's

        ! the factorization that is kept (see the module documentation):
        logical  :: current = .false.                !! whether there is one, for the current matrices
        logical, dimension(:), allocatable :: in_set !! its working set `dimension(m+n)`
        real(wp) :: diag = 0.0_wp                    !! what was added to its Hessian's diagonal
        real(wp) :: reg  = 0.0_wp                    !! its regularization \( \epsilon \)
        integer  :: rank = 0                         !! the rank `r` of its quasi-Newton Hessian's low-rank part
                                                     !! (`0` if it has none)
        real(wp), dimension(:,:), allocatable :: t_lu  !! the LU factors of its matrix \( T \) `dimension(r,r)`
        integer,  dimension(:),   allocatable :: t_piv !! the row pivots of that factorization `dimension(r)`

        contains

        procedure, public :: initialize      => kkt_initialize
        procedure, public :: new_matrices    => kkt_new_matrices
        procedure, public :: factor          => kkt_factor
        procedure, public :: factor_identity => kkt_factor_identity
        procedure, public :: solve           => kkt_solve
        procedure, public :: destroy         => kkt_destroy

    end type sqpopt_kkt_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  set up the KKT matrix of a problem with `n` variables, `m` constraints,
!  and the Jacobian sparsity pattern `jac_irow`/`jac_icol`. The Hessian's
!  pattern is `hess_irow`/`hess_icol` (each off-diagonal element once, as
!  for [[set_hessian_sparsity]]); without it, the Hessian is diagonal. The
!  pattern of \( K \) is its whole diagonal, then the Hessian's elements,
!  then the Jacobian's. `ok` is false, and the matrix is left disabled, if
!  the library was built without MUMPS or the solver couldn't be started.

    subroutine kkt_initialize(me, n, m, jac_irow, jac_icol, ok, hess_irow, hess_icol, threads)

    class(sqpopt_kkt_type), intent(inout) :: me
    integer,               intent(in)  :: n         !! number of variables
    integer,               intent(in)  :: m         !! number of constraints
    integer, dimension(:), intent(in)  :: jac_irow  !! row indices of the Jacobian's nonzeros
    integer, dimension(:), intent(in)  :: jac_icol  !! column indices of the Jacobian's nonzeros
    logical,               intent(out) :: ok        !! whether the matrix can now be factored
    integer, dimension(:), intent(in), optional :: hess_irow !! row indices of the Hessian's nonzeros
    integer, dimension(:), intent(in), optional :: hess_icol !! column indices of the Hessian's nonzeros
    integer,               intent(in), optional :: threads   !! number of OpenMP threads of the sparse solver (see
                                                             !! [[symmetric_solver_initialize]])

    integer :: k
    logical :: oom

    call me%destroy()
    me%n = n
    me%m = m
    if (present(hess_irow) .and. present(hess_icol)) then
        me%h_irow = hess_irow
        me%h_icol = hess_icol
    else
        allocate(me%h_irow(0), me%h_icol(0))
    end if

    call me%solver%initialize(n + m, &
                              [(k, k=1, n+m), max(me%h_irow, me%h_icol), n + jac_irow], &
                              [(k, k=1, n+m), min(me%h_irow, me%h_icol), jac_icol], ok, threads=threads)
    if (.not. ok) then
        ! (keep the reason, for the caller)
        oom = me%solver%out_of_memory
        call me%destroy()
        me%solver%out_of_memory = oom
        return
    end if
    allocate(me%val(n + m + size(me%h_irow) + size(jac_irow)), me%in_set(m + n))
    me%enabled = .true.

    end subroutine kkt_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  note that the Hessian's or the Jacobian's values have changed (a new
!  major iteration), so the factorization that is kept no longer applies.

    subroutine kkt_new_matrices(me)

    class(sqpopt_kkt_type), intent(inout) :: me

    me%current = .false.

    end subroutine kkt_new_matrices
!*******************************************************************************

!*******************************************************************************
!>
!  factor the KKT matrix of the working set `status`, with the Hessian
!  `hessian` (the exact one, with its current values and shift, or a
!  quasi-Newton one, with its current pairs: see the module documentation).
!  On success, `n_negative` and `singular` describe the working set's face,
!  and [[kkt_solve]] can be called. If this is the factorization that is
!  kept, nothing is computed. The matrix must have been initialized with
!  the Hessian's sparsity pattern in the first case, and without one in the
!  second.
!
!  `reg` is the regularization \( \epsilon \) of the module documentation
!  (default 0). A small positive value makes the matrix of a working set
!  with dependent rows nonsingular: the solve then satisfies those rows in
!  a (heavily weighted) least-squares sense instead of exactly.
!
!  `ok` is false if the matrix isn't enabled, or the factorization failed
!  (it is then disabled for the rest of the solve).

    subroutine kkt_factor(me, hessian, jac, status, ok, reg)

    class(sqpopt_kkt_type),     intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! the Hessian (see [[sqpopt_hessian_module]])
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! the constraint Jacobian `dimension(m,n)`
    integer, dimension(:),      intent(in)    :: status  !! the working set: nonzero for each general row, then each
                                                         !! variable bound, that is in it `dimension(m+n)`
    logical,                    intent(out)   :: ok      !! whether the factorization is available
    real(wp), optional,         intent(in)    :: reg     !! the regularization \( \epsilon \ge 0 \) (default 0)

    real(wp) :: eps
    real(wp), dimension(0) :: no_values
    logical :: kept

    eps = 0.0_wp
    if (present(reg)) eps = reg
    if (hessian%exact) then
        call factor_values(me, hessian%shift, hessian%h_val, jac, status, eps, ok)
    else
        call factor_values(me, 1.0_wp/hessian%gamma + hessian%shift, no_values, jac, status, eps, ok, kept)
        if (ok .and. .not. kept) call low_rank_update(me, hessian, ok)
    end if

    end subroutine kkt_factor
!*******************************************************************************

!*******************************************************************************
!>
!  factor the KKT matrix of the working set `status` with the Hessian
!  `scale` \( \times I \) and the regularization `reg` (\( \epsilon \)),
!  as [[kkt_factor]] does for the exact Hessian. The matrix must have been
!  initialized without a Hessian pattern.

    subroutine kkt_factor_identity(me, scale, jac, status, reg, ok)

    class(sqpopt_kkt_type),     intent(inout) :: me
    real(wp),                   intent(in)    :: scale  !! the Hessian is `scale` times the identity
    type(sqpopt_sparse_matrix), intent(in)    :: jac    !! the constraint Jacobian `dimension(m,n)`
    integer, dimension(:),      intent(in)    :: status !! the working set (see [[kkt_factor]]) `dimension(m+n)`
    real(wp),                   intent(in)    :: reg    !! the regularization \( \epsilon \ge 0 \)
    logical,                    intent(out)   :: ok     !! whether the factorization is available

    real(wp), dimension(0) :: no_values

    call factor_values(me, scale, no_values, jac, status, reg, ok)

    end subroutine kkt_factor_identity
!*******************************************************************************

!*******************************************************************************
!>
!  assemble and factor \( K \) (see the module documentation) for the
!  Hessian with the values `h_val` (in the pattern given to
!  [[kkt_initialize]]) plus `diag` on its diagonal, unless that is the
!  factorization that is kept (`kept`).

    subroutine factor_values(me, diag, h_val, jac, status, reg, ok, kept)

    class(sqpopt_kkt_type),     intent(inout) :: me
    real(wp),                   intent(in)    :: diag   !! added to the Hessian's diagonal
    real(wp), dimension(:),     intent(in)    :: h_val  !! the Hessian's values, in its pattern's order
    type(sqpopt_sparse_matrix), intent(in)    :: jac    !! the constraint Jacobian `dimension(m,n)`
    integer, dimension(:),      intent(in)    :: status !! the working set (see [[kkt_factor]]) `dimension(m+n)`
    real(wp),                   intent(in)    :: reg    !! the regularization \( \epsilon \ge 0 \)
    logical,                    intent(out)   :: ok     !! whether the factorization is available
    logical, optional,          intent(out)   :: kept   !! whether it is the one that was kept (nothing was computed)

    integer :: n, m, nh, k, i, j

    ok = .false.
    if (present(kept)) kept = .false.
    if (.not. me%enabled) return

    if (me%current) then
        if (diag == me%diag .and. reg == me%reg .and. all((status /= 0) .eqv. me%in_set)) then
            ok = .true.
            if (present(kept)) kept = .true.
            return
        end if
    end if

    n  = me%n
    m  = me%m
    nh = size(me%h_irow)

    ! the diagonal: `diag`, or 1 for a variable at a bound of the working
    ! set; and `-reg`, or -1 for a constraint outside it:
    do j = 1, n
        me%val(j) = merge(1.0_wp, diag, status(m+j) /= 0)
    end do
    do i = 1, m
        me%val(n+i) = merge(-reg, -1.0_wp, status(i) /= 0)
    end do
    ! the Hessian's elements between free variables:
    do k = 1, nh
        i = me%h_irow(k)
        j = me%h_icol(k)
        me%val(n+m+k) = merge(0.0_wp, h_val(k), status(m+i) /= 0 .or. status(m+j) /= 0)
    end do
    ! the Jacobian's elements of the working set's rows, in the free variables:
    do k = 1, jac%nnz
        me%val(n+m+nh+k) = merge(jac%val(k), 0.0_wp, status(jac%irow(k)) /= 0 .and. status(m+jac%icol(k)) == 0)
    end do

    me%current = .false.
    call me%solver%factor(me%val, ok)
    if (.not. ok) then
        ! (fall back on the matrix-free methods from here on)
        me%enabled = .false.
        return
    end if
    me%n_negative = max(me%solver%n_negative - m, 0)
    me%singular   = me%solver%n_null > 0
    me%current    = .true.
    me%in_set     = status /= 0
    me%diag       = diag
    me%reg        = reg
    me%rank       = 0

    end subroutine factor_values
!*******************************************************************************

!*******************************************************************************
!>
!  complete the factorization of \( K_0 \) just made for the quasi-Newton
!  Hessian `hessian` (see the module documentation): form the matrix
!  \( T \) (one solve with \( K_0 \) for each of its `r` columns) and
!  its LU factorization, and correct `n_negative` and `singular` for the
!  low-rank part. (A BFGS matrix is positive definite, so its face has no
!  negative curvature, and the inertias are only computed for SR1.)

    subroutine low_rank_update(me, hessian, ok)

    class(sqpopt_kkt_type),    intent(inout) :: me
    type(sqpopt_hessian_type), intent(inout) :: hessian !! the quasi-Newton Hessian
    logical,                   intent(out)   :: ok      !! whether the solves succeeded

    real(wp), dimension(:,:), allocatable :: t, sm
    real(wp), dimension(me%n + me%m) :: v
    real(wp) :: sigma
    integer  :: r, j, n, pos_t, pos_m, n_neg, n_zero
    logical  :: nonsingular

    ok = .true.
    n  = me%n
    r  = hessian%low_rank_size()
    me%rank = r
    if (r == 0) return

    allocate(t(r,r), sm(r,r))
    call hessian%low_rank_middle(sm, sigma)
    sm = sigma*sm

    ! T = sigma*M + E^T K0^{-1} E, a column at a time:
    do j = 1, r
        call hessian%low_rank_column(j, v(1:n))
        where (me%in_set(me%m+1:)) v(1:n) = 0.0_wp
        v(n+1:) = 0.0_wp
        call me%solver%solve(v, ok, refine=.false.)
        if (.not. ok) then
            me%enabled = .false.
            me%current = .false.
            return
        end if
        where (me%in_set(me%m+1:)) v(1:n) = 0.0_wp
        call hessian%low_rank_multiply_transpose(v(1:n), t(:,j))
    end do
    t = sm + 0.5_wp*(t + transpose(t))

    ! the inertia of K, from those of K0, T, and sigma*M:
    if (sigma > 0.0_wp) then
        call dense_symmetric_inertia(sm, pos_m, n_neg, n_zero)
        call dense_symmetric_inertia(t, pos_t, n_neg, n_zero)
        me%n_negative = max(me%solver%n_negative + pos_t - pos_m - me%m, 0)
    end if

    if (allocated(me%t_lu)) deallocate(me%t_lu, me%t_piv)
    allocate(me%t_piv(r))
    call move_alloc(t, me%t_lu)
    call dense_lu_factor(me%t_lu, me%t_piv, nonsingular)
    if (.not. nonsingular) me%singular = .true.

    end subroutine low_rank_update
!*******************************************************************************

!*******************************************************************************
!>
!  solve \( K v = r \) with the factorization that is kept: `v` holds `r`
!  on entry and the solution on exit. Its first `n` elements belong to the
!  variables, and the other `m` to the constraints. Because of the fixed
!  pattern (see the module documentation), the solution's element is the
!  right-hand side's own for a variable at a bound of the working set, and
!  minus that for a constraint outside it.
!
!  `hessian` must be given if the factorization was made for a
!  quasi-Newton Hessian (the same one): its low-rank part is then applied
!  by the Sherman-Morrison-Woodbury formula (see the module documentation).
!
!  `ok` is false if there is no factorization, or no finite solution.

    subroutine kkt_solve(me, v, ok, hessian)

    class(sqpopt_kkt_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: v  !! the right-hand side, overwritten by the solution `dimension(n+m)`
    logical,                intent(out)   :: ok !! whether a solution was found
    type(sqpopt_hessian_type), optional, intent(in) :: hessian !! the quasi-Newton Hessian that was factored

    real(wp), dimension(:), allocatable :: w, e
    integer :: n

    ok = .false.
    if (.not. (me%enabled .and. me%current)) return
    call me%solver%solve(v, ok)
    if (.not. ok .or. me%rank == 0) return

    ! v = v - K0^{-1} E T^{-1} E^T v:
    ok = .false.
    if (.not. present(hessian)) return
    n = me%n
    allocate(w(me%rank), e(n + me%m))
    e(1:n) = merge(0.0_wp, v(1:n), me%in_set(me%m+1:))
    call hessian%low_rank_multiply_transpose(e(1:n), w)
    call dense_lu_solve(me%t_lu, me%t_piv, w)
    call hessian%low_rank_multiply(w, e(1:n))
    where (me%in_set(me%m+1:)) e(1:n) = 0.0_wp
    e(n+1:) = 0.0_wp
    call me%solver%solve(e, ok)
    if (ok) v = v - e

    end subroutine kkt_solve
!*******************************************************************************

!*******************************************************************************
!>
!  free the solver and the arrays, and return to the initial (disabled)
!  state. It does nothing if the matrix was never initialized.

    subroutine kkt_destroy(me)

    class(sqpopt_kkt_type), intent(inout) :: me

    call me%solver%destroy()
    if (allocated(me%h_irow)) deallocate(me%h_irow)
    if (allocated(me%h_icol)) deallocate(me%h_icol)
    if (allocated(me%val))    deallocate(me%val)
    if (allocated(me%in_set)) deallocate(me%in_set)
    if (allocated(me%t_lu))   deallocate(me%t_lu)
    if (allocated(me%t_piv))  deallocate(me%t_piv)
    me%rank       = 0
    me%enabled    = .false.
    me%n_negative = 0
    me%singular   = .false.
    me%n          = 0
    me%m          = 0
    me%current    = .false.
    me%diag       = 0.0_wp
    me%reg        = 0.0_wp

    end subroutine kkt_destroy
!*******************************************************************************

    end module sqpopt_kkt_module
!*******************************************************************************
