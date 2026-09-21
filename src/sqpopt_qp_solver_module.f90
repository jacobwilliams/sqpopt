!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Solves the quadratic programming (QP) subproblem generated at each
!  major SQP iteration:
!
!  $$ \min_{p} \; \frac{1}{2} p^T H p + g^T p $$
!
!  subject to the linearized constraints and variable bounds:
!
!  $$ c_l - c(x) \le J(x) p \le c_u - c(x) $$
!  $$ x_l - x \le p \le x_u - x $$
!
!  Neither \( H \) nor \( J \) is ever formed as a dense matrix: \( H \)
!  is a matrix-free limited-memory operator (see [[sqpopt_hessian_module]])
!  and \( J \) is stored as sparse COO triplets (see [[sqpopt_types_module]]).
!
!  **Composite-step algorithm** (the default, chosen so that
!  neither `H` nor `J` ever needs to be factorized or formed as a dense
!  `n x n`/`m x n` array):
!
!  1. A Lagrange multiplier estimate is obtained from the least-squares
!     stationarity condition \( J_A^T \lambda_A \approx g \), solved with
!     `LSQR`, where \( J_A \) is the sub-matrix of rows of \( J \) belonging
!     to the *active set* (equality constraints, plus inequality
!     constraints currently at or beyond one of their bounds -- see
!     `active_tol`). Multipliers for inactive constraints are fixed at zero
!     (a simple stand-in for full complementarity, since this method does not
!     maintain a proper working set across iterations).
!  2. A *normal step* \( p_n \) is computed as the minimum-norm solution of
!     \( J p_n = \text{viol} \), where `viol` is the linearized constraint
!     violation (using the *full* `J`; rows with zero violation contribute
!     nothing, so inactive constraints do not need to be filtered here),
!     again solved with `LSQR`.
!  3. A *tangential step* is computed as the (approximate) unconstrained
!     quasi-Newton direction \( v = H^{-1} g \) (matrix-free two-loop
!     recursion, see [[sqpopt_hessian_module]]), **projected onto the null
!     space of \( J_A \)** (the same active-set sub-matrix as step 1) so
!     that it does not reintroduce infeasibility in the active constraints:
!     \( p_t = v - J_A^T z \), where `z` minimizes \( \lVert J_A^T z - v
!     \rVert_2 \) (again solved with `LSQR`).
!  4. \( p = p_n - p_t \) is rescaled if \( \lVert p \rVert_2 \) exceeds
!     `max_step` (a simple trust-region-style safeguard against the
!     composite step occasionally overshooting), then adjusted so that
!     \( x+p \) respects the variable bounds, using the strategy selected
!     by `bound_enforcement` (`sqpopt_bounds_scalar`, the default, clips
!     only the violating components; `sqpopt_bounds_vector` rescales the
!     whole step uniformly instead, preserving its direction).
!
!  This composite-step method (the default, `sqpopt_qp_composite`) does
!  not enforce the linearized general-constraint bounds exactly, relying
!  on the outer major SQP iterations to converge to feasibility instead;
!  two other modes are available for problems that need the linearized
!  constraints solved exactly: `sqpopt_qp_dense` (a dense active-set QP,
!  see [[sqpopt_qp_dense_module]]) and `sqpopt_qp_reduced_hessian` (a
!  sparse/matrix-free active-set QP, see
!  [[sqpopt_qp_reduced_hessian_module]]). Both of those enforce variable
!  bounds exactly as part of the QP solve itself, so `bound_enforcement`
!  does not apply to them.

    module sqpopt_qp_solver_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_success
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_linalg_module,  only: sparse_matvec_transpose
    use sqpopt_qp_dense_module, only: sqpopt_dense_qp_type
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_reduced_hessian_qp_type
    use lsqr_module,           only: lsqr_solver_ez

    implicit none

    private

    integer, parameter, public :: sqpopt_qp_composite       = 1  !! composite-step heuristic (default, see module docs)
    integer, parameter, public :: sqpopt_qp_dense           = 2  !! opt-in dense active-set QP solver (see [[sqpopt_qp_dense_module]])
    integer, parameter, public :: sqpopt_qp_reduced_hessian = 3  !! opt-in sparse (projected-CG) active-set QP solver (see [[sqpopt_qp_reduced_hessian_module]])

    integer, parameter, public :: sqpopt_bounds_vector = 1  !! rescale the *entire* step `p` by the same factor so that
                                                             !! `x+p` just touches the first bound it would otherwise
                                                             !! violate, preserving `p`'s direction exactly
    integer, parameter, public :: sqpopt_bounds_scalar = 2  !! (default) clip only the violating components of `x+p`
                                                             !! to their bound; the other components of `p` are left
                                                             !! unchanged

    type, public :: sqpopt_qp_solver_type
        !! workspace and options for the QP subproblem solver.

        integer  :: mode                = sqpopt_qp_composite  !! which QP algorithm to use (see the `sqpopt_qp_*` constants)
        real(wp) :: max_step           = 2.0_wp                 !! trust-region-style cap on \( \lVert p \rVert_2 \);
                                                                 !! the step is rescaled if it is exceeded (safeguards
                                                                 !! against the composite step occasionally
                                                                 !! overshooting -- see [[sqpopt_qp_solver_module]])
        real(wp) :: active_tol         = 1.0e-6_wp              !! an inequality constraint is considered part of the
                                                                 !! active set if it is within `active_tol` of (or beyond)
                                                                 !! one of its bounds (equality constraints are always active)
        integer  :: bound_enforcement  = sqpopt_bounds_scalar   !! how `mode==sqpopt_qp_composite` enforces the variable
                                                                 !! bounds `x_lb<=x+p<=x_ub` on its computed step (see the
                                                                 !! `sqpopt_bounds_*` constants); not used by `sqpopt_qp_dense`/
                                                                 !! `sqpopt_qp_reduced_hessian`, which enforce bounds exactly
                                                                 !! as part of the QP solve itself
        type(sqpopt_dense_qp_type)           :: dense_qp    !! the dense QP solver (used only when `mode==sqpopt_qp_dense`)
        type(sqpopt_reduced_hessian_qp_type) :: sparse_qp   !! the sparse QP solver (used only when `mode==sqpopt_qp_reduced_hessian`)

        contains

        procedure, public :: solve => solve_qp_subproblem

    end type sqpopt_qp_solver_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda`, dispatching to the
!  algorithm selected by `me%mode`: the composite-step heuristic
!  (default, see the module-level documentation) or one of the opt-in
!  active-set solvers (see [[sqpopt_qp_dense_module]],
!  [[sqpopt_qp_reduced_hessian_module]]).

    subroutine solve_qp_subproblem(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    class(sqpopt_qp_solver_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! matrix-free Hessian approximation (never a dense `n x n` matrix)
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! sparse constraint Jacobian, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x       !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g       !! objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c       !! current constraint values `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: x_lb    !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: x_ub    !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c_lb    !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: c_ub    !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:),     intent(out)   :: p       !! computed search direction `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: lambda  !! Lagrange multipliers for the linearized constraints `dimension(m)`
    integer,                    intent(out)   :: istat   !! status code (see [[sqpopt_types_module]])

    select case (me%mode)
    case (sqpopt_qp_dense)
        call me%dense_qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)
        ! the dense solver enforces bounds/constraints exactly, but still apply the
        ! same trust-region cap as the composite step, for a consistent step-size safeguard:
        if (norm2(p) > me%max_step) p = p*(me%max_step/norm2(p))
    case (sqpopt_qp_reduced_hessian)
        call me%sparse_qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)
        if (norm2(p) > me%max_step) p = p*(me%max_step/norm2(p))
    case default
        call solve_composite_step(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)
    end select

    end subroutine solve_qp_subproblem
!*******************************************************************************

!*******************************************************************************
!>
!  the composite-step algorithm (see the module-level documentation).

    subroutine solve_composite_step(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    class(sqpopt_qp_solver_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! matrix-free Hessian approximation (never a dense `n x n` matrix)
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! sparse constraint Jacobian, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x       !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g       !! objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c       !! current constraint values `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: x_lb    !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: x_ub    !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c_lb    !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: c_ub    !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:),     intent(out)   :: p       !! computed search direction `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: lambda  !! Lagrange multipliers for the linearized constraints `dimension(m)`
    integer,                    intent(out)   :: istat   !! status code (see [[sqpopt_types_module]])

    integer :: n, m, k, istop, m_active
    real(wp), dimension(size(c)) :: viol
    real(wp), dimension(size(g)) :: v, jtz
    type(lsqr_solver_ez) :: lsqr
    logical, dimension(size(c)) :: active
    integer, dimension(:), allocatable :: active_rows
    real(wp), dimension(:), allocatable :: lambda_active, z
    type(sqpopt_sparse_matrix) :: jac_a

    n = size(g)
    m = size(c)

    if (m == 0) then
        ! no general constraints: the tangential step is just the
        ! unconstrained quasi-Newton step (no null-space projection needed):
        call hessian%inverse_vector_product(g, v)
        p = -v
    else

        ! determine the active set: equality constraints, plus inequality
        ! constraints currently at (or beyond) one of their bounds:
        do k = 1, m
            active(k) = (c_ub(k)-c_lb(k) <= me%active_tol) .or. &
                        (c(k)-c_lb(k) <= me%active_tol) .or. (c_ub(k)-c(k) <= me%active_tol)
        end do
        active_rows = pack([(k, k=1,m)], active)
        m_active = size(active_rows)
        lambda = 0.0_wp

        if (m_active == 0) then

            ! no active constraints: same as the unconstrained case:
            call hessian%inverse_vector_product(g, v)
            p = -v

        else

            call select_active_rows(jac, active_rows, jac_a)
            allocate(lambda_active(m_active), z(m_active))

            ! (1) Lagrange multiplier estimate: least-squares solve of J_A^T*lambda = g,
            !     i.e. A*lambda = g with A = J_A^T (an n x m_active matrix, stored by
            !     transposing the (irow,icol) pattern of the active-set sub-Jacobian):
            call lsqr%initialize(n, m_active, jac_a%val, jac_a%icol, jac_a%irow)
            call lsqr%solve(g, 0.0_wp, lambda_active, istop)
            lambda(active_rows) = lambda_active

            ! (2) normal step: minimum-norm solution of J*p_n = viol, where `viol`
            !     is the change in c(x) needed to satisfy the linearized bounds
            !     (using the full J; inactive rows have zero violation so they
            !     do not need to be filtered out here):
            do k = 1, m
                viol(k) = min(max(c(k), c_lb(k)), c_ub(k)) - c(k)
            end do
            call lsqr%initialize(m, n, jac%val, jac%irow, jac%icol)
            call lsqr%solve(viol, 0.0_wp, p, istop)

            ! (3) tangential step: the (approximate) unconstrained quasi-Newton step
            !     v = H^{-1} g, projected onto the null space of J_A so that it does
            !     not reintroduce infeasibility in the active constraints. The
            !     null-space projection of v is `v - J_A^T z`, where z minimizes
            !     ||J_A^T z - v||_2 (again solved with LSQR -- the same
            !     least-squares problem shape as the multiplier estimate in step 1,
            !     just with a different right-hand side):
            call hessian%inverse_vector_product(g, v)
            call lsqr%initialize(n, m_active, jac_a%val, jac_a%icol, jac_a%irow)
            call lsqr%solve(v, 0.0_wp, z, istop)
            call sparse_matvec_transpose(jac_a, z, jtz)
            p = p - (v - jtz)

        end if

    end if

    ! (4) trust-region-style safeguard: rescale the step if it is
    !     unreasonably large (the composite step is only an
    !     approximate QP solution and can occasionally overshoot):
    if (norm2(p) > me%max_step) p = p*(me%max_step/norm2(p))

    ! (5) enforce the variable bounds x_lb<=x+p<=x_ub on the (possibly
    !     rescaled) step, using the strategy selected by `me%bound_enforcement`:
    select case (me%bound_enforcement)
    case (sqpopt_bounds_vector)
        call rescale_step_to_bounds(x, x_lb, x_ub, p)
    case default ! sqpopt_bounds_scalar
        do k = 1, n
            p(k) = min(max(x(k)+p(k), x_lb(k)), x_ub(k)) - x(k)
        end do
    end select

    istat = sqpopt_success

    end subroutine solve_composite_step
!*******************************************************************************

!*******************************************************************************
!>
!  `sqpopt_bounds_vector` bound enforcement: find the largest \( \alpha \in
!  [0,1] \) such that \( x+\alpha p \) satisfies every variable bound, then
!  rescale the whole step `p := \alpha p`. Unlike component-wise clipping,
!  this preserves the step's direction exactly (only its length changes).

    subroutine rescale_step_to_bounds(x, x_lb, x_ub, p)

    real(wp), dimension(:), intent(in)    :: x, x_lb, x_ub
    real(wp), dimension(:), intent(inout) :: p

    real(wp) :: alpha, alpha_k
    integer  :: k

    alpha = 1.0_wp
    do k = 1, size(x)
        if (p(k) > 0.0_wp .and. x(k)+p(k) > x_ub(k)) then
            alpha_k = (x_ub(k) - x(k))/p(k)
            alpha = min(alpha, alpha_k)
        else if (p(k) < 0.0_wp .and. x(k)+p(k) < x_lb(k)) then
            alpha_k = (x_lb(k) - x(k))/p(k)
            alpha = min(alpha, alpha_k)
        end if
    end do
    alpha = max(alpha, 0.0_wp)

    p = alpha*p

    end subroutine rescale_step_to_bounds
!*******************************************************************************

!*******************************************************************************
!>
!  build the sub-matrix of `jac` containing only the rows listed in
!  `active_rows` (renumbered `1..size(active_rows)`), used to restrict the
!  multiplier estimate and null-space projection to the active set.

    subroutine select_active_rows(jac, active_rows, jac_a)

    type(sqpopt_sparse_matrix), intent(in)  :: jac         !! full constraint Jacobian, `dimension(m,n)`
    integer, dimension(:),      intent(in)  :: active_rows !! original row indices to keep, `dimension(m_active)`
    type(sqpopt_sparse_matrix), intent(out) :: jac_a       !! resulting sub-matrix, `dimension(size(active_rows),n)`

    integer, dimension(:), allocatable :: row_map  !! original row -> new row index (0 if not active)
    integer :: k, idx, nnz_a

    allocate(row_map(jac%nrows))
    row_map = 0
    do k = 1, size(active_rows)
        row_map(active_rows(k)) = k
    end do

    nnz_a = count(row_map(jac%irow(1:jac%nnz)) > 0)
    jac_a%nrows = size(active_rows)
    jac_a%ncols = jac%ncols
    jac_a%nnz   = nnz_a
    allocate(jac_a%irow(nnz_a), jac_a%icol(nnz_a), jac_a%val(nnz_a))

    idx = 0
    do k = 1, jac%nnz
        if (row_map(jac%irow(k)) > 0) then
            idx = idx + 1
            jac_a%irow(idx) = row_map(jac%irow(k))
            jac_a%icol(idx) = jac%icol(k)
            jac_a%val(idx)  = jac%val(k)
        end if
    end do

    end subroutine select_active_rows
!*******************************************************************************

    end module sqpopt_qp_solver_module
!*******************************************************************************
