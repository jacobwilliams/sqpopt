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
!  **v1 algorithm** (a simplified composite-step method, chosen so that
!  neither `H` nor `J` ever needs to be factorized or formed as a dense
!  `n x n`/`m x n` array):
!
!  1. A Lagrange multiplier estimate is obtained from the least-squares
!     stationarity condition \( J_A^T \lambda_A \approx g \), solved with
!     `LSQR`, where \( J_A \) is the sub-matrix of rows of \( J \) belonging
!     to the *active set* (equality constraints, plus inequality
!     constraints currently at or beyond one of their bounds -- see
!     `active_tol`). Multipliers for inactive constraints are fixed at zero
!     (a simple stand-in for full complementarity, since v1 does not
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
!     composite step occasionally overshooting), then clipped
!     component-wise so that \( x+p \) respects the variable bounds.
!
!  This is a deliberate v1 simplification of a full active-set QP solve
!  (it does not enforce the linearized general-constraint bounds exactly,
!  relying on the outer major SQP iterations to converge to feasibility);
!  a more rigorous active-set/interior-point QP solver is a natural future
!  enhancement (see `PLAN.md`).

    module sqpopt_qp_solver_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_success
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_linalg_module,  only: sqpopt_linsolve_lusol, sparse_matvec_transpose
    use lsqr_module,           only: lsqr_solver_ez

    implicit none

    private

    type, public :: sqpopt_qp_solver_type
        !! workspace and options for the QP subproblem solver.

        integer  :: max_iter           = 0                     !! maximum number of iterations allowed for the QP solver
        integer  :: linear_solver_mode = sqpopt_linsolve_lusol  !! sparse linear solver used for the KKT system
        real(wp) :: max_step           = 2.0_wp                 !! trust-region-style cap on \( \lVert p \rVert_2 \);
                                                                 !! the step is rescaled if it is exceeded (safeguards
                                                                 !! against the v1 composite step occasionally
                                                                 !! overshooting -- see [[sqpopt_qp_solver_module]])
        real(wp) :: active_tol         = 1.0e-6_wp              !! an inequality constraint is considered part of the
                                                                 !! active set if it is within `active_tol` of (or beyond)
                                                                 !! one of its bounds (equality constraints are always active)

        contains

        procedure, public :: solve => solve_qp_subproblem

    end type sqpopt_qp_solver_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda` (see the module-level
!  documentation for the v1 algorithm used).

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
    !     unreasonably large (the v1 composite step is only an
    !     approximate QP solution and can occasionally overshoot):
    if (norm2(p) > me%max_step) p = p*(me%max_step/norm2(p))

    ! (5) clip the (possibly rescaled) step so that x+p respects the variable bounds:
    do k = 1, n
        p(k) = min(max(x(k)+p(k), x_lb(k)), x_ub(k)) - x(k)
    end do

    istat = sqpopt_success

    end subroutine solve_qp_subproblem
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
