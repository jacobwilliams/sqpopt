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
!     stationarity condition \( J^T \lambda \approx g \), solved with `LSQR`.
!  2. A *normal step* \( p_n \) is computed as the minimum-norm solution of
!     \( J p_n = \text{viol} \), where `viol` is the linearized constraint
!     violation, again solved with `LSQR`.
!  3. A *tangential step* is computed as the (approximate) unconstrained
!     quasi-Newton direction \( v = H^{-1} g \) (matrix-free two-loop
!     recursion, see [[sqpopt_hessian_module]]), **projected onto the null
!     space of \( J \)** so that it does not reintroduce constraint
!     infeasibility: \( p_t = v - J^T z \), where `z` minimizes
!     \( \lVert J^T z - v \rVert_2 \) (again solved with `LSQR`).
!  4. \( p = p_n - p_t \) is clipped component-wise so that \( x+p \)
!     respects the variable bounds.
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

        integer :: max_iter           = 0                     !! maximum number of iterations allowed for the QP solver
        integer :: linear_solver_mode = sqpopt_linsolve_lusol  !! sparse linear solver used for the KKT system

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

    integer :: n, m, k, istop
    real(wp), dimension(size(c)) :: viol, z
    real(wp), dimension(size(g)) :: v, jtz
    type(lsqr_solver_ez) :: lsqr

    n = size(g)
    m = size(c)

    if (m == 0) then
        ! no general constraints: the tangential step is just the
        ! unconstrained quasi-Newton step (no null-space projection needed):
        call hessian%inverse_vector_product(g, v)
        p = -v
    else

        ! (1) Lagrange multiplier estimate: least-squares solve of J^T*lambda = g,
        !     i.e. A*lambda = g with A = J^T (an n x m matrix, stored by
        !     transposing the (irow,icol) pattern of the sparse Jacobian):
        call lsqr%initialize(n, m, jac%val, jac%icol, jac%irow)
        call lsqr%solve(g, 0.0_wp, lambda, istop)

        ! (2) normal step: minimum-norm solution of J*p_n = viol, where `viol`
        !     is the change in c(x) needed to satisfy the linearized bounds:
        do k = 1, m
            viol(k) = min(max(c(k), c_lb(k)), c_ub(k)) - c(k)
        end do
        call lsqr%initialize(m, n, jac%val, jac%irow, jac%icol)
        call lsqr%solve(viol, 0.0_wp, p, istop)

        ! (3) tangential step: the (approximate) unconstrained quasi-Newton step
        !     v = H^{-1} g, projected onto the null space of J so that it does
        !     not reintroduce constraint infeasibility. The null-space
        !     projection of v is `v - J^T z`, where z minimizes ||J^T z - v||_2
        !     (again solved with LSQR -- the same least-squares problem shape
        !     as the multiplier estimate in step 1, just with a different
        !     right-hand side):
        call hessian%inverse_vector_product(g, v)
        call lsqr%initialize(n, m, jac%val, jac%icol, jac%irow)
        call lsqr%solve(v, 0.0_wp, z, istop)
        call sparse_matvec_transpose(jac, z, jtz)
        p = p - (v - jtz)

    end if

    ! (4) clip the combined step so that x+p respects the variable bounds:
    do k = 1, n
        p(k) = min(max(x(k)+p(k), x_lb(k)), x_ub(k)) - x(k)
    end do

    istat = sqpopt_success

    end subroutine solve_qp_subproblem
!*******************************************************************************

    end module sqpopt_qp_solver_module
!*******************************************************************************
