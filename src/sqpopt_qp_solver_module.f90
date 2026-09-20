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
!  The resulting QP is solved via a sparse active-set method whose KKT
!  system is factorized/solved using the sparse solver requested by
!  `linear_solver_mode` (see [[sqpopt_linalg_module]]).

    module sqpopt_qp_solver_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_linalg_module,  only: sqpopt_linsolve_lusol

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
!  the associated Lagrange multipliers `lambda`.

    subroutine solve_qp_subproblem(me, hessian, jac, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    class(sqpopt_qp_solver_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! matrix-free Hessian approximation (never a dense `n x n` matrix)
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! sparse constraint Jacobian, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: g       !! objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c       !! current constraint values `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: x_lb    !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: x_ub    !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c_lb    !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: c_ub    !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:),     intent(out)   :: p       !! computed search direction `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: lambda  !! Lagrange multipliers for the linearized constraints `dimension(m)`
    integer,                    intent(out)   :: istat   !! status code (see [[sqpopt_types_module]])

    ! TODO: implement

    end subroutine solve_qp_subproblem
!*******************************************************************************

    end module sqpopt_qp_solver_module
!*******************************************************************************
