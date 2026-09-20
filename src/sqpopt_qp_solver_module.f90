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

    module sqpopt_qp_solver_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    type, public :: sqpopt_qp_solver_type
        !! workspace and options for the QP subproblem solver.

        integer :: max_iter = 0  !! maximum number of iterations allowed for the QP solver

        contains

        procedure, public :: solve => solve_qp_subproblem

    end type sqpopt_qp_solver_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda`.

    subroutine solve_qp_subproblem(me, h, g, jac, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    class(sqpopt_qp_solver_type), intent(inout) :: me
    real(wp), dimension(:,:), intent(in)  :: h       !! Hessian approximation `dimension(n,n)`
    real(wp), dimension(:),   intent(in)  :: g       !! objective gradient `dimension(n)`
    real(wp), dimension(:,:), intent(in)  :: jac     !! constraint Jacobian `dimension(m,n)`
    real(wp), dimension(:),   intent(in)  :: c       !! current constraint values `dimension(m)`
    real(wp), dimension(:),   intent(in)  :: x_lb    !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),   intent(in)  :: x_ub    !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),   intent(in)  :: c_lb    !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),   intent(in)  :: c_ub    !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:),   intent(out) :: p       !! computed search direction `dimension(n)`
    real(wp), dimension(:),   intent(out) :: lambda  !! Lagrange multipliers for the linearized constraints `dimension(m)`
    integer,                   intent(out) :: istat  !! status code (see [[sqpopt_types_module]])

    ! TODO: implement

    end subroutine solve_qp_subproblem
!*******************************************************************************

    end module sqpopt_qp_solver_module
!*******************************************************************************
