!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The core SQP major iteration: evaluates the problem functions, updates
!  the (limited-memory, matrix-free) Hessian approximation, solves the
!  sparse QP subproblem for the search direction, and performs a line
!  search to update the current point. No dense `n x n` or `m x n`
!  matrix is ever formed.

    module sqpopt_iterate_module

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type
    use sqpopt_hessian_module,    only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type

    implicit none

    private

    public :: sqpopt_iterate

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  perform one major SQP iteration, updating `x` and `lambda` in place.

    subroutine sqpopt_iterate(problem, options, hessian, qp_solver, linesearch, x, lambda, istat)

    type(sqpopt_problem_type),    intent(inout) :: problem     !! problem definition
    type(sqpopt_options_type),    intent(in)    :: options     !! solver options
    type(sqpopt_hessian_type),    intent(inout) :: hessian     !! Hessian of the Lagrangian approximation
    type(sqpopt_qp_solver_type),  intent(inout) :: qp_solver   !! QP subproblem solver
    type(sqpopt_linesearch_type), intent(inout) :: linesearch  !! merit function / line search
    real(wp), dimension(:), intent(inout) :: x       !! current point, updated on exit `dimension(n)`
    real(wp), dimension(:), intent(inout) :: lambda  !! current Lagrange multipliers, updated on exit `dimension(m)`
    integer,                 intent(out)   :: istat   !! status code (see [[sqpopt_types_module]])

    ! TODO: implement

    end subroutine sqpopt_iterate
!*******************************************************************************

    end module sqpopt_iterate_module
!*******************************************************************************
