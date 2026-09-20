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
    use sqpopt_types_module,      only: sqpopt_sparse_matrix, sqpopt_success
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type, sqpopt_hessian_sr1
    use sqpopt_hessian_module,    only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type
    use sqpopt_linalg_module,     only: sparse_matvec_transpose
    use sqpopt_convergence_module, only: check_convergence

    implicit none

    private

    public :: sqpopt_iterate

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  perform one major SQP iteration, updating `x` and `lambda` in place.
!  `x_prev`/`gl_prev` hold the previous point / Lagrangian gradient used
!  to form the quasi-Newton `(s,y)` pair; they should be passed in
!  *unallocated* before the first call (the Hessian update is skipped on
!  the first iteration, since no previous point is yet available).

    subroutine sqpopt_iterate(problem, options, hessian, qp_solver, linesearch, &
                               x, lambda, x_prev, gl_prev, converged, istat)

    type(sqpopt_problem_type),    intent(inout) :: problem     !! problem definition
    type(sqpopt_options_type),    intent(in)    :: options     !! solver options
    type(sqpopt_hessian_type),    intent(inout) :: hessian     !! Hessian of the Lagrangian approximation
    type(sqpopt_qp_solver_type),  intent(inout) :: qp_solver   !! QP subproblem solver
    type(sqpopt_linesearch_type), intent(inout) :: linesearch  !! merit function / line search
    real(wp), dimension(:), intent(inout) :: x       !! current point, updated on exit `dimension(n)`
    real(wp), dimension(:), intent(inout) :: lambda  !! current Lagrange multipliers, updated on exit `dimension(m)`
    real(wp), dimension(:), allocatable, intent(inout) :: x_prev  !! previous point (unallocated before the 1st call)
    real(wp), dimension(:), allocatable, intent(inout) :: gl_prev !! previous Lagrangian gradient (unallocated before the 1st call)
    logical,                 intent(out)   :: converged !! true if `x` (on entry) already satisfies the convergence criteria
    integer,                 intent(out)   :: istat     !! status code (see [[sqpopt_types_module]])

    real(wp) :: f
    real(wp), dimension(problem%n) :: g, gl, p
    real(wp), dimension(problem%m) :: c, new_lambda
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: alpha
    integer :: qp_istat, ls_istat, conv_istat

    ! evaluate the problem functions and the sparse Jacobian at the current point:
    call problem%eval_f(x, f)
    call problem%eval_g(x, g)
    call problem%eval_c(x, c)
    jac%nrows = problem%m
    jac%ncols = problem%n
    jac%nnz   = problem%jac_nnz
    jac%irow  = problem%jac_irow
    jac%icol  = problem%jac_icol
    allocate(jac%val(problem%jac_nnz))
    call problem%eval_jac(x, jac%val)

    ! check convergence at the current point before taking a step:
    call check_convergence(x, g, jac, c, problem%x_lb, problem%x_ub, problem%c_lb, problem%c_ub, &
                            lambda, options%ktol, options%ctol, converged, conv_istat)
    if (converged) then
        istat = sqpopt_success
        return
    end if

    ! Lagrangian gradient at the current point (used for the next quasi-Newton update):
    gl = g
    block
        real(wp), dimension(problem%n) :: jtlam
        call sparse_matvec_transpose(jac, lambda, jtlam)
        gl = g - jtlam
    end block

    ! update the quasi-Newton Hessian approximation using the previous step
    ! (skipped on the very first iteration, since there is no previous point):
    if (allocated(x_prev)) then
        if (options%hessian_mode == sqpopt_hessian_sr1) then
            call hessian%update_sr1(x - x_prev, gl - gl_prev)
        else
            ! `sqpopt_hessian_exact` is not yet supported in this v1 driver;
            ! falls back to the BFGS update (see PLAN.md for future work).
            call hessian%update_bfgs(x - x_prev, gl - gl_prev)
        end if
    end if

    ! solve the linearized QP subproblem for the search direction and multipliers:
    call qp_solver%solve(hessian, jac, x, g, c, problem%x_lb, problem%x_ub, &
                          problem%c_lb, problem%c_ub, p, new_lambda, qp_istat)

    ! line search along `p` to (approximately) minimize the merit function:
    call linesearch%search(problem%eval_f, problem%eval_c, x, p, problem%c_lb, problem%c_ub, alpha, ls_istat)

    ! save the current point/gradient for the next quasi-Newton update:
    x_prev  = x
    gl_prev = gl

    ! update the point and multipliers:
    x      = x + alpha*p
    lambda = new_lambda

    istat = sqpopt_success

    end subroutine sqpopt_iterate
!*******************************************************************************

    end module sqpopt_iterate_module
!*******************************************************************************
