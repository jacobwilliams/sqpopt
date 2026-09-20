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
    use sqpopt_linalg_module,     only: sparse_matvec_transpose, sparse_matvec
    use sqpopt_convergence_module, only: check_convergence
    use lsqr_module,              only: lsqr_solver_ez

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

    ! update the merit function's penalty parameter so that it dominates the
    ! current multiplier estimates (as in slsqp): this is required for the
    ! l1 exact penalty function's minimizer to coincide with the true
    ! constrained optimum (Han/Powell); without it, the merit function can
    ! prefer a "compromise" infeasible point over the true solution:
    if (size(new_lambda) > 0) linesearch%penalty = max(linesearch%penalty, maxval(abs(new_lambda)) + 1.0_wp)

    ! safeguard (as in slsqp): if `p` is not a descent direction for the
    ! merit function (can happen since the v1 composite step is only an
    ! approximate QP solution, so it lacks the usual guarantee that the
    ! *optimal* QP solution is a descent direction), reset the Hessian
    ! approximation to the identity and recompute `p` once from scratch:
    if (dot_product(g,p) - linesearch%penalty*constraint_violation(c, problem%c_lb, problem%c_ub) >= 0.0_wp) then
        call hessian%reset()
        call qp_solver%solve(hessian, jac, x, g, c, problem%x_lb, problem%x_ub, &
                              problem%c_lb, problem%c_ub, p, new_lambda, qp_istat)
    end if

    ! second-order correction (SOC): for strongly nonlinear constraints, the
    ! full step `p` can be rejected by the merit function even when it is a
    ! genuinely good step, because the *linearized* constraint prediction
    ! differs from the true (nonlinear) constraint value at `x+p` (the
    ! "Maratos effect"). Correct for this by solving for an additional small
    ! step that accounts for the true constraint residual at `x+p`, and use
    ! the corrected step if it has a better merit function value:
    if (problem%m > 0) then
        call second_order_correction(problem, linesearch, jac, x, c, p)
    end if

    ! line search along `p` to (approximately) minimize the merit function:
    call linesearch%search(problem%eval_f, problem%eval_c, x, p, f, g, c, problem%c_lb, problem%c_ub, alpha, ls_istat)

    ! save the current point/gradient for the next quasi-Newton update:
    x_prev  = x
    gl_prev = gl

    ! update the point and multipliers:
    x      = x + alpha*p
    lambda = new_lambda

    istat = sqpopt_success

    end subroutine sqpopt_iterate
!*******************************************************************************

!*******************************************************************************
!>
!  the \( \ell_1 \) constraint violation measure \( \lVert \max(c_l-c,0,c-c_u)
!  \rVert_1 \), used by the descent-direction safeguard above.

    pure function constraint_violation(c, c_lb, c_ub) result(v)

    real(wp), dimension(:), intent(in) :: c, c_lb, c_ub
    real(wp) :: v

    v = sum(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))

    end function constraint_violation
!*******************************************************************************

!*******************************************************************************
!>
!  second-order correction: refine `p` using the true (nonlinear)
!  constraint residual at `x+p` rather than its linear prediction, and
!  replace `p` with the corrected step if it improves the merit function
!  (a standard remedy for the Maratos effect near nonlinear constraints).

    subroutine second_order_correction(problem, linesearch, jac, x, c, p)

    type(sqpopt_problem_type),    intent(inout) :: problem
    type(sqpopt_linesearch_type), intent(inout) :: linesearch
    type(sqpopt_sparse_matrix),   intent(in)    :: jac
    real(wp), dimension(:),       intent(in)    :: x
    real(wp), dimension(:),       intent(in)    :: c
    real(wp), dimension(:),       intent(inout) :: p

    real(wp), dimension(size(c)) :: jp, c_p, resid, c_soc
    real(wp), dimension(size(p)) :: p_corr, p_soc
    real(wp) :: f_p, f_soc, phi_p, phi_soc
    type(lsqr_solver_ez) :: lsqr
    integer :: istop

    call sparse_matvec(jac, p, jp)
    call problem%eval_c(x+p, c_p)
    resid = c_p - (c+jp)  !! nonlinear residual left uncorrected by the linear model

    call lsqr%initialize(problem%m, problem%n, jac%val, jac%irow, jac%icol)
    call lsqr%solve(-resid, 0.0_wp, p_corr, istop)
    p_soc = p + p_corr

    call problem%eval_f(x+p, f_p)
    call linesearch%eval_merit(f_p, c_p, problem%c_lb, problem%c_ub, phi_p)

    call problem%eval_f(x+p_soc, f_soc)
    call problem%eval_c(x+p_soc, c_soc)
    call linesearch%eval_merit(f_soc, c_soc, problem%c_lb, problem%c_ub, phi_soc)

    if (phi_soc < phi_p) p = p_soc

    end subroutine second_order_correction
!*******************************************************************************

    end module sqpopt_iterate_module
!*******************************************************************************
