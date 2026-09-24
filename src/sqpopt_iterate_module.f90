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
    use sqpopt_types_module,      only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_user_requested_stop, sqpopt_report_func
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type, sqpopt_hessian_sr1
    use sqpopt_hessian_module,    only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_filter
    use sqpopt_linalg_module,     only: sparse_matvec_transpose
    use sqpopt_convergence_module, only: check_convergence
    use sqpopt_soc_module,        only: second_order_correction
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type

    implicit none

    private

    public :: sqpopt_iterate

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  perform one major SQP iteration, updating `x` and `lambda` in place.
!  `x_prev`/`gl_prev`/`f_prev` hold the previous point / Lagrangian
!  gradient / objective value, used to form the quasi-Newton `(s,y)` pair
!  and to check `options%ftol`/`xtol`; they should be passed in
!  *unallocated* before the first call (the Hessian update and the
!  stalled-progress convergence test are both skipped on the first
!  iteration, since no previous point is yet available).

    subroutine sqpopt_iterate(problem, options, hessian, qp_solver, linesearch, trust_region, &
                               x, lambda, x_prev, gl_prev, f_prev, iter, report, converged, istat)

    type(sqpopt_problem_type),    intent(inout) :: problem     !! problem definition
    type(sqpopt_options_type),    intent(in)    :: options     !! solver options
    type(sqpopt_hessian_type),    intent(inout) :: hessian     !! Hessian of the Lagrangian approximation
    type(sqpopt_qp_solver_type),  intent(inout) :: qp_solver   !! QP subproblem solver
    type(sqpopt_linesearch_type), intent(inout) :: linesearch  !! merit function / line search
    type(sqpopt_trust_region_type), intent(inout) :: trust_region !! trust-region globalization (used instead of
                                                                   !! `linesearch` when `trust_region%enabled`)
    real(wp), dimension(:), intent(inout) :: x       !! current point, updated on exit `dimension(n)`
    real(wp), dimension(:), intent(inout) :: lambda  !! current Lagrange multipliers, updated on exit `dimension(m)`
    real(wp), dimension(:), allocatable, intent(inout) :: x_prev  !! previous point (unallocated before the 1st call)
    real(wp), dimension(:), allocatable, intent(inout) :: gl_prev !! previous Lagrangian gradient (unallocated before the 1st call)
    real(wp),               allocatable, intent(inout) :: f_prev  !! previous objective value (unallocated before the 1st call)
    integer,                intent(in)    :: iter      !! major iteration number (starts at 1), passed to `report`
    procedure(sqpopt_report_func), optional, pointer :: report !! optional user progress-reporting callback (see [[sqpopt_types_module]])
    logical,                 intent(out)   :: converged !! true if `x` (on entry) already satisfies the convergence criteria
    integer,                 intent(out)   :: istat     !! status code (see [[sqpopt_types_module]])

    real(wp) :: f !! current objective function value
    real(wp), dimension(problem%n) :: g, gl, p, x_new
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

    ! report progress on the current iterate, if the user has supplied a
    ! callback, before doing any further work this iteration -- this
    ! reports every major iterate, including the initial guess (iter=1,
    ! before any step has been taken) and the final, converged point:
    if (present(report)) then
        if (associated(report)) then
            block
                logical :: user_stop
                call report(iter, x, f, c, lambda, user_stop)
                if (user_stop) then
                    istat     = sqpopt_user_requested_stop
                    converged = .false.
                    return
                end if
            end block
        end if
    end if

    ! check convergence at the current point before taking a step (the
    ! stalled-progress test also needs `f_prev`/`x_prev`, so is only
    ! available from the 2nd iteration onward):
    if (allocated(x_prev)) then
        call check_convergence(x, g, jac, c, problem%x_lb, problem%x_ub, problem%c_lb, problem%c_ub, &
                                lambda, options%ktol, options%ctol, converged, conv_istat, &
                                f=f, f_prev=f_prev, x_prev=x_prev, ftol=options%ftol, xtol=options%xtol)
    else
        call check_convergence(x, g, jac, c, problem%x_lb, problem%x_ub, problem%c_lb, problem%c_ub, &
                                lambda, options%ktol, options%ctol, converged, conv_istat)
    end if
    if (options%print_level >= 1) write(*,'(A,I5,A,ES13.5,A,L1)') ' sqpopt iter ', iter, ': f = ', f, ', converged = ', converged
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
            ! `sqpopt_hessian_exact` is not yet supported; falls back to the BFGS update:
            call hessian%update_bfgs(x - x_prev, gl - gl_prev)
        end if
    end if

    if (trust_region%enabled) then

        ! trust-region globalization: re-solves the QP as needed with a
        ! shrinking radius and its own accept/reject test (merit-ratio or
        ! filter, depending on `linesearch%mode`) instead of a line search
        ! along one fixed `p` -- see [[sqpopt_trust_region_module]]:
        call trust_region%step(problem, hessian, qp_solver, linesearch, x, g, f, c, jac, x_new, new_lambda, alpha, istat)

    else

        ! solve the linearized QP subproblem for the search direction and multipliers:
        call qp_solver%solve(hessian, jac, x, g, c, problem%x_lb, problem%x_ub, &
                              problem%c_lb, problem%c_ub, p, new_lambda, qp_istat)

        ! update the merit function's penalty parameter so that it dominates the
        ! current multiplier estimates (as in `slsqp`): for `sqpopt_merit_l1` this
        ! is required for the exact penalty function's minimizer to coincide with
        ! the true constrained optimum (Han/Powell); for
        ! `sqpopt_merit_augmented_lagrangian` the same rule is used as a simple
        ! (if not exactly optimal) substitute for the theoretically-correct
        ! closed-form threshold, which would require tracking the QP's own
        ! multiplier separately from `lambda`. Without a large-enough penalty,
        ! the merit function can prefer a "compromise" infeasible point over
        ! the true solution:
        if (size(new_lambda) > 0) linesearch%penalty = max(linesearch%penalty, maxval(abs(new_lambda)) + 1.0_wp)

        ! safeguard (as in `slsqp`): if `p` is not a descent direction for the
        ! merit function (the linearized QP solve is not always guaranteed to
        ! produce one), reset the Hessian approximation to the identity and
        ! recompute `p` once from scratch:
        block
            real(wp) :: dphi0
            call linesearch%directional_derivative(jac, g, p, c, problem%c_lb, problem%c_ub, new_lambda, dphi0)
            if (dphi0 >= 0.0_wp) then
                call hessian%reset()
                call qp_solver%solve(hessian, jac, x, g, c, problem%x_lb, problem%x_ub, &
                                      problem%c_lb, problem%c_ub, p, new_lambda, istat)
                if (istat /= sqpopt_success) return
            end if
        end block

        ! second-order correction (SOC): for strongly nonlinear constraints, the
        ! full step `p` can be rejected by the merit function even when it is a
        ! genuinely good step, because the *linearized* constraint prediction
        ! differs from the true (nonlinear) constraint value at `x+p` (the
        ! "Maratos effect"). Correct for this by solving for an additional small
        ! step that accounts for the true constraint residual at `x+p`, and use
        ! the corrected step if it has a better merit function value:
        if (problem%m > 0) then
            call second_order_correction(problem, linesearch, jac, x, c, new_lambda, p)
        end if

        ! line search along `p` to (approximately) minimize the merit function
        ! (or, in `sqpopt_linesearch_filter` mode, to find a point acceptable
        ! to the filter -- that mode needs the QP's own predicted decrease in
        ! `f`, `q = -(g^Tp + 0.5*p^THp)`, computed here since only this routine
        ! has access to `hessian`; skipped for the other modes, which ignore
        ! `q`, since `hv_product` isn't free):
        block
            real(wp) :: q
            q = 0.0_wp
            if (linesearch%mode == sqpopt_linesearch_filter) then
                block
                    real(wp), dimension(problem%n) :: hp
                    call hessian%hv_product(p, hp)
                    q = -(dot_product(g, p) + 0.5_wp*dot_product(p, hp))
                end block
            end if
            call linesearch%search(problem%eval_f, problem%eval_c, x, p, f, g, c, jac, new_lambda, &
                                    problem%c_lb, problem%c_ub, q, alpha, x_new, istat)
        end block

    end if

    ! save the current point/gradient/objective for the next quasi-Newton
    ! update and stalled-progress convergence test:
    x_prev  = x
    gl_prev = gl
    f_prev  = f

    ! update the point and multipliers: `x_new` is always well-defined here
    ! (even when `istat==sqpopt_line_search_failed`, e.g. the `alpha_min`
    ! floor is deliberately still accepted -- see `armijo_line_search` --
    ! or, in `sqpopt_linesearch_watchdog` mode, `x_new` may instead be an
    ! earlier best point on backtrack). This update is never skipped based
    ! on `istat`, since `alpha`/`x_new` are always meaningful regardless of
    ! whether the sufficient-decrease test was satisfied:
    x      = x_new
    lambda = new_lambda

    if (options%print_level >= 1) write(*,'(A,I5,A,ES13.5,A,ES10.2,A,I0)') &
        ' sqpopt iter ', iter, ': f = ', f, ', alpha = ', alpha, ', qp_istat = ', qp_istat

    istat = sqpopt_success

    end subroutine sqpopt_iterate
!*******************************************************************************

    end module sqpopt_iterate_module
!*******************************************************************************
