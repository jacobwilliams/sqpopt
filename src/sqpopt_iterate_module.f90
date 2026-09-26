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
    use sqpopt_types_module,      only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_user_requested_stop, sqpopt_report_func, &
                                         sqpopt_infeasible, sqpopt_function_error, sqpopt_all_finite, sqpopt_unbounded, &
                                         sqpopt_acceptable, sqpopt_infinity
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type, sqpopt_hessian_sr1
    use sqpopt_hessian_module,    only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_filter, l1_violation
    use sqpopt_linalg_module,     only: sparse_matvec_transpose
    use sqpopt_convergence_module, only: check_convergence
    use sqpopt_soc_module,        only: soc_step
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_restoration_module,  only: restoration_step

    implicit none

    private

    public :: sqpopt_iterate, sqpopt_evaluate_point

    type, public :: sqpopt_iter_info
        !! information about one major iteration, for the iteration log
        real(wp) :: f         = 0.0_wp  !! objective at the start of the iteration (of the scaled problem)
        real(wp) :: kkt       = 0.0_wp  !! KKT error there (see [[check_convergence]])
        real(wp) :: feas      = 0.0_wp  !! feasibility error there
        real(wp) :: alpha     = 0.0_wp  !! step length taken
        real(wp) :: step_norm = 0.0_wp  !! \( \lVert x_{k+1}-x_k \rVert_2 \)
        real(wp) :: penalty   = 0.0_wp  !! merit function penalty parameter
        integer  :: qp_istat  = 0       !! status of the QP solve
        integer  :: qp_iter   = 0       !! active-set iterations of the QP solve
        logical  :: restoration = .false. !! whether a feasibility-restoration step was taken
        logical  :: stepped   = .false. !! whether the iteration got as far as computing a step
    end type sqpopt_iter_info

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  perform one major SQP iteration, updating `x` and `lambda` in place.
!  `x_prev`/`gl_prev`/`f_prev` hold the previous point, the Lagrangian
!  gradient there (evaluated with the multipliers passed back out in
!  `lambda`, as the quasi-Newton update requires), and the objective value.
!  They are used to form the quasi-Newton `(s,y)` pair and for the
!  stalled-progress and infeasibility tests. Pass them in *unallocated*
!  before the first call. The Hessian update and those two tests are
!  skipped on the first iteration, since there is no previous point yet.
!
!  On exit, `done` is true if the solver should stop at the (unchanged)
!  input point `x`, with the reason in `istat`: `sqpopt_success`,
!  `sqpopt_stalled`, or `sqpopt_infeasible` (from [[check_convergence]]),
!  `sqpopt_acceptable` (the acceptable-level test held for
!  `options%acceptable_iter` consecutive iterations, counted in
!  `n_acceptable`), `sqpopt_unbounded` (the objective is below
!  `options%obj_lower_limit` at a feasible point),
!  `sqpopt_user_requested_stop` (from the `report` callback, or a user
!  function returning `status<0`), or `sqpopt_function_error` (a problem
!  function returned a non-finite value, or failed, at `x`). Otherwise `istat` reports
!  how the step went: `sqpopt_success`, `sqpopt_qp_solve_failed` (the QP
!  solver hit its iteration limit, and its last step was used anyway), or
!  `sqpopt_line_search_failed` (no acceptable step was found, so `x` is
!  unchanged). After a failed step the Hessian approximation is reset, so
!  the next iteration tries a different direction, and `f_prev` is
!  deallocated, so the stalled-progress test (which would otherwise see
!  "no change") is skipped on the next iteration.

    subroutine sqpopt_iterate(problem, options, hessian, qp_solver, linesearch, trust_region, &
                               x, lambda, x_prev, gl_prev, f_prev, jac, n_acceptable, iter, report, done, istat, info)

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
    type(sqpopt_sparse_matrix),          intent(inout) :: jac     !! workspace for the constraint Jacobian: its sparsity
                                                                  !! structure is set on the first call (when `jac%val` is
                                                                  !! unallocated) and reused, and its values are updated
    integer,                intent(inout) :: n_acceptable !! number of consecutive iterations (so far) at which the
                                                          !! acceptable-level test has held (`0` before the 1st call)
    integer,                intent(in)    :: iter      !! major iteration number (starts at 1), passed to `report`
    procedure(sqpopt_report_func), optional, pointer :: report !! optional user progress-reporting callback (see [[sqpopt_types_module]])
    logical,                 intent(out)   :: done      !! true if the solver should stop at `x` (see `istat` for why)
    integer,                 intent(out)   :: istat     !! status code (see above and [[sqpopt_types_module]])
    type(sqpopt_iter_info),  intent(out)   :: info      !! information about this iteration, for the log

    real(wp) :: f !! current objective function value
    real(wp), dimension(problem%n) :: g, gl, p, x_new
    real(wp), dimension(problem%m) :: c, new_lambda
    real(wp) :: alpha
    integer :: qp_istat, step_istat
    logical :: restore

    done = .false.

    ! evaluate the problem functions and the sparse Jacobian at the current point:
    call problem%f(x, f)
    call problem%g(x, g)
    call problem%c(x, c)
    if (.not. allocated(jac%val)) then
        jac%nrows = problem%m
        jac%ncols = problem%n
        jac%nnz   = problem%jac_nnz
        jac%irow  = problem%jac_irow
        jac%icol  = problem%jac_icol
        allocate(jac%val(problem%jac_nnz))
    end if
    call problem%jac(x, jac%val)
    info%f = f

    ! a user function asked to stop:
    if (problem%stop_requested) then
        istat = sqpopt_user_requested_stop
        done  = .true.
        return
    end if

    ! every accepted trial point had finite `f` and `c`, so a non-finite
    ! value (or a failed evaluation) here is either at the starting point or
    ! in the derivatives:
    if (.not. (sqpopt_all_finite([f]) .and. sqpopt_all_finite(g) .and. &
               sqpopt_all_finite(c) .and. sqpopt_all_finite(jac%val))) then
        istat = sqpopt_function_error
        done  = .true.
        return
    end if

    ! report progress on the current iterate, if the user has supplied a
    ! callback, before doing any further work this iteration -- this
    ! reports every major iterate, including the initial guess (iter=1,
    ! before any step has been taken) and the final, converged point.
    ! The values are converted back to the original (unscaled) problem:
    if (present(report)) then
        if (associated(report)) then
            block
                logical :: user_stop
                user_stop = .false.
                if (associated(problem%user_data)) then
                    call report(iter, x, f/problem%f_scale, c/problem%c_scale, &
                                lambda*problem%c_scale/problem%f_scale, user_stop, problem%user_data)
                else
                    call report(iter, x, f/problem%f_scale, c/problem%c_scale, &
                                lambda*problem%c_scale/problem%f_scale, user_stop)
                end if
                if (user_stop) then
                    istat = sqpopt_user_requested_stop
                    done  = .true.
                    return
                end if
            end block
        end if
    end if

    ! check convergence at the current point before taking a step. The
    ! stalled-progress test also needs `f_prev`/`x_prev`, and the
    ! infeasibility test `x_prev`, so neither is done on the first iteration;
    ! the stalled-progress test is also skipped right after a failed step,
    ! when `f_prev` is deallocated (an unallocated actual argument counts as
    ! absent for an optional dummy argument):
    call check_convergence(x, g, jac, c, problem%x_lb, problem%x_ub, problem%c_lb, problem%c_ub, &
                            lambda, options%ktol, options%ctol, done, istat, &
                            f=f, f_prev=f_prev, x_prev=x_prev, ftol=options%ftol, xtol=options%xtol, &
                            kkt_error=info%kkt, feas_error=info%feas)
    if (done) return

    ! unbounded: the objective is below its limit at a feasible point:
    if (options%obj_lower_limit > -sqpopt_infinity .and. info%feas <= options%ctol) then
        if (f/problem%f_scale < options%obj_lower_limit) then
            istat = sqpopt_unbounded
            done  = .true.
            return
        end if
    end if

    ! acceptable-level convergence (as in IPOPT): the looser tolerances have
    ! held for `acceptable_iter` consecutive iterations:
    if (options%acceptable_iter > 0) then
        if (info%kkt <= options%acceptable_ktol .and. info%feas <= options%acceptable_ctol) then
            n_acceptable = n_acceptable + 1
        else
            n_acceptable = 0
        end if
        if (n_acceptable >= options%acceptable_iter) then
            istat = sqpopt_acceptable
            done  = .true.
            return
        end if
    end if

    ! Lagrangian gradient at the current point, with the current multipliers
    ! (the multipliers the previous iteration's `gl_prev` was also formed with):
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

    qp_istat = sqpopt_success
    restore  = .false.

    if (trust_region%enabled) then

        ! trust-region globalization: re-solves the QP as needed with a
        ! shrinking radius and its own accept/reject test (merit-ratio or
        ! filter, depending on `linesearch%mode`) instead of a line search
        ! along one fixed `p` -- see [[sqpopt_trust_region_module]]:
        call trust_region%step(problem, hessian, qp_solver, linesearch, x, g, f, c, jac, x_new, new_lambda, alpha, step_istat)

    else

        ! solve the linearized QP subproblem for the search direction and multipliers:
        call qp_solver%solve(hessian, jac, x, g, c, problem%x_lb, problem%x_ub, &
                              problem%c_lb, problem%c_ub, p, new_lambda, qp_istat)
        restore = qp_istat == sqpopt_infeasible

        if (.not. restore) then

            call update_penalty(linesearch, new_lambda)

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
                                          problem%c_lb, problem%c_ub, p, new_lambda, qp_istat)
                    restore = qp_istat == sqpopt_infeasible
                    if (.not. restore) call update_penalty(linesearch, new_lambda)
                end if
            end block

        end if

        if (restore) then

            ! the linearized constraints are inconsistent, so the QP step can't
            ! be trusted: take a step toward feasibility instead, keeping the
            ! current multipliers (see [[sqpopt_restoration_module]]):
            new_lambda = lambda
            call restoration_step(problem, jac, x, c, qp_solver%max_step, x_new, alpha, step_istat)

        else

            ! line search along `p` to (approximately) minimize the merit function
            ! (or, in `sqpopt_linesearch_filter` mode, to find a point acceptable
            ! to the filter). If the full step is rejected because of constraint
            ! curvature (the Maratos effect), the line search also tries its
            ! second-order correction (see `soc` below):
            if (problem%m > 0) then
                call linesearch%search(eval_f_cached, eval_c_cached, x, p, f, g, c, jac, new_lambda, &
                                        problem%c_lb, problem%c_ub, alpha, x_new, step_istat, soc=soc)
            else
                call linesearch%search(eval_f_cached, eval_c_cached, x, p, f, g, c, jac, new_lambda, &
                                        problem%c_lb, problem%c_ub, alpha, x_new, step_istat)
            end if

            if (step_istat /= sqpopt_success .and. linesearch%mode == sqpopt_linesearch_filter) then
                ! the filter line search failed (no acceptable step length): as
                ! in Wächter & Biegler's method, add the current point to the
                ! filter and, if infeasible, take a feasibility restoration
                ! step instead (keeping the current multipliers):
                block
                    real(wp) :: theta
                    theta = l1_violation(c, problem%c_lb, problem%c_ub)
                    call linesearch%filter_record(theta, f)
                    if (theta > 0.0_wp) then
                        restore    = .true.
                        new_lambda = lambda
                        call restoration_step(problem, jac, x, c, qp_solver%max_step, x_new, alpha, step_istat)
                    end if
                end block
            end if

        end if

    end if

    ! a user function asked to stop during the step (the point is left unchanged):
    if (problem%stop_requested) then
        istat = sqpopt_user_requested_stop
        done  = .true.
        return
    end if

    ! save the current point/objective for the next stalled-progress test,
    ! and the Lagrangian gradient at the current point *evaluated with the
    ! new multipliers*, so that the next quasi-Newton pair is
    ! y = grad L(x_new, lambda_new) - grad L(x, lambda_new):
    x_prev = x
    f_prev = f
    block
        real(wp), dimension(problem%n) :: jtlam
        call sparse_matvec_transpose(jac, new_lambda, jtlam)
        gl_prev = g - jtlam
    end block

    ! update the point and multipliers: `x_new` is always well-defined here
    ! (even when `step_istat==sqpopt_line_search_failed`, e.g. the `alpha_min`
    ! floor is deliberately still accepted -- see `armijo_line_search` --
    ! or, in `sqpopt_linesearch_watchdog` mode, `x_new` may instead be an
    ! earlier best point on backtrack). This update is never skipped based
    ! on the status, since `alpha`/`x_new` are always meaningful regardless of
    ! whether the sufficient-decrease test was satisfied:
    info%stepped     = .true.
    info%alpha       = alpha
    info%step_norm   = norm2(x_new - x)
    info%penalty     = linesearch%penalty
    info%qp_istat    = qp_istat
    info%qp_iter     = qp_solver%n_iter
    info%restoration = restore

    x      = x_new
    lambda = new_lambda

    if (step_istat /= sqpopt_success) then
        ! no (acceptable) step was found along this direction: start the
        ! next iteration from a fresh Hessian approximation, so it computes
        ! a different direction, and don't let the stalled-progress test
        ! mistake "no step taken" for convergence:
        call hessian%reset()
        if (allocated(f_prev)) deallocate(f_prev)
    end if

    ! report the first failure, if any (a QP that stopped at its iteration
    ! limit still produced a usable step, so it's reported below a line-search
    ! failure):
    if (step_istat /= sqpopt_success) then
        istat = step_istat
    else if (qp_istat /= sqpopt_success .and. .not. restore) then
        istat = qp_istat
    else
        istat = sqpopt_success
    end if

    contains

        subroutine eval_f_cached(xx, ff)
        !! `f`, through the problem's evaluation cache (see [[sqpopt_problem_module]])
        real(wp), dimension(:), intent(in)  :: xx
        real(wp),               intent(out) :: ff
        call problem%f(xx, ff)
        end subroutine eval_f_cached

        subroutine eval_c_cached(xx, cc)
        !! `c`, through the problem's evaluation cache (see [[sqpopt_problem_module]])
        real(wp), dimension(:), intent(in)  :: xx
        real(wp), dimension(:), intent(out) :: cc
        call problem%c(xx, cc)
        end subroutine eval_c_cached

        subroutine soc(p_trial, c_trial, p_soc, ok)
        !! the second-order correction of a rejected trial step, given to
        !! the line search (see [[sqpopt_soc_module]])
        real(wp), dimension(:), intent(in)  :: p_trial !! the rejected trial step
        real(wp), dimension(:), intent(in)  :: c_trial !! constraint values at `x+p_trial`
        real(wp), dimension(:), intent(out) :: p_soc   !! the corrected step
        logical,                intent(out) :: ok      !! true if `p_soc` is usable
        call soc_step(jac, x, p_trial, c, c_trial, problem%c_lb, problem%c_ub, problem%x_lb, problem%x_ub, p_soc, ok)
        end subroutine soc

    end subroutine sqpopt_iterate
!*******************************************************************************

!*******************************************************************************
!>
!  evaluate the (scaled) problem at `x`, and the resulting KKT and
!  feasibility errors and variable-bound multipliers `z` (with the
!  constraint multipliers `lambda`): \( z = g - J^T\lambda \) for a variable
!  at one of its bounds, zero otherwise (for the Lagrangian
!  \( f - \lambda^Tc - z^Tx \)). Used to report the final state of a solve.

    subroutine sqpopt_evaluate_point(problem, options, x, lambda, jac, f, c, kkt, feas, z)

    type(sqpopt_problem_type),  intent(inout) :: problem !! problem definition
    type(sqpopt_options_type),  intent(in)    :: options !! solver options
    real(wp), dimension(:),     intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: lambda  !! constraint multipliers `dimension(m)`
    type(sqpopt_sparse_matrix), intent(inout) :: jac     !! Jacobian workspace (as for [[sqpopt_iterate]])
    real(wp),                   intent(out)   :: f       !! objective at `x`
    real(wp), dimension(:),     intent(out)   :: c       !! constraints at `x` `dimension(m)`
    real(wp),                   intent(out)   :: kkt     !! KKT error at `x`
    real(wp),                   intent(out)   :: feas    !! feasibility error at `x`
    real(wp), dimension(:),     intent(out)   :: z       !! variable-bound multipliers `dimension(n)`

    real(wp), dimension(size(x)) :: g, jtlam
    logical :: converged
    integer :: istat, j

    call problem%f(x, f)
    call problem%g(x, g)
    call problem%c(x, c)
    if (.not. allocated(jac%val)) then
        jac%nrows = problem%m
        jac%ncols = problem%n
        jac%nnz   = problem%jac_nnz
        jac%irow  = problem%jac_irow
        jac%icol  = problem%jac_icol
        allocate(jac%val(problem%jac_nnz))
    end if
    call problem%jac(x, jac%val)

    call check_convergence(x, g, jac, c, problem%x_lb, problem%x_ub, problem%c_lb, problem%c_ub, &
                            lambda, options%ktol, options%ctol, converged, istat, &
                            kkt_error=kkt, feas_error=feas)

    call sparse_matvec_transpose(jac, lambda, jtlam)
    z = 0.0_wp
    do j = 1, size(x)
        if (x(j) - problem%x_lb(j) <= options%ctol .or. problem%x_ub(j) - x(j) <= options%ctol) &
            z(j) = g(j) - jtlam(j)
    end do

    end subroutine sqpopt_evaluate_point
!*******************************************************************************

!*******************************************************************************
!>
!  update the merit function's penalty parameter so that it dominates the
!  current multiplier estimates (as in `slsqp`): for `sqpopt_merit_l1` this
!  is required for the exact penalty function's minimizer to coincide with
!  the true constrained optimum (Han/Powell); for
!  `sqpopt_merit_augmented_lagrangian` the same rule is used as a simple
!  (if not exactly optimal) substitute for the theoretically-correct
!  closed-form threshold, which would require tracking the QP's own
!  multiplier separately from `lambda`. Without a large-enough penalty,
!  the merit function can prefer a "compromise" infeasible point over
!  the true solution.

    subroutine update_penalty(linesearch, lambda)

    type(sqpopt_linesearch_type), intent(inout) :: linesearch !! merit function / line search
    real(wp), dimension(:),       intent(in)    :: lambda     !! current multiplier estimates `dimension(m)`

    if (size(lambda) > 0) linesearch%penalty = max(linesearch%penalty, maxval(abs(lambda)) + 1.0_wp)

    end subroutine update_penalty
!*******************************************************************************

    end module sqpopt_iterate_module
!*******************************************************************************
