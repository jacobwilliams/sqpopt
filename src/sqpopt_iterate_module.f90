!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The core SQP major iteration ([[sqpopt_iterate]]): evaluates the problem
!  functions, tests for convergence, updates the Hessian approximation
!  (limited-memory quasi-Newton, or the user's exact Hessian with an
!  inertia-correcting shift, found from the QP solver's tests or, with
!  `options%inertia_control`, from a factorization, which also corrects
!  SR1: see [[sqpopt_inertia_module]]; before the exact Hessian is
!  evaluated, multipliers that came from a QP with a large shift are
!  re-estimated by least squares), solves the QP subproblem for the
!  search direction (directly, if it can, with `options%direct_qp`: see
!  [[sqpopt_qp_direct_module]]; and re-solving it with diverging-multiplier
!  constraints elastic, see `options%elastic_multiplier_limit`), and takes
!  the step: by a line
!  search (with second-order corrections), or a trust-region step, or, when
!  no acceptable step is found at an infeasible point or the QP is
!  inconsistent, a feasibility restoration step or phase (see
!  [[sqpopt_restoration_module]]). What happened is returned in a
!  [[sqpopt_iter_info]], for the iteration log (see `options%print_level`);
!  at `print_level >= 3`, the details are also written to the detailed log
!  (see [[sqpopt_log_module]]). No dense `n x n` or `m x n` matrix is ever
!  formed (except by the dense QP solver).

    module sqpopt_iterate_module

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_types_module,      only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_user_requested_stop, sqpopt_report_func, &
                                         sqpopt_infeasible, sqpopt_function_error, sqpopt_all_finite, sqpopt_unbounded, &
                                         sqpopt_acceptable, sqpopt_infinity, sqpopt_stalled, sqpopt_qp_solve_failed, &
                                         sqpopt_out_of_memory
    use sqpopt_problem_module,    only: sqpopt_problem_type, sqpopt_derivatives_fast, sqpopt_derivatives_accurate
    use sqpopt_options_module,    only: sqpopt_options_type
    use sqpopt_hessian_module,    only: sqpopt_hessian_type, sqpopt_hessian_sr1, sqpopt_hessian_exact
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_filter, sqpopt_linesearch_funnel, &
                                         l1_violation
    use sqpopt_linalg_module,     only: sparse_matvec_transpose
    use sqpopt_convergence_module, only: check_convergence
    use sqpopt_soc_module,        only: soc_step
    use sqpopt_log_module,        only: sqpopt_log_type, sqpopt_log_detail, fmt_e, fmt_i, plural, qp_status_text
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_restoration_module,  only: restoration_step, escape_step, sqpopt_restoration_type, &
                                          sqpopt_restoration_phase
    use sqpopt_inertia_module,      only: sqpopt_inertia_type
    use sqpopt_kkt_module,          only: sqpopt_kkt_type
    use sqpopt_least_squares_module, only: sqpopt_least_squares_type, multiplier_estimate
    use sqpopt_qp_direct_module,    only: direct_outcome_text

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
        logical  :: restoration = .false. !! whether a feasibility-restoration step was taken (a single step, or
                                          !! an iteration of a restoration phase: see `phase`)
        logical  :: phase     = .false. !! whether the step was an iteration of a restoration phase
        logical  :: stepped   = .false. !! whether the iteration got as far as computing a step
        logical  :: soc       = .false. !! whether the accepted step was second-order corrected
        logical  :: hess_reset = .false. !! whether the Hessian approximation was reset (or its shift increased)
        logical  :: elastic   = .false. !! whether the QP was re-solved with diverging-multiplier constraints elastic
        logical  :: escape    = .false. !! whether an escape step (from a stationary point of the violation) was taken
        logical  :: nonmonotone = .false. !! whether the step came from the line search's non-monotone retry
        logical  :: relaxed   = .false. !! whether the step was a watchdog relaxed step
        real(wp) :: stat_unscaled = 0.0_wp !! stationarity error of the *unscaled* problem at the start of the iteration
        real(wp) :: lam_max   = 0.0_wp  !! largest multiplier magnitude of the unscaled problem, after the step
        real(wp) :: glob      = 0.0_wp  !! the globalization's state after the step: the merit penalty, the number
                                        !! of filter entries, the funnel width, or the trust-region radius
        real(wp) :: hess_measure = 0.0_wp !! the number of stored quasi-Newton pairs, or the exact Hessian's shift
        integer  :: n_fc      = 0       !! calls of `fc` during the iteration (set by the caller)
        logical  :: derivatives = .false. !! whether the solver switched from fast to accurate derivatives (see
                                          !! `options%derivative_accuracy`)
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
!  function returning `status<0`), `sqpopt_function_error` (a problem
!  function returned a non-finite value, or failed, at `x`), or
!  `sqpopt_out_of_memory` (a QP solve couldn't allocate its matrices, or a
!  factorization ran out of memory).
!  Otherwise `istat` reports
!  how the step went: `sqpopt_success`, `sqpopt_qp_solve_failed` (the QP
!  solver hit its iteration limit, and its last step was used anyway), or
!  `sqpopt_line_search_failed` (no acceptable step was found, so `x` is
!  unchanged). After a failed step the Hessian approximation is reset, so
!  the next iteration tries a different direction, and `f_prev` is
!  deallocated, so the stalled-progress test (which would otherwise see
!  "no change") is skipped on the next iteration. `info` describes the
!  iteration (its measures and events), for the iteration log.

    subroutine sqpopt_iterate(problem, options, hessian, qp_solver, linesearch, trust_region, &
                               x, lambda, x_prev, gl_prev, f_prev, viol_prev, jac, n_acceptable, n_stalled, n_escape, &
                               restoration, inertia, kkt, least_squares, iter, report, done, &
                               istat, info)

    type(sqpopt_problem_type),    intent(inout)   :: problem      !! problem definition
    type(sqpopt_options_type),    intent(in)      :: options      !! solver options
    type(sqpopt_hessian_type),    intent(inout)   :: hessian      !! Hessian of the Lagrangian approximation
    type(sqpopt_qp_solver_type),  intent(inout)   :: qp_solver    !! QP subproblem solver
    type(sqpopt_linesearch_type), intent(inout)   :: linesearch   !! merit function / line search
    type(sqpopt_trust_region_type), intent(inout) :: trust_region !! trust-region globalization (used instead of
                                                                   !! `linesearch` when `trust_region%enabled`)
    real(wp), dimension(:), intent(inout) :: x       !! current point, updated on exit `dimension(n)`
    real(wp), dimension(:), intent(inout) :: lambda  !! current Lagrange multipliers, updated on exit `dimension(m)`
    real(wp), dimension(:), allocatable, intent(inout) :: x_prev  !! previous point (unallocated before the 1st call)
    real(wp), dimension(:), allocatable, intent(inout) :: gl_prev !! previous Lagrangian gradient (unallocated before the 1st call)
    real(wp),               allocatable, intent(inout) :: f_prev  !! previous objective value (unallocated before the 1st call)
    real(wp),               allocatable, intent(inout) :: viol_prev !! previous constraint violation (unallocated before the 1st call)
    type(sqpopt_sparse_matrix),          intent(inout) :: jac     !! workspace for the constraint Jacobian: its sparsity
                                                                  !! structure is set on the first call (when `jac%val` is
                                                                  !! unallocated) and reused, and its values are updated
    integer,                intent(inout) :: n_stalled    !! number of consecutive iterations (so far) at which the
                                                          !! stalled-progress test has held (`0` before the 1st call)
    integer,                intent(inout) :: n_acceptable !! number of consecutive iterations (so far) at which the
                                                          !! acceptable-level test has held (`0` before the 1st call)
    integer,                intent(inout) :: n_escape     !! number of second-order escapes from a stationary point of
                                                          !! the violation taken so far (`0` before the 1st call)
    type(sqpopt_restoration_type), intent(inout) :: restoration !! feasibility restoration phase state (see
                                                                !! [[sqpopt_restoration_module]])
    type(sqpopt_inertia_type), intent(inout) :: inertia !! inertia control of the exact or SR1 Hessian (used if
                                                        !! `inertia%enabled`, see [[sqpopt_inertia_module]])
    type(sqpopt_kkt_type),     intent(inout) :: kkt     !! the KKT matrix of the QP's working set, for the inertia
                                                        !! control and the QP solver's direct method (used if
                                                        !! `kkt%enabled`, see [[sqpopt_kkt_module]])
    type(sqpopt_least_squares_type), intent(inout) :: least_squares !! direct least-squares solver, for the
                                                        !! restoration steps and second-order corrections (used
                                                        !! if `least_squares%enabled`, see
                                                        !! [[sqpopt_least_squares_module]])
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

    integer, parameter :: max_escape = 3 !! maximum number of second-order escapes (see [[escape_step]])
    type(sqpopt_log_type) :: lg !! the detailed log
    real(wp) :: stat_err
    integer :: n_stalled0 !! `n_stalled` on entry
    logical :: slow       !! whether the objective changed little from the previous iterate (see
                          !! `options%acceptable_obj_change_tol`)
    real(wp) :: shift_floor !! with inertia control: the Hessian's shift at the start of the iteration (the part
                            !! of it that is carried over from failed steps)
    logical :: keep_shift   !! with inertia control: whether the step failed, so the shift it was computed with
                            !! is carried over (increased) to the next iteration

    done = .false.
    lg = linesearch%log   ! (the detailed log, set up by `solve`)

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

    ! (the factorizations kept from the previous iteration no longer apply:
    ! `estimate_multipliers`, below, may use the least-squares one)
    call kkt%new_matrices()
    call least_squares%new_matrices()

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

    ! with the exact Hessian: if the last step's shift dominated it, the QP's
    ! multipliers are mostly an artifact of the shift, and the Hessian that
    ! is about to be evaluated with them would be too (see `estimate_multipliers`)
    if (hessian%shift_dominant) call estimate_multipliers()
    hessian%shift_dominant = .false.

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
    ! With fast derivatives (see `options%derivative_accuracy`), switch to
    ! accurate ones, for the rest of the solve, as soon as the point is near
    ! a solution (the KKT and feasibility errors are below
    ! `derivative_switch_tol`), or the solver would stop here (convergence,
    ! acceptable level, or infeasibility, which inaccurate derivatives may
    ! fake), or progress has stalled; then redo the tests with them:
    n_stalled0 = n_stalled
    do
        call convergence_tests()
        if (problem%derivative_accuracy /= sqpopt_derivatives_fast) exit
        if (.not. (done .or. n_stalled > 0 .or. &
                   (info%kkt <= options%derivative_switch_tol .and. info%feas <= options%derivative_switch_tol) .or. &
                   (options%acceptable_iter > 0 .and. n_acceptable + 1 >= options%acceptable_iter .and. &
                    info%kkt <= options%acceptable_ktol .and. info%feas <= options%acceptable_ctol))) exit
        call problem%set_derivative_accuracy(sqpopt_derivatives_accurate)
        info%derivatives = .true.
        call lg%put(sqpopt_log_detail, 'switched to accurate derivatives (KKT error '//fmt_e(info%kkt)// &
                    ', feasibility error '//fmt_e(info%feas)//'): gradient and Jacobian re-evaluated')
        call problem%g(x, g)
        call problem%jac(x, jac%val)
        call kkt%new_matrices()
        call least_squares%new_matrices()
        if (problem%stop_requested) then
            istat = sqpopt_user_requested_stop
            done  = .true.
            return
        end if
        if (.not. (sqpopt_all_finite(g) .and. sqpopt_all_finite(jac%val))) then
            istat = sqpopt_function_error
            done  = .true.
            return
        end if
        done      = .false.
        n_stalled = n_stalled0
    end do

    ! a point that is stationary for the violation may still be a saddle of
    ! it (e.g. on a symmetry plane of the problem, which exactly computed
    ! steps never leave): before declaring the problem infeasible, look for a
    ! second-order decrease of the violation (see [[escape_step]]), and if
    ! one is found, continue from there (a limited number of times):
    if (done .and. istat == sqpopt_infeasible .and. n_escape < max_escape) then
        call escape_step(problem, jac, x, c, options%ktol, x_new, step_istat)
        if (step_istat == sqpopt_success) then
            call lg%put(sqpopt_log_detail, 'escape step from a stationary point of the violation '// &
                        '(the Hessian approximation is reset)')
            info%escape     = .true.
            info%hess_reset = .true.
            n_escape = n_escape + 1
            done  = .false.
            istat = sqpopt_success
            x_prev = x
            if (allocated(f_prev)) deallocate(f_prev)
            viol_prev = maxval(max(problem%c_lb-c, 0.0_wp) + max(c-problem%c_ub, 0.0_wp))
            block
                real(wp), dimension(problem%n) :: jtlam
                call sparse_matvec_transpose(jac, lambda, jtlam)
                gl_prev = g - jtlam
            end block
            call hessian%reset()
            info%stepped     = .true.
            info%alpha       = 1.0_wp
            info%step_norm   = norm2(x_new - x)
            info%penalty     = linesearch%merit%penalty
            info%restoration = .true.
            x = x_new
            return
        end if
    end if
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
    ! held for `acceptable_iter` consecutive iterations, during which the
    ! objective (of the original problem) no longer changed much:
    if (options%acceptable_iter > 0) then
        slow = .true.
        if (allocated(f_prev)) slow = abs(f - f_prev) <= &
                                      options%acceptable_obj_change_tol*max(problem%f_scale, abs(f))
        if (info%kkt <= options%acceptable_ktol .and. info%feas <= options%acceptable_ctol .and. slow) then
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

    if (options%hessian_mode == sqpopt_hessian_exact) then
        ! the user's exact Hessian of the Lagrangian, at the current point
        ! with the current multipliers (see [[sqpopt_hessian_module]]):
        block
            real(wp), dimension(problem%hess_nnz) :: hval
            call problem%hess(x, lambda, hval)
            if (problem%stop_requested) then
                istat = sqpopt_user_requested_stop
                done  = .true.
                return
            end if
            if (.not. sqpopt_all_finite(hval)) then
                istat = sqpopt_function_error
                done  = .true.
                return
            end if
            ! (the shift is decreased only after a good step: not after a
            ! failed one, or during a run of very short ones)
            call hessian%set_values(hval, decay=qp_solver%n_short == 0 .and. allocated(f_prev))
        end block
    else if (allocated(x_prev) .and. .not. info%derivatives) then
        ! update the quasi-Newton Hessian approximation using the previous step
        ! (skipped on the very first iteration, since there is no previous point,
        ! and when the derivatives were just switched to accurate ones, since
        ! `gl_prev` was formed with the fast ones):
        if (options%hessian_mode == sqpopt_hessian_sr1) then
            call hessian%update_sr1(x - x_prev, gl - gl_prev)
        else
            call hessian%update_bfgs(x - x_prev, gl - gl_prev)
        end if
    end if

    ! (the Hessian has changed, so a factorization with it no longer applies)
    call kkt%new_matrices()

    qp_istat = sqpopt_success
    restore  = .false.
    shift_floor = hessian%shift
    keep_shift  = .false.

    if (restoration%active) then

        ! in a feasibility restoration phase: minimize the violation, keeping
        ! the current multipliers (see [[sqpopt_restoration_module]]):
        restore    = .true.
        new_lambda = lambda
        linesearch%merit%joint_active = .false.
        info%phase = .true.
        call restoration_phase_iteration()

    else if (trust_region%enabled) then

        ! with inertia control, first shift the exact Hessian as the working
        ! set that the QP starts from needs (the trust region's QPs are not
        ! re-solved for the working set they end with: its bounds on the step
        ! keep them bounded):
        if (inertia%enabled) then
            block
                logical :: shifted, inertia_ok
                call correct_inertia(.true., shifted, inertia_ok)
            end block
        end if

        ! trust-region globalization: re-solves the QP as needed with a
        ! shrinking radius and its own accept/reject test (merit-ratio,
        ! filter, or funnel, depending on `linesearch%mode`) instead of a line
        ! search along one fixed `p` -- see [[sqpopt_trust_region_module]]:
        call trust_region%step(problem, hessian, qp_solver, linesearch, x, g, f, c, jac, x_new, new_lambda, alpha, &
                               step_istat, kkt=kkt, inertia=inertia, least_squares=least_squares)
        if (hessian%shift > shift_floor) info%hess_reset = .true.
        info%soc = trust_region%used_soc

        ! if no step was acceptable at an infeasible point, start a
        ! restoration phase (as in Fletcher & Leyffer's trust-region filter SQP):
        if (step_istat /= sqpopt_success .and. options%restoration_mode == sqpopt_restoration_phase) then
            ! (only if infeasible beyond `ctol`, as the convergence test measures it:
            ! a round-off violation has nothing for restoration to reduce)
            if (info%feas > options%ctol) then
                restore    = .true.
                new_lambda = lambda
                call start_restoration_phase(record=.true.)
            end if
        end if

    else

        ! with inertia control, first shift the Hessian as the working set
        ! that the QP will most likely end with needs (the one it starts
        ! from), so that the QP usually has to be solved only once:
        if (inertia%enabled) then
            block
                logical :: shifted, inertia_ok
                call correct_inertia(.true., shifted, inertia_ok)
            end block
        end if

        ! solve the linearized QP subproblem for the search direction and multipliers:
        call solve_qp()
        if (qp_istat == sqpopt_out_of_memory) then
            ! (there is no step to search along: stop here, see below)
            istat = sqpopt_out_of_memory
            done  = .true.
            return
        end if
        restore = qp_istat == sqpopt_infeasible
        if (.not. restore) call elastic_resolve()

        if (.not. restore) then

            call update_penalty()

            ! safeguard (as in `slsqp`): if `p` is not a descent direction for the
            ! merit function (the linearized QP solve is not always guaranteed to
            ! produce one), reset the Hessian approximation to the identity and
            ! recompute `p` once from scratch. With the exact Hessian, which may
            ! be indefinite, "reset" increases its shift instead (see
            ! [[sqpopt_hessian_module]]), and the QP is re-solved as often as
            ! needed, also when it failed or found negative curvature (a
            ! nonconvex QP: this is the inertia correction). With inertia
            ! control (for the exact or the SR1 Hessian), the test for a
            ! nonconvex QP is instead the inertia of the KKT matrix of the QP's
            ! final working set, which also gives the shift (see
            ! [[sqpopt_inertia_module]]):
            block
                real(wp) :: dphi0
                integer :: n_shift
                logical :: shifted, inertia_ok, unusable
                integer, parameter :: max_shift = 15 !! (from the smallest shift to the largest, x10 each time)
                n_shift = 0
                do
                    call linesearch%merit%directional_derivative(jac, g, p, c, problem%c_lb, problem%c_ub, new_lambda, dphi0)
                    shifted = .false.
                    if (options%hessian_mode == sqpopt_hessian_exact .or. inertia%enabled) then
                        inertia_ok = .false.
                        if (inertia%enabled .and. n_shift < max_shift) call correct_inertia(.false., shifted, inertia_ok)
                        unusable = shifted .or. dphi0 >= 0.0_wp .or. qp_istat == sqpopt_qp_solve_failed
                        ! (negative curvature that the QP met before its final working
                        ! set is not a reason to shift, if the inertia there is right)
                        if (.not. inertia_ok) unusable = unusable .or. qp_solver%negative_curvature
                        if (.not. unusable .or. n_shift >= max_shift) exit
                        n_shift = n_shift + 1
                    else
                        if (.not. (dphi0 >= 0.0_wp) .or. n_shift >= 1) exit
                        n_shift = 1
                    end if
                    info%hess_reset = .true.
                    if (shifted) then
                        call lg%put(sqpopt_log_detail, 'QP re-solved with the shifted Hessian')
                    else if (inertia%enabled) then
                        call inertia%raise(hessian)
                        call lg%put(sqpopt_log_detail, 'QP step not usable (failed, or not a descent '// &
                                    'direction): Hessian shift increased to '//fmt_e(hessian%shift)//', QP re-solved')
                    else if (options%hessian_mode == sqpopt_hessian_exact) then
                        call hessian%reset()
                        call lg%put(sqpopt_log_detail, 'QP step not usable (nonconvex, failed, or not a descent '// &
                                    'direction): Hessian shift increased to '//fmt_e(hessian%shift)//', QP re-solved')
                    else
                        call hessian%reset()
                        call kkt%new_matrices()   ! (the factorization kept was for the pairs just discarded)
                        call lg%put(sqpopt_log_detail, 'QP step is not a descent direction: Hessian approximation '// &
                                    'reset, QP re-solved')
                    end if
                    call solve_qp()
                    restore = qp_istat == sqpopt_infeasible
                    if (restore) exit
                    call update_penalty()
                end do
            end block

        end if

        if (restore) then

            ! the linearized constraints are inconsistent, so the QP step can't
            ! be trusted: start a restoration phase, or take a step toward
            ! feasibility, keeping the current multipliers (see
            ! [[sqpopt_restoration_module]]):
            new_lambda = lambda
            linesearch%merit%joint_active = .false.
            call lg%put(sqpopt_log_detail, 'linearized constraints are inconsistent: restoration step')
            call restoration_step(problem, jac, x, c, qp_solver%max_step*qp_solver%step_scale, x_new, alpha, step_istat, &
                                  least_squares=least_squares)
            call note_restoration_step('Gauss-Newton')
            if (step_istat /= sqpopt_success .and. norm2(p) > 0.0_wp) then
                ! no first-order decrease of the violation is possible from `x`
                ! (it is stationary for the violation, e.g. `J=0` at a maximum
                ! of it): try the QP's elastic step instead, along which the
                ! violation may still decrease to second order. If this fails
                ! too, the next iteration's infeasibility test stops at `x`.
                call restoration_step(problem, jac, x, c, qp_solver%max_step*qp_solver%step_scale, &
                                      x_new, alpha, step_istat, direction=p)
                call note_restoration_step('along the elastic QP step')
            end if

        else

            ! line search along `p` to (approximately) minimize the merit function
            ! (or, in `sqpopt_linesearch_filter`/`_funnel` mode, to find a point acceptable
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

            info%soc         = linesearch%used_soc
            info%nonmonotone = linesearch%used_nonmonotone
            info%relaxed     = linesearch%used_relaxed

            ! with the augmented Lagrangian's joint step, the new multipliers are
            ! those along the step, not the QP's (see [[update_penalty_parameter]]):
            if (linesearch%merit%joint_active) new_lambda = lambda + alpha*(new_lambda - lambda)

            ! adapt the step-length cap like a trust radius (see [[sqpopt_qp_solver_module]]):
            if (step_istat == sqpopt_success .and. qp_solver%capped .and. alpha >= 1.0_wp) then
                qp_solver%step_scale = min(2.0_wp*qp_solver%step_scale, 1.0e10_wp)
            else if (alpha < 1.0_wp) then
                qp_solver%step_scale = max(1.0_wp, 0.5_wp*qp_solver%step_scale)
            end if
            ! a run of very short steps means the quasi-Newton model is poor
            ! (e.g. the search is crawling along a curved constraint): start
            ! it afresh, as after a failed step
            if (step_istat == sqpopt_success .and. alpha < 1.0e-2_wp) then
                qp_solver%n_short = qp_solver%n_short + 1
            else
                qp_solver%n_short = 0
            end if
            if (qp_solver%n_short >= 3) then
                call hessian%reset()
                keep_shift = .true.
                info%hess_reset = .true.
                call lg%put(sqpopt_log_detail, '3 very short steps in a row: Hessian approximation reset')
                qp_solver%n_short = 0
            end if

            if (step_istat /= sqpopt_success .and. (linesearch%mode == sqpopt_linesearch_filter .or. &
                                                     linesearch%mode == sqpopt_linesearch_funnel)) then
                ! the filter (or funnel) line search failed (no acceptable step
                ! length): as in Wächter & Biegler's method, add the current
                ! point to the filter (or shrink the funnel toward it) and, if
                ! infeasible, take a feasibility restoration step instead
                ! (keeping the current multipliers):
                block
                    real(wp) :: theta
                    theta = l1_violation(c, problem%c_lb, problem%c_ub)
                    if (linesearch%mode == sqpopt_linesearch_funnel) then
                        call linesearch%funnel%restoration(theta)
                        call lg%put(sqpopt_log_detail, 'no acceptable step: funnel tightened to width '// &
                                    fmt_e(linesearch%funnel%width))
                    else
                        call linesearch%filter%record(theta, f)
                        call lg%put(sqpopt_log_detail, 'no acceptable step: the current point is added to the filter')
                    end if
                    ! (restoration only if infeasible beyond `ctol`, as the convergence
                    ! test measures it: a round-off violation has nothing to reduce)
                    if (info%feas > options%ctol) then
                        restore    = .true.
                        new_lambda = lambda
                        if (options%restoration_mode == sqpopt_restoration_phase) then
                            call start_restoration_phase(record=.false.)
                        else
                            call restoration_step(problem, jac, x, c, qp_solver%max_step*qp_solver%step_scale, &
                                                  x_new, alpha, step_istat, least_squares=least_squares)
                            call note_restoration_step('Gauss-Newton')
                        end if
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

    ! a QP solve of the step (any of them: a re-solve, or one of the trust
    ! region's or the restoration phase's) couldn't allocate its matrices, or
    ! a factorization ran out of memory (the point is left unchanged):
    if (qp_solver%out_of_memory .or. kkt%solver%out_of_memory .or. least_squares%kkt%solver%out_of_memory) then
        istat = sqpopt_out_of_memory
        done  = .true.
        return
    end if

    ! save the current point/objective for the next stalled-progress test,
    ! and the Lagrangian gradient at the current point *evaluated with the
    ! new multipliers*, so that the next quasi-Newton pair is
    ! y = grad L(x_new, lambda_new) - grad L(x, lambda_new):
    x_prev = x
    f_prev = f
    viol_prev = 0.0_wp
    if (problem%m > 0) viol_prev = maxval(max(problem%c_lb-c, 0.0_wp) + max(c-problem%c_ub, 0.0_wp))
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
    info%penalty     = linesearch%merit%penalty
    info%qp_istat    = qp_istat
    info%qp_iter     = qp_solver%n_iter
    info%restoration = restore
    if (problem%m > 0) info%lam_max = maxval(abs(new_lambda*problem%c_scale))/problem%f_scale
    info%glob         = glob_value()
    info%hess_measure = merge(hessian%shift, real(hessian%n_history, wp), options%hessian_mode == sqpopt_hessian_exact)

    ! (see `estimate_multipliers`: the shift's term of the QP's stationarity
    ! condition, against the gradient's)
    hessian%shift_dominant = options%hessian_mode == sqpopt_hessian_exact .and. .not. restore .and. problem%m > 0 .and. &
                             hessian%shift*maxval(abs(x_new - x)) >= max(1.0_wp, maxval(abs(g)))

    x      = x_new
    lambda = new_lambda

    if (step_istat /= sqpopt_success) then
        ! no (acceptable) step was found along this direction: start the
        ! next iteration from a fresh Hessian approximation, so it computes
        ! a different direction, and don't let the stalled-progress test
        ! mistake "no step taken" for convergence:
        call hessian%reset()
        keep_shift = .true.
        info%hess_reset = .true.
        call lg%put(sqpopt_log_detail, 'no acceptable step: Hessian approximation reset')
        if (allocated(f_prev)) deallocate(f_prev)
        if (problem%derivative_accuracy == sqpopt_derivatives_fast) then
            ! (the failure may be due to the inaccurate derivatives: use accurate
            ! ones from the next iteration on; `x_prev` is dropped, since
            ! `gl_prev` was formed with the fast ones)
            call problem%set_derivative_accuracy(sqpopt_derivatives_accurate)
            info%derivatives = .true.
            call lg%put(sqpopt_log_detail, 'switched to accurate derivatives (from the next iteration)')
            if (allocated(x_prev)) deallocate(x_prev)
        end if
    end if

    ! with inertia control, the shift that this step needed is found again at
    ! the next iteration (starting from a third of it): only the increases
    ! after a failed step, or a run of very short ones, are carried over
    if (inertia%enabled) then
        if (hessian%shift > 0.0_wp) inertia%shift_last = hessian%shift
        if (.not. keep_shift) hessian%shift = shift_floor
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

        subroutine estimate_multipliers()
        !! replace the constraint multipliers by their least-squares estimate at
        !! the current point (see [[multiplier_estimate]]), over the equality
        !! constraints and the inequalities whose multiplier is nonzero, in the
        !! variables that are not at a bound. An inequality's multiplier is only
        !! replaced if the estimate has the same sign.
        !!
        !! This is for the exact Hessian, after a step whose shift
        !! \( \delta \) dominated it
        !! (\( \delta \lVert p \rVert_\infty \ge \max(1, \lVert g \rVert_\infty) \)).
        !! The QP's multipliers satisfy \( g + (H + \delta I)p = J^T\lambda \),
        !! so they then mostly balance the shift's term \( \delta p \), and can
        !! be orders of magnitude too large. The exact Hessian is evaluated with
        !! the multipliers, so it would be too, and need a still larger shift:
        !! on a hanging-chain problem this fed on itself until the multipliers
        !! were `1e13`, the steps `1e-9`, and the solver stopped as stalled far
        !! from the solution. The estimate doesn't depend on the shift.
        logical,  dimension(problem%m) :: rows
        real(wp), dimension(problem%m) :: estimate
        logical :: ok
        integer :: i
        rows = (problem%c_ub - problem%c_lb <= 0.0_wp) .or. lambda /= 0.0_wp
        estimate = lambda
        call multiplier_estimate(jac, rows, x > problem%x_lb .and. x < problem%x_ub, g, estimate, ok, &
                                 least_squares=least_squares)
        if (.not. ok) return
        do i = 1, problem%m
            if (problem%c_ub(i) - problem%c_lb(i) <= 0.0_wp .or. estimate(i)*lambda(i) > 0.0_wp) lambda(i) = estimate(i)
        end do
        call lg%put(sqpopt_log_detail, 'multipliers re-estimated by least squares (the Hessian''s shift dominated '// &
                    'the last step): largest '//fmt_e(maxval(abs(lambda))))
        end subroutine estimate_multipliers

        subroutine convergence_tests()
        !! the convergence tests at the current point, setting `done` and `istat`
        !! (see [[check_convergence]]), and the count of consecutive stalled iterations
        call check_convergence(x, g, jac, c, problem%x_lb, problem%x_ub, problem%c_lb, problem%c_ub, &
                                lambda, options%ktol, options%ctol, done, istat, &
                                f=f, f_prev=f_prev, x_prev=x_prev, ftol=options%ftol, xtol=options%xtol, &
                                kkt_error=info%kkt, feas_error=info%feas, viol_prev=viol_prev, &
                                dual_inf_tol=options%dual_inf_tol, f_scale=problem%f_scale, stat_error=stat_err)
        info%stat_unscaled = stat_err/problem%f_scale
        ! the stalled-progress test must hold for `stall_iter` consecutive
        ! iterations: a single negligible step (e.g. a short line-search step on a
        ! badly scaled problem) is not a stall, and stopping on it made results
        ! depend on last-bit differences between platforms
        if (done .and. istat == sqpopt_stalled) then
            n_stalled = n_stalled + 1
            if (n_stalled < options%stall_iter) then
                done  = .false.
                istat = sqpopt_success
            end if
        else
            n_stalled = 0
        end if
        end subroutine convergence_tests

        subroutine solve_qp()
        !! solve the linearized QP subproblem for the search direction `p` and
        !! the multipliers `new_lambda` (see [[solve_qp_subproblem]]). Its
        !! direct method, if in use, may raise the exact Hessian's shift.
        real(wp) :: shift0
        real(wp), dimension(size(x)) :: lb, ub
        shift0 = hessian%shift
        ! (the variables' bounds, and the limits on their steps if any: see `set_max_step`)
        call problem%step_bounds(x, lb, ub)
        call qp_solver%solve(hessian, jac, x, g, c, lb, ub, &
                              problem%c_lb, problem%c_ub, p, new_lambda, qp_istat, kkt=kkt, inertia=inertia)
        if (hessian%shift > shift0) then
            info%hess_reset = .true.
            call lg%put(sqpopt_log_detail, 'direct QP: nonconvex face, Hessian shift '//fmt_e(shift0)//' -> '// &
                        fmt_e(hessian%shift))
        end if
        call note_qp()
        end subroutine solve_qp

        subroutine note_qp()
        !! the detailed log's line for the QP solve just done
        if (.not. lg%on(sqpopt_log_detail)) return
        if (qp_solver%unconstrained_used) then
            call lg%put(sqpopt_log_detail, 'QP: the unconstrained quasi-Newton step is feasible (no QP solver run)')
            return
        end if
        if (qp_solver%direct_used) then
            call lg%put(sqpopt_log_detail, 'direct QP: '//plural(qp_solver%n_iter, 'change', 'changes')// &
                        ' of the working set, working set '//fmt_i(qp_solver%n_working)//', ok')
            return
        end if
        if (qp_solver%direct_outcome >= 0) then
            call lg%put(sqpopt_log_detail, 'direct QP: gave up after '// &
                        plural(qp_solver%direct_changes, 'change', 'changes')//' of the working set ('// &
                        direct_outcome_text(qp_solver%direct_outcome)//')')
        end if
        call lg%put(sqpopt_log_detail, qp_solver%mode_name(problem%n)//': '// &
                    plural(qp_solver%n_iter, 'iteration', 'iterations')//', working set '//fmt_i(qp_solver%n_working)// &
                    ', '//plural(qp_solver%n_slacks, 'elastic slack', 'elastic slacks')//', '//qp_status_text(qp_istat)// &
                    trim(merge(', negative curvature', '                    ', qp_solver%negative_curvature)))
        end subroutine note_qp

        subroutine correct_inertia(predicted, shifted, ok)
        !! inertia control: increase the Hessian's shift, if necessary,
        !! until it has no negative curvature on the null space of the QP
        !! solver's working set (see [[inertia_correct]]): before the
        !! iteration's QP solve, the one that solve will start from (see
        !! [[starting_working_set]]), and after it, the one it ended with.
        logical, intent(in)  :: predicted !! whether this is before the iteration's QP solve (else after it)
        logical, intent(out) :: shifted   !! whether the shift was increased
        logical, intent(out) :: ok        !! whether no negative curvature is left
        integer, dimension(:), allocatable :: status
        real(wp) :: shift0
        integer  :: n_negative, n_factor0
        shifted = .false.
        ok      = .false.
        if (predicted) then
            call qp_solver%starting_working_set(problem%x_lb, problem%x_ub, problem%c_lb, problem%c_ub, status)
        else
            call qp_solver%working_set(problem%n, problem%m, status)
            if (.not. allocated(status)) return
        end if
        shift0    = hessian%shift
        n_factor0 = kkt%solver%n_factor
        call inertia%correct(kkt, hessian, jac, status, shifted, ok, n_negative)
        if (shifted) info%hess_reset = .true.
        if (.not. lg%on(sqpopt_log_detail)) return
        if (.not. inertia%enabled) then
            call lg%put(sqpopt_log_detail, 'inertia control: the factorization failed, continuing without it')
        else if (shifted) then
            call lg%put(sqpopt_log_detail, 'inertia control ('//trim(merge('starting', 'final   ', predicted))// &
                        ' working set): '//plural(n_negative, 'direction', 'directions')// &
                        ' of negative curvature, Hessian shift '//fmt_e(shift0)//' -> '//fmt_e(hessian%shift)//' ('// &
                        plural(kkt%solver%n_factor - n_factor0, 'factorization', 'factorizations')//')')
            if (.not. ok) call lg%put(sqpopt_log_detail, 'inertia control: negative curvature is left at the largest shift')
        end if
        end subroutine correct_inertia

        subroutine note_restoration_step(how)
        !! the detailed log's line for a single restoration step just taken (or not)
        character(len=*), intent(in) :: how !! which kind of step it was, for the message
        real(wp), dimension(problem%m) :: c_new
        if (.not. lg%on(sqpopt_log_detail)) return
        if (step_istat == sqpopt_success) then
            call problem%c(x_new, c_new)
            call lg%put(sqpopt_log_detail, 'restoration step ('//how//'): violation '// &
                        fmt_e(l1_violation(c, problem%c_lb, problem%c_ub))//' -> '// &
                        fmt_e(l1_violation(c_new, problem%c_lb, problem%c_ub)))
        else
            call lg%put(sqpopt_log_detail, 'restoration step ('//how//'): no decrease of the violation found')
        end if
        end subroutine note_restoration_step

        function glob_value() result(v)
        !! the globalization's state (see `sqpopt_iter_info%glob`)
        real(wp) :: v
        if (trust_region%enabled) then
            v = trust_region%radius
        else if (linesearch%mode == sqpopt_linesearch_filter) then
            v = 0.0_wp
            if (allocated(linesearch%filter%theta)) v = real(size(linesearch%filter%theta), wp)
        else if (linesearch%mode == sqpopt_linesearch_funnel) then
            v = linesearch%funnel%width
        else
            v = linesearch%merit%penalty
        end if
        end function glob_value

        pure function glob_name() result(name)
        !! what the restoration phase's exit test uses
        character(len=6) :: name
        name = merge('funnel', 'filter', linesearch%mode == sqpopt_linesearch_funnel)
        end function glob_name

        subroutine eval_f_cached(xx, ff)
        !! `f`, through the problem's evaluation cache (see [[sqpopt_problem_module]])
        real(wp), dimension(:), intent(in)  :: xx !! point `dimension(n)`
        real(wp),               intent(out) :: ff !! scaled objective at `xx`
        call problem%f(xx, ff)
        end subroutine eval_f_cached

        subroutine eval_c_cached(xx, cc)
        !! `c`, through the problem's evaluation cache (see [[sqpopt_problem_module]])
        real(wp), dimension(:), intent(in)  :: xx !! point `dimension(n)`
        real(wp), dimension(:), intent(out) :: cc !! scaled constraints at `xx` `dimension(m)`
        call problem%c(xx, cc)
        end subroutine eval_c_cached

        subroutine soc(p_trial, c_trial, p_soc, ok)
        !! the second-order correction of a rejected trial step, given to
        !! the line search (see [[sqpopt_soc_module]])
        real(wp), dimension(:), intent(in)  :: p_trial !! the rejected trial step
        real(wp), dimension(:), intent(in)  :: c_trial !! constraint values at `x+p_trial`
        real(wp), dimension(:), intent(out) :: p_soc   !! the corrected step
        logical,                intent(out) :: ok      !! true if `p_soc` is usable
        real(wp), dimension(size(x)) :: lb, ub
        call problem%step_bounds(x, lb, ub)
        call soc_step(jac, x, p_trial, c, c_trial, problem%c_lb, problem%c_ub, lb, ub, p_soc, ok, &
                      least_squares=least_squares)
        end subroutine soc

        subroutine start_restoration_phase(record)
        !! start a feasibility restoration phase at `x`, and take its first
        !! iteration. If `record`, first add `x` to the filter (or tighten the
        !! funnel toward it), so the iterations can't cycle back to it (a
        !! failed filter or funnel line search has already done so).
        logical, intent(in) :: record !! whether to first add `x` to the filter (or tighten the funnel toward it)
        real(wp) :: theta
        theta = l1_violation(c, problem%c_lb, problem%c_ub)
        if (record) then
            if (linesearch%mode == sqpopt_linesearch_funnel) then
                call linesearch%funnel%restoration(theta)
            else if (linesearch%mode == sqpopt_linesearch_filter) then
                call linesearch%filter%record(theta, f)
            end if
        end if
        call restoration%enter(x, theta, qp_solver)
        call lg%put(sqpopt_log_detail, 'restoration phase started, at violation '//fmt_e(theta))
        info%phase = .true.
        call restoration_phase_iteration()
        end subroutine start_restoration_phase

        subroutine restoration_phase_iteration()
        !! one iteration of the restoration phase (falling back to a
        !! Gauss-Newton step if it fails), and end the phase if the new point
        !! is good enough (see [[restoration_phase_done]]). If both steps
        !! fail, the phase also ends, so that the next iteration tries the
        !! optimality QP again instead of retrying the phase from the same point.
        real(wp) :: f_new, theta0, theta_new
        real(wp), dimension(problem%m) :: c_new
        logical :: ended
        character(len=:), allocatable :: how
        theta0 = l1_violation(c, problem%c_lb, problem%c_ub)
        how = 'feasibility QP'
        call restoration%step(problem, jac, x, c, x_new, alpha, step_istat)
        if (step_istat /= sqpopt_success) then
            how = 'Gauss-Newton fallback'
            call restoration_step(problem, jac, x, c, qp_solver%max_step*qp_solver%step_scale, x_new, alpha, step_istat, &
                                  least_squares=least_squares)
        end if
        if (step_istat == sqpopt_success) then
            call problem%f(x_new, f_new)
            call problem%c(x_new, c_new)
            theta_new = l1_violation(c_new, problem%c_lb, problem%c_ub)
            ended = restoration%done(linesearch, theta_new, f_new, &
                                     options%restoration_exit_factor, options%ctol, options%restoration_max_iter)
            if (lg%on(sqpopt_log_detail)) then
                call lg%put(sqpopt_log_detail, 'restoration phase iteration '//fmt_i(restoration%n_iter)//' ('//how// &
                            '): violation '//fmt_e(theta0)//' -> '//fmt_e(theta_new))
                if (ended) then
                    if (theta_new <= options%ctol) then
                        how = 'feasible'
                    else if (restoration%n_iter >= options%restoration_max_iter) then
                        how = 'iteration limit'
                    else
                        how = 'violation reduced, and acceptable to the '//trim(glob_name())
                    end if
                    call lg%put(sqpopt_log_detail, 'restoration phase ended ('//how//')')
                end if
            end if
        else
            restoration%active = .false.
            call lg%put(sqpopt_log_detail, 'restoration phase ended (no step reduces the violation)')
        end if
        end subroutine restoration_phase_iteration

        subroutine elastic_resolve()
        !! if a constraint's multiplier is diverging, re-solve the QP with that
        !! constraint elastic (SNOPT-style elastic mode). At a point where the
        !! constraint qualification fails (e.g. a constraint tangent to a
        !! bound), the linearization lets each step cover only a fraction of
        !! the distance to it, so the iterates creep toward the point with a
        !! multiplier that grows every iteration, and may converge there to a
        !! point that is not a minimizer. An \( \ell_1 \) penalty with a
        !! bounded weight on the constraint (the elastic slack) caps its
        !! multiplier, and lets the step relax the linearization and leave.
        !!
        !! A constraint `i` triggers it when its multiplier's "push"
        !! \( |\lambda_i| \lVert \nabla c_i \rVert_\infty \) exceeds
        !! `options%elastic_multiplier_limit` \( \times \max(1, \lVert g \rVert_\infty) \)
        !! both now and at the previous iterate, and has grown by more than
        !! `growth` since then. The re-solve is done at most `max_resolves`
        !! times per solve: at a solution that itself has an unbounded
        !! multiplier (e.g. a cusp), it would otherwise keep pushing the
        !! iterates away.
        real(wp), parameter :: growth       = 1.5_wp !! the multiplier growth that indicates divergence
        integer,  parameter :: max_resolves = 3      !! maximum elastic re-solves per solve
        real(wp) :: lim, rownorm(problem%m), wmax
        real(wp), dimension(problem%n) :: lb_step, ub_step
        integer  :: sgn(problem%m), k, i
        if (options%elastic_multiplier_limit <= 0.0_wp .or. problem%m == 0 .or. &
            qp_solver%n_elastic >= max_resolves) return
        lim = options%elastic_multiplier_limit*max(1.0_wp, maxval(abs(g)))
        rownorm = 0.0_wp
        do k = 1, jac%nnz
            rownorm(jac%irow(k)) = max(rownorm(jac%irow(k)), abs(jac%val(k)))
        end do
        sgn  = 0
        wmax = 0.0_wp
        do i = 1, problem%m
            if (abs(new_lambda(i))*rownorm(i) > lim .and. abs(lambda(i))*rownorm(i) > lim .and. &
                abs(new_lambda(i)) > growth*abs(lambda(i))) then
                sgn(i) = merge(1, -1, new_lambda(i) > 0.0_wp)   ! (the side of the row that is active)
                wmax   = max(wmax, rownorm(i))
            end if
        end do
        if (all(sgn == 0)) return
        qp_solver%n_elastic = qp_solver%n_elastic + 1
        info%elastic = .true.
        call lg%put(sqpopt_log_detail, 'elastic re-solve: '// &
                    plural(count(sgn /= 0), 'constraint', 'constraints')//' with a diverging multiplier, weight '// &
                    fmt_e(lim/wmax))
        ! (the weight caps each elastic row's push at about the limit)
        call problem%step_bounds(x, lb_step, ub_step)
        call qp_solver%solve(hessian, jac, x, g, c, lb_step, ub_step, &
                              problem%c_lb, problem%c_ub, p, new_lambda, qp_istat, &
                              elastic_sign=sgn, elastic_weight=lim/wmax)
        call note_qp()
        restore = qp_istat == sqpopt_infeasible
        end subroutine elastic_resolve

        subroutine update_penalty()
        !! update the merit function's penalty parameter for the QP step `p`
        !! and multipliers `new_lambda` (see [[update_penalty_parameter]])
        real(wp), dimension(problem%n) :: hp
        call hessian%hv_product(p, hp)
        call linesearch%merit%update_penalty(jac, g, p, dot_product(p, hp), c, problem%c_lb, problem%c_ub, &
                                       lambda, new_lambda)
        end subroutine update_penalty

    end subroutine sqpopt_iterate
!*******************************************************************************

!*******************************************************************************
!>
!  evaluate the (scaled) problem at `x`, and the resulting KKT and
!  feasibility errors and variable-bound multipliers `z` (with the
!  constraint multipliers `lambda`): \( z = g - J^T\lambda \) for a variable
!  at one of its bounds, zero otherwise (for the Lagrangian
!  \( f - \lambda^Tc - z^Tx \)). Used to report the final state of a solve.

    subroutine sqpopt_evaluate_point(problem, options, x, lambda, jac, f, c, kkt, feas, z, stat_error)

    type(sqpopt_problem_type),  intent(inout) :: problem    !! problem definition
    type(sqpopt_options_type),  intent(in)    :: options    !! solver options
    real(wp), dimension(:),     intent(in)    :: x          !! point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: lambda     !! constraint multipliers `dimension(m)`
    type(sqpopt_sparse_matrix), intent(inout) :: jac        !! Jacobian workspace (as for [[sqpopt_iterate]])
    real(wp),                   intent(out)   :: f          !! objective at `x`
    real(wp), dimension(:),     intent(out)   :: c          !! constraints at `x` `dimension(m)`
    real(wp),                   intent(out)   :: kkt        !! KKT error at `x`
    real(wp),                   intent(out)   :: feas       !! feasibility error at `x`
    real(wp), dimension(:),     intent(out)   :: z          !! variable-bound multipliers `dimension(n)`
    real(wp), optional,         intent(out)   :: stat_error !! the stationarity residual (of the scaled problem,
                                                             !! without the multiplier scaling; see [[check_convergence]])

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
                            kkt_error=kkt, feas_error=feas, stat_error=stat_error)

    call sparse_matvec_transpose(jac, lambda, jtlam)
    z = 0.0_wp
    do j = 1, size(x)
        if (x(j) - problem%x_lb(j) <= options%ctol .or. problem%x_ub(j) - x(j) <= options%ctol) then
            z(j) = g(j) - jtlam(j)
        end if
    end do

    end subroutine sqpopt_evaluate_point
!*******************************************************************************


    end module sqpopt_iterate_module
!*******************************************************************************
