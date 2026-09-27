!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Module containing the main object-oriented interface to the `sqpopt`
!  Sequential Quadratic Programming solver. It is used via the
!  [[sqpopt_type]] class, which is the only public entity in this module.
!  The other components of the algorithm (problem definition, options,
!  Hessian approximation, QP subproblem solver, line search, and
!  convergence checking) are each implemented in their own module so
!  that they may be developed, tested, and swapped out independently.
!  Internally, sparse (COO) storage is used by default for the
!  constraint Jacobian and the Lagrangian Hessian is a matrix-free
!  limited-memory operator -- dense `n x n`/`m x n` arrays are never formed.

    module sqpopt_module

    use, intrinsic :: iso_fortran_env, only: int64
    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_types_module,      only: sqpopt_success, sqpopt_max_iter_reached, &
                                         sqpopt_user_requested_stop, sqpopt_report_func, &
                                         sqpopt_invalid_input, sqpopt_status_message, sqpopt_sparse_matrix, &
                                         sqpopt_results_type, sqpopt_max_evals_reached, sqpopt_time_limit_reached, &
                                         sqpopt_qp_solve_failed, sqpopt_infinity, sqpopt_all_finite
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type, sqpopt_hessian_bfgs, sqpopt_hessian_sr1, sqpopt_hessian_exact
    use sqpopt_hessian_module,    only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type, sqpopt_qp_auto, sqpopt_qp_dense, &
                                         sqpopt_qp_reduced_hessian
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_null_space_lu, sqpopt_null_space_lsqr
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_armijo, sqpopt_linesearch_funnel, &
                                         sqpopt_merit_l1, sqpopt_merit_augmented_lagrangian, &
                                         sqpopt_penalty_multipliers, sqpopt_penalty_model
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_restoration_module,  only: sqpopt_restoration_type, sqpopt_restoration_phase, &
                                          sqpopt_restoration_gauss_newton
    use sqpopt_iterate_module,    only: sqpopt_iterate, sqpopt_evaluate_point, sqpopt_iter_info

    implicit none

    private

    type, public :: sqpopt_type
        !! main class for the SQP optimizer.

        private

        type(sqpopt_problem_type)    :: problem      !! the nonlinear program to be solved (working copy)
        type(sqpopt_options_type)    :: options      !! solver options
        type(sqpopt_hessian_type)    :: hessian      !! Hessian of the Lagrangian approximation
        type(sqpopt_qp_solver_type)  :: qp_solver    !! QP subproblem solver
        type(sqpopt_linesearch_type) :: linesearch   !! merit function / line search
        type(sqpopt_trust_region_type) :: trust_region !! trust-region globalization (opt-in alternative to `linesearch`)
        type(sqpopt_restoration_type)  :: restoration  !! feasibility restoration phase state (reset on each `solve`)

        ! the components as given to `initialize`, restored at the start of every
        ! `solve`, so that no state (scaling, caches, penalty, filter, watchdog,
        ! trust-region radius, QP warm start) carries over from a previous solve:
        type(sqpopt_problem_type)      :: problem0
        type(sqpopt_linesearch_type)   :: linesearch0
        type(sqpopt_trust_region_type) :: trust_region0
        type(sqpopt_qp_solver_type)    :: qp_solver0

        real(wp), dimension(:), allocatable :: x       !! current/final optimization variables
        real(wp), dimension(:), allocatable :: lambda  !! current/final Lagrange multipliers (of the scaled problem)

        procedure(sqpopt_report_func), pointer, nopass :: report => null() !! optional user progress-reporting
                                                                           !! callback, invoked once per major
                                                                           !! iteration (see [[sqpopt_types_module]])

        type(sqpopt_results_type) :: results !! the outcome of the last `solve`

        contains

        private

        procedure, public :: initialize     => sqpopt_initialize
        procedure, public :: solve          => sqpopt_solve
        procedure, public :: get_solution   => sqpopt_get_solution
        procedure, public :: get_results    => sqpopt_get_results
        procedure, public :: status_message => sqpopt_get_status_message
        procedure, public :: destroy        => sqpopt_destroy

    end type sqpopt_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize (or reinitialize) an [[sqpopt_type]] solver instance.

    subroutine sqpopt_initialize(me, problem, options, hessian, qp_solver, linesearch, trust_region, report)

    class(sqpopt_type), intent(inout) :: me
    type(sqpopt_problem_type),optional,intent(in)    :: problem      !! the nonlinear program to be solved
    type(sqpopt_options_type),optional,intent(in)    :: options      !! solver options
    type(sqpopt_hessian_type),optional,intent(in)    :: hessian      !! Hessian of the Lagrangian approximation
    type(sqpopt_qp_solver_type),optional,intent(in)  :: qp_solver    !! QP subproblem solver
    type(sqpopt_linesearch_type),optional,intent(in) :: linesearch   !! merit function / line search
    type(sqpopt_trust_region_type),optional,intent(in) :: trust_region !! trust-region globalization (opt-in
                                                                        !! alternative to `linesearch`, see
                                                                        !! [[sqpopt_trust_region_module]])
    procedure(sqpopt_report_func), optional :: report               !! optional user progress-reporting callback,
                                                                     !! called once per major iteration with the
                                                                     !! current iterate; set its `user_stop` output
                                                                     !! to request the solver stop early
                                                                     !! (see [[sqpopt_types_module]])

    ! use inputs if present, else use defaults:
    if (present(problem))    then; me%problem0 = problem; else; me%problem0 = sqpopt_problem_type(); end if
    if (present(options))    then; me%options = options; else; me%options = sqpopt_options_type(); end if
    if (present(hessian))    then; me%hessian = hessian; else; me%hessian = sqpopt_hessian_type(); end if
    if (present(qp_solver))  then; me%qp_solver0 = qp_solver; else; me%qp_solver0 = sqpopt_qp_solver_type(); end if
    if (present(linesearch)) then; me%linesearch0 = linesearch; else; me%linesearch0 = sqpopt_linesearch_type(); end if
    if (present(trust_region)) then
        me%trust_region0 = trust_region
    else
        me%trust_region0 = sqpopt_trust_region_type()
    end if
    me%report => null()
    if (present(report)) me%report => report

    if (allocated(me%x))      deallocate(me%x)
    if (allocated(me%lambda)) deallocate(me%lambda)
    me%results = sqpopt_results_type()
    me%results%message = ''

    end subroutine sqpopt_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  solve the nonlinear program starting from the initial guess `x0` (and,
!  optionally, initial constraint multipliers `lambda0`, with the sign
!  convention of [[sqpopt_results_type]]), running major SQP iterations
!  until convergence or a stopping criterion is met. The outcome is
!  available from `get_solution`, `get_results`, and `status_message`.

    subroutine sqpopt_solve(me, x0, istat, lambda0)

    class(sqpopt_type),     intent(inout) :: me
    real(wp), dimension(:), intent(in)    :: x0      !! initial guess for the optimization variables `dimension(n)`
    integer,                intent(out)   :: istat   !! status code (see [[sqpopt_types_module]])
    real(wp), dimension(:), intent(in), optional :: lambda0 !! initial guess for the constraint multipliers `dimension(m)`

    real(wp), dimension(:), allocatable :: x_prev, gl_prev  !! quasi-Newton state (unallocated until the 2nd iteration)
    real(wp), allocatable :: f_prev  !! previous objective value, for the `options%ftol` stalled-progress test (unallocated until the 2nd iteration)
    real(wp), allocatable :: viol_prev !! previous constraint violation, for the infeasibility test (unallocated until the 2nd iteration)
    type(sqpopt_sparse_matrix) :: jac !! Jacobian workspace (structure set once, values updated each iteration)
    type(sqpopt_iter_info) :: info
    logical :: done, valid
    integer :: iter_istat, iter, n_fail, n_acceptable, n_stalled, n_escape
    integer(int64) :: t_start, t_now, t_rate
    character(len=:), allocatable :: msg
    type(sqpopt_restoration_type) :: fresh_restoration !! (default-initialized)

    call system_clock(t_start, t_rate)
    valid = .false.

    ! start every solve from the components exactly as configured:
    me%problem      = me%problem0
    me%linesearch   = me%linesearch0
    me%trust_region = me%trust_region0
    me%restoration  = fresh_restoration
    me%qp_solver    = me%qp_solver0

    me%x = x0
    if (allocated(me%lambda)) deallocate(me%lambda)
    allocate(me%lambda(max(me%problem%m,0)))
    me%lambda = 0.0_wp
    me%results%iterations = 0

    ! check the inputs before doing anything else:
    call me%problem%validate(istat, msg)
    if (istat == sqpopt_success) call validate_options(me, istat, msg)
    if (istat == sqpopt_success .and. size(x0) /= me%problem%n) then
        istat = sqpopt_invalid_input
        msg   = 'x0 must have size n'
    end if
    if (istat == sqpopt_success .and. .not. sqpopt_all_finite(x0)) then
        ! (checked before `x0` is projected onto the bounds below, which would
        ! silently turn a NaN into a bound)
        istat = sqpopt_invalid_input
        msg   = 'x0 must be finite'
    end if
    if (istat == sqpopt_success .and. present(lambda0)) then
        if (size(lambda0) /= me%problem%m) then
            istat = sqpopt_invalid_input
            msg   = 'lambda0 must have size m'
        else if (.not. sqpopt_all_finite(lambda0)) then
            istat = sqpopt_invalid_input
            msg   = 'lambda0 must be finite'
        end if
    end if
    if (istat /= sqpopt_success) then
        call finish(istat, msg)
        return
    end if
    valid = .true.

    ! start from a point that satisfies the variable bounds, so the user
    ! functions are never evaluated outside them:
    me%x = min(max(x0, me%problem%x_lb), me%problem%x_ub)

    ! the QP step-length cap starts relative to the size of the variables (so
    ! a solution far away isn't approached by doubling the cap from
    ! `max_step`; see [[sqpopt_qp_solver_module]]):
    me%qp_solver%step_scale = max(1.0_wp, maxval(abs(me%x))/me%qp_solver%max_step)

    ! empty evaluation caches, and (optionally) gradient-based scaling:
    call me%problem%reset_evaluations()
    if (me%options%scaling) call me%problem%compute_scaling(me%x, me%options%scaling_max_gradient)
    if (present(lambda0)) me%lambda = lambda0*me%problem%f_scale/me%problem%c_scale

    call me%hessian%initialize(me%problem%n, lbfgs_memory(me%options%lbfgs_memory, me%problem%n), &
                                use_sr1=(me%options%hessian_mode == sqpopt_hessian_sr1), &
                                scale0=me%options%hessian_scale0)
    if (me%options%hessian_mode == sqpopt_hessian_exact) then
        call me%hessian%set_exact(me%problem%hess_irow, me%problem%hess_icol)
    end if
    me%qp_solver%mode        = me%options%qp_solver_mode
    me%linesearch%mode       = me%options%linesearch_mode
    me%linesearch%merit%mode           = me%options%merit_mode
    me%linesearch%merit%penalty_update = me%options%penalty_update

    if (me%options%print_level >= 1) call print_header()

    n_fail = 0
    n_acceptable = 0
    n_stalled    = 0
    n_escape     = 0
    do iter = 1, me%options%max_iter
        me%results%iterations = iter
        call sqpopt_iterate(me%problem, me%options, me%hessian, me%qp_solver, me%linesearch, me%trust_region, &
                             me%x, me%lambda, x_prev, gl_prev, f_prev, viol_prev, jac, n_acceptable, n_stalled, n_escape, &
                             me%restoration, iter, me%report, &
                             done, iter_istat, info)
        if (me%options%print_level >= 1) call print_iteration(iter, info, iter_istat)
        if (done) then
            ! converged, stalled, acceptable, infeasible, unbounded, function error, or user stop:
            call finish(iter_istat)
            return
        end if
        if (iter_istat == sqpopt_success) then
            n_fail = 0
        else
            ! a failed QP solve or line search: tolerate a few in a row, since
            ! the next iteration's re-linearization and Hessian update often
            ! recover, but don't loop until `max_iter` on a stuck iteration:
            n_fail = n_fail + 1
            if (n_fail >= me%options%max_consecutive_failures) then
                call finish(iter_istat)
                return
            end if
        end if
        ! evaluation and time limits:
        if (me%options%max_evals > 0 .and. me%problem%n_eval_fc >= me%options%max_evals) then
            call finish(sqpopt_max_evals_reached)
            return
        end if
        if (me%options%max_time > 0.0_wp) then
            call system_clock(t_now)
            if (real(t_now-t_start, wp)/real(t_rate, wp) >= me%options%max_time) then
                call finish(sqpopt_time_limit_reached)
                return
            end if
        end if
    end do

    call finish(sqpopt_max_iter_reached)

    contains

        subroutine finish(stat, detail)
        !! set the final status and message, fill in the results (for the
        !! original, unscaled problem), and print the summary
        integer, intent(in) :: stat
        character(len=*), intent(in), optional :: detail !! extra detail appended to the status message
        real(wp) :: fs
        real(wp), dimension(me%problem%m) :: cs
        real(wp), dimension(size(me%x)) :: zs

        istat = stat
        me%results%istat   = stat
        me%results%message = sqpopt_status_message(stat)
        if (present(detail)) then
            if (len(detail) > 0) me%results%message = me%results%message//': '//detail
        end if

        me%results%x = me%x
        if (valid) then
            call sqpopt_evaluate_point(me%problem, me%options, me%x, me%lambda, jac, fs, cs, &
                                       me%results%kkt_error, me%results%feasibility_error, zs)
            me%results%f      = fs/me%problem%f_scale
            me%results%c      = cs/me%problem%c_scale
            me%results%lambda = me%lambda*me%problem%c_scale/me%problem%f_scale
            me%results%z      = zs/me%problem%f_scale
            ! (the feasibility error of the original problem:)
            me%results%feasibility_error = 0.0_wp
            if (me%problem%m > 0) me%results%feasibility_error = maxval(max(me%problem%c_lb/me%problem%c_scale &
                - me%results%c, 0.0_wp) + max(me%results%c - me%problem%c_ub/me%problem%c_scale, 0.0_wp))
            me%results%n_eval_fc   = me%problem%n_eval_fc
            me%results%n_eval_gjac = me%problem%n_eval_gjac
            me%results%n_eval_hess = me%problem%n_eval_hess
        else
            me%results%f = 0.0_wp
            me%results%c = spread(0.0_wp, 1, max(me%problem%m,0))
            me%results%lambda = me%lambda
            me%results%z = spread(0.0_wp, 1, size(x0))
            me%results%n_eval_fc = 0; me%results%n_eval_gjac = 0; me%results%n_eval_hess = 0
        end if
        call system_clock(t_now)
        me%results%time = real(t_now-t_start, wp)/real(t_rate, wp)

        ! (not for invalid inputs, which may include `output_unit` itself)
        if (valid .and. me%options%print_level >= 1) call print_summary()
        end subroutine finish

        subroutine print_header()
        !! the iteration log's column headings
        integer :: u
        u = me%options%output_unit
        write(u,'(A)') ''
        write(u,'(A,I0,A,I0,A,ES9.2,A,ES9.2)') ' sqpopt: n = ', me%problem%n, ', m = ', me%problem%m, &
            ', objective scale = ', me%problem%f_scale, ', min constraint scale = ', &
            minval([1.0_wp, me%problem%c_scale])
        if (me%options%print_level >= 2) then
            write(u,'(A6,A18,3A11,A11,A11,A7,A6)') 'iter', 'objective', 'infeas', 'kkt', 'alpha', &
                'penalty', '|step|', 'qp_it', 'flags'
        else
            write(u,'(A6,A18,3A11,A6)') 'iter', 'objective', 'infeas', 'kkt', 'alpha', 'flags'
        end if
        end subroutine print_header

        subroutine print_iteration(iter, info, iter_istat)
        !! one line of the iteration log. Flags: `R` = feasibility
        !! restoration step, `Q` = QP solve failed, `F` = no acceptable step
        !! found (the Hessian approximation is reset).
        integer,                intent(in) :: iter, iter_istat
        type(sqpopt_iter_info), intent(in) :: info
        character(len=3) :: flags
        integer :: u
        u = me%options%output_unit
        flags = ''
        if (info%restoration) flags = trim(flags)//'R'
        if (info%qp_istat == sqpopt_qp_solve_failed) flags = trim(flags)//'Q'
        if (info%stepped .and. iter_istat /= sqpopt_success .and. iter_istat /= sqpopt_qp_solve_failed) then
            flags = trim(flags)//'F'
        end if
        if (.not. info%stepped) then
            write(u,'(I6,ES18.9,2ES11.2)') iter, info%f/me%problem%f_scale, info%feas, info%kkt
        else if (me%options%print_level >= 2) then
            write(u,'(I6,ES18.9,5ES11.2,I7,2X,A)') iter, info%f/me%problem%f_scale, info%feas, info%kkt, &
                info%alpha, info%penalty, info%step_norm, info%qp_iter, flags
        else
            write(u,'(I6,ES18.9,3ES11.2,2X,A)') iter, info%f/me%problem%f_scale, info%feas, info%kkt, &
                info%alpha, flags
        end if
        end subroutine print_iteration

        subroutine print_summary()
        !! the final summary
        integer :: u
        u = me%options%output_unit
        write(u,'(A)') ''
        write(u,'(A,I0,2A)')   ' sqpopt: status ', me%results%istat, ': ', me%results%message
        write(u,'(A,ES18.10)') '   objective          = ', me%results%f
        write(u,'(A,ES10.2)')  '   feasibility error  = ', me%results%feasibility_error
        write(u,'(A,ES10.2)')  '   KKT error (scaled) = ', me%results%kkt_error
        write(u,'(A,I0)')      '   iterations         = ', me%results%iterations
        if (me%results%n_eval_hess > 0) then
            write(u,'(A,3(I0,A))') '   evaluations        = ', me%results%n_eval_fc, ' fc, ', me%results%n_eval_gjac, &
                                   ' gjac, ', me%results%n_eval_hess, ' hess'
        else
            write(u,'(A,2(I0,A))') '   evaluations        = ', me%results%n_eval_fc, ' fc, ', me%results%n_eval_gjac, ' gjac'
        end if
        write(u,'(A,F0.3,A)')  '   time               = ', me%results%time, ' s'
        write(u,'(A)') ''
        end subroutine print_summary

    end subroutine sqpopt_solve
!*******************************************************************************

!*******************************************************************************
!>
!  check the options and every component's settings for invalid values.

    subroutine validate_options(me, istat, msg)

    class(sqpopt_type),            intent(in)  :: me
    integer,                       intent(out) :: istat !! `sqpopt_success` or `sqpopt_invalid_input`
    character(len=:), allocatable, intent(out) :: msg   !! description of the problem found (empty if none)

    logical :: opened

    istat = sqpopt_invalid_input
    msg   = ''

    associate (o => me%options, ls => me%linesearch, tr => me%trust_region, qp => me%qp_solver)

    ! ---- options ----
    if (o%max_iter < 0) then
        msg = 'options%max_iter must be >= 0'
        return
    end if
    if (.not. (o%elastic_multiplier_limit >= 0.0_wp)) then
        msg = 'options%elastic_multiplier_limit must be >= 0 (0: disabled)'; return
    end if
    if (o%lbfgs_memory < 0) then
        msg = 'options%lbfgs_memory must be >= 0 (0: automatic)'
        return
    end if
    if (o%max_consecutive_failures < 1) then
        msg = 'options%max_consecutive_failures must be >= 1'
        return
    end if
    if (o%hessian_mode < sqpopt_hessian_bfgs .or. o%hessian_mode > sqpopt_hessian_exact) then
        msg = 'options%hessian_mode is not a valid sqpopt_hessian_* value'
        return
    end if
    if (o%hessian_mode == sqpopt_hessian_exact .and. .not. associated(me%problem%eval_hess)) then
        msg = 'options%hessian_mode = sqpopt_hessian_exact requires the hess function (set_functions) '// &
              'and its sparsity pattern (set_hessian_sparsity)'
        return
    end if
    if (all(o%qp_solver_mode /= [sqpopt_qp_auto, sqpopt_qp_dense, sqpopt_qp_reduced_hessian])) then
        msg = 'options%qp_solver_mode is not a valid sqpopt_qp_* value'
        return
    end if
    if (o%linesearch_mode < sqpopt_linesearch_armijo .or. o%linesearch_mode > sqpopt_linesearch_funnel) then
        msg = 'options%linesearch_mode is not a valid sqpopt_linesearch_* value'
        return
    end if
    if (o%merit_mode < sqpopt_merit_l1 .or. o%merit_mode > sqpopt_merit_augmented_lagrangian) then
        msg = 'options%merit_mode is not a valid sqpopt_merit_* value'
        return
    end if
    if (all(o%penalty_update /= [sqpopt_penalty_multipliers, sqpopt_penalty_model])) then
        msg = 'options%penalty_update is not a valid sqpopt_penalty_* value'
        return
    end if
    if (.not. (o%ktol > 0.0_wp .and. o%ctol > 0.0_wp)) then
        msg = 'options%ktol and options%ctol must be > 0'
        return
    end if
    if (.not. (o%dual_inf_tol > 0.0_wp)) then
        msg = 'options%dual_inf_tol must be > 0'
        return
    end if
    if (.not. (o%ftol >= 0.0_wp .and. o%xtol >= 0.0_wp)) then
        msg = 'options%ftol and options%xtol must be >= 0'
        return
    end if
    if (.not. (o%acceptable_ktol > 0.0_wp .and. o%acceptable_ctol > 0.0_wp) .or. o%acceptable_iter < 0 &
        .or. o%stall_iter < 1) then
        msg = 'options%acceptable_ktol/acceptable_ctol must be > 0, acceptable_iter >= 0, and stall_iter >= 1'
        return
    end if
    if (all(o%restoration_mode /= [sqpopt_restoration_phase, sqpopt_restoration_gauss_newton]) .or. &
        .not. (o%restoration_exit_factor > 0.0_wp .and. o%restoration_exit_factor < 1.0_wp) .or. &
        o%restoration_max_iter < 1) then
        msg = 'options%restoration_mode must be a sqpopt_restoration_* value, restoration_exit_factor in (0,1), '// &
              'and restoration_max_iter >= 1'
        return
    end if
    if (o%max_evals < 0 .or. .not. (o%max_time >= 0.0_wp)) then
        msg = 'options%max_evals and options%max_time must be >= 0'
        return
    end if
    if (.not. (o%scaling_max_gradient > 0.0_wp .and. o%hessian_scale0 > 0.0_wp)) then
        msg = 'options%scaling_max_gradient and options%hessian_scale0 must be > 0'
        return
    end if
    if (o%print_level > 0) then
        inquire(unit=o%output_unit, opened=opened)
        if (.not. opened) then
            msg = 'options%output_unit is not an open unit'
            return
        end if
    end if

    ! ---- line search ----
    if (.not. (ls%sigma > 0.0_wp .and. ls%sigma < 1.0_wp .and. ls%backtrack > 0.0_wp .and. ls%backtrack < 1.0_wp)) then
        msg = 'linesearch%sigma and linesearch%backtrack must be in (0,1)'
        return
    end if
    if (.not. (ls%alpha_min > 0.0_wp .and. ls%tol > 0.0_wp .and. ls%major_step_limit > 0.0_wp) &
        .or. ls%max_ls_iter < 1) then
        msg = 'linesearch%alpha_min, tol, and major_step_limit must be > 0, and max_ls_iter >= 1'
        return
    end if
    if (.not. (ls%merit%penalty >= 0.0_wp) .or. ls%watchdog_relaxed_len < 0 .or. ls%watchdog_cooldown_len < 0) then
        msg = 'linesearch%merit%penalty, watchdog_relaxed_len, and watchdog_cooldown_len must be >= 0'
        return
    end if
    if (ls%nonmonotone_len < 0 .or. .not. (ls%merit%penalty_rho > 0.0_wp .and. ls%merit%penalty_rho < 1.0_wp)) then
        msg = 'linesearch%nonmonotone_len must be >= 0, and merit%penalty_rho in (0,1)'
        return
    end if
    if (.not. (ls%filter%gamma_theta > 0.0_wp .and. ls%filter%gamma_theta < 1.0_wp .and. &
               ls%filter%gamma_phi > 0.0_wp .and. ls%filter%delta > 0.0_wp .and. &
               ls%filter%s_theta > 1.0_wp .and. ls%filter%s_phi > 1.0_wp .and. &
               ls%filter%eta_phi > 0.0_wp .and. ls%filter%eta_phi < 0.5_wp .and. &
               ls%filter%theta_min_fact > 0.0_wp .and. ls%filter%theta_max_fact > ls%filter%theta_min_fact .and. &
               ls%filter%gamma_alpha > 0.0_wp .and. ls%filter%gamma_alpha <= 1.0_wp)) then
        msg = 'a linesearch%filter parameter is out of range'
        return
    end if
    if (.not. (ls%funnel%width_min > 0.0_wp .and. ls%funnel%width_fact >= 1.0_wp .and. &
               ls%funnel%beta > 0.0_wp .and. ls%funnel%beta < 1.0_wp .and. &
               ls%funnel%kappa > 0.0_wp .and. ls%funnel%kappa < 1.0_wp .and. &
               ls%funnel%delta > 0.0_wp .and. ls%funnel%s_theta > 1.0_wp .and. &
               ls%funnel%eta > 0.0_wp .and. ls%funnel%eta < 0.5_wp .and. ls%funnel%gamma > 0.0_wp) &
        .or. all(ls%funnel%update /= [1, 2])) then
        msg = 'a linesearch%funnel parameter is out of range'
        return
    end if

    ! ---- trust region ----
    if (tr%enabled) then
        if (.not. (tr%radius0 > 0.0_wp .and. tr%radius_min > 0.0_wp .and. tr%radius_max >= tr%radius0 .and. &
                   tr%radius0 >= tr%radius_min)) then
            msg = 'trust_region radii must satisfy 0 < radius_min <= radius0 <= radius_max'
            return
        end if
        if (.not. (tr%eta1 > 0.0_wp .and. tr%eta1 <= tr%eta2 .and. tr%eta2 < 1.0_wp .and. &
                   tr%shrink_factor > 0.0_wp .and. tr%shrink_factor < 1.0_wp .and. tr%expand_factor >= 1.0_wp) &
            .or. tr%max_retries < 1) then
            msg = 'trust_region: need 0 < eta1 <= eta2 < 1, 0 < shrink_factor < 1, expand_factor >= 1, '// &
                  'and max_retries >= 1'
            return
        end if
    end if

    ! ---- QP solver ----
    if (.not. qp%max_step > 0.0_wp .or. qp%auto_dense_max_n < 0) then
        msg = 'qp_solver%max_step must be > 0, and auto_dense_max_n >= 0'
        return
    end if
    associate (d => qp%dense_qp, r => qp%sparse_qp)
    if (.not. (d%active_tol > 0.0_wp .and. d%opt_tol > 0.0_wp .and. d%feas_tol > 0.0_wp .and. &
               d%elastic_weight > 0.0_wp .and. d%elastic_weight_max >= d%elastic_weight) .or. d%max_iter < 1) then
        msg = 'a qp_solver%dense_qp setting is out of range'
        return
    end if
    if (.not. (r%active_tol > 0.0_wp .and. r%opt_tol > 0.0_wp .and. r%feas_tol > 0.0_wp .and. &
               r%pcg_rtol > 0.0_wp .and. r%elastic_weight > 0.0_wp .and. r%elastic_weight_max >= r%elastic_weight &
               .and. r%lsqr_atol >= 0.0_wp .and. r%lsqr_btol >= 0.0_wp .and. r%lsqr_conlim >= 0.0_wp) &
        .or. r%max_iter < 1 .or. r%dense_max_ns < 0 .or. all(r%null_space /= [sqpopt_null_space_lu, sqpopt_null_space_lsqr])) then
        msg = 'a qp_solver%sparse_qp setting is out of range'
        return
    end if
    end associate

    end associate

    istat = sqpopt_success

    end subroutine validate_options
!*******************************************************************************

!*******************************************************************************
!>
!  a description of the status of the last `solve` (including, for
!  `sqpopt_invalid_input`, what was invalid).

    function sqpopt_get_status_message(me) result(msg)

    class(sqpopt_type), intent(in) :: me
    character(len=:), allocatable :: msg

    if (allocated(me%results%message)) then
        msg = me%results%message
    else
        msg = ''
    end if

    end function sqpopt_get_status_message
!*******************************************************************************

!*******************************************************************************
!>
!  return the final solution of the last `solve`, the constraint
!  multipliers, and (optionally) the variable-bound multipliers (for the
!  original problem; see [[sqpopt_results_type]] for the sign convention).

    subroutine sqpopt_get_solution(me, x, lambda, z)

    class(sqpopt_type),     intent(in)  :: me
    real(wp), dimension(:), intent(out) :: x       !! optimization variables `dimension(n)`
    real(wp), dimension(:), intent(out) :: lambda  !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out), optional :: z !! variable-bound multipliers `dimension(n)`

    x      = me%results%x
    lambda = me%results%lambda
    if (present(z)) z = me%results%z

    end subroutine sqpopt_get_solution
!*******************************************************************************

!*******************************************************************************
!>
!  return the full results of the last `solve` (status, counts, final
!  point, multipliers, and errors; see [[sqpopt_results_type]]).

    subroutine sqpopt_get_results(me, results)

    class(sqpopt_type),        intent(in)  :: me
    type(sqpopt_results_type), intent(out) :: results

    results = me%results

    end subroutine sqpopt_get_results
!*******************************************************************************

!*******************************************************************************
!>
!  destroy the solver instance, deallocating all internal arrays.

    subroutine sqpopt_destroy(me)

    class(sqpopt_type), intent(inout) :: me

    call me%initialize()

    end subroutine sqpopt_destroy
!*******************************************************************************

!*******************************************************************************
!>
!  the number of `(s,y)` pairs the limited-memory Hessian keeps, for the
!  option value `memory` (see `sqpopt_options_type%lbfgs_memory`) and `n`
!  variables: `memory` itself, or, if it is `0` (automatic),
!  \( \max(10, \min(n, 100)) \).

    pure integer function lbfgs_memory(memory, n)

    integer, intent(in) :: memory !! the option value
    integer, intent(in) :: n      !! number of variables

    if (memory > 0) then
        lbfgs_memory = memory
    else
        lbfgs_memory = max(10, min(n, 100))
    end if

    end function lbfgs_memory
!*******************************************************************************

    end module sqpopt_module
!*******************************************************************************
