!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Module containing the main object-oriented interface to the `sqpopt`
!  Sequential Quadratic Programming solver. It is used via the
!  [[sqpopt_type]] class, which is the only public entity in this module.
!  The other components of the algorithm (problem definition, options,
!  Hessian approximation, QP subproblem solvers, line search, merit
!  function, filter and funnel, trust region, feasibility restoration,
!  second-order correction, convergence checking, the optional sparse
!  factorizations (inertia control, direct QP steps, and direct
!  least-squares solves), the diagnostics, and the detailed log)
!  are each implemented in their own module so that they may be
!  developed, tested, and swapped out independently. Internally, sparse
!  (COO) storage is used for the constraint Jacobian, and the Hessian of
!  the Lagrangian is a matrix-free limited-memory operator (or the user's
!  sparse exact Hessian), so dense `n x n`/`m x n` arrays are only formed
!  by the dense QP solver (used for small problems, see
!  [[sqpopt_qp_solver_module]]).

    module sqpopt_module

    use, intrinsic :: iso_fortran_env, only: int64
    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_types_module,      only: sqpopt_success, sqpopt_max_iter_reached, &
                                         sqpopt_user_requested_stop, sqpopt_report_func, &
                                         sqpopt_invalid_input, sqpopt_status_message, sqpopt_sparse_matrix, &
                                         sqpopt_results_type, sqpopt_max_evals_reached, sqpopt_time_limit_reached, &
                                         sqpopt_qp_solve_failed, sqpopt_infinity, sqpopt_all_finite, &
                                         sqpopt_iter_info, sqpopt_iteration_flags
    use sqpopt_problem_module,    only: sqpopt_problem_type, sqpopt_derivatives_fast, sqpopt_derivatives_accurate
    use sqpopt_options_module,    only: sqpopt_options_type
    use sqpopt_hessian_module,    only: sqpopt_hessian_type, sqpopt_hessian_bfgs, sqpopt_hessian_sr1, sqpopt_hessian_exact
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type, sqpopt_qp_auto, sqpopt_qp_dense, &
                                         sqpopt_qp_reduced_hessian, sqpopt_qp_daqp
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_null_space_lu, sqpopt_null_space_lsqr
    use sqpopt_eigen_module, only: sqpopt_eigen_auto, sqpopt_eigen_jacobi, sqpopt_eigen_lapack, sqpopt_eigen_ql, &
                                   sqpopt_eigen_has_lapack
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_armijo, sqpopt_linesearch_funnel, &
                                         sqpopt_linesearch_filter, sqpopt_linesearch_watchdog, &
                                         sqpopt_merit_l1, sqpopt_merit_augmented_lagrangian, &
                                         sqpopt_penalty_multipliers, sqpopt_penalty_model
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_restoration_module,  only: sqpopt_restoration_type, sqpopt_restoration_phase, &
                                          sqpopt_restoration_gauss_newton
    use sqpopt_iterate_module,    only: sqpopt_iterate, sqpopt_evaluate_point
    use sqpopt_log_module,        only: sqpopt_log_type, sqpopt_log_detail, fmt_e, fmt_i, plural
    use sqpopt_inertia_module,    only: sqpopt_inertia_type
    use sqpopt_kkt_module,        only: sqpopt_kkt_type
    use sqpopt_least_squares_module,    only: sqpopt_least_squares_type
    use sqpopt_symmetric_solver_module, only: sqpopt_linear_solver_available, sqpopt_linear_solver_name, &
                                              sqpopt_linear_solver_max_order, &
                                              sqpopt_linear_solver_mumps, sqpopt_linear_solver_lapack
    use sqpopt_diagnostics_module,      only: sqpopt_diagnostics_type, sqpopt_diagnostics_write

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
    type(sqpopt_problem_type),optional,intent(in)      :: problem      !! the nonlinear program to be solved
    type(sqpopt_options_type),optional,intent(in)      :: options      !! solver options
    type(sqpopt_hessian_type),optional,intent(in)      :: hessian      !! Hessian of the Lagrangian approximation
    type(sqpopt_qp_solver_type),optional,intent(in)    :: qp_solver    !! QP subproblem solver
    type(sqpopt_linesearch_type),optional,intent(in)   :: linesearch   !! merit function / line search
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
    integer :: iter_istat, iter, n_fail, n_acceptable, n_stalled, n_escape, n_fc0
    logical :: in_phase !! whether the previous iteration was in a restoration phase (to count the phases)
    integer :: ios !! (for the I/O statements of the printed output, which must never stop the solver)
    integer :: detail_unit !! scratch file for the detail lines of the current iteration (`print_level >= 3`;
                           !! `-1` if none)
    integer(int64) :: t_start, t_now, t_rate
    character(len=:), allocatable :: msg
    type(sqpopt_restoration_type) :: fresh_restoration !! (default-initialized)
    ! the optional sparse factorizations (see `options%linear_solver`). They live for one solve, and are freed by
    ! `finish`:
    type(sqpopt_kkt_type)     :: kkt     !! the KKT matrix of the QP's working set (see `options%inertia_control`
                                         !! and `options%direct_qp`)
    type(sqpopt_inertia_type) :: inertia !! inertia control of the exact Hessian (see `options%inertia_control`)
    type(sqpopt_least_squares_type) :: least_squares !! direct least-squares solver (see
                                                     !! `options%direct_least_squares`)
    type(sqpopt_diagnostics_type) :: diagnostics !! the diagnostics of this solve (see `options%diagnostic_level`)
    logical :: started

    call system_clock(t_start, t_rate)
    valid = .false.
    detail_unit = -1

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
    me%results = sqpopt_results_type()

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
    call me%problem%set_derivative_accuracy(me%options%derivative_accuracy)
    if (me%options%scaling) call me%problem%compute_scaling(me%x, me%options%scaling_max_gradient, &
                                                            me%options%scaling_min_value)
    if (present(lambda0)) me%lambda = lambda0*me%problem%f_scale/me%problem%c_scale

    ! the diagnostics (from level 2: the report on the starting point, from the
    ! derivatives that the scaling, or else the first iteration, evaluates there):
    call diagnostics%start(me%problem, me%options, me%x)
    me%qp_solver%dense_qp%keep_slacks  = me%options%diagnostic_level >= 1
    me%qp_solver%sparse_qp%keep_slacks = me%options%diagnostic_level >= 1

    call me%hessian%initialize(me%problem%n, &
                                lbfgs_memory(me%options%lbfgs_memory, me%problem%n, me%options%direct_qp), &
                                use_sr1=(me%options%hessian_mode == sqpopt_hessian_sr1), &
                                scale0=me%options%hessian_scale0)
    if (me%options%hessian_mode == sqpopt_hessian_exact) then
        call me%hessian%set_exact(me%problem%hess_irow, me%problem%hess_icol)
    end if
    ! the optional sparse factorizations (if the sparse solver can't be
    ! started, the solve continues without them):
    if (me%options%hessian_mode == sqpopt_hessian_exact) then
        if (me%options%inertia_control .or. me%options%direct_qp) then
            call kkt%initialize(me%problem%n, me%problem%m, me%problem%jac_irow, me%problem%jac_icol, started, &
                                hess_irow=me%problem%hess_irow, hess_icol=me%problem%hess_icol, &
                                threads=me%options%factorization_threads, solver=me%options%linear_solver)
        end if
        inertia%enabled = me%options%inertia_control .and. kkt%enabled
    else if (me%options%direct_qp .or. (me%options%inertia_control .and. me%options%hessian_mode == sqpopt_hessian_sr1)) then
        ! (a quasi-Newton Hessian has no sparsity pattern: see [[sqpopt_kkt_module]])
        call kkt%initialize(me%problem%n, me%problem%m, me%problem%jac_irow, me%problem%jac_icol, started, &
                            threads=me%options%factorization_threads, solver=me%options%linear_solver)
        ! (the BFGS matrix is positive definite: only SR1 needs the inertia control)
        inertia%enabled = me%options%inertia_control .and. me%options%hessian_mode == sqpopt_hessian_sr1 .and. kkt%enabled
    end if
    if (me%options%direct_least_squares .and. me%problem%m > 0) then
        call least_squares%initialize(me%problem%n, me%problem%m, me%problem%jac_irow, me%problem%jac_icol, started, &
                                      threads=me%options%factorization_threads, solver=me%options%linear_solver)
    end if
    me%qp_solver%mode        = me%options%qp_solver_mode
    me%qp_solver%direct      = me%options%direct_qp .and. kkt%enabled
    me%linesearch%mode       = me%options%linesearch_mode
    me%linesearch%merit%mode           = me%options%merit_mode
    me%linesearch%merit%penalty_update = me%options%penalty_update
    ! (the components that write detail lines to the log)
    me%linesearch%log   = sqpopt_log_type(unit=me%options%output_unit, level=me%options%print_level)
    if (me%options%print_level >= sqpopt_log_detail) then
        ! (the detail lines of an iteration go to a scratch file, and are copied
        ! out after the iteration's line of the log, see `print_details`)
        open(newunit=detail_unit, status='scratch', action='readwrite', form='formatted', iostat=ios)
        if (ios == 0) then
            me%linesearch%log%unit = detail_unit
        else
            detail_unit = -1   ! (then the detail lines go straight to `output_unit`, before their iteration's line)
        end if
    end if
    me%trust_region%log = me%linesearch%log

    if (me%options%print_level >= 1) call print_header()
    in_phase = .false.

    n_fail = 0
    n_acceptable = 0
    n_stalled    = 0
    n_escape     = 0
    do iter = 1, me%options%max_iter
        me%results%iterations = iter
        n_fc0 = me%problem%n_eval_fc
        call sqpopt_iterate(me%problem, me%options, me%hessian, me%qp_solver, me%linesearch, me%trust_region, &
                             me%x, me%lambda, x_prev, gl_prev, f_prev, viol_prev, jac, n_acceptable, n_stalled, n_escape, &
                             me%restoration, inertia, kkt, least_squares, diagnostics, iter, me%report, &
                             done, iter_istat, info)
        info%n_fc = me%problem%n_eval_fc - n_fc0
        call count_events()
        call diagnostics%record(iter, info, iter_istat, me%problem, me%qp_solver, me%options, &
                                kkt%solver%time + least_squares%kkt%solver%time)
        if (me%options%print_level >= 1) call print_iteration(iter, info, iter_istat)
        if (me%options%print_level >= sqpopt_log_detail) call print_details()
        if (done) then
            ! converged, stalled, acceptable, infeasible, unbounded, function error, out of memory, or user stop:
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
        integer, intent(in)                    :: stat   !! the final status code
        character(len=*), intent(in), optional :: detail !! extra detail appended to the status message
        real(wp) :: fs
        real(wp), dimension(:), allocatable :: cs
        real(wp), dimension(:), allocatable :: zs
        integer :: ios

        allocate(cs(me%problem%m), zs(size(me%x)))

        istat = stat
        me%results%istat   = stat
        me%results%message = sqpopt_status_message(stat)
        if (present(detail)) then
            if (len(detail) > 0) me%results%message = me%results%message//': '//detail
        end if

        me%results%x = me%x
        if (valid) then
            call sqpopt_evaluate_point(me%problem, me%options, me%x, me%lambda, jac, fs, cs, &
                                       me%results%kkt_error, me%results%feasibility_error, zs, &
                                       stat_error=me%results%stationarity_error)
            me%results%stationarity_error = me%results%stationarity_error/me%problem%f_scale
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
            if (allocated(me%results%c)) deallocate(me%results%c)
            allocate(me%results%c(max(me%problem%m,0)), source=0.0_wp)
            me%results%lambda = me%lambda
            if (allocated(me%results%z)) deallocate(me%results%z)
            allocate(me%results%z(size(x0)), source=0.0_wp)
            me%results%n_eval_fc = 0; me%results%n_eval_gjac = 0; me%results%n_eval_hess = 0
        end if
        call system_clock(t_now)
        me%results%time = real(t_now-t_start, wp)/real(t_rate, wp)
        me%results%time_functions = me%problem%time_user
        me%results%time_qp        = me%qp_solver%time
        me%results%n_qp_solves    = me%qp_solver%n_solves
        me%results%n_direct_qp    = me%qp_solver%n_direct
        me%results%n_unconstrained_qp = me%qp_solver%n_unconstrained
        me%results%n_daqp_fallbacks   = me%qp_solver%n_daqp_fallbacks
        me%results%n_factorizations   = kkt%solver%n_factor + least_squares%kkt%solver%n_factor
        me%results%time_factorization = kkt%solver%time + least_squares%kkt%solver%time
        if (valid) call diagnostics%finish(me%problem, me%options, me%results, me%x, me%lambda, cs, jac, &
                                           me%qp_solver, me%hessian%shift, kkt%enabled .and. kkt%singular)
        call kkt%destroy()
        call least_squares%destroy()

        ! (not for invalid inputs, which may include `output_unit` itself)
        if (valid .and. me%options%print_level >= 1) call print_summary()
        if (detail_unit /= -1) then
            close(detail_unit, iostat=ios)
            detail_unit = -1
        end if
        end subroutine finish

        subroutine print_details()
        !! copy the current iteration's detail lines from the scratch file to the log, and empty it
        character(len=1024) :: line
        integer :: ios
        if (detail_unit == -1) return
        rewind(detail_unit, iostat=ios)
        do
            read(detail_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            write(me%options%output_unit, '(A)', iostat=ios) trim(line)
        end do
        rewind(detail_unit, iostat=ios)
        endfile(detail_unit, iostat=ios)
        rewind(detail_unit, iostat=ios)
        end subroutine print_details

        subroutine count_events()
        !! add this iteration's events to the counts in the results
        me%results%n_qp_iterations = me%results%n_qp_iterations + info%qp_iter
        if (info%soc)         me%results%n_soc            = me%results%n_soc + 1
        if (info%hess_reset)  me%results%n_hessian_resets = me%results%n_hessian_resets + 1
        if (info%elastic)     me%results%n_elastic        = me%results%n_elastic + 1
        if (info%escape)      me%results%n_escape         = me%results%n_escape + 1
        if (info%derivatives) me%results%derivative_switch_iteration = iter
        if (info%restoration) me%results%n_restoration_steps = me%results%n_restoration_steps + 1
        if (info%phase .and. .not. in_phase) me%results%n_restoration_phases = me%results%n_restoration_phases + 1
        in_phase = info%phase .and. me%restoration%active
        end subroutine count_events

        subroutine print_header()
        !! the problem, the method, the scaling, the tolerances, and the
        !! iteration log's column headings (`print_level >= 1`; the legend is
        !! printed with the summary, see `print_legend`), at
        !! `print_level >= 3` the smallest constraint scale factors, and with
        !! `diagnostic_level >= 2` the diagnostics' report on the starting point
        integer :: u, n_eq
        integer :: ios
        u = me%options%output_unit
        n_eq = 0
        if (me%problem%m > 0) n_eq = count(me%problem%c_ub - me%problem%c_lb <= 0.0_wp)
        write(u, '(A)', iostat=ios) ''
        write(u, '(A)', iostat=ios) ' sqpopt: '//plural(me%problem%n, 'variable', 'variables')//', '// &
                       plural(me%problem%m, 'constraint', 'constraints')//' ('//plural(n_eq, 'equality', 'equalities')// &
                       '), '//plural(me%problem%jac_nnz, 'Jacobian nonzero', 'Jacobian nonzeros')
        write(u, '(A)', iostat=ios) '   method:     '//method_text()
        if (me%options%scaling) then
            write(u, '(A)', iostat=ios) '   scaling:    objective x '//fmt_e(me%problem%f_scale)//constraint_scale_text()
        else
            write(u, '(A)', iostat=ios) '   scaling:    off'
        end if
        if (allocated(me%problem%x_max_step)) then
            write(u, '(A)', iostat=ios) '   max step:   '// &
                fmt_i(count(me%problem%x_max_step < sqpopt_infinity))//' of '// &
                plural(me%problem%n, 'variable', 'variables')//' limited (set_max_step), smallest limit '// &
                fmt_e(minval(me%problem%x_max_step))
        end if
        write(u, '(A)', iostat=ios) '   tolerances: ktol '//fmt_e(me%options%ktol)//', ctol '//fmt_e(me%options%ctol)// &
                       ', dual_inf_tol '//fmt_e(me%options%dual_inf_tol)//', max_iter '//fmt_i(me%options%max_iter)
        if (me%options%print_level >= sqpopt_log_detail .and. me%problem%m > 0) call print_scale_factors()
        call sqpopt_diagnostics_write(u, diagnostics%problem_report())
        write(u, '(A)', iostat=ios) ''
        if (me%options%print_level >= 2) then
            write(u, '(A6,A17,3A10,A8,A10,2A6,2A10,2A10,2X,A)', iostat=ios) 'iter', 'objective', 'infeas*', 'kkt*', 'alpha', &
                'fc', '|step|', 'qp_it', 'ls_fc', '|lambda|', 'stat', glob_heading(), hess_heading(), 'flags'
        else
            write(u, '(A6,A17,3A10,A8,2X,A)', iostat=ios) 'iter', 'objective', 'infeas*', 'kkt*', 'alpha', 'fc', 'flags'
        end if
        end subroutine print_header

        subroutine print_legend()
        !! what the columns and flags mean (printed with the summary, so the log stays compact)
        integer :: u
        integer :: ios
        u = me%options%output_unit
        if (me%options%print_level >= 2) then
            write(u, '(A)', iostat=ios) '   columns: * = of the scaled problem (the convergence test); objective, |lambda| (the'
            write(u, '(A)', iostat=ios) '            largest multiplier) and stat (the stationarity error) are of the original'
            write(u, '(A)', iostat=ios) '            problem; fc = calls of fc so far, ls_fc = in this iteration; qp_it = QP'
            write(u, '(A)', iostat=ios) '            iterations; '//trim(adjustl(glob_heading()))//' = '//glob_meaning()//'; '// &
                           trim(adjustl(hess_heading()))//' = '//hess_meaning()
        else
            write(u, '(A)', iostat=ios) '   columns: * = of the scaled problem (the convergence test); the objective is of the'
            write(u, '(A)', iostat=ios) '            original problem; fc = calls of fc so far'
        end if
        write(u, '(A)', iostat=ios) '   flags:   R restoration step, P restoration phase, S second-order correction, H Hessian'
        write(u, '(A)', iostat=ios) '            reset, E elastic QP re-solve, X escape step, N non-monotone step, W watchdog'
        write(u, '(A)', iostat=ios) '            relaxed step, D switched to accurate derivatives, Q QP failed, F no'
        write(u, '(A)', iostat=ios) '            acceptable step'
        end subroutine print_legend

        subroutine print_scale_factors()
        !! the smallest constraint scale factors (at most 10, `print_level >= 3`)
        integer :: u, i, k, n_show
        integer, dimension(:), allocatable :: order
        integer :: ios
        allocate(order(me%problem%m))
        u = me%options%output_unit
        if (all(me%problem%c_scale == 1.0_wp)) return
        ! (the smallest factors first: the constraints scaled down the most)
        do i = 1, me%problem%m
            order(i) = i
        end do
        call sort_by_scale(order)
        n_show = min(me%problem%m, 10)
        write(u, '(A)', iostat=ios) '   constraint scale factors (smallest '//fmt_i(n_show)//' of '//fmt_i(me%problem%m)//'):'
        do k = 1, n_show
            i = order(k)
            write(u, '(A)', iostat=ios) '     c('//fmt_i(i)//') x '//fmt_e(me%problem%c_scale(i))
        end do
        end subroutine print_scale_factors

        subroutine sort_by_scale(order)
        !! sort constraint indices by increasing scale factor (insertion sort: `m` is only printed at level 3)
        integer, dimension(:), intent(inout) :: order !! the constraint indices, sorted in place
        integer :: i, j, t
        do i = 2, size(order)
            t = order(i)
            j = i - 1
            do while (j >= 1)
                if (me%problem%c_scale(order(j)) <= me%problem%c_scale(t)) exit
                order(j+1) = order(j)
                j = j - 1
            end do
            order(j+1) = t
        end do
        end subroutine sort_by_scale

        function method_text() result(str)
        !! the algorithms used, e.g. "filter line search, dense QP (auto), L-BFGS Hessian (10 pairs)"
        character(len=:), allocatable :: str
        character(len=:), allocatable :: glob, hess, merit
        merit = ''
        select case (me%options%merit_mode)
        case (sqpopt_merit_augmented_lagrangian); merit = 'augmented Lagrangian merit'
        case default;                              merit = 'l1 merit'
        end select
        if (me%options%penalty_update == sqpopt_penalty_model) then
            merit = merit//', model penalty'
        else
            merit = merit//', multiplier penalty'
        end if
        select case (me%options%linesearch_mode)
        case (sqpopt_linesearch_filter);   glob = 'filter'
        case (sqpopt_linesearch_funnel);   glob = 'funnel'
        case (sqpopt_linesearch_armijo);   glob = 'Armijo line search ('//merit//')'
        case (sqpopt_linesearch_watchdog); glob = 'watchdog line search ('//merit//')'
        case default;                      glob = 'exact line search ('//merit//')'
        end select
        if (me%trust_region%enabled) then
            select case (me%options%linesearch_mode)
            case (sqpopt_linesearch_filter, sqpopt_linesearch_funnel); glob = 'trust region with the '//glob
            case default; glob = 'trust region (ratio test, '//merit//')'
            end select
        else if (me%options%linesearch_mode == sqpopt_linesearch_filter .or. &
                 me%options%linesearch_mode == sqpopt_linesearch_funnel) then
            glob = glob//' line search'
        end if
        select case (me%options%hessian_mode)
        case (sqpopt_hessian_exact)
            hess = 'exact Hessian'
            if (inertia%enabled .and. me%qp_solver%direct) then
                hess = hess//' (inertia control, direct QP)'
            else if (inertia%enabled) then
                hess = hess//' (inertia control)'
            else if (me%qp_solver%direct) then
                hess = hess//' (direct QP)'
            end if
        case (sqpopt_hessian_sr1)
            hess = 'L-SR1 Hessian ('//fmt_i(me%hessian%max_history)//' pairs'
            if (me%hessian%convexify) hess = hess//', convexified'
            if (inertia%enabled) hess = hess//', inertia control'
            if (me%qp_solver%direct) hess = hess//', direct QP'
            hess = hess//')'
        case default
            hess = 'L-BFGS Hessian ('//fmt_i(me%hessian%max_history)//' pairs'
            if (me%qp_solver%direct) hess = hess//', direct QP'
            hess = hess//')'
        end select
        str = glob//', '//me%qp_solver%mode_name(me%problem%n)//', '//hess
        if (me%problem%m > 0) then
            if (me%options%restoration_mode == sqpopt_restoration_phase) then
                str = str//', restoration phases'
            else
                str = str//', Gauss-Newton restoration steps'
            end if
            if (least_squares%enabled) str = str//', direct least squares'
        end if
        if (kkt%enabled .or. least_squares%enabled) then
            str = str//', '//sqpopt_linear_solver_name(me%options%linear_solver)
            if (me%options%linear_solver /= sqpopt_linear_solver_mumps) then
                continue   ! (QDLDL is single-threaded)
            else if (me%options%factorization_threads == 0) then
                str = str//', factorizations on the OpenMP threads'
            else if (me%options%factorization_threads > 1) then
                str = str//', factorizations on '//fmt_i(me%options%factorization_threads)//' threads'
            end if
        end if
        end function method_text

        function constraint_scale_text() result(str)
        !! the range of the constraint scale factors
        character(len=:), allocatable :: str
        if (me%problem%m == 0) then
            str = ''
        else if (all(me%problem%c_scale == 1.0_wp)) then
            str = ', constraints unscaled'
        else
            str = ', constraints x ['//fmt_e(minval(me%problem%c_scale))//', '//fmt_e(maxval(me%problem%c_scale))// &
                  '] ('//fmt_i(count(me%problem%c_scale /= 1.0_wp))//' of '//fmt_i(me%problem%m)//' scaled)'
        end if
        end function constraint_scale_text

        pure function glob_heading() result(str)
        !! the heading of the globalization's column (`print_level >= 2`)
        character(len=10) :: str
        if (me%trust_region%enabled) then
            str = 'radius'
        else if (me%options%linesearch_mode == sqpopt_linesearch_filter) then
            str = 'filter'
        else if (me%options%linesearch_mode == sqpopt_linesearch_funnel) then
            str = 'funnel'
        else
            str = 'penalty'
        end if
        str = adjustr(str)
        end function glob_heading

        function glob_meaning() result(str)
        !! what the globalization's column of the log (`print_level >= 2`) shows, for the legend
        character(len=:), allocatable :: str
        if (me%trust_region%enabled) then
            str = 'trust-region radius'
        else if (me%options%linesearch_mode == sqpopt_linesearch_filter) then
            str = 'filter entries'
        else if (me%options%linesearch_mode == sqpopt_linesearch_funnel) then
            str = 'funnel width'
        else
            str = 'merit penalty parameter'
        end if
        end function glob_meaning

        pure function hess_heading() result(str)
        !! the heading of the Hessian's column (`print_level >= 2`)
        character(len=10) :: str
        str = merge('shift', 'pairs', me%options%hessian_mode == sqpopt_hessian_exact)
        str = adjustr(str)
        end function hess_heading

        function hess_meaning() result(str)
        !! what the Hessian's column of the log shows, for the legend
        character(len=:), allocatable :: str
        if (me%options%hessian_mode == sqpopt_hessian_exact) then
            str = 'Hessian shift (inertia correction)'
        else
            str = 'stored quasi-Newton pairs'
        end if
        end function hess_meaning

        subroutine print_iteration(iter, info, iter_istat)
        !! one line of the iteration log (see `print_legend` for the flags)
        integer,                intent(in) :: iter       !! major iteration number
        integer,                intent(in) :: iter_istat !! the iteration's status code
        type(sqpopt_iter_info), intent(in) :: info       !! what happened in the iteration
        character(len=12) :: flags
        character(len=10) :: gcol, hcol
        integer :: u
        integer :: ios
        u = me%options%output_unit
        flags = sqpopt_iteration_flags(info, iter_istat)
        if (.not. info%stepped) then
            ! (the final point: no step was taken from it)
            if (me%options%print_level >= 2) then
                write(u, '(I6,ES17.9,2ES10.2,A10,I8,A10,2A6,A10,ES10.2)', iostat=ios) iter, info%f/me%problem%f_scale, &
                    info%feas, info%kkt, '', me%problem%n_eval_fc, '', '', '', '', info%stat_unscaled
            else
                write(u, '(I6,ES17.9,2ES10.2,A10,I8)', iostat=ios) iter, info%f/me%problem%f_scale, info%feas, info%kkt, '', &
                    me%problem%n_eval_fc
            end if
        else if (me%options%print_level >= 2) then
            if (.not. me%trust_region%enabled .and. me%options%linesearch_mode == sqpopt_linesearch_filter) then
                write(gcol,'(I10)', iostat=ios) nint(info%glob)
            else
                write(gcol,'(ES10.2)', iostat=ios) info%glob
            end if
            if (ios /= 0) gcol = '      ****'
            if (me%options%hessian_mode == sqpopt_hessian_exact) then
                write(hcol,'(ES10.2)', iostat=ios) info%hess_measure
            else
                write(hcol,'(I10)', iostat=ios) nint(info%hess_measure)
            end if
            if (ios /= 0) hcol = '      ****'
            write(u, '(I6,ES17.9,3ES10.2,I8,ES10.2,2I6,2ES10.2,2A10,2X,A)', iostat=ios) iter, info%f/me%problem%f_scale, &
                info%feas, info%kkt, info%alpha, me%problem%n_eval_fc, info%step_norm, info%qp_iter, info%n_fc, &
                info%lam_max, info%stat_unscaled, gcol, hcol, trim(flags)
        else
            write(u, '(I6,ES17.9,3ES10.2,I8,2X,A)', iostat=ios) iter, info%f/me%problem%f_scale, info%feas, info%kkt, &
                info%alpha, me%problem%n_eval_fc, trim(flags)
        end if
        end subroutine print_iteration

        subroutine print_summary()
        !! the final summary (and, at `print_level >= 3`, the solution; and with
        !! `diagnostic_level >= 1`, the diagnosis)
        integer :: u
        real(wp) :: t_other
        character(len=:), allocatable :: events
        integer :: ios
        u = me%options%output_unit
        write(u, '(A)', iostat=ios) ''
        call print_legend()
        write(u, '(A)', iostat=ios) ''
        write(u, '(A,I0,2A)', iostat=ios)   ' sqpopt: status ', me%results%istat, ': ', me%results%message
        write(u, '(A,ES18.10)', iostat=ios) '   objective           = ', me%results%f
        write(u, '(A)', iostat=ios)         '   feasibility error   = '//fmt_e(me%results%feasibility_error)// &
                               '  (original problem)'
        write(u, '(A)', iostat=ios)         '   KKT error           = '//fmt_e(me%results%kkt_error)// &
                               '  (scaled problem; ktol = '//fmt_e(me%options%ktol)//')'
        write(u, '(A)', iostat=ios)         '   stationarity error  = '//fmt_e(me%results%stationarity_error)// &
                               '  (original problem; dual_inf_tol = '//fmt_e(me%options%dual_inf_tol)//')'
        if (me%problem%m > 0) then
            write(u, '(A)', iostat=ios)     '   largest multiplier  = '//fmt_e(maxval(abs(me%results%lambda)))
        end if
        write(u, '(A)', iostat=ios)         '   active              = '//active_text()
        write(u, '(A,I0)', iostat=ios)      '   iterations          = ', me%results%iterations
        if (me%results%n_eval_hess > 0) then
            write(u, '(A,3(I0,A))', iostat=ios) '   evaluations         = ', me%results%n_eval_fc, &
                                   ' fc, ', me%results%n_eval_gjac, &
                                   ' gjac, ', me%results%n_eval_hess, ' hess'
        else
            write(u, '(A,2(I0,A))', iostat=ios) '   evaluations         = ', &
                  me%results%n_eval_fc, ' fc, ', me%results%n_eval_gjac, ' gjac'
        end if
        if (me%results%n_direct_qp > 0) then
            write(u, '(A,3(I0,A))', iostat=ios) '   QP iterations       = ', me%results%n_qp_iterations, ' (', &
                                   me%results%n_direct_qp, ' of ', me%results%n_qp_solves, ' QPs solved directly)'
        else if (me%results%n_unconstrained_qp > 0) then
            write(u, '(A,3(I0,A))', iostat=ios) '   QP iterations       = ', me%results%n_qp_iterations, ' (', &
                                   me%results%n_unconstrained_qp, ' of ', me%results%n_qp_solves, &
                                   ' QPs solved by the unconstrained step)'
        else
            write(u, '(A,I0)', iostat=ios)  '   QP iterations       = ', me%results%n_qp_iterations
        end if
        if (me%results%n_daqp_fallbacks > 0) then
            write(u, '(A,2(I0,A))', iostat=ios) '   DAQP fallbacks      = ', me%results%n_daqp_fallbacks, ' of ', &
                                   me%results%n_qp_solves, ' QPs solved by the dense QP solver instead'
        end if
        if (me%results%n_factorizations > 0) then
            write(u, '(A,I0)', iostat=ios)  '   factorizations      = ', me%results%n_factorizations
        end if
        events = ''
        call add_event(events, me%results%n_soc, 'second-order correction', 'second-order corrections')
        call add_event(events, me%results%n_restoration_phases, 'restoration phase', 'restoration phases')
        call add_event(events, me%results%n_restoration_steps, 'restoration step', 'restoration steps')
        call add_event(events, me%results%n_hessian_resets, 'Hessian reset', 'Hessian resets')
        call add_event(events, me%results%n_elastic, 'elastic re-solve', 'elastic re-solves')
        call add_event(events, me%results%n_escape, 'escape step', 'escape steps')
        if (me%results%derivative_switch_iteration > 0) then
            if (len(events) > 0) events = events//', '
            events = events//'accurate derivatives from iteration '//fmt_i(me%results%derivative_switch_iteration)
        end if
        if (len(events) == 0) events = 'none'
        write(u, '(A)', iostat=ios)         '   events              = '//events
        t_other = max(0.0_wp, me%results%time - me%results%time_functions - me%results%time_qp &
                              - me%results%time_factorization)
        if (me%results%n_factorizations > 0) then
            write(u, '(A)', iostat=ios)     '   time                = '//fmt_f(me%results%time)//' s (user functions '// &
                                   fmt_f(me%results%time_functions)//' s, QP '//fmt_f(me%results%time_qp)// &
                                   ' s, factorizations '//fmt_f(me%results%time_factorization)// &
                                   ' s, other '//fmt_f(t_other)//' s)'
        else
            write(u, '(A)', iostat=ios)     '   time                = '//fmt_f(me%results%time)//' s (user functions '// &
                                   fmt_f(me%results%time_functions)//' s, QP '//fmt_f(me%results%time_qp)// &
                                   ' s, other '//fmt_f(t_other)//' s)'
        end if
        if (me%options%print_level >= sqpopt_log_detail) call print_solution()
        if (allocated(me%results%diagnosis%report)) then
            write(u, '(A)', iostat=ios) ''
            call sqpopt_diagnostics_write(u, me%results%diagnosis%report)
        end if
        write(u, '(A)', iostat=ios) ''
        end subroutine print_summary

        subroutine add_event(events, n, one, many)
        !! add a count to the summary's list of events (e.g. "3 Hessian resets")
        character(len=:), allocatable, intent(inout) :: events    !! the list so far
        integer,          intent(in)                 :: n         !! the count
        character(len=*), intent(in)                 :: one, many !! the event's name, singular and plural
        if (n == 0) return
        if (len(events) > 0) events = events//', '
        if (n == 1) then
            events = events//'1 '//one
        else
            events = events//fmt_i(n)//' '//many
        end if
        end subroutine add_event

        function fmt_f(t) result(str)
        !! a time in seconds
        real(wp), intent(in) :: t !! the time, in seconds
        character(len=:), allocatable :: str
        character(len=32) :: buf
        integer :: ios
        write(buf,'(F0.3)', iostat=ios) t
        if (ios /= 0) then
            str = '****'
            return
        end if
        str = trim(buf)
        if (str(1:1) == '.') str = '0'//str
        end function fmt_f

        function active_text() result(str)
        !! how many constraints and variable bounds are active at the solution
        character(len=:), allocatable :: str
        integer :: n_c, n_x
        n_c = count(constraint_side() /= 0)
        n_x = count(bound_side() /= 0)
        str = fmt_i(n_c)//' of '//fmt_i(me%problem%m)//' constraints, '//fmt_i(n_x)//' of '// &
              fmt_i(me%problem%n)//' variable bounds'
        end function active_text

        function constraint_side() result(side)
        !! for each constraint: 0 = inactive, -1 = at its lower bound, +1 = at its upper bound, 2 = equality
        integer, dimension(:), allocatable :: side
        real(wp) :: lb, ub, tol
        integer :: i
        allocate(side(me%problem%m))
        side = 0
        do i = 1, me%problem%m
            lb = me%problem%c_lb(i)/me%problem%c_scale(i)
            ub = me%problem%c_ub(i)/me%problem%c_scale(i)
            tol = max(me%options%ctol, 1.0e-8_wp)
            if (ub - lb <= 0.0_wp) then
                side(i) = 2
            else if (abs(me%results%c(i) - lb) <= tol*max(1.0_wp, abs(lb))) then
                side(i) = -1
            else if (abs(me%results%c(i) - ub) <= tol*max(1.0_wp, abs(ub))) then
                side(i) = 1
            end if
        end do
        end function constraint_side

        function bound_side() result(side)
        !! for each variable: 0 = free, -1 = at its lower bound, +1 = at its upper bound, 2 = fixed
        integer, dimension(:), allocatable :: side
        real(wp) :: tol
        integer :: j
        allocate(side(me%problem%n))
        side = 0
        tol = max(me%options%ctol, 1.0e-8_wp)
        do j = 1, me%problem%n
            if (me%problem%x_ub(j) - me%problem%x_lb(j) <= 0.0_wp) then
                side(j) = 2
            else if (abs(me%results%x(j) - me%problem%x_lb(j)) <= tol*max(1.0_wp, abs(me%problem%x_lb(j)))) then
                side(j) = -1
            else if (abs(me%results%x(j) - me%problem%x_ub(j)) <= tol*max(1.0_wp, abs(me%problem%x_ub(j)))) then
                side(j) = 1
            end if
        end do
        end function bound_side

        subroutine print_solution()
        !! the solution: the variables and constraints, with their bounds, multipliers, and which are active
        !! (`print_level >= 3`; at most `max_rows` of each)
        integer, parameter :: max_rows = 100
        integer :: u, i
        integer, dimension(:), allocatable :: cside
        integer, dimension(:), allocatable :: xside
        integer :: ios
        allocate(cside(me%problem%m), xside(me%problem%n))
        u = me%options%output_unit
        cside = constraint_side()
        xside = bound_side()
        write(u, '(A)', iostat=ios) ''
        write(u, '(A)', iostat=ios) '   variables:'
        write(u, '(A8,5A17)', iostat=ios) 'j', 'x', 'lower', 'upper', 'z', 'active'
        do i = 1, min(me%problem%n, max_rows)
            write(u, '(I8,4A17,A17)', iostat=ios) i, num(me%results%x(i)), num(me%problem%x_lb(i)), num(me%problem%x_ub(i)), &
                num(me%results%z(i)), side_text(xside(i), 'fixed   ')
        end do
        if (me%problem%n > max_rows) write(u, '(A)', iostat=ios) '     ... ('//fmt_i(me%problem%n - max_rows)//' more)'
        if (me%problem%m > 0) then
            write(u, '(A)', iostat=ios) ''
            write(u, '(A)', iostat=ios) '   constraints:'
            write(u, '(A8,5A17)', iostat=ios) 'i', 'c', 'lower', 'upper', 'lambda', 'active'
            do i = 1, min(me%problem%m, max_rows)
                write(u, '(I8,4A17,A17)', iostat=ios) i, num(me%results%c(i)), num(me%problem%c_lb(i)/me%problem%c_scale(i)), &
                    num(me%problem%c_ub(i)/me%problem%c_scale(i)), num(me%results%lambda(i)), &
                    side_text(cside(i), 'equality')
            end do
            if (me%problem%m > max_rows) write(u, '(A)', iostat=ios) '     ... ('//fmt_i(me%problem%m - max_rows)//' more)'
        end if
        end subroutine print_solution

        pure function side_text(side, both) result(str)
        !! the "active" column of the solution tables
        integer,          intent(in) :: side !! `-1` at the lower bound, `+1` at the upper bound, `2` fixed (or an equality), `0` neither
        character(len=*), intent(in) :: both !! the text for `side==2`
        character(len=17) :: str
        select case (side)
        case (-1); str = 'lower'
        case (1);  str = 'upper'
        case (2);  str = both
        case default; str = ''
        end select
        str = adjustr(str)
        end function side_text

        function num(v) result(str)
        !! a value for the solution tables (infinite bounds as `-inf`/`inf`)
        real(wp), intent(in) :: v !! the value
        character(len=17) :: str
        integer :: ios
        if (v <= -sqpopt_infinity) then
            str = '-inf'
            str = adjustr(str)
        else if (v >= sqpopt_infinity) then
            str = 'inf'
            str = adjustr(str)
        else
            write(str,'(ES17.8)', iostat=ios) v
            if (ios /= 0) str = '             ****'
        end if
        end function num

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
    if (.not. sqpopt_linear_solver_available(o%linear_solver)) then
        if (o%linear_solver == sqpopt_linear_solver_mumps) then
            msg = 'options%linear_solver = sqpopt_linear_solver_mumps requires a library built with MUMPS '// &
                  '(the HAS_MUMPS preprocessor directive)'
        else if (o%linear_solver == sqpopt_linear_solver_lapack) then
            msg = 'options%linear_solver = sqpopt_linear_solver_lapack requires a library built with LAPACK '// &
                  '(the HAS_LAPACK preprocessor directive)'
        else
            msg = 'options%linear_solver is not a valid sqpopt_linear_solver_* value'
        end if
        return
    end if
    ! (the factorizations' matrices have the order n+m: a solver with a size limit, the
    ! dense ones, is refused here rather than left unused)
    if (((o%inertia_control .and. o%hessian_mode /= sqpopt_hessian_bfgs) .or. o%direct_qp .or. &
         (o%direct_least_squares .and. me%problem%m > 0)) .and. &
        me%problem%n + me%problem%m > sqpopt_linear_solver_max_order(o%linear_solver)) then
        msg = 'options%linear_solver: the '//sqpopt_linear_solver_name(o%linear_solver)// &
              ' solver is limited to matrices of order '//fmt_i(sqpopt_linear_solver_max_order(o%linear_solver))// &
              ' (this problem''s are of order n+m = '//fmt_i(me%problem%n + me%problem%m)// &
              '): use sqpopt_linear_solver_qdldl'
        return
    end if
    if (o%factorization_threads < 0) then
        msg = 'options%factorization_threads must be >= 0 (0: as the OpenMP environment says)'
        return
    end if
    if (all(o%qp_solver_mode /= [sqpopt_qp_auto, sqpopt_qp_dense, sqpopt_qp_reduced_hessian, sqpopt_qp_daqp])) then
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
    if (all(o%derivative_accuracy /= [sqpopt_derivatives_fast, sqpopt_derivatives_accurate]) .or. &
        .not. (o%derivative_switch_tol >= 0.0_wp)) then
        msg = 'options%derivative_accuracy must be a sqpopt_derivatives_* value, and derivative_switch_tol >= 0'
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
    if (.not. o%acceptable_obj_change_tol >= 0.0_wp) then
        msg = 'options%acceptable_obj_change_tol must be >= 0'
        return
    end if
    if (.not. (o%scaling_min_value >= 0.0_wp .and. o%scaling_min_value <= 1.0_wp)) then
        msg = 'options%scaling_min_value must be in [0, 1]'
        return
    end if
    if (o%print_level > 0) then
        inquire(unit=o%output_unit, opened=opened)
        if (.not. opened) then
            msg = 'options%output_unit is not an open unit'
            return
        end if
    end if
    if (o%diagnostic_level < 0 .or. o%diagnostic_level > 3) then
        msg = 'options%diagnostic_level must be 0, 1, 2, or 3'
        return
    end if
    if (o%diagnostic_level >= 2 .and. o%diagnostics_unit /= -1) then
        inquire(unit=o%diagnostics_unit, opened=opened)
        if (.not. opened) then
            msg = 'options%diagnostics_unit is not an open unit'
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
    if (qp%direct_max_changes < 0 .or. .not. qp%direct_tol > 0.0_wp) then
        msg = 'qp_solver%direct_max_changes must be >= 0, and direct_tol > 0'
        return
    end if
    associate (d => qp%dense_qp, r => qp%sparse_qp, a => qp%daqp_qp)
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
    if (.not. (a%primal_tol > 0.0_wp .and. a%dual_tol > 0.0_wp) .or. a%max_iter < 1) then
        msg = 'a qp_solver%daqp_qp setting is out of range'
        return
    end if
    end associate

    ! ---- Hessian ----
    if (all(me%hessian%eigen_solver /= [sqpopt_eigen_auto, sqpopt_eigen_jacobi, sqpopt_eigen_lapack, &
                                             sqpopt_eigen_ql])) then
        msg = 'hessian%eigen_solver is not a valid sqpopt_eigen_* value'
        return
    end if
    if (me%hessian%eigen_solver == sqpopt_eigen_lapack .and. .not. sqpopt_eigen_has_lapack) then
        msg = 'hessian%eigen_solver = sqpopt_eigen_lapack needs a library built with LAPACK (HAS_LAPACK)'
        return
    end if

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
    type(sqpopt_results_type), intent(out) :: results !! the results of the last `solve`

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
!  \( \max(10, \min(n, 100)) \), or 10 with the direct QP method (`direct`),
!  whose cost grows with the square of the number of pairs.

    pure integer function lbfgs_memory(memory, n, direct)

    integer, intent(in) :: memory !! the option value
    integer, intent(in) :: n      !! number of variables
    logical, intent(in) :: direct !! whether the direct QP method is in use (`options%direct_qp`)

    integer, parameter :: direct_memory = 10 !! the automatic memory with the direct QP method

    if (memory > 0) then
        lbfgs_memory = memory
    else if (direct) then
        lbfgs_memory = direct_memory
    else
        lbfgs_memory = max(10, min(n, 100))
    end if

    end function lbfgs_memory
!*******************************************************************************

    end module sqpopt_module
!*******************************************************************************
