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

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_types_module,      only: sqpopt_success, sqpopt_max_iter_reached, &
                                         sqpopt_user_requested_stop, sqpopt_report_func, &
                                         sqpopt_invalid_input, sqpopt_status_message, sqpopt_sparse_matrix
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type, sqpopt_hessian_bfgs, sqpopt_hessian_sr1, sqpopt_hessian_exact
    use sqpopt_hessian_module,    only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type, sqpopt_qp_auto, sqpopt_qp_reduced_hessian
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_armijo, sqpopt_linesearch_filter, &
                                         sqpopt_merit_l1, sqpopt_merit_augmented_lagrangian
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_iterate_module,    only: sqpopt_iterate

    implicit none

    private

    type, public :: sqpopt_type
        !! main class for the SQP optimizer.

        private

        type(sqpopt_problem_type)    :: problem      !! the nonlinear program to be solved
        type(sqpopt_options_type)    :: options      !! solver options
        type(sqpopt_hessian_type)    :: hessian      !! Hessian of the Lagrangian approximation
        type(sqpopt_qp_solver_type)  :: qp_solver    !! QP subproblem solver
        type(sqpopt_linesearch_type) :: linesearch   !! merit function / line search
        type(sqpopt_trust_region_type) :: trust_region !! trust-region globalization (opt-in alternative to `linesearch`)
        type(sqpopt_linesearch_type)   :: linesearch0   !! `linesearch` as given to `initialize`: restored at the start of
                                                        !! every `solve`, so no state (penalty, filter, watchdog) carries over
        type(sqpopt_trust_region_type) :: trust_region0 !! `trust_region` as given to `initialize` (same reason)
        type(sqpopt_qp_solver_type)    :: qp_solver0    !! `qp_solver` as given to `initialize` (same reason: no QP
                                                        !! warm-start state carries over)

        real(wp), dimension(:), allocatable :: x       !! current/final optimization variables
        real(wp), dimension(:), allocatable :: lambda  !! current/final Lagrange multipliers

        procedure(sqpopt_report_func), pointer, nopass :: report => null() !! optional user progress-reporting
                                                                           !! callback, invoked once per major
                                                                           !! iteration (see [[sqpopt_types_module]])

        integer :: iter  = 0  !! number of major iterations performed
        integer :: istat = 0  !! solver status code (see [[sqpopt_types_module]])
        character(len=:), allocatable :: message !! description of the final status (see `status_message`)

        contains

        private

        procedure, public :: initialize   => sqpopt_initialize
        procedure, public :: solve        => sqpopt_solve
        procedure, public :: get_solution => sqpopt_get_solution
        procedure, public :: destroy      => sqpopt_destroy
        procedure, public :: status_message => sqpopt_get_status_message

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
    procedure(sqpopt_report_func), optional, pointer :: report      !! optional user progress-reporting callback,
                                                                     !! called once per major iteration with the
                                                                     !! current iterate; set its `user_stop` output
                                                                     !! to request the solver stop early
                                                                     !! (see [[sqpopt_types_module]])

    ! use inputs if present, else use defaults:
    if (present(problem))    then; me%problem = problem; else; me%problem = sqpopt_problem_type(); end if
    if (present(options))    then; me%options = options; else; me%options = sqpopt_options_type(); end if
    if (present(hessian))    then; me%hessian = hessian; else; me%hessian = sqpopt_hessian_type(); end if
    if (present(qp_solver))  then; me%qp_solver = qp_solver; else; me%qp_solver = sqpopt_qp_solver_type(); end if
    if (present(linesearch)) then; me%linesearch = linesearch; else; me%linesearch = sqpopt_linesearch_type(); end if
    if (present(trust_region)) then; me%trust_region = trust_region; else; me%trust_region = sqpopt_trust_region_type(); end if
    me%linesearch0   = me%linesearch
    me%trust_region0 = me%trust_region
    me%qp_solver0    = me%qp_solver
    me%report => null()
    if (present(report)) then
        if (associated(report)) me%report => report
    end if

    if (allocated(me%x))      deallocate(me%x)
    if (allocated(me%lambda)) deallocate(me%lambda)
    me%iter  = 0
    me%istat = 0
    me%message = ''

    end subroutine sqpopt_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  solve the nonlinear program starting from the initial guess `x0`,
!  running major SQP iterations until convergence or a stopping
!  criterion is met.

    subroutine sqpopt_solve(me, x0, istat)

    class(sqpopt_type),     intent(inout) :: me
    real(wp), dimension(:), intent(in)    :: x0     !! initial guess for the optimization variables `dimension(n)`
    integer,                intent(out)   :: istat  !! status code (see [[sqpopt_types_module]])

    real(wp), dimension(:), allocatable :: x_prev, gl_prev  !! quasi-Newton state (unallocated until the 2nd iteration)
    real(wp), allocatable :: f_prev  !! previous objective value, for the `options%ftol` stalled-progress test (unallocated until the 2nd iteration)
    type(sqpopt_sparse_matrix) :: jac !! Jacobian workspace (structure set once, values updated each iteration)
    logical :: done
    integer :: iter_istat, iter, n_fail
    character(len=:), allocatable :: msg

    me%x = x0
    if (allocated(me%lambda)) deallocate(me%lambda)
    allocate(me%lambda(max(me%problem%m,0)))
    me%lambda = 0.0_wp
    me%iter = 0

    ! check the inputs before doing anything else:
    call me%problem%validate(istat, msg)
    if (istat == sqpopt_success) call validate_options(me%options, istat, msg)
    if (istat == sqpopt_success .and. size(x0) /= me%problem%n) then
        istat = sqpopt_invalid_input
        msg   = 'x0 must have size n'
    end if
    if (istat /= sqpopt_success) then
        call finish(istat, msg)
        return
    end if

    ! start from a point that satisfies the variable bounds, so the user
    ! functions are never evaluated outside them:
    me%x = min(max(x0, me%problem%x_lb), me%problem%x_ub)

    ! start every solve from the components exactly as configured (no
    ! state from a previous solve carries over), with an empty evaluation cache:
    call me%problem%reset_cache()
    me%linesearch   = me%linesearch0
    me%trust_region = me%trust_region0
    me%qp_solver    = me%qp_solver0
    call me%hessian%initialize(me%problem%n, me%options%lbfgs_memory, &
                                use_sr1=(me%options%hessian_mode == sqpopt_hessian_sr1))
    me%qp_solver%mode               = me%options%qp_solver_mode
    me%linesearch%mode              = me%options%linesearch_mode
    me%linesearch%merit_mode        = me%options%merit_mode

    n_fail = 0
    do iter = 1, me%options%max_iter
        me%iter = iter
        call sqpopt_iterate(me%problem, me%options, me%hessian, me%qp_solver, me%linesearch, me%trust_region, &
                             me%x, me%lambda, x_prev, gl_prev, f_prev, jac, iter, me%report, done, iter_istat)
        if (done) then
            ! converged, stalled, infeasible, or user stop:
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
    end do

    call finish(sqpopt_max_iter_reached)

    contains

        subroutine finish(stat, detail)
        !! set the final status (and its message) on both `istat` and `me`
        integer, intent(in) :: stat
        character(len=*), intent(in), optional :: detail !! extra detail appended to the status message
        istat    = stat
        me%istat = stat
        me%message = sqpopt_status_message(stat)
        if (present(detail)) then
            if (len(detail) > 0) me%message = me%message//': '//detail
        end if
        end subroutine finish

    end subroutine sqpopt_solve
!*******************************************************************************

!*******************************************************************************
!>
!  check the options for invalid values.

    subroutine validate_options(options, istat, msg)

    type(sqpopt_options_type),     intent(in)  :: options
    integer,                       intent(out) :: istat !! `sqpopt_success` or `sqpopt_invalid_input`
    character(len=:), allocatable, intent(out) :: msg   !! description of the problem found (empty if none)

    istat = sqpopt_invalid_input
    msg   = ''

    if (options%max_iter < 0) then
        msg = 'options%max_iter must be >= 0'
    else if (options%lbfgs_memory < 1) then
        msg = 'options%lbfgs_memory must be >= 1'
    else if (options%max_consecutive_failures < 1) then
        msg = 'options%max_consecutive_failures must be >= 1'
    else if (options%hessian_mode < sqpopt_hessian_bfgs .or. options%hessian_mode > sqpopt_hessian_exact) then
        msg = 'options%hessian_mode is not a valid sqpopt_hessian_* value'
    else if (options%qp_solver_mode < sqpopt_qp_auto .or. options%qp_solver_mode > sqpopt_qp_reduced_hessian) then
        msg = 'options%qp_solver_mode is not a valid sqpopt_qp_* value'
    else if (options%linesearch_mode < sqpopt_linesearch_armijo .or. &
             options%linesearch_mode > sqpopt_linesearch_filter) then
        msg = 'options%linesearch_mode is not a valid sqpopt_linesearch_* value'
    else if (options%merit_mode < sqpopt_merit_l1 .or. options%merit_mode > sqpopt_merit_augmented_lagrangian) then
        msg = 'options%merit_mode is not a valid sqpopt_merit_* value'
    else if (.not. (options%ktol > 0.0_wp .and. options%ctol > 0.0_wp)) then
        msg = 'options%ktol and options%ctol must be > 0'
    else if (.not. (options%ftol >= 0.0_wp .and. options%xtol >= 0.0_wp)) then
        msg = 'options%ftol and options%xtol must be >= 0'
    else
        istat = sqpopt_success
    end if

    end subroutine validate_options
!*******************************************************************************

!*******************************************************************************
!>
!  a description of the status of the last `solve` (including, for
!  `sqpopt_invalid_input`, what was invalid).

    function sqpopt_get_status_message(me) result(msg)

    class(sqpopt_type), intent(in) :: me
    character(len=:), allocatable :: msg

    if (allocated(me%message)) then
        msg = me%message
    else
        msg = ''
    end if

    end function sqpopt_get_status_message
!*******************************************************************************

!*******************************************************************************
!>
!  return the current (or final) solution and associated Lagrange multipliers.

    subroutine sqpopt_get_solution(me, x, lambda)

    class(sqpopt_type),     intent(in)  :: me
    real(wp), dimension(:), intent(out) :: x       !! optimization variables `dimension(n)`
    real(wp), dimension(:), intent(out) :: lambda  !! Lagrange multipliers `dimension(m)`

    x      = me%x
    lambda = me%lambda

    end subroutine sqpopt_get_solution
!*******************************************************************************

!*******************************************************************************
!>
!  destroy the solver instance, deallocating all internal arrays.

    subroutine sqpopt_destroy(me)

    class(sqpopt_type), intent(inout) :: me

    call me%initialize()

    end subroutine sqpopt_destroy
!*******************************************************************************

    end module sqpopt_module
!*******************************************************************************