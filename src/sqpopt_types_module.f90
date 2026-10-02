!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Common types and status codes used throughout `sqpopt`. Sparse matrices
!  are always stored in coordinate (COO) format with 1-based indices --
!  the same convention used by the `lusol` and `LSQR` dependencies
!  -- so that they may be passed to those solvers without conversion.

    module sqpopt_types_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    private

    ! solver status/exit codes:
    integer, parameter, public :: sqpopt_success              = 0  !! converged successfully
    integer, parameter, public :: sqpopt_acceptable           = 1  !! converged to the "acceptable" (looser) tolerances, for
                                                                   !! several consecutive iterations, but not to the normal ones
    integer, parameter, public :: sqpopt_stalled              = 2  !! stopped -- the point is feasible but the objective and variables
                                                                   !! have stopped changing (see `ftol`/`xtol`) before the KKT test
                                                                   !! was satisfied -- usually an acceptable, if imprecise, solution

    integer, parameter, public :: sqpopt_max_iter_reached     = 10 !! stopped -- maximum number of iterations reached
    integer, parameter, public :: sqpopt_user_requested_stop  = 11 !! stopped -- user requested stop
    integer, parameter, public :: sqpopt_max_evals_reached    = 12 !! stopped -- maximum number of function evaluations reached
    integer, parameter, public :: sqpopt_time_limit_reached   = 13 !! stopped -- time limit reached

    integer, parameter, public :: sqpopt_infeasible           = 21  !! stopped -- problem appears to be infeasible
    integer, parameter, public :: sqpopt_line_search_failed   = 22  !! stopped -- line search failed to find an acceptable step
    integer, parameter, public :: sqpopt_qp_solve_failed      = 23  !! stopped -- QP subproblem solver failed
    integer, parameter, public :: sqpopt_invalid_input        = 24  !! stopped -- the problem definition or options are invalid
    integer, parameter, public :: sqpopt_function_error       = 25  !! stopped -- the problem functions returned a non-finite value
                                                                    !! (NaN or Inf) at the current point (trial points with
                                                                    !! non-finite values are rejected instead)
    integer, parameter, public :: sqpopt_unbounded            = 26 !! stopped -- the objective fell below its lower limit at a
                                                                   !! feasible point (the problem appears to be unbounded)
    integer, parameter, public :: sqpopt_out_of_memory        = 27 !! stopped -- an array could not be allocated (the dense
                                                                   !! matrices of the dense QP solver)

    real(wp), parameter, public :: sqpopt_infinity = 1.0e20_wp !! any bound with magnitude `>= sqpopt_infinity` is treated
                                                                !! as absent (bounds are clamped to `[-sqpopt_infinity,
                                                                !! sqpopt_infinity]` so that e.g. `-huge(1.0_wp)` is safe to use)

    public :: sqpopt_status_message
    public :: sqpopt_all_finite
    public :: sqpopt_iteration_flags
    public :: l1_violation
    public :: merit_slack

    type, public :: sqpopt_sparse_matrix
        !! a sparse matrix stored in coordinate (COO) format
        integer :: nrows = 0  !! number of rows
        integer :: ncols = 0  !! number of columns
        integer :: nnz   = 0  !! number of nonzero elements
        integer,  dimension(:), allocatable :: irow  !! row indices `dimension(nnz)`
        integer,  dimension(:), allocatable :: icol  !! column indices `dimension(nnz)`
        real(wp), dimension(:), allocatable :: val   !! nonzero values `dimension(nnz)`
    end type sqpopt_sparse_matrix

    type, public :: sqpopt_diagnosis_type
        !! the diagnostics of a solve (see `options%diagnostic_level`): which
        !! constraints and variables are responsible for the errors at the
        !! final point, what the active set looks like there, and, from level
        !! `2`, what the history of the iterations and the starting point
        !! say. Values are for the original (unscaled) problem. Everything
        !! here is also in the text of `report` and `problem_report`, which
        !! is what the solver prints.
        integer  :: level = 0 !! the diagnostic level these were computed at (`0`: none, and nothing below is set)

        ! ---- the final point (level >= 1) ----
        integer,  dimension(:), allocatable :: violated_constraints !! the constraints with the largest violation of
                                                                    !! their bounds, largest first (at most 5)
        real(wp), dimension(:), allocatable :: violations           !! their violations
        integer,  dimension(:), allocatable :: stationarity_variables !! the variables with the largest stationarity
                                                                      !! residual (the gradient of the Lagrangian,
                                                                      !! projected onto the bounds), largest first
                                                                      !! (at most 5)
        real(wp), dimension(:), allocatable :: stationarity_residuals !! their residuals
        integer,  dimension(:), allocatable :: multiplier_constraints !! the constraints with the largest multipliers,
                                                                      !! largest first (at most 5)
        real(wp), dimension(:), allocatable :: multipliers          !! their multipliers
        integer  :: n_active_constraints = 0 !! constraints at one of their bounds (including the equalities)
        integer  :: n_active_bounds      = 0 !! variables at one of their bounds (including the fixed ones)
        integer  :: n_weakly_active      = 0 !! active inequality constraints and bounds whose multiplier is
                                             !! (nearly) zero
        integer  :: n_wrong_sign         = 0 !! active inequality constraints and bounds whose multiplier has the
                                             !! wrong sign for the bound they are at
        integer  :: n_dependent          = -1 !! active constraints whose gradients are linearly dependent on those
                                              !! of the others, in the variables that are not at a bound (`-1`: not
                                              !! computed)
        integer,  dimension(:), allocatable :: elastic_constraints !! the constraints that were elastic in the last
                                                                   !! QP subproblem (its linearized constraints were
                                                                   !! inconsistent, or it was an elastic re-solve),
                                                                   !! largest slack first (at most 5)
        real(wp), dimension(:), allocatable :: elastic_slacks      !! their slacks: by how much the QP's step violates
                                                                   !! each linearized constraint

        ! ---- the history of the iterations (level >= 2) ----
        real(wp) :: convergence_rate = 0.0_wp !! the factor by which the error (the larger of the KKT error over
                                              !! `ktol` and the feasibility error over `ctol`) changed per
                                              !! iteration, fitted over the last iterations (at most 20; `< 1`:
                                              !! converging; `0`: not enough iterations to tell)
        integer  :: iterations_needed = -1    !! the further iterations it would take to reach the tolerances at
                                              !! that rate (`-1`: unknown, or not converging)
        integer  :: n_active_set_flips = 0    !! the consecutive iterations, up to the last one, in which the QP's
                                              !! working set went back to what it was two iterations before
        logical  :: objective_derivative_suspect = .false. !! whether the objective's changes along the steps
                                                           !! disagreed with its gradient (see `derivative_suspects`)
        integer,  dimension(:), allocatable :: derivative_suspects !! the constraints whose changes along the steps
                                                                   !! disagreed with their Jacobian rows in most of
                                                                   !! the steps (the worst first, at most 5): their
                                                                   !! derivatives may be wrong
        integer  :: slowest_iteration = 0     !! the iteration that took the most time
        real(wp) :: slowest_iteration_time = 0.0_wp !! its time (seconds)

        ! ---- the starting point (level >= 2) ----
        real(wp) :: constraint_gradient_ratio = 0.0_wp !! the largest over the smallest size of a constraint's
                                                       !! gradient at the starting point (`0`: no constraints)
        real(wp) :: variable_gradient_ratio   = 0.0_wp !! the largest over the smallest size of a variable's
                                                       !! column of the objective gradient and the Jacobian there
        integer  :: n_constant_constraints    = 0      !! constraints with no entry in the Jacobian's pattern
        integer  :: n_single_variable_constraints = 0  !! constraints with a single entry there (a bound on the
                                                       !! variable, if they are linear)
        integer  :: n_dependent_equalities    = -1     !! equality constraints whose gradients at the starting
                                                       !! point are linearly dependent on those of the others
                                                       !! (`-1`: not computed)
        logical  :: probe_not_finite = .false. !! level `3`: whether `f` or `c` is not finite at one of the two
                                               !! points probed next to the starting point

        character(len=:), allocatable :: problem_report !! the report on the starting point (level `>= 2`), as
                                                        !! printed before the iteration log: its lines, each
                                                        !! ended by a newline character
        character(len=:), allocatable :: report         !! the diagnosis of the solve, as printed after the
                                                        !! summary: its lines, each ended by a newline character
    end type sqpopt_diagnosis_type

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

    type, public :: sqpopt_results_type
        !! the outcome of a solve (see `sqpopt_type%get_results`). Values are
        !! for the original (unscaled) problem.
        integer  :: istat      = 0       !! status code
        character(len=:), allocatable :: message !! description of the status
        integer  :: iterations = 0       !! number of major iterations performed
        integer  :: n_eval_fc   = 0      !! number of calls of the objective and constraint function (`fc`)
        integer  :: n_eval_gjac = 0      !! number of calls of the gradient and Jacobian function (`gjac`)
        integer  :: n_eval_hess = 0      !! number of calls of the Lagrangian Hessian function (`hess`)
        real(wp) :: f          = 0.0_wp  !! objective function value at `x`
        real(wp) :: kkt_error  = 0.0_wp  !! KKT (stationarity/complementarity) error at `x` (of the scaled problem,
                                          !! as used by the convergence test)
        real(wp) :: feasibility_error = 0.0_wp !! largest violation of a constraint or variable bound at `x`
        real(wp) :: stationarity_error = 0.0_wp !! stationarity error at `x` of the original (unscaled) problem (see
                                                 !! `options%dual_inf_tol`)
        real(wp) :: time       = 0.0_wp  !! wall-clock time of the solve (seconds)
        real(wp) :: time_functions = 0.0_wp !! of which, in the user's functions
        real(wp) :: time_qp        = 0.0_wp !! of which, in the QP subproblem solver (without the factorizations)
        real(wp) :: time_factorization = 0.0_wp !! of which, in the sparse factorizations and their solves (see
                                                !! `n_factorizations`)
        integer  :: n_qp_iterations      = 0 !! total active-set iterations of the QP subproblem solves (and, for
                                             !! those solved directly, changes of the working set)
        integer  :: n_qp_solves          = 0 !! QP subproblems solved (without those of restoration phases)
        integer  :: n_direct_qp          = 0 !! of which, solved directly, without the active-set solver (see
                                             !! `options%direct_qp`)
        integer  :: n_unconstrained_qp   = 0 !! of which, solved by the unconstrained step of the limited-memory
                                             !! BFGS Hessian, without a QP solver (see
                                             !! `qp_solver%unconstrained_step`)
        integer  :: n_factorizations     = 0 !! sparse factorizations, by the inertia control, the direct QP
                                             !! method, and the direct least-squares solves (see
                                             !! `options%inertia_control`, `direct_qp`, and
                                             !! `direct_least_squares`; `0` without them)
        integer  :: n_soc                = 0 !! accepted second-order-corrected steps
        integer  :: n_hessian_resets     = 0 !! iterations in which the Hessian approximation was reset (or, with the
                                             !! exact Hessian, its shift increased)
        integer  :: n_restoration_steps  = 0 !! iterations that took a feasibility-restoration step (including those
                                             !! of restoration phases)
        integer  :: n_restoration_phases = 0 !! restoration phases entered
        integer  :: n_elastic            = 0 !! QP re-solves with constraints elastic (see
                                             !! `options%elastic_multiplier_limit`)
        integer  :: n_escape             = 0 !! escape steps from a stationary point of the violation
        integer  :: derivative_switch_iteration = 0 !! the iteration at which the solver switched from fast to
                                                    !! accurate derivatives (`0` if it didn't: see
                                                    !! `options%derivative_accuracy`)
        real(wp), dimension(:), allocatable :: x      !! final point `dimension(n)`
        real(wp), dimension(:), allocatable :: c      !! constraint values at `x` `dimension(m)`
        real(wp), dimension(:), allocatable :: lambda !! constraint multipliers `dimension(m)` (for the Lagrangian
                                                      !! \( f - \lambda^T c - z^T x \): \( \lambda_i \ge 0 \) at a lower
                                                      !! bound, \( \le 0 \) at an upper bound)
        real(wp), dimension(:), allocatable :: z      !! variable-bound multipliers `dimension(n)` (same sign convention;
                                                      !! zero for a variable not at a bound)
        type(sqpopt_diagnosis_type) :: diagnosis      !! the diagnostics of the solve (only with
                                                      !! `options%diagnostic_level >= 1`)
    end type sqpopt_results_type

    abstract interface
        subroutine sqpopt_report_func(iter, x, f, c, lambda, user_stop, data)
            !! user-supplied callback invoked once per major SQP iteration
            !! for progress monitoring (see `sqpopt_type%initialize`'s
            !! `report` argument). Set `user_stop=.true.` to have the
            !! solver stop after the current iteration
            !! (`istat=sqpopt_user_requested_stop`). The values are for the
            !! original (unscaled) problem.
            import :: wp
            implicit none
            integer,                intent(in)    :: iter      !! major iteration number (starts at 1)
            real(wp), dimension(:), intent(in)    :: x         !! current optimization variables `dimension(n)`
            real(wp),               intent(in)    :: f         !! current objective function value
            real(wp), dimension(:), intent(in)    :: c         !! current constraint values `dimension(m)`
            real(wp), dimension(:), intent(in)    :: lambda    !! current Lagrange multiplier estimate `dimension(m)`
            logical,                intent(out)   :: user_stop !! set `.true.` to request the solver stop
            class(*), optional,     intent(inout) :: data      !! the user data given to `set_functions`
        end subroutine sqpopt_report_func
    end interface

    public :: sqpopt_report_func

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  a human-readable description of a solver status code.

    pure function sqpopt_status_message(istat) result(msg)

    integer, intent(in) :: istat !! status code (one of the `sqpopt_*` status parameters)
    character(len=:), allocatable :: msg

    select case (istat)
    case (sqpopt_success);             msg = 'converged successfully'
    case (sqpopt_acceptable);          msg = 'converged to an acceptable level'
    case (sqpopt_stalled);             msg = 'feasible, but no further progress is being made'
    case (sqpopt_max_iter_reached);    msg = 'maximum number of iterations reached'
    case (sqpopt_user_requested_stop); msg = 'user requested stop'
    case (sqpopt_max_evals_reached);   msg = 'maximum number of function evaluations reached'
    case (sqpopt_time_limit_reached);  msg = 'time limit reached'
    case (sqpopt_infeasible);          msg = 'problem appears to be (locally) infeasible'
    case (sqpopt_line_search_failed);  msg = 'line search failed to find an acceptable step'
    case (sqpopt_qp_solve_failed);     msg = 'QP subproblem solver failed'
    case (sqpopt_invalid_input);       msg = 'invalid problem definition or options'
    case (sqpopt_function_error);      msg = 'the problem functions returned a non-finite value (NaN or Inf), '// &
                                               'or could not be evaluated, at the current point'
    case (sqpopt_unbounded);           msg = 'the objective fell below its lower limit (the problem appears to be unbounded)'
    case (sqpopt_out_of_memory);       msg = 'out of memory (an array could not be allocated)'
    case default;                      msg = 'unknown status code'
    end select

    end function sqpopt_status_message
!*******************************************************************************

!*******************************************************************************
!>
!  the flags of an iteration, as in the last column of the iteration log
!  (see `options%print_level`): `R` restoration step, `P` restoration
!  phase, `S` second-order correction, `H` Hessian reset, `E` elastic QP
!  re-solve, `X` escape step, `N` non-monotone step, `W` watchdog relaxed
!  step, `D` switched to accurate derivatives, `Q` QP failed, `F` no
!  acceptable step.

    pure function sqpopt_iteration_flags(info, iter_istat) result(flags)

    type(sqpopt_iter_info), intent(in) :: info       !! what happened in the iteration
    integer,                intent(in) :: iter_istat !! the iteration's status code
    character(len=12) :: flags

    flags = ''
    if (info%restoration .and. .not. info%phase) flags = trim(flags)//'R'
    if (info%phase)       flags = trim(flags)//'P'
    if (info%soc)         flags = trim(flags)//'S'
    if (info%hess_reset)  flags = trim(flags)//'H'
    if (info%elastic)     flags = trim(flags)//'E'
    if (info%escape)      flags = trim(flags)//'X'
    if (info%nonmonotone) flags = trim(flags)//'N'
    if (info%relaxed)     flags = trim(flags)//'W'
    if (info%derivatives) flags = trim(flags)//'D'
    if (info%qp_istat == sqpopt_qp_solve_failed) flags = trim(flags)//'Q'
    if (info%stepped .and. iter_istat /= sqpopt_success .and. iter_istat /= sqpopt_qp_solve_failed) then
        flags = trim(flags)//'F'
    end if

    end function sqpopt_iteration_flags
!*******************************************************************************

!*******************************************************************************
!>
!  true if every element of `v` is finite (not NaN or +/-Inf). Uses
!  `ieee_is_finite`, which (unlike an ordinary comparison) does not
!  raise an IEEE invalid exception on a NaN.

    pure logical function sqpopt_all_finite(v)

    real(wp), dimension(:), intent(in) :: v !! the values to check

    sqpopt_all_finite = all(ieee_is_finite(v))

    end function sqpopt_all_finite
!*******************************************************************************

!*******************************************************************************
!>
!  the \( \ell_1 \) constraint violation \( h(x) = \lVert \max(c_l-c,0,c-c_u)
!  \rVert_1 \) (the filter's and the funnel's \( \theta \), and the violation
!  term of the `sqpopt_merit_l1` merit function).

    pure function l1_violation(c, c_lb, c_ub) result(h)

    real(wp), dimension(:), intent(in) :: c    !! constraint values `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_lb !! their lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_ub !! their upper bounds `dimension(m)`
    real(wp) :: h

    h = sum(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))

    end function l1_violation
!*******************************************************************************

!*******************************************************************************
!>
!  the roundoff-level slack allowed when comparing a trial merit function
!  (or objective) value against the current one,
!  \( 10 \epsilon \max(1,|\phi_0|) \) (as in IPOPT's `Compare_le`).

    pure function merit_slack(phi0) result(slack)

    real(wp), intent(in) :: phi0 !! merit function (or objective) value at the current point
    real(wp) :: slack

    slack = 10.0_wp*epsilon(1.0_wp)*max(1.0_wp, abs(phi0))

    end function merit_slack
!*******************************************************************************

    end module sqpopt_types_module
!*******************************************************************************
