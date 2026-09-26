!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Merit function evaluation and line search used to globalize the SQP
!  iterations (ensures progress towards both optimality and feasibility).
!
!  Two merit functions are available (`sqpopt_linesearch_type%merit_mode`):
!
!  * `sqpopt_merit_l1` (**default**) -- the standard non-smooth \( \ell_1 \)
!    exact penalty function (as in `slsqp`).
!  * `sqpopt_merit_augmented_lagrangian` -- a smooth augmented Lagrangian
!    merit function (Gill, Murray, Saunders & Wright, *"Some Theoretical
!    Properties of an Augmented Lagrangian Merit Function"*, SOL 86-6R --
!    the merit function used in NPSOL/NPSQP and, in spirit, SNOPT). Unlike
!    the \( \ell_1 \) function, it is twice continuously differentiable,
!    which is the reason SNOPT-family solvers do not need a second-order
!    correction (see [[sqpopt_soc_module]]) to avoid the Maratos effect.
!
!  Three line search strategies are available (`sqpopt_linesearch_type%mode`):
!
!  * `sqpopt_linesearch_armijo` (**default**) -- a standard backtracking
!    line search with an Armijo-type sufficient-decrease test on the merit
!    function (as used by default in `slsqp`). An exact line search is
!    usually overkill (it requires many more function evaluations for a
!    marginal benefit), so this is the recommended/default mode.
!  * `sqpopt_linesearch_exact` -- (approximately) minimizes the merit
!    function along the search direction using the derivative-free [[fmin]]
!    routine (from the `fmin` dependency), rather than a hand-written
!    exact-search implementation.
!  * `sqpopt_linesearch_watchdog` -- Powell's watchdog technique
!    (Chamberlain, Lemarechal, Pedersen & Powell, *Math. Prog. Study 16*
!    (1982)): tracks the best point found so
!    far and, for a short window after a genuine improvement, *relaxes* the
!    sufficient-decrease test to allow the merit function to temporarily
!    get worse (accepting the full quasi-Newton step outright) rather than
!    stalling near a curved/simultaneously-active constraint boundary (the
!    Maratos effect). If the relaxed window is used up without a new best
!    point, it backtracks all the way to the best point found so far and
!    disables relaxed acceptance for `watchdog_cooldown_len` iterations.
!    This targets the same failure mode as the second-order correction and
!    the augmented Lagrangian merit function, via a different mechanism.
!  * `sqpopt_linesearch_filter` -- a **filter** line search (Fletcher &
!    Leyffer, *"Nonlinear programming without a penalty function"*, Math.
!    Program. 91 (2002)), with the globally convergent line-search rules of
!    Wächter & Biegler (*"Line search filter methods for nonlinear
!    programming: motivation and global convergence"*, SIAM J. Optim. 16
!    (2005), and the IPOPT paper, Math. Program. 106 (2006)): no merit
!    function or penalty parameter; a trial point is judged by the pair
!    \( (\theta, \varphi) \) of \( \ell_1 \) constraint violation and
!    objective value, and must not be dominated by the filter. A
!    *switching condition* decides whether a step must reduce the objective
!    (an Armijo test, "f-type" step) or may trade objective for feasibility;
!    only non-f-type steps enlarge the filter. See [[filter_line_search]].
!
!  All three of the merit-function-based modes start each search from an
!  initial trial step length
!  \( \alpha_0 \le 1 \) (see [[initial_step_length]]) rather than always
!  `1`, capped by `major_step_limit` (SNOPT's "Major step limit" option)
!  so that no variable can change by more than a factor of
!  `major_step_limit` relative to its current magnitude in a single major
!  iteration -- a QP step that is technically feasible can still be
!  unreasonably large (e.g. from a poorly-scaled problem or an early,
!  inaccurate Hessian approximation), and this guards against the
!  resulting merit-function evaluations diverging or becoming undefined.


    module sqpopt_linesearch_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_success, sqpopt_line_search_failed, sqpopt_sparse_matrix, sqpopt_all_finite
    use sqpopt_linalg_module,  only: sparse_matvec, sparse_matvec_transpose
    use fmin_module,           only: fmin

    implicit none

    private

    public :: l1_violation

    integer, parameter, public :: sqpopt_linesearch_armijo   = 1  !! backtracking Armijo-type line search (default)
    integer, parameter, public :: sqpopt_linesearch_exact    = 2  !! (approximate) exact 1-D minimization of the merit function, via [[fmin]]
    integer, parameter, public :: sqpopt_linesearch_watchdog = 3  !! Powell's watchdog technique (relaxed acceptance + backtracking, see module docs)
    integer, parameter, public :: sqpopt_linesearch_filter   = 4  !! Fletcher & Leyffer's filter method (no merit function/penalty parameter, see module docs)

    integer, parameter, public :: sqpopt_merit_l1                   = 1  !! non-smooth \( \ell_1 \) exact penalty merit function (default)
    integer, parameter, public :: sqpopt_merit_augmented_lagrangian = 2  !! smooth augmented Lagrangian merit function (NPSOL/SNOPT-style)

    abstract interface
        subroutine sqpopt_ls_objective_func(x, f)
            !! evaluates the (scaled) objective at `x` (NaN if it can't be evaluated)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)  :: x
            real(wp),               intent(out) :: f
        end subroutine sqpopt_ls_objective_func
        subroutine sqpopt_ls_constraint_func(x, c)
            !! evaluates the (scaled) constraints at `x` (NaN if they can't be evaluated)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)  :: x
            real(wp), dimension(:), intent(out) :: c
        end subroutine sqpopt_ls_constraint_func
        subroutine sqpopt_soc_func(p, c_trial, p_soc, ok)
            !! computes the second-order-corrected version `p_soc` of a
            !! rejected trial step `p`, given the constraint values
            !! `c_trial` at `x+p` (see [[sqpopt_soc_module]]); `ok` is false
            !! if no usable correction is available
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)  :: p       !! the rejected trial step `dimension(n)`
            real(wp), dimension(:), intent(in)  :: c_trial !! constraint values at `x+p` `dimension(m)`
            real(wp), dimension(:), intent(out) :: p_soc   !! the corrected step `dimension(n)`
            logical,                intent(out) :: ok      !! true if `p_soc` is usable
        end subroutine sqpopt_soc_func
    end interface
    public :: sqpopt_soc_func

    type, public :: sqpopt_linesearch_type
        !! options and state for the merit function and line search.

        integer  :: mode        = sqpopt_linesearch_armijo !! line search strategy to use
        integer  :: merit_mode  = sqpopt_merit_l1           !! merit function to use
        real(wp) :: penalty     = 1.0_wp    !! current penalty parameter used in the merit function
                                            !! (called \( \mu \) for `sqpopt_merit_l1`, \( \rho \) for
                                            !! `sqpopt_merit_augmented_lagrangian`)
        real(wp) :: tol         = 1.0e-4_wp !! desired tolerance on the minimizer (`sqpopt_linesearch_exact` mode)
        real(wp) :: sigma       = 0.1_wp    !! Armijo sufficient-decrease parameter, \( 0 < \sigma < 1 \) (`sqpopt_linesearch_armijo` mode)
        real(wp) :: backtrack   = 0.5_wp    !! step-length reduction factor at each backtracking step (`sqpopt_linesearch_armijo` mode)
        real(wp) :: alpha_min   = 1.0e-10_wp !! minimum step length (`armijo`/`watchdog`/`filter` modes): the search
                                             !! fails (no step is taken, `istat=sqpopt_line_search_failed`) if no
                                             !! acceptable step is found before `alpha` would drop below this
        integer  :: max_ls_iter = 40        !! maximum number of trial step lengths per search (`armijo`/`watchdog`/`filter` modes)
        real(wp) :: major_step_limit = 2.0_wp !! caps the *initial* trial step length (before any backtracking) so that
                                              !! no variable changes by more than this fraction of \( \max(1,|x_j|) \)
                                              !! (SNOPT's "Major step limit" option): the search starts from
                                              !! \( \alpha_0 = \min\!\left(1, \dfrac{\texttt{major\_step\_limit}}
                                              !! {\max_j |p_j|/\max(1,|x_j|)}\right) \) instead of always \( \alpha_0=1 \).
                                              !! Guards against divergence from a QP step that is technically feasible
                                              !! but unreasonably large relative to the current point; set to `huge(1.0_wp)`
                                              !! to disable

        integer  :: watchdog_relaxed_len    = 2     !! number of relaxed steps tolerated before requiring a new best point (`sqpopt_linesearch_watchdog` mode)
        integer  :: watchdog_cooldown_len   = 10    !! number of iterations relaxed acceptance is disabled for after a backtrack (`sqpopt_linesearch_watchdog` mode)

        ! filter parameters (`sqpopt_linesearch_filter` mode, and trust region with the filter), with the
        ! default values from Wächter & Biegler (see [[filter_line_search]]):
        real(wp) :: filter_gamma_theta    = 1.0e-5_wp !! margin \( \gamma_\theta \): a step must reduce \( \theta \) by this fraction...
        real(wp) :: filter_gamma_phi      = 1.0e-5_wp !! ...or reduce \( \varphi \) by \( \gamma_\varphi\theta \)
        real(wp) :: filter_delta          = 1.0_wp    !! switching condition constant \( \delta \)
        real(wp) :: filter_s_theta        = 1.1_wp    !! switching condition exponent \( s_\theta \)
        real(wp) :: filter_s_phi          = 2.3_wp    !! switching condition exponent \( s_\varphi \)
        real(wp) :: filter_eta_phi        = 1.0e-4_wp !! Armijo constant \( \eta_\varphi \) for f-type steps
        real(wp) :: filter_theta_max_fact = 1.0e4_wp  !! \( \theta_{max} = \) this \( \times \max(1,\theta_0) \): no point with a larger violation is accepted
        real(wp) :: filter_theta_min_fact = 1.0e-4_wp !! \( \theta_{min} = \) this \( \times \max(1,\theta_0) \): below it, f-type steps are allowed
        real(wp) :: filter_gamma_alpha    = 0.05_wp   !! safety factor \( \gamma_\alpha \) in the minimum step length before restoration

        ! internal state for the filter (not user options -- persists across major iterations):
        logical  :: filter_ready = .false.        !! whether the filter below has been initialized
        real(wp) :: filter_theta_max = 0.0_wp     !! \( \theta_{max} \)
        real(wp) :: filter_theta_min = 0.0_wp     !! \( \theta_{min} \)
        real(wp), dimension(:), allocatable :: filter_theta !! constraint-violation value of each filter entry
        real(wp), dimension(:), allocatable :: filter_phi   !! objective value of each filter entry

        ! internal state for `sqpopt_linesearch_watchdog` mode (not user options -- persists across major iterations):
        logical  :: watchdog_ready               = .false. !! whether the best-point tracking below has been initialized
        integer  :: watchdog_relaxed_remaining   = 0        !! iterations left in the current relaxed window
        integer  :: watchdog_cooldown_remaining  = 0        !! iterations left before relaxed acceptance may reactivate
        real(wp) :: watchdog_w_opt               = 0.0_wp   !! best merit value found so far
        real(wp), dimension(:), allocatable :: watchdog_x_opt !! best point found so far `dimension(n)`

        contains

        procedure, public :: eval_merit              => eval_merit_function
        procedure, public :: directional_derivative  => merit_directional_derivative
        procedure, public :: search                  => line_search
        procedure, public :: filter_prepare          => filter_prepare_state
        procedure, public :: filter_accept           => filter_step_acceptable
        procedure, public :: filter_record           => filter_augment

    end type sqpopt_linesearch_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  evaluate the merit function used to measure progress, dispatching on
!  `me%merit_mode`:
!
!  * `sqpopt_merit_l1`:
!    $$ \phi(x) = f(x) + \mu \lVert \max(c_l - c(x), 0, c(x) - c_u) \rVert_1 $$
!  * `sqpopt_merit_augmented_lagrangian`:
!    $$ \phi(x,\lambda,\rho) = f(x) - \lambda^T\!\left(c(x)-s\right) + \tfrac{1}{2}\rho \lVert c(x)-s \rVert_2^2 $$
!    where the slack `s` is the closed-form minimizer of \( \phi \) subject
!    to \( c_l \le s \le c_u \) (see [[augmented_lagrangian_slacks]]).

    subroutine eval_merit_function(me, f, c, c_lb, c_ub, lambda, phi)

    class(sqpopt_linesearch_type), intent(inout) :: me
    real(wp),                intent(in)  :: f      !! objective function value
    real(wp), dimension(:), intent(in)  :: c      !! constraint values `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: lambda !! Lagrange multiplier estimate `dimension(m)`
                                                  !! (only used by `sqpopt_merit_augmented_lagrangian`)
    real(wp),                intent(out) :: phi    !! value of the merit function

    real(wp), dimension(size(c)) :: s, r

    select case (me%merit_mode)
    case (sqpopt_merit_augmented_lagrangian)
        call augmented_lagrangian_slacks(me, c, c_lb, c_ub, lambda, s)
        r   = c - s
        phi = f - dot_product(lambda, r) + 0.5_wp*me%penalty*dot_product(r, r)
    case default
        phi = f + me%penalty*sum(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))
    end select

    end subroutine eval_merit_function
!*******************************************************************************

!*******************************************************************************
!>
!  the closed-form slack `s` that minimizes the augmented Lagrangian merit
!  function subject to \( c_l \le s \le c_u \): this generalizes the
!  non-negative slack \( s_i=\max(0,c_i-\lambda_i/\rho) \) of the original
!  (single-sided) formulation to `sqpopt`'s two-sided constraint bounds
!  (equality rows, where `c_lb=c_ub`, are handled automatically since `s`
!  is then clipped to that single value regardless of \( \lambda,\rho \)).

    subroutine augmented_lagrangian_slacks(me, c, c_lb, c_ub, lambda, s)

    class(sqpopt_linesearch_type), intent(in)  :: me
    real(wp), dimension(:), intent(in)  :: c !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb !! lower bounds on the constraints `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub !! upper bounds on the constraints `dimension(m)`
    real(wp), dimension(:), intent(in)  :: lambda !! Lagrange multipliers for the constraints `dimension(m)`
    real(wp), dimension(:), intent(out) :: s !! closed-form slack that minimizes the augmented Lagrangian merit function

    if (me%penalty <= 0.0_wp) then
        s = min(max(c, c_lb), c_ub)
    else
        s = min(max(c - lambda/me%penalty, c_lb), c_ub)
    end if

    end subroutine augmented_lagrangian_slacks
!*******************************************************************************

!*******************************************************************************
!>
!  the (approximate) directional derivative \( D(\phi;p) \) of the merit
!  function along `p`, used by the Armijo sufficient-decrease test and by
!  the descent-direction safeguard in [[sqpopt_iterate_module]]. Dispatches
!  on `me%merit_mode`:
!
!  * `sqpopt_merit_l1`: \( D(\phi;p) = g^Tp + \mu \, D(\lVert \text{viol}(c + \alpha Jp)
!    \rVert_1; \alpha=0^+) \), the exact one-sided derivative of the
!    linearized violation. For a step that satisfies the linearized
!    constraints this is the usual \( g^Tp - \mu \lVert \text{viol}(x)
!    \rVert_1 \); for one that doesn't (e.g. shortened by the step-length
!    cap, or an elastic step), that formula would overstate the decrease,
!    so that no step length could pass the sufficient-decrease test.
!  * `sqpopt_merit_augmented_lagrangian`: \( D(\phi;p) = (g - J^T\lambda +
!    \rho J^T(c-s))^Tp \), holding \( \lambda \) and `s` fixed at their
!    current values (a simplification of the full NPSQP theory, which
!    differentiates along a joint `(x,lambda,s)` step; adequate since
!    `sqpopt` only updates `lambda` once per major iteration).

    subroutine merit_directional_derivative(me, jac, g, p, c, c_lb, c_ub, lambda, dphi0)

    class(sqpopt_linesearch_type), intent(in) :: me
    type(sqpopt_sparse_matrix), intent(in) :: jac  !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:), intent(in) :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(in) :: p        !! search direction `dimension(n)`
    real(wp), dimension(:), intent(in) :: c        !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_lb     !! lower bounds on the constraints `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_ub     !! upper bounds on the constraints `dimension(m)`
    real(wp), dimension(:), intent(in) :: lambda   !! Lagrange multipliers for the constraints `dimension(m)`
    real(wp), intent(out) :: dphi0 !! directional derivative of the merit function along `p`

    real(wp), dimension(size(c)) :: s !! slack variables for the augmented Lagrangian
    real(wp), dimension(size(g)) :: jtlam, jtr

    select case (me%merit_mode)
    case (sqpopt_merit_augmented_lagrangian)
        call augmented_lagrangian_slacks(me, c, c_lb, c_ub, lambda, s)
        call sparse_matvec_transpose(jac, lambda, jtlam)
        call sparse_matvec_transpose(jac, c-s, jtr)
        dphi0 = dot_product(g - jtlam + me%penalty*jtr, p)
    case default
        block
            real(wp), dimension(size(c)) :: jp
            real(wp) :: rate
            integer :: i
            call sparse_matvec(jac, p, jp)
            rate = 0.0_wp
            do i = 1, size(c)
                if (c(i) < c_lb(i)) then
                    rate = rate - jp(i)
                else if (c(i) > c_ub(i)) then
                    rate = rate + jp(i)
                else if (c(i) == c_lb(i) .and. jp(i) < 0.0_wp) then
                    rate = rate - jp(i)
                else if (c(i) == c_ub(i) .and. jp(i) > 0.0_wp) then
                    rate = rate + jp(i)
                end if
            end do
            dphi0 = dot_product(g, p) + me%penalty*rate
        end block
    end select

    end subroutine merit_directional_derivative
!*******************************************************************************

!*******************************************************************************
!>
!  perform a line search along the direction `p` to find an accepted step
!  length `alpha`, dispatching to the strategy selected by `me%mode`.
!
!  If no acceptable step is found, `istat=sqpopt_line_search_failed` and
!  `x_new=x` (no step is taken; in `sqpopt_linesearch_watchdog` mode a
!  failed relaxed window instead returns the best point found so far).
!  Trial points where `f` or `c` is not finite (NaN or Inf) are always
!  rejected, so the search backtracks away from them.
!
!  If `soc` is present, a second-order correction of the *first* trial
!  step is also tried when that step is rejected and did not reduce the
!  constraint violation (see [[sqpopt_soc_module]]); in
!  `sqpopt_linesearch_exact` mode, when the minimizer along `p` falls
!  short of the full step and the full step did not reduce the violation.

    subroutine line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat, soc)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_ls_objective_func)  :: eval_f  !! evaluates \( f(x) \)
    procedure(sqpopt_ls_constraint_func) :: eval_c  !! evaluates \( c(x) \)
    real(wp), dimension(:), intent(in)  :: x      !! current point `dimension(n)`
    real(wp), dimension(:), intent(in)  :: p      !! search direction `dimension(n)`
    real(wp),                intent(in)  :: f      !! objective function value at `x`
    real(wp), dimension(:), intent(in)  :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(in)  :: c      !! constraint values at `x` `dimension(m)`
    type(sqpopt_sparse_matrix), intent(in) :: jac !! constraint Jacobian at `x`, `dimension(m,n)`
                                                  !! (only used by `sqpopt_merit_augmented_lagrangian`)
    real(wp), dimension(:), intent(in)  :: lambda !! Lagrange multiplier estimate `dimension(m)`
                                                  !! (only used by `sqpopt_merit_augmented_lagrangian`)
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp),                intent(out) :: alpha  !! accepted step length (`0` if no step was taken)
    real(wp), dimension(:), intent(out) :: x_new   !! the accepted new point `dimension(n)` (normally
                                                    !! `x + alpha*p`; see above for the exceptions)
    integer,                  intent(out) :: istat  !! status code (see [[sqpopt_types_module]])
    procedure(sqpopt_soc_func), optional :: soc     !! computes a second-order-corrected step

    select case (me%mode)
    case (sqpopt_linesearch_exact)
        call exact_line_search(me, eval_f, eval_c, x, p, f, c, lambda, c_lb, c_ub, alpha, x_new, istat, soc)
    case (sqpopt_linesearch_watchdog)
        call watchdog_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat, soc)
    case (sqpopt_linesearch_filter)
        call filter_line_search(me, eval_f, eval_c, x, p, f, g, c, c_lb, c_ub, alpha, x_new, istat, soc)
    case default
        call armijo_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat, soc)
    end select

    end subroutine line_search
!*******************************************************************************

!*******************************************************************************
!>
!  the initial trial step length used by every line search mode before any
!  backtracking (SNOPT's "Major step limit" option):
!  $$ \alpha_0 = \min\!\left(1, \frac{\texttt{step\_limit}}
!  {\max_j |p_j|/\max(1,|x_j|)}\right) $$
!  i.e. `p` is scaled down (never up) so that no variable changes by more
!  than a factor of `step_limit` relative to \( \max(1,|x_j|) \). Returns
!  `1` unchanged if `p` is zero.

    function initial_step_length(x, p, step_limit) result(alpha0)

    real(wp), dimension(:), intent(in) :: x, p !! current point and search direction, `dimension(n)`
    real(wp),                intent(in) :: step_limit !! `sqpopt_linesearch_type%major_step_limit`
    real(wp) :: alpha0 !! initial trial step length

    real(wp) :: rmax
    integer  :: k

    rmax = 0.0_wp
    do k = 1, size(x)
        rmax = max(rmax, abs(p(k))/max(1.0_wp, abs(x(k))))
    end do

    ! (written so that `step_limit=huge(1.0_wp)` can't overflow)
    if (rmax > step_limit) then
        alpha0 = step_limit/rmax
    else
        alpha0 = 1.0_wp
    end if

    end function initial_step_length
!*******************************************************************************

!*******************************************************************************
!>
!  the roundoff-level slack allowed when comparing a trial merit function
!  value against the current one, \( 10 \epsilon \max(1,|\phi_0|) \) (as
!  in IPOPT's `Compare_le`).

    pure function merit_slack(phi0) result(slack)

    real(wp), intent(in) :: phi0 !! merit function value at the current point
    real(wp) :: slack

    slack = 10.0_wp*epsilon(1.0_wp)*max(1.0_wp, abs(phi0))

    end function merit_slack
!*******************************************************************************

!*******************************************************************************
!>
!  evaluate `f` and `c` at a trial point; `ok` is false if any value is
!  not finite (NaN or Inf), in which case the point must be rejected.

    subroutine eval_fc(eval_f, eval_c, x_trial, f_trial, c_trial, ok)

    procedure(sqpopt_ls_objective_func)  :: eval_f
    procedure(sqpopt_ls_constraint_func) :: eval_c
    real(wp), dimension(:), intent(in)  :: x_trial !! trial point `dimension(n)`
    real(wp),               intent(out) :: f_trial !! objective function value at `x_trial`
    real(wp), dimension(:), intent(out) :: c_trial !! constraint values at `x_trial` `dimension(m)`
    logical,                intent(out) :: ok      !! true if `f_trial` and `c_trial` are all finite

    call eval_f(x_trial, f_trial)
    call eval_c(x_trial, c_trial)
    ok = sqpopt_all_finite([f_trial]) .and. sqpopt_all_finite(c_trial)

    end subroutine eval_fc
!*******************************************************************************

!*******************************************************************************
!>
!  evaluate `f`, `c`, and the merit function at a trial point; `ok` is
!  false (and `phi_trial=huge`) if any value is not finite.

    subroutine eval_trial(me, eval_f, eval_c, x_trial, c_lb, c_ub, lambda, c_trial, phi_trial, ok)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_ls_objective_func)  :: eval_f
    procedure(sqpopt_ls_constraint_func) :: eval_c
    real(wp), dimension(:), intent(in)  :: x_trial   !! trial point `dimension(n)`
    real(wp), dimension(:), intent(in)  :: c_lb      !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub      !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: lambda    !! Lagrange multiplier estimate `dimension(m)`
    real(wp), dimension(:), intent(out) :: c_trial   !! constraint values at `x_trial` `dimension(m)`
    real(wp),               intent(out) :: phi_trial !! merit function value at `x_trial`
    logical,                intent(out) :: ok        !! true if everything is finite

    real(wp) :: f_trial

    call eval_fc(eval_f, eval_c, x_trial, f_trial, c_trial, ok)
    if (ok) then
        call me%eval_merit(f_trial, c_trial, c_lb, c_ub, lambda, phi_trial)
        ok = sqpopt_all_finite([phi_trial])
    end if
    if (.not. ok) phi_trial = huge(1.0_wp)

    end subroutine eval_trial
!*******************************************************************************

!*******************************************************************************
!>
!  the backtracking Armijo search shared by [[armijo_line_search]] and
!  [[watchdog_line_search]]: starting from `alpha0`, `alpha` is reduced by
!  `backtrack` until
!  $$ \phi(x+\alpha p) \le \phi(x) + \sigma \alpha \min(D(\phi;p),0) $$
!  (the `min` keeps a non-descent `p` from loosening the test), for at most
!  `max_ls_iter` trials and not below `alpha_min`. The comparison allows a
!  roundoff-level slack (see [[merit_slack]]), so that near a solution,
!  where the change in the merit function is at the level of rounding
!  error, a good step isn't rejected (which would otherwise backtrack
!  `alpha` to a uselessly tiny value). If the first trial is
!  rejected without reducing the constraint violation and `soc` is
!  present, the second-order-corrected step is also tried at that point.

    subroutine backtrack_search(me, eval_f, eval_c, x, p, alpha0, phi0, dphi0, c, c_lb, c_ub, lambda, &
                                alpha, x_new, phi_new, accepted, soc)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_ls_objective_func)  :: eval_f
    procedure(sqpopt_ls_constraint_func) :: eval_c
    real(wp), dimension(:), intent(in)  :: x        !! current point `dimension(n)`
    real(wp), dimension(:), intent(in)  :: p        !! search direction `dimension(n)`
    real(wp),               intent(in)  :: alpha0   !! initial trial step length
    real(wp),               intent(in)  :: phi0     !! merit function value at `x`
    real(wp),               intent(in)  :: dphi0    !! directional derivative of the merit function along `p`
    real(wp), dimension(:), intent(in)  :: c        !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb     !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub     !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: lambda   !! Lagrange multiplier estimate `dimension(m)`
    real(wp),               intent(out) :: alpha    !! accepted step length (meaningful only if `accepted`)
    real(wp), dimension(:), intent(out) :: x_new    !! accepted point `dimension(n)` (meaningful only if `accepted`)
    real(wp),               intent(out) :: phi_new  !! merit function value at `x_new` (meaningful only if `accepted`)
    logical,                intent(out) :: accepted !! true if a step satisfying the Armijo test was found
    procedure(sqpopt_soc_func), optional :: soc     !! computes a second-order-corrected step

    real(wp), dimension(size(x)) :: x_trial, p_soc
    real(wp), dimension(size(c)) :: c_trial
    real(wp) :: phi_trial, slope, slack
    logical  :: ok, soc_ok
    integer  :: it

    slope    = min(dphi0, 0.0_wp)
    slack    = merit_slack(phi0)
    accepted = .false.
    alpha    = alpha0

    do it = 1, me%max_ls_iter

        x_trial = x + alpha*p
        call eval_trial(me, eval_f, eval_c, x_trial, c_lb, c_ub, lambda, c_trial, phi_trial, ok)
        if (ok .and. phi_trial <= phi0 + me%sigma*alpha*slope + slack) then
            accepted = .true.
            exit
        end if

        if (it == 1 .and. ok .and. present(soc)) then
            if (l1_violation(c_trial, c_lb, c_ub) >= l1_violation(c, c_lb, c_ub)) then
                ! the full step didn't reduce the constraint violation, so its
                ! rejection may be due to constraint curvature (the Maratos
                ! effect): try the second-order-corrected step before backtracking:
                call soc(alpha*p, c_trial, p_soc, soc_ok)
                if (soc_ok) then
                    x_trial = x + p_soc
                    call eval_trial(me, eval_f, eval_c, x_trial, c_lb, c_ub, lambda, c_trial, phi_trial, ok)
                    if (ok .and. phi_trial <= phi0 + me%sigma*alpha*slope + slack) then
                        accepted = .true.
                        exit
                    end if
                end if
            end if
        end if

        alpha = me%backtrack*alpha
        if (alpha < me%alpha_min) exit

    end do

    if (accepted) then
        x_new   = x_trial
        phi_new = phi_trial
    end if

    end subroutine backtrack_search
!*******************************************************************************

!*******************************************************************************
!>
!  backtracking line search with an Armijo-type sufficient-decrease test
!  on the merit function (see [[backtrack_search]]), starting from
!  \( \alpha_0 \) (see [[initial_step_length]], normally `1` unless capped
!  by `major_step_limit`). If no acceptable step is found, no step is taken
!  (`x_new=x`, `istat=sqpopt_line_search_failed`).

    subroutine armijo_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat, soc)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_ls_objective_func)  :: eval_f
    procedure(sqpopt_ls_constraint_func) :: eval_c
    real(wp), dimension(:), intent(in) :: x !! current point `x`
    real(wp), dimension(:), intent(in) :: p !! search direction `p`
    real(wp), dimension(:), intent(in) :: g !! gradient of the objective at `x`
    real(wp), dimension(:), intent(in) :: c !! constraint values at `x`
    real(wp), dimension(:), intent(in) :: lambda !! Lagrange multipliers at `x`
    real(wp), dimension(:), intent(in) :: c_lb   !! lower bounds on the constraints
    real(wp), dimension(:), intent(in) :: c_ub   !! upper bounds on the constraints
    type(sqpopt_sparse_matrix), intent(in) :: jac    !! constraint Jacobian at `x` (`dimension(m,n)`)
    real(wp),                   intent(in)  :: f     !! objective function value at `x`
    real(wp),                   intent(out) :: alpha !! step length along `p` (`0` if no step was taken)
    real(wp), dimension(:),     intent(out) :: x_new !! new point `x + alpha*p` (or the corrected step)
    integer,                    intent(out) :: istat !! status of the line search (success or failure)
    procedure(sqpopt_soc_func), optional :: soc      !! computes a second-order-corrected step

    real(wp) :: phi0, dphi0, phi_new
    logical  :: accepted

    call me%eval_merit(f, c, c_lb, c_ub, lambda, phi0)
    call me%directional_derivative(jac, g, p, c, c_lb, c_ub, lambda, dphi0)

    call backtrack_search(me, eval_f, eval_c, x, p, initial_step_length(x, p, me%major_step_limit), &
                          phi0, dphi0, c, c_lb, c_ub, lambda, alpha, x_new, phi_new, accepted, soc)

    if (accepted) then
        istat = sqpopt_success
    else
        alpha = 0.0_wp
        x_new = x
        istat = sqpopt_line_search_failed
    end if

    end subroutine armijo_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  (approximately) minimize the merit function along `p` using the
!  derivative-free 1-D minimizer [[fmin]]. The step is only accepted if it
!  does not increase the merit function by more than roundoff (see
!  [[merit_slack]]; otherwise no step is taken and
!  `istat=sqpopt_line_search_failed`); non-finite trial values are treated
!  as `huge` by the 1-D minimization. If the minimizer falls short of the
!  (initial) full step, and the full step did not reduce the constraint
!  violation (the signature of the Maratos effect, which slows the exact
!  search to linear convergence), the second-order-corrected full step is
!  also tried (if `soc` is present), and taken if its merit is lower.

    subroutine exact_line_search(me, eval_f, eval_c, x, p, f, c, lambda, c_lb, c_ub, alpha, x_new, istat, soc)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_ls_objective_func)  :: eval_f  !! evaluates \( f(x) \)
    procedure(sqpopt_ls_constraint_func) :: eval_c  !! evaluates \( c(x) \)
    real(wp), dimension(:), intent(in)  :: x      !! current point `dimension(n)`
    real(wp), dimension(:), intent(in)  :: p      !! search direction `dimension(n)`
    real(wp),               intent(in)  :: f      !! objective function value at `x`
    real(wp), dimension(:), intent(in)  :: c      !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in)  :: lambda !! Lagrange multiplier estimate `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp),                intent(out) :: alpha  !! accepted step length (`0` if no step was taken)
    real(wp), dimension(:), intent(out) :: x_new   !! the accepted new point `dimension(n)`
    integer,                  intent(out) :: istat  !! status code (see [[sqpopt_types_module]])
    procedure(sqpopt_soc_func), optional :: soc     !! computes a second-order-corrected step

    real(wp) :: phi0, phi, alpha0, phi_full, phi_soc
    real(wp), dimension(size(c)) :: c_full, c_soc
    real(wp), dimension(size(x)) :: p_soc
    logical :: ok, soc_ok

    call me%eval_merit(f, c, c_lb, c_ub, lambda, phi0)
    alpha0 = initial_step_length(x, p, me%major_step_limit)
    alpha  = fmin(merit_along_direction, 0.0_wp, alpha0, me%tol)
    phi    = merit_along_direction(alpha)
    x_new  = x + alpha*p

    if (present(soc) .and. alpha < alpha0*(1.0_wp - sqrt(me%tol))) then
        call eval_trial(me, eval_f, eval_c, x + alpha0*p, c_lb, c_ub, lambda, c_full, phi_full, ok)
        if (ok) then
            if (l1_violation(c_full, c_lb, c_ub) >= l1_violation(c, c_lb, c_ub)) then
                call soc(alpha0*p, c_full, p_soc, soc_ok)
                if (soc_ok) then
                    call eval_trial(me, eval_f, eval_c, x + p_soc, c_lb, c_ub, lambda, c_soc, phi_soc, ok)
                    if (ok .and. phi_soc < phi) then
                        phi   = phi_soc
                        alpha = alpha0
                        x_new = x + p_soc
                    end if
                end if
            end if
        end if
    end if

    if (phi < phi0 + merit_slack(phi0)) then
        istat = sqpopt_success
    else
        alpha = 0.0_wp
        x_new = x
        istat = sqpopt_line_search_failed
    end if

    contains

    !*******************************************************************************
    !>
    !  the merit function \( \phi(x + \alpha p) \) along the search direction,
    !  in the form required by [[fmin]] (`huge` at a non-finite trial point).
    !  Uses the host-associated variables set by [[exact_line_search]].

        function merit_along_direction(alpha) result(phi)

        real(wp), intent(in) :: alpha !! step length along the search direction
        real(wp) :: phi !! merit function value at the trial point

        real(wp), dimension(size(c_lb)) :: c_trial
        logical :: ok

        call eval_trial(me, eval_f, eval_c, x + alpha*p, c_lb, c_ub, lambda, c_trial, phi, ok)

        end function merit_along_direction
    !*******************************************************************************

    end subroutine exact_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  Powell's watchdog technique (Chamberlain, Lemarechal, Pedersen & Powell,
!  *Math. Prog. Study 16* (1982)): a variant of
!  [[armijo_line_search]] that, once a genuine improvement has been made,
!  allows a bounded number of subsequent *relaxed* steps -- accepting the
!  full step `x+p` outright even if it does not satisfy the sufficient-
!  decrease test (but never a step to a non-finite point) -- rather than
!  stalling near a curved or simultaneously-active constraint boundary (the
!  Maratos effect). If none of those relaxed steps beats the best point
!  found so far, the search backtracks all the way to that best point and
!  disables relaxed acceptance for `watchdog_cooldown_len` further calls.
!  If the standard search fails with no relaxed window available, no step
!  is taken.

    subroutine watchdog_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat, soc)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_ls_objective_func)  :: eval_f
    procedure(sqpopt_ls_constraint_func) :: eval_c
    real(wp), dimension(:), intent(in)  :: x !! current point in the search space
    real(wp), dimension(:), intent(in)  :: p !! search direction
    real(wp), dimension(:), intent(in)  :: g !! gradient of the objective function at `x`
    real(wp), dimension(:), intent(in)  :: c !! constraint function values at `x`
    real(wp), dimension(:), intent(in)  :: lambda !! Lagrange multipliers at `x`
    real(wp), dimension(:), intent(in)  :: c_lb !! lower bounds on the constraints
    real(wp), dimension(:), intent(in)  :: c_ub !! upper bounds on the constraints
    type(sqpopt_sparse_matrix), intent(in)  :: jac !! Jacobian of the constraints at `x`
    real(wp),                   intent(in)  :: f !! objective function value at `x`
    real(wp),                   intent(out) :: alpha !! step length found by the line search
    real(wp), dimension(:),     intent(out) :: x_new !! new point after the line search
    integer,                    intent(out) :: istat !! status of the line search (0 if successful)
    procedure(sqpopt_soc_func), optional :: soc      !! computes a second-order-corrected step

    real(wp), dimension(size(c)) :: c_trial !! constraint values at the trial point
    real(wp) :: phi0, dphi0, phi_trial, alpha0 !! merit function values, directional derivative, initial step length
    logical :: standard_ok, relaxed_used, ok !! flags indicating if standard or relaxed line search succeeded

    call me%eval_merit(f, c, c_lb, c_ub, lambda, phi0)
    call me%directional_derivative(jac, g, p, c, c_lb, c_ub, lambda, dphi0)

    ! initialize the best-point-so-far tracking the first time this is called:
    if (.not. me%watchdog_ready) then
        me%watchdog_x_opt = x
        me%watchdog_w_opt = phi0
        me%watchdog_ready = .true.
    end if
    if (me%watchdog_cooldown_remaining > 0) me%watchdog_cooldown_remaining = me%watchdog_cooldown_remaining - 1

    ! standard backtracking Armijo search, exactly as in [[armijo_line_search]]:
    alpha0 = initial_step_length(x, p, me%major_step_limit)
    call backtrack_search(me, eval_f, eval_c, x, p, alpha0, phi0, dphi0, c, c_lb, c_ub, lambda, &
                          alpha, x_new, phi_trial, standard_ok, soc)

    relaxed_used = .false.
    if (.not. standard_ok .and. me%watchdog_relaxed_remaining > 0 .and. me%watchdog_cooldown_remaining == 0) then
        ! the standard sufficient-decrease test failed; the watchdog
        ! technique allows a relaxed step here instead of stalling (the
        ! merit function may temporarily get worse), still capped by the
        ! major step limit like the initial standard-search step:
        alpha = alpha0
        x_new = x + alpha*p
        call eval_trial(me, eval_f, eval_c, x_new, c_lb, c_ub, lambda, c_trial, phi_trial, ok)
        relaxed_used = ok
    end if

    if (.not. (standard_ok .or. relaxed_used)) then
        ! plain Armijo search failed, with no (usable) relaxed window: no step
        alpha = 0.0_wp
        x_new = x
        istat = sqpopt_line_search_failed
        return
    end if

    istat = sqpopt_success

    block
        logical :: is_new_best
        is_new_best = phi_trial < me%watchdog_w_opt
        if (is_new_best) then
            ! update the safety net whenever there is genuine improvement:
            me%watchdog_x_opt = x_new
            me%watchdog_w_opt = phi_trial
        end if

        if (standard_ok .and. alpha >= 0.99_wp*alpha0) then
            ! a good, (nearly) full accepted step: "reward" the next few
            ! calls with a fresh window of relaxed acceptance (only a step
            ! that met the strict sufficient-decrease test at close to full
            ! length earns this, not just any improvement):
            me%watchdog_relaxed_remaining = me%watchdog_relaxed_len
        else if (relaxed_used) then
            if (is_new_best) then
                ! the relaxed step paid off with real progress: keep the
                ! window open for another attempt:
                me%watchdog_relaxed_remaining = me%watchdog_relaxed_len
            else
                ! a relaxed step was taken but did not beat the best point
                ! so far; consume one attempt from the relaxed budget:
                me%watchdog_relaxed_remaining = me%watchdog_relaxed_remaining - 1
                if (me%watchdog_relaxed_remaining <= 0) then
                    ! the relaxed attempts are exhausted without improvement:
                    ! back-track all the way to the best point found so far,
                    ! and disable relaxed acceptance for
                    ! `watchdog_cooldown_len` further calls:
                    x_new = me%watchdog_x_opt
                    alpha = 0.0_wp
                    me%watchdog_cooldown_remaining = me%watchdog_cooldown_len
                    istat = sqpopt_line_search_failed
                end if
            end if
        end if
    end block

    end subroutine watchdog_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  the filter line search (see the module-level documentation), with the
!  rules of Wächter & Biegler (2005, 2006). With
!  \( \theta \) the \( \ell_1 \) constraint violation and \( \varphi=f \),
!  a trial point \( x+\alpha p \) is accepted if it is acceptable to the
!  filter (see [[filter_step_acceptable]]) and either
!
!  * (f-type step) the *switching condition*
!    \( g^Tp<0 \) and \( \alpha(-g^Tp)^{s_\varphi} > \delta\theta_k^{s_\theta} \)
!    holds and \( \theta_k \le \theta_{min} \), and the Armijo condition
!    \( \varphi(x+\alpha p) \le \varphi_k + \eta_\varphi \alpha g^Tp \)
!    holds; or
!  * otherwise, it sufficiently reduces the violation or the objective:
!    \( \theta \le (1-\gamma_\theta)\theta_k \) or
!    \( \varphi \le \varphi_k - \gamma_\varphi\theta_k \).
!
!  Only a step that is *not* f-type adds \( ((1-\gamma_\theta)\theta_k,
!  \varphi_k-\gamma_\varphi\theta_k) \) to the filter. `alpha` is
!  backtracked until acceptance, or until it falls below
!  \( \alpha_{min} \) (Wächter & Biegler's formula, the smallest step at
!  which acceptance is still possible), in which case the search fails
!  (`istat=sqpopt_line_search_failed`, `x_new=x`); [[sqpopt_iterate_module]]
!  then adds the current point to the filter and takes a feasibility
!  restoration step. If the first trial is rejected without reducing
!  \( \theta \), its second-order correction is also tried (if `soc` is
!  present), with the same acceptance test.

    subroutine filter_line_search(me, eval_f, eval_c, x, p, f, g, c, c_lb, c_ub, alpha, x_new, istat, soc)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_ls_objective_func)  :: eval_f
    procedure(sqpopt_ls_constraint_func) :: eval_c
    real(wp), dimension(:), intent(in)  :: x      !! current point `dimension(n)`
    real(wp), dimension(:), intent(in)  :: p      !! search direction `dimension(n)`
    real(wp),               intent(in)  :: f      !! objective function value at `x`
    real(wp), dimension(:), intent(in)  :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(in)  :: c      !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp),               intent(out) :: alpha  !! accepted step length (`0` if no step was taken)
    real(wp), dimension(:), intent(out) :: x_new  !! the accepted new point `dimension(n)`
    integer,                intent(out) :: istat  !! status code (see [[sqpopt_types_module]])
    procedure(sqpopt_soc_func), optional :: soc   !! computes a second-order-corrected step

    real(wp), dimension(size(x)) :: x_trial, p_soc
    real(wp), dimension(size(c)) :: c_trial
    real(wp) :: f_trial, theta0, theta_t, gtp, alpha_lim
    logical :: ok, soc_ok, f_type
    integer :: it

    theta0 = l1_violation(c, c_lb, c_ub)
    call me%filter_prepare(theta0)
    gtp = dot_product(g, p)
    alpha_lim = max(me%alpha_min, filter_min_step(me, theta0, gtp))

    alpha = initial_step_length(x, p, me%major_step_limit)
    do it = 1, me%max_ls_iter

        x_trial = x + alpha*p
        call eval_fc(eval_f, eval_c, x_trial, f_trial, c_trial, ok)
        if (ok) then
            theta_t = l1_violation(c_trial, c_lb, c_ub)
            if (me%filter_accept(theta0, f, gtp, alpha, theta_t, f_trial, f_type)) then
                call accept()
                return
            end if
        end if

        if (it == 1 .and. ok .and. present(soc)) then
            if (theta_t >= theta0) then
                ! the full step didn't reduce the constraint violation: try
                ! the second-order-corrected step before backtracking:
                call soc(alpha*p, c_trial, p_soc, soc_ok)
                if (soc_ok) then
                    x_trial = x + p_soc
                    call eval_fc(eval_f, eval_c, x_trial, f_trial, c_trial, ok)
                    if (ok) then
                        theta_t = l1_violation(c_trial, c_lb, c_ub)
                        if (me%filter_accept(theta0, f, gtp, alpha, theta_t, f_trial, f_type)) then
                            call accept()
                            return
                        end if
                    end if
                end if
            end if
        end if

        alpha = me%backtrack*alpha
        if (alpha < alpha_lim) exit

    end do

    ! no acceptable point was found: no step is taken
    alpha = 0.0_wp
    x_new = x
    istat = sqpopt_line_search_failed

    contains

        subroutine accept()
        !! accept `x_trial`, augmenting the filter unless it was an f-type step
        if (.not. f_type) call me%filter_record(theta0, f)
        x_new = x_trial
        istat = sqpopt_success
        end subroutine accept

    end subroutine filter_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the filter (empty) and its bounds \( \theta_{max} \),
!  \( \theta_{min} \) from the constraint violation `theta0` at the
!  starting point, the first time it is used; a no-op after that. Shared
!  by [[filter_line_search]] and the trust region's filter acceptance, so
!  both use the same filter.

    subroutine filter_prepare_state(me, theta0)

    class(sqpopt_linesearch_type), intent(inout) :: me
    real(wp),                      intent(in)    :: theta0 !! constraint violation at the starting point

    if (me%filter_ready) return

    me%filter_theta_max = me%filter_theta_max_fact*max(1.0_wp, theta0)
    me%filter_theta_min = me%filter_theta_min_fact*max(1.0_wp, theta0)
    if (allocated(me%filter_theta)) deallocate(me%filter_theta)
    if (allocated(me%filter_phi))   deallocate(me%filter_phi)
    allocate(me%filter_theta(0), me%filter_phi(0))
    me%filter_ready = .true.

    end subroutine filter_prepare_state
!*******************************************************************************

!*******************************************************************************
!>
!  the \( \ell_1 \) constraint violation \( h(x) = \lVert \max(c_l-c,0,c-c_u)
!  \rVert_1 \) (the filter's \( \theta \), and the violation term of the
!  `sqpopt_merit_l1` merit function).

    pure function l1_violation(c, c_lb, c_ub) result(h)

    real(wp), dimension(:), intent(in) :: c, c_lb, c_ub
    real(wp) :: h

    h = sum(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))

    end function l1_violation
!*******************************************************************************

!*******************************************************************************
!>
!  whether a trial point with violation `theta_t` and objective `phi_t`,
!  reached with step length `alpha` from the current point
!  (`theta_k`, `phi_k`), is acceptable (see [[filter_line_search]] for the
!  rules). `gtp` is the predicted change in the objective along the full
!  step (\( g^Tp \) for the line search; the trust region passes the
!  quadratic model's \( -q \), with `alpha=1`). `f_type` is set if the
!  switching condition held (so the step was judged on the objective alone,
!  and must not enlarge the filter). The trial point must also be
!  acceptable to the filter itself: \( \theta_t < \theta_{max} \), and not
!  dominated by any filter entry (\( \theta_t < \theta_j \) or
!  \( \varphi_t < \varphi_j \) for every entry `j`).

    function filter_step_acceptable(me, theta_k, phi_k, gtp, alpha, theta_t, phi_t, f_type) result(ok)

    class(sqpopt_linesearch_type), intent(in) :: me
    real(wp), intent(in)  :: theta_k, phi_k   !! violation and objective at the current point
    real(wp), intent(in)  :: gtp              !! predicted change in the objective along the full step
    real(wp), intent(in)  :: alpha            !! step length
    real(wp), intent(in)  :: theta_t, phi_t   !! violation and objective at the trial point
    logical,  intent(out) :: f_type           !! true if the switching condition held
    logical :: ok

    integer :: j

    ! switching condition, alpha*(-gtp)^s_phi > delta*theta_k^s_theta
    ! (compared in log space, so the powers can't overflow):
    f_type = gtp < 0.0_wp .and. theta_k <= me%filter_theta_min
    if (f_type .and. theta_k > 0.0_wp) &
        f_type = log(alpha) + me%filter_s_phi*log(-gtp) > log(me%filter_delta) + me%filter_s_theta*log(theta_k)

    if (f_type) then
        ! Armijo condition on the objective (with a roundoff-level slack):
        ok = phi_t <= phi_k + me%filter_eta_phi*alpha*gtp + merit_slack(phi_k)
    else
        ! sufficient reduction of the violation or the objective:
        ok = theta_t <= (1.0_wp - me%filter_gamma_theta)*theta_k .or. &
             phi_t   <= phi_k - me%filter_gamma_phi*theta_k
    end if
    if (.not. ok) return

    ! acceptable to the filter:
    ok = theta_t < me%filter_theta_max
    if (.not. ok) return
    do j = 1, size(me%filter_theta)
        if (theta_t >= me%filter_theta(j) .and. phi_t >= me%filter_phi(j)) then
            ok = .false.
            return
        end if
    end do

    end function filter_step_acceptable
!*******************************************************************************

!*******************************************************************************
!>
!  augment the filter with the current point (`theta_k`, `phi_k`), with
!  margins: the entry \( ((1-\gamma_\theta)\theta_k,
!  \varphi_k-\gamma_\varphi\theta_k) \) is added, and any existing entries
!  it dominates are removed.

    subroutine filter_augment(me, theta_k, phi_k)

    class(sqpopt_linesearch_type), intent(inout) :: me
    real(wp),                      intent(in)    :: theta_k, phi_k !! violation and objective at the current point

    real(wp) :: theta_new, phi_new
    logical, dimension(:), allocatable :: keep

    if (.not. me%filter_ready) call me%filter_prepare(theta_k)

    theta_new = (1.0_wp - me%filter_gamma_theta)*theta_k
    phi_new   = phi_k - me%filter_gamma_phi*theta_k

    if (size(me%filter_theta) > 0) then
        keep = .not. (me%filter_theta >= theta_new .and. me%filter_phi >= phi_new)
        me%filter_theta = pack(me%filter_theta, keep)
        me%filter_phi   = pack(me%filter_phi,   keep)
    end if
    me%filter_theta = [me%filter_theta, theta_new]
    me%filter_phi   = [me%filter_phi,   phi_new]

    end subroutine filter_augment
!*******************************************************************************

!*******************************************************************************
!>
!  Wächter & Biegler's minimum step length \( \alpha_{min} \): below it,
!  no trial step can be acceptable (so the line search should give up and
!  go to feasibility restoration).

    pure function filter_min_step(me, theta_k, gtp) result(alpha_min)

    class(sqpopt_linesearch_type), intent(in) :: me
    real(wp), intent(in) :: theta_k !! violation at the current point
    real(wp), intent(in) :: gtp     !! \( g^Tp \)
    real(wp) :: alpha_min

    if (gtp < 0.0_wp) then
        alpha_min = min(me%filter_gamma_theta, me%filter_gamma_phi*theta_k/(-gtp))
        if (theta_k <= me%filter_theta_min) then
            if (theta_k > 0.0_wp) then
                alpha_min = min(alpha_min, exp(log(me%filter_delta) + me%filter_s_theta*log(theta_k) &
                                               - me%filter_s_phi*log(-gtp)))
            else
                alpha_min = 0.0_wp
            end if
        end if
    else
        alpha_min = me%filter_gamma_theta
    end if
    alpha_min = me%filter_gamma_alpha*alpha_min

    end function filter_min_step
!*******************************************************************************

    end module sqpopt_linesearch_module
!*******************************************************************************
