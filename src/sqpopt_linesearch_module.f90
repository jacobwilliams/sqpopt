!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The line searches used to globalize the SQP iterations (ensures
!  progress towards both optimality and feasibility). The acceptance tests
!  they use are in their own modules, held here as components: the merit
!  function and its penalty parameter (`merit`, see
!  [[sqpopt_merit_module]]), the filter (`filter`, see
!  [[sqpopt_filter_module]]), and the funnel (`funnel`, see
!  [[sqpopt_funnel_module]]). The trust region (see
!  [[sqpopt_trust_region_module]]) uses the same components, so the two
!  globalizations share the filter or funnel.
!
!  Five line search strategies are available (`sqpopt_linesearch_type%mode`):
!
!  * `sqpopt_linesearch_armijo` -- a standard backtracking line search
!    with an Armijo-type sufficient-decrease test on the merit function (as
!    used by default in `slsqp`).
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
!  * `sqpopt_linesearch_filter` (**default**) -- a **filter** line search (Fletcher &
!    Leyffer, *"Nonlinear programming without a penalty function"*, Math.
!    Program. 91 (2002)), with the globally convergent line-search rules of
!    Wächter & Biegler (*"Line search filter methods for nonlinear
!    programming: motivation and global convergence"*, SIAM J. Optim. 16
!    (2005), and the IPOPT paper, Math. Program. 106 (2006)): no merit
!    function or penalty parameter; a trial point is judged by the pair
!    \( (\theta, \varphi) \) of \( \ell_1 \) constraint violation and
!    objective value, and must not be dominated by the filter (see
!    [[sqpopt_filter_module]]). See [[filter_line_search]].
!    It is the default: on the Hock-Schittkowski test set (see
!    `test/test_hs_suite.f90`) it solves the most problems, with the
!    fewest function evaluations, of all the line search / merit function /
!    penalty combinations (and it has no penalty parameter to tune).
!  * `sqpopt_linesearch_funnel` -- the **funnel** method (Kiessling, Leyffer
!    & Vanaret, *"A unified funnel restoration SQP algorithm"*, Math.
!    Program. (2025); as implemented in the Uno solver): like the filter, no
!    merit function or penalty parameter, but the filter's list of points is
!    replaced by a single number, the funnel width, a bound on the
!    \( \ell_1 \) constraint violation that shrinks as the iterations
!    progress. A trial point must be inside the funnel; the same switching
!    condition then decides whether it must reduce the objective (an
!    Armijo test, "f-type") or sufficiently reduce the violation relative
!    to the funnel ("h-type", which shrinks the funnel). See
!    [[funnel_line_search]].
!
!  Two options from NLPQLP (Schittkowski) apply to the backtracking
!  searches: `interpolate` (on by default) chooses each backtracking step
!  length by safeguarded quadratic interpolation (see
!  [[next_step_length]]) rather than a fixed factor; and `nonmonotone_len`
!  (off by default) retries a failed search non-monotonically, against the
!  worst of the recent iterates (not used by the funnel search). On the Hock-Schittkowski test set the
!  interpolation solves one more problem with the filter search, with
!  fewer function evaluations (and a third fewer with the \( \ell_1 \)
!  merit function); the non-monotone retry doesn't help the filter search,
!  but does help the merit-function ones (e.g. with the augmented
!  Lagrangian and interpolation, 273 solved instead of 270).
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
    use sqpopt_types_module,   only: sqpopt_success, sqpopt_line_search_failed, sqpopt_sparse_matrix, sqpopt_all_finite, &
                                     l1_violation, merit_slack
    use sqpopt_merit_module,   only: sqpopt_merit_type, sqpopt_merit_l1, sqpopt_merit_augmented_lagrangian, &
                                     sqpopt_penalty_multipliers, sqpopt_penalty_model
    use sqpopt_filter_module,  only: sqpopt_filter_type
    use sqpopt_funnel_module,  only: sqpopt_funnel_type
    use fmin_module,           only: fmin
    use sqpopt_log_module,     only: sqpopt_log_type, sqpopt_log_detail, fmt_e, fmt_g

    implicit none

    private

    ! (re-exported, so that users and the other modules can get them from here:)
    public :: l1_violation
    public :: sqpopt_merit_l1, sqpopt_merit_augmented_lagrangian, sqpopt_penalty_multipliers, sqpopt_penalty_model

    integer, parameter, public :: sqpopt_linesearch_armijo   = 1  !! backtracking Armijo-type line search on a merit function
    integer, parameter, public :: sqpopt_linesearch_exact    = 2  !! (approximate) exact 1-D minimization of the merit function, via [[fmin]]
    integer, parameter, public :: sqpopt_linesearch_watchdog = 3  !! Powell's watchdog technique (relaxed acceptance + backtracking, see module docs)
    integer, parameter, public :: sqpopt_linesearch_filter   = 4  !! (default) Fletcher & Leyffer's filter method (no merit
                                                                  !! function/penalty parameter, see module docs)
    integer, parameter, public :: sqpopt_linesearch_funnel   = 5  !! the funnel method of Kiessling, Leyffer & Vanaret (no
                                                                  !! merit function/penalty parameter, see module docs)

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

        integer  :: mode        = sqpopt_linesearch_filter !! line search strategy to use (set from `options%linesearch_mode`)
        real(wp) :: tol         = 1.0e-4_wp !! desired tolerance on the minimizer (`sqpopt_linesearch_exact` mode)
        real(wp) :: sigma       = 0.1_wp    !! Armijo sufficient-decrease parameter, \( 0 < \sigma < 1 \) (`sqpopt_linesearch_armijo` mode)
        real(wp) :: backtrack   = 0.5_wp    !! step-length reduction factor at each backtracking step (`sqpopt_linesearch_armijo` mode)
        real(wp) :: alpha_min   = 1.0e-10_wp !! minimum step length (`armijo`/`watchdog`/`filter` modes): the search
                                             !! fails (no step is taken, `istat=sqpopt_line_search_failed`) if no
                                             !! acceptable step is found before `alpha` would drop below this
        integer  :: max_ls_iter = 40        !! maximum number of trial step lengths per search (`armijo`/`watchdog`/`filter` modes)
        logical  :: interpolate = .true.    !! choose each backtracking step length by safeguarded quadratic
                                            !! interpolation (as in NLPQLP) instead of the fixed factor `backtrack`
                                            !! (see [[next_step_length]]; `armijo`/`watchdog`/`filter` modes)
        integer  :: nonmonotone_len = 0     !! if `>0`, when a search fails it is retried non-monotonically
                                            !! (as in NLPQLP): against the worst merit value (`filter` mode: the
                                            !! worst violation and objective) of the last `nonmonotone_len`
                                            !! iterates instead of the current one (`armijo`/`filter` modes)
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

        ! the acceptance tests (options and state; see their modules):
        type(sqpopt_merit_type)  :: merit  !! the merit function and its penalty parameter (`armijo`, `exact`,
                                           !! and `watchdog` modes, and the trust region's merit-ratio test)
        type(sqpopt_filter_type) :: filter !! the filter (`filter` mode, and the trust region with the filter)
        type(sqpopt_funnel_type) :: funnel !! the funnel (`funnel` mode, and the trust region with the funnel)

        ! internal state for the non-monotone fallback (not user options): the merit values (`armijo`), or
        ! the violations and objective values (`filter`), at the most recent iterates
        integer  :: nm_count = 0
        real(wp), dimension(:), allocatable :: nm_phi, nm_theta, nm_f

        ! internal state for `sqpopt_linesearch_watchdog` mode (not user options -- persists across major iterations):
        logical  :: watchdog_ready               = .false. !! whether the best-point tracking below has been initialized
        integer  :: watchdog_relaxed_remaining   = 0        !! iterations left in the current relaxed window
        integer  :: watchdog_cooldown_remaining  = 0        !! iterations left before relaxed acceptance may reactivate
        real(wp) :: watchdog_w_opt               = 0.0_wp   !! best merit value found so far
        real(wp), dimension(:), allocatable :: watchdog_x_opt !! best point found so far `dimension(n)`

        ! the detailed log (set by `solve`), and what the last search used (outputs, for the log):
        type(sqpopt_log_type) :: log
        logical :: used_soc         = .false. !! the accepted step was a second-order-corrected one
        logical :: used_nonmonotone = .false. !! the accepted step came from the non-monotone retry
        logical :: used_relaxed     = .false. !! the step was a watchdog relaxed step

        contains

        procedure, public :: search                   => line_search
        procedure, public :: globalization_acceptable => point_acceptable

    end type sqpopt_linesearch_type

    contains
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

    me%used_soc         = .false.
    me%used_nonmonotone = .false.
    me%used_relaxed     = .false.

    select case (me%mode)
    case (sqpopt_linesearch_exact)
        call exact_line_search(me, eval_f, eval_c, x, p, f, c, lambda, c_lb, c_ub, alpha, x_new, istat, soc)
    case (sqpopt_linesearch_watchdog)
        call watchdog_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat, soc)
    case (sqpopt_linesearch_filter)
        call filter_line_search(me, eval_f, eval_c, x, p, f, g, c, c_lb, c_ub, alpha, x_new, istat, soc)
    case (sqpopt_linesearch_funnel)
        call funnel_line_search(me, eval_f, eval_c, x, p, f, g, c, c_lb, c_ub, alpha, x_new, istat, soc)
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

    subroutine eval_trial(me, eval_f, eval_c, x_trial, c_lb, c_ub, lambda, c_trial, phi_trial, ok, alpha)

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
    real(wp), optional,     intent(in)  :: alpha     !! step length (for the joint step, see [[eval_merit_function]])

    real(wp) :: f_trial

    call eval_fc(eval_f, eval_c, x_trial, f_trial, c_trial, ok)
    if (ok) then
        if (present(alpha)) then
            call me%merit%eval(f_trial, c_trial, c_lb, c_ub, lambda, phi_trial, alpha)
        else
            call me%merit%eval(f_trial, c_trial, c_lb, c_ub, lambda, phi_trial, 1.0_wp)
        end if
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
                                alpha, x_new, phi_new, accepted, soc, phi_ref)

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
    real(wp), optional,     intent(in)  :: phi_ref  !! reference value for the sufficient-decrease test, instead of
                                                    !! `phi0` (the non-monotone retry, see `nonmonotone_len`)

    real(wp), dimension(size(x)) :: x_trial, p_soc
    real(wp), dimension(size(c)) :: c_trial
    real(wp) :: phi_trial, slope, slack, phi_r
    logical  :: ok, soc_ok
    integer  :: it

    slope    = min(dphi0, 0.0_wp)
    phi_r    = phi0
    if (present(phi_ref)) phi_r = max(phi0, phi_ref)
    slack    = merit_slack(phi_r)
    accepted = .false.
    alpha    = alpha0

    do it = 1, me%max_ls_iter

        x_trial = x + alpha*p
        call eval_trial(me, eval_f, eval_c, x_trial, c_lb, c_ub, lambda, c_trial, phi_trial, ok, alpha)
        call log_merit_trial(me, alpha, phi_trial, phi_r + me%sigma*alpha*slope + slack, ok, .false.)
        if (ok .and. phi_trial <= phi_r + me%sigma*alpha*slope + slack) then
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
                    call eval_trial(me, eval_f, eval_c, x_trial, c_lb, c_ub, lambda, c_trial, phi_trial, ok, alpha)
                    call log_merit_trial(me, alpha, phi_trial, phi_r + me%sigma*alpha*slope + slack, ok, .true.)
                    if (ok .and. phi_trial <= phi_r + me%sigma*alpha*slope + slack) then
                        accepted = .true.
                        me%used_soc = .true.
                        exit
                    end if
                end if
            end if
        end if

        if (ok) then
            alpha = next_step_length(me, alpha, phi_r, slope, phi_trial)  ! (through the reference value tested against)
        else
            alpha = me%backtrack*alpha
        end if
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

    call me%merit%eval(f, c, c_lb, c_ub, lambda, phi0)
    call me%merit%directional_derivative(jac, g, p, c, c_lb, c_ub, lambda, dphi0)

    call backtrack_search(me, eval_f, eval_c, x, p, initial_step_length(x, p, me%major_step_limit), &
                          phi0, dphi0, c, c_lb, c_ub, lambda, alpha, x_new, phi_new, accepted, soc)

    if (.not. accepted .and. me%nonmonotone_len > 0 .and. me%nm_count > 0) then
        ! the non-monotone retry: against the worst merit value of the recent iterates
        if (maxval(me%nm_phi(1:me%nm_count)) > phi0) then
            call me%log%put(sqpopt_log_detail, 'ls  non-monotone retry, against merit '// &
                            fmt_g(maxval(me%nm_phi(1:me%nm_count))))
            call backtrack_search(me, eval_f, eval_c, x, p, initial_step_length(x, p, me%major_step_limit), &
                                  phi0, dphi0, c, c_lb, c_ub, lambda, alpha, x_new, phi_new, accepted, soc, &
                                  phi_ref=maxval(me%nm_phi(1:me%nm_count)))
            me%used_nonmonotone = accepted
        end if
    end if
    call nonmonotone_push(me, phi0, 0.0_wp, f)

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

    call me%merit%eval(f, c, c_lb, c_ub, lambda, phi0)
    alpha0 = initial_step_length(x, p, me%major_step_limit)
    alpha  = fmin(merit_along_direction, 0.0_wp, alpha0, me%tol)
    phi    = merit_along_direction(alpha)
    x_new  = x + alpha*p

    if (present(soc) .and. alpha < alpha0*(1.0_wp - sqrt(me%tol))) then
        call eval_trial(me, eval_f, eval_c, x + alpha0*p, c_lb, c_ub, lambda, c_full, phi_full, ok, alpha0)
        if (ok) then
            if (l1_violation(c_full, c_lb, c_ub) >= l1_violation(c, c_lb, c_ub)) then
                call soc(alpha0*p, c_full, p_soc, soc_ok)
                if (soc_ok) then
                    call eval_trial(me, eval_f, eval_c, x + p_soc, c_lb, c_ub, lambda, c_soc, phi_soc, ok, alpha0)
                    call me%log%put(sqpopt_log_detail, 'ls  second-order correction: merit '//fmt_g(phi_soc)// &
                                    merge(' (better, taken)   ', ' (not better)      ', ok .and. phi_soc < phi))
                    if (ok .and. phi_soc < phi) then
                        phi   = phi_soc
                        alpha = alpha0
                        x_new = x + p_soc
                        me%used_soc = .true.
                    end if
                end if
            end if
        end if
    end if

    call me%log%put(sqpopt_log_detail, 'ls  exact search: alpha '//fmt_e(alpha)//', merit '//fmt_g(phi)// &
                    ' (from '//fmt_g(phi0)//')')
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

        call eval_trial(me, eval_f, eval_c, x + alpha*p, c_lb, c_ub, lambda, c_trial, phi, ok, alpha)

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

    call me%merit%eval(f, c, c_lb, c_ub, lambda, phi0)
    call me%merit%directional_derivative(jac, g, p, c, c_lb, c_ub, lambda, dphi0)

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
        call eval_trial(me, eval_f, eval_c, x_new, c_lb, c_ub, lambda, c_trial, phi_trial, ok, alpha0)
        relaxed_used = ok
        me%used_relaxed = ok
        if (ok) call me%log%put(sqpopt_log_detail, 'ls  watchdog: relaxed step taken, merit '//fmt_g(phi_trial))
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
    logical :: ok, soc_ok, f_type, acc
    integer :: it

    theta0 = l1_violation(c, c_lb, c_ub)
    call me%filter%prepare(theta0)
    gtp = dot_product(g, p)
    alpha_lim = max(me%alpha_min, me%filter%min_step(theta0, gtp))

    alpha = initial_step_length(x, p, me%major_step_limit)
    do it = 1, me%max_ls_iter

        x_trial = x + alpha*p
        call eval_fc(eval_f, eval_c, x_trial, f_trial, c_trial, ok)
        if (ok) then
            theta_t = l1_violation(c_trial, c_lb, c_ub)
            acc = me%filter%accept(theta0, f, gtp, alpha, theta_t, f_trial, f_type)
            call log_theta_trial(me, alpha, theta_t, f_trial, .true., acc, f_type, .false.)
            if (acc) then
                call accept()
                return
            end if
        else
            call log_theta_trial(me, alpha, 0.0_wp, 0.0_wp, .false., .false., .false., .false.)
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
                        acc = me%filter%accept(theta0, f, gtp, alpha, theta_t, f_trial, f_type)
                        call log_theta_trial(me, alpha, theta_t, f_trial, .true., acc, f_type, .true.)
                        if (acc) then
                            me%used_soc = .true.
                            call accept()
                            return
                        end if
                    end if
                end if
            end if
        end if

        if (ok .and. me%interpolate) then
            ! interpolate the violation if the trial made it worse (the usual reason
            ! for rejection; its slope along a QP step is -theta0), else the objective
            if (theta_t > theta0 .and. theta0 > 0.0_wp) then
                alpha = next_step_length(me, alpha, theta0, -theta0, theta_t)
            else
                alpha = next_step_length(me, alpha, f, min(gtp, 0.0_wp), f_trial)
            end if
        else
            alpha = me%backtrack*alpha
        end if
        if (alpha < alpha_lim) exit

    end do

    if (me%nonmonotone_len > 0 .and. me%nm_count > 0) then
        ! the non-monotone retry: accept a point no worse, in both the violation
        ! and the objective (with an Armijo margin), than the worst of the
        ! recent iterates (and within theta_max)
        block
            real(wp) :: theta_ref, f_ref
            theta_ref = max(theta0, maxval(me%nm_theta(1:me%nm_count)))
            f_ref     = max(f, maxval(me%nm_f(1:me%nm_count)))
            if (theta_ref > theta0 .or. f_ref > f) then
                call me%log%put(sqpopt_log_detail, 'ls  non-monotone retry, against violation '//fmt_e(theta_ref)// &
                                ' and objective '//fmt_g(f_ref))
                alpha = initial_step_length(x, p, me%major_step_limit)
                do it = 1, me%max_ls_iter
                    x_trial = x + alpha*p
                    call eval_fc(eval_f, eval_c, x_trial, f_trial, c_trial, ok)
                    if (ok) then
                        theta_t = l1_violation(c_trial, c_lb, c_ub)
                        if (theta_t <= theta_ref .and. theta_t <= me%filter%theta_max .and. &
                            f_trial <= f_ref + me%filter%eta_phi*alpha*min(gtp, 0.0_wp)) then
                            call nonmonotone_push(me, 0.0_wp, theta0, f)
                            call me%filter%record(theta0, f)
                            me%used_nonmonotone = .true.
                            x_new = x_trial
                            istat = sqpopt_success
                            return
                        end if
                    end if
                    alpha = me%backtrack*alpha
                    if (alpha < me%alpha_min) exit
                end do
            end if
        end block
    end if
    call nonmonotone_push(me, 0.0_wp, theta0, f)

    ! no acceptable point was found: no step is taken
    call me%log%put(sqpopt_log_detail, 'ls  no acceptable step length (the next one would be below '// &
                    fmt_e(alpha_lim)//')')
    alpha = 0.0_wp
    x_new = x
    istat = sqpopt_line_search_failed

    contains

        subroutine accept()
        !! accept `x_trial`, augmenting the filter unless it was an f-type step
        call nonmonotone_push(me, 0.0_wp, theta0, f)
        if (.not. f_type) call me%filter%record(theta0, f)
        x_new = x_trial
        istat = sqpopt_success
        end subroutine accept

    end subroutine filter_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  the next trial step length after `alpha` was rejected: `backtrack*alpha`,
!  or, if `interpolate`, the minimizer of the quadratic interpolating
!  \( \phi(0) \), \( \phi'(0) \), and \( \phi(\alpha) \) (as in NLPQLP),
!  $$ \bar\alpha = \frac{-\tfrac12\alpha^2\phi'(0)}{\phi(\alpha)-\phi(0)-\alpha\phi'(0)} $$
!  safeguarded to \( [0.1\alpha, 0.5\alpha] \) (and `backtrack*alpha` if the
!  quadratic has no minimizer, e.g. \( \phi'(0) \ge 0 \)).

    pure function next_step_length(me, alpha, phi0, dphi0, phi_a) result(alpha_new)

    class(sqpopt_linesearch_type), intent(in) :: me
    real(wp), intent(in) :: alpha  !! the rejected step length
    real(wp), intent(in) :: phi0   !! \( \phi(0) \)
    real(wp), intent(in) :: dphi0  !! \( \phi'(0) \)
    real(wp), intent(in) :: phi_a  !! \( \phi(\alpha) \)
    real(wp) :: alpha_new

    real(wp) :: curv

    alpha_new = me%backtrack*alpha
    if (.not. me%interpolate .or. dphi0 >= 0.0_wp) return
    curv = phi_a - phi0 - alpha*dphi0
    if (curv <= 0.0_wp) return
    alpha_new = min(max(-0.5_wp*alpha**2*dphi0/curv, 0.1_wp*alpha), 0.5_wp*alpha)

    end function next_step_length
!*******************************************************************************

!*******************************************************************************
!>
!  add the current iterate's merit value (`armijo` mode), or violation and
!  objective value (`filter` mode), to the queue used by the non-monotone
!  retry (keeping the most recent `nonmonotone_len`).

    subroutine nonmonotone_push(me, phi, theta, f)

    class(sqpopt_linesearch_type), intent(inout) :: me
    real(wp), intent(in) :: phi, theta, f

    integer :: n

    n = me%nonmonotone_len
    if (n <= 0) return
    if (.not. allocated(me%nm_phi)) then
        allocate(me%nm_phi(n), me%nm_theta(n), me%nm_f(n))
        me%nm_count = 0
    end if
    if (me%nm_count == n) then
        me%nm_phi(1:n-1)   = me%nm_phi(2:n)
        me%nm_theta(1:n-1) = me%nm_theta(2:n)
        me%nm_f(1:n-1)     = me%nm_f(2:n)
    else
        me%nm_count = me%nm_count + 1
    end if
    me%nm_phi(me%nm_count)   = phi
    me%nm_theta(me%nm_count) = theta
    me%nm_f(me%nm_count)     = f

    end subroutine nonmonotone_push
!*******************************************************************************

!*******************************************************************************
!>
!  the funnel line search (see the module-level documentation), with the
!  rules of Kiessling, Leyffer & Vanaret, as implemented in the Uno solver.
!  With \( \theta \) the \( \ell_1 \) constraint violation, \( \varphi=f \),
!  and \( \tau \) the funnel width, a trial point \( x+\alpha p \) is
!  accepted if it is inside the funnel, \( \theta_t \le \tau \), and either
!
!  * (f-type step) the *switching condition*
!    \( \alpha(-g^Tp) > \delta\theta_k^{s_\theta} \) holds, and the Armijo
!    condition \( \varphi(x+\alpha p) \le \varphi_k + \eta\alpha g^Tp \)
!    holds; or
!  * (h-type step) otherwise, the violation is sufficiently inside the
!    funnel: \( \theta_t \le \beta\tau \). The funnel then shrinks (see
!    [[funnel_shrink]]).
!
!  (With `funnel%require_current`, the trial point must also be acceptable
!  with respect to the current point, see [[funnel_step_acceptable]].)
!  `alpha` is backtracked until acceptance, or until it falls below
!  `alpha_min`, in which case the search fails
!  (`istat=sqpopt_line_search_failed`, `x_new=x`); [[sqpopt_iterate_module]]
!  then shrinks the funnel toward the current violation (so the iterations
!  can't cycle back to it, see [[funnel_shrink_restoration]]) and takes a
!  feasibility restoration step. If the first trial is rejected without
!  reducing \( \theta \), its second-order correction is also tried (if
!  `soc` is present), with the same acceptance test. The non-monotone retry
!  (`nonmonotone_len`) is not used.

    subroutine funnel_line_search(me, eval_f, eval_c, x, p, f, g, c, c_lb, c_ub, alpha, x_new, istat, soc)

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
    real(wp) :: f_trial, theta0, theta_t, gtp
    logical :: ok, soc_ok, f_type, acc
    integer :: it

    theta0 = l1_violation(c, c_lb, c_ub)
    call me%funnel%prepare(theta0)
    gtp = dot_product(g, p)

    alpha = initial_step_length(x, p, me%major_step_limit)
    do it = 1, me%max_ls_iter

        x_trial = x + alpha*p
        call eval_fc(eval_f, eval_c, x_trial, f_trial, c_trial, ok)
        if (ok) then
            theta_t = l1_violation(c_trial, c_lb, c_ub)
            acc = me%funnel%accept(theta0, f, -alpha*gtp, theta_t, f_trial, f_type)
            call log_theta_trial(me, alpha, theta_t, f_trial, .true., acc, f_type, .false.)
            if (acc) then
                call accept()
                return
            end if
        else
            call log_theta_trial(me, alpha, 0.0_wp, 0.0_wp, .false., .false., .false., .false.)
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
                        acc = me%funnel%accept(theta0, f, -alpha*gtp, theta_t, f_trial, f_type)
                        call log_theta_trial(me, alpha, theta_t, f_trial, .true., acc, f_type, .true.)
                        if (acc) then
                            me%used_soc = .true.
                            call accept()
                            return
                        end if
                    end if
                end if
            end if
        end if

        if (ok .and. me%interpolate) then
            ! interpolate the violation if the trial made it worse, else the
            ! objective (as in [[filter_line_search]])
            if (theta_t > theta0 .and. theta0 > 0.0_wp) then
                alpha = next_step_length(me, alpha, theta0, -theta0, theta_t)
            else
                alpha = next_step_length(me, alpha, f, min(gtp, 0.0_wp), f_trial)
            end if
        else
            alpha = me%backtrack*alpha
        end if
        if (alpha < me%alpha_min) exit

    end do

    ! no acceptable point was found: no step is taken
    call me%log%put(sqpopt_log_detail, 'ls  no acceptable step length (the next one would be below '// &
                    fmt_e(me%alpha_min)//')')
    alpha = 0.0_wp
    x_new = x
    istat = sqpopt_line_search_failed

    contains

        subroutine accept()
        !! accept `x_trial`, shrinking the funnel after an h-type step
        if (.not. f_type) call me%funnel%record(theta0, theta_t)
        x_new = x_trial
        istat = sqpopt_success
        end subroutine accept

    end subroutine funnel_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  the detailed log's line for a trial point of a merit-function search:
!  its step length, merit value, and the value it had to reach.

    subroutine log_merit_trial(me, alpha, phi, target, ok, is_soc)

    class(sqpopt_linesearch_type), intent(in) :: me
    real(wp), intent(in) :: alpha, phi, target
    logical,  intent(in) :: ok     !! whether the functions were finite there
    logical,  intent(in) :: is_soc !! whether it was the second-order-corrected step

    character(len=:), allocatable :: what

    if (.not. me%log%on(sqpopt_log_detail)) return
    what = merge('ls  SOC   alpha ', 'ls  trial alpha ', is_soc)
    if (.not. ok) then
        call me%log%put(sqpopt_log_detail, what//fmt_e(alpha)//': non-finite function value, rejected')
    else
        call me%log%put(sqpopt_log_detail, what//fmt_e(alpha)//': merit '//fmt_g(phi)//', needed <= '// &
                        fmt_g(target)//merge(', accepted', ', rejected', phi <= target))
    end if

    end subroutine log_merit_trial
!*******************************************************************************

!*******************************************************************************
!>
!  the detailed log's line for a trial point of a filter or funnel search:
!  its step length, violation, objective, and whether (and as which type of
!  step) it was accepted.

    subroutine log_theta_trial(me, alpha, theta, f, ok, accepted, f_type, is_soc)

    class(sqpopt_linesearch_type), intent(in) :: me
    real(wp), intent(in) :: alpha, theta, f
    logical,  intent(in) :: ok       !! whether the functions were finite there
    logical,  intent(in) :: accepted
    logical,  intent(in) :: f_type   !! whether the switching condition held (an f-type step)
    logical,  intent(in) :: is_soc   !! whether it was the second-order-corrected step

    character(len=:), allocatable :: what, verdict

    if (.not. me%log%on(sqpopt_log_detail)) return
    what = merge('ls  SOC   alpha ', 'ls  trial alpha ', is_soc)
    if (.not. ok) then
        call me%log%put(sqpopt_log_detail, what//fmt_e(alpha)//': non-finite function value, rejected')
        return
    end if
    if (accepted) then
        if (f_type) then
            verdict = 'accepted (f-type: objective decrease)'
        else
            verdict = 'accepted (h-type: violation/objective)'
        end if
    else if (me%mode == sqpopt_linesearch_funnel) then
        verdict = 'rejected (funnel width '//fmt_e(me%funnel%width)//')'
    else
        verdict = 'rejected (filter: '//trim(adjustl(itoa(size(me%filter%theta))))//' entries)'
    end if
    call me%log%put(sqpopt_log_detail, what//fmt_e(alpha)//': violation '//fmt_e(theta)//', objective '// &
                    fmt_g(f)//', '//verdict)

    contains

        pure function itoa(i) result(str)
        integer, intent(in) :: i
        character(len=16) :: str
        write(str, '(I0)') i
        end function itoa

    end subroutine log_theta_trial
!*******************************************************************************

!*******************************************************************************
!>
!  whether a point with violation `theta` and objective `phi` is acceptable
!  to the filter (`sqpopt_linesearch_filter` mode) or to the funnel
!  (`sqpopt_linesearch_funnel` mode); always true in the other modes. Used
!  to decide when a feasibility restoration phase can end (see
!  [[sqpopt_restoration_module]]).

    function point_acceptable(me, theta, phi) result(ok)

    class(sqpopt_linesearch_type), intent(in) :: me
    real(wp), intent(in) :: theta !! constraint violation
    real(wp), intent(in) :: phi   !! objective
    logical :: ok

    select case (me%mode)
    case (sqpopt_linesearch_filter)
        ok = me%filter%acceptable(theta, phi)
    case (sqpopt_linesearch_funnel)
        ok = me%funnel%acceptable(theta)
    case default
        ok = .true.
    end select

    end function point_acceptable
!*******************************************************************************

    end module sqpopt_linesearch_module
!*******************************************************************************
