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
!    correction (see [[sqpopt_iterate_module]]) to avoid the Maratos effect.
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
!  * `sqpopt_linesearch_filter` -- a line-search adaptation of Fletcher &
!    Leyffer's **filter** method (*"Nonlinear programming without a penalty
!    function"*, Math. Program. 91 (2002), see `references/fletcher.pdf`):
!    dispenses with the merit function/`penalty` parameter entirely and
!    instead accepts a trial point `x+alpha*p` if the pair `(f,h)` of
!    objective value and \( \ell_1 \) constraint violation is not
!    dominated by any `(f,h)` pair from a previously-accepted iterate (the
!    "filter"), using the paper's sufficient-reduction envelope (their eqs.
!    3-4) to exclude points arbitrarily close to the filter. The original
!    paper's algorithm is trust-region-based, with an explicit feasibility
!    restoration phase for infeasible QPs and NW/SE filter corner rules;
!    this port instead backtracks `alpha` (like the other three modes)
!    until an acceptable point is found, and omits the restoration phase
!    and corner rules (the paper itself notes the corner rules can be
!    dispensed with -- §3.6). See [[filter_line_search]] for details.
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
    use sqpopt_types_module,   only: sqpopt_success, sqpopt_line_search_failed, sqpopt_sparse_matrix
    use sqpopt_problem_module, only: sqpopt_objective_func, sqpopt_constraint_func
    use sqpopt_linalg_module,  only: sparse_matvec_transpose
    use fmin_module,           only: fmin

    implicit none

    private

    public :: l1_violation, filter_penalty_estimate

    integer, parameter, public :: sqpopt_linesearch_armijo   = 1  !! backtracking Armijo-type line search (default)
    integer, parameter, public :: sqpopt_linesearch_exact    = 2  !! (approximate) exact 1-D minimization of the merit function, via [[fmin]]
    integer, parameter, public :: sqpopt_linesearch_watchdog = 3  !! Powell's watchdog technique (relaxed acceptance + backtracking, see module docs)
    integer, parameter, public :: sqpopt_linesearch_filter   = 4  !! Fletcher & Leyffer's filter method (no merit function/penalty parameter, see module docs)

    integer, parameter, public :: sqpopt_merit_l1                   = 1  !! non-smooth \( \ell_1 \) exact penalty merit function (default)
    integer, parameter, public :: sqpopt_merit_augmented_lagrangian = 2  !! smooth augmented Lagrangian merit function (NPSOL/SNOPT-style)

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
        real(wp) :: alpha_min   = 0.1_wp    !! minimum step length, \( 0 < \alpha_{min} < 1 \) (`sqpopt_linesearch_armijo` mode):
                                            !! backtracking never goes below this; if the floor is reached without
                                            !! satisfying the sufficient-decrease test, `alpha_min` is accepted anyway
                                            !! (this avoids ever taking a useless near-zero step)
        integer  :: max_ls_iter = 20        !! maximum number of Armijo backtracking steps (`sqpopt_linesearch_armijo` mode)
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

        real(wp) :: filter_beta     = 0.99_wp    !! envelope constant \( \beta \) in the filter's sufficient-reduction test, eq. (3) (`sqpopt_linesearch_filter` mode)
        real(wp) :: filter_alpha1   = 0.25_wp    !! envelope constant \( \alpha_1 \) (weight on the QP-predicted decrease `q`), eq. (4) (`sqpopt_linesearch_filter` mode)
        real(wp) :: filter_alpha2   = 1.0e-4_wp  !! envelope constant \( \alpha_2 \) (weight on \( h \mu \)), eq. (4) (`sqpopt_linesearch_filter` mode)
        real(wp) :: filter_ubd      = 100.0_wp   !! default floor on the upper bound `u` on the constraint violation, §3.2 (`sqpopt_linesearch_filter` mode)
        real(wp) :: filter_tt       = 1.25_wp    !! multiplier on the initial constraint violation used to set `u`, §3.2 (`sqpopt_linesearch_filter` mode)
        real(wp) :: filter_feas_tol = 1.0e-8_wp  !! below this constraint violation, a trial point is treated as "feasible": if both the
                                                  !! current point and the trial point are feasible, plain sufficient decrease in `f` is
                                                  !! also required (the filter test alone is vacuous when `h` stays at/near zero, e.g. for
                                                  !! problems with no nonlinear constraints) (`sqpopt_linesearch_filter` mode)

        ! internal state for `sqpopt_linesearch_filter` mode (not user options -- persists across major iterations):
        logical  :: filter_ready = .false.  !! whether the filter/upper-bound below has been initialized
        real(wp) :: filter_u     = 0.0_wp    !! current upper bound `u` on the constraint violation, §3.2
        real(wp), dimension(:), allocatable :: filter_f  !! objective values of the filter's `(f,h,q,mu)` entries
        real(wp), dimension(:), allocatable :: filter_h  !! constraint-violation values of the filter's entries
        real(wp), dimension(:), allocatable :: filter_q  !! QP-predicted objective decrease at each entry, used by eq. (4)
        real(wp), dimension(:), allocatable :: filter_mu !! penalty-parameter estimate at each entry, used by eq. (4)

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
        procedure, public :: filter_test             => filter_acceptable
        procedure, public :: filter_record           => filter_add

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
!  * `sqpopt_merit_l1`: \( D(\phi;p) = g^Tp - \mu \lVert \text{viol}(x) \rVert_1 \)
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
        dphi0 = dot_product(g, p) - me%penalty*sum(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))
    end select

    end subroutine merit_directional_derivative
!*******************************************************************************

!*******************************************************************************
!>
!  perform a line search along the direction `p` to find an accepted step
!  length `alpha`, dispatching to the strategy selected by `me%mode`.

    subroutine line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, q, alpha, x_new, istat)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: eval_f  !! evaluates \( f(x) \)
    procedure(sqpopt_constraint_func) :: eval_c  !! evaluates \( c(x) \)
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
    real(wp),                intent(in)  :: q      !! the QP's predicted decrease in `f` along `p`,
                                                    !! \( q = -\left(g^Tp + \tfrac{1}{2}p^THp\right) \)
                                                    !! (only used by `sqpopt_linesearch_filter`)
    real(wp),                intent(out) :: alpha  !! accepted step length
    real(wp), dimension(:), intent(out) :: x_new   !! the accepted new point `dimension(n)` (normally
                                                    !! `x + alpha*p`, except `sqpopt_linesearch_watchdog`
                                                    !! may instead return an earlier best point on backtrack)
    integer,                  intent(out) :: istat  !! status code (see [[sqpopt_types_module]])

    select case (me%mode)
    case (sqpopt_linesearch_exact)
        call exact_line_search(me, eval_f, eval_c, x, p, lambda, c_lb, c_ub, alpha, x_new, istat)
    case (sqpopt_linesearch_watchdog)
        call watchdog_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat)
    case (sqpopt_linesearch_filter)
        call filter_line_search(me, eval_f, eval_c, x, p, f, c, lambda, c_lb, c_ub, q, alpha, x_new, istat)
    case default
        call armijo_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat)
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

    if (rmax > 0.0_wp) then
        alpha0 = min(1.0_wp, step_limit/rmax)
    else
        alpha0 = 1.0_wp
    end if

    end function initial_step_length
!*******************************************************************************

!*******************************************************************************
!>
!  backtracking line search with an Armijo-type sufficient-decrease test
!  on the merit function (as used by default in `slsqp`):
!  starting from \( \alpha_0 \) (see [[initial_step_length]], normally `1`
!  unless capped by `major_step_limit`), `alpha` is repeatedly reduced by
!  `backtrack` until
!  $$ \phi(x+\alpha p) \le \phi(x) + \sigma \alpha D(\phi;p) $$
!  `alpha` never goes below `alpha_min`: if the floor is reached without
!  satisfying the test (which can happen since `p` is only an approximate
!  QP solution, so is not guaranteed to be a descent direction for `phi` in
!  every case), `alpha_min` is accepted anyway -- this is much cheaper than
!  falling back to an exact line search, and still guarantees the step
!  never shrinks to a useless near-zero value.

    subroutine armijo_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: eval_f
    procedure(sqpopt_constraint_func) :: eval_c
    real(wp), dimension(:), intent(in) :: x !! current point `x`
    real(wp), dimension(:), intent(in) :: p !! search direction `p`
    real(wp), dimension(:), intent(in) :: g !! gradient of the objective at `x`
    real(wp), dimension(:), intent(in) :: c !! constraint values at `x`
    real(wp), dimension(:), intent(in) :: lambda !! Lagrange multipliers at `x`
    real(wp), dimension(:), intent(in) :: c_lb   !! lower bounds on the constraints
    real(wp), dimension(:), intent(in) :: c_ub   !! upper bounds on the constraints
    type(sqpopt_sparse_matrix), intent(in) :: jac    !! constraint Jacobian at `x` (`dimension(m,n)`)
    real(wp),                   intent(in)  :: f     !! objective function value at `x`
    real(wp),                   intent(out) :: alpha !! step length along `p`
    real(wp), dimension(:),     intent(out) :: x_new !! new point `x + alpha*p`
    integer,                    intent(out) :: istat !! status of the line search (success or failure)

    real(wp), dimension(size(x)) :: x_trial !! trial point `x + alpha*p`
    real(wp), dimension(size(c)) :: c_trial !! constraint values at the trial point
    real(wp) :: phi0, dphi0, phi_trial, f_trial !! merit function values and directional derivative
    integer :: it !! iteration counter for the line search loop

    call me%eval_merit(f, c, c_lb, c_ub, lambda, phi0)
    call me%directional_derivative(jac, g, p, c, c_lb, c_ub, lambda, dphi0)

    alpha = initial_step_length(x, p, me%major_step_limit)
    do it = 1, me%max_ls_iter
        x_trial = x + alpha*p
        call eval_f(x_trial, f_trial)
        call eval_c(x_trial, c_trial)
        call me%eval_merit(f_trial, c_trial, c_lb, c_ub, lambda, phi_trial)
        if (phi_trial <= phi0 + me%sigma*alpha*dphi0) then
            istat = sqpopt_success
            x_new = x_trial
            return
        end if
        if (alpha <= me%alpha_min) exit
        alpha = max(me%backtrack*alpha, me%alpha_min)
    end do

    ! backtracking reached the floor without satisfying the Armijo test;
    ! accept `alpha_min` anyway rather than continuing to shrink toward zero:
    alpha = me%alpha_min
    x_new = x_trial
    istat = sqpopt_line_search_failed

    end subroutine armijo_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  (approximately) minimize the merit function along `p` using the
!  derivative-free 1-D minimizer [[fmin]].

    subroutine exact_line_search(me, eval_f, eval_c, x, p, lambda, c_lb, c_ub, alpha, x_new, istat)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: eval_f  !! evaluates \( f(x) \)
    procedure(sqpopt_constraint_func) :: eval_c  !! evaluates \( c(x) \)
    real(wp), dimension(:), intent(in)  :: x      !! current point `dimension(n)`
    real(wp), dimension(:), intent(in)  :: p      !! search direction `dimension(n)`
    real(wp), dimension(:), intent(in)  :: lambda !! Lagrange multiplier estimate `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp),                intent(out) :: alpha  !! accepted step length
    real(wp), dimension(:), intent(out) :: x_new   !! the accepted new point `dimension(n)`
    integer,                  intent(out) :: istat  !! status code (see [[sqpopt_types_module]])

    alpha = fmin(merit_along_direction, 0.0_wp, initial_step_length(x, p, me%major_step_limit), me%tol)
    x_new = x + alpha*p
    istat = sqpopt_success

    contains

    !*******************************************************************************
    !>
    !  the merit function \( \phi(x + \alpha p) \) along the search direction,
    !  in the form required by [[fmin]]. Uses the host-associated variables
    !  set by [[exact_line_search]].

        function merit_along_direction(alpha) result(phi)

        real(wp), intent(in) :: alpha !! step length along the search direction
        real(wp) :: phi !! merit function value at the trial point

        real(wp), dimension(size(x)) :: x_trial
        real(wp), dimension(size(c_lb)) :: c_trial
        real(wp) :: f_trial

        x_trial = x + alpha*p
        call eval_f(x_trial, f_trial)
        call eval_c(x_trial, c_trial)
        call me%eval_merit(f_trial, c_trial, c_lb, c_ub, lambda, phi)

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
!  decrease test -- rather than stalling near a curved or simultaneously-
!  active constraint boundary (the Maratos effect). If none of those
!  relaxed steps beats the best point found so far, the search backtracks
!  all the way to that best point and disables relaxed acceptance for
!  `watchdog_cooldown_len` further calls.

    subroutine watchdog_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: eval_f
    procedure(sqpopt_constraint_func) :: eval_c
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

    real(wp), dimension(size(x)) :: x_trial !! trial point during the line search
    real(wp), dimension(size(c)) :: c_trial !! constraint values at the trial point
    real(wp) :: phi0, dphi0, phi_trial, f_trial, alpha0 !! merit function values, directional derivative, trial objective, initial step length
    logical :: standard_ok, relaxed_used !! flags indicating if standard or relaxed line search succeeded
    integer :: it !! iteration counter for the line search loop

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
    alpha = alpha0
    standard_ok = .false.
    do it = 1, me%max_ls_iter
        x_trial = x + alpha*p
        call eval_f(x_trial, f_trial)
        call eval_c(x_trial, c_trial)
        call me%eval_merit(f_trial, c_trial, c_lb, c_ub, lambda, phi_trial)
        if (phi_trial <= phi0 + me%sigma*alpha*dphi0) then
            standard_ok = .true.
            exit
        end if
        if (alpha <= me%alpha_min) exit
        alpha = max(me%backtrack*alpha, me%alpha_min)
    end do

    relaxed_used = .false.
    if (.not. standard_ok .and. me%watchdog_relaxed_remaining > 0 .and. me%watchdog_cooldown_remaining == 0) then
        ! the standard sufficient-decrease test failed even at `alpha_min`;
        ! the watchdog technique allows a relaxed step here instead of
        ! stalling (the merit function may temporarily get worse), still
        ! capped by the major step limit like the initial standard-search step:
        alpha = alpha0
        x_trial = x + alpha*p
        call eval_f(x_trial, f_trial)
        call eval_c(x_trial, c_trial)
        call me%eval_merit(f_trial, c_trial, c_lb, c_ub, lambda, phi_trial)
        relaxed_used = .true.
    end if

    x_new = x_trial
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
        else if (.not. standard_ok) then
            ! plain Armijo floor reached, with no relaxed window available:
            istat = sqpopt_line_search_failed
        end if
    end block

    end subroutine watchdog_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  Fletcher & Leyffer's filter method (see the module-level documentation
!  and `references/fletcher.pdf`), adapted to a backtracking line search:
!  a trial point `x+alpha*p` is accepted if its `(f,h)` pair -- objective
!  value and \( \ell_1 \) constraint violation \( h = \lVert
!  \max(c_l-c,0,c-c_u) \rVert_1 \) -- is not dominated by any prior
!  accepted iterate's `(f,h)` pair (the "filter"), using the paper's
!  eqs. (3)-(4) sufficient-reduction envelope so that points arbitrarily
!  close to an existing filter entry are excluded (see [[filter_acceptable]]).
!  No merit function or penalty parameter is used. If no trial `alpha`
!  is accepted before `alpha_min`, that floor is accepted anyway (as in
!  [[armijo_line_search]]) and still recorded in the filter.

    subroutine filter_line_search(me, eval_f, eval_c, x, p, f, c, lambda, c_lb, c_ub, q, alpha, x_new, istat)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: eval_f
    procedure(sqpopt_constraint_func) :: eval_c
    real(wp), dimension(:), intent(in)  :: x, p, c, lambda, c_lb, c_ub
    real(wp),                intent(in)  :: f      !! objective function value at `x`
    real(wp),                intent(in)  :: q      !! the QP's predicted decrease in `f` along `p`
    real(wp),                intent(out) :: alpha
    real(wp), dimension(:), intent(out) :: x_new
    integer,                  intent(out) :: istat

    real(wp), dimension(size(x)) :: x_trial
    real(wp), dimension(size(c)) :: c_trial
    real(wp) :: f_trial, h_trial, h0, mu, alpha0
    logical :: both_feasible, accept
    integer :: it

    h0 = l1_violation(c, c_lb, c_ub)
    both_feasible = h0 <= me%filter_feas_tol

    call me%filter_prepare(h0)

    mu = filter_penalty_estimate(lambda)

    alpha0 = initial_step_length(x, p, me%major_step_limit)
    alpha = alpha0
    do it = 1, me%max_ls_iter
        x_trial = x + alpha*p
        call eval_f(x_trial, f_trial)
        call eval_c(x_trial, c_trial)
        h_trial = l1_violation(c_trial, c_lb, c_ub)

        accept = filter_acceptable(me, f_trial, h_trial)
        if (accept .and. both_feasible .and. h_trial <= me%filter_feas_tol) then
            ! both the current and trial points are essentially feasible, so
            ! h stays at/near zero and the filter test alone is vacuous (see
            ! the module docs) -- also require plain descent in `f`:
            accept = f_trial < f
        end if

        if (accept) then
            call filter_add(me, f_trial, h_trial, q, mu)
            x_new = x_trial
            istat = sqpopt_success
            return
        end if
        if (alpha <= me%alpha_min) exit
        alpha = max(me%backtrack*alpha, me%alpha_min)
    end do

    ! backtracking reached the floor without an acceptable point; accept
    ! `alpha_min` anyway (as in [[armijo_line_search]]) and still record it
    ! in the filter so later iterations don't keep proposing the same point:
    alpha = me%alpha_min
    x_new = x_trial
    call filter_add(me, f_trial, h_trial, q, mu)
    istat = sqpopt_line_search_failed

    end subroutine filter_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the filter (empty list) and its upper bound `filter_u`
!  (§3.2) the first time it is used, given the constraint violation `h0`
!  at the starting point; a no-op on every subsequent call. Shared by
!  [[filter_line_search]] and `sqpopt_trust_region_module`'s filter-based
!  acceptance test, so both start from the same, single filter.

    subroutine filter_prepare_state(me, h0)

    class(sqpopt_linesearch_type), intent(inout) :: me
    real(wp),                      intent(in)    :: h0

    if (me%filter_ready) return

    me%filter_u = max(me%filter_ubd, me%filter_tt*h0)
    if (allocated(me%filter_f)) deallocate(me%filter_f, me%filter_h, me%filter_q, me%filter_mu)
    allocate(me%filter_f(0), me%filter_h(0), me%filter_q(0), me%filter_mu(0))
    me%filter_ready = .true.

    end subroutine filter_prepare_state
!*******************************************************************************

!*******************************************************************************
!>
!  the \( \ell_1 \) constraint violation \( h(x) = \lVert \max(c_l-c,0,c-c_u)
!  \rVert_1 \) used by `sqpopt_linesearch_filter` (the same measure used,
!  with a `penalty` multiplier, by the `sqpopt_merit_l1` merit function).

    pure function l1_violation(c, c_lb, c_ub) result(h)

    real(wp), dimension(:), intent(in) :: c, c_lb, c_ub
    real(wp) :: h

    h = sum(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))

    end function l1_violation
!*******************************************************************************

!*******************************************************************************
!>
!  an estimate of the filter's penalty parameter \( \mu \) for the current
!  iterate (Fletcher & Leyffer, \u00a73.5): the least power of ten larger than
!  \( \lVert \lambda \rVert_\infty \), clipped to \( [10^{-6},10^6] \).

    pure function filter_penalty_estimate(lambda) result(mu)

    real(wp), dimension(:), intent(in) :: lambda
    real(wp) :: mu, lam_inf

    lam_inf = 0.0_wp
    if (size(lambda) > 0) lam_inf = maxval(abs(lambda))

    if (lam_inf <= 0.0_wp) then
        mu = 1.0e-6_wp
    else
        mu = min(max(10.0_wp**ceiling(log10(lam_inf)), 1.0e-6_wp), 1.0e6_wp)
    end if

    end function filter_penalty_estimate
!*******************************************************************************

!*******************************************************************************
!>
!  whether the pair `(f_trial,h_trial)` is acceptable to the current
!  filter: not exceeding the upper bound `filter_u` on the constraint
!  violation (\u00a73.2), and, for every existing filter entry `l`, satisfying
!  the sufficient-reduction envelope (eqs. 3-4):
!  $$ h_{trial} \le \beta h^{(l)} \quad \text{or} \quad
!  f_{trial} \le f^{(l)} - \max\!\left(\alpha_1 q^{(l)}, \alpha_2 h^{(l)}
!  \mu^{(l)}\right) $$

    function filter_acceptable(me, f_trial, h_trial) result(ok)

    class(sqpopt_linesearch_type), intent(in) :: me
    real(wp),                      intent(in) :: f_trial, h_trial
    logical :: ok

    integer :: l

    ok = h_trial <= me%filter_u
    if (.not. ok) return

    do l = 1, size(me%filter_f)
        if (.not. (h_trial <= me%filter_beta*me%filter_h(l) .or. &
                   f_trial <= me%filter_f(l) - max(me%filter_alpha1*me%filter_q(l), &
                                                    me%filter_alpha2*me%filter_h(l)*me%filter_mu(l)))) then
            ok = .false.
            return
        end if
    end do

    end function filter_acceptable
!*******************************************************************************

!*******************************************************************************
!>
!  add `(f_new,h_new,q_new,mu_new)` to the filter, first removing any
!  existing entries that it dominates (an entry `l` is dominated by the
!  new point if `f_new<=f^(l)` and `h_new<=h^(l)`).

    subroutine filter_add(me, f_new, h_new, q_new, mu_new)

    class(sqpopt_linesearch_type), intent(inout) :: me
    real(wp),                      intent(in)    :: f_new, h_new, q_new, mu_new

    logical, dimension(:), allocatable :: keep

    if (size(me%filter_f) > 0) then
        allocate(keep(size(me%filter_f)))
        keep = .not. (f_new <= me%filter_f .and. h_new <= me%filter_h)
        me%filter_f  = pack(me%filter_f,  keep)
        me%filter_h  = pack(me%filter_h,  keep)
        me%filter_q  = pack(me%filter_q,  keep)
        me%filter_mu = pack(me%filter_mu, keep)
    end if

    me%filter_f  = [me%filter_f,  f_new]
    me%filter_h  = [me%filter_h,  h_new]
    me%filter_q  = [me%filter_q,  q_new]
    me%filter_mu = [me%filter_mu, mu_new]

    end subroutine filter_add
!*******************************************************************************

    end module sqpopt_linesearch_module
!*******************************************************************************
