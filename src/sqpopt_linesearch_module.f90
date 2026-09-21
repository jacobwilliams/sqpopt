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
!    (1982), as used in `references/vf13`): tracks the best point found so
!    far and, for a short window after a genuine improvement, *relaxes* the
!    sufficient-decrease test to allow the merit function to temporarily
!    get worse (accepting the full quasi-Newton step outright) rather than
!    stalling near a curved/simultaneously-active constraint boundary (the
!    Maratos effect). If the relaxed window is used up without a new best
!    point, it backtracks all the way to the best point found so far and
!    disables relaxed acceptance for `watchdog_cooldown_len` iterations.
!    This targets the same failure mode as the second-order correction and
!    the augmented Lagrangian merit function, via a different mechanism.


    module sqpopt_linesearch_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_success, sqpopt_line_search_failed, sqpopt_sparse_matrix
    use sqpopt_problem_module, only: sqpopt_objective_func, sqpopt_constraint_func
    use sqpopt_linalg_module,  only: sparse_matvec_transpose
    use fmin_module,           only: fmin

    implicit none

    private

    integer, parameter, public :: sqpopt_linesearch_armijo   = 1  !! backtracking Armijo-type line search (default)
    integer, parameter, public :: sqpopt_linesearch_exact    = 2  !! (approximate) exact 1-D minimization of the merit function, via [[fmin]]
    integer, parameter, public :: sqpopt_linesearch_watchdog = 3  !! Powell's watchdog technique (relaxed acceptance + backtracking, see module docs)

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

        integer  :: watchdog_relaxed_len    = 2     !! number of relaxed steps tolerated before requiring a new best point (`sqpopt_linesearch_watchdog` mode)
        integer  :: watchdog_cooldown_len   = 10    !! number of iterations relaxed acceptance is disabled for after a backtrack (`sqpopt_linesearch_watchdog` mode)

        ! internal state for `sqpopt_linesearch_watchdog` mode (not user options -- persists across major iterations):
        logical  :: watchdog_ready               = .false. !! whether the best-point tracking below has been initialized
        integer  :: watchdog_relaxed_remaining   = 0        !! iterations left in the current relaxed window
        integer  :: watchdog_cooldown_remaining  = 0        !! iterations left before relaxed acceptance may reactivate
        real(wp) :: watchdog_w_opt               = 0.0_wp   !! best merit value found so far
        real(wp), dimension(:), allocatable :: watchdog_x_opt !! best point found so far `dimension(n)`

        contains

        procedure, public :: eval_merit             => eval_merit_function
        procedure, public :: directional_derivative  => merit_directional_derivative
        procedure, public :: search                  => line_search

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
    real(wp), dimension(:), intent(in)  :: c, c_lb, c_ub, lambda
    real(wp), dimension(:), intent(out) :: s

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
    type(sqpopt_sparse_matrix), intent(in) :: jac     !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:), intent(in) :: g           !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(in) :: p           !! search direction `dimension(n)`
    real(wp), dimension(:), intent(in) :: c, c_lb, c_ub, lambda
    real(wp), intent(out) :: dphi0

    real(wp), dimension(size(c)) :: s
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

    subroutine line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat)

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
    case default
        call armijo_line_search(me, eval_f, eval_c, x, p, f, g, c, jac, lambda, c_lb, c_ub, alpha, x_new, istat)
    end select

    end subroutine line_search
!*******************************************************************************

!*******************************************************************************
!>
!  backtracking line search with an Armijo-type sufficient-decrease test
!  on the merit function (as used by default in `slsqp`):
!  starting from `alpha=1`, `alpha` is repeatedly reduced by `backtrack`
!  until
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
    real(wp), dimension(:), intent(in)  :: x, p, g, c, lambda, c_lb, c_ub
    type(sqpopt_sparse_matrix), intent(in) :: jac
    real(wp),                intent(in)  :: f
    real(wp),                intent(out) :: alpha
    real(wp), dimension(:), intent(out) :: x_new
    integer,                  intent(out) :: istat

    real(wp), dimension(size(x)) :: x_trial
    real(wp), dimension(size(c)) :: c_trial
    real(wp) :: phi0, dphi0, phi_trial, f_trial
    integer :: it

    call me%eval_merit(f, c, c_lb, c_ub, lambda, phi0)
    call me%directional_derivative(jac, g, p, c, c_lb, c_ub, lambda, dphi0)

    alpha = 1.0_wp
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

    alpha = fmin(merit_along_direction, 0.0_wp, 1.0_wp, me%tol)
    x_new = x + alpha*p
    istat = sqpopt_success

    contains

    !*******************************************************************************
    !>
    !  the merit function \( \phi(x + \alpha p) \) along the search direction,
    !  in the form required by [[fmin]]. Uses the host-associated variables
    !  set by [[exact_line_search]].

        function merit_along_direction(alpha) result(phi)

        real(wp), intent(in) :: alpha
        real(wp) :: phi

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
!  *Math. Prog. Study 16* (1982); see `references/vf13`): a variant of
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
    real(wp), dimension(:), intent(in)  :: x, p, g, c, lambda, c_lb, c_ub
    type(sqpopt_sparse_matrix), intent(in) :: jac
    real(wp),                intent(in)  :: f
    real(wp),                intent(out) :: alpha
    real(wp), dimension(:), intent(out) :: x_new
    integer,                  intent(out) :: istat

    real(wp), dimension(size(x)) :: x_trial
    real(wp), dimension(size(c)) :: c_trial
    real(wp) :: phi0, dphi0, phi_trial, f_trial
    logical :: standard_ok, relaxed_used
    integer :: it

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
    alpha = 1.0_wp
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
        ! the watchdog technique allows a relaxed *full* step here instead
        ! of stalling (the merit function may temporarily get worse):
        alpha = 1.0_wp
        x_trial = x + p
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

        if (standard_ok .and. alpha >= 0.99_wp) then
            ! a good, (nearly) full accepted step: "reward" the next few
            ! calls with a fresh window of relaxed acceptance (mirroring
            ! VF13's ISWDOG logic, which opens the relaxed window after an
            ! iteration whose merit reduction was large enough):
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

    end module sqpopt_linesearch_module
!*******************************************************************************
