!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Merit function evaluation and line search used to globalize the SQP
!  iterations (ensures progress towards both optimality and feasibility).
!  Two line search strategies are available (`sqpopt_linesearch_type%mode`):
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

    module sqpopt_linesearch_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_success, sqpopt_line_search_failed
    use sqpopt_problem_module, only: sqpopt_objective_func, sqpopt_constraint_func
    use fmin_module,           only: fmin

    implicit none

    private

    integer, parameter, public :: sqpopt_linesearch_armijo = 1  !! backtracking Armijo-type line search (default)
    integer, parameter, public :: sqpopt_linesearch_exact  = 2  !! (approximate) exact 1-D minimization of the merit function, via [[fmin]]

    type, public :: sqpopt_linesearch_type
        !! options and state for the merit function and line search.

        integer  :: mode        = sqpopt_linesearch_armijo !! line search strategy to use
        real(wp) :: penalty     = 1.0_wp    !! current penalty parameter used in the merit function
        real(wp) :: tol         = 1.0e-4_wp !! desired tolerance on the minimizer (`sqpopt_linesearch_exact` mode)
        real(wp) :: sigma       = 0.1_wp    !! Armijo sufficient-decrease parameter, \( 0 < \sigma < 1 \) (`sqpopt_linesearch_armijo` mode)
        real(wp) :: backtrack   = 0.5_wp    !! step-length reduction factor at each backtracking step (`sqpopt_linesearch_armijo` mode)
        real(wp) :: alpha_min   = 0.1_wp    !! minimum step length, \( 0 < \alpha_{min} < 1 \) (`sqpopt_linesearch_armijo` mode):
                                            !! backtracking never goes below this; if the floor is reached without
                                            !! satisfying the sufficient-decrease test, `alpha_min` is accepted anyway
                                            !! (this avoids ever taking a useless near-zero step)
        integer  :: max_ls_iter = 20        !! maximum number of Armijo backtracking steps (`sqpopt_linesearch_armijo` mode)

        contains

        procedure, public :: eval_merit => eval_merit_function
        procedure, public :: search     => line_search

    end type sqpopt_linesearch_type

    ! module-level state used to pass context into the fmin callback below,
    ! since `fmin` requires a plain `function(alpha) result(phi)` interface
    ! with no way to pass extra context through its argument list:
    procedure(sqpopt_objective_func),  pointer, save, private :: ls_eval_f => null()
    procedure(sqpopt_constraint_func), pointer, save, private :: ls_eval_c => null()
    real(wp), dimension(:), allocatable, save, private :: ls_x, ls_p, ls_c_lb, ls_c_ub
    real(wp), save, private :: ls_penalty = 1.0_wp

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  evaluate the exact penalty (merit) function used to measure progress:
!
!  $$ \phi(x) = f(x) + \mu \lVert \max(c_l - c(x), 0, c(x) - c_u) \rVert_1 $$

    subroutine eval_merit_function(me, f, c, c_lb, c_ub, phi)

    class(sqpopt_linesearch_type), intent(inout) :: me
    real(wp),                intent(in)  :: f      !! objective function value
    real(wp), dimension(:), intent(in)  :: c      !! constraint values `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp),                intent(out) :: phi    !! value of the merit function

    phi = f + me%penalty*sum(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))

    end subroutine eval_merit_function
!*******************************************************************************

!*******************************************************************************
!>
!  perform a line search along the direction `p` to find an accepted step
!  length `alpha`, dispatching to the strategy selected by `me%mode`.

    subroutine line_search(me, eval_f, eval_c, x, p, f, g, c, c_lb, c_ub, alpha, istat)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: eval_f  !! evaluates \( f(x) \)
    procedure(sqpopt_constraint_func) :: eval_c  !! evaluates \( c(x) \)
    real(wp), dimension(:), intent(in)  :: x      !! current point `dimension(n)`
    real(wp), dimension(:), intent(in)  :: p      !! search direction `dimension(n)`
    real(wp),                intent(in)  :: f      !! objective function value at `x`
    real(wp), dimension(:), intent(in)  :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(in)  :: c      !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp),                intent(out) :: alpha  !! accepted step length
    integer,                  intent(out) :: istat  !! status code (see [[sqpopt_types_module]])

    select case (me%mode)
    case (sqpopt_linesearch_exact)
        call exact_line_search(me, eval_f, eval_c, x, p, c_lb, c_ub, alpha, istat)
    case default
        call armijo_line_search(me, eval_f, eval_c, x, p, f, g, c, c_lb, c_ub, alpha, istat)
    end select

    end subroutine line_search
!*******************************************************************************

!*******************************************************************************
!>
!  backtracking line search with an Armijo-type sufficient-decrease test
!  on the \( \ell_1 \) merit function (as used by default in `slsqp`):
!  starting from `alpha=1`, `alpha` is repeatedly reduced by `backtrack`
!  until
!  $$ \phi(x+\alpha p) \le \phi(x) + \sigma \alpha D(\phi;p) $$
!  where \( D(\phi;p) = g^T p - \mu \lVert \text{viol}(x) \rVert_1 \) is
!  the (approximate) directional derivative of the merit function along
!  `p`. `alpha` never goes below `alpha_min`: if the floor is reached
!  without satisfying the test (which can happen since `p` is only an
!  approximate QP solution, so is not guaranteed to be a descent direction
!  for `phi` in every case), `alpha_min` is accepted anyway -- this is much
!  cheaper than falling back to an exact line search, and still guarantees
!  the step never shrinks to a useless near-zero value.

    subroutine armijo_line_search(me, eval_f, eval_c, x, p, f, g, c, c_lb, c_ub, alpha, istat)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: eval_f
    procedure(sqpopt_constraint_func) :: eval_c
    real(wp), dimension(:), intent(in)  :: x, p, g, c, c_lb, c_ub
    real(wp),                intent(in)  :: f
    real(wp),                intent(out) :: alpha
    integer,                  intent(out) :: istat

    real(wp), dimension(size(x)) :: x_trial
    real(wp), dimension(size(c)) :: c_trial
    real(wp) :: phi0, dphi0, v0, phi_trial, f_trial
    integer :: it

    v0    = sum(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))
    phi0  = f + me%penalty*v0
    dphi0 = dot_product(g, p) - me%penalty*v0

    alpha = 1.0_wp
    do it = 1, me%max_ls_iter
        x_trial = x + alpha*p
        call eval_f(x_trial, f_trial)
        call eval_c(x_trial, c_trial)
        phi_trial = f_trial + me%penalty*sum(max(c_lb-c_trial, 0.0_wp) + max(c_trial-c_ub, 0.0_wp))
        if (phi_trial <= phi0 + me%sigma*alpha*dphi0) then
            istat = sqpopt_success
            return
        end if
        if (alpha <= me%alpha_min) exit
        alpha = max(me%backtrack*alpha, me%alpha_min)
    end do

    ! backtracking reached the floor without satisfying the Armijo test;
    ! accept `alpha_min` anyway rather than continuing to shrink toward zero:
    alpha = me%alpha_min
    istat = sqpopt_line_search_failed

    end subroutine armijo_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  (approximately) minimize the merit function along `p` using the
!  derivative-free 1-D minimizer [[fmin]].

    subroutine exact_line_search(me, eval_f, eval_c, x, p, c_lb, c_ub, alpha, istat)

    class(sqpopt_linesearch_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: eval_f  !! evaluates \( f(x) \)
    procedure(sqpopt_constraint_func) :: eval_c  !! evaluates \( c(x) \)
    real(wp), dimension(:), intent(in)  :: x      !! current point `dimension(n)`
    real(wp), dimension(:), intent(in)  :: p      !! search direction `dimension(n)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp),                intent(out) :: alpha  !! accepted step length
    integer,                  intent(out) :: istat  !! status code (see [[sqpopt_types_module]])

    ls_eval_f => eval_f
    ls_eval_c => eval_c
    ls_x      = x
    ls_p      = p
    ls_c_lb   = c_lb
    ls_c_ub   = c_ub
    ls_penalty = me%penalty

    alpha = fmin(merit_along_direction, 0.0_wp, 1.0_wp, me%tol)
    istat = sqpopt_success

    end subroutine exact_line_search
!*******************************************************************************

!*******************************************************************************
!>
!  the merit function \( \phi(x + \alpha p) \) along the search direction,
!  in the form required by [[fmin]]. Uses the module-level state set by
!  [[exact_line_search]] just before calling `fmin`.

    function merit_along_direction(alpha) result(phi)

    real(wp), intent(in) :: alpha
    real(wp) :: phi

    real(wp), dimension(size(ls_x)) :: x_trial, c_trial
    real(wp) :: f_trial

    x_trial = ls_x + alpha*ls_p
    call ls_eval_f(x_trial, f_trial)
    call ls_eval_c(x_trial, c_trial)
    phi = f_trial + ls_penalty*sum(max(ls_c_lb-c_trial, 0.0_wp) + max(c_trial-ls_c_ub, 0.0_wp))

    end function merit_along_direction
!*******************************************************************************

    end module sqpopt_linesearch_module
!*******************************************************************************
