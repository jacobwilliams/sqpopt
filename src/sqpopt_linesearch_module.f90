!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Merit function evaluation and line search used to globalize the SQP
!  iterations (ensures progress towards both optimality and feasibility).
!  The 1-D minimization of the merit function along the search direction
!  is done using the derivative-free [[fmin]] routine (from the `fmin`
!  dependency), rather than a hand-written backtracking search.

    module sqpopt_linesearch_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_success
    use sqpopt_problem_module, only: sqpopt_objective_func, sqpopt_constraint_func
    use fmin_module,           only: fmin

    implicit none

    private

    type, public :: sqpopt_linesearch_type
        !! options and state for the merit function and line search.

        real(wp) :: penalty = 1.0_wp  !! current penalty parameter used in the merit function
        real(wp) :: tol     = 1.0e-4_wp !! desired tolerance on the 1-D line search minimizer

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
!  perform a line search along the direction `p` to find a step length
!  `alpha` \( \in (0,1] \) that (approximately) minimizes the merit
!  function, using the derivative-free 1-D minimizer [[fmin]].

    subroutine line_search(me, eval_f, eval_c, x, p, c_lb, c_ub, alpha, istat)

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

    end subroutine line_search
!*******************************************************************************

!*******************************************************************************
!>
!  the merit function \( \phi(x + \alpha p) \) along the search direction,
!  in the form required by [[fmin]]. Uses the module-level state set by
!  [[line_search]] just before calling `fmin`.

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
