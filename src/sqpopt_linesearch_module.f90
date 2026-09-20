!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Merit function evaluation and line search used to globalize the SQP
!  iterations (ensures progress towards both optimality and feasibility).

    module sqpopt_linesearch_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    type, public :: sqpopt_linesearch_type
        !! options and state for the merit function and line search.

        real(wp) :: penalty = 0.0_wp  !! current penalty parameter used in the merit function

        contains

        procedure, public :: eval_merit => eval_merit_function
        procedure, public :: search     => line_search

    end type sqpopt_linesearch_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  evaluate the exact penalty (merit) function used to measure progress:
!
!  $$ \phi(x) = f(x) + \mu \lVert \max(c_l - c(x), 0, c(x) - c_u) \rVert $$

    subroutine eval_merit_function(me, f, c, c_lb, c_ub, phi)

    class(sqpopt_linesearch_type), intent(inout) :: me
    real(wp),                intent(in)  :: f      !! objective function value
    real(wp), dimension(:), intent(in)  :: c      !! constraint values `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp),                intent(out) :: phi    !! value of the merit function

    ! TODO: implement

    end subroutine eval_merit_function
!*******************************************************************************

!*******************************************************************************
!>
!  perform a line search along the direction `p` to find a step length
!  `alpha` \( \in (0,1] \) that gives sufficient decrease in the merit function.

    subroutine line_search(me, x, p, f, g, c, c_lb, c_ub, alpha, istat)

    class(sqpopt_linesearch_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: x      !! current point `dimension(n)`
    real(wp), dimension(:), intent(in)  :: p      !! search direction `dimension(n)`
    real(wp),                intent(in)  :: f      !! objective function value at `x`
    real(wp), dimension(:), intent(in)  :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(in)  :: c      !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp),                intent(out) :: alpha  !! accepted step length
    integer,                  intent(out) :: istat  !! status code (see [[sqpopt_types_module]])

    ! TODO: implement

    end subroutine line_search
!*******************************************************************************

    end module sqpopt_linesearch_module
!*******************************************************************************
