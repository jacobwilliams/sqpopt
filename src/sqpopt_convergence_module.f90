!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Stopping criteria for the SQP algorithm, based on the Karush-Kuhn-Tucker
!  (KKT) optimality conditions and feasibility of the constraints and bounds.

    module sqpopt_convergence_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    public :: check_convergence

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  check whether the current iterate satisfies the convergence criteria
!  (optimality, feasibility, and lack of progress tests).

    subroutine check_convergence(x, f, g, c, c_lb, c_ub, lambda, ktol, ctol, converged, istat)

    real(wp), dimension(:), intent(in)  :: x         !! current point `dimension(n)`
    real(wp),                intent(in)  :: f         !! objective function value at `x`
    real(wp), dimension(:), intent(in)  :: g         !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(in)  :: c         !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb      !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub      !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: lambda    !! Lagrange multipliers `dimension(m)`
    real(wp),                intent(in)  :: ktol      !! KKT optimality tolerance
    real(wp),                intent(in)  :: ctol      !! feasibility tolerance
    logical,                  intent(out) :: converged !! true if the convergence criteria are satisfied
    integer,                  intent(out) :: istat     !! status code (see [[sqpopt_types_module]])

    ! TODO: implement

    end subroutine check_convergence
!*******************************************************************************

    end module sqpopt_convergence_module
!*******************************************************************************
