!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Stopping criteria for the SQP algorithm, based on the Karush-Kuhn-Tucker
!  (KKT) optimality conditions and feasibility of the constraints and bounds.

    module sqpopt_convergence_module

    use sqpopt_kinds,         only: wp => sqpopt_module_wp
    use sqpopt_types_module,  only: sqpopt_sparse_matrix, sqpopt_success
    use sqpopt_linalg_module, only: sparse_matvec_transpose

    implicit none

    private

    public :: check_convergence

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  check whether the current iterate satisfies the convergence criteria:
!  dual feasibility / stationarity of the *projected* Lagrangian gradient
!  \( r = g - J^T \lambda \) (projected onto the active variable bounds,
!  as in the standard bound-constrained KKT/projected-gradient test), and
!  primal feasibility of the constraints and variable bounds.

    subroutine check_convergence(x, g, jac, c, x_lb, x_ub, c_lb, c_ub, lambda, ktol, ctol, converged, istat)

    real(wp), dimension(:),     intent(in)  :: x         !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)  :: g         !! objective gradient at `x` `dimension(n)`
    type(sqpopt_sparse_matrix), intent(in)  :: jac       !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)  :: c         !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:),     intent(in)  :: x_lb      !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)  :: x_ub      !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)  :: c_lb      !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)  :: c_ub      !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)  :: lambda    !! Lagrange multipliers `dimension(m)`
    real(wp),                   intent(in)  :: ktol      !! KKT optimality tolerance
    real(wp),                   intent(in)  :: ctol      !! feasibility tolerance
    logical,                     intent(out) :: converged !! true if the convergence criteria are satisfied
    integer,                     intent(out) :: istat     !! status code (see [[sqpopt_types_module]])

    real(wp), dimension(size(g)) :: jtlam, r
    real(wp) :: kkt_res, c_viol, x_viol, ri
    integer :: i

    call sparse_matvec_transpose(jac, lambda, jtlam)
    r = g - jtlam

    ! projected-gradient stationarity test: a nonzero reduced gradient
    ! component is only a violation if it points into the feasible region
    ! at an active bound (see e.g. the `pgtol` test used by L-BFGS-B):
    kkt_res = 0.0_wp
    do i = 1, size(x)
        if (x_ub(i)-x_lb(i) <= ctol) then
            ri = 0.0_wp                    !! fixed variable: no stationarity requirement
        else if (x(i)-x_lb(i) <= ctol) then
            ri = min(r(i), 0.0_wp)         !! at lower bound: only r(i)<0 is a violation
        else if (x_ub(i)-x(i) <= ctol) then
            ri = max(r(i), 0.0_wp)         !! at upper bound: only r(i)>0 is a violation
        else
            ri = r(i)                      !! free variable: full stationarity required
        end if
        kkt_res = max(kkt_res, abs(ri))
    end do

    c_viol  = maxval(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))
    x_viol  = maxval(max(x_lb-x, 0.0_wp) + max(x-x_ub, 0.0_wp))

    converged = kkt_res <= ktol .and. c_viol <= ctol .and. x_viol <= ctol
    istat = sqpopt_success

    end subroutine check_convergence
!*******************************************************************************

    end module sqpopt_convergence_module
!*******************************************************************************
