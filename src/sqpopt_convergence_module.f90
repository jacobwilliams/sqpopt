!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Stopping criteria for the SQP algorithm, based on the Karush-Kuhn-Tucker
!  (KKT) optimality conditions and feasibility of the constraints and
!  bounds, plus the stalled-progress and local-infeasibility tests.

    module sqpopt_convergence_module

    use sqpopt_kinds,         only: wp => sqpopt_module_wp
    use sqpopt_types_module,  only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_stalled, sqpopt_infeasible
    use sqpopt_linalg_module, only: sparse_matvec_transpose

    implicit none

    private

    public :: check_convergence

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  check whether the solver should stop at the current iterate, and why:
!
!  * `istat=sqpopt_success`: the KKT conditions hold. That is:
!    - stationarity of the *projected* Lagrangian gradient
!      \( r = g - J^T \lambda \) (projected onto the active variable
!      bounds, as in the standard bound-constrained projected-gradient
!      test);
!    - dual feasibility and complementarity of the constraint multipliers
!      (for the Lagrangian \( f - \lambda^T c \), \( \lambda_i>0 \) is
!      only allowed at a lower bound \( c_i=c_{l,i} \) and
!      \( \lambda_i<0 \) only at an upper bound \( c_i=c_{u,i} \), measured
!      as \( \lambda_i (c_i-c_{l,i}) \) resp. \( -\lambda_i (c_{u,i}-c_i) \);
!      equality constraints are free);
!    - primal feasibility of the constraints and variable bounds.
!
!    The stationarity and complementarity residuals are measured against
!    `ktol*s_d`, where \( s_d = \max(s_{max}, \lVert\lambda\rVert_1/m)/s_{max} \)
!    with \( s_{max}=100 \) (the scaling used by IPOPT): the tolerance is
!    only loosened when the multipliers are large on average, and then only
!    in proportion, so diverging multipliers (e.g. at a degenerate point)
!    can't mask a non-KKT point.
!  * `istat=sqpopt_stalled`: `f`, `f_prev`, `x_prev`, `ftol`, and `xtol` are
!    all supplied, the point is feasible, and the objective and the
!    variables have both stopped changing (by less than `ftol`/`xtol`,
!    relatively) since the previous iterate, without the KKT test having
!    been satisfied.
!  * `istat=sqpopt_infeasible`: only tested when `x_prev` is supplied (so
!    never at the initial point). The constraints are violated, but the
!    point is stationary for the constraint violation
!    \( \tfrac12 \lVert r_c \rVert_2^2 \) (\( r_c \) the signed violation of
!    each constraint's bounds) subject to the variable bounds, i.e. the
!    projected gradient \( J^T r_c \) is below `ktol*|r_c|`. So no
!    first-order progress toward feasibility is possible from here. If
!    `viol_prev` (the violation at `x_prev`) is supplied, the violation
!    must also have stopped decreasing (by less than 1% since `x_prev`): a
!    point can be stationary for the violation without being a minimum of
!    it (e.g. `J=0` at a maximum, or where the linearization is
!    degenerate), and the solver can still make progress from there to
!    second order, so it only stops once it has failed to.
!
!  `converged` is true in all three cases (i.e. it means "stop here").

    subroutine check_convergence(x, g, jac, c, x_lb, x_ub, c_lb, c_ub, lambda, ktol, ctol, converged, istat, &
                                  f, f_prev, x_prev, ftol, xtol, kkt_error, feas_error, viol_prev, &
                                  dual_inf_tol, f_scale)

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
    logical,                     intent(out) :: converged !! true if the solver should stop at `x` (see `istat` for why)
    integer,                     intent(out) :: istat     !! why: `sqpopt_success`, `sqpopt_stalled`, or `sqpopt_infeasible`
                                                          !! if `converged`, else `sqpopt_success` (see [[sqpopt_types_module]])
    real(wp),                     optional, intent(in) :: f      !! objective value at `x` (enables the stalled-progress test with `f_prev`/`ftol`)
    real(wp),                     optional, intent(in) :: f_prev !! objective value at the previous iterate
    real(wp), dimension(:),       optional, intent(in) :: x_prev !! the previous iterate `dimension(n)` (enables the
                                                                 !! stalled-progress test with `xtol`, and the infeasibility test)
    real(wp),                     optional, intent(in) :: ftol   !! relative objective-change tolerance for the stalled-progress test
    real(wp),                     optional, intent(in) :: xtol   !! relative variable-change tolerance for the stalled-progress test
    real(wp),                     optional, intent(out) :: kkt_error  !! the KKT error: the larger of the stationarity and
                                                                      !! complementarity residuals, divided by their scaling
                                                                      !! (so it is directly comparable to `ktol`)
    real(wp),                     optional, intent(out) :: feas_error !! the largest violation of a constraint or variable bound
    real(wp),                     optional, intent(in)  :: viol_prev  !! the largest constraint violation at `x_prev`
                                                                      !! (enables the no-progress condition of the
                                                                      !! infeasibility test)
    real(wp),                     optional, intent(in)  :: dual_inf_tol !! if present (with `f_scale`), the KKT test
                                                                        !! also requires the stationarity residual of
                                                                        !! the *unscaled* problem to be at most this
                                                                        !! (as IPOPT's `dual_inf_tol`)
    real(wp),                     optional, intent(in)  :: f_scale      !! the objective's scale factor: the unscaled
                                                                        !! stationarity residual is the scaled one
                                                                        !! divided by it

    real(wp), parameter :: s_max = 100.0_wp !! multiplier-scaling threshold (see above)

    real(wp), dimension(size(g)) :: jtlam, r
    real(wp), dimension(size(c)) :: rc
    real(wp) :: kkt_res, dual_res, c_viol, x_viol, rel_f, rel_x, lam_scale
    integer :: i

    istat     = sqpopt_success
    converged = .false.

    call sparse_matvec_transpose(jac, lambda, jtlam)
    r = g - jtlam

    ! projected-gradient stationarity test: a nonzero reduced gradient
    ! component is only a violation if it points into the feasible region
    ! at an active bound (see e.g. the `pgtol` test used by L-BFGS-B):
    kkt_res = projected_inf_norm(r, x, x_lb, x_ub, ctol)

    ! dual feasibility & complementarity of the constraint multipliers:
    dual_res = 0.0_wp
    do i = 1, size(c)
        if (c_ub(i)-c_lb(i) <= ctol) cycle  ! equality constraint: free multiplier
        if (lambda(i) > 0.0_wp) then
            dual_res = max(dual_res, lambda(i)*max(c(i)-c_lb(i), 0.0_wp))
        else if (lambda(i) < 0.0_wp) then
            dual_res = max(dual_res, -lambda(i)*max(c_ub(i)-c(i), 0.0_wp))
        end if
    end do

    lam_scale = 1.0_wp
    if (size(lambda) > 0) lam_scale = max(s_max, sum(abs(lambda))/size(lambda))/s_max

    c_viol = 0.0_wp
    if (size(c) > 0) c_viol = maxval(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))
    x_viol = maxval(max(x_lb-x, 0.0_wp) + max(x-x_ub, 0.0_wp))
    if (present(kkt_error))  kkt_error  = max(kkt_res, dual_res)/lam_scale
    if (present(feas_error)) feas_error = max(c_viol, x_viol)

    if (kkt_res <= ktol*lam_scale .and. dual_res <= ktol*lam_scale .and. &
        c_viol <= ctol .and. x_viol <= ctol .and. unscaled_ok()) then
        converged = .true.
        return
    end if

    if (c_viol <= ctol .and. x_viol <= ctol .and. &
        present(f) .and. present(f_prev) .and. present(x_prev) .and. present(ftol) .and. present(xtol)) then
        ! feasible, but the KKT test hasn't yet reached `ktol`: stop anyway
        ! if the objective and the variables have both stalled (relative to
        ! the previous iterate), rather than looping until `max_iter` on
        ! marginal, ever-shrinking steps:
        rel_f = abs(f-f_prev)/max(1.0_wp, abs(f))
        rel_x = norm2(x-x_prev)/max(1.0_wp, norm2(x))
        if (rel_f <= ftol .and. rel_x <= xtol) then
            converged = .true.
            istat     = sqpopt_stalled
            return
        end if
    end if

    if (c_viol > ctol .and. present(x_prev)) then
        if (present(viol_prev)) then
            if (c_viol < 0.99_wp*viol_prev) return   ! (still making progress toward feasibility)
        end if
        ! infeasible: is `x` stationary for the (squared, l2) constraint
        ! violation, subject to the variable bounds? Then no first-order
        ! progress toward feasibility is possible from here. The gradient
        ! `J^T r_c` is measured relative to `|r_c|` (it is proportional to
        ! it), so a nearly-feasible point isn't mistaken for a stationary one:
        rc = min(c-c_lb, 0.0_wp) + max(c-c_ub, 0.0_wp)
        call sparse_matvec_transpose(jac, rc, r)
        if (projected_inf_norm(r, x, x_lb, x_ub, ctol) <= ktol*norm2(rc)) then
            converged = .true.
            istat     = sqpopt_infeasible
        end if
    end if

    contains

        logical function unscaled_ok()
        !! the unscaled stationarity test (see `dual_inf_tol`): with the
        !! objective scaled by `f_scale` and each constraint by its own
        !! factor, the multipliers scale by their ratio, so the gradient of
        !! the unscaled Lagrangian is the scaled one divided by `f_scale`
        unscaled_ok = .true.
        if (present(dual_inf_tol) .and. present(f_scale)) &
            unscaled_ok = kkt_res <= dual_inf_tol*f_scale
        end function unscaled_ok

    end subroutine check_convergence
!*******************************************************************************

!*******************************************************************************
!>
!  the infinity norm of the gradient `r`, projected onto the variable
!  bounds: a component only counts if moving along `-r` would stay feasible
!  (free variable, or at a bound but pointing into the interior).

    pure function projected_inf_norm(r, x, x_lb, x_ub, tol) result(res)

    real(wp), dimension(:), intent(in) :: r     !! gradient `dimension(n)`
    real(wp), dimension(:), intent(in) :: x     !! current point `dimension(n)`
    real(wp), dimension(:), intent(in) :: x_lb  !! variable lower bounds `dimension(n)`
    real(wp), dimension(:), intent(in) :: x_ub  !! variable upper bounds `dimension(n)`
    real(wp),               intent(in) :: tol   !! distance within which a bound counts as active
    real(wp) :: res

    real(wp) :: ri
    integer :: i

    res = 0.0_wp
    do i = 1, size(x)
        if (x_ub(i)-x_lb(i) <= tol) then
            ri = 0.0_wp                    !! fixed variable: no stationarity requirement
        else if (x(i)-x_lb(i) <= tol) then
            ri = min(r(i), 0.0_wp)         !! at lower bound: only r(i)<0 is a violation
        else if (x_ub(i)-x(i) <= tol) then
            ri = max(r(i), 0.0_wp)         !! at upper bound: only r(i)>0 is a violation
        else
            ri = r(i)                      !! free variable: full stationarity required
        end if
        res = max(res, abs(ri))
    end do

    end function projected_inf_norm
!*******************************************************************************

    end module sqpopt_convergence_module
!*******************************************************************************
