!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Feasibility restoration, used by [[sqpopt_iterate_module]] when the QP
!  subproblem reports that the linearized constraints are inconsistent
!  (`sqpopt_infeasible`), so that no QP step can be trusted, or when the
!  filter or funnel line search, or the trust region, finds no acceptable
!  step at an infeasible point (one whose largest constraint violation, as
!  the convergence test measures it, exceeds `options%ctol`: a round-off
!  violation has nothing to restore). Two strategies are available
!  (`options%restoration_mode`):
!
!  * `sqpopt_restoration_phase` (**default**): when the filter or funnel
!    line search, or the trust region, finds no acceptable step, a
!    **feasibility restoration phase**, as in filter-SQP methods (Fletcher
!    & Leyffer; Wächter & Biegler; and the Uno solver's
!    `FeasibilityRestoration`). (An inconsistent QP is still handled by the
!    Gauss-Newton step below, falling back to the QP's elastic step, which
!    also accounts for the objective: on the Hock-Schittkowski problems
!    this finds better local solutions than a pure feasibility phase.) The
!    solver switches to minimizing the constraint violation for as many
!    major iterations as needed (see [[restoration_phase_step]]): each one
!    solves a *feasibility QP* with the regular QP solvers, whose objective
!    is a proximal term toward the point where the phase started, and
!    whose linearized constraints are enforced (or, if they are
!    inconsistent, their \( \ell_1 \) violation is minimized by the QP's
!    elastic mode), followed by a backtracking search on the \( \ell_1 \)
!    violation. The phase ends (see [[restoration_phase_done]]) once the
!    violation has dropped below `restoration_exit_factor` times its value
!    at the start of the phase *and* the point is acceptable to the filter
!    (or funnel), when the point is feasible, after
!    `restoration_max_iter` iterations, or when a phase step (and its
!    Gauss-Newton fallback) fails, so that the next iteration tries the
!    optimality QP again. The point where the phase started
!    is added to the filter (or the funnel is tightened toward it), so the
!    iterations can't cycle back to it.
!  * `sqpopt_restoration_gauss_newton`: a single Gauss-Newton step on the
!    constraint violation each time (the original, lightweight strategy;
!    with the filter or funnel line search, the current point is added to
!    the filter, or the funnel tightened, first):
!
!  $$ \min_x \; \tfrac12 \lVert r_c(x) \rVert_2^2 \quad \text{s.t.} \quad x_l \le x \le x_u $$
!
!  where \( r_c \) is the signed violation of each constraint's bounds.
!  Repeated restoration steps converge to either a feasible point (after
!  which the normal SQP iterations resume) or a point that is stationary
!  for the violation, which [[check_convergence]] then reports as
!  `sqpopt_infeasible`. This step is also the fallback when a restoration
!  phase step fails.
!
!  In both cases, before a point that is stationary for the violation is
!  reported as infeasible, [[escape_step]] looks for a second-order
!  decrease of the violation.

    module sqpopt_restoration_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_line_search_failed, sqpopt_all_finite, &
                                     sqpopt_infeasible, sqpopt_qp_solve_failed
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, l1_violation
    use sqpopt_linalg_module,     only: sparse_matvec
    use lsqr_module,           only: lsqr_solver_ez
    use sqpopt_least_squares_module, only: sqpopt_least_squares_type

    implicit none

    private

    public :: restoration_step, escape_step

    integer, parameter, public :: sqpopt_restoration_phase        = 1 !! (default) a feasibility restoration phase
                                                                      !! (see the module documentation)
    integer, parameter, public :: sqpopt_restoration_gauss_newton = 2 !! a single Gauss-Newton step on the violation

    type, public :: sqpopt_restoration_type
        !! the state of the feasibility restoration phase (internal; reset on each `solve`)
        logical  :: active    = .false.  !! whether the solver is in the restoration phase
        integer  :: n_iter    = 0        !! number of iterations of the current phase
        integer  :: n_phases  = 0        !! number of restoration phases so far
        real(wp) :: theta_ref = 0.0_wp   !! \( \ell_1 \) violation where the phase started
        real(wp), dimension(:), allocatable :: x_ref !! the point where the phase started (the proximal center)
        type(sqpopt_qp_solver_type) :: qp   !! the QP solver for the feasibility QPs (a copy of the main one,
                                            !! so that the latter's warm start is kept for the optimality QPs)
        type(sqpopt_hessian_type)   :: hess !! the feasibility QP's Hessian (the proximal term, \( \zeta I \))
        contains
        procedure, public :: enter => restoration_phase_enter
        procedure, public :: step  => restoration_phase_step
        procedure, public :: done  => restoration_phase_done
    end type sqpopt_restoration_type

    real(wp), parameter :: zeta = 1.0_wp !! weight \( \zeta \) of the feasibility QP's proximal term (only its ratio
                                         !! to the QP's elastic weight matters: with consistent linearized constraints
                                         !! the step is the one closest to the proximal center)

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  take one Gauss-Newton step toward feasibility of the constraints (see
!  the module-level documentation), with a backtracking Armijo search on
!  \( \tfrac12 \lVert r_c \rVert_2^2 \). If no decrease is found, `x_new=x`
!  and `istat=sqpopt_line_search_failed` (also if `LSQR` stops at its
!  iteration limit before finding the Gauss-Newton step). With
!  `least_squares` (if it is enabled), the Gauss-Newton step is computed by
!  a direct solve, and by `LSQR` only if that fails.
!
!  If `direction` is present, it is searched along instead of the
!  Gauss-Newton direction (for when that can't make progress: at a point
!  that is stationary for the violation, a decrease may still be possible
!  to second order, e.g. along the QP's elastic step where `J=0`). Along a
!  direction without first-order decrease, any decrease of the violation
!  is accepted.

    subroutine restoration_step(problem, jac, x, c, max_step, x_new, alpha, istat, direction, least_squares)

    type(sqpopt_problem_type),  intent(inout) :: problem  !! problem definition
    type(sqpopt_sparse_matrix), intent(in)    :: jac      !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x        !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c        !! constraint values at `x` `dimension(m)`
    real(wp),                   intent(in)    :: max_step !! cap on \( \lVert p \rVert_2 \)
    real(wp), dimension(:),     intent(out)   :: x_new    !! new point `dimension(n)`
    real(wp),                   intent(out)   :: alpha    !! accepted step length (`0` if none)
    integer,                    intent(out)   :: istat    !! status code (see [[sqpopt_types_module]])
    real(wp), dimension(:), optional, intent(in) :: direction !! search direction to use instead of Gauss-Newton `dimension(n)`
    type(sqpopt_least_squares_type), optional, intent(inout) :: least_squares !! the direct least-squares solver
                                                                              !! (see [[sqpopt_least_squares_module]])

    real(wp), parameter :: sigma     = 1.0e-4_wp !! Armijo sufficient-decrease parameter
    real(wp), parameter :: backtrack = 0.5_wp    !! step-length reduction factor
    integer,  parameter :: max_ls    = 30        !! maximum number of backtracking steps
    integer,  parameter :: lsqr_itnlim_stop = 5  !! `LSQR`'s `istop` for "iteration limit reached"

    real(wp), dimension(size(c)) :: rc, c_trial
    real(wp), dimension(size(x)) :: p, x_trial
    real(wp) :: h0, h_trial, dh0
    type(lsqr_solver_ez) :: lsqr
    integer :: istop, it
    logical :: solved

    rc = violation(c, problem%c_lb, problem%c_ub)
    h0 = 0.5_wp*dot_product(rc, rc)

    ! Gauss-Newton step: minimum-norm solution of J*p = -r_c, then made to
    ! respect the variable bounds and the step-length cap:
    if (present(direction)) then
        p = direction
    else
        solved = .false.
        if (present(least_squares)) then
            block
                logical, dimension(size(c)) :: all_rows
                all_rows = .true.
                call least_squares%min_norm(jac, all_rows, -rc, p, solved)
            end block
        end if
        if (.not. solved) then
            call lsqr%initialize(problem%m, problem%n, jac%val, jac%irow, jac%icol, &
                                 itnlim=2*(problem%m+problem%n)+10)
            call lsqr%solve(-rc, 0.0_wp, p, istop)
            if (istop == lsqr_itnlim_stop .or. .not. sqpopt_all_finite(p)) then
                ! LSQR didn't converge: no usable Gauss-Newton step
                alpha = 0.0_wp
                x_new = x
                istat = sqpopt_line_search_failed
                return
            end if
        end if
    end if
    p = min(max(x+p, problem%x_lb), problem%x_ub) - x
    if (norm2(p) > max_step) p = p*(max_step/norm2(p))

    ! directional derivative of 0.5*|r_c|^2 along p: r_c^T J p
    block
        real(wp), dimension(size(c)) :: jp
        integer :: k
        jp = 0.0_wp
        do k = 1, jac%nnz
            jp(jac%irow(k)) = jp(jac%irow(k)) + jac%val(k)*p(jac%icol(k))
        end do
        dh0 = dot_product(rc, jp)
    end block

    alpha = 1.0_wp
    if (present(direction)) dh0 = min(dh0, 0.0_wp)   ! (then any decrease is accepted)
    if (dh0 < 0.0_wp .or. present(direction)) then
        do it = 1, max_ls
            x_trial = x + alpha*p
            call problem%c(x_trial, c_trial)
            if (sqpopt_all_finite(c_trial)) then
                rc = violation(c_trial, problem%c_lb, problem%c_ub)
                h_trial = 0.5_wp*dot_product(rc, rc)
            else
                h_trial = huge(1.0_wp)  ! a non-finite trial point is always rejected
            end if
            if (h_trial <= h0 + sigma*alpha*dh0 .and. h_trial < h0) then
                x_new = x_trial
                istat = sqpopt_success
                return
            end if
            alpha = backtrack*alpha
        end do
    end if

    ! no decrease in the violation along the Gauss-Newton direction:
    alpha = 0.0_wp
    x_new = x
    istat = sqpopt_line_search_failed

    end subroutine restoration_step
!*******************************************************************************

!*******************************************************************************
!>
!  look for a second-order decrease of the constraint violation from a point
!  `x` that is stationary for it (where [[check_convergence]] would report
!  `sqpopt_infeasible`). Such a point may be a saddle of the violation rather
!  than a minimum: e.g. on a symmetry plane of the problem (`x_j=0`, with
!  every function even in `x_j`), where the gradients' `x_j` components all
!  vanish, so no first-order method, and no exactly computed step, ever
!  leaves the plane.
!
!  Each variable whose column of the Jacobian is negligible in the violated
!  rows (up to `max_probe` of them) is perturbed by `t*max(1,|x_j|)` in each
!  direction, for `t` = `1e-3`, `1e-2`, `1e-1` (within its bounds). The
!  first trial point with a lower (squared, l2) violation is returned in
!  `x_new`, with `istat=sqpopt_success`; if there is none, `x_new=x` and
!  `istat=sqpopt_line_search_failed`. This costs at most `6*max_probe`
!  evaluations of the constraints.
!
!  A column is negligible if its largest element (in the violated rows) is
!  below `tol`, the tolerance at which [[check_convergence]] takes the
!  gradient of the violation to be zero, or below `sqrt(epsilon)` times the
!  largest element of all the columns. The first test doesn't depend on the
!  working precision: the iterates approach such a plane without landing on
!  it exactly (in `real64`, round-off usually pushes them off it instead;
!  in `real128` it doesn't, and e.g. TP88 stops next to the plane `x_2=0`
!  with a column of `1e-8`, far above `sqrt(epsilon)`).

    subroutine escape_step(problem, jac, x, c, tol, x_new, istat)

    type(sqpopt_problem_type),  intent(inout) :: problem  !! problem definition
    type(sqpopt_sparse_matrix), intent(in)    :: jac      !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x        !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c        !! constraint values at `x` `dimension(m)`
    real(wp),                   intent(in)    :: tol      !! a Jacobian column is negligible below this (the convergence
                                                          !! test's `ktol`)
    real(wp), dimension(:),     intent(out)   :: x_new    !! new point `dimension(n)`
    integer,                    intent(out)   :: istat    !! status code (see [[sqpopt_types_module]])

    integer,  parameter :: max_probe = 10 !! maximum number of variables probed
    real(wp), parameter :: steps(3) = [1.0e-3_wp, 1.0e-2_wp, 1.0e-1_wp] !! relative perturbations tried
    real(wp), parameter :: col_tol = sqrt(epsilon(1.0_wp)) !! a column is also negligible below this (relative)

    real(wp), dimension(size(c)) :: rc, c_trial
    real(wp), dimension(size(x)) :: colmax, x_trial
    real(wp) :: h0, h_trial, jmax, xj
    integer :: j, k, i, s, n_probe

    x_new = x
    istat = sqpopt_line_search_failed

    rc = violation(c, problem%c_lb, problem%c_ub)
    h0 = 0.5_wp*dot_product(rc, rc)
    if (h0 <= 0.0_wp) return

    ! largest Jacobian element of each variable in the violated rows:
    colmax = 0.0_wp
    do k = 1, jac%nnz
        if (rc(jac%irow(k)) /= 0.0_wp) colmax(jac%icol(k)) = max(colmax(jac%icol(k)), abs(jac%val(k)))
    end do
    jmax = maxval(colmax)

    n_probe = 0
    do j = 1, size(x)
        if (colmax(j) > max(tol, col_tol*jmax)) cycle
        n_probe = n_probe + 1
        if (n_probe > max_probe) exit
        do i = 1, size(steps)
            do s = 1, -1, -2
                xj = min(max(x(j) + s*steps(i)*max(1.0_wp, abs(x(j))), problem%x_lb(j)), problem%x_ub(j))
                if (xj == x(j)) cycle
                x_trial    = x
                x_trial(j) = xj
                call problem%c(x_trial, c_trial)
                if (.not. sqpopt_all_finite(c_trial)) cycle
                rc = violation(c_trial, problem%c_lb, problem%c_ub)
                h_trial = 0.5_wp*dot_product(rc, rc)
                if (h_trial < h0) then
                    x_new = x_trial
                    istat = sqpopt_success
                    return
                end if
            end do
        end do
    end do

    end subroutine escape_step
!*******************************************************************************

!*******************************************************************************
!>
!  start a feasibility restoration phase at `x` (with \( \ell_1 \)
!  violation `theta`): record where it started, and set up the feasibility
!  QP's solver (a copy of `qp_solver`) and Hessian.

    subroutine restoration_phase_enter(me, x, theta, qp_solver)

    class(sqpopt_restoration_type), intent(inout) :: me
    real(wp), dimension(:),         intent(in)    :: x         !! the current point `dimension(n)`
    real(wp),                       intent(in)    :: theta     !! its \( \ell_1 \) constraint violation
    type(sqpopt_qp_solver_type),    intent(in)    :: qp_solver !! the main QP solver (copied)

    me%active    = .true.
    me%n_iter    = 0
    me%n_phases  = me%n_phases + 1
    me%theta_ref = theta
    me%x_ref     = x
    me%qp        = qp_solver
    call me%hess%initialize(size(x), 1, scale0=zeta)

    end subroutine restoration_phase_enter
!*******************************************************************************

!*******************************************************************************
!>
!  one iteration of the feasibility restoration phase from `x`: solve the
!  feasibility QP
!
!  $$ \min_p \; \zeta (x-x_{ref})^T p + \tfrac12 \zeta \lVert p \rVert^2
!     \quad \text{s.t.} \quad c_l \le c + J p \le c_u, \quad x_l \le x+p \le x_u $$
!
!  (whose elastic mode minimizes the \( \ell_1 \) violation of the
!  linearized constraints if they are inconsistent), then backtrack along
!  `p` until the \( \ell_1 \) violation \( \theta \) satisfies the Armijo
!  condition \( \theta(x+\alpha p) \le \theta(x) - \eta\alpha\,\text{pred} \),
!  where \( \text{pred} = \theta(c) - \theta(c+Jp) \) is the decrease
!  predicted by the linearization. If the QP predicts no decrease, or no
!  step length is accepted, `x_new=x` and
!  `istat=sqpopt_line_search_failed`.

    subroutine restoration_phase_step(me, problem, jac, x, c, x_new, alpha, istat)

    class(sqpopt_restoration_type), intent(inout) :: me
    type(sqpopt_problem_type),  intent(inout) :: problem  !! problem definition
    type(sqpopt_sparse_matrix), intent(in)    :: jac      !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x        !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c        !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:),     intent(out)   :: x_new    !! new point `dimension(n)`
    real(wp),                   intent(out)   :: alpha    !! accepted step length (`0` if none)
    integer,                    intent(out)   :: istat    !! status code (see [[sqpopt_types_module]])

    real(wp), parameter :: eta       = 1.0e-4_wp !! Armijo constant
    real(wp), parameter :: alpha_min = 1.0e-8_wp !! smallest step length tried

    real(wp), dimension(size(x)) :: g_r, p, x_trial
    real(wp), dimension(size(c)) :: lambda_r, jp, c_trial
    real(wp) :: theta0, theta_t, pred
    integer :: qp_istat

    me%n_iter = me%n_iter + 1
    x_new = x
    alpha = 0.0_wp
    istat = sqpopt_line_search_failed

    theta0 = l1_violation(c, problem%c_lb, problem%c_ub)
    g_r = zeta*(x - me%x_ref)
    call me%qp%solve(me%hess, jac, x, g_r, c, problem%x_lb, problem%x_ub, problem%c_lb, problem%c_ub, &
                     p, lambda_r, qp_istat)
    if (qp_istat /= sqpopt_success .and. qp_istat /= sqpopt_infeasible .and. qp_istat /= sqpopt_qp_solve_failed) return

    ! the decrease in the violation predicted by the linearization:
    call sparse_matvec(jac, p, jp)
    pred = theta0 - l1_violation(c + jp, problem%c_lb, problem%c_ub)
    if (.not. (pred > 1.0e-12_wp*max(1.0_wp, theta0))) return

    alpha = 1.0_wp
    do
        x_trial = x + alpha*p
        call problem%c(x_trial, c_trial)
        if (sqpopt_all_finite(c_trial)) then
            theta_t = l1_violation(c_trial, problem%c_lb, problem%c_ub)
            if (theta_t <= theta0 - eta*alpha*pred) then
                x_new = x_trial
                istat = sqpopt_success
                return
            end if
        end if
        alpha = 0.5_wp*alpha
        if (alpha < alpha_min) exit
    end do
    alpha = 0.0_wp

    end subroutine restoration_phase_step
!*******************************************************************************

!*******************************************************************************
!>
!  whether the restoration phase can end at a point with \( \ell_1 \)
!  violation `theta` and objective `f`: the violation is at most
!  `exit_factor` times its value where the phase started and the point is
!  acceptable to the filter or funnel (see [[sqpopt_linesearch_type]]'s
!  `globalization_acceptable`); or the point is feasible
!  (\( \theta \le \) `feas_tol`); or the phase has taken `max_iter`
!  iterations. Ends the phase (`active=.false.`) if so.

    function restoration_phase_done(me, linesearch, theta, f, exit_factor, feas_tol, max_iter) result(done)

    class(sqpopt_restoration_type), intent(inout) :: me
    type(sqpopt_linesearch_type),   intent(in)    :: linesearch  !! supplies the filter or funnel
    real(wp),                       intent(in)    :: theta       !! \( \ell_1 \) violation at the new point
    real(wp),                       intent(in)    :: f           !! objective at the new point
    real(wp),                       intent(in)    :: exit_factor !! required reduction of the violation
    real(wp),                       intent(in)    :: feas_tol    !! a violation this small ends the phase regardless
    integer,                        intent(in)    :: max_iter    !! maximum number of iterations of a phase
    logical :: done

    done = (theta <= exit_factor*me%theta_ref .and. linesearch%globalization_acceptable(theta, f)) &
           .or. theta <= feas_tol .or. me%n_iter >= max_iter
    if (done) me%active = .false.

    end function restoration_phase_done
!*******************************************************************************

!*******************************************************************************
!>
!  the signed violation of the constraint bounds: `c-c_lb` below the lower
!  bound, `c-c_ub` above the upper bound, and `0` in between.

    pure function violation(c, c_lb, c_ub) result(rc)

    real(wp), dimension(:), intent(in) :: c    !! constraint values `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_lb !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_ub !! constraint upper bounds `dimension(m)`
    real(wp), dimension(size(c)) :: rc

    rc = min(c-c_lb, 0.0_wp) + max(c-c_ub, 0.0_wp)

    end function violation
!*******************************************************************************

    end module sqpopt_restoration_module
!*******************************************************************************
