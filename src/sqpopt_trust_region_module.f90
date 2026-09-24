!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Opt-in **trust-region** globalization strategy for the SQP major
!  iteration -- an alternative to the line-search-based approach in
!  [[sqpopt_linesearch_module]] (see `plan/TRUST_REGION_PLAN.md` for the
!  full design rationale). Disabled by default (`sqpopt_trust_region_type
!  %enabled = .false.`), in which case `sqpopt_iterate_module` behaves
!  exactly as if this module did not exist.
!
!  When enabled, a major iteration no longer solves the linearized QP
!  once and then searches for a step length `alpha` along the resulting
!  `p`. Instead, the QP is re-solved as needed with a **trust-region
!  radius**: the variable bounds passed to `qp_solver%solve` are
!  temporarily tightened to `max(x_lb,x-radius)`/`min(x_ub,x+radius)`
!  (an \( \ell_\infty \) ball around `x`, intersected with the problem's
!  own bounds -- exactly the box constraint Fletcher & Leyffer's own
!  trust-region QP uses, and for the same reason: it keeps the QP a
!  simple bound-constrained problem, not a genuinely quadratically-
!  constrained one). No changes are needed to any of the three
!  `qp_solver_mode`s to support this, since they all already enforce
!  variable bounds.
!
!  The resulting step `p` (after the same second-order correction used by
!  the line-search path, see [[sqpopt_soc_module]]) is tested for
!  acceptance, and the radius is grown or shrunk accordingly:
!
!  * if `linesearch%mode == sqpopt_linesearch_filter`, acceptance reuses
!    the filter's own `(f,h)` domination test (`linesearch%filter_test`)
!    -- this combination is the *literal* Fletcher & Leyffer filter-SQP
!    algorithm (`references/fletcher.pdf`), of which
!    `sqpopt_linesearch_filter` on its own is only a line-search
!    adaptation.
!  * otherwise, acceptance uses the classical trust-region-SQP ratio test
!    (Nocedal & Wright, *Numerical Optimization*, Ch. 18): \( \rho =
!    \text{ared}/\text{pred} \), the ratio of the actual to the
!    predicted decrease in the merit function selected by
!    `linesearch%merit_mode`, accepted when \( \rho \ge \eta_1 \).
!
!  If a step is rejected, the radius is shrunk (`shrink_factor`) and the
!  QP is re-solved (bounded by `max_retries` retries per major iteration)
!  -- there is no `alpha` to backtrack along, since a genuinely different
!  `p` is computed each time. If every retry is rejected, the last
!  (smallest-radius) trial point is accepted anyway
!  (`istat=sqpopt_line_search_failed`), the same "accept it anyway rather
!  than freeze" fallback [[sqpopt_linesearch_module]]'s `alpha_min` floor uses.
!
!  Because this bypasses `linesearch%search` entirely, the line-search-
!  specific options (`major_step_limit`, `alpha_min`/`sigma`/`backtrack`/
!  `max_ls_iter`, and the watchdog fields) are not used when trust region
!  is enabled -- there is no `alpha`, so nothing to cap or backtrack.

    module sqpopt_trust_region_module

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_types_module,      only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_line_search_failed
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_hessian_module,    only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_filter, &
                                         l1_violation, filter_penalty_estimate
    use sqpopt_linalg_module,     only: sparse_matvec
    use sqpopt_soc_module,        only: second_order_correction

    implicit none

    private

    type, public :: sqpopt_trust_region_type
        !! options and state for trust-region-based globalization
        !! (an alternative to `linesearch%search`, see module docs).

        logical  :: enabled       = .false.    !! if true, use trust-region radius management instead of a line search
        real(wp) :: radius0       = 1.0_wp      !! initial trust-region radius
        real(wp) :: radius_min    = 1.0e-8_wp   !! below this, a major iteration's retries give up (no step taken)
        real(wp) :: radius_max    = 1.0e3_wp    !! ceiling on the radius
        real(wp) :: eta1          = 0.1_wp      !! ratio threshold to accept a step (merit-ratio acceptance only)
        real(wp) :: eta2          = 0.75_wp     !! ratio threshold to also grow the radius (merit-ratio acceptance only)
        real(wp) :: shrink_factor = 0.5_wp      !! `radius *= shrink_factor` on a rejected step
        real(wp) :: expand_factor = 2.0_wp      !! `radius *= expand_factor` on a step that used the full radius and was accepted
        integer  :: max_retries   = 20          !! maximum QP re-solves (with a shrinking radius) per major iteration

        ! internal state (persists across major iterations, like the watchdog/filter state on `sqpopt_linesearch_type`):
        logical  :: ready  = .false. !! whether `radius` has been initialized from `radius0` yet
        real(wp) :: radius = 0.0_wp  !! current trust-region radius

        contains

        procedure, public :: step => trust_region_step

    end type sqpopt_trust_region_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  attempt to find and accept a trust-region step, re-solving the QP
!  subproblem with a shrinking radius as needed (see the module-level
!  documentation for the full algorithm). Updates `me%radius` in place.

    subroutine trust_region_step(me, problem, hessian, qp_solver, linesearch, &
                                  x, g, f, c, jac, x_new, new_lambda, alpha, istat)

    class(sqpopt_trust_region_type), intent(inout) :: me
    type(sqpopt_problem_type),      intent(inout) :: problem
    type(sqpopt_hessian_type),      intent(inout) :: hessian     !! matrix-free Hessian approximation
    type(sqpopt_qp_solver_type),    intent(inout) :: qp_solver
    type(sqpopt_linesearch_type),   intent(inout) :: linesearch  !! supplies `merit_mode`/`eval_merit` (ratio test) or
                                                                  !! the filter (`mode==sqpopt_linesearch_filter`)
    real(wp), dimension(:),         intent(in)  :: x       !! current point `dimension(n)`
    real(wp), dimension(:),         intent(in)  :: g       !! objective gradient at `x` `dimension(n)`
    real(wp),                       intent(in)  :: f       !! objective value at `x`
    real(wp), dimension(:),         intent(in)  :: c       !! constraint values at `x` `dimension(m)`
    type(sqpopt_sparse_matrix),     intent(in)  :: jac     !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:),         intent(out) :: x_new      !! new point `dimension(n)` (the last, smallest-radius trial point if no retry was clearly accepted)
    real(wp), dimension(:),         intent(out) :: new_lambda !! new multipliers `dimension(m)`, from the QP solve that produced `x_new`
    real(wp),                       intent(out) :: alpha      !! `1` if a step was accepted, `0` otherwise (informational only -- there
                                                               !! is no line-search step length under trust-region globalization)
    integer,                        intent(out) :: istat      !! status code (see [[sqpopt_types_module]])

    integer :: retry, qp_istat
    logical :: use_filter, accept, both_feasible
    real(wp), dimension(size(x)) :: p, x_lb2, x_ub2, hp
    real(wp), dimension(size(c)) :: jp, c_trial, c_lin
    real(wp) :: f_trial, h0, h_trial, h_lin, q, pred, ared, ratio, mu, phi0, phi_trial

    if (.not. me%ready) then
        me%radius = me%radius0
        me%ready  = .true.
    end if

    use_filter    = (linesearch%mode == sqpopt_linesearch_filter)
    h0            = l1_violation(c, problem%c_lb, problem%c_ub)
    both_feasible = h0 <= linesearch%filter_feas_tol
    if (use_filter) call linesearch%filter_prepare(h0)

    do retry = 1, me%max_retries

        x_lb2 = max(problem%x_lb, x - me%radius)
        x_ub2 = min(problem%x_ub, x + me%radius)

        call qp_solver%solve(hessian, jac, x, g, c, x_lb2, x_ub2, problem%c_lb, problem%c_ub, p, new_lambda, qp_istat)

        if (problem%m > 0) call second_order_correction(problem, linesearch, jac, x, c, new_lambda, p)

        ! keep the merit function's penalty parameter dominating the current
        ! multiplier estimates, same rule as the line-search path (needed for
        ! `eval_merit`'s ared/pred to be meaningful in the merit-ratio branch below):
        if (size(new_lambda) > 0) linesearch%penalty = max(linesearch%penalty, maxval(abs(new_lambda)) + 1.0_wp)

        call hessian%hv_product(p, hp)
        q = -(dot_product(g, p) + 0.5_wp*dot_product(p, hp))

        call problem%eval_f(x+p, f_trial)
        call problem%eval_c(x+p, c_trial)
        h_trial = l1_violation(c_trial, problem%c_lb, problem%c_ub)

        if (use_filter) then

            accept = linesearch%filter_test(f_trial, h_trial)
            if (accept .and. both_feasible .and. h_trial <= linesearch%filter_feas_tol) accept = f_trial < f
            ratio = 1.0_wp !! not used for the ratio test in this branch, only for the "grow radius" gate below

        else

            call sparse_matvec(jac, p, jp)
            c_lin = c + jp
            h_lin = l1_violation(c_lin, problem%c_lb, problem%c_ub)
            pred  = q + linesearch%penalty*(h0 - h_lin) !! predicted decrease in the merit function

            call linesearch%eval_merit(f, c, problem%c_lb, problem%c_ub, new_lambda, phi0)
            call linesearch%eval_merit(f_trial, c_trial, problem%c_lb, problem%c_ub, new_lambda, phi_trial)
            ared = phi0 - phi_trial !! actual decrease in the merit function

            if (pred > 1.0e-12_wp) then
                ratio = ared/pred
            else
                ! the model doesn't even predict improvement: fall back to
                ! whether the true merit function improved anyway:
                ratio = merge(1.0_wp, -1.0_wp, ared > 0.0_wp)
            end if
            accept = ratio >= me%eta1

        end if

        if (accept) then

            if (use_filter) then
                mu = filter_penalty_estimate(new_lambda)
                call linesearch%filter_record(f_trial, h_trial, q, mu)
            end if

            if (maxval(abs(p)) >= 0.99_wp*me%radius) then
                ! the step used (approximately) the full trust region --
                ! grow it for next time, gated on ratio quality for the
                ! merit-ratio branch (Nocedal & Wright's textbook rule),
                ! or unconditionally for the filter branch (Fletcher &
                ! Leyffer's simpler "accepted implies possibly grow" rule):
                if (use_filter .or. ratio >= me%eta2) me%radius = min(me%radius*me%expand_factor, me%radius_max)
            end if

            x_new = x + p
            alpha = 1.0_wp
            istat = sqpopt_success
            return

        else

            me%radius = max(me%radius*me%shrink_factor, me%radius_min)
            if (me%radius <= me%radius_min) exit

        end if

    end do

    ! every retry was rejected (or the radius floor was reached): accept
    ! the last (smallest-radius) trial point anyway, mirroring the line-
    ! search modes' own "accept `alpha_min` anyway" fallback -- `x` is
    ! never left completely frozen, which would otherwise risk repeating
    ! the exact same failed retry sequence on the next major iteration:
    if (use_filter) then
        mu = filter_penalty_estimate(new_lambda)
        call linesearch%filter_record(f_trial, h_trial, q, mu)
    end if
    x_new = x + p
    alpha = 0.0_wp
    istat = sqpopt_line_search_failed

    end subroutine trust_region_step
!*******************************************************************************

    end module sqpopt_trust_region_module
!*******************************************************************************
