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
!  The resulting step `p` is tested for acceptance (and, if it is rejected
!  without reducing the constraint violation, so is its second-order
!  correction, see [[sqpopt_soc_module]]), and the radius is grown or
!  shrunk accordingly:
!
!  * if `linesearch%mode == sqpopt_linesearch_filter`, acceptance uses the
!    same filter test as the filter line search (`linesearch%filter%accept`,
!    with the quadratic model's predicted decrease `q` in place of
!    \( g^Tp \) in the switching condition) -- a trust-region filter-SQP
!    method in the spirit of Fletcher & Leyffer (`references/fletcher.pdf`).
!  * if `linesearch%mode == sqpopt_linesearch_funnel`, acceptance uses the
!    funnel test of the funnel line search (`linesearch%funnel%accept`,
!    with `q` as the predicted decrease) -- the trust-region funnel SQP
!    method of Kiessling, Leyffer & Vanaret (Uno's `funnelsqp` preset).
!  * otherwise, acceptance uses the classical trust-region-SQP ratio test
!    (Nocedal & Wright, *Numerical Optimization*, Ch. 18): \( \rho =
!    \text{ared}/\text{pred} \), the ratio of the actual to the
!    predicted decrease in the merit function selected by
!    `linesearch%merit%mode`, accepted when \( \rho \ge \eta_1 \).
!
!  If a step is rejected, the radius is shrunk (`shrink_factor`) and the
!  QP is re-solved (bounded by `max_retries` retries per major iteration)
!  -- there is no `alpha` to backtrack along, since a genuinely different
!  `p` is computed each time. Trial points where `f` or `c` is not finite
!  are always rejected. If every retry is rejected, no step is taken
!  (`istat=sqpopt_line_search_failed`), as for a failed line search.
!
!  Because this bypasses `linesearch%search` entirely, the line-search-
!  specific options (`major_step_limit`, `alpha_min`/`sigma`/`backtrack`/
!  `max_ls_iter`, and the watchdog fields) are not used when trust region
!  is enabled -- there is no `alpha`, so nothing to cap or backtrack.
!
!  At `options%print_level >= 3`, every QP re-solve (radius, QP status, and
!  the trial point's ratio, or violation and objective, and whether it was
!  accepted) is written to the detailed log (`log`, see
!  [[sqpopt_log_module]]).

    module sqpopt_trust_region_module

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_types_module,      only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_line_search_failed, sqpopt_all_finite
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_hessian_module,    only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_filter, sqpopt_linesearch_funnel, &
                                         l1_violation
    use sqpopt_linalg_module,     only: sparse_matvec
    use sqpopt_soc_module,        only: soc_step
    use sqpopt_kkt_module,           only: sqpopt_kkt_type
    use sqpopt_inertia_module,       only: sqpopt_inertia_type
    use sqpopt_least_squares_module, only: sqpopt_least_squares_type
    use sqpopt_log_module,        only: sqpopt_log_type, sqpopt_log_detail, fmt_e, fmt_g, qp_status_text

    implicit none

    private

    type, public :: sqpopt_trust_region_type
        !! options and state for trust-region-based globalization
        !! (an alternative to `linesearch%search`, see module docs).

        logical  :: enabled       = .false.    !! if true, use trust-region radius management instead of a line search
        real(wp) :: radius0       = 1.0_wp      !! initial trust-region radius
        real(wp) :: radius_min    = 1.0e-8_wp   !! below this, a major iteration's retries give up (no step is taken)
        real(wp) :: radius_max    = 1.0e3_wp    !! ceiling on the radius
        real(wp) :: eta1          = 0.1_wp      !! ratio threshold to accept a step (merit-ratio acceptance only)
        real(wp) :: eta2          = 0.75_wp     !! ratio threshold to also grow the radius (merit-ratio acceptance only)
        real(wp) :: shrink_factor = 0.5_wp      !! `radius *= shrink_factor` on a rejected step
        real(wp) :: expand_factor = 2.0_wp      !! `radius *= expand_factor` on a step that used the full radius and was accepted
        integer  :: max_retries   = 20          !! maximum QP re-solves (with a shrinking radius) per major iteration

        ! internal state (persists across major iterations, like the watchdog/filter state on `sqpopt_linesearch_type`):
        logical  :: ready  = .false. !! whether `radius` has been initialized from `radius0` yet
        real(wp) :: radius = 0.0_wp  !! current trust-region radius

        ! the detailed log (set by `solve`), and whether the last step was second-order corrected (output):
        type(sqpopt_log_type) :: log
        logical :: used_soc = .false.

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
                                  x, g, f, c, jac, x_new, new_lambda, alpha, istat, kkt, inertia, least_squares)

    class(sqpopt_trust_region_type), intent(inout) :: me
    type(sqpopt_problem_type),      intent(inout) :: problem    !! problem definition
    type(sqpopt_hessian_type),      intent(inout) :: hessian    !! matrix-free Hessian approximation
    type(sqpopt_qp_solver_type),    intent(inout) :: qp_solver  !! QP subproblem solver
    type(sqpopt_linesearch_type),   intent(inout) :: linesearch !! supplies the merit function (`merit`, ratio test) or
                                                                  !! the filter (`filter`, `mode==sqpopt_linesearch_filter`)
                                                                  !! or funnel (`funnel`, `mode==sqpopt_linesearch_funnel`)
    real(wp), dimension(:),         intent(in)  :: x       !! current point `dimension(n)`
    real(wp), dimension(:),         intent(in)  :: g       !! objective gradient at `x` `dimension(n)`
    real(wp),                       intent(in)  :: f       !! objective value at `x`
    real(wp), dimension(:),         intent(in)  :: c       !! constraint values at `x` `dimension(m)`
    type(sqpopt_sparse_matrix),     intent(in)  :: jac     !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:),         intent(out) :: x_new      !! new point `dimension(n)` (`x` if no retry was accepted)
    real(wp), dimension(:),         intent(out) :: new_lambda !! new multipliers `dimension(m)`, from the QP solve that produced `x_new`
    real(wp),                       intent(out) :: alpha      !! `1` if a step was accepted, `0` otherwise (informational only -- there
                                                               !! is no line-search step length under trust-region globalization)
    integer,                        intent(out) :: istat      !! status code (see [[sqpopt_types_module]])
    type(sqpopt_kkt_type), optional, intent(inout) :: kkt     !! the KKT matrix, for the QP solver's direct method
                                                              !! (see [[solve_qp_subproblem]])
    type(sqpopt_inertia_type), optional, intent(inout) :: inertia !! the inertia control, for the same
    type(sqpopt_least_squares_type), optional, intent(inout) :: least_squares !! the direct least-squares solver, for
                                                              !! the second-order correction (see [[soc_step]])

    integer :: retry, qp_istat
    logical :: use_filter, use_funnel, accept, ok, soc_ok, f_type
    real(wp), dimension(:), allocatable :: p, p_soc, x_lb2, x_ub2, hp
    real(wp), dimension(:), allocatable :: jp, c_trial, c_lin
    real(wp) :: f_trial, h0, h_trial, q, pred, ratio, phi0, phi_model

    allocate(p(size(x)), p_soc(size(x)), x_lb2(size(x)), x_ub2(size(x)), hp(size(x)))
    allocate(jp(size(c)), c_trial(size(c)), c_lin(size(c)))

    if (.not. me%ready) then
        me%radius = me%radius0
        me%ready  = .true.
    end if

    me%used_soc = .false.
    use_funnel = (linesearch%mode == sqpopt_linesearch_funnel)
    use_filter = (linesearch%mode == sqpopt_linesearch_filter) .or. use_funnel !! (no merit function in either)
    h0         = l1_violation(c, problem%c_lb, problem%c_ub)
    if (use_funnel) then
        call linesearch%funnel%prepare(h0)
    else if (use_filter) then
        call linesearch%filter%prepare(h0)
    end if

    do retry = 1, me%max_retries

        ! (the box, within the variables' bounds and their step limits, if any)
        call problem%step_bounds(x, x_lb2, x_ub2)
        x_lb2 = max(x_lb2, x - me%radius)
        x_ub2 = min(x_ub2, x + me%radius)
        ! the trust region controls the step length, so the QP solver's own
        ! (line-search) cap on it, `max_step*step_scale`, must not bind inside
        ! the box (whose steps are at most `sqrt(n)*radius` long): otherwise
        ! every step is cut to `max_step`, never uses the full radius, and the
        ! radius can't grow either (e.g. TP220, which starts 25000 away from
        ! its solution, crawled 2 per iteration)
        qp_solver%step_scale = max(qp_solver%step_scale, 1.01_wp*sqrt(real(size(x), wp))*me%radius/qp_solver%max_step)

        call qp_solver%solve(hessian, jac, x, g, c, x_lb2, x_ub2, problem%c_lb, problem%c_ub, p, new_lambda, qp_istat, &
                             kkt=kkt, inertia=inertia)

        ! keep the merit function's penalty parameter dominating the current
        ! multiplier estimates, same rule as the line-search path (needed for
        ! `merit%eval`'s ared/pred to be meaningful in the merit-ratio test):
        if (size(new_lambda) > 0) linesearch%merit%penalty = max(linesearch%merit%penalty, maxval(abs(new_lambda)) + 1.0_wp)

        ! the model's predicted decrease in `f`, and in the merit function: the
        ! merit function (whichever `merit_mode` is selected) evaluated at the
        ! QP model's objective `f-q` and linearized constraints `c+Jp` (for
        ! `sqpopt_merit_l1` this is the classical `q + penalty*(h(c)-h(c+Jp))`):
        call hessian%hv_product(p, hp)
        q = -(dot_product(g, p) + 0.5_wp*dot_product(p, hp))
        if (.not. use_filter) then
            call sparse_matvec(jac, p, jp)
            c_lin = c + jp
            call linesearch%merit%eval(f, c, problem%c_lb, problem%c_ub, new_lambda, phi0)
            call linesearch%merit%eval(f-q, c_lin, problem%c_lb, problem%c_ub, new_lambda, phi_model)
            pred = phi0 - phi_model
        end if

        x_new = x + p   ! (the trial point)
        call evaluate(x_new, accept)
        call log_trial('tr  radius '//fmt_e(me%radius)//', QP '//qp_status_text(qp_istat)//': ')

        if (.not. accept .and. ok .and. problem%m > 0) then
            if (h_trial >= h0) then
                ! the step didn't reduce the constraint violation, so its
                ! rejection may be due to constraint curvature (the Maratos
                ! effect): try the second-order-corrected step (kept inside
                ! the trust region) before shrinking the radius:
                call soc_step(jac, x, p, c, c_trial, problem%c_lb, problem%c_ub, x_lb2, x_ub2, p_soc, soc_ok, &
                              least_squares=least_squares)
                if (soc_ok) then
                    x_new = x + p_soc
                    call evaluate(x_new, accept)
                    call log_trial('tr  second-order correction: ')
                    if (accept) then
                        p = p_soc
                        me%used_soc = .true.
                    end if
                end if
            end if
        end if

        if (accept) then

            ! (a step that isn't f-type adds the current point to the filter,
            ! or shrinks the funnel)
            if (use_funnel .and. .not. f_type) then
                call linesearch%funnel%record(h0, h_trial)
            else if (use_filter .and. .not. f_type) then
                call linesearch%filter%record(h0, f)
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

    ! every retry was rejected (or the radius floor was reached): no step
    ! is taken. The next major iteration starts from the (small) radius
    ! reached here, with a reset Hessian approximation (see
    ! [[sqpopt_iterate_module]]), so it won't repeat the same retries:
    call me%log%put(sqpopt_log_detail, 'tr  no acceptable step (radius '//fmt_e(me%radius)//')')
    x_new = x
    alpha = 0.0_wp
    istat = sqpopt_line_search_failed

    contains

        subroutine log_trial(prefix)
        !! the detailed log's line for the last evaluated trial point
        character(len=*), intent(in) :: prefix !! the start of the line (what was tried)
        character(len=:), allocatable :: line
        if (.not. me%log%on(sqpopt_log_detail)) return
        if (.not. ok) then
            line = prefix//'non-finite function value, rejected'
        else if (use_filter) then
            line = prefix//'violation '//fmt_e(h_trial)//', objective '//fmt_g(f_trial)
        else
            line = prefix//'ratio '//fmt_e(ratio)//' (actual/predicted merit decrease)'
        end if
        if (ok) then
            if (accept) then
                line = line//', accepted'
            else
                line = line//', rejected'
            end if
        end if
        call me%log%put(sqpopt_log_detail, line)
        end subroutine log_trial


        subroutine evaluate(x_trial, accept)
        !! evaluate the trial point `x_trial` (setting `f_trial`, `c_trial`,
        !! `h_trial`, `ok`, and `ratio`), and whether it is acceptable. A
        !! trial point where `f` or `c` is not finite is never acceptable.
        real(wp), dimension(:), intent(in)  :: x_trial !! the trial point `dimension(n)`
        logical,                intent(out) :: accept  !! whether it is acceptable
        real(wp) :: phi_trial, ared

        call problem%f(x_trial, f_trial)
        call problem%c(x_trial, c_trial)
        ok = sqpopt_all_finite([f_trial]) .and. sqpopt_all_finite(c_trial)
        accept = .false.
        f_type = .false.
        ratio  = -1.0_wp
        if (.not. ok) return
        h_trial = l1_violation(c_trial, problem%c_lb, problem%c_ub)

        if (use_filter) then

            if (use_funnel) then
                accept = linesearch%funnel%accept(h0, f, q, h_trial, f_trial, f_type)
            else
                accept = linesearch%filter%accept(h0, f, -q, 1.0_wp, h_trial, f_trial, f_type)
            end if
            ratio = 1.0_wp !! not used for the ratio test in this branch, only for the "grow radius" gate

        else

            call linesearch%merit%eval(f_trial, c_trial, problem%c_lb, problem%c_ub, new_lambda, phi_trial)
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

        end subroutine evaluate

    end subroutine trust_region_step
!*******************************************************************************

    end module sqpopt_trust_region_module
!*******************************************************************************
