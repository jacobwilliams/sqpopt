!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The funnel (Kiessling, Leyffer & Vanaret, *"A unified funnel
!  restoration SQP algorithm"*, Math. Program. (2025); as implemented in
!  the Uno solver): the acceptance test of the funnel line search (see
!  [[sqpopt_linesearch_module]]'s `funnel_line_search`) and of the trust
!  region with the funnel (see [[sqpopt_trust_region_module]]), and the
!  funnel width, which persists across major iterations.
!
!  Like the filter (see [[sqpopt_filter_module]]), there is no merit
!  function or penalty parameter, but the filter's list of points is
!  replaced by a single number, the funnel width \( \tau \), a bound on the
!  \( \ell_1 \) constraint violation that shrinks as the iterations
!  progress. A trial point must be inside the funnel; a switching
!  condition then decides whether it must reduce the objective (an Armijo
!  test, "f-type") or sufficiently reduce the violation relative to the
!  funnel ("h-type", which shrinks the funnel). See
!  [[funnel_step_acceptable]].

    module sqpopt_funnel_module

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: merit_slack

    implicit none

    private

    type, public :: sqpopt_funnel_type
        !! options and state of the funnel (the `funnel` component of
        !! [[sqpopt_linesearch_type]]), with the default values of the Uno
        !! solver.

        real(wp) :: width_min  = 1.0_wp    !! initial funnel width \( \tau_0 = \max( \) this, `width_fact`
                                           !! \( \times\,\theta_0) \)
        real(wp) :: width_fact = 1.5_wp    !! see `width_min`
        real(wp) :: beta       = 0.9999_wp !! an h-type step must reach \( \theta \le \beta\tau \) (\( 0<\beta<1 \))
        real(wp) :: kappa      = 0.5_wp    !! after an h-type step, the new width is (at least) the convex
                                           !! combination \( \kappa\theta_k + (1-\kappa)\theta_t \) (see `update`)
        integer  :: update     = 1         !! funnel width update after an h-type step: `1` =
                                           !! \( \max(\beta\tau, \kappa\theta_k+(1-\kappa)\theta_t) \) if the
                                           !! violation decreased, else \( \beta\tau \); `2` = \( \kappa\tau +
                                           !! (1-\kappa)\theta_t \)
        real(wp) :: delta      = 0.999_wp  !! switching condition constant \( \delta \): an f-type step needs
                                           !! \( \alpha(-g^Tp) > \delta\theta_k^{s_\theta} \)
        real(wp) :: s_theta    = 2.0_wp    !! switching condition exponent \( s_\theta \)
        real(wp) :: eta        = 1.0e-4_wp !! Armijo constant \( \eta \) for f-type steps
        logical  :: require_current = .false. !! also require every trial point to be acceptable with respect to
                                              !! the current point: \( \theta_t < \beta\theta_k \) or
                                              !! \( \varphi_t \le \varphi_k - \gamma\theta_t \)
        real(wp) :: gamma      = 1.0e-3_wp !! \( \gamma \) in `require_current`

        ! internal state (not user options -- persists across major iterations):
        logical  :: ready = .false.  !! whether the funnel width below has been initialized
        real(wp) :: width = 0.0_wp   !! the funnel width \( \tau \)

        contains

        procedure, public :: prepare     => funnel_prepare_state
        procedure, public :: accept      => funnel_step_acceptable
        procedure, public :: record      => funnel_shrink
        procedure, public :: restoration => funnel_shrink_restoration
        procedure, public :: acceptable  => funnel_point_acceptable

    end type sqpopt_funnel_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the funnel width from the constraint violation `theta0` at
!  the starting point, the first time it is used:
!  \( \tau_0 = \max(\) `width_min`, `width_fact`
!  \( \times\,\theta_0) \). After that, it only makes sure the current
!  point is inside the funnel, \( \tau \ge \theta_k \) (a restoration or
!  escape step reduces a different measure of the violation, so it can
!  occasionally end slightly outside). Shared by [[funnel_line_search]] and
!  the trust region's funnel acceptance.

    subroutine funnel_prepare_state(me, theta_k)

    class(sqpopt_funnel_type), intent(inout) :: me
    real(wp),                      intent(in)    :: theta_k !! constraint violation at the current point

    if (.not. me%ready) then
        me%width = max(me%width_min, me%width_fact*theta_k)
        me%ready = .true.
    end if
    me%width = max(me%width, theta_k)

    end subroutine funnel_prepare_state
!*******************************************************************************
!*******************************************************************************
!>
!  whether a trial point with violation `theta_t` and objective `phi_t` is
!  acceptable to the funnel, from the current point (`theta_k`, `phi_k`)
!  (see [[funnel_line_search]] for the rules). `pred` is the predicted
!  decrease in the objective (\( -\alpha g^Tp \) for the line search; the
!  trust region passes the quadratic model's decrease `q`). `f_type` is set
!  if the switching condition held (so the step was judged on the objective
!  alone, and must not shrink the funnel). Doesn't change the funnel: the
!  caller shrinks it after accepting an h-type step (see [[funnel_shrink]]).

    function funnel_step_acceptable(me, theta_k, phi_k, pred, theta_t, phi_t, f_type) result(ok)

    class(sqpopt_funnel_type), intent(in) :: me
    real(wp), intent(in)  :: theta_k, phi_k   !! violation and objective at the current point
    real(wp), intent(in)  :: pred             !! predicted decrease in the objective
    real(wp), intent(in)  :: theta_t, phi_t   !! violation and objective at the trial point
    logical,  intent(out) :: f_type           !! true if the switching condition held
    logical :: ok

    f_type = .false.

    ! inside the funnel:
    ok = theta_t <= me%width
    if (.not. ok) return

    ! (optionally) acceptable with respect to the current point:
    if (me%require_current) then
        ok = theta_t < me%beta*theta_k .or. phi_t <= phi_k - me%gamma*theta_t
        if (.not. ok) return
    end if

    ! switching condition, pred > delta*theta_k^s_theta (compared in log
    ! space, so the power can't underflow or overflow):
    f_type = pred > 0.0_wp
    if (f_type .and. theta_k > 0.0_wp) &
        f_type = log(pred) > log(me%delta) + me%s_theta*log(theta_k)

    if (f_type) then
        ! Armijo condition on the objective (with a roundoff-level slack):
        ok = phi_t <= phi_k - me%eta*pred + merit_slack(phi_k)
    else
        ! h-type: sufficiently inside the funnel:
        ok = theta_t <= me%beta*me%width
    end if

    end function funnel_step_acceptable
!*******************************************************************************
!*******************************************************************************
!>
!  shrink the funnel after an accepted h-type step from violation `theta_k`
!  to `theta_t` (see `update` for the two rules).

    subroutine funnel_shrink(me, theta_k, theta_t)

    class(sqpopt_funnel_type), intent(inout) :: me
    real(wp),                      intent(in)    :: theta_k !! violation at the current point
    real(wp),                      intent(in)    :: theta_t !! violation at the accepted point

    if (me%update == 2) then
        me%width = me%kappa*me%width + (1.0_wp - me%kappa)*theta_t
    else if (theta_t <= theta_k) then
        me%width = max(me%beta*me%width, &
                              me%kappa*theta_k + (1.0_wp - me%kappa)*theta_t)
    else
        me%width = me%beta*me%width
    end if

    end subroutine funnel_shrink
!*******************************************************************************
!*******************************************************************************
!>
!  shrink the funnel toward the current violation `theta_k` before a
!  feasibility restoration step, so the iterations can't cycle back to the
!  current point: \( \tau = \kappa\tau + (1-\kappa)\theta_k \).

    subroutine funnel_shrink_restoration(me, theta_k)

    class(sqpopt_funnel_type), intent(inout) :: me
    real(wp),                      intent(in)    :: theta_k !! violation at the current point

    if (.not. me%ready) call me%prepare(theta_k)
    me%width = me%kappa*me%width + (1.0_wp - me%kappa)*theta_k

    end subroutine funnel_shrink_restoration
!*******************************************************************************

!*******************************************************************************
!>
!  whether a point with violation `theta` is inside the funnel,
!  \( \theta \le \tau \) (always true before the funnel is initialized).

    function funnel_point_acceptable(me, theta) result(ok)

    class(sqpopt_funnel_type), intent(in) :: me
    real(wp), intent(in) :: theta !! constraint violation
    logical :: ok

    ok = .true.
    if (me%ready) ok = theta <= me%width

    end function funnel_point_acceptable
!*******************************************************************************

    end module sqpopt_funnel_module
!*******************************************************************************
