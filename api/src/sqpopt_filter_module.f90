!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The filter (Fletcher & Leyffer, *"Nonlinear programming without a
!  penalty function"*, Math. Program. 91 (2002)), with the globally
!  convergent rules of Wächter & Biegler (*"Line search filter methods for
!  nonlinear programming: motivation and global convergence"*, SIAM J.
!  Optim. 16 (2005), and the IPOPT paper, Math. Program. 106 (2006)): the
!  acceptance test of the filter line search (see
!  [[sqpopt_linesearch_module]]'s `filter_line_search`) and of the trust
!  region with the filter (see [[sqpopt_trust_region_module]]), and the
!  filter itself, which persists across major iterations.
!
!  A trial point is judged by the pair \( (\theta, \varphi) \) of
!  \( \ell_1 \) constraint violation and objective value, and must not be
!  dominated by the filter. A *switching condition* decides whether a step
!  must reduce the objective (an Armijo test, "f-type" step) or may trade
!  objective for feasibility; only non-f-type steps enlarge the filter
!  (see [[filter_step_acceptable]]).

    module sqpopt_filter_module

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: merit_slack

    implicit none

    private

    type, public :: sqpopt_filter_type
        !! options and state of the filter (the `filter` component of
        !! [[sqpopt_linesearch_type]]), with the default values from
        !! Wächter & Biegler.

        real(wp) :: gamma_theta    = 1.0e-5_wp !! margin \( \gamma_\theta \): a step must reduce \( \theta \) by this fraction...
        real(wp) :: gamma_phi      = 1.0e-5_wp !! ...or reduce \( \varphi \) by \( \gamma_\varphi\theta \)
        real(wp) :: delta          = 1.0_wp    !! switching condition constant \( \delta \)
        real(wp) :: s_theta        = 1.1_wp    !! switching condition exponent \( s_\theta \)
        real(wp) :: s_phi          = 2.3_wp    !! switching condition exponent \( s_\varphi \)
        real(wp) :: eta_phi        = 1.0e-4_wp !! Armijo constant \( \eta_\varphi \) for f-type steps
        real(wp) :: theta_max_fact = 1.0e4_wp  !! \( \theta_{max} = \) this \( \times \max(1,\theta_0) \): no point with a larger violation is accepted
        real(wp) :: theta_min_fact = 1.0e-4_wp !! \( \theta_{min} = \) this \( \times \max(1,\theta_0) \): below it, f-type steps are allowed
        real(wp) :: gamma_alpha    = 0.05_wp   !! safety factor \( \gamma_\alpha \) in the minimum step length before restoration

        ! internal state (not user options -- persists across major iterations):
        logical  :: ready = .false.        !! whether the filter below has been initialized
        real(wp) :: theta_max = 0.0_wp     !! \( \theta_{max} \)
        real(wp) :: theta_min = 0.0_wp     !! \( \theta_{min} \)
        real(wp), dimension(:), allocatable :: theta !! constraint-violation value of each filter entry
        real(wp), dimension(:), allocatable :: phi   !! objective value of each filter entry

        contains

        procedure, public :: prepare    => filter_prepare_state
        procedure, public :: accept     => filter_step_acceptable
        procedure, public :: record     => filter_augment
        procedure, public :: min_step   => filter_min_step
        procedure, public :: acceptable => filter_point_acceptable

    end type sqpopt_filter_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the filter (empty) and its bounds \( \theta_{max} \),
!  \( \theta_{min} \) from the constraint violation `theta0` at the
!  starting point, the first time it is used; a no-op after that. Shared
!  by [[filter_line_search]] and the trust region's filter acceptance, so
!  both use the same filter.

    subroutine filter_prepare_state(me, theta0)

    class(sqpopt_filter_type), intent(inout) :: me
    real(wp),                      intent(in)    :: theta0 !! constraint violation at the starting point

    if (me%ready) return

    me%theta_max = me%theta_max_fact*max(1.0_wp, theta0)
    me%theta_min = me%theta_min_fact*max(1.0_wp, theta0)
    if (allocated(me%theta)) deallocate(me%theta)
    if (allocated(me%phi))   deallocate(me%phi)
    allocate(me%theta(0), me%phi(0))
    me%ready = .true.

    end subroutine filter_prepare_state
!*******************************************************************************
!*******************************************************************************
!>
!  whether a trial point with violation `theta_t` and objective `phi_t`,
!  reached with step length `alpha` from the current point
!  (`theta_k`, `phi_k`), is acceptable (see [[filter_line_search]] for the
!  rules). `gtp` is the predicted change in the objective along the full
!  step (\( g^Tp \) for the line search; the trust region passes the
!  quadratic model's \( -q \), with `alpha=1`). `f_type` is set if the
!  switching condition held (so the step was judged on the objective alone,
!  and must not enlarge the filter). The trial point must also be
!  acceptable to the filter itself: \( \theta_t < \theta_{max} \), and not
!  dominated by any filter entry (\( \theta_t < \theta_j \) or
!  \( \varphi_t < \varphi_j \) for every entry `j`).

    function filter_step_acceptable(me, theta_k, phi_k, gtp, alpha, theta_t, phi_t, f_type) result(ok)

    class(sqpopt_filter_type), intent(in) :: me
    real(wp), intent(in)  :: theta_k, phi_k   !! violation and objective at the current point
    real(wp), intent(in)  :: gtp              !! predicted change in the objective along the full step
    real(wp), intent(in)  :: alpha            !! step length
    real(wp), intent(in)  :: theta_t, phi_t   !! violation and objective at the trial point
    logical,  intent(out) :: f_type           !! true if the switching condition held
    logical :: ok

    integer :: j

    ! switching condition, alpha*(-gtp)^s_phi > delta*theta_k^s_theta
    ! (compared in log space, so the powers can't overflow):
    f_type = gtp < 0.0_wp .and. theta_k <= me%theta_min
    if (f_type .and. theta_k > 0.0_wp) then
        f_type = log(alpha) + me%s_phi*log(-gtp) > log(me%delta) + me%s_theta*log(theta_k)
    end if

    if (f_type) then
        ! Armijo condition on the objective (with a roundoff-level slack):
        ok = phi_t <= phi_k + me%eta_phi*alpha*gtp + merit_slack(phi_k)
    else
        ! sufficient reduction of the violation or the objective:
        ok = theta_t <= (1.0_wp - me%gamma_theta)*theta_k .or. &
             phi_t   <= phi_k - me%gamma_phi*theta_k
    end if
    if (.not. ok) return

    ! acceptable to the filter:
    ok = theta_t < me%theta_max
    if (.not. ok) return
    do j = 1, size(me%theta)
        if (theta_t >= me%theta(j) .and. phi_t >= me%phi(j)) then
            ok = .false.
            return
        end if
    end do

    end function filter_step_acceptable
!*******************************************************************************
!*******************************************************************************
!>
!  augment the filter with the current point (`theta_k`, `phi_k`), with
!  margins: the entry \( ((1-\gamma_\theta)\theta_k,
!  \varphi_k-\gamma_\varphi\theta_k) \) is added, and any existing entries
!  it dominates are removed.

    subroutine filter_augment(me, theta_k, phi_k)

    class(sqpopt_filter_type), intent(inout) :: me
    real(wp),                      intent(in)    :: theta_k, phi_k !! violation and objective at the current point

    real(wp) :: theta_new, phi_new
    logical, dimension(:), allocatable :: keep

    if (.not. me%ready) call me%prepare(theta_k)

    theta_new = (1.0_wp - me%gamma_theta)*theta_k
    phi_new   = phi_k - me%gamma_phi*theta_k

    if (size(me%theta) > 0) then
        keep = .not. (me%theta >= theta_new .and. me%phi >= phi_new)
        me%theta = pack(me%theta, keep)
        me%phi   = pack(me%phi,   keep)
    end if
    me%theta = [me%theta, theta_new]
    me%phi   = [me%phi,   phi_new]

    end subroutine filter_augment
!*******************************************************************************
!*******************************************************************************
!>
!  Wächter & Biegler's minimum step length \( \alpha_{min} \): below it,
!  no trial step can be acceptable (so the line search should give up and
!  go to feasibility restoration).

    pure function filter_min_step(me, theta_k, gtp) result(alpha_min)

    class(sqpopt_filter_type), intent(in) :: me
    real(wp), intent(in) :: theta_k !! violation at the current point
    real(wp), intent(in) :: gtp     !! \( g^Tp \)
    real(wp) :: alpha_min

    if (gtp < 0.0_wp) then
        alpha_min = min(me%gamma_theta, me%gamma_phi*theta_k/(-gtp))
        if (theta_k <= me%theta_min) then
            if (theta_k > 0.0_wp) then
                alpha_min = min(alpha_min, exp(log(me%delta) + me%s_theta*log(theta_k) &
                                               - me%s_phi*log(-gtp)))
            else
                alpha_min = 0.0_wp
            end if
        end if
    else
        alpha_min = me%gamma_theta
    end if
    alpha_min = me%gamma_alpha*alpha_min

    end function filter_min_step
!*******************************************************************************

!*******************************************************************************
!>
!  whether a point with violation `theta` and objective `phi` is acceptable
!  to the filter: \( \theta < \theta_{max} \), and not dominated by any
!  filter entry (always true before the filter is initialized).

    function filter_point_acceptable(me, theta, phi) result(ok)

    class(sqpopt_filter_type), intent(in) :: me
    real(wp), intent(in) :: theta !! constraint violation
    real(wp), intent(in) :: phi   !! objective
    logical :: ok

    ok = .true.
    if (.not. me%ready) return
    ok = theta < me%theta_max
    if (ok .and. size(me%theta) > 0) ok = .not. any(theta >= me%theta .and. phi >= me%phi)

    end function filter_point_acceptable
!*******************************************************************************

    end module sqpopt_filter_module
!*******************************************************************************
