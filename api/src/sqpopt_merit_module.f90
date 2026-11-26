!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The merit functions, and their penalty parameter, used by the
!  merit-function line searches (`armijo`, `exact`, and `watchdog`, see
!  [[sqpopt_linesearch_module]]) and by the trust region's merit-ratio
!  test (see [[sqpopt_trust_region_module]]).
!
!  Two merit functions are available (`sqpopt_merit_type%mode`, set from
!  `options%merit_mode`):
!
!  * `sqpopt_merit_l1` (**default**) -- the standard non-smooth \( \ell_1 \)
!    exact penalty function (as in `slsqp`).
!  * `sqpopt_merit_augmented_lagrangian` -- a smooth augmented Lagrangian
!    merit function (Gill, Murray, Saunders & Wright, *"Some Theoretical
!    Properties of an Augmented Lagrangian Merit Function"*, SOL 86-6R --
!    the merit function used in NPSOL/NPSQP and, in spirit, SNOPT). Unlike
!    the \( \ell_1 \) function, it is twice continuously differentiable,
!    which is the reason SNOPT-family solvers do not need a second-order
!    correction (see [[sqpopt_soc_module]]) to avoid the Maratos effect.
!
!  Their penalty parameter is updated by one of two rules
!  (`penalty_update`, set from `options%penalty_update`, see
!  [[update_penalty_parameter]]): `sqpopt_penalty_multipliers`
!  (**default**; keeps it above the multiplier estimates, never
!  decreasing) or `sqpopt_penalty_model` (each merit's own principled
!  rule: Byrd-Nocedal model reduction for \( \ell_1 \), and for the
!  augmented Lagrangian the Gill-Murray-Saunders-Wright rule, with a joint
!  step in the variables, multipliers, and slacks, and a penalty that can
!  decrease).

    module sqpopt_merit_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, l1_violation
    use sqpopt_linalg_module,  only: sparse_matvec, sparse_matvec_transpose

    implicit none

    private

    integer, parameter, public :: sqpopt_merit_l1                   = 1  !! non-smooth \( \ell_1 \) exact penalty merit function (default)
    integer, parameter, public :: sqpopt_merit_augmented_lagrangian = 2  !! smooth augmented Lagrangian merit function (NPSOL/SNOPT-style)

    integer, parameter, public :: sqpopt_penalty_multipliers = 1  !! (default) penalty kept above the multipliers,
                                                                  !! \( \mu \ge \lVert\lambda\rVert_\infty + 1 \); never decreases
    integer, parameter, public :: sqpopt_penalty_model       = 2  !! each merit function's own principled rule (see
                                                                  !! [[update_penalty_parameter]]): Byrd-Nocedal
                                                                  !! model reduction (`sqpopt_merit_l1`), or
                                                                  !! Gill-Murray-Saunders-Wright (`sqpopt_merit_augmented_lagrangian`)

    type, public :: sqpopt_merit_type
        !! options and state of the merit function (the `merit` component of
        !! [[sqpopt_linesearch_type]]).

        integer  :: mode           = sqpopt_merit_l1            !! merit function to use (set from `options%merit_mode`)
        real(wp) :: penalty        = 1.0_wp                     !! current penalty parameter (called \( \mu \) for
                                                                 !! `sqpopt_merit_l1`, \( \rho \) for
                                                                 !! `sqpopt_merit_augmented_lagrangian`)
        integer  :: penalty_update = sqpopt_penalty_multipliers !! how the penalty parameter is updated (see the
                                                                 !! `sqpopt_penalty_*` constants; set from
                                                                 !! `options%penalty_update`)
        real(wp) :: penalty_rho    = 0.1_wp                     !! `sqpopt_penalty_model` with `sqpopt_merit_l1`: the
                                                                 !! fraction \( \rho \) of the linearized violation
                                                                 !! reduction the penalty must credit (Nocedal & Wright
                                                                 !! eq. 18.36), \( 0 < \rho < 1 \)

        ! internal state for `sqpopt_penalty_model` with `sqpopt_merit_augmented_lagrangian` (not user
        ! options): the joint step in the multipliers, and the floor limiting how often the penalty decreases
        logical  :: joint_active  = .false.  !! whether the merit's multipliers move along the step (see
                                             !! [[update_penalty_parameter]])
        real(wp), dimension(:), allocatable :: lambda0 !! the multipliers at the start of the step `dimension(m)`
        real(wp), dimension(:), allocatable :: s0      !! the slacks at the start of the step `dimension(m)`
        real(wp), dimension(:), allocatable :: q       !! the slacks' change per unit step `dimension(m)`
        real(wp) :: penalty_floor = 1.0e-2_wp !! the penalty can only decrease above this, which doubles each time

        contains

        procedure, public :: eval                   => eval_merit_function
        procedure, public :: directional_derivative => merit_directional_derivative
        procedure, public :: update_penalty         => update_penalty_parameter

    end type sqpopt_merit_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  evaluate the merit function used to measure progress, dispatching on
!  `me%mode`:
!
!  * `sqpopt_merit_l1`:
!    $$ \phi(x) = f(x) + \mu \lVert \max(c_l - c(x), 0, c(x) - c_u) \rVert_1 $$
!  * `sqpopt_merit_augmented_lagrangian`:
!    $$ \phi(x,\lambda,\rho) = f(x) - \lambda^T\!\left(c(x)-s\right) + \tfrac{1}{2}\rho \lVert c(x)-s \rVert_2^2 $$
!    where the slack `s` is the closed-form minimizer of \( \phi \) subject
!    to \( c_l \le s \le c_u \) (see [[augmented_lagrangian_slacks]]).
!    With the joint step (`joint_active`), the multipliers and the slacks
!    are those at step length `alpha` along it, \( \lambda_0 +
!    \alpha(\lambda - \lambda_0) \) and \( s_0 + \alpha q \) (`alpha=0` if
!    absent), instead of the slacks' closed-form minimizer (see
!    [[update_penalty_parameter]]).

    subroutine eval_merit_function(me, f, c, c_lb, c_ub, lambda, phi, alpha)

    class(sqpopt_merit_type), intent(inout) :: me
    real(wp),                intent(in)  :: f      !! objective function value
    real(wp), dimension(:), intent(in)  :: c      !! constraint values `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:), intent(in)  :: lambda !! Lagrange multiplier estimate `dimension(m)`
                                                  !! (only used by `sqpopt_merit_augmented_lagrangian`)
    real(wp),                intent(out) :: phi    !! value of the merit function
    real(wp), optional,      intent(in)  :: alpha  !! step length along the joint step (see above)

    real(wp), dimension(:), allocatable :: s, r, lam

    allocate(s(size(c)), r(size(c)), lam(size(c)))

    select case (me%mode)
    case (sqpopt_merit_augmented_lagrangian)
        if (me%joint_active) then
            ! (the multipliers and the slacks both move along the step)
            if (present(alpha)) then
                lam = me%lambda0 + alpha*(lambda - me%lambda0)
                s   = me%s0 + alpha*me%q
            else
                lam = me%lambda0
                s   = me%s0
            end if
        else
            lam = lambda
            call augmented_lagrangian_slacks(me, c, c_lb, c_ub, lam, s)
        end if
        r   = c - s
        phi = f - dot_product(lam, r) + 0.5_wp*me%penalty*dot_product(r, r)
    case default
        phi = f + me%penalty*sum(max(c_lb-c, 0.0_wp) + max(c-c_ub, 0.0_wp))
    end select

    end subroutine eval_merit_function
!*******************************************************************************
!*******************************************************************************
!>
!  the closed-form slack `s` that minimizes the augmented Lagrangian merit
!  function subject to \( c_l \le s \le c_u \): this generalizes the
!  non-negative slack \( s_i=\max(0,c_i-\lambda_i/\rho) \) of the original
!  (single-sided) formulation to `sqpopt`'s two-sided constraint bounds
!  (equality rows, where `c_lb=c_ub`, are handled automatically since `s`
!  is then clipped to that single value regardless of \( \lambda,\rho \)).

    subroutine augmented_lagrangian_slacks(me, c, c_lb, c_ub, lambda, s)

    class(sqpopt_merit_type), intent(in)  :: me
    real(wp), dimension(:), intent(in)  :: c !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_lb !! lower bounds on the constraints `dimension(m)`
    real(wp), dimension(:), intent(in)  :: c_ub !! upper bounds on the constraints `dimension(m)`
    real(wp), dimension(:), intent(in)  :: lambda !! Lagrange multipliers for the constraints `dimension(m)`
    real(wp), dimension(:), intent(out) :: s !! closed-form slack that minimizes the augmented Lagrangian merit function

    if (me%penalty <= 0.0_wp) then
        s = min(max(c, c_lb), c_ub)
    else
        s = min(max(c - lambda/me%penalty, c_lb), c_ub)
    end if

    end subroutine augmented_lagrangian_slacks
!*******************************************************************************
!*******************************************************************************
!>
!  the (approximate) directional derivative \( D(\phi;p) \) of the merit
!  function along `p`, used by the Armijo sufficient-decrease test and by
!  the descent-direction safeguard in [[sqpopt_iterate_module]]. Dispatches
!  on `me%mode`:
!
!  * `sqpopt_merit_l1`: \( D(\phi;p) = g^Tp + \mu \, D(\lVert \text{viol}(c + \alpha Jp)
!    \rVert_1; \alpha=0^+) \), the exact one-sided derivative of the
!    linearized violation. For a step that satisfies the linearized
!    constraints this is the usual \( g^Tp - \mu \lVert \text{viol}(x)
!    \rVert_1 \); for one that doesn't (e.g. shortened by the step-length
!    cap, or an elastic step), that formula would overstate the decrease,
!    so that no step length could pass the sufficient-decrease test.
!  * `sqpopt_merit_augmented_lagrangian`: \( D(\phi;p) = (g - J^T\lambda +
!    \rho J^T(c-s))^Tp \), holding \( \lambda \) and `s` fixed at their
!    current values; with the joint step (`joint_active`, see
!    [[update_penalty_parameter]]), where \( r = c-s \) changes by
!    \( d = Jp - q \) and \( \lambda \) by \( \xi \) per unit step, it is
!    \( g^Tp - \xi^Tr_0 - \lambda_0^Td + \rho\, r_0^Td \).

    subroutine merit_directional_derivative(me, jac, g, p, c, c_lb, c_ub, lambda, dphi0)

    class(sqpopt_merit_type), intent(in) :: me
    type(sqpopt_sparse_matrix), intent(in) :: jac  !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:), intent(in) :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(in) :: p        !! search direction `dimension(n)`
    real(wp), dimension(:), intent(in) :: c        !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_lb     !! lower bounds on the constraints `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_ub     !! upper bounds on the constraints `dimension(m)`
    real(wp), dimension(:), intent(in) :: lambda   !! Lagrange multipliers for the constraints `dimension(m)`
    real(wp), intent(out) :: dphi0 !! directional derivative of the merit function along `p`

    real(wp), dimension(:), allocatable :: s !! slack variables for the augmented Lagrangian
    real(wp), dimension(:), allocatable :: jtlam, jtr

    allocate(s(size(c)), jtlam(size(g)), jtr(size(g)))

    select case (me%mode)
    case (sqpopt_merit_augmented_lagrangian)
        if (me%joint_active) then
            ! (the multipliers and the slacks move too: by `xi = lambda - lambda0`
            ! and `q` per unit step, so `r = c-s` changes by `d = Jp - q`)
            block
                real(wp), dimension(:), allocatable :: jp, d
                allocate(jp(size(c)), d(size(c)))
                call sparse_matvec(jac, p, jp)
                d = jp - me%q
                s = c - me%s0     ! (= r0)
                dphi0 = dot_product(g, p) - dot_product(lambda - me%lambda0, s) - dot_product(me%lambda0, d) &
                        + me%penalty*dot_product(s, d)
            end block
        else
            call augmented_lagrangian_slacks(me, c, c_lb, c_ub, lambda, s)
            call sparse_matvec_transpose(jac, lambda, jtlam)
            s = c - s   ! (= r)
            call sparse_matvec_transpose(jac, s, jtr)
            dphi0 = dot_product(g - jtlam + me%penalty*jtr, p)
        end if
    case default
        block
            real(wp), dimension(:), allocatable :: jp
            real(wp) :: rate
            integer :: i
            allocate(jp(size(c)))
            call sparse_matvec(jac, p, jp)
            rate = 0.0_wp
            do i = 1, size(c)
                if (c(i) < c_lb(i)) then
                    rate = rate - jp(i)
                else if (c(i) > c_ub(i)) then
                    rate = rate + jp(i)
                else if (c(i) == c_lb(i) .and. jp(i) < 0.0_wp) then
                    rate = rate - jp(i)
                else if (c(i) == c_ub(i) .and. jp(i) > 0.0_wp) then
                    rate = rate + jp(i)
                end if
            end do
            dphi0 = dot_product(g, p) + me%penalty*rate
        end block
    end select

    end subroutine merit_directional_derivative
!*******************************************************************************
!*******************************************************************************
!>
!  update the merit function's penalty parameter after the QP solve, for
!  the step `p` (with \( p^THp \) = `php`), according to `penalty_update`:
!
!  * `sqpopt_penalty_multipliers`: \( \mu \ge \lVert\lambda_{QP}\rVert_\infty + 1 \)
!    (Han/Powell, as in `slsqp`), which never decreases. Simple, but at a
!    degenerate point, where the multipliers are huge, the penalty becomes
!    huge too, and the line search then only accepts tiny steps.
!  * `sqpopt_penalty_model`, with `sqpopt_merit_l1`: Byrd-Nocedal's
!    model-reduction rule (Nocedal & Wright, *Numerical Optimization*,
!    eq. 18.36): the penalty only increases if the step's predicted merit
!    reduction would not credit at least a fraction `penalty_rho` of the
!    linearized violation reduction \( \Delta v = v(c) - v(c+Jp) \):
!    $$ \mu \ge \frac{g^Tp + \tfrac12 \max(p^THp, 0)}{(1-\rho)\,\Delta v} $$
!    This depends on the step itself rather than on the multiplier
!    estimates, so it stays moderate at degenerate points.
!  * `sqpopt_penalty_model`, with `sqpopt_merit_augmented_lagrangian`:
!    Gill, Murray, Saunders & Wright (SOL 86-6R; NPSOL). The line search
!    moves the multipliers and the slacks too (`joint_active`): from the
!    current \( \lambda_0 \) toward the QP's \( \lambda_{QP} \)
!    (\( \xi = \lambda_{QP}-\lambda_0 \)), and from the merit's minimizer
!    \( s_0 \) toward the QP's linearized constraint values
!    (\( q = \text{clip}(c+Jp) - s_0 \)). The penalty is set so that the
!    merit's slope is at most \( -\tfrac12 p^THp \): with
!    \( r = c-s_0 \) and \( d = Jp-q \) (\( = -r \) for a consistent QP),
!    the slope is \( A + \rho B \), \( A = g^Tp - \xi^Tr - \lambda_0^Td \),
!    \( B = r^Td \), so (if \( B<0 \)) \( \hat\rho = (A + \tfrac12
!    p^THp)/(-B) \). If \( \rho < \hat\rho \), \( \rho \) increases to
!    \( \max(\hat\rho, 2\rho) \); if \( \rho > 4\max(\hat\rho, \rho_f) \), it
!    *decreases* to \( \max(\hat\rho, \rho_f, \sqrt{\rho\max(\hat\rho,\rho_f)}) \),
!    where the floor \( \rho_f \) doubles after each decrease, so that it
!    can only decrease finitely often (as the theory requires).

    subroutine update_penalty_parameter(me, jac, g, p, php, c, c_lb, c_ub, lambda, lambda_qp)

    class(sqpopt_merit_type), intent(inout) :: me
    type(sqpopt_sparse_matrix), intent(in) :: jac       !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:),     intent(in) :: g         !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:),     intent(in) :: p         !! the step `dimension(n)`
    real(wp),                   intent(in) :: php       !! \( p^THp \)
    real(wp), dimension(:),     intent(in) :: c         !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:),     intent(in) :: c_lb      !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),     intent(in) :: c_ub      !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:),     intent(in) :: lambda    !! the current multipliers `dimension(m)`
    real(wp), dimension(:),     intent(in) :: lambda_qp !! the QP's multipliers `dimension(m)`

    real(wp), dimension(:), allocatable :: jp, s, r
    real(wp) :: dv, req, a, b, rho_hat, target

    allocate(jp(size(c)), s(size(c)), r(size(c)))

    me%joint_active = .false.
    if (size(c) == 0) return

    select case (me%penalty_update)

    case (sqpopt_penalty_model)

        call sparse_matvec(jac, p, jp)

        select case (me%mode)

        case (sqpopt_merit_augmented_lagrangian)
            ! the joint step: the slacks start at the merit's minimizer (resetting
            ! them can only decrease the merit) and move toward the QP's
            ! linearized constraint values, so that `r = c-s` decreases like
            ! the violation (`d = Jp - q = -r0` if the QP is consistent):
            me%joint_active = .true.
            me%lambda0 = lambda
            call augmented_lagrangian_slacks(me, c, c_lb, c_ub, lambda, s)
            me%s0 = s
            me%q  = min(max(c + jp, c_lb), c_ub) - s
            r = c - s
            a = dot_product(g, p) - dot_product(lambda_qp - lambda, r) - dot_product(lambda, jp - me%q)
            b = dot_product(r, jp - me%q)
            if (b < -tiny(1.0_wp)) then
                rho_hat = (a + 0.5_wp*max(php, 0.0_wp))/(-b)
                target  = max(rho_hat, me%penalty_floor)
                if (me%penalty < rho_hat) then
                    me%penalty = max(rho_hat, 2.0_wp*me%penalty)
                else if (me%penalty > 4.0_wp*target) then
                    me%penalty = max(target, sqrt(me%penalty*target))
                    me%penalty_floor = 2.0_wp*me%penalty_floor
                end if
            end if

        case default   ! (l1)
            r  = c + jp   ! (the linearized constraints)
            dv = l1_violation(c, c_lb, c_ub) - l1_violation(r, c_lb, c_ub)
            if (dv > 0.0_wp) then
                req = (dot_product(g, p) + 0.5_wp*max(php, 0.0_wp))/((1.0_wp - me%penalty_rho)*dv)
                if (me%penalty < req) me%penalty = 1.1_wp*req
            end if

        end select

    case default   ! (sqpopt_penalty_multipliers)
        me%penalty = max(me%penalty, maxval(abs(lambda_qp)) + 1.0_wp)
    end select

    end subroutine update_penalty_parameter
!*******************************************************************************

    end module sqpopt_merit_module
!*******************************************************************************
