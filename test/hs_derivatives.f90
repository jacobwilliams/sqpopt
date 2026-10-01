!*******************************************************************************
!>
!  Derivative checking and finite-difference derivatives for the
!  Hock-Schittkowski test harnesses (`test_hs_suite.f90`, `test_hs_slsqp.f90`),
!  so that they make the same decisions: a problem's gradient or Jacobian
!  that doesn't match central differences (see [[check_derivatives]]) is
!  replaced by central differences (one-sided at a variable bound, see
!  [[fd_step]]).
!
!  The routines work on the problem in `hs_current`, which the harness sets
!  up (with `hs_setup`) before calling them.

    module hs_derivatives_module

    use sqpopt_kinds, only: dp => sqpopt_module_wp
    use hs_problems_module, only: hs_problem, hs_f, hs_g, hs_c, hs_jac
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    private

    type(hs_problem), public :: hs_current !! the problem being solved (its bounds are used by [[fd_step]])

    public :: check_derivatives, fd_gradient, fd_jacobian, fd_step

    contains

    subroutine check_derivatives(id, n, m, fd_g, fd_jac)
    !! compare the problem's gradient and Jacobian with finite differences at
    !! the starting point and at a nearby point, and switch to finite
    !! differences for whichever doesn't match. Each is compared at three
    !! difference steps (`h`, `h/100`, `10h`), and only counts as a mismatch
    !! if even the closest of the three disagrees (by more than 1e-3,
    !! relatively): finite differences of badly scaled functions can be
    !! inaccurate at any one step (e.g. TP376, 383), but a wrong derivative
    !! disagrees at all of them.
    integer, intent(in)    :: id, n, m      !! the problem (`hs_current` must be set up for it) and its size
    logical, intent(inout) :: fd_g, fd_jac  !! set if the gradient/Jacobian doesn't match
    real(dp), dimension(n) :: x, g
    real(dp), dimension(n,3) :: gd
    real(dp), dimension(m,n) :: jac
    real(dp), dimension(m,n,3) :: jd
    real(dp), dimension(3), parameter :: hf = [1.0_dp, 1.0e-2_dp, 1.0e1_dp]
    integer :: t, j_, k
    do t = 1, 2
        x = hs_current%x0
        if (t == 2) x = x + 0.01_dp*max(1.0_dp, abs(x))*[(sin(real(3*t+j_, dp)), j_=1,n)]
        x = min(max(x, hs_current%x_lb), hs_current%x_ub)
        call hs_g(id, x, g)
        do k = 1, 3
            call fd_gradient(id, x, gd(:,k), hf(k))
        end do
        if (mismatch(g, gd)) fd_g = .true.
        if (m > 0) then
            call hs_jac(id, x, jac)
            do k = 1, 3
                call fd_jacobian(id, x, jd(:,:,k), hf(k))
            end do
            if (mismatch(reshape(jac, [size(jac)]), reshape(jd, [size(jac), 3]))) fd_jac = .true.
        end if
    end do
    end subroutine check_derivatives

    logical function mismatch(a, b)
    !! whether analytic values `a` differ from the finite-difference values
    !! `b(:,k)` (at three steps `k`) by more than 1e-3 (relatively) at every
    !! step. A non-finite analytic value where a finite difference is finite
    !! is a mismatch (e.g. the original TP25's gradient took a fractional
    !! power of a negative number); non-finite finite differences are skipped.
    real(dp), dimension(:),   intent(in) :: a !! the analytic values
    real(dp), dimension(:,:), intent(in) :: b !! the finite-difference values, at three steps (columns)
    integer :: q, k
    real(dp) :: best
    mismatch = .false.
    do q = 1, size(a)
        best = huge(1.0_dp)
        do k = 1, size(b,2)
            if (.not. ieee_is_finite(b(q,k))) cycle
            if (.not. ieee_is_finite(a(q))) then
                best = huge(1.0_dp)
                exit
            end if
            best = min(best, abs(a(q) - b(q,k))/max(1.0_dp, abs(a(q)), abs(b(q,k))))
        end do
        if (all(.not. ieee_is_finite(b(q,:)))) cycle
        if (best > 1.0e-3_dp) mismatch = .true.
    end do
    end function mismatch

    subroutine fd_gradient(id, x, g, hfac, forward)
    !! finite-difference gradient of problem `id`'s objective (see [[fd_step]])
    integer,                intent(in)  :: id   !! problem number
    real(dp), dimension(:), intent(in)  :: x    !! point `dimension(n)`
    real(dp), dimension(:), intent(out) :: g    !! the finite-difference gradient `dimension(n)`
    real(dp), optional,     intent(in)  :: hfac !! factor on the default difference step
    logical,  optional,     intent(in)  :: forward !! use (first-order) forward differences instead (see
                                                   !! [[forward_step]]): half the function evaluations,
                                                   !! but less accurate
    real(dp), dimension(size(x)) :: xp
    real(dp) :: f0, f1, f2, h
    integer :: j, side
    call hs_f(id, x, f0)
    if (present(forward)) then
        if (forward) then
            do j = 1, size(x)
                h = forward_step(x, j)
                xp = x; xp(j) = x(j) + h;  call hs_f(id, xp, f1)
                g(j) = (f1 - f0)/h
            end do
            return
        end if
    end if
    do j = 1, size(x)
        call fd_step(x, j, h, side, hfac)
        xp = x; xp(j) = x(j) + h;          call hs_f(id, xp, f1)
        xp(j) = x(j) + merge(-h, 2*h, side == 0); call hs_f(id, xp, f2)
        if (side == 0) then
            g(j) = (f1 - f2)/(2.0_dp*h)
        else
            g(j) = (-3.0_dp*f0 + 4.0_dp*f1 - f2)/(2.0_dp*h)
        end if
    end do
    end subroutine fd_gradient

    subroutine fd_jacobian(id, x, jac, hfac, forward)
    !! finite-difference Jacobian of problem `id`'s constraints (see [[fd_step]])
    integer,                  intent(in)  :: id   !! problem number
    real(dp), dimension(:),   intent(in)  :: x    !! point `dimension(n)`
    real(dp), dimension(:,:), intent(out) :: jac  !! the finite-difference Jacobian `dimension(m,n)`
    real(dp), optional,       intent(in)  :: hfac !! factor on the default difference step
    logical,  optional,       intent(in)  :: forward !! use (first-order) forward differences instead (see
                                                     !! [[forward_step]])
    real(dp), dimension(size(x)) :: xp
    real(dp), dimension(size(jac,1)) :: c0, c1, c2
    real(dp) :: h
    integer :: j, side
    call hs_c(id, x, c0)
    if (present(forward)) then
        if (forward) then
            do j = 1, size(x)
                h = forward_step(x, j)
                xp = x; xp(j) = x(j) + h;  call hs_c(id, xp, c1)
                jac(:,j) = (c1 - c0)/h
            end do
            return
        end if
    end if
    do j = 1, size(x)
        call fd_step(x, j, h, side, hfac)
        xp = x; xp(j) = x(j) + h;          call hs_c(id, xp, c1)
        xp(j) = x(j) + merge(-h, 2*h, side == 0); call hs_c(id, xp, c2)
        if (side == 0) then
            jac(:,j) = (c1 - c2)/(2.0_dp*h)
        else
            jac(:,j) = (-3.0_dp*c0 + 4.0_dp*c1 - c2)/(2.0_dp*h)
        end if
    end do
    end subroutine fd_jacobian

    subroutine fd_step(x, j, h, side, hfac)
    !! the step for differencing along `x(j)`: central (`side=0`, points
    !! `x(j)+h` and `x(j)-h`), or, if that would cross one of the problem's
    !! variable bounds, a second-order one-sided difference into the
    !! interior (`side=1`, points `x(j)+h` and `x(j)+2h`, with `h<0` at an
    !! upper bound). Some problems' functions have a kink at a bound (e.g.
    !! TP358 clips `x` to its bounds), where a central difference would
    !! average across it (TP331, 358, 376, 383 start on such a bound).
    real(dp), dimension(:), intent(in)  :: x    !! point `dimension(n)`
    integer,                intent(in)  :: j    !! the variable to difference along
    real(dp),               intent(out) :: h    !! the step
    integer,                intent(out) :: side !! `0` central, `+1`/`-1` one-sided forward/backward
    real(dp), optional,     intent(in)  :: hfac !! factor on the default step (default 1)
    h = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp, abs(x(j)))
    if (present(hfac)) h = hfac*h
    side = 0
    if (x(j) - h < hs_current%x_lb(j)) then
        side = 1
    else if (x(j) + h > hs_current%x_ub(j)) then
        side = 1
        h = -h
    end if
    end subroutine fd_step

    real(dp) function forward_step(x, j) result(h)
    !! the step for a forward difference along `x(j)`: `sqrt(eps)` relative,
    !! backward instead if a forward step would cross the upper bound
    real(dp), dimension(:), intent(in) :: x !! point `dimension(n)`
    integer,                intent(in) :: j !! the variable to difference along
    h = sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(x(j)))
    if (x(j) + h > hs_current%x_ub(j)) h = -h
    end function forward_step

    end module hs_derivatives_module
!*******************************************************************************
