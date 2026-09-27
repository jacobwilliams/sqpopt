program test_acceptance

    !! Unit test of the acceptance tests of the globalization, each on its
    !! own (without the solver): [[sqpopt_filter_module]],
    !! [[sqpopt_funnel_module]], and [[sqpopt_merit_module]].
    !!
    !! * Filter: `theta_max`/`theta_min` from the starting violation; the
    !!   sufficient-reduction test of a non-f-type step, and the Armijo test
    !!   of an f-type one (switching condition); an entry added by `record`
    !!   rejects the points it dominates, and is removed when a later entry
    !!   dominates it; the minimum step length is positive.
    !! * Funnel: the initial width; the h-type test (inside `beta*width`) and
    !!   the f-type Armijo test; the width after an h-type step (both update
    !!   rules) and before a restoration step.
    !! * Merit: the l1 and augmented Lagrangian values; the l1 directional
    !!   derivative against a finite difference of the linearized merit; the
    !!   `sqpopt_penalty_multipliers` update.
    !! * The line search's dispatcher `globalization_acceptable` uses the
    !!   filter or the funnel according to the mode.

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_types_module,      only: sqpopt_sparse_matrix, l1_violation
    use sqpopt_filter_module,     only: sqpopt_filter_type
    use sqpopt_funnel_module,     only: sqpopt_funnel_type
    use sqpopt_merit_module,      only: sqpopt_merit_type, sqpopt_merit_l1, sqpopt_merit_augmented_lagrangian
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_filter, sqpopt_linesearch_funnel, &
                                        sqpopt_linesearch_armijo

    implicit none

    real(wp), parameter :: tol = 1.0e-12_wp

    write(*,*) '----------------------------'
    write(*,*) 'test_acceptance'
    write(*,*) '----------------------------'

    call test_filter()
    call test_funnel()
    call test_merit()
    call test_dispatch()

    write(*,*) 'test_acceptance PASSED'

contains

    subroutine check(ok, what)
    logical,          intent(in) :: ok
    character(len=*), intent(in) :: what
    if (.not. ok) error stop 'test_acceptance FAILED: '//what
    end subroutine check

    subroutine test_filter()
    type(sqpopt_filter_type) :: flt
    logical :: f_type

    call flt%prepare(2.0_wp)
    call check(flt%ready, 'filter not initialized')
    call check(abs(flt%theta_max - 2.0e4_wp) < tol .and. abs(flt%theta_min - 2.0e-4_wp) < tol, 'filter theta_max/min')
    call flt%prepare(100.0_wp)   ! (a no-op once initialized)
    call check(abs(flt%theta_max - 2.0e4_wp) < tol, 'filter re-initialized')

    ! non-f-type (theta_k > theta_min): sufficient reduction of the violation or the objective
    call check(flt%accept(1.0_wp, 5.0_wp, -1.0_wp, 1.0_wp, 0.5_wp, 6.0_wp, f_type), 'filter: violation reduction rejected')
    call check(.not. f_type, 'filter: f-type with theta_k > theta_min')
    call check(flt%accept(1.0_wp, 5.0_wp, -1.0_wp, 1.0_wp, 1.5_wp, 4.0_wp, f_type), 'filter: objective reduction rejected')
    call check(.not. flt%accept(1.0_wp, 5.0_wp, -1.0_wp, 1.0_wp, 1.0_wp, 5.0_wp, f_type), 'filter: no reduction accepted')
    call check(.not. flt%accept(1.0_wp, 5.0_wp, -1.0_wp, 1.0_wp, 3.0e4_wp, -1.0e10_wp, f_type), &
               'filter: point above theta_max accepted')

    ! f-type (theta_k = 0 <= theta_min, g^Tp < 0): the Armijo test on the objective
    call check(flt%accept(0.0_wp, 5.0_wp, -1.0_wp, 1.0_wp, 0.0_wp, 4.0_wp, f_type) .and. f_type, 'filter: f-type step')
    call check(.not. flt%accept(0.0_wp, 5.0_wp, -1.0_wp, 1.0_wp, 0.0_wp, 5.0_wp, f_type), 'filter: f-type without decrease')

    ! an entry rejects the points it dominates:
    call flt%record(1.0_wp, 5.0_wp)
    call check(size(flt%theta) == 1, 'filter: entry not added')
    call check(.not. flt%acceptable(2.0_wp, 6.0_wp), 'filter: dominated point acceptable')
    call check(flt%acceptable(0.5_wp, 10.0_wp) .and. flt%acceptable(2.0_wp, 1.0_wp), 'filter: non-dominated point rejected')
    call check(.not. flt%accept(1.0_wp, 5.0_wp, -1.0_wp, 1.0_wp, 1.0_wp, 5.0_wp, f_type), 'filter: dominated trial accepted')

    ! ... and is removed when a new entry dominates it:
    call flt%record(0.5_wp, 3.0_wp)
    call check(size(flt%theta) == 1 .and. abs(flt%theta(1) - 0.5_wp*(1.0_wp - flt%gamma_theta)) < tol, &
               'filter: dominated entry not removed')
    call flt%record(2.0_wp, 1.0_wp)
    call check(size(flt%theta) == 2, 'filter: non-dominated entry not kept')

    call check(flt%min_step(1.0_wp, -1.0_wp) > 0.0_wp, 'filter: minimum step length')
    print '(A)', 'filter: OK'
    end subroutine test_filter

    subroutine test_funnel()
    type(sqpopt_funnel_type) :: fun
    logical :: f_type
    real(wp) :: w

    call check(fun%acceptable(1.0e10_wp), 'funnel: point rejected before initialization')
    call fun%prepare(2.0_wp)
    call check(abs(fun%width - 3.0_wp) < tol, 'funnel: initial width')   ! max(width_min, width_fact*theta0)
    call check(fun%acceptable(3.0_wp) .and. .not. fun%acceptable(3.1_wp), 'funnel: acceptable')

    ! h-type (no predicted decrease): sufficiently inside the funnel
    call check(fun%accept(2.0_wp, 5.0_wp, 0.0_wp, 2.5_wp, 9.0_wp, f_type) .and. .not. f_type, 'funnel: h-type step')
    call check(.not. fun%accept(2.0_wp, 5.0_wp, 0.0_wp, 3.5_wp, 1.0_wp, f_type), 'funnel: point outside accepted')

    ! f-type (switching condition): the Armijo test on the objective
    call check(fun%accept(1.0e-3_wp, 5.0_wp, 1.0_wp, 1.0e-3_wp, 4.0_wp, f_type) .and. f_type, 'funnel: f-type step')
    call check(.not. fun%accept(1.0e-3_wp, 5.0_wp, 1.0_wp, 1.0e-3_wp, 5.0_wp, f_type), 'funnel: f-type without decrease')

    ! shrinking, rule 1: max(beta*width, kappa*theta_k + (1-kappa)*theta_t)
    w = fun%width
    call fun%record(2.0_wp, 1.0_wp)
    call check(abs(fun%width - max(fun%beta*w, 1.5_wp)) < tol, 'funnel: update rule 1')
    ! rule 2: kappa*width + (1-kappa)*theta_t
    fun%update = 2
    w = fun%width
    call fun%record(2.0_wp, 1.0_wp)
    call check(abs(fun%width - (fun%kappa*w + (1.0_wp - fun%kappa)*1.0_wp)) < tol, 'funnel: update rule 2')
    ! before a restoration step: toward the current violation
    w = fun%width
    call fun%restoration(0.5_wp)
    call check(abs(fun%width - (fun%kappa*w + (1.0_wp - fun%kappa)*0.5_wp)) < tol, 'funnel: restoration')
    print '(A)', 'funnel: OK'
    end subroutine test_funnel

    subroutine test_merit()
    type(sqpopt_merit_type) :: mer
    type(sqpopt_sparse_matrix) :: jac
    real(wp), dimension(2) :: c, c_lb, c_ub, lambda, lambda_qp
    real(wp), dimension(3) :: g, p
    real(wp) :: phi, phi_h, dphi, h

    ! c1 = 3 violates c1 <= 2; c2 = 0 is at its lower bound 0
    c = [3.0_wp, 0.0_wp];  c_lb = [-1.0_wp, 0.0_wp];  c_ub = [2.0_wp, 5.0_wp]
    lambda = 0.0_wp
    mer%mode = sqpopt_merit_l1
    mer%penalty = 10.0_wp
    call mer%eval(1.0_wp, c, c_lb, c_ub, lambda, phi)
    call check(abs(phi - (1.0_wp + 10.0_wp*1.0_wp)) < tol, 'merit: l1 value')

    ! the directional derivative against a one-sided difference of the linearized merit
    jac%nrows = 2; jac%ncols = 3; jac%nnz = 4
    jac%irow = [1, 1, 2, 2];  jac%icol = [1, 2, 2, 3];  jac%val = [1.0_wp, -2.0_wp, 1.0_wp, -1.0_wp]
    g = [1.0_wp, 0.5_wp, -1.0_wp];  p = [-1.0_wp, 0.5_wp, 0.25_wp]
    call mer%directional_derivative(jac, g, p, c, c_lb, c_ub, lambda, dphi)
    h = 1.0e-7_wp
    phi_h = 1.0_wp + h*dot_product(g, p) + mer%penalty*l1_violation(c + h*[p(1)-2.0_wp*p(2), p(2)-p(3)], c_lb, c_ub)
    call check(abs((phi_h - phi)/h - dphi) < 1.0e-6_wp, 'merit: l1 directional derivative')

    ! augmented Lagrangian, lambda = 0: the slacks are clip(c), so phi = f + rho/2*|c - clip(c)|^2
    mer%mode = sqpopt_merit_augmented_lagrangian
    call mer%eval(1.0_wp, c, c_lb, c_ub, lambda, phi)
    call check(abs(phi - (1.0_wp + 0.5_wp*10.0_wp*1.0_wp)) < tol, 'merit: augmented Lagrangian value')

    ! penalty update (multipliers rule): at least max|lambda_qp| + 1, never decreasing
    lambda_qp = [-20.0_wp, 3.0_wp]
    call mer%update_penalty(jac, g, p, 1.0_wp, c, c_lb, c_ub, lambda, lambda_qp)
    call check(abs(mer%penalty - 21.0_wp) < tol, 'merit: penalty increase')
    lambda_qp = [1.0_wp, 1.0_wp]
    call mer%update_penalty(jac, g, p, 1.0_wp, c, c_lb, c_ub, lambda, lambda_qp)
    call check(abs(mer%penalty - 21.0_wp) < tol, 'merit: penalty decreased')
    print '(A)', 'merit: OK'
    end subroutine test_merit

    subroutine test_dispatch()
    type(sqpopt_linesearch_type) :: ls
    call ls%filter%record(1.0_wp, 5.0_wp)
    call ls%funnel%prepare(0.1_wp)            ! (width 1)
    ls%mode = sqpopt_linesearch_filter
    call check(.not. ls%globalization_acceptable(2.0_wp, 6.0_wp), 'dispatch: filter')
    ls%mode = sqpopt_linesearch_funnel
    call check(ls%globalization_acceptable(0.5_wp, 100.0_wp) .and. .not. ls%globalization_acceptable(2.0_wp, 0.0_wp), &
               'dispatch: funnel')
    ls%mode = sqpopt_linesearch_armijo
    call check(ls%globalization_acceptable(1.0e10_wp, 1.0e10_wp), 'dispatch: armijo')
    print '(A)', 'dispatch: OK'
    end subroutine test_dispatch

end program test_acceptance
