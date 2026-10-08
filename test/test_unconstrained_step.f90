program test_unconstrained_step

    !! The unconstrained step of the QP front end (`qp_solver%unconstrained_step`,
    !! see [[sqpopt_qp_solver_module]]): with the limited-memory BFGS Hessian, a
    !! QP whose unconstrained minimizer \( -H^{-1} g \) satisfies the bounds
    !! and the linearized constraints is solved by that step, without a QP
    !! solver. The problems are
    !!
    !!     minimize   sum_i (x_i - 2)^2 + sum_{i<n} (x_{i+1} - x_i)^2/2
    !!     subject to x_lb <= x <= x_ub,   c_lb <= sum_i x_i <= c_ub
    !!
    !! from \( x = 0 \), with 40 variables (the dense QP solver) and 300 (the
    !! sparse one), and the bounds of five cases:
    !!
    !! * no bounds, and a constraint that is never active: every QP must be
    !!   solved by the unconstrained step, and the solution, \( x = 2 \), must
    !!   be the one found with the step turned off, by the QP solver;
    !! * the same, with the limited-memory SR1 Hessian, whose step is the
    !!   minimizer within the step cap, from its spectral decomposition (see
    !!   [[sqpopt_spectral_module]]): the problem is a convex quadratic, so the
    !!   matrix stays positive definite, and every QP must again be solved by
    !!   the step;
    !! * upper bounds of 1 on the odd variables, which are active at the
    !!   solution: the step can't solve the QPs near the solution, and the
    !!   solution must again be the one found without it;
    !! * an active inequality constraint (`sum x <= n`): likewise, with the
    !!   solution \( x = 1 \) and a multiplier of \( -2 \);
    !! * an equality constraint (`sum x = n`): the step can only be used
    !!   by a QP whose unconstrained minimizer happens to satisfy it, and
    !!   the solution must be the same.

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_sr1
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type
    use sqpopt_types_module,     only: sqpopt_results_type, sqpopt_success, sqpopt_acceptable, sqpopt_stalled
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp   !! "no bound"
    real(wp), parameter :: tol = 1.0e-4_wp   !! tolerance on the solution

    integer, parameter :: sizes(2) = [40, 300]
    integer :: k

    write(*,*) '----------------------------'
    write(*,*) 'test_unconstrained_step'
    write(*,*) '----------------------------'

    do k = 1, size(sizes)
        call test('free', sizes(k))
        call test('sr1', sizes(k))
        call test('bounds', sizes(k))
        call test('inequality', sizes(k))
        call test('equality', sizes(k))
    end do

    print '(A)', 'test_unconstrained_step PASSED'

    contains

    subroutine test(name, n)
    !! solve one case with and without the unconstrained step, and check the results
    character(len=*), intent(in) :: name !! the case (see the program's documentation)
    integer,          intent(in) :: n    !! number of variables
    real(wp), dimension(n) :: x_on, x_off, x_star
    real(wp), dimension(1) :: lambda_on, lambda_off
    type(sqpopt_results_type) :: r_on, r_off
    integer :: i

    call run(name, n, .true., x_on, lambda_on, r_on)
    call run(name, n, .false., x_off, lambda_off, r_off)
    print '(A12,A,I4,5(A,I0),A,ES9.2)', name, ': n =', n, ', QPs ', r_on%n_qp_solves, ', by the unconstrained step ', &
        r_on%n_unconstrained_qp, ', iterations ', r_on%iterations, ' (', r_off%iterations, &
        ' without it), QP iterations ', r_on%n_qp_iterations, ', |x_on - x_off| = ', maxval(abs(x_on - x_off))

    if (r_off%n_unconstrained_qp /= 0) error stop 'test_unconstrained_step FAILED: the step was used when turned off'
    if (maxval(abs(x_on - x_off)) > tol) error stop 'test_unconstrained_step FAILED: '//name//': a different solution'
    if (abs(lambda_on(1) - lambda_off(1)) > tol) then
        error stop 'test_unconstrained_step FAILED: '//name//': a different multiplier'
    end if

    select case (name)
    case ('free')
        x_star = 2.0_wp
        if (r_on%n_unconstrained_qp /= r_on%n_qp_solves .or. r_on%n_qp_solves == 0) then
            error stop 'test_unconstrained_step FAILED: free: not every QP was solved by the step'
        end if
        if (r_on%n_qp_iterations /= 0) error stop 'test_unconstrained_step FAILED: free: a QP solver was run'
        if (lambda_on(1) /= 0.0_wp) error stop 'test_unconstrained_step FAILED: free: a nonzero multiplier'
    case ('sr1')
        x_star = 2.0_wp
        if (r_on%n_unconstrained_qp /= r_on%n_qp_solves .or. r_on%n_qp_solves == 0) then
            error stop 'test_unconstrained_step FAILED: sr1: not every QP was solved by the step'
        end if
        if (r_on%n_qp_iterations /= 0) error stop 'test_unconstrained_step FAILED: sr1: a QP solver was run'
    case ('bounds')
        ! (the solution of the bound-constrained problem isn't known in closed form: compare with
        ! the run without the step, above, and check the bounds)
        x_star = x_on
        if (any(x_on(1:n:2) > 1.0_wp)) error stop 'test_unconstrained_step FAILED: bounds: a bound is violated'
        if (.not. all(x_on(1:n:2) > 1.0_wp - tol)) error stop 'test_unconstrained_step FAILED: bounds: not active'
        if (r_on%n_unconstrained_qp >= r_on%n_qp_solves) then
            error stop 'test_unconstrained_step FAILED: bounds: every QP was solved by the step'
        end if
    case ('inequality')
        x_star = 1.0_wp
        if (abs(lambda_on(1) + 2.0_wp) > tol) error stop 'test_unconstrained_step FAILED: inequality: wrong multiplier'
        if (r_on%n_unconstrained_qp >= r_on%n_qp_solves) then
            error stop 'test_unconstrained_step FAILED: inequality: every QP was solved by the step'
        end if
    case default
        x_star = 1.0_wp
        if (abs(lambda_on(1) + 2.0_wp) > tol) error stop 'test_unconstrained_step FAILED: equality: wrong multiplier'
    end select
    do i = 1, n
        if (abs(x_on(i) - x_star(i)) > tol) error stop 'test_unconstrained_step FAILED: '//name//': wrong solution'
    end do

    end subroutine test

    subroutine run(name, n, step, x, lambda, r)
    !! solve one case
    character(len=*),          intent(in)  :: name   !! the case (see the program's documentation)
    integer,                   intent(in)  :: n      !! number of variables
    logical,                   intent(in)  :: step   !! `qp_solver%unconstrained_step`
    real(wp), dimension(:),    intent(out) :: x      !! the solution `dimension(n)`
    real(wp), dimension(:),    intent(out) :: lambda !! the constraint's multiplier `dimension(1)`
    type(sqpopt_results_type), intent(out) :: r      !! the results

    type(sqpopt_type)           :: solver
    type(sqpopt_problem_type)   :: problem
    type(sqpopt_options_type)   :: options
    type(sqpopt_qp_solver_type) :: qp
    real(wp), dimension(n) :: x_lb, x_ub
    real(wp) :: c_lb, c_ub
    integer :: i, istat

    x_lb = -big
    x_ub = big
    c_lb = -big
    c_ub = 10.0_wp*real(n, wp)
    select case (name)
    case ('sr1');        options%hessian_mode = sqpopt_hessian_sr1
    case ('bounds');     x_ub(1:n:2) = 1.0_wp
    case ('inequality'); c_ub = real(n, wp)
    case ('equality');   c_lb = real(n, wp); c_ub = real(n, wp)
    case default
    end select

    call problem%set_problem_size(n=n, m=1)
    call problem%set_bounds(x_lb, x_ub, [c_lb], [c_ub])
    call problem%set_jacobian_sparsity(n, [(1, i=1,n)], [(i, i=1,n)])
    call problem%set_functions(fc=fc, gjac=gjac)

    options%max_iter = 500
    qp%unconstrained_step = step
    call solver%initialize(problem=problem, options=options, qp_solver=qp)
    call solver%solve([(0.0_wp, i=1,n)], istat)
    call solver%get_solution(x, lambda)
    call solver%get_results(r)
    if (istat /= sqpopt_success .and. istat /= sqpopt_acceptable .and. istat /= sqpopt_stalled) then
        print '(A,A,I0,A,L1,A,I0)', name, ': n = ', n, ', step = ', step, ', istat = ', istat
        error stop 'test_unconstrained_step FAILED: a solve did not converge'
    end if

    end subroutine run

    subroutine fc(x, f, c, status, data)
    !! the objective and the constraint
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(1)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (none)
    integer :: n
    n = size(x)
    f = sum((x - 2.0_wp)**2) + 0.5_wp*sum((x(2:n) - x(1:n-1))**2)
    c(1) = sum(x)
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient and the constraint Jacobian
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` `dimension(n)`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    integer :: n
    n = size(x)
    g = 2.0_wp*(x - 2.0_wp)
    g(2:n)   = g(2:n)   + (x(2:n) - x(1:n-1))
    g(1:n-1) = g(1:n-1) - (x(2:n) - x(1:n-1))
    jac_val = 1.0_wp
    end subroutine gjac

end program test_unconstrained_step
