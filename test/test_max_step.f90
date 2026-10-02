program test_max_step

    !! The limit on the change of each variable in one major iteration
    !! (`problem%set_max_step`). The Rosenbrock function on the unit disk (see
    !! `test_rosenbrock_disk`),
    !!
    !!     minimize   100 (y - x^2)^2 + (1 - x)^2
    !!     subject to x^2 + y^2 <= 1
    !!
    !! is solved from (-1.2, 1) with the limits (0.2, 0.05) on the steps of
    !! `x` and `y`, with every globalization: the filter and funnel line
    !! searches, a merit-function line search, the watchdog, and the trust
    !! region, and with the exact Hessian. Every solve must converge to the
    !! solution, (0.78642, 0.61770), and no iterate may differ from the
    !! previous one by more than the limits (the iterates are recorded by the
    !! `report` callback). Without limits the first step of the default solve
    !! is 2.5 long, so the limits must also make the solve take more
    !! iterations.
    !!
    !! A limit of `sqpopt_infinity` is no limit: with (`sqpopt_infinity`,
    !! 0.05) the steps of `x` are free and those of `y` limited. A problem
    !! without constraints is limited too (its QPs would otherwise be solved
    !! by the unconstrained step of the QP front end). And invalid limits
    !! (the wrong size, or not positive) must be reported as invalid input.

    use sqpopt_module,              only: sqpopt_type
    use sqpopt_problem_module,      only: sqpopt_problem_type
    use sqpopt_options_module,      only: sqpopt_options_type
    use sqpopt_hessian_module,      only: sqpopt_hessian_exact
    use sqpopt_linesearch_module,   only: sqpopt_linesearch_funnel, sqpopt_linesearch_armijo, sqpopt_linesearch_watchdog
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_types_module,        only: sqpopt_results_type, sqpopt_success, sqpopt_acceptable, sqpopt_stalled, &
                                          sqpopt_invalid_input, sqpopt_infinity
    use sqpopt_kinds,               only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: x_star(2) = [0.7864152_wp, 0.6176983_wp] !! the solution on the disk
    real(wp), parameter :: tol = 1.0e-5_wp                          !! tolerance on it

    integer, parameter :: max_path = 2000
    real(wp) :: path(2, max_path) !! the iterates
    integer  :: n_path            !! number of iterates recorded
    logical  :: constrained       !! whether the constraint is part of the problem being solved

    integer :: iter_free, iter_limited

    write(*,*) '----------------------------'
    write(*,*) 'test_max_step'
    write(*,*) '----------------------------'

    constrained = .true.
    call run('no limits', [sqpopt_infinity, sqpopt_infinity], x_star, iter_free)
    call run('filter (default)', [0.2_wp, 0.05_wp], x_star, iter_limited)
    if (iter_limited <= iter_free) error stop 'test_max_step FAILED: the limits did not lengthen the solve'
    call run('funnel', [0.2_wp, 0.05_wp], x_star, iter_limited)
    call run('armijo', [0.2_wp, 0.05_wp], x_star, iter_limited)
    call run('watchdog', [0.2_wp, 0.05_wp], x_star, iter_limited)
    call run('trust region', [0.2_wp, 0.05_wp], x_star, iter_limited)
    call run('exact Hessian', [0.2_wp, 0.05_wp], x_star, iter_limited)
    call run('y only', [sqpopt_infinity, 0.05_wp], x_star, iter_limited)

    ! without the constraint: the minimizer is (1, 1)
    constrained = .false.
    call run('unconstrained', [0.2_wp, 0.05_wp], [1.0_wp, 1.0_wp], iter_limited)

    call invalid('wrong size', [0.2_wp])
    call invalid('zero', [0.2_wp, 0.0_wp])
    call invalid('negative', [-1.0_wp, 0.05_wp])

    print '(A)', 'test_max_step PASSED'

    contains

    subroutine setup(problem, exact)
    !! define the problem
    type(sqpopt_problem_type), intent(out) :: problem !! the problem
    logical,                   intent(in)  :: exact   !! whether to supply the Hessian
    integer, dimension(0) :: none
    if (constrained) then
        call problem%set_problem_size(n=2, m=1)
        call problem%set_bounds(x_lb=[-1.0e20_wp, -1.0e20_wp], x_ub=[1.0e20_wp, 1.0e20_wp], &
                                c_lb=[-1.0e20_wp], c_ub=[1.0_wp])
        call problem%set_jacobian_sparsity(nnz=2, irow=[1, 1], icol=[1, 2])
    else
        call problem%set_problem_size(n=2, m=0)
        call problem%set_bounds(x_lb=[-1.0e20_wp, -1.0e20_wp], x_ub=[1.0e20_wp, 1.0e20_wp], &
                                c_lb=[real(wp) ::], c_ub=[real(wp) ::])
        call problem%set_jacobian_sparsity(nnz=0, irow=none, icol=none)
    end if
    if (exact) then
        call problem%set_functions(fc=fc, gjac=gjac, hess=hess)
        call problem%set_hessian_sparsity(3, [1, 2, 2], [1, 1, 2])
    else
        call problem%set_functions(fc=fc, gjac=gjac)
    end if
    end subroutine setup

    subroutine run(label, max_step, x_expected, iterations)
    !! solve the problem with the step limits `max_step`, and check the iterates and the solution
    character(len=*),       intent(in)  :: label      !! the configuration
    real(wp), dimension(:), intent(in)  :: max_step   !! the limits `dimension(2)`
    real(wp), dimension(:), intent(in)  :: x_expected !! the solution `dimension(2)`
    integer,                intent(out) :: iterations !! the major iterations taken
    type(sqpopt_type)              :: solver
    type(sqpopt_problem_type)      :: problem
    type(sqpopt_options_type)      :: options
    type(sqpopt_trust_region_type) :: trust_region
    type(sqpopt_results_type)      :: r
    real(wp) :: x(2), step(2), largest(2)
    real(wp), dimension(:), allocatable :: lambda
    integer  :: istat, k

    call setup(problem, exact=label == 'exact Hessian')
    if (label /= 'no limits') call problem%set_max_step(max_step)
    options%max_iter = 1000
    select case (label)
    case ('funnel');        options%linesearch_mode = sqpopt_linesearch_funnel
    case ('armijo');        options%linesearch_mode = sqpopt_linesearch_armijo
    case ('watchdog');      options%linesearch_mode = sqpopt_linesearch_watchdog
    case ('trust region');  trust_region%enabled = .true.
    case ('exact Hessian'); options%hessian_mode = sqpopt_hessian_exact
    case default
    end select

    n_path = 0
    call solver%initialize(problem=problem, options=options, trust_region=trust_region, report=report)
    call solver%solve([-1.2_wp, 1.0_wp], istat)
    allocate(lambda(problem%m))
    call solver%get_solution(x, lambda)
    call solver%get_results(r)
    iterations = r%iterations

    largest = 0.0_wp
    do k = 2, n_path
        step = abs(path(:, k) - path(:, k-1))
        largest = max(largest, step)
    end do
    print '(A18,A,I0,A,I0,A,2F10.6,A,2F9.5)', label, ': istat = ', istat, ', iterations = ', iterations, &
        ', x = ', x, ', largest steps = ', largest

    if (istat /= sqpopt_success .and. istat /= sqpopt_acceptable .and. istat /= sqpopt_stalled) then
        error stop 'test_max_step FAILED: '//label//': the solve did not converge'
    end if
    if (maxval(abs(x - x_expected)) > tol) error stop 'test_max_step FAILED: '//label//': wrong solution'
    if (n_path < 2) error stop 'test_max_step FAILED: the iterates were not reported'
    if (label /= 'no limits') then
        if (any(largest > max_step*(1.0_wp + 1.0e-12_wp))) then
            error stop 'test_max_step FAILED: '//label//': a step is longer than its limit'
        end if
    else if (.not. largest(1) > 1.0_wp) then
        error stop 'test_max_step FAILED: the unlimited solve has no long step (the test needs one)'
    end if

    end subroutine run

    subroutine invalid(label, max_step)
    !! the limits `max_step` must be reported as invalid input
    character(len=*),       intent(in) :: label    !! the case
    real(wp), dimension(:), intent(in) :: max_step !! the limits
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    integer :: istat
    call setup(problem, exact=.false.)
    call problem%set_max_step(max_step)
    call solver%initialize(problem=problem)
    call solver%solve([-1.2_wp, 1.0_wp], istat)
    print '(A18,A,I0,2A)', label, ': istat = ', istat, '  ', solver%status_message()
    if (istat /= sqpopt_invalid_input) error stop 'test_max_step FAILED: invalid limits were accepted'
    end subroutine invalid

    subroutine fc(x, f, c, status, data)
    !! the objective and the constraint
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(2)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint value at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (none)
    f = 100.0_wp*(x(2) - x(1)**2)**2 + (1.0_wp - x(1))**2
    if (size(c) > 0) c(1) = x(1)**2 + x(2)**2
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient and the constraint Jacobian
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(2)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    g = [-400.0_wp*x(1)*(x(2) - x(1)**2) - 2.0_wp*(1.0_wp - x(1)), 200.0_wp*(x(2) - x(1)**2)]
    if (size(jac_val) > 0) jac_val = 2.0_wp*x
    end subroutine gjac

    subroutine hess(x, lambda, hess_val, status, data)
    !! the Hessian of the Lagrangian `f - lambda*c`, lower triangle: (1,1), (2,1), (2,2)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(2)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multiplier `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x`
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    real(wp) :: lam
    lam = 0.0_wp
    if (size(lambda) > 0) lam = lambda(1)
    hess_val = [1200.0_wp*x(1)**2 - 400.0_wp*x(2) + 2.0_wp - 2.0_wp*lam, -400.0_wp*x(1), 200.0_wp - 2.0_wp*lam]
    end subroutine hess

    subroutine report(iter, x, f, c, lambda, user_stop, data)
    !! record the iterate
    integer,                intent(in)    :: iter      !! major iteration number (starts at 1)
    real(wp), dimension(:), intent(in)    :: x         !! current point `dimension(2)`
    real(wp),               intent(in)    :: f         !! current objective value
    real(wp), dimension(:), intent(in)    :: c         !! current constraint value `dimension(m)`
    real(wp), dimension(:), intent(in)    :: lambda    !! current multiplier estimate `dimension(m)`
    logical,                intent(out)   :: user_stop !! set `.true.` to request the solver stop
    class(*), optional,     intent(inout) :: data      !! the user data given to `set_functions` (none)
    user_stop = .false.
    if (n_path < max_path) then
        n_path = n_path + 1
        path(:, n_path) = x
    end if
    end subroutine report

end program test_max_step
