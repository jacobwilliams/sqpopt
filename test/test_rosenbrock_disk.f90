program test_rosenbrock_disk

    !! A problem with two variables and one constraint, small enough to draw:
    !! the Rosenbrock function on the unit disk,
    !!
    !!     minimize   100 (y - x^2)^2 + (1 - x)^2
    !!     subject to x^2 + y^2 <= 1
    !!
    !! from the usual starting point of the Rosenbrock function, (-1.2, 1),
    !! which is outside the disk. The unconstrained minimizer, (1, 1), is
    !! outside it too, so the constraint is active at the solution,
    !! (0.78642, 0.61770), where the objective is 0.045675 and the
    !! constraint's multiplier is -0.12150 (negative: an upper bound).
    !!
    !! The solve must converge to that point, with the default options and
    !! with the exact Hessian.
    !!
    !! With the option `--path=FILE`, the iterates of the default solve (the
    !! iteration, `x`, `y`, the objective, the constraint's value, and the
    !! solver's KKT error there, one line each) are also written to `FILE`,
    !! for the figure of the guide's Performance page (see
    !! `tools/rosenbrock_disk_figure.py`). The KKT error is the one of the
    !! solver's convergence test (of the scaled problem), from the per-iteration
    !! table of its diagnostics (`diagnostic_level = 2`, which doesn't change
    !! the iterates):
    !!
    !!     fpm test test_rosenbrock_disk -- --path=build/rosenbrock_disk_path.txt

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_hessian_module, only: sqpopt_hessian_exact
    use sqpopt_types_module,   only: sqpopt_results_type, sqpopt_success
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: x_star(2)   = [0.7864152_wp, 0.6176983_wp] !! the solution
    real(wp), parameter :: f_star      = 0.0456748_wp                 !! its objective
    real(wp), parameter :: lambda_star = -0.1214966_wp                !! its multiplier
    real(wp), parameter :: tol         = 1.0e-5_wp                    !! tolerance on them

    integer, parameter :: max_path = 200
    real(wp) :: path(4, max_path) !! `x`, `y`, the objective, and the constraint at each iterate
    real(wp) :: path_kkt(max_path) !! the solver's KKT error at each iterate (with `--path`)
    integer  :: n_path            !! number of iterates recorded

    character(len=512) :: arg, path_file
    integer :: i, u

    write(*,*) '----------------------------'
    write(*,*) 'test_rosenbrock_disk'
    write(*,*) '----------------------------'

    path_file = ''
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        if (arg(1:7) == '--path=') then
            path_file = arg(8:)
        else
            error stop 'test_rosenbrock_disk: unknown option '//trim(arg)
        end if
    end do

    call test('L-BFGS (default)', .false.)
    if (path_file /= '') then
        open(newunit=u, file=trim(path_file), status='replace', action='write')
        write(u, '(A)') '# iterates of test_rosenbrock_disk (default options): iteration, x, y, objective, x^2 + y^2, '// &
                        'KKT error'
        do i = 1, n_path
            write(u, '(I4,5ES24.15)') i - 1, path(:, i), path_kkt(i)
        end do
        close(u)
        print '(A,I0,2A)', 'wrote ', n_path, ' iterates to ', trim(path_file)
    end if
    call test('exact Hessian', .true.)

    print '(A)', 'test_rosenbrock_disk PASSED'

    contains

    subroutine test(label, exact)
    !! solve the problem, and check the solution
    character(len=*), intent(in) :: label !! the configuration, for the output
    logical,          intent(in) :: exact !! whether to use the exact Hessian
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    real(wp) :: x(2), lambda(1)
    integer  :: istat, u_diag
    logical  :: record_kkt

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-1.0e20_wp, -1.0e20_wp], x_ub=[1.0e20_wp, 1.0e20_wp], &
                            c_lb=[-1.0e20_wp], c_ub=[1.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1, 1], icol=[1, 2])
    if (exact) then
        call problem%set_functions(fc=fc, gjac=gjac, hess=hess)
        call problem%set_hessian_sparsity(3, [1, 2, 2], [1, 1, 2])
        options%hessian_mode = sqpopt_hessian_exact
    else
        call problem%set_functions(fc=fc, gjac=gjac)
    end if

    ! (with `--path`, the default solve also records the solver's KKT error at each iterate)
    record_kkt = path_file /= '' .and. .not. exact
    if (record_kkt) then
        open(newunit=u_diag, status='scratch', action='readwrite', form='formatted')
        options%diagnostic_level = 2
        options%diagnostics_unit = u_diag
    end if

    n_path = 0
    path_kkt = -1.0_wp
    call solver%initialize(problem=problem, options=options, report=report)
    call solver%solve([-1.2_wp, 1.0_wp], istat)
    if (record_kkt) then
        call read_kkt(u_diag)
        close(u_diag)
    end if
    call solver%get_solution(x, lambda)
    call solver%get_results(r)

    print '(A18,A,I0,A,I0,A,I0,A,2F11.7,A,F11.7,A,F11.7)', label, ': istat = ', istat, ', iterations = ', &
        r%iterations, ', fc = ', r%n_eval_fc, ', x = ', x, ', f = ', r%f, ', lambda = ', lambda
    if (istat /= sqpopt_success) error stop 'test_rosenbrock_disk FAILED: the solve did not converge'
    if (maxval(abs(x - x_star)) > tol) error stop 'test_rosenbrock_disk FAILED: wrong solution'
    if (abs(r%f - f_star) > tol) error stop 'test_rosenbrock_disk FAILED: wrong objective'
    if (abs(lambda(1) - lambda_star) > tol) error stop 'test_rosenbrock_disk FAILED: wrong multiplier'
    if (sum(x**2) > 1.0_wp + 1.0e-7_wp) error stop 'test_rosenbrock_disk FAILED: the solution is outside the disk'
    if (n_path < 2) error stop 'test_rosenbrock_disk FAILED: the iterates were not reported'
    if (maxval(abs(path(1:2, n_path) - x)) > tol) then
        error stop 'test_rosenbrock_disk FAILED: the last iterate reported is not the solution'
    end if

    end subroutine test

    subroutine fc(x, f, c, status, data)
    !! the objective and the constraint
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(2)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint value at `x` `dimension(1)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (none)
    f = 100.0_wp*(x(2) - x(1)**2)**2 + (1.0_wp - x(1))**2
    c(1) = x(1)**2 + x(2)**2
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient and the constraint Jacobian
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(2)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` `dimension(2)`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    g = [-400.0_wp*x(1)*(x(2) - x(1)**2) - 2.0_wp*(1.0_wp - x(1)), 200.0_wp*(x(2) - x(1)**2)]
    jac_val = 2.0_wp*x
    end subroutine gjac

    subroutine hess(x, lambda, hess_val, status, data)
    !! the Hessian of the Lagrangian `f - lambda*c`, lower triangle: (1,1), (2,1), (2,2)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(2)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multiplier `dimension(1)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x`
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    hess_val = [1200.0_wp*x(1)**2 - 400.0_wp*x(2) + 2.0_wp - 2.0_wp*lambda(1), -400.0_wp*x(1), &
                200.0_wp - 2.0_wp*lambda(1)]
    end subroutine hess

    subroutine read_kkt(u)
    !! read the KKT error of each iteration from the diagnostics' table (see
    !! [[sqpopt_diagnostics_module]]): its row `k` is for the point at the start of
    !! iteration `k`, the iterate `k-1` of the path
    integer, intent(in) :: u !! the unit the diagnostics wrote to (rewound here)
    character(len=1024) :: line
    real(wp) :: f, feas, kkt
    integer :: k, ios
    rewind(u)
    read(u, '(A)', iostat=ios) line   ! (the header)
    do
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) exit
        read(line, *, iostat=ios) k, f, feas, kkt
        if (ios /= 0) cycle
        if (k >= 1 .and. k <= max_path) path_kkt(k) = kkt
    end do
    end subroutine read_kkt

    subroutine report(iter, x, f, c, lambda, user_stop, data)
    !! record the iterate
    integer,                intent(in)    :: iter      !! major iteration number (starts at 1)
    real(wp), dimension(:), intent(in)    :: x         !! current point `dimension(2)`
    real(wp),               intent(in)    :: f         !! current objective value
    real(wp), dimension(:), intent(in)    :: c         !! current constraint value `dimension(1)`
    real(wp), dimension(:), intent(in)    :: lambda    !! current multiplier estimate `dimension(1)`
    logical,                intent(out)   :: user_stop !! set `.true.` to request the solver stop
    class(*), optional,     intent(inout) :: data      !! the user data given to `set_functions` (none)
    user_stop = .false.
    if (n_path < max_path) then
        n_path = n_path + 1
        path(:, n_path) = [x(1), x(2), f, c(1)]
    end if
    end subroutine report

end program test_rosenbrock_disk
