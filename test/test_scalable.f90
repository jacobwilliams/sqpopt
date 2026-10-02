program test_scalable

    !! The scalable test functions of [[scalable_functions_module]] (converted
    !! from the Julia package NonlinearOptimizationTestFunctions.jl): bound
    !! constrained problems, with no other constraints, for any number of
    !! variables.
    !!
    !! First the derivatives of every function are checked by central
    !! differences (the gradient against the function, and the Hessian against
    !! the gradient), at a point near its starting point, for 8 variables.
    !!
    !! Then every function is solved from its starting point, for each number
    !! of variables, with the limited-memory BFGS Hessian, and, if it has a
    !! sparse Hessian, with that; in a build with MUMPS, also with the exact
    !! Hessian, inertia control, and the direct QP method. The solver must
    !! converge (`sqpopt_success`, `sqpopt_acceptable`, or `sqpopt_stalled`)
    !! to a point within the bounds, with a smaller objective than the
    !! starting point's, where the gradient, projected on the bounds, is small
    !! compared with the starting point's. The result is then one of:
    !!
    !! * `global`: the objective is within `f_tol` of the known minimum;
    !! * `loose`: it isn't, but the distance to the minimum was reduced by a
    !!   factor of `gap_tol`. The solver scales the objective by its gradient
    !!   at the starting point, and tests convergence on the scaled problem,
    !!   so from a starting point with a very large gradient (`zakharov`,
    !!   `schwefel12`) it stops at a point that is not close to the minimizer
    !!   in absolute terms;
    !! * `local`: another stationary point, for a function that has other
    !!   local minimizers (`scalable_multimodal`);
    !! * `FAILED`: anything else.
    !!
    !! The test fails if a solve fails that isn't in the list of the known
    !! failures (`known_failures`), or one in the list passes.
    !!
    !! By default the numbers of variables are 20 and 1000. The iteration
    !! limit is high because `rosenbrock` needs about 5 iterations per
    !! variable with the limited-memory BFGS Hessian, and 1.7 with the exact
    !! one. These problems have no constraints, so with the limited-memory
    !! BFGS Hessian their QPs are solved by the unconstrained step (see
    !! [[sqpopt_qp_solver_module]]) wherever the bounds aren't active.
    !!
    !! Command-line options (any of them makes the run a non-default
    !! configuration, which runs every solve and reports the results, but
    !! doesn't fail the test):
    !!
    !! * `--n=N1,N2,...`: the numbers of variables (multiples of 4, for
    !!   `powell_singular`)
    !! * `--function=NAME`: only that function
    !! * `--hessian=bfgs|exact|direct`: only that configuration (`direct` is
    !!   the exact Hessian with inertia control and the direct QP method, and
    !!   needs a build with MUMPS)
    !! * `--max-iter=K`: `options%max_iter` (default 10000)
    !! * `--print=L`: `options%print_level`

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_exact
    use sqpopt_inertia_module,   only: sqpopt_has_mumps
    use sqpopt_types_module,     only: sqpopt_results_type, sqpopt_success, sqpopt_acceptable, sqpopt_stalled
    use sqpopt_kinds,            only: wp => sqpopt_module_wp
    use scalable_functions_module

    implicit none

    integer, parameter :: cfg_bfgs = 1, cfg_exact = 2, cfg_direct = 3 !! the configurations
    character(len=*), parameter :: cfg_name(3) = [character(len=6) :: 'L-BFGS', 'exact', 'direct']

    real(wp), parameter :: f_tol   = 1.0e-5_wp !! tolerance on the minimum, relative to `max(1, |f_min|)`
    real(wp), parameter :: gap_tol = 1.0e-8_wp !! reduction of `f - f_min` from the starting point for a `loose` result
    real(wp), parameter :: pg_tol = 1.0e-3_wp !! tolerance on the projected gradient, relative to `max(1, |g(x0)|)`
    real(wp), parameter :: fd_tol = 1.0e-5_wp !! relative tolerance of the derivative checks

    type :: known_failure
        !! a solve that is known to fail
        character(len=20) :: name !! the function
        integer :: n              !! the number of variables
        integer :: cfg            !! the configuration
    end type known_failure

    !> the solves of the default configuration that are known to fail (see the comments at each)
    type(known_failure), parameter :: known_failures(0) = [known_failure ::]

    integer, dimension(:), allocatable :: sizes
    character(len=:), allocatable :: only_function
    integer :: only_cfg, max_iter, print_level
    logical :: default_run
    integer :: id, k, cfg, n_global, n_loose, n_local, n_failed, n_unexpected
    type(scalable_function_type) :: fun

    write(*,*) '----------------------------'
    write(*,*) 'test_scalable'
    write(*,*) '----------------------------'

    call read_arguments()

    do id = 1, n_scalable_functions
        call check_derivatives(id)
    end do
    print '(A)', 'test_scalable [derivatives] PASSED'

    n_global = 0
    n_loose  = 0
    n_local  = 0
    n_failed = 0
    n_unexpected = 0
    print '(/,A)', '          function       n  config  istat   iter     fc                f            f_min'// &
                   '   |proj g|   time(s)  result'
    do id = 1, n_scalable_functions
        if (only_function /= '' .and. only_function /= scalable_function_name(id)) cycle
        do k = 1, size(sizes)
            call scalable_function_setup(id, sizes(k), fun)
            do cfg = cfg_bfgs, cfg_direct
                if (only_cfg /= 0 .and. cfg /= only_cfg) cycle
                if (cfg /= cfg_bfgs .and. .not. fun%has_hessian) cycle
                if (cfg == cfg_direct .and. .not. sqpopt_has_mumps) cycle
                call solve(fun, cfg)
            end do
        end do
    end do

    print '(/,4(A,I0))', 'global minimum: ', n_global, ', loose: ', n_loose, ', another stationary point: ', n_local, &
                         ', failed: ', n_failed
    if (.not. default_run) then
        print '(A)', 'test_scalable: non-default configuration, the results are not checked'
    else if (n_unexpected > 0) then
        error stop 'test_scalable FAILED: a solve failed that is not a known failure, or a known failure passed'
    else
        print '(A)', 'test_scalable PASSED'
    end if

    contains

    subroutine read_arguments()
    !! read the command-line options (see the program's documentation)
    character(len=256) :: arg
    integer :: i, j, ios, count
    sizes = [20, 1000]
    only_function = ''
    only_cfg    = 0
    max_iter    = 10000
    print_level = 0
    default_run = command_argument_count() == 0
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        if (arg(1:4) == '--n=') then
            count = 1
            do j = 5, len_trim(arg)
                if (arg(j:j) == ',') count = count + 1
            end do
            deallocate(sizes)
            allocate(sizes(count))
            read(arg(5:), *, iostat=ios) sizes
            if (ios /= 0 .or. any(sizes < 4) .or. any(mod(sizes, 4) /= 0)) then
                error stop 'test_scalable: --n needs numbers that are multiples of 4'
            end if
        else if (arg(1:11) == '--function=') then
            only_function = trim(arg(12:))
        else if (arg == '--hessian=bfgs') then
            only_cfg = cfg_bfgs
        else if (arg == '--hessian=exact') then
            only_cfg = cfg_exact
        else if (arg == '--hessian=direct') then
            only_cfg = cfg_direct
        else if (arg(1:11) == '--max-iter=') then
            read(arg(12:), *, iostat=ios) max_iter
            if (ios /= 0) error stop 'test_scalable: bad --max-iter'
        else if (arg(1:8) == '--print=') then
            read(arg(9:), *, iostat=ios) print_level
            if (ios /= 0) error stop 'test_scalable: bad --print'
        else
            error stop 'test_scalable: unknown option '//trim(arg)
        end if
    end do
    end subroutine read_arguments

    subroutine check_derivatives(id)
    !! check the gradient and the Hessian of function `id` by central differences, for 8 variables
    integer, intent(in) :: id !! which function
    integer, parameter :: n = 8
    type(scalable_function_type) :: fun
    real(wp) :: x(n), xp(n), g(n), gp(n), gm(n), g_fd(n), h_fd(n,n), h(n,n), step, err
    real(wp), dimension(:), allocatable :: hval
    integer :: i, k

    call scalable_function_setup(id, n, fun)
    ! (a point without symmetries, inside the bounds)
    x = fun%x0 + [(0.3_wp*sin(1.7_wp*real(i, wp)) + 0.05_wp*real(i, wp), i=1,n)]
    x = min(max(x, fun%x_lb + 1.0e-3_wp), fun%x_ub - 1.0e-3_wp)

    call fun%g(x, g)
    do i = 1, n
        step = 1.0e-6_wp*max(1.0_wp, abs(x(i)))
        xp = x
        xp(i) = x(i) + step
        g_fd(i) = fun%f(xp)
        call fun%g(xp, gp)
        xp(i) = x(i) - step
        g_fd(i) = (g_fd(i) - fun%f(xp))/(2.0_wp*step)
        call fun%g(xp, gm)
        h_fd(:,i) = (gp - gm)/(2.0_wp*step)
    end do
    err = maxval(abs(g - g_fd))/max(1.0_wp, maxval(abs(g)))
    if (.not. err <= fd_tol) then
        print '(A,A,ES10.2)', fun%name, ': gradient error ', err
        error stop 'test_scalable FAILED: a gradient differs from its finite differences'
    end if

    if (fun%has_hessian) then
        allocate(hval(size(fun%hess_irow)))
        call fun%h(x, hval)
        h = 0.0_wp
        do k = 1, size(hval)
            h(fun%hess_irow(k), fun%hess_icol(k)) = hval(k)
            h(fun%hess_icol(k), fun%hess_irow(k)) = hval(k)
        end do
        err = maxval(abs(h - h_fd))/max(1.0_wp, maxval(abs(h)))
        if (.not. err <= fd_tol) then
            print '(A,A,ES10.2)', fun%name, ': Hessian error ', err
            error stop 'test_scalable FAILED: a Hessian differs from its finite differences'
        end if
    end if
    end subroutine check_derivatives

    subroutine solve(fun, cfg)
    !! solve the problem of `fun` with configuration `cfg`, and report and count the result
    type(scalable_function_type), intent(inout) :: fun !! the function
    integer,                      intent(in)    :: cfg !! the configuration (`cfg_bfgs`, `cfg_exact`, or `cfg_direct`)

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    real(wp), dimension(fun%n) :: x, g
    real(wp), dimension(0) :: lambda
    integer, dimension(0) :: no_rows
    real(wp) :: f0, g0, pg
    integer :: istat, n, j
    logical :: converged, reached, known
    character(len=6) :: verdict

    n = fun%n
    call problem%set_problem_size(n=n, m=0)
    call problem%set_bounds(x_lb=fun%x_lb, x_ub=fun%x_ub, c_lb=[real(wp) ::], c_ub=[real(wp) ::])
    call problem%set_jacobian_sparsity(nnz=0, irow=no_rows, icol=no_rows)
    if (cfg == cfg_bfgs) then
        call problem%set_functions(fc=fc, gjac=gjac, data=fun)
    else
        call problem%set_functions(fc=fc, gjac=gjac, hess=hess, data=fun)
        call problem%set_hessian_sparsity(size(fun%hess_irow), fun%hess_irow, fun%hess_icol)
        options%hessian_mode = sqpopt_hessian_exact
    end if
    if (cfg == cfg_direct) then
        options%inertia_control = .true.
        options%direct_qp       = .true.
    end if
    options%max_iter    = max_iter
    options%print_level = print_level

    f0 = fun%f(fun%x0)
    call fun%g(fun%x0, g)
    g0 = maxval(abs(g))

    call solver%initialize(problem=problem, options=options)
    call solver%solve(fun%x0, istat)
    call solver%get_solution(x, lambda)
    call solver%get_results(r)

    ! the gradient, projected on the bounds:
    call fun%g(x, g)
    where (x <= fun%x_lb) g = min(g, 0.0_wp)
    where (x >= fun%x_ub) g = max(g, 0.0_wp)
    pg = maxval(abs(g))

    converged = (istat == sqpopt_success .or. istat == sqpopt_acceptable .or. istat == sqpopt_stalled) .and. &
                all(x >= fun%x_lb) .and. all(x <= fun%x_ub) .and. r%f <= f0 .and. pg <= pg_tol*max(1.0_wp, g0)
    reached = r%f - fun%f_min <= f_tol*max(1.0_wp, abs(fun%f_min))
    if (converged .and. reached) then
        verdict = 'global'
        n_global = n_global + 1
    else if (converged .and. r%f - fun%f_min <= gap_tol*(f0 - fun%f_min)) then
        verdict = 'loose'
        n_loose = n_loose + 1
    else if (converged .and. fun%kind == scalable_multimodal) then
        verdict = 'local'
        n_local = n_local + 1
    else
        verdict = 'FAILED'
        n_failed = n_failed + 1
    end if

    known = .false.
    do j = 1, size(known_failures)
        if (known_failures(j)%name == fun%name .and. known_failures(j)%n == n .and. known_failures(j)%cfg == cfg) known = .true.
    end do
    if (known .neqv. (verdict == 'FAILED')) n_unexpected = n_unexpected + 1

    print '(A18,I8,2X,A6,I7,I7,I7,2ES17.8,ES11.2,F10.3,2X,A,A)', fun%name, n, cfg_name(cfg), istat, r%iterations, &
        r%n_eval_fc, r%f, fun%f_min, pg, r%time, trim(verdict), merge(' (known)', '        ', known)

    end subroutine solve

    subroutine fc(x, f, c, status, data)
    !! the objective (there are no constraints)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(0)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the function (a `scalable_function_type`)
    f = 0.0_wp
    select type (data)
    type is (scalable_function_type)
        f = data%f(x)
    class default
        status = -1
    end select
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient (the Jacobian is empty)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` `dimension(0)`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the function (a `scalable_function_type`)
    g = 0.0_wp
    select type (data)
    type is (scalable_function_type)
        call data%g(x, g)
    class default
        status = -1
    end select
    end subroutine gjac

    subroutine hess(x, lambda, hess_val, status, data)
    !! the Hessian of the Lagrangian, which is the objective's
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(0)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the function (a `scalable_function_type`)
    hess_val = 0.0_wp
    select type (data)
    type is (scalable_function_type)
        call data%h(x, hess_val)
    class default
        status = -1
    end select
    end subroutine hess

end program test_scalable
