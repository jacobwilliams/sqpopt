program benchmark_large

    !! Large sparse problems with analytic second derivatives, for measuring
    !! what the exact Hessian and the options that use a sparse factorization
    !! (`options%inertia_control`, `direct_qp`, and `direct_least_squares`)
    !! do to the run time, and where the time goes:
    !!
    !!    fpm run --example benchmark_large --profile release -- [--scale=S] [--problem=NAME] [--config=NAME]
    !!        [--no-bfgs] [--no-active-set] [--least-squares] [--no-least-squares] [--memory=K] [--threads=T] [--print=L]
    !!        [--linear-solver=qdldl|mumps] [--sr1]
    !!
    !! With MUMPS (see the README):
    !!
    !!    pixi run run-mumps --example benchmark_large --profile release -- --scale=10
    !!
    !! `--scale=S` multiplies the problem sizes below by `S` (default 1).
    !! `--problem=NAME` runs only that problem (`control`, `rosenbrock`, `wells`,
    !! `circles`, or `hyperbolas`), and `--config=NAME` only that configuration
    !! (`bfgs`, `exact`, `bfgs-direct`, `inertia`, `direct`, or one of the
    !! SR1 ones: see below). `--sr1` runs the SR1 configurations instead of
    !! the others.
    !! `--no-bfgs` leaves out the L-BFGS runs, and `--no-active-set` every run
    !! without `direct_qp` (at large sizes those take most of the time).
    !! `direct_least_squares` is on in the runs with `direct_qp`:
    !! `--least-squares` turns it on in every run, and `--no-least-squares`
    !! off in every run. `--memory=K` sets `options%lbfgs_memory`, `--threads=T`
    !! sets `options%factorization_threads`, `--print=L` sets
    !! `options%print_level`, and `--linear-solver=` sets
    !! `options%linear_solver` (QDLDL by default; `mumps` needs a build with
    !! MUMPS).
    !!
    !! Problems (sizes for `S = 1`):
    !!
    !! * `control`: a discretized nonlinear optimal-control problem with
    !!   `N = 5000` steps (`n = 2N+1` variables, `m = N+1` equality
    !!   constraints, bounds on the controls; many degrees of freedom):
    !!
    !!       minimize   sum_k h*(y_k^2 + u_k^2)/2 + (y_N)^2
    !!       subject to y_0 = 1,
    !!                  y_{k+1} = y_k + h*(u_k - y_k^3),   k = 0..N-1
    !!                  -0.3 <= u_k <= 0.3
    !!
    !! * `rosenbrock`: the chained Rosenbrock function with `n = 5000`
    !!   variables, `n/2` circle constraints, and bounds:
    !!
    !!       minimize   sum_i 100*(x_{i+1}-x_i^2)^2 + (1-x_i)^2
    !!       subject to x_i^2 + x_{i+1}^2 <= 1.5,   i = 1,3,5,...
    !!                  -2 <= x <= 2
    !!
    !! * `wells`: a chain of double-well potentials with `n = 10000`
    !!   variables and `n/4` nonlinear equality constraints. It starts near
    !!   `x = 0`, where every well is concave, so the Hessian of the
    !!   Lagrangian is indefinite there (a test of the Hessian's shift):
    !!
    !!       minimize   sum_i (x_i^4/4 - x_i^2/2 + x_i/10) + sum_i (x_{i+1}-x_i)^2/2
    !!       subject to x_a*x_b + x_c - x_d^2 = 0   for each block (a,b,c,d) of 4 variables
    !!                  -1.5 <= x <= 0.8
    !!
    !! * `circles`: a chain of `n - 1` coupled circle constraints on
    !!   `n = 10000` variables, with the objective of the Maratos example (see
    !!   `test/test_maratos.f90`). It starts on the constraints, where full
    !!   steps increase both the objective and the violation, so the line
    !!   search needs second-order corrections; and each constraint shares a
    !!   variable with the next, so the corrections are ill-conditioned
    !!   least-squares problems (a test of `direct_least_squares`). The
    !!   solution is `x = (1,0,1,0,...)`, with `f = -n/2`:
    !!
    !!       minimize   2*sum_i (x_i^2 + x_{i+1}^2 - 1) - sum_{i odd} x_i
    !!       subject to x_i^2 + x_{i+1}^2 = 1,   i = 1..n-1
    !!
    !! * `hyperbolas` (only with `--problem=hyperbolas`): a chain of `n - 1`
    !!   hyperbola constraints on `n = 10000` variables, from a start where
    !!   the linearized constraints can't be satisfied within the bounds, so
    !!   the first iterations are Gauss-Newton restoration steps (the other
    !!   use of `direct_least_squares`). Most of its time goes to the QPs that
    !!   find the linearization inconsistent, which is why it is not in the
    !!   default set:
    !!
    !!       minimize   sum_i (x_i - a_i)^2/2,   a_i = 1 + sin(0.37 i)/2
    !!       subject to x_i*x_{i+1} >= 1,   i = 1..n-1
    !!                  0.1 <= x <= 2
    !!
    !! Each is solved with L-BFGS (the default; configuration `bfgs`), with the
    !! exact Hessian (`exact`), with L-BFGS and the
    !! direct QP method (`bfgs-direct`), with the exact Hessian and inertia
    !! control (`inertia`), and with those and the direct QP method (`direct`;
    !! the direct runs also use direct least-squares solves). With `--sr1`, it
    !! is solved with the limited-memory SR1 Hessian instead: plain (`sr1`),
    !! shifted by its smallest eigenvalue before each QP (`hessian%convexify`,
    !! `sr1-convexify`), with inertia control (`sr1-inertia`), and the last two
    !! with the direct QP method (`sr1-convexify-direct`,
    !! `sr1-inertia-direct`). The columns are
    !! the status, major iterations, calls of `fc`, the objective, the total
    !! time and the parts of it in the QP solver and in the factorizations,
    !! the number of QPs solved directly out of all QPs, the number of
    !! factorizations, and the numbers of second-order corrections and of
    !! restoration steps.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_hessian_module, only: sqpopt_hessian_exact, sqpopt_hessian_sr1, sqpopt_hessian_bfgs, sqpopt_hessian_type
    use sqpopt_symmetric_solver_module, only: sqpopt_linear_solver_mumps, sqpopt_linear_solver_qdldl
    use sqpopt_types_module,   only: sqpopt_results_type
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    integer  :: nsteps   !! `control`: number of steps `N`
    real(wp) :: h        !! `control`: step size
    real(wp) :: scale    !! `--scale`
    logical  :: with_bfgs, with_active_set
    logical  :: with_sr1 !! `--sr1`: the SR1 configurations instead of the others
    integer  :: least_squares !! `direct_least_squares`: `1` in every run, `-1` in none, `0` in the direct ones
    character(len=:), allocatable :: only_config !! `--config` (empty: all of them)
    integer  :: print_level !! `--print`
    integer  :: memory      !! `--memory` (`0`: automatic)
    integer  :: threads     !! `--threads` (`0`: as the OpenMP environment says)
    integer  :: linear_solver !! `--linear-solver`
    integer  :: i, ios
    character(len=64) :: arg
    character(len=:), allocatable :: only !! `--problem` (empty: all of them)

    scale = 1.0_wp
    with_bfgs = .true.
    with_active_set = .true.
    with_sr1 = .false.
    least_squares = 0
    only_config = ''
    print_level = 0
    memory = 0
    threads = 1
    linear_solver = sqpopt_linear_solver_qdldl
    only = ''
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        if (arg(1:8) == '--scale=') then
            read(arg(9:), *, iostat=ios) scale
            if (ios /= 0 .or. .not. scale > 0.0_wp) error stop 'benchmark_large: bad --scale value'
        else if (arg(1:10) == '--problem=') then
            only = trim(arg(11:))
            if (all(only /= [character(len=10) :: 'control', 'rosenbrock', 'wells', 'circles', 'hyperbolas'])) then
                error stop 'benchmark_large: bad --problem value'
            end if
        else if (arg == '--no-bfgs') then
            with_bfgs = .false.
        else if (arg == '--no-active-set') then
            with_active_set = .false.
        else if (arg == '--sr1') then
            with_sr1 = .true.
        else if (arg == '--no-least-squares') then
            least_squares = -1
        else if (arg == '--least-squares') then
            least_squares = 1
        else if (arg(1:9) == '--config=') then
            only_config = trim(arg(10:))
            if (all(only_config /= [character(len=20) :: 'bfgs', 'exact', 'bfgs-direct', 'inertia', 'direct', 'sr1', &
                                    'sr1-convexify', 'sr1-inertia', 'sr1-convexify-direct', 'sr1-inertia-direct'])) then
                error stop 'benchmark_large: bad --config value'
            end if
        else if (arg(1:9) == '--memory=') then
            read(arg(10:), *, iostat=ios) memory
            if (ios /= 0 .or. memory < 0) error stop 'benchmark_large: bad --memory value'
        else if (arg(1:10) == '--threads=') then
            read(arg(11:), *, iostat=ios) threads
            if (ios /= 0 .or. threads < 0) error stop 'benchmark_large: bad --threads value'
        else if (arg == '--linear-solver=mumps') then
            linear_solver = sqpopt_linear_solver_mumps
        else if (arg == '--linear-solver=qdldl') then
            linear_solver = sqpopt_linear_solver_qdldl
        else if (arg(1:8) == '--print=') then
            read(arg(9:), *, iostat=ios) print_level
            if (ios /= 0) error stop 'benchmark_large: bad --print value'
        else
            error stop 'benchmark_large: unknown option (see the header of example/benchmark_large.f90)'
        end if
    end do

    write(*,'(A)') ''
    write(*,'(A11,A8,A8,2X,A24,A6,A6,A6,A16,3A9,A10,3A6)') 'problem', 'n', 'm', 'configuration          ', 'istat', 'iter', &
        'fc', 'f', 'time', 'QP', 'factor', 'direct', 'fact', 'soc', 'rest'

    if (only == '' .or. only == 'control')    call run_all('control')
    if (only == '' .or. only == 'rosenbrock') call run_all('rosenbrock')
    if (only == '' .or. only == 'wells')      call run_all('wells')
    if (only == '' .or. only == 'circles')    call run_all('circles')
    ! (not in the default set: its inconsistent QPs take most of the time)
    if (only == 'hyperbolas') call run_all('hyperbolas')

    contains

    subroutine run_all(name)
    !! solve one problem with every configuration
    character(len=*), intent(in) :: name !! the problem
    if (with_sr1) then
        if (with_active_set .and. wanted('sr1')) call run(name, 'L-SR1', sqpopt_hessian_sr1, .false., .false.)
        if (with_active_set .and. wanted('sr1-convexify')) then
            call run(name, 'L-SR1, convexify', sqpopt_hessian_sr1, .false., .false., convexify=.true.)
        end if
        if (with_active_set .and. wanted('sr1-inertia')) then
            call run(name, 'L-SR1, inertia', sqpopt_hessian_sr1, .true., .false.)
        end if
        if (wanted('sr1-convexify-direct')) then
            call run(name, 'L-SR1, convexify, direct', sqpopt_hessian_sr1, .false., .true., convexify=.true.)
        end if
        if (wanted('sr1-inertia-direct')) call run(name, 'L-SR1, inertia, direct', sqpopt_hessian_sr1, .true., .true.)
    else
        if (with_bfgs .and. with_active_set .and. wanted('bfgs')) then
            call run(name, 'L-BFGS', sqpopt_hessian_bfgs, .false., .false.)
        end if
        if (with_active_set .and. wanted('exact')) call run(name, 'exact Hessian', sqpopt_hessian_exact, .false., .false.)
        if (with_bfgs .and. wanted('bfgs-direct')) call run(name, 'L-BFGS, direct', sqpopt_hessian_bfgs, .false., .true.)
        if (with_active_set .and. wanted('inertia')) call run(name, 'exact, inertia', sqpopt_hessian_exact, .true., .false.)
        if (wanted('direct')) call run(name, 'exact, inertia, direct', sqpopt_hessian_exact, .true., .true.)
    end if
    write(*,'(A)') ''
    end subroutine run_all

    logical function wanted(config)
    !! whether the configuration `config` is to be run (see `--config`)
    character(len=*), intent(in) :: config !! the configuration's name
    wanted = only_config == '' .or. only_config == config
    end function wanted

    subroutine run(name, config, hessian_mode, inertia, direct, convexify)
    !! solve one problem with one configuration, and print a line of results
    character(len=*), intent(in) :: name         !! the problem
    character(len=*), intent(in) :: config       !! the configuration, for the output
    integer,          intent(in) :: hessian_mode !! `options%hessian_mode`
    logical,          intent(in) :: inertia      !! `options%inertia_control`
    logical,          intent(in) :: direct       !! `options%direct_qp` and `options%direct_least_squares`
    logical, optional, intent(in) :: convexify   !! `hessian%convexify` (default `.false.`)

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_hessian_type) :: hessian
    type(sqpopt_results_type) :: r
    real(wp), dimension(:), allocatable :: x0
    character(len=10) :: share
    integer :: istat

    select case (name)
    case ('control');    call setup_control(nint(5000*scale), problem, x0)
    case ('rosenbrock'); call setup_rosenbrock(2*nint(2500*scale), problem, x0)
    case ('hyperbolas'); call setup_hyperbolas(nint(10000*scale), problem, x0)
    case ('circles');    call setup_circles(nint(10000*scale), problem, x0)
    case default;        call setup_wells(4*nint(2500*scale), problem, x0)
    end select

    options%max_iter = 2000
    options%print_level  = print_level
    options%lbfgs_memory = memory
    options%factorization_threads = threads
    options%hessian_mode = hessian_mode
    if (present(convexify)) hessian%convexify = convexify
    options%inertia_control      = inertia
    options%direct_qp            = direct
    options%direct_least_squares = least_squares == 1 .or. (direct .and. least_squares == 0)
    options%linear_solver        = linear_solver

    call solver%initialize(problem=problem, options=options, hessian=hessian)
    call solver%solve(x0, istat)
    call solver%get_results(r)

    write(share, '(I0,A,I0)') r%n_direct_qp, '/', r%n_qp_solves
    write(*,'(A11,I8,I8,2X,A24,I6,I6,I6,ES16.8,3F9.3,A10,3I6)') name, problem%n, problem%m, config, r%istat, r%iterations, &
        r%n_eval_fc, r%f, r%time, r%time_qp, r%time_factorization, adjustr(share), r%n_factorizations, r%n_soc, &
        r%n_restoration_steps

    end subroutine run

    !------------------------------------------------------------------------
    ! control problem
    !------------------------------------------------------------------------

    subroutine set_uniform_bounds(problem, x_lb, x_ub, c_lb, c_ub)
    !! give every variable the bounds `x_lb` and `x_ub`, and every constraint `c_lb` and `c_ub` (after
    !! `set_problem_size`). The arrays are allocated here: an array expression given as an argument
    !! (`spread(...)`, an array constructor) is a temporary, which some compilers put on the stack, and
    !! these problems are large.
    type(sqpopt_problem_type), intent(inout) :: problem !! the problem definition
    real(wp),                  intent(in)    :: x_lb    !! lower bound of every variable
    real(wp),                  intent(in)    :: x_ub    !! upper bound of every variable
    real(wp),                  intent(in)    :: c_lb    !! lower bound of every constraint
    real(wp),                  intent(in)    :: c_ub    !! upper bound of every constraint
    real(wp), dimension(:), allocatable :: xl, xu, cl, cu
    allocate(xl(problem%n), xu(problem%n), cl(problem%m), cu(problem%m))
    xl = x_lb
    xu = x_ub
    cl = c_lb
    cu = c_ub
    call problem%set_bounds(xl, xu, cl, cu)
    end subroutine set_uniform_bounds

    subroutine set_banded_hessian(problem, tridiagonal)
    !! set the sparsity pattern of the Hessian of the Lagrangian: the diagonal and, if `tridiagonal`, the
    !! subdiagonal after it (the pattern's arrays are allocated here, see `set_uniform_bounds`)
    type(sqpopt_problem_type), intent(inout) :: problem     !! the problem definition
    logical,                   intent(in)    :: tridiagonal !! whether the subdiagonal is in the pattern
    integer, dimension(:), allocatable :: irow, icol
    integer :: n, nnz, i
    n   = problem%n
    nnz = n
    if (tridiagonal) nnz = 2*n - 1
    allocate(irow(nnz), icol(nnz))
    do i = 1, n
        irow(i) = i
        icol(i) = i
    end do
    do i = n + 1, nnz
        irow(i) = i - n + 1
        icol(i) = i - n
    end do
    call problem%set_hessian_sparsity(nnz, irow, icol)
    end subroutine set_banded_hessian

    subroutine setup_control(nn, problem, x0)
    !! the discretized optimal-control problem with `nn` steps
    integer,                             intent(in)  :: nn      !! number of time steps `N` (`n = 2N+1` variables)
    type(sqpopt_problem_type),           intent(out) :: problem !! the problem definition
    real(wp), dimension(:), allocatable, intent(out) :: x0      !! the starting point
    real(wp), dimension(:), allocatable :: x_lb, x_ub, c_b
    integer, dimension(:), allocatable :: irow, icol
    integer :: n, m, k, nnz
    nsteps = nn
    h = 5.0_wp/nsteps
    n = 2*nsteps + 1
    m = nsteps + 1

    ! variables: y_0..y_N are 1..N+1, u_0..u_{N-1} are N+2..2N+1
    allocate(x_lb(n), x_ub(n))
    x_lb = -1.0e20_wp; x_ub = 1.0e20_wp
    x_lb(nsteps+2:n) = -0.3_wp; x_ub(nsteps+2:n) = 0.3_wp

    ! Jacobian pattern: row 1 = y_0; row k+2 = (y_{k+1}, y_k, u_k)
    nnz = 1 + 3*nsteps
    allocate(irow(nnz), icol(nnz))
    irow(1) = 1; icol(1) = 1
    do k = 0, nsteps-1
        irow(2+3*k:4+3*k) = k+2
        icol(2+3*k) = k+2          ! y_{k+1}
        icol(3+3*k) = k+1          ! y_k
        icol(4+3*k) = nsteps+2+k   ! u_k
    end do

    call problem%set_problem_size(n=n, m=m)
    allocate(c_b(m))
    c_b    = 0.0_wp
    c_b(1) = 1.0_wp
    call problem%set_bounds(x_lb, x_ub, c_b, c_b)
    call problem%set_jacobian_sparsity(nnz, irow, icol)
    call problem%set_functions(fc=fc_control, gjac=gjac_control, hess=hess_control)
    ! the Hessian of the Lagrangian is diagonal:
    call set_banded_hessian(problem, tridiagonal=.false.)

    allocate(x0(n))
    x0 = 0.0_wp
    x0(1:nsteps+1) = 1.0_wp
    end subroutine setup_control

    subroutine fc_control(x, f, c, status, data)
    !! objective and constraints of the control problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: k
    f = 0.5_wp*h*(sum(x(1:nsteps)**2) + sum(x(nsteps+2:)**2)) + x(nsteps+1)**2
    c(1) = x(1)
    do k = 0, nsteps-1
        c(k+2) = x(k+2) - x(k+1) - h*(x(nsteps+2+k) - x(k+1)**3)
    end do
    end subroutine fc_control

    subroutine gjac_control(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian values of the control problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: k
    g = h*x
    g(nsteps+1) = 2.0_wp*x(nsteps+1)
    jac(1) = 1.0_wp
    do k = 0, nsteps-1
        jac(2+3*k) = 1.0_wp
        jac(3+3*k) = -1.0_wp + 3.0_wp*h*x(k+1)**2
        jac(4+3*k) = -h
    end do
    end subroutine gjac_control

    subroutine hess_control(x, lambda, hess_val, status, data)
    !! the (diagonal) Hessian of the Lagrangian of the control problem: the
    !! objective's, minus `lambda(k+2)*6*h*y_k` from constraint `k+2`
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: k
    hess_val = h
    hess_val(nsteps+1) = 2.0_wp
    do k = 0, nsteps-1
        hess_val(k+1) = hess_val(k+1) - lambda(k+2)*6.0_wp*h*x(k+1)
    end do
    end subroutine hess_control

    !------------------------------------------------------------------------
    ! chained Rosenbrock with circle constraints
    !------------------------------------------------------------------------

    subroutine setup_rosenbrock(n, problem, x0)
    !! the constrained chained-Rosenbrock problem with `n` variables
    integer,                             intent(in)  :: n       !! number of variables (even)
    type(sqpopt_problem_type),           intent(out) :: problem !! the problem definition
    real(wp), dimension(:), allocatable, intent(out) :: x0      !! the starting point
    integer, dimension(:), allocatable :: irow, icol
    integer :: m, i
    m = n/2
    allocate(irow(2*m), icol(2*m))
    do i = 1, m
        irow(2*i-1:2*i) = i
        icol(2*i-1) = 2*i-1
        icol(2*i)   = 2*i
    end do
    call problem%set_problem_size(n=n, m=m)
    call set_uniform_bounds(problem, -2.0_wp, 2.0_wp, -1.0e20_wp, 1.5_wp)
    call problem%set_jacobian_sparsity(2*m, irow, icol)
    call problem%set_functions(fc=fc_rosen, gjac=gjac_rosen, hess=hess_rosen)
    ! the Hessian of the Lagrangian is tridiagonal: the diagonal, then the subdiagonal
    call set_banded_hessian(problem, tridiagonal=.true.)
    allocate(x0(n))
    do i = 1, n
        x0(i) = merge(-1.2_wp, 1.0_wp, mod(i,2) == 1)
    end do
    end subroutine setup_rosenbrock

    subroutine fc_rosen(x, f, c, status, data)
    !! objective and constraints of the chained Rosenbrock problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n, i
    n = size(x)
    f = sum(100.0_wp*(x(2:n)-x(1:n-1)**2)**2 + (1.0_wp-x(1:n-1))**2)
    do i = 1, size(c)
        c(i) = x(2*i-1)**2 + x(2*i)**2
    end do
    end subroutine fc_rosen

    subroutine gjac_rosen(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian values of the chained Rosenbrock problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n, i
    n = size(x)
    g = 0.0_wp
    g(1:n-1) = -400.0_wp*x(1:n-1)*(x(2:n)-x(1:n-1)**2) - 2.0_wp*(1.0_wp-x(1:n-1))
    g(2:n)   = g(2:n) + 200.0_wp*(x(2:n)-x(1:n-1)**2)
    do i = 1, size(jac)/2
        jac(2*i-1) = 2.0_wp*x(2*i-1)
        jac(2*i)   = 2.0_wp*x(2*i)
    end do
    end subroutine gjac_rosen

    subroutine hess_rosen(x, lambda, hess_val, status, data)
    !! the (tridiagonal) Hessian of the Lagrangian of the chained Rosenbrock
    !! problem: the diagonal (`n` values), then the subdiagonal (`n-1`)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: n, i
    n = size(x)
    hess_val = 0.0_wp
    hess_val(1:n-1) = 1200.0_wp*x(1:n-1)**2 - 400.0_wp*x(2:n) + 2.0_wp
    hess_val(2:n)   = hess_val(2:n) + 200.0_wp
    hess_val(n+1:2*n-1) = -400.0_wp*x(1:n-1)
    ! the circle constraints x_{2i-1}^2 + x_{2i}^2, with Hessian 2*I on their pair:
    do i = 1, size(lambda)
        hess_val(2*i-1) = hess_val(2*i-1) - 2.0_wp*lambda(i)
        hess_val(2*i)   = hess_val(2*i)   - 2.0_wp*lambda(i)
    end do
    end subroutine hess_rosen

    !------------------------------------------------------------------------
    ! chain of double wells (nonconvex)
    !------------------------------------------------------------------------

    subroutine setup_wells(n, problem, x0)
    !! the chain of double wells with `n` variables (a multiple of 4)
    integer,                             intent(in)  :: n       !! number of variables
    type(sqpopt_problem_type),           intent(out) :: problem !! the problem definition
    real(wp), dimension(:), allocatable, intent(out) :: x0      !! the starting point
    integer, dimension(:), allocatable :: irow, icol
    integer :: m, k, i
    m = n/4
    ! constraint k has the variables 4k-3 .. 4k (a, b, c, d):
    allocate(irow(4*m), icol(4*m))
    do k = 1, m
        irow(4*k-3:4*k) = k
        icol(4*k-3:4*k) = [4*k-3, 4*k-2, 4*k-1, 4*k]
    end do
    call problem%set_problem_size(n=n, m=m)
    call set_uniform_bounds(problem, -1.5_wp, 0.8_wp, 0.0_wp, 0.0_wp)
    call problem%set_jacobian_sparsity(4*m, irow, icol)
    call problem%set_functions(fc=fc_wells, gjac=gjac_wells, hess=hess_wells)
    ! the Hessian of the Lagrangian is tridiagonal: the diagonal, then the subdiagonal
    call set_banded_hessian(problem, tridiagonal=.true.)
    allocate(x0(n))
    do i = 1, n
        x0(i) = 0.1_wp*sin(real(i, wp))
    end do
    end subroutine setup_wells

    subroutine fc_wells(x, f, c, status, data)
    !! objective and constraints of the chain of double wells
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n, k
    n = size(x)
    f = sum(0.25_wp*x**4 - 0.5_wp*x**2 + 0.1_wp*x) + 0.5_wp*sum((x(2:n) - x(1:n-1))**2)
    do k = 1, size(c)
        c(k) = x(4*k-3)*x(4*k-2) + x(4*k-1) - x(4*k)**2
    end do
    end subroutine fc_wells

    subroutine gjac_wells(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian values of the chain of double wells
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n, k
    n = size(x)
    g = x**3 - x + 0.1_wp
    g(1:n-1) = g(1:n-1) - (x(2:n) - x(1:n-1))
    g(2:n)   = g(2:n)   + (x(2:n) - x(1:n-1))
    do k = 1, size(jac)/4
        jac(4*k-3:4*k) = [x(4*k-2), x(4*k-3), 1.0_wp, -2.0_wp*x(4*k)]
    end do
    end subroutine gjac_wells

    subroutine hess_wells(x, lambda, hess_val, status, data)
    !! the (tridiagonal) Hessian of the Lagrangian of the chain of double
    !! wells: the diagonal (`n` values), then the subdiagonal (`n-1`)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: n, k
    n = size(x)
    hess_val(1:n) = 3.0_wp*x**2 - 1.0_wp + 2.0_wp
    hess_val(1)   = hess_val(1) - 1.0_wp
    hess_val(n)   = hess_val(n) - 1.0_wp
    hess_val(n+1:2*n-1) = -1.0_wp
    ! constraint k: x_a*x_b + x_c - x_d^2, with Hessian elements (b,a) = 1 and (d,d) = -2:
    do k = 1, size(lambda)
        hess_val(n+4*k-3) = hess_val(n+4*k-3) - lambda(k)
        hess_val(4*k)     = hess_val(4*k) + 2.0_wp*lambda(k)
    end do
    end subroutine hess_wells

    !------------------------------------------------------------------------
    ! chain of hyperbolas
    !------------------------------------------------------------------------

    subroutine setup_hyperbolas(n, problem, x0)
    !! the chain of hyperbola constraints with `n` variables
    integer,                             intent(in)  :: n       !! number of variables
    type(sqpopt_problem_type),           intent(out) :: problem !! the problem definition
    real(wp), dimension(:), allocatable, intent(out) :: x0      !! the starting point
    integer, dimension(:), allocatable :: irow, icol
    integer :: m, i
    m = n - 1
    allocate(irow(2*m), icol(2*m))
    do i = 1, m
        irow(2*i-1:2*i) = i
        icol(2*i-1:2*i) = [i, i+1]
    end do
    call problem%set_problem_size(n=n, m=m)
    call set_uniform_bounds(problem, 0.1_wp, 2.0_wp, 1.0_wp, 1.0e20_wp)
    call problem%set_jacobian_sparsity(2*m, irow, icol)
    call problem%set_functions(fc=fc_hyperbolas, gjac=gjac_hyperbolas, hess=hess_hyperbolas)
    ! the Hessian of the Lagrangian is tridiagonal: the diagonal, then the subdiagonal
    call set_banded_hessian(problem, tridiagonal=.true.)
    allocate(x0(n))
    x0 = 0.1_wp
    end subroutine setup_hyperbolas

    pure function hyperbola_target(i) result(a)
    !! the value the objective pulls variable `i` toward
    integer, intent(in) :: i !! index of the variable
    real(wp) :: a
    a = 1.0_wp + 0.5_wp*sin(0.37_wp*real(i, wp))
    end function hyperbola_target

    subroutine fc_hyperbolas(x, f, c, status, data)
    !! objective and constraints of the chain of hyperbolas
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n, i
    n = size(x)
    f = 0.0_wp
    do i = 1, n
        f = f + (x(i) - hyperbola_target(i))**2
    end do
    f = 0.5_wp*f
    c = x(1:n-1)*x(2:n)
    end subroutine fc_hyperbolas

    subroutine gjac_hyperbolas(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian values of the chain of hyperbolas
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n, i
    n = size(x)
    do i = 1, n
        g(i) = x(i) - hyperbola_target(i)
    end do
    do i = 1, n-1
        jac(2*i-1:2*i) = [x(i+1), x(i)]
    end do
    end subroutine gjac_hyperbolas

    subroutine hess_hyperbolas(x, lambda, hess_val, status, data)
    !! the (tridiagonal) Hessian of the Lagrangian of the chain of hyperbolas:
    !! the identity, minus `lambda_i` between the two variables of constraint `i`
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: n
    n = size(x)
    hess_val(1:n)  = 1.0_wp
    hess_val(n+1:) = -lambda
    end subroutine hess_hyperbolas

    !------------------------------------------------------------------------
    ! chain of circles
    !------------------------------------------------------------------------

    subroutine setup_circles(n, problem, x0)
    !! the chain of circle constraints with `n` variables
    integer,                             intent(in)  :: n       !! number of variables
    type(sqpopt_problem_type),           intent(out) :: problem !! the problem definition
    real(wp), dimension(:), allocatable, intent(out) :: x0      !! the starting point
    integer, dimension(:), allocatable :: irow, icol
    integer :: m, i
    m = n - 1
    allocate(irow(2*m), icol(2*m))
    do i = 1, m
        irow(2*i-1:2*i) = i
        icol(2*i-1:2*i) = [i, i+1]
    end do
    call problem%set_problem_size(n=n, m=m)
    call set_uniform_bounds(problem, -1.0e20_wp, 1.0e20_wp, 1.0_wp, 1.0_wp)
    call problem%set_jacobian_sparsity(2*m, irow, icol)
    call problem%set_functions(fc=fc_circles, gjac=gjac_circles, hess=hess_circles)
    call set_banded_hessian(problem, tridiagonal=.false.)
    allocate(x0(n))
    do i = 1, n
        x0(i) = merge(cos(0.3_wp), sin(0.3_wp), mod(i,2) == 1)
    end do
    end subroutine setup_circles

    subroutine fc_circles(x, f, c, status, data)
    !! objective and constraints of the chain of circles
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n
    n = size(x)
    c = x(1:n-1)**2 + x(2:n)**2
    f = 2.0_wp*sum(c - 1.0_wp) - sum(x(1:n:2))
    end subroutine fc_circles

    subroutine gjac_circles(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian values of the chain of circles
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n, i
    n = size(x)
    g = 8.0_wp*x
    g(1) = 4.0_wp*x(1)
    g(n) = 4.0_wp*x(n)
    g(1:n:2) = g(1:n:2) - 1.0_wp
    do i = 1, n-1
        jac(2*i-1:2*i) = [2.0_wp*x(i), 2.0_wp*x(i+1)]
    end do
    end subroutine gjac_circles

    subroutine hess_circles(x, lambda, hess_val, status, data)
    !! the (diagonal) Hessian of the Lagrangian of the chain of circles
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: n
    n = size(x)
    hess_val = 8.0_wp
    hess_val(1) = 4.0_wp
    hess_val(n) = 4.0_wp
    hess_val(1:n-1) = hess_val(1:n-1) - 2.0_wp*lambda
    hess_val(2:n)   = hess_val(2:n)   - 2.0_wp*lambda
    end subroutine hess_circles

end program benchmark_large
