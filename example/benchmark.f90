program benchmark

    !! Scalable benchmark problems for measuring `sqpopt`'s efficiency
    !! (function evaluations and run time), at sizes where the default
    !! `sqpopt_qp_auto` mode picks the dense QP solver and at sizes where it
    !! picks the sparse reduced-Hessian one. Run with:
    !!
    !!    fpm run --example benchmark --profile release
    !!
    !! Command-line options: `--qp=auto|dense|sparse|daqp` (the QP solver,
    !! `options%qp_solver_mode`; the default is `auto`), and `--max-n=N` (skip
    !! the problems with more than `N` variables: the dense QP solvers form
    !! dense matrices of order `n`).
    !!
    !! Problems:
    !!
    !! * `control`: a discretized nonlinear optimal-control problem with `N`
    !!   steps (`n = 2N+1` variables: states `y_0..y_N` and controls
    !!   `u_0..u_{N-1}`; `m = N+1` equality constraints):
    !!
    !!       minimize   sum_k h*(y_k^2 + u_k^2)/2 + (y_N)^2
    !!       subject to y_0 = 1,
    !!                  y_{k+1} = y_k + h*(u_k - y_k^3),   k = 0..N-1
    !!                  -0.3 <= u_k <= 0.3
    !!
    !! * `rosenbrock`: the chained Rosenbrock function with `n/2` coupled
    !!   circle constraints (`m = n/2` nonlinear inequalities) and bounds:
    !!
    !!       minimize   sum_i 100*(x_{i+1}-x_i^2)^2 + (1-x_i)^2
    !!       subject to x_i^2 + x_{i+1}^2 <= 1.5,   i = 1,3,5,...
    !!                  -2 <= x <= 2

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_types_module,   only: sqpopt_results_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_auto, sqpopt_qp_dense, sqpopt_qp_reduced_hessian, sqpopt_qp_daqp
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    integer :: nsteps               !! `control`: number of steps `N`
    real(wp) :: h                   !! `control`: step size
    integer :: qp_mode              !! `--qp=`: `options%qp_solver_mode`
    integer :: max_n                !! `--max-n=`: the largest number of variables to solve

    call parse_arguments()
    write(*,'(A)') ''
    write(*,'(A12,A8,A6,A6,A6,A8,A8,A10,A16)') 'problem', 'size', 'n', 'm', 'istat', &
        'n_fc', 'n_gjac', 'time (s)', 'f'

    if (2*50+1 <= max_n)  call run_control(50)
    if (2*150+1 <= max_n) call run_control(150)
    if (2*500+1 <= max_n) call run_control(500)
    if (40 <= max_n)      call run_rosenbrock(40)
    if (300 <= max_n)     call run_rosenbrock(300)
    if (2000 <= max_n)    call run_rosenbrock(2000)

    contains

    subroutine parse_arguments()
    !! read the command-line options (see the program's documentation)
    character(len=64) :: arg
    integer :: i, ios
    qp_mode = sqpopt_qp_auto
    max_n   = huge(1)
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        select case (trim(arg))
        case ('--qp=auto');   qp_mode = sqpopt_qp_auto
        case ('--qp=dense');  qp_mode = sqpopt_qp_dense
        case ('--qp=sparse'); qp_mode = sqpopt_qp_reduced_hessian
        case ('--qp=daqp');   qp_mode = sqpopt_qp_daqp
        case default
            if (arg(1:8) == '--max-n=') then
                read(arg(9:), *, iostat=ios) max_n
                if (ios /= 0) error stop 'benchmark: bad --max-n value'
            else
                error stop 'benchmark: unknown option '//trim(arg)
            end if
        end select
    end do
    end subroutine parse_arguments

    subroutine report(name, size_param, n, m, solver)
    !! print one line of results (from the solver's results object)
    character(len=*), intent(in)  :: name       !! the problem
    integer, intent(in)           :: size_param !! its size parameter (time steps, or variables)
    integer, intent(in)           :: n          !! number of variables
    integer, intent(in)           :: m          !! number of constraints
    type(sqpopt_type), intent(in) :: solver     !! the solver, after the solve
    type(sqpopt_results_type) :: r
    call solver%get_results(r)
    write(*,'(A12,I8,I6,I6,I6,I8,I8,F10.3,ES16.8)') name, size_param, n, m, r%istat, &
        r%n_eval_fc, r%n_eval_gjac, r%time, r%f
    end subroutine report

    !------------------------------------------------------------------------
    ! control problem
    !------------------------------------------------------------------------

    subroutine run_control(nn)
    !! solve the discretized optimal-control problem with `nn` time steps, and print the results
    integer, intent(in) :: nn !! number of time steps `N` (`n = 2N+1` variables)
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp), dimension(:), allocatable :: x, x_lb, x_ub
    integer, dimension(:), allocatable :: irow, icol
    integer :: n, m, k, nnz, istat
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
    call problem%set_bounds(x_lb, x_ub, spread(0.0_wp,1,m), spread(0.0_wp,1,m))
    problem%c_lb(1) = 1.0_wp; problem%c_ub(1) = 1.0_wp
    call problem%set_jacobian_sparsity(nnz, irow, icol)
    call problem%set_functions(fc=fc_control, gjac=gjac_control)

    options%max_iter = 2000
    options%qp_solver_mode = qp_mode
    allocate(x(n))
    x = 0.0_wp
    x(1:nsteps+1) = 1.0_wp

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x, istat)
    call report('control', nsteps, n, m, solver)

    end subroutine run_control

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

    !------------------------------------------------------------------------
    ! chained Rosenbrock with circle constraints
    !------------------------------------------------------------------------

    subroutine run_rosenbrock(n)
    !! solve the constrained chained-Rosenbrock problem with `n` variables, and print the results
    integer, intent(in) :: n !! number of variables
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp), dimension(:), allocatable :: x
    integer, dimension(:), allocatable :: irow, icol
    integer :: m, i, istat
    m = n/2
    allocate(irow(2*m), icol(2*m))
    do i = 1, m
        irow(2*i-1:2*i) = i
        icol(2*i-1) = 2*i-1
        icol(2*i)   = 2*i
    end do
    call problem%set_problem_size(n=n, m=m)
    call problem%set_bounds(spread(-2.0_wp,1,n), spread(2.0_wp,1,n), spread(-1.0e20_wp,1,m), spread(1.5_wp,1,m))
    call problem%set_jacobian_sparsity(2*m, irow, icol)
    call problem%set_functions(fc=fc_rosen, gjac=gjac_rosen)

    options%max_iter = 2000
    options%qp_solver_mode = qp_mode
    allocate(x(n))
    do i = 1, n
        x(i) = merge(-1.2_wp, 1.0_wp, mod(i,2) == 1)
    end do

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x, istat)
    call report('rosenbrock', n, n, m, solver)

    end subroutine run_rosenbrock

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

end program benchmark
