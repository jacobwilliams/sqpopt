program settings_study

    !! Measures how the QP solver and Hessian settings perform as the size and
    !! the number of degrees of freedom of a problem change (the numbers
    !! behind the guide's "Choosing settings" page). Run with:
    !!
    !!    fpm run --example settings_study --profile release
    !!
    !! Problems:
    !!
    !! * `control`: a discretized nonlinear optimal-control problem with `N`
    !!   steps and `K` piecewise-constant controls (`u_j` applies on `N/K`
    !!   consecutive steps), so `n = N+1+K`, `m = N+1` equality constraints,
    !!   and `K` degrees of freedom:
    !!
    !!       minimize   sum_k h*y_k^2/2 + (y_N)^2 + (N/K)*sum_j h*u_j^2/2
    !!       subject to y_0 = 1,
    !!                  y_{k+1} = y_k + h*(u_{j(k)} - y_k^3),   k = 0..N-1
    !!                  -0.3 <= u_j <= 0.3
    !!
    !! * `elliptic`: a semilinear elliptic control problem on an `N x N` grid
    !!   (`n = 2N^2`, `m = N^2` equality constraints with the 5-point
    !!   Laplacian stencil; about `N^2` degrees of freedom, less the active
    !!   control bounds), as in `test/test_sparse_options.f90`:
    !!
    !!       minimize   sum (y - y_d)^2/2 + alpha*sum u^2/2
    !!       subject to A*y + y^3 - u = 0,   u_l <= u <= u_u
    !!
    !! * `dense`: `n` variables and `n/2` constraints (half equalities, half
    !!   inequalities), each depending on every variable (see `setup_dense`).
    !!
    !! Studies:
    !!
    !! 1. QP solver against problem size: the dense QP and the sparse QP (with
    !!    `LU` and with `LSQR` null spaces), on `control` and `elliptic` with
    !!    many degrees of freedom (`K = N`), from `n` of about 50 to 1600.
    !! 2. Degrees of freedom: the `control` problem with `N = 1000` steps and
    !!    `K` from 5 to 1000 controls, with the L-BFGS memory, the exact
    !!    Hessian, and `dense_max_ns`.
    !! 3. QP solver against problem size, on `dense`, from `n = 20` to 400.
    !!
    !! (The dense QP is skipped above `n = 500`, and `LSQR` above `n = 1000`:
    !! each would take minutes.)

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_bfgs, sqpopt_hessian_exact
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type, sqpopt_qp_dense, sqpopt_qp_reduced_hessian
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_null_space_lu, sqpopt_null_space_lsqr
    use sqpopt_types_module,     only: sqpopt_results_type, sqpopt_infinity
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: pi = acos(-1.0_wp)
    integer,  parameter :: dense_max_n = 500 !! the dense QP is skipped above this many variables (it takes minutes)
    integer,  parameter :: lsqr_max_n = 1000 !! and the `LSQR` null space above this many

    type :: config_type
        !! a configuration of the QP solver and the Hessian
        character(len=26) :: name = ''                     !! label, for the output
        integer :: qp_mode      = sqpopt_qp_reduced_hessian !! `options%qp_solver_mode`
        integer :: hessian_mode = sqpopt_hessian_bfgs      !! `options%hessian_mode`
        integer :: lbfgs_memory = 0                        !! `options%lbfgs_memory` (`0`: the default)
        integer :: null_space   = sqpopt_null_space_lu     !! `sparse_qp%null_space`
        integer :: dense_max_ns = -1                       !! `sparse_qp%dense_max_ns` (`< 0`: the default)
    end type config_type

    ! problem data (for the user functions):
    integer  :: nsteps                           !! `control`: number of steps `N`
    integer  :: nctrl                            !! `control`: number of controls `K`
    real(wp) :: h                                !! `control`: step size
    integer  :: ngrid                            !! `elliptic`: grid size `N`
    real(wp), parameter :: alpha = 0.01_wp       !! `elliptic`: weight of the controls
    real(wp), dimension(:), allocatable :: yd    !! `elliptic`: target state `dimension(N^2)`
    real(wp), dimension(:), allocatable :: jac0  !! `elliptic`: the constant Jacobian values
    integer,  dimension(:), allocatable :: diag_pos !! `elliptic`: position of each row's `y_k` entry in the Jacobian
    real(wp), dimension(:,:), allocatable :: amat   !! `dense`: the linear part of the constraints `dimension(m,n)`

    type(config_type), dimension(*), parameter :: qp_configs = [ &
        config_type(name='dense QP', qp_mode=sqpopt_qp_dense), &
        config_type(name='sparse QP (LU)'), &
        config_type(name='sparse QP (LSQR)', null_space=sqpopt_null_space_lsqr) ]

    type(config_type), dimension(*), parameter :: dof_configs = [ &
        config_type(name='default (L-BFGS, 100 pairs)'), &
        config_type(name='L-BFGS, 10 pairs', lbfgs_memory=10), &
        config_type(name='L-BFGS, 500 pairs', lbfgs_memory=500), &
        config_type(name='dense_max_ns = 0', dense_max_ns=0), &
        config_type(name='dense_max_ns = 1000', dense_max_ns=1000), &
        config_type(name='exact Hessian', hessian_mode=sqpopt_hessian_exact) ]

    integer, dimension(*), parameter :: control_sizes  = [25, 50, 100, 200, 400, 800] !! `N` (with `K = N`)
    integer, dimension(*), parameter :: elliptic_sizes = [5, 7, 10, 14, 20, 28]       !! grid sizes
    integer, dimension(*), parameter :: dof_sizes = [5, 20, 100, 400, 1000]           !! `K` (with `N = 1000`)
    integer, dimension(*), parameter :: dense_sizes = [20, 50, 100, 200, 400]         !! `n` of the `dense` problem

    type(sqpopt_problem_type) :: problem
    real(wp), dimension(:), allocatable :: x0
    integer :: i, j

    write(*,'(/A)') '1. QP solver against problem size (many degrees of freedom)'
    call header()
    do i = 1, size(control_sizes)
        call setup_control(control_sizes(i), control_sizes(i), problem, x0)
        do j = 1, size(qp_configs)
            call run('control', qp_configs(j), problem, x0)
        end do
    end do
    do i = 1, size(elliptic_sizes)
        call setup_elliptic(elliptic_sizes(i), problem, x0)
        do j = 1, size(qp_configs)
            call run('elliptic', qp_configs(j), problem, x0)
        end do
    end do

    write(*,'(/A)') '2. Degrees of freedom: control, N = 1000 steps, K controls'
    call header()
    do i = 1, size(dof_sizes)
        call setup_control(1000, dof_sizes(i), problem, x0)
        do j = 1, size(dof_configs)
            ! (the dense reduced Hessian only differs from the default with more than 50 superbasics)
            if (dof_configs(j)%dense_max_ns > 0 .and. dof_sizes(i) <= 50) cycle
            call run('control', dof_configs(j), problem, x0)
        end do
    end do

    write(*,'(/A)') '3. QP solver against problem size (dense Jacobian)'
    call header()
    do i = 1, size(dense_sizes)
        call setup_dense(dense_sizes(i), problem, x0)
        do j = 1, 2
            call run('dense', qp_configs(j), problem, x0)
        end do
    end do

    contains

    subroutine header()
    !! print the column headings of a study's table
    write(*,'(A10,A6,A6,A5,1X,A28,A6,A6,A7,A9,A10,A10,A16)') 'problem', 'n', 'm', 'dof', 'configuration', &
        'istat', 'iter', 'fc', 'qp_iter', 'time (s)', 'QP (s)', 'f'
    end subroutine header

    subroutine run(name, cfg, problem, x0)
    !! solve a problem with one configuration, and print one line of results
    character(len=*),          intent(in) :: name    !! the problem, for the output
    type(config_type),         intent(in) :: cfg     !! the configuration
    type(sqpopt_problem_type), intent(in) :: problem !! the problem
    real(wp), dimension(:),    intent(in) :: x0      !! starting point `dimension(n)`
    type(sqpopt_type)           :: solver
    type(sqpopt_options_type)   :: options
    type(sqpopt_qp_solver_type) :: qp
    type(sqpopt_results_type)   :: r
    integer :: istat

    if ((cfg%qp_mode == sqpopt_qp_dense .and. problem%n > dense_max_n) .or. &
        (cfg%null_space == sqpopt_null_space_lsqr .and. problem%n > lsqr_max_n)) then
        write(*,'(A10,I6,I6,I5,1X,A28,A)') name, problem%n, problem%m, problem%n - problem%m, cfg%name, &
            '  (skipped: too slow)'
        return
    end if

    options%max_iter       = 3000
    options%qp_solver_mode = cfg%qp_mode
    options%hessian_mode   = cfg%hessian_mode
    options%lbfgs_memory   = cfg%lbfgs_memory
    qp%sparse_qp%null_space = cfg%null_space
    if (cfg%dense_max_ns >= 0) qp%sparse_qp%dense_max_ns = cfg%dense_max_ns

    call solver%initialize(problem=problem, options=options, qp_solver=qp)
    call solver%solve(x0, istat)
    call solver%get_results(r)

    write(*,'(A10,I6,I6,I5,1X,A28,I6,I6,I7,I9,F10.3,F10.3,ES16.8)') name, problem%n, problem%m, &
        problem%n - problem%m, cfg%name, r%istat, r%iterations, r%n_eval_fc, r%n_qp_iterations, &
        r%time, r%time_qp, r%f
    flush(6)

    end subroutine run

    !------------------------------------------------------------------------
    ! control: N steps, K piecewise-constant controls
    !------------------------------------------------------------------------

    pure integer function ctrl(k)
    !! the column of the control that applies on step `k` (`0..N-1`)
    integer, intent(in) :: k !! the step
    ctrl = nsteps + 2 + (k*nctrl)/nsteps
    end function ctrl

    subroutine setup_control(nn, kk, problem, x0)
    !! set up the control problem
    integer,                   intent(in)  :: nn      !! number of steps `N`
    integer,                   intent(in)  :: kk      !! number of controls `K` (a divisor of `N`)
    type(sqpopt_problem_type), intent(out) :: problem !! the problem
    real(wp), dimension(:), allocatable, intent(out) :: x0 !! starting point `dimension(n)`
    real(wp), dimension(:), allocatable :: x_lb, x_ub
    integer, dimension(:), allocatable :: irow, icol
    integer :: n, m, k, nnz
    nsteps = nn
    nctrl  = kk
    h = 5.0_wp/nsteps
    n = nsteps + 1 + nctrl
    m = nsteps + 1

    ! variables: y_0..y_N are 1..N+1, u_1..u_K are N+2..N+1+K
    allocate(x_lb(n), x_ub(n))
    x_lb = -sqpopt_infinity; x_ub = sqpopt_infinity
    x_lb(nsteps+2:n) = -0.3_wp; x_ub(nsteps+2:n) = 0.3_wp

    ! Jacobian pattern: row 1 = y_0; row k+2 = (y_{k+1}, y_k, u_{j(k)})
    nnz = 1 + 3*nsteps
    allocate(irow(nnz), icol(nnz))
    irow(1) = 1; icol(1) = 1
    do k = 0, nsteps-1
        irow(2+3*k:4+3*k) = k+2
        icol(2+3*k) = k+2
        icol(3+3*k) = k+1
        icol(4+3*k) = ctrl(k)
    end do

    call problem%set_problem_size(n=n, m=m)
    call problem%set_bounds(x_lb, x_ub, [1.0_wp, spread(0.0_wp,1,m-1)], [1.0_wp, spread(0.0_wp,1,m-1)])
    call problem%set_jacobian_sparsity(nnz, irow, icol)
    call problem%set_functions(fc=fc_control, gjac=gjac_control, hess=hess_control)
    call problem%set_hessian_sparsity(n, [(k, k=1,n)], [(k, k=1,n)])

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
    f = 0.5_wp*h*sum(x(1:nsteps)**2) + x(nsteps+1)**2 + 0.5_wp*h*(real(nsteps, wp)/nctrl)*sum(x(nsteps+2:)**2)
    c(1) = x(1)
    do k = 0, nsteps-1
        c(k+2) = x(k+2) - x(k+1) - h*(x(ctrl(k)) - x(k+1)**3)
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
    g(1:nsteps)  = h*x(1:nsteps)
    g(nsteps+1)  = 2.0_wp*x(nsteps+1)
    g(nsteps+2:) = h*(real(nsteps, wp)/nctrl)*x(nsteps+2:)
    jac(1) = 1.0_wp
    do k = 0, nsteps-1
        jac(2+3*k) = 1.0_wp
        jac(3+3*k) = -1.0_wp + 3.0_wp*h*x(k+1)**2
        jac(4+3*k) = -h
    end do
    end subroutine gjac_control

    subroutine hess_control(x, lambda, hess_val, status, data)
    !! the (diagonal) Hessian of the Lagrangian of the control problem
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: k
    hess_val(1:nsteps)  = h
    hess_val(nsteps+1)  = 2.0_wp
    hess_val(nsteps+2:) = h*real(nsteps, wp)/nctrl
    do k = 0, nsteps-1
        hess_val(k+1) = hess_val(k+1) - lambda(k+2)*6.0_wp*h*x(k+1)
    end do
    end subroutine hess_control

    !------------------------------------------------------------------------
    ! elliptic: semilinear elliptic control on an N x N grid
    !------------------------------------------------------------------------

    pure function laplacian(v) result(av)
    !! the 5-point Laplacian stencil (without the `1/h^2`) on the `N x N` grid,
    !! with zero boundary values
    real(wp), dimension(:), intent(in) :: v  !! values at the grid points, by columns `dimension(N^2)`
    real(wp), dimension(size(v)) :: av
    integer :: i, j, k
    do j = 1, ngrid
        do i = 1, ngrid
            k = (j-1)*ngrid + i
            av(k) = 4.0_wp*v(k)
            if (i > 1)     av(k) = av(k) - v(k-1)
            if (i < ngrid) av(k) = av(k) - v(k+1)
            if (j > 1)     av(k) = av(k) - v(k-ngrid)
            if (j < ngrid) av(k) = av(k) - v(k+ngrid)
        end do
    end do
    end function laplacian

    subroutine setup_elliptic(nn, problem, x0)
    !! set up the elliptic control problem (its data manufactured from a
    !! chosen solution, with one in seven of the control bounds active)
    integer,                   intent(in)  :: nn      !! grid size `N` (`n = 2N^2` variables)
    type(sqpopt_problem_type), intent(out) :: problem !! the problem
    real(wp), dimension(:), allocatable, intent(out) :: x0 !! starting point `dimension(n)`
    real(wp), dimension(:), allocatable :: y, u, lam, u_lb, u_ub
    integer,  dimension(:), allocatable :: irow, icol
    integer :: ny, n, nnz, nk, i, j, k
    logical, dimension(6) :: keep
    ngrid = nn
    ny = ngrid**2
    n  = 2*ny

    allocate(y(ny), lam(ny), u_lb(ny), u_ub(ny))
    do j = 1, ngrid
        do i = 1, ngrid
            k = (j-1)*ngrid + i
            y(k) = 0.5_wp*sin(2.0_wp*pi*i/(ngrid+1))*cos(3.0_wp*pi*j/(ngrid+1)) + 0.2_wp*cos(real(k, wp))
        end do
    end do
    u = laplacian(y) + y**3
    lam = alpha*u
    u_lb = -10.0_wp
    u_ub =  10.0_wp
    do k = 1, ny
        select case (mod(k, 7))
        case (0)
            u_ub(k) = u(k)
            lam(k)  = lam(k) + 0.05_wp
        case (3)
            u_lb(k) = u(k)
            lam(k)  = lam(k) - 0.05_wp
        case default
        end select
    end do
    yd = y + laplacian(lam) + 3.0_wp*y**2*lam

    nnz = 2*ny + 4*ngrid*(ngrid-1)
    if (allocated(jac0)) deallocate(jac0, diag_pos)
    allocate(irow(nnz), icol(nnz), jac0(nnz), diag_pos(ny))
    jac0 = -1.0_wp
    nnz = 0
    do j = 1, ngrid
        do i = 1, ngrid
            k = (j-1)*ngrid + i
            keep = [.true., i > 1, i < ngrid, j > 1, j < ngrid, .true.]
            nk = count(keep)
            irow(nnz+1:nnz+nk) = k
            icol(nnz+1:nnz+nk) = pack([k, k-1, k+1, k-ngrid, k+ngrid, ny+k], keep)
            diag_pos(k) = nnz + 1
            nnz = nnz + nk
        end do
    end do

    call problem%set_problem_size(n=n, m=ny)
    call problem%set_bounds([spread(-sqpopt_infinity, 1, ny), u_lb], [spread(sqpopt_infinity, 1, ny), u_ub], &
                            spread(0.0_wp, 1, ny), spread(0.0_wp, 1, ny))
    call problem%set_jacobian_sparsity(nnz, irow, icol)
    call problem%set_functions(fc=fc_elliptic, gjac=gjac_elliptic, hess=hess_elliptic)
    call problem%set_hessian_sparsity(n, [(k, k=1,n)], [(k, k=1,n)])

    allocate(x0(n))
    x0 = 0.0_wp

    end subroutine setup_elliptic

    subroutine fc_elliptic(x, f, c, status, data)
    !! objective and constraints of the elliptic control problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: ny
    ny = ngrid**2
    f = 0.5_wp*sum((x(1:ny) - yd)**2) + 0.5_wp*alpha*sum(x(ny+1:)**2)
    c = laplacian(x(1:ny)) + x(1:ny)**3 - x(ny+1:)
    end subroutine fc_elliptic

    subroutine gjac_elliptic(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian values of the elliptic control problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: ny
    ny = ngrid**2
    g(1:ny)  = x(1:ny) - yd
    g(ny+1:) = alpha*x(ny+1:)
    jac = jac0
    jac(diag_pos) = 4.0_wp + 3.0_wp*x(1:ny)**2
    end subroutine gjac_elliptic

    subroutine hess_elliptic(x, lambda, hess_val, status, data)
    !! the (diagonal) Hessian of the Lagrangian of the elliptic control problem
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: ny
    ny = ngrid**2
    hess_val(1:ny)  = 1.0_wp - 6.0_wp*lambda*x(1:ny)
    hess_val(ny+1:) = alpha
    end subroutine hess_elliptic

    !------------------------------------------------------------------------
    ! dense: every constraint depends on every variable
    !------------------------------------------------------------------------

    subroutine setup_dense(n, problem, x0)
    !! set up the dense problem: `n` variables, `m = n/2` constraints (the odd
    !! ones equalities, the even ones inequalities), with a dense Jacobian:
    !!
    !!     minimize   sum_j (1 + mod(j,5))*(x_j - 1)^2/2
    !!     subject to sum_j a_ij*x_j + 0.1*x_i^2 (= or <=) b_i,   -2 <= x <= 2
    !!
    !! with `a_ij = sin(0.7*i*j)/sqrt(n) + 2*delta_ij`, and `b` such that
    !! `x_j = cos(j)/2` is feasible
    integer,                   intent(in)  :: n       !! number of variables (even)
    type(sqpopt_problem_type), intent(out) :: problem !! the problem
    real(wp), dimension(:), allocatable, intent(out) :: x0 !! starting point `dimension(n)`
    real(wp), dimension(:), allocatable :: xf, b, c_lb
    integer :: m, i, j
    m = n/2
    if (allocated(amat)) deallocate(amat)
    allocate(amat(m,n))
    do j = 1, n
        do i = 1, m
            amat(i,j) = sin(0.7_wp*i*j)/sqrt(real(n, wp))
        end do
    end do
    do i = 1, m
        amat(i,i) = amat(i,i) + 2.0_wp
    end do
    xf = [(0.5_wp*cos(real(j, wp)), j=1,n)]
    b = matmul(amat, xf) + 0.1_wp*xf(1:m)**2
    c_lb = b
    c_lb(2::2) = -sqpopt_infinity

    call problem%set_problem_size(n=n, m=m)
    call problem%set_bounds(spread(-2.0_wp, 1, n), spread(2.0_wp, 1, n), c_lb, b)
    call problem%set_jacobian_sparsity(m*n, [((i, j=1,n), i=1,m)], [((j, j=1,n), i=1,m)])
    call problem%set_functions(fc=fc_dense, gjac=gjac_dense, hess=hess_dense)
    call problem%set_hessian_sparsity(n, [(j, j=1,n)], [(j, j=1,n)])

    allocate(x0(n))
    x0 = 0.0_wp

    end subroutine setup_dense

    subroutine fc_dense(x, f, c, status, data)
    !! objective and constraints of the dense problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: j
    f = 0.5_wp*sum([(1 + mod(j,5), j=1,size(x))]*(x - 1.0_wp)**2)
    c = matmul(amat, x) + 0.1_wp*x(1:size(c))**2
    end subroutine fc_dense

    subroutine gjac_dense(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian values of the dense problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n, i, j
    n = size(x)
    g = [(1 + mod(j,5), j=1,n)]*(x - 1.0_wp)
    jac = reshape(transpose(amat), [size(jac)])
    do i = 1, size(amat, 1)
        jac((i-1)*n + i) = jac((i-1)*n + i) + 0.2_wp*x(i)
    end do
    end subroutine gjac_dense

    subroutine hess_dense(x, lambda, hess_val, status, data)
    !! the (diagonal) Hessian of the Lagrangian of the dense problem
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: j
    hess_val = [(real(1 + mod(j,5), wp), j=1,size(x))]
    hess_val(1:size(lambda)) = hess_val(1:size(lambda)) - 0.2_wp*lambda
    end subroutine hess_dense

end program settings_study
