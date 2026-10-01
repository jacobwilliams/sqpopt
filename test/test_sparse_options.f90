program test_sparse_options

    !! Large sparse problems, each solved with every combination of the
    !! options that matter for the sparse QP (`sqpopt_qp_reduced_hessian`,
    !! which `sqpopt_qp_auto` picks at these sizes) and the Hessian:
    !!
    !! * the null-space method: `sqpopt_null_space_lu` (the default) or
    !!   `sqpopt_null_space_lsqr`;
    !! * `dense_max_ns = 250` (a dense Cholesky factorization of the reduced
    !!   Hessian for faces with up to 250 superbasics, instead of CG: with
    !!   the default 50, these problems' faces are all solved by CG);
    !! * `warm_start = .false.` (every QP from a crash start);
    !! * the Hessian: L-BFGS (default and short memory), L-SR1, or the
    !!   user's exact sparse Hessian.
    !!
    !! (The sizes, about 300-1000 variables, keep the `LSQR` solves, whose
    !! cost grows quickly with `n`, to a few seconds.)
    !!
    !! Every problem is manufactured from a chosen solution and multipliers
    !! (with strict complementarity and a positive definite Hessian of the
    !! Lagrangian), so the solution `x*` and the optimal objective are known
    !! exactly, and every solve must reach them:
    !!
    !! * `elliptic`: a semilinear elliptic control problem on an `N x N` grid
    !!   (`n = 2N^2`, `m = N^2` equality constraints with the 5-point
    !!   Laplacian stencil, whose basis LU factors fill in), with bounds on
    !!   the controls, one in seven of them active at the solution:
    !!
    !!       minimize   sum (y - y_d)^2/2 + alpha*sum u^2/2
    !!       subject to A*y + y^3 - u = 0,   u_l <= u <= u_u
    !!
    !! * `ball`: a linear objective on a ball, with `n` variables in one dense
    !!   constraint row and a third of the bounds active (`n-1` superbasics
    !!   at most, so the reduced Hessian is never small):
    !!
    !!       minimize   w^T x
    !!       subject to sum x^2 <= r^2,   -1 <= x <= 1
    !!
    !! * `chain`: a projection onto a chain of `n-1` convex inequality
    !!   constraints, two in five of them active at the solution, and one in
    !!   five of the lower bounds, from a start that violates every
    !!   constraint (so the first QPs have `n-1` elastic slacks):
    !!
    !!       minimize   sum (x - a)^2/2
    !!       subject to x_i^2 + x_{i+1}^2 <= 1,   i = 1..n-1
    !!                  0 <= x <= 2

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_bfgs, sqpopt_hessian_sr1, sqpopt_hessian_exact
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_null_space_lu, sqpopt_null_space_lsqr
    use sqpopt_types_module,     only: sqpopt_results_type, sqpopt_success, sqpopt_acceptable, sqpopt_stalled, &
                                       sqpopt_infinity
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: rel_tol  = 1.0e-6_wp !! relative tolerance on the optimal objective
    real(wp), parameter :: x_tol    = 1.0e-4_wp !! tolerance on the distance to `x*` (relative to \( \max(1,|x^*|_\infty) \))
    real(wp), parameter :: feas_tol = 1.0e-6_wp !! tolerance on the constraint violation
    real(wp), parameter :: pi = acos(-1.0_wp)

    type :: config_type
        !! a configuration of the sparse QP and the Hessian
        character(len=24) :: name = ''                  !! label, for the output
        integer :: hessian_mode = sqpopt_hessian_bfgs   !! `options%hessian_mode`
        integer :: lbfgs_memory = 0                     !! `options%lbfgs_memory` (`0`: the default)
        integer :: null_space   = sqpopt_null_space_lu  !! `sparse_qp%null_space`
        integer :: dense_max_ns = -1                    !! `sparse_qp%dense_max_ns` (`< 0`: the default)
        logical :: warm_start   = .true.                !! `sparse_qp%warm_start`
    end type config_type

    type(config_type), dimension(*), parameter :: configs = [ &
        config_type(name='default'), &
        config_type(name='LU, dense reduced Hess.', dense_max_ns=250), &
        config_type(name='LU, cold start', warm_start=.false.), &
        config_type(name='LSQR', null_space=sqpopt_null_space_lsqr), &
        config_type(name='L-BFGS memory 5', lbfgs_memory=5), &
        config_type(name='L-SR1', hessian_mode=sqpopt_hessian_sr1), &
        config_type(name='exact Hessian', hessian_mode=sqpopt_hessian_exact), &
        config_type(name='exact Hess., dense red.', hessian_mode=sqpopt_hessian_exact, dense_max_ns=250), &
        config_type(name='exact Hessian, LSQR', hessian_mode=sqpopt_hessian_exact, null_space=sqpopt_null_space_lsqr) ]

    ! problem data (for the user functions):
    integer  :: ngrid                            !! `elliptic`: grid size `N`
    real(wp), parameter :: alpha = 0.01_wp       !! `elliptic`: weight of the controls
    real(wp), dimension(:), allocatable :: yd    !! `elliptic`: target state `dimension(N^2)`
    real(wp), dimension(:), allocatable :: jac0  !! `elliptic`: the constant Jacobian values
    integer,  dimension(:), allocatable :: diag_pos !! `elliptic`: position of each row's `y_k` entry in the Jacobian
    real(wp), dimension(:), allocatable :: w     !! `ball`: objective coefficients `dimension(n)`
    real(wp), dimension(:), allocatable :: a     !! `chain`: the point projected `dimension(n)`

    type(sqpopt_problem_type) :: problem
    real(wp), dimension(:), allocatable :: x0, x_star
    real(wp) :: f_star
    integer  :: i

    write(*,*) '----------------------------'
    write(*,*) 'test_sparse_options'
    write(*,*) '----------------------------'

    call setup_elliptic(16, problem, x0, x_star, f_star)
    do i = 1, size(configs)
        call run('elliptic', configs(i), problem, x0, x_star, f_star)
    end do

    call setup_ball(1002, problem, x0, x_star, f_star)
    do i = 1, size(configs)
        call run('ball', configs(i), problem, x0, x_star, f_star)
    end do

    call setup_chain(300, problem, x0, x_star, f_star)
    do i = 1, size(configs)
        call run('chain', configs(i), problem, x0, x_star, f_star)
    end do

    print '(A)', 'test_sparse_options PASSED'

    contains

    subroutine run(name, cfg, problem, x0, x_star, f_star)
    !! solve a problem with one configuration, and check the solution
    character(len=*),          intent(in) :: name    !! the problem, for the output
    type(config_type),         intent(in) :: cfg     !! the configuration
    type(sqpopt_problem_type), intent(in) :: problem !! the problem
    real(wp), dimension(:),    intent(in) :: x0      !! starting point `dimension(n)`
    real(wp), dimension(:),    intent(in) :: x_star  !! the known solution `dimension(n)`
    real(wp),                  intent(in) :: f_star  !! the known optimal objective
    type(sqpopt_type)           :: solver
    type(sqpopt_options_type)   :: options
    type(sqpopt_qp_solver_type) :: qp
    type(sqpopt_results_type)   :: r
    character(len=:), allocatable :: label
    real(wp) :: x_err
    integer  :: istat

    options%max_iter     = 3000
    options%hessian_mode = cfg%hessian_mode
    options%lbfgs_memory = cfg%lbfgs_memory
    qp%sparse_qp%null_space = cfg%null_space
    qp%sparse_qp%warm_start = cfg%warm_start
    if (cfg%dense_max_ns >= 0) qp%sparse_qp%dense_max_ns = cfg%dense_max_ns

    call solver%initialize(problem=problem, options=options, qp_solver=qp)
    call solver%solve(x0, istat)
    call solver%get_results(r)

    label = name//' ('//trim(cfg%name)//')'
    x_err = maxval(abs(r%x - x_star))/max(1.0_wp, maxval(abs(x_star)))
    print '(A36,A,I6,A,I3,A,ES16.8,A,ES9.2,A,ES9.2,A,I5,A,I5,A,I7,A,F7.3,A)', label, ': n=', problem%n, &
        ' istat=', r%istat, ' f=', r%f, ' |x-x*|=', x_err, ' viol=', r%feasibility_error, &
        ' iter=', r%iterations, ' fc=', r%n_eval_fc, ' qp_iter=', r%n_qp_iterations, ' time=', r%time, ' s'

    if (problem%n <= qp%auto_dense_max_n) error stop 'test_sparse_options FAILED: problem too small for the sparse QP'
    if (cfg%hessian_mode == sqpopt_hessian_exact .and. r%n_eval_hess == 0) then
        error stop 'test_sparse_options FAILED: '//label//': the Hessian was not used'
    end if
    if (r%istat /= sqpopt_success .and. r%istat /= sqpopt_acceptable .and. r%istat /= sqpopt_stalled) then
        error stop 'test_sparse_options FAILED: '//label//': did not converge'
    end if
    if (abs(r%f - f_star) > rel_tol*max(1.0_wp, abs(f_star))) then
        error stop 'test_sparse_options FAILED: '//label//': wrong objective'
    end if
    if (x_err > x_tol) error stop 'test_sparse_options FAILED: '//label//': wrong solution'
    if (r%feasibility_error > feas_tol) error stop 'test_sparse_options FAILED: '//label//': infeasible'

    end subroutine run

    !------------------------------------------------------------------------
    ! elliptic: semilinear elliptic control on an N x N grid
    !------------------------------------------------------------------------

    pure function laplacian(v) result(av)
    !! the 5-point Laplacian stencil (without the `1/h^2`) on the `N x N` grid,
    !! with zero boundary values: `4*v_k` minus the (up to four) neighbours
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

    subroutine setup_elliptic(nn, problem, x0, x_star, f_star)
    !! set up the elliptic control problem, from its manufactured solution
    integer,                   intent(in)  :: nn      !! grid size `N` (`n = 2N^2` variables)
    type(sqpopt_problem_type), intent(out) :: problem !! the problem
    real(wp), dimension(:), allocatable, intent(out) :: x0     !! starting point `dimension(n)`
    real(wp), dimension(:), allocatable, intent(out) :: x_star !! the solution `dimension(n)`
    real(wp),                  intent(out) :: f_star  !! the optimal objective
    real(wp), dimension(:), allocatable :: y, u, lam, u_lb, u_ub
    integer,  dimension(:), allocatable :: irow, icol
    integer :: ny, n, nnz, nk, i, j, k
    logical, dimension(6) :: keep
    ngrid = nn
    ny = ngrid**2
    n  = 2*ny

    ! the solution: a smooth state plus a rough part, its control, and
    ! multipliers (`L = f + lambda^T c`), with the control at its upper bound
    ! where mod(k,7) = 0 and at its lower bound where mod(k,7) = 3:
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
    ! (stationarity in `y`: y - y_d + (A + 3*diag(y^2))*lambda = 0)
    yd = y + laplacian(lam) + 3.0_wp*y**2*lam
    x_star = [y, u]
    f_star = 0.5_wp*sum((y - yd)**2) + 0.5_wp*alpha*sum(u**2)

    ! Jacobian pattern, row by row: y_k, its neighbours, then u_k
    nnz = 2*ny + 4*ngrid*(ngrid-1)
    allocate(irow(nnz), icol(nnz), jac0(nnz), diag_pos(ny))
    jac0 = -1.0_wp
    nnz = 0
    do j = 1, ngrid
        do i = 1, ngrid
            k = (j-1)*ngrid + i
            ! (only the neighbours inside the grid)
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
    ! ball: linear objective on a ball, with one dense constraint row
    !------------------------------------------------------------------------

    subroutine setup_ball(n, problem, x0, x_star, f_star)
    !! set up the ball problem, from its manufactured solution
    integer,                   intent(in)  :: n       !! number of variables (a multiple of 6)
    type(sqpopt_problem_type), intent(out) :: problem !! the problem
    real(wp), dimension(:), allocatable, intent(out) :: x0     !! starting point `dimension(n)`
    real(wp), dimension(:), allocatable, intent(out) :: x_star !! the solution `dimension(n)`
    real(wp),                  intent(out) :: f_star  !! the optimal objective
    real(wp), dimension(6), parameter :: w_cycle = [3.0_wp, -1.0_wp, 0.5_wp, -2.5_wp, 1.5_wp, -0.25_wp]
    integer :: i

    ! with the constraint's multiplier 1 (`w + 2*x = 0` where no bound is
    ! active), the solution is x = -w/2 clipped to [-1,1] (strictly, since
    ! no |w_i| is 2), and the radius is whatever that point has:
    w = [(w_cycle(mod(i-1, 6) + 1), i=1,n)]
    x_star = min(max(-0.5_wp*w, -1.0_wp), 1.0_wp)
    f_star = dot_product(w, x_star)

    call problem%set_problem_size(n=n, m=1)
    call problem%set_bounds(spread(-1.0_wp, 1, n), spread(1.0_wp, 1, n), [-sqpopt_infinity], [sum(x_star**2)])
    call problem%set_jacobian_sparsity(n, spread(1, 1, n), [(i, i=1,n)])
    call problem%set_functions(fc=fc_ball, gjac=gjac_ball, hess=hess_ball)
    call problem%set_hessian_sparsity(n, [(i, i=1,n)], [(i, i=1,n)])

    allocate(x0(n))
    x0 = 0.1_wp

    end subroutine setup_ball

    subroutine fc_ball(x, f, c, status, data)
    !! objective and constraint of the ball problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = dot_product(w, x)
    c(1) = sum(x**2)
    end subroutine fc_ball

    subroutine gjac_ball(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian values of the ball problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    g = w
    jac = 2.0_wp*x
    end subroutine gjac_ball

    subroutine hess_ball(x, lambda, hess_val, status, data)
    !! the (diagonal) Hessian of the Lagrangian of the ball problem
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    hess_val = -2.0_wp*lambda(1)
    end subroutine hess_ball

    !------------------------------------------------------------------------
    ! chain: projection onto a chain of convex inequality constraints
    !------------------------------------------------------------------------

    subroutine setup_chain(n, problem, x0, x_star, f_star)
    !! set up the chain problem, from its manufactured solution
    integer,                   intent(in)  :: n       !! number of variables (a multiple of 5)
    type(sqpopt_problem_type), intent(out) :: problem !! the problem
    real(wp), dimension(:), allocatable, intent(out) :: x0     !! starting point `dimension(n)`
    real(wp), dimension(:), allocatable, intent(out) :: x_star !! the solution `dimension(n)`
    real(wp),                  intent(out) :: f_star  !! the optimal objective
    real(wp), dimension(5), parameter :: x_cycle = [0.6_wp, 0.8_wp, 0.6_wp, 0.2_wp, 0.0_wp]
    real(wp), dimension(:), allocatable :: mu
    integer, dimension(:), allocatable :: irow, icol
    integer :: m, i

    ! the solution: rows 5k+1 and 5k+2 active (0.6^2 + 0.8^2 = 1), with
    ! multipliers `mu >= 0` (`L = f + mu^T c`), and every fifth variable at
    ! its lower bound 0, with bound multiplier 0.5:
    m = n - 1
    x_star = [(x_cycle(mod(i-1, 5) + 1), i=1,n)]
    allocate(mu(0:n))
    mu = 0.0_wp
    do i = 1, m
        if (mod(i, 5) == 1 .or. mod(i, 5) == 2) mu(i) = 0.5_wp + 0.25_wp*mod(i, 3)
    end do
    ! (stationarity: x - a + J^T*mu - z = 0)
    a = x_star + 2.0_wp*x_star*(mu(0:n-1) + mu(1:n))
    where (x_star == 0.0_wp) a = a - 0.5_wp
    f_star = 0.5_wp*sum((x_star - a)**2)

    allocate(irow(2*m), icol(2*m))
    do i = 1, m
        irow(2*i-1:2*i) = i
        icol(2*i-1) = i
        icol(2*i)   = i + 1
    end do
    call problem%set_problem_size(n=n, m=m)
    call problem%set_bounds(spread(0.0_wp, 1, n), spread(2.0_wp, 1, n), spread(-sqpopt_infinity, 1, m), &
                            spread(1.0_wp, 1, m))
    call problem%set_jacobian_sparsity(2*m, irow, icol)
    call problem%set_functions(fc=fc_chain, gjac=gjac_chain, hess=hess_chain)
    call problem%set_hessian_sparsity(n, [(i, i=1,n)], [(i, i=1,n)])

    allocate(x0(n))
    x0 = 1.0_wp

    end subroutine setup_chain

    subroutine fc_chain(x, f, c, status, data)
    !! objective and constraints of the chain problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n
    n = size(x)
    f = 0.5_wp*sum((x - a)**2)
    c = x(1:n-1)**2 + x(2:n)**2
    end subroutine fc_chain

    subroutine gjac_chain(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian values of the chain problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: n
    n = size(x)
    g = x - a
    jac(1::2) = 2.0_wp*x(1:n-1)
    jac(2::2) = 2.0_wp*x(2:n)
    end subroutine gjac_chain

    subroutine hess_chain(x, lambda, hess_val, status, data)
    !! the (diagonal) Hessian of the Lagrangian of the chain problem
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    hess_val = 1.0_wp
    hess_val(1:size(lambda))   = hess_val(1:size(lambda))   - 2.0_wp*lambda
    hess_val(2:size(lambda)+1) = hess_val(2:size(lambda)+1) - 2.0_wp*lambda
    end subroutine hess_chain

end program test_sparse_options
