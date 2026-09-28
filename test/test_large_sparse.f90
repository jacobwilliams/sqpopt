program test_large_sparse

    !! Large sparse problems, to test the solver with the sparse QP
    !! (`sqpopt_qp_reduced_hessian`, the default for problems with more than
    !! `auto_dense_max_n` variables) at scale. Each large problem is solved
    !! with the default options (`sqpopt_qp_auto`, which picks the sparse QP
    !! here) and with the sparse QP forced. For comparison, a small instance
    !! of each (`n` about 50; the dense QP is impractical at the large sizes)
    !! is solved with the dense QP (`sqpopt_qp_dense`) and with the sparse
    !! QP. Finally, the large problems are solved with the user's exact
    !! (sparse) Hessian of the Lagrangian (`sqpopt_hessian_exact`) instead of
    !! the quasi-Newton approximation. Every solve must converge to the known
    !! optimum:
    !!
    !! * `control`: a discretized nonlinear optimal-control problem with
    !!   `N = 500` steps (`n = 2N+1 = 1001` variables, `m = N+1 = 501`
    !!   equality constraints, a banded Jacobian with 3 nonzeros per row, and
    !!   bounds on the controls; small size: `N = 25`, `n = 51`):
    !!
    !!       minimize   sum_k h*(y_k^2 + u_k^2)/2 + (y_N)^2
    !!       subject to y_0 = 1,
    !!                  y_{k+1} = y_k + h*(u_k - y_k^3),   k = 0..N-1
    !!                  -0.3 <= u_k <= 0.3
    !!
    !! * `rosenbrock`: the chained Rosenbrock function with `n = 2000`
    !!   variables and `n/2 = 1000` circle constraints (2 nonzeros per row),
    !!   and bounds (small size: `n = 50`):
    !!
    !!       minimize   sum_i 100*(x_{i+1}-x_i^2)^2 + (1-x_i)^2
    !!       subject to x_i^2 + x_{i+1}^2 <= 1.5,   i = 1,3,5,...
    !!                  -2 <= x <= 2
    !!
    !! (the same problems as `example/benchmark.f90`).

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type, sqpopt_hessian_exact
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type, sqpopt_qp_auto, sqpopt_qp_reduced_hessian, &
                                       sqpopt_qp_dense
    use sqpopt_types_module,     only: sqpopt_results_type, sqpopt_success, sqpopt_acceptable, sqpopt_stalled
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: rel_tol = 1.0e-5_wp  !! relative tolerance on the optimal objective
    real(wp), parameter :: feas_tol = 1.0e-6_wp !! tolerance on the constraint violation

    integer  :: nsteps   !! `control`: number of steps `N`
    real(wp) :: h        !! `control`: step size
    integer  :: mode
    integer, parameter :: large_modes(2) = [sqpopt_qp_auto,  sqpopt_qp_reduced_hessian] !! QP modes for the large problems
    integer, parameter :: small_modes(2) = [sqpopt_qp_dense, sqpopt_qp_reduced_hessian] !! QP modes for the small ones

    write(*,*) '----------------------------'
    write(*,*) 'test_large_sparse'
    write(*,*) '----------------------------'

    ! large problems: the default mode (which picks the sparse QP) and the sparse QP forced
    do mode = 1, size(large_modes)
        call run_control(500, large_modes(mode), 3.30614092e-1_wp)
        call run_rosenbrock(2000, large_modes(mode), 1.97466586e3_wp)
    end do

    ! small problems: the dense QP vs. the sparse QP
    do mode = 1, size(small_modes)
        call run_control(25, small_modes(mode), 3.36087467e-1_wp)
        call run_rosenbrock(50, small_modes(mode), 4.43669595e1_wp)
    end do

    ! large problems with the exact Hessian (and the default QP mode):
    call run_control(500, sqpopt_qp_auto, 3.30614092e-1_wp, exact=.true.)
    ! (the chained Rosenbrock problem is nonconvex, with many local solutions:
    ! the Newton iterations converge, in far fewer iterations, to a different
    ! one than the quasi-Newton ones)
    call run_rosenbrock(2000, sqpopt_qp_auto, 1.97826254e3_wp, exact=.true.)

    print '(A)', 'test_large_sparse PASSED'

    contains

    subroutine check(name, qp_mode, n, solver, f_star, exact)
    !! check the solve: converged (or stopped at an acceptable/stalled point),
    !! to the known optimum, and feasible (and, in the `sqpopt_qp_auto` mode,
    !! with the sparse QP)
    character(len=*),  intent(in) :: name    !! the problem, for the output
    integer,           intent(in) :: qp_mode !! `options%qp_solver_mode` used
    integer,           intent(in) :: n       !! number of variables
    type(sqpopt_type), intent(in) :: solver  !! the solver, after the solve
    real(wp),          intent(in) :: f_star  !! the known optimal objective
    logical,           intent(in) :: exact   !! whether the exact Hessian was used
    type(sqpopt_results_type) :: r
    type(sqpopt_qp_solver_type) :: qp
    character(len=:), allocatable :: label
    call solver%get_results(r)
    select case (qp_mode)
    case (sqpopt_qp_auto);            label = name//' (auto QP)  '
    case (sqpopt_qp_reduced_hessian); label = name//' (sparse QP)'
    case default;                     label = name//' (dense QP) '
    end select
    if (exact) label = name//' (exact Hessian)'
    print '(A32,A,I6,A,I3,A,ES16.8,A,ES9.2,A,I5,A,I5,A,F7.3,A)', label, ': n=', n, ' istat=', r%istat, &
        ' f=', r%f, ' viol=', r%feasibility_error, ' fc=', r%n_eval_fc, ' gjac=', r%n_eval_gjac, &
        ' time=', r%time, ' s'
    if (exact .and. r%n_eval_hess == 0) error stop 'test_large_sparse FAILED: '//label//': the Hessian was not used'
    ! (the default mode picks the sparse QP for problems this size)
    if (qp_mode == sqpopt_qp_auto .and. n <= qp%auto_dense_max_n) then
        error stop 'test_large_sparse FAILED: problem too small for the sparse QP'
    end if
    if (r%istat /= sqpopt_success .and. r%istat /= sqpopt_acceptable .and. r%istat /= sqpopt_stalled) then
        error stop 'test_large_sparse FAILED: '//label//': did not converge'
    end if
    if (abs(r%f - f_star) > rel_tol*max(1.0_wp, abs(f_star))) then
        error stop 'test_large_sparse FAILED: '//label//': wrong objective'
    end if
    if (r%feasibility_error > feas_tol) then
        error stop 'test_large_sparse FAILED: '//label//': infeasible'
    end if
    end subroutine check

    !------------------------------------------------------------------------
    ! control problem
    !------------------------------------------------------------------------

    subroutine run_control(nn, qp_mode, f_star, exact)
    !! solve the discretized optimal-control problem, and check the solution
    integer,  intent(in)           :: nn      !! number of time steps `N` (`n = 2N+1` variables)
    integer,  intent(in)           :: qp_mode !! `options%qp_solver_mode`
    real(wp), intent(in)           :: f_star  !! the known optimal objective
    logical,  intent(in), optional :: exact   !! use the exact Hessian (default `.false.`)
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
    call problem%set_functions(fc=fc_control, gjac=gjac_control, hess=hess_control)
    ! the Hessian of the Lagrangian is diagonal:
    call problem%set_hessian_sparsity(n, [(k, k=1,n)], [(k, k=1,n)])

    options%max_iter       = 2000
    options%qp_solver_mode = qp_mode
    if (present(exact)) then
        if (exact) options%hessian_mode = sqpopt_hessian_exact
    end if
    allocate(x(n))
    x = 0.0_wp
    x(1:nsteps+1) = 1.0_wp

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x, istat)
    call check('control', qp_mode, n, solver, f_star, options%hessian_mode == sqpopt_hessian_exact)

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

    subroutine gjac_control(x, g, jac, status, data)
    !! objective gradient and constraint Jacobian values of the control problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
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

    subroutine run_rosenbrock(n, qp_mode, f_star, exact)
    !! solve the constrained chained-Rosenbrock problem, and check the solution
    integer,  intent(in)           :: n       !! number of variables
    integer,  intent(in)           :: qp_mode !! `options%qp_solver_mode`
    real(wp), intent(in)           :: f_star  !! the known optimal objective
    logical,  intent(in), optional :: exact   !! use the exact Hessian (default `.false.`)
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
    call problem%set_functions(fc=fc_rosen, gjac=gjac_rosen, hess=hess_rosen)
    ! the Hessian of the Lagrangian is tridiagonal: the diagonal, then the subdiagonal
    call problem%set_hessian_sparsity(2*n-1, [(i, i=1,n), (i+1, i=1,n-1)], [(i, i=1,n), (i, i=1,n-1)])

    options%max_iter       = 2000
    options%qp_solver_mode = qp_mode
    if (present(exact)) then
        if (exact) options%hessian_mode = sqpopt_hessian_exact
    end if
    allocate(x(n))
    do i = 1, n
        x(i) = merge(-1.2_wp, 1.0_wp, mod(i,2) == 1)
    end do

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x, istat)
    call check('rosenbrock', qp_mode, n, solver, f_star, options%hessian_mode == sqpopt_hessian_exact)

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

    subroutine gjac_rosen(x, g, jac, status, data)
    !! objective gradient and constraint Jacobian values of the chained Rosenbrock problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac    !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
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

end program test_large_sparse
