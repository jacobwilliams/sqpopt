program test_sr1_convexify

    !! The convexified limited-memory SR1 Hessian (`hessian%convexify`, see
    !! [[sqpopt_spectral_module]] and [[sqpopt_iterate_module]]): before each
    !! QP, an indefinite SR1 matrix is shifted by its smallest eigenvalue, so
    !! that every QP is convex. The problem is a sum of double wells,
    !!
    !!     minimize   sum_i ((x_i^2 - 1)^2 + 0.1*x_i) + 0.05*sum_{i<n} (x_{i+1} - x_i)^2
    !!     subject to -2 <= x <= 2
    !!
    !! from a point near the humps (`x_i = 0.1`), where the curvature is
    !! negative, so the SR1 matrix becomes indefinite. With 20 variables (the
    !! dense QP solver) and 300 (the sparse one), it checks that
    !!
    !! * the solve converges to a minimum (every `|x_i|` near 1), and the
    !!   detailed log has the "convexified" detail line and the method line's
    !!   "convexified";
    !! * the detailed log (`print_level=3`) doesn't change the results;
    !! * a second `solve` gives the same results (no state carries over);
    !! * the eigensolvers (`hessian%eigen_solver`: QL, Jacobi, and LAPACK if
    !!   the library has it) give the same solution.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_hessian_module, only: sqpopt_hessian_sr1, sqpopt_hessian_type
    use sqpopt_eigen_module,   only: sqpopt_eigen_auto, sqpopt_eigen_ql, sqpopt_eigen_jacobi, sqpopt_eigen_lapack, &
                                     sqpopt_eigen_has_lapack
    use sqpopt_types_module,   only: sqpopt_results_type, sqpopt_success, sqpopt_acceptable
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: sizes(2) = [20, 300]
    integer :: k

    write(*,*) '----------------------------'
    write(*,*) 'test_sr1_convexify'
    write(*,*) '----------------------------'

    do k = 1, size(sizes)
        call test(sizes(k))
    end do

    print '(A)', 'test_sr1_convexify PASSED'

    contains

    subroutine test(n)
    !! all the checks (see the program's documentation) with `n` variables
    integer, intent(in) :: n !! number of variables

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_hessian_type) :: hessian
    type(sqpopt_results_type) :: r0, r1, r2, re
    real(wp), dimension(:), allocatable :: x0, x_lb, x_ub
    integer, dimension(0) :: no_rows
    integer, dimension(3) :: solvers
    character(len=256) :: line
    logical :: found_detail, found_method
    integer :: u, ios, istat, i

    allocate(x0(n), x_lb(n), x_ub(n))
    x0   = 0.1_wp
    x_lb = -2.0_wp
    x_ub =  2.0_wp
    call problem%set_problem_size(n=n, m=0)
    call problem%set_bounds(x_lb=x_lb, x_ub=x_ub, c_lb=[real(wp) ::], c_ub=[real(wp) ::])
    call problem%set_jacobian_sparsity(nnz=0, irow=no_rows, icol=no_rows)
    call problem%set_functions(fc=fc, gjac=gjac)

    options%hessian_mode = sqpopt_hessian_sr1
    options%max_iter     = 500
    hessian%convexify    = .true.

    ! ---- the detailed log, which must have the convexification's lines ----
    open(newunit=u, status='scratch', action='readwrite', form='formatted')
    options%print_level = 3
    options%output_unit = u
    call solver%initialize(problem=problem, options=options, hessian=hessian)
    call solver%solve(x0, istat)
    call solver%get_results(r0)
    rewind(u)
    found_detail = .false.
    found_method = .false.
    do
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line, 'SR1 Hessian convexified') > 0) found_detail = .true.
        if (index(line, 'method:') > 0 .and. index(line, 'convexified') > 0) found_method = .true.
    end do
    close(u)
    call check_solution('detailed log', istat, r0)
    print '(A,I4,A,I0,A,I0,A,2L2)', 'n =', n, ': iterations ', r0%iterations, ', fc ', r0%n_eval_fc, &
        ', detail line and method found:', found_detail, found_method
    if (.not. found_detail) error stop 'test_sr1_convexify FAILED: the matrix was never convexified'
    if (.not. found_method) error stop 'test_sr1_convexify FAILED: the method line does not say "convexified"'

    ! ---- without printing, twice on the same object: the same results ----
    options%print_level = 0
    call solver%initialize(problem=problem, options=options, hessian=hessian)
    call solver%solve(x0, istat)
    call solver%get_results(r1)
    call check_solution('print_level=0', istat, r1)
    call solver%solve(x0, istat)
    call solver%get_results(r2)
    call check_solution('second solve', istat, r2)
    if (any(r1%x /= r0%x) .or. r1%iterations /= r0%iterations .or. r1%n_eval_fc /= r0%n_eval_fc) &
        error stop 'test_sr1_convexify FAILED: printing changed the result'
    if (any(r2%x /= r1%x) .or. r2%iterations /= r1%iterations .or. r2%n_eval_fc /= r1%n_eval_fc) &
        error stop 'test_sr1_convexify FAILED: the second solve differs'

    ! ---- every eigensolver: the same solution ----
    solvers = [sqpopt_eigen_ql, sqpopt_eigen_jacobi, sqpopt_eigen_lapack]
    do i = 1, size(solvers)
        if (solvers(i) == sqpopt_eigen_lapack .and. .not. sqpopt_eigen_has_lapack) cycle
        hessian%eigen_solver = solvers(i)
        call solver%initialize(problem=problem, options=options, hessian=hessian)
        call solver%solve(x0, istat)
        call solver%get_results(re)
        call check_solution('eigen_solver', istat, re)
        print '(A,I0,A,I0,A,ES9.2)', '  eigen_solver ', solvers(i), ': iterations ', re%iterations, &
            ', |x - x_auto| = ', maxval(abs(re%x - r1%x))
        if (maxval(abs(re%x - r1%x)) > 1.0e-5_wp) error stop 'test_sr1_convexify FAILED: eigensolvers differ'
    end do
    hessian%eigen_solver = sqpopt_eigen_auto

    end subroutine test

    subroutine check_solution(label, istat, r)
    !! check that a solve converged to a minimum: every `|x_i|` near 1
    character(len=*),          intent(in) :: label !! the case, for the output
    integer,                   intent(in) :: istat !! the status returned by `solve`
    type(sqpopt_results_type), intent(in) :: r     !! its results
    if (istat /= sqpopt_success .and. istat /= sqpopt_acceptable) then
        print '(A,A,I0)', label, ': istat = ', istat
        error stop 'test_sr1_convexify FAILED: a solve did not converge'
    end if
    if (r%istat /= istat) error stop 'test_sr1_convexify FAILED: results%istat differs'
    if (any(abs(abs(r%x) - 1.0_wp) > 0.1_wp)) then
        print '(A,A,2ES12.4)', label, ': min, max |x| ', minval(abs(r%x)), maxval(abs(r%x))
        error stop 'test_sr1_convexify FAILED: not a minimum'
    end if
    end subroutine check_solution

    subroutine fc(x, f, c, status, data)
    !! the objective (no constraints)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(0)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (none)
    integer :: n
    n = size(x)
    f = sum((x**2 - 1.0_wp)**2 + 0.1_wp*x) + 0.05_wp*sum((x(2:n) - x(1:n-1))**2)
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient (no constraints)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` `dimension(0)`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    integer :: n
    n = size(x)
    g = 4.0_wp*x*(x**2 - 1.0_wp) + 0.1_wp
    g(2:n)   = g(2:n)   + 0.1_wp*(x(2:n) - x(1:n-1))
    g(1:n-1) = g(1:n-1) - 0.1_wp*(x(2:n) - x(1:n-1))
    end subroutine gjac

end program test_sr1_convexify
