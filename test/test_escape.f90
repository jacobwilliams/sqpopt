program test_escape

    !! Unit test of [[escape_step]]: a point *next to* a symmetry plane of the
    !! problem, where the constraint violation is stationary to within the
    !! convergence tolerance, must still be probed across the plane.
    !!
    !!   subject to  c = x2**2 - x1**2 - 1 >= 0
    !!
    !! `c` is even in `x2`, and on the plane `x2 = 0` the problem is
    !! infeasible, with the violation `1 + x1**2` smallest at `x1 = 0`. At
    !! `x = (1e-7, 1e-9)` the Jacobian is `(-2e-7, 2e-9)`: both columns are
    !! below `ktol = 1e-6` (so the point is stationary for the violation, as
    !! [[check_convergence]] measures it), but the `x2` column is not
    !! negligible relative to the `x1` column. The escape step must find the
    !! decrease of the violation along `x2` (no step along `x1` has one).
    !! Iterates like this arise in `real128` (TP88 of the HS suite), where
    !! round-off doesn't push them off the plane as it does in `real64`.
    !!
    !! With `tol = 0`, only the relative test is left, which finds no
    !! negligible column here, so nothing is probed.

    use sqpopt_restoration_module, only: escape_step
    use sqpopt_problem_module,     only: sqpopt_problem_type
    use sqpopt_types_module,       only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_line_search_failed
    use sqpopt_kinds,              only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: ktol = 1.0e-6_wp !! the convergence test's default tolerance
    real(wp), parameter :: x0(2) = [1.0e-7_wp, 1.0e-9_wp] !! the point next to the plane `x2 = 0`

    type(sqpopt_problem_type)  :: problem
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: x_new(2), c(1), c_new(1)
    integer  :: istat

    write(*,*) '----------------------------'
    write(*,*) 'test_escape'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-10.0_wp,-10.0_wp], x_ub=[10.0_wp,10.0_wp], c_lb=[0.0_wp], c_ub=[1.0e20_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc, gjac=gjac)
    call problem%reset_evaluations()

    jac%nrows = 1
    jac%ncols = 2
    jac%nnz  = 2
    jac%irow = [1, 1]
    jac%icol = [1, 2]
    jac%val  = [-2.0_wp*x0(1), 2.0_wp*x0(2)]
    call problem%c(x0, c)

    call escape_step(problem, jac, x0, c, ktol, x_new, istat)
    call problem%c(x_new, c_new)
    print '(A,2ES12.4,A,ES12.4,A,ES12.4,A,I0)', 'tol=ktol: x_new=', x_new, '  c: ', c(1), ' -> ', c_new(1), '  istat=', istat
    if (istat /= sqpopt_success) error stop 'test_escape FAILED: no escape step found'
    if (x_new(1) /= x0(1) .or. abs(x_new(2)) < 1.0e-4_wp) error stop 'test_escape FAILED: the step is not along x2'
    if (c_new(1) <= c(1)) error stop 'test_escape FAILED: the violation did not decrease'

    call escape_step(problem, jac, x0, c, 0.0_wp, x_new, istat)
    print '(A,2ES12.4,A,I0)', 'tol=0:    x_new=', x_new, '  istat=', istat
    if (istat /= sqpopt_line_search_failed .or. any(x_new /= x0)) &
        error stop 'test_escape FAILED: a column above the relative tolerance was probed'

    print '(A)', 'test_escape PASSED'

    contains

    subroutine fc(x, f, c, status, data)
    !! the objective (unused here) and the constraint
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    associate(unused => status); end associate
    if (present(data)) continue
    f = x(1)**2 + x(2)**2
    c(1) = x(2)**2 - x(1)**2 - 1.0_wp
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient (unused here) and the constraint Jacobian
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the constraint Jacobian at `x`
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    associate(unused => accuracy, unused2 => status); end associate
    if (present(data)) continue
    g = 2.0_wp*x
    jac_val = [-2.0_wp*x(1), 2.0_wp*x(2)]
    end subroutine gjac

end program test_escape
