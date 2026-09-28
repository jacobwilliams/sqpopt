program test_derivative_accuracy

    !! Test of `options%derivative_accuracy` (see [[sqpopt_options_module]]):
    !! with `sqpopt_derivatives_fast`, the solver asks `gjac` for fast
    !! derivatives first and switches, once and for good, to accurate ones
    !! before it stops. Here the "fast" derivatives are forward differences
    !! with a deliberately large step (so they are far off), and the
    !! solution must still match the one found with accurate (analytic)
    !! derivatives throughout. With the default options, `gjac` is only ever
    !! asked for accurate derivatives.
    !!
    !!   minimize   (1-x1)^2 + 100*(x2-x1^2)^2    (Rosenbrock)
    !!   subject to x1^2 + x2^2 <= 1.5

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type, sqpopt_derivatives_fast, sqpopt_derivatives_accurate
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_types_module,   only: sqpopt_success, sqpopt_results_type
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    real(wp) :: xa(2), xf(2), lam(1)
    integer  :: istat
    integer  :: n_fast, n_accurate !! `gjac` calls asking for fast/accurate derivatives
    logical  :: back_to_fast       !! whether a fast call came after an accurate one

    write(*,*) '----------------------------'
    write(*,*) 'test_derivative_accuracy'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-5.0_wp,-5.0_wp], x_ub=[5.0_wp,5.0_wp], c_lb=[-1.0e20_wp], c_ub=[1.5_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_fd_or_exact)

    ! ---- the default: accurate derivatives throughout ----
    call reset_counts()
    call solver%initialize(problem=problem, options=options)
    call solver%solve([-1.2_wp, 1.0_wp], istat)
    call solver%get_solution(xa, lam)
    call solver%get_results(r)
    print '(A,2F12.8,A,I0,A,2I4)', 'accurate: x=', xa, ' istat=', istat, ' gjac calls fast/accurate:', n_fast, n_accurate
    if (istat /= sqpopt_success) error stop 'test_derivative_accuracy FAILED: accurate solve did not converge'
    if (n_fast /= 0) error stop 'test_derivative_accuracy FAILED: fast derivatives asked for by default'
    if (r%derivative_switch_iteration /= 0) error stop 'test_derivative_accuracy FAILED: switch reported by default'

    ! ---- fast, then accurate ----
    call reset_counts()
    options%derivative_accuracy = sqpopt_derivatives_fast
    call solver%initialize(problem=problem, options=options)
    call solver%solve([-1.2_wp, 1.0_wp], istat)
    call solver%get_solution(xf, lam)
    call solver%get_results(r)
    print '(A,2F12.8,A,I0,A,2I4,A,I0)', 'fast:     x=', xf, ' istat=', istat, ' gjac calls fast/accurate:', &
        n_fast, n_accurate, ' switched at iteration ', r%derivative_switch_iteration
    if (istat /= sqpopt_success) error stop 'test_derivative_accuracy FAILED: fast solve did not converge'
    if (n_fast == 0 .or. n_accurate == 0) error stop 'test_derivative_accuracy FAILED: no switch'
    if (back_to_fast) error stop 'test_derivative_accuracy FAILED: switched back to fast derivatives'
    if (r%derivative_switch_iteration <= 0) error stop 'test_derivative_accuracy FAILED: switch not reported'
    if (maxval(abs(xf - xa)) > 1.0e-6_wp) error stop 'test_derivative_accuracy FAILED: inaccurate solution'

    ! ---- an invalid value ----
    options%derivative_accuracy = 0
    call solver%initialize(problem=problem, options=options)
    call solver%solve([-1.2_wp, 1.0_wp], istat)
    if (istat == sqpopt_success) error stop 'test_derivative_accuracy FAILED: invalid derivative_accuracy accepted'

    print '(A)', 'test_derivative_accuracy PASSED'

    contains

    subroutine reset_counts()
    !! zero the `gjac` call counts
    n_fast = 0
    n_accurate = 0
    back_to_fast = .false.
    end subroutine reset_counts

    subroutine fc_obj_cons(x, f, c, status, data)
    !! `fc` for `set_functions`: the objective and the constraint
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = (1.0_wp - x(1))**2 + 100.0_wp*(x(2) - x(1)**2)**2
    c(1) = x(1)**2 + x(2)**2
    associate(unused => status); end associate
    end subroutine fc_obj_cons

    subroutine gjac_fd_or_exact(x, g, jac_val, accuracy, status, data)
    !! `gjac` for `set_functions`: coarse forward differences of `fc` when
    !! fast derivatives are asked for, else the analytic derivatives
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    real(wp), parameter :: h = 1.0e-4_wp !! (a deliberately large step, so the fast derivatives are poor)
    real(wp) :: f0, f1, c0(1), c1(1), xp(2)
    integer :: j
    select case (accuracy)
    case (sqpopt_derivatives_fast)
        n_fast = n_fast + 1
        if (n_accurate > 0) back_to_fast = .true.
        call fc_obj_cons(x, f0, c0, status, data)
        do j = 1, 2
            xp = x
            xp(j) = x(j) + h
            call fc_obj_cons(xp, f1, c1, status, data)
            g(j) = (f1 - f0)/h
            jac_val(j) = (c1(1) - c0(1))/h
        end do
    case (sqpopt_derivatives_accurate)
        n_accurate = n_accurate + 1
        g = [-2.0_wp*(1.0_wp - x(1)) - 400.0_wp*x(1)*(x(2) - x(1)**2), 200.0_wp*(x(2) - x(1)**2)]
        jac_val = [2.0_wp*x(1), 2.0_wp*x(2)]
    case default
        error stop 'test_derivative_accuracy FAILED: invalid accuracy argument'
    end select
    end subroutine gjac_fd_or_exact

end program test_derivative_accuracy
