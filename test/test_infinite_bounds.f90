program test_infinite_bounds

    !! Regression test: `+/-huge(1.0_wp)` may be used for absent bounds
    !! (they are treated as infinite, see `sqpopt_infinity`) without any
    !! floating-point overflow. Floating-point halting is enabled for
    !! overflow, invalid operations, and division by zero, so any such
    !! exception anywhere in the solver stops the test.
    !!
    !!   minimize   (x1-2)^2 + (x2-3)^2
    !!   subject to x1 + x2 <= 3   (no lower bound)
    !!              x unbounded
    !!
    !! known solution: x* = (1, 2)

    use, intrinsic :: ieee_exceptions, only: ieee_set_halting_mode, ieee_overflow, ieee_invalid, ieee_divide_by_zero
    use sqpopt_module,            only: sqpopt_type
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_auto, sqpopt_qp_dense, sqpopt_qp_reduced_hessian
    use sqpopt_linesearch_module, only: sqpopt_linesearch_armijo, sqpopt_linesearch_filter
    use sqpopt_types_module,      only: sqpopt_success
    use sqpopt_kinds,             only: wp => sqpopt_module_wp

    implicit none

    integer,  parameter :: modes(3) = [sqpopt_qp_auto, sqpopt_qp_dense, sqpopt_qp_reduced_hessian]
    real(wp), parameter :: inf = huge(1.0_wp)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: xsol(2), lam(1)
    integer  :: istat, i, ls

    write(*,*) '----------------------------'
    write(*,*) 'test_infinite_bounds'
    write(*,*) '----------------------------'

    call ieee_set_halting_mode(ieee_overflow,       .true.)
    call ieee_set_halting_mode(ieee_invalid,        .true.)
    call ieee_set_halting_mode(ieee_divide_by_zero, .true.)

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-inf,-inf], x_ub=[inf,inf], c_lb=[-inf], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv)

    do ls = 1, 2
        do i = 1, size(modes)
            options = sqpopt_options_type()
            options%qp_solver_mode  = modes(i)
            options%linesearch_mode = merge(sqpopt_linesearch_armijo, sqpopt_linesearch_filter, ls == 1)
            call solver%initialize(problem=problem, options=options)
            call solver%solve([0.0_wp, 0.0_wp], istat)
            call solver%get_solution(xsol, lam)
            print '(A,I0,A,I0,A,2F10.6,A,I0)', 'linesearch_mode=', options%linesearch_mode, &
                ' qp_solver_mode=', modes(i), ': x=', xsol, '  istat=', istat
            if (istat /= sqpopt_success) error stop 'test_infinite_bounds FAILED: did not reach sqpopt_success'
            if (maxval(abs(xsol-xexpect)) > 1.0e-5_wp) error stop 'test_infinite_bounds FAILED: wrong solution'
        end do
    end do

    print '(A)', 'test_infinite_bounds PASSED'

    contains

    subroutine obj(x, f, status, data)
    !! the objective
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = (x(1)-2.0_wp)**2 + (x(2)-3.0_wp)**2
    end subroutine obj

    subroutine grad(x, g, status, data)
    !! the objective's gradient
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    g = [2.0_wp*(x(1)-2.0_wp), 2.0_wp*(x(2)-3.0_wp)]
    end subroutine grad

    subroutine cons(x, c, status, data)
    !! the constraints
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    c(1) = x(1) + x(2)
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    !! the nonzero values of the constraint Jacobian (in the sparsity pattern's order)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    associate(unused => x); end associate
    jac_val = 1.0_wp
    end subroutine jacv

    subroutine fc_obj_cons(x, f, c, status, data)
    !! `fc` for `set_functions`: the objective (`obj`) and the constraints (`cons`)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    call obj(x, f, status, data)
    if (status == 0) call cons(x, c, status, data)
    end subroutine fc_obj_cons

    subroutine gjac_grad_jacv(x, g, jac_val, status, data)
    !! `gjac` for `set_functions`: the gradient (`grad`) and the Jacobian values (`jacv`)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g       !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    call grad(x, g, status, data)
    if (status == 0) call jacv(x, jac_val, status, data)
    end subroutine gjac_grad_jacv


end program test_infinite_bounds
