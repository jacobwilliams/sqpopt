program test_input_validation

    !! Regression test: invalid problem definitions and options are
    !! rejected up front with `sqpopt_invalid_input` (and a message saying
    !! what is wrong), rather than crashing or silently misbehaving. Also
    !! checks that a starting point outside the variable bounds is moved
    !! inside them before any function is evaluated.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type, sqpopt_hessian_exact
    use sqpopt_types_module,   only: sqpopt_invalid_input, sqpopt_success
    use sqpopt_linesearch_module,   only: sqpopt_linesearch_type
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_qp_solver_module,    only: sqpopt_qp_solver_type
    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf

    implicit none

    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    logical :: outside_bounds  !! set if a function is ever evaluated outside the variable bounds

    write(*,*) '----------------------------'
    write(*,*) 'test_input_validation'
    write(*,*) '----------------------------'

    ! problem size never set:
    problem = sqpopt_problem_type()
    call expect_invalid('problem size not set', problem, sqpopt_options_type(), [0.0_wp, 0.0_wp])

    ! functions never set:
    problem = sqpopt_problem_type()
    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds([-1.0_wp,-1.0_wp], [1.0_wp,1.0_wp], [1.0_wp], [1.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call expect_invalid('functions not set', problem, sqpopt_options_type(), [0.0_wp, 0.0_wp])

    ! lower bound above upper bound:
    call valid_problem(problem)
    call problem%set_bounds([-1.0_wp, 2.0_wp], [1.0_wp,1.0_wp], [1.0_wp], [1.0_wp])
    call expect_invalid('x_lb > x_ub', problem, sqpopt_options_type(), [0.0_wp, 0.0_wp])

    ! Jacobian index out of range:
    call valid_problem(problem)
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,2], icol=[1,2])
    call expect_invalid('Jacobian row index > m', problem, sqpopt_options_type(), [0.0_wp, 0.0_wp])

    ! wrong-size starting point:
    call valid_problem(problem)
    call expect_invalid('size(x0) /= n', problem, sqpopt_options_type(), [0.0_wp, 0.0_wp, 0.0_wp])

    ! non-finite starting point (projecting a NaN onto the bounds would silently
    ! turn it into a bound) and multipliers:
    call valid_problem(problem)
    call expect_invalid('NaN in x0', problem, sqpopt_options_type(), [ieee_value(1.0_wp, ieee_quiet_nan), 0.0_wp])
    call expect_invalid('Inf in x0', problem, sqpopt_options_type(), [0.0_wp, ieee_value(1.0_wp, ieee_positive_inf)])
    block
        type(sqpopt_type) :: solver
        integer :: istat
        call solver%initialize(problem=problem)
        call solver%solve([0.0_wp, 0.0_wp], istat, lambda0=[ieee_value(1.0_wp, ieee_quiet_nan)])
        print '(A,I0,2A)', 'NaN in lambda0: istat=', istat, '  ', solver%status_message()
        if (istat /= sqpopt_invalid_input) error stop 'test_input_validation FAILED: NaN in lambda0'
    end block

    ! invalid options:
    call valid_problem(problem)
    options = sqpopt_options_type()
    options%lbfgs_memory = -1
    call expect_invalid('lbfgs_memory = -1', problem, options, [0.0_wp, 0.0_wp])
    options = sqpopt_options_type()
    options%qp_solver_mode = 99
    call expect_invalid('qp_solver_mode = 99', problem, options, [0.0_wp, 0.0_wp])
    options = sqpopt_options_type()
    options%ktol = 0.0_wp
    call expect_invalid('ktol = 0', problem, options, [0.0_wp, 0.0_wp])
    options = sqpopt_options_type()
    options%acceptable_iter = -1
    call expect_invalid('acceptable_iter = -1', problem, options, [0.0_wp, 0.0_wp])
    options = sqpopt_options_type()
    options%print_level = 1
    options%output_unit = 98765   ! (not an open unit)
    call expect_invalid('output_unit not open', problem, options, [0.0_wp, 0.0_wp])
    options = sqpopt_options_type()
    options%hessian_mode = sqpopt_hessian_exact   ! (the problem has no hess function)
    call expect_invalid('exact Hessian without hess', problem, options, [0.0_wp, 0.0_wp])
    block
        type(sqpopt_problem_type) :: p2
        call valid_problem(p2)
        call p2%set_hessian_sparsity(nnz=1, irow=[3], icol=[1])
        call expect_invalid('Hessian row index > n', p2, sqpopt_options_type(), [0.0_wp, 0.0_wp])
    end block

    ! invalid component settings:
    block
        type(sqpopt_linesearch_type)   :: ls
        type(sqpopt_trust_region_type) :: tr
        type(sqpopt_qp_solver_type)    :: qp
        ls%sigma = 1.5_wp
        call expect_invalid('linesearch%sigma = 1.5', problem, sqpopt_options_type(), [0.0_wp, 0.0_wp], linesearch=ls)
        tr%enabled = .true.
        tr%eta1 = 0.9_wp
        tr%eta2 = 0.5_wp
        call expect_invalid('trust_region eta1 > eta2', problem, sqpopt_options_type(), [0.0_wp, 0.0_wp], trust_region=tr)
        qp%max_step = 0.0_wp
        call expect_invalid('qp_solver%max_step = 0', problem, sqpopt_options_type(), [0.0_wp, 0.0_wp], qp_solver=qp)
    end block

    ! wrong-size initial multipliers:
    block
        type(sqpopt_type) :: solver
        integer :: istat
        call solver%initialize(problem=problem)
        call solver%solve([0.0_wp, 0.0_wp], istat, lambda0=[1.0_wp, 2.0_wp])
        print '(A,I0,2A)', 'size(lambda0) /= m: istat=', istat, '  ', solver%status_message()
        if (istat /= sqpopt_invalid_input) error stop 'test_input_validation FAILED: size(lambda0) /= m'
    end block

    ! a valid problem, started outside the bounds [-2,2]: x0 is projected
    ! onto the bounds before anything is evaluated:
    block
        type(sqpopt_type) :: solver
        real(wp) :: xsol(2), lam(1)
        integer :: istat
        call valid_problem(problem)
        outside_bounds = .false.
        call solver%initialize(problem=problem)
        call solver%solve([5.0_wp, -0.5_wp], istat)
        call solver%get_solution(xsol, lam)
        print '(A,2F10.6,A,I0)', 'x0 outside the bounds: x=', xsol, '  istat=', istat
        if (istat /= sqpopt_success) error stop 'test_input_validation FAILED: valid problem did not converge'
        if (maxval(abs(xsol + 1.0_wp/sqrt(2.0_wp))) > 1.0e-6_wp) then
            error stop 'test_input_validation FAILED: valid problem converged to the wrong point'
        end if
        if (outside_bounds) error stop 'test_input_validation FAILED: a function was evaluated outside the bounds'
    end block

    print '(A)', 'test_input_validation PASSED'

    contains

    subroutine expect_invalid(label, problem, options, x0, linesearch, trust_region, qp_solver)
    !! check that `solve` rejects the inputs with `sqpopt_invalid_input`
    character(len=*),          intent(in)                :: label        !! the case, for the output
    type(sqpopt_problem_type), intent(in)                :: problem      !! problem definition
    type(sqpopt_options_type), intent(in)                :: options      !! solver options
    real(wp), dimension(:),    intent(in)                :: x0           !! starting point
    type(sqpopt_linesearch_type),   intent(in), optional :: linesearch   !! line search settings (default if absent)
    type(sqpopt_trust_region_type), intent(in), optional :: trust_region !! trust-region settings (default if absent)
    type(sqpopt_qp_solver_type),    intent(in), optional :: qp_solver    !! QP solver settings (default if absent)
    type(sqpopt_type) :: solver
    integer :: istat
    call solver%initialize(problem=problem, options=options, linesearch=linesearch, trust_region=trust_region, &
                           qp_solver=qp_solver)
    call solver%solve(x0, istat)
    print '(A,A,I0,2A)', label, ': istat=', istat, '  ', solver%status_message()
    if (istat /= sqpopt_invalid_input) error stop 'test_input_validation FAILED: '//label
    end subroutine expect_invalid

    !> minimize x1 + x2 s.t. x1^2 + x2^2 = 1, -2 <= x <= 2
    !! (solution x = -(1,1)/sqrt(2))
    subroutine valid_problem(problem)
    type(sqpopt_problem_type), intent(out) :: problem !! the (valid) problem definition
    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds([-2.0_wp,-2.0_wp], [2.0_wp,2.0_wp], [1.0_wp], [1.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv)
    end subroutine valid_problem

    subroutine note_bounds(x)
    !! record whether a user function was called outside the variable bounds
    real(wp), dimension(:), intent(in) :: x !! the point a user function was called at
    if (any(abs(x) > 2.0_wp)) outside_bounds = .true.
    end subroutine note_bounds

    subroutine obj(x, f, status, data)
    !! the objective
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    call note_bounds(x)
    f = x(1) + x(2)
    end subroutine obj

    subroutine grad(x, g, status, data)
    !! the objective's gradient
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    call note_bounds(x)
    g = 1.0_wp
    end subroutine grad

    subroutine cons(x, c, status, data)
    !! the constraints
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    call note_bounds(x)
    c(1) = x(1)**2 + x(2)**2
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    !! the nonzero values of the constraint Jacobian (in the sparsity pattern's order)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    call note_bounds(x)
    jac_val = 2.0_wp*x
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


end program test_input_validation
