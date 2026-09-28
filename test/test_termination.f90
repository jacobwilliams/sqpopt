program test_termination

    !! Tests of the optional termination criteria, and of writing the
    !! iteration log to a given unit:
    !!
    !! * `max_evals`: stop after that many objective evaluations;
    !! * `max_time`: stop once the time limit is reached;
    !! * `obj_lower_limit`: stop (`sqpopt_unbounded`) if the objective falls
    !!   below it at a feasible point (on an unbounded problem);
    !! * `acceptable_*`: with an unreachable `ktol`, stop at the acceptable
    !!   level (here, the first iterate with a KKT error below 1e-3) instead
    !!   of continuing;
    !! * `print_level=2` with `output_unit` set to a scratch file: the log
    !!   and summary are written there;
    !! * `print_level=3`: the detailed log has each of its parts (the method,
    !!   the column headings, detail lines, the events, and the solution
    !!   tables), the results' event counts and times are consistent, and the
    !!   solution is the same as with `print_level=0` (printing must not
    !!   change the result).
    !!
    !! Problem (bounded case): the Rosenbrock function subject to
    !! x1^2 + x2^2 <= 1.5 (solution near (0.9072, 0.8228)).

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_types_module,   only: sqpopt_max_evals_reached, sqpopt_time_limit_reached, sqpopt_unbounded, &
                                     sqpopt_acceptable, sqpopt_results_type
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem, unbounded
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat, u, ios
    character(len=200) :: line
    logical :: found_summary

    write(*,*) '----------------------------'
    write(*,*) 'test_termination'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-5.0_wp,-5.0_wp], x_ub=[5.0_wp,5.0_wp], c_lb=[-1.0e20_wp], c_ub=[1.5_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv)

    ! ---- max_evals ----
    options = sqpopt_options_type()
    options%max_evals = 5
    call solver%initialize(problem=problem, options=options)
    call solver%solve([-1.2_wp, 1.0_wp], istat)
    call solver%get_results(r)
    print '(A,I0,A,I0)', 'max_evals: istat=', istat, ' n_eval_fc=', r%n_eval_fc
    if (istat /= sqpopt_max_evals_reached .or. r%n_eval_fc < 5) error stop 'test_termination FAILED: max_evals'

    ! ---- max_time ----
    options = sqpopt_options_type()
    options%max_time = 1.0e-9_wp
    call solver%initialize(problem=problem, options=options)
    call solver%solve([-1.2_wp, 1.0_wp], istat)
    call solver%get_results(r)
    print '(A,I0,A,I0)', 'max_time: istat=', istat, ' iterations=', r%iterations
    if (istat /= sqpopt_time_limit_reached .or. r%iterations /= 1) error stop 'test_termination FAILED: max_time'

    ! ---- obj_lower_limit (minimize -x1-x2 with x unbounded) ----
    call unbounded%set_problem_size(n=2, m=0)
    call unbounded%set_bounds(x_lb=[-1.0e20_wp,-1.0e20_wp], x_ub=[1.0e20_wp,1.0e20_wp], &
                              c_lb=[real(wp)::], c_ub=[real(wp)::])
    call unbounded%set_jacobian_sparsity(nnz=0, irow=[integer::], icol=[integer::])
    call unbounded%set_functions(fc=fc_obj_lin_cons0, gjac=gjac_grad_lin_jac0)
    options = sqpopt_options_type()
    options%obj_lower_limit = -10.0_wp
    call solver%initialize(problem=unbounded, options=options)
    call solver%solve([0.0_wp, 0.0_wp], istat)
    call solver%get_results(r)
    print '(A,I0,A,F10.4,2A)', 'unbounded: istat=', istat, ' f=', r%f, '  ', solver%status_message()
    if (istat /= sqpopt_unbounded .or. r%f >= -10.0_wp) error stop 'test_termination FAILED: obj_lower_limit'

    ! ---- acceptable level, with the log written to a scratch file ----
    open(newunit=u, status='scratch', action='readwrite', form='formatted')
    options = sqpopt_options_type()
    options%ktol = 1.0e-300_wp            ! (effectively unreachable)
    options%acceptable_ktol = 1.0e-3_wp
    options%acceptable_iter = 1
    options%print_level = 2
    options%output_unit = u
    call solver%initialize(problem=problem, options=options)
    call solver%solve([-1.2_wp, 1.0_wp], istat)
    call solver%get_results(r)
    print '(A,I0,A,2F10.6,A,I0)', 'acceptable: istat=', istat, ' x=', r%x, ' iterations=', r%iterations
    if (istat /= sqpopt_acceptable) error stop 'test_termination FAILED: acceptable'
    if (r%kkt_error > 1.0e-3_wp) error stop 'test_termination FAILED: acceptable KKT error'
    if (maxval(abs(r%x - [0.9072340_wp, 0.8227555_wp])) > 1.0e-2_wp) error stop 'test_termination FAILED: acceptable x'
    rewind(u)
    found_summary = .false.
    do
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line, 'sqpopt: status') > 0) found_summary = .true.
    end do
    close(u)
    if (.not. found_summary) error stop 'test_termination FAILED: no summary written to output_unit'

    ! ---- the detailed log (print_level=3), which must not change the result ----
    block
        type(sqpopt_results_type) :: r0
        logical, dimension(6) :: found
        character(len=*), dimension(6), parameter :: parts = [character(len=24) :: '   method:', '  iter ', &
            '        . ', '   events ', '   variables:', '   constraints:']
        integer :: k
        options = sqpopt_options_type()
        call solver%initialize(problem=problem, options=options)
        call solver%solve([-1.2_wp, 1.0_wp], istat)
        call solver%get_results(r0)
        open(newunit=u, status='scratch', action='readwrite', form='formatted')
        options%print_level = 3
        options%output_unit = u
        call solver%initialize(problem=problem, options=options)
        call solver%solve([-1.2_wp, 1.0_wp], istat)
        call solver%get_results(r)
        rewind(u)
        found = .false.
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            do k = 1, size(parts)
                if (index(line, trim(parts(k))) == 1) found(k) = .true.
            end do
        end do
        close(u)
        print '(A,6L2,A,I0,A,I0)', 'print_level=3: parts found', found, '  QP iterations ', r%n_qp_iterations, &
            '  iterations ', r%iterations
        if (.not. all(found)) error stop 'test_termination FAILED: a part of the detailed log is missing'
        if (any(r%x /= r0%x) .or. r%iterations /= r0%iterations .or. r%n_eval_fc /= r0%n_eval_fc) &
            error stop 'test_termination FAILED: printing changed the result'
        if (r%n_qp_iterations < r%iterations - 1 .or. r%time_functions < 0.0_wp .or. r%time_qp < 0.0_wp .or. &
            r%time_functions + r%time_qp > r%time + 1.0e-3_wp) error stop 'test_termination FAILED: counts or times'
    end block

    print '(A)', 'test_termination PASSED'

    contains

    subroutine obj(x, f, status, data)
    !! the objective
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = (1.0_wp-x(1))**2 + 100.0_wp*(x(2)-x(1)**2)**2
    end subroutine obj

    subroutine grad(x, g, status, data)
    !! the objective's gradient
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    g(1) = -2.0_wp*(1.0_wp-x(1)) - 400.0_wp*x(1)*(x(2)-x(1)**2)
    g(2) = 200.0_wp*(x(2)-x(1)**2)
    end subroutine grad

    subroutine cons(x, c, status, data)
    !! the constraints
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    c(1) = x(1)**2 + x(2)**2
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    !! the nonzero values of the constraint Jacobian (in the sparsity pattern's order)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    jac_val = 2.0_wp*x
    end subroutine jacv

    subroutine obj_lin(x, f, status, data)
    !! the objective of the unbounded case, \( -x_1-x_2 \)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = -x(1) - x(2)
    end subroutine obj_lin

    subroutine grad_lin(x, g, status, data)
    !! the gradient of `obj_lin`
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    g = -1.0_wp
    end subroutine grad_lin

    subroutine cons0(x, c, status, data)
    !! the constraints of the unbounded case (it has none)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    end subroutine cons0

    subroutine jac0(x, jac_val, status, data)
    !! the Jacobian of the unbounded case (it has no constraints)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    end subroutine jac0

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

    subroutine fc_obj_lin_cons0(x, f, c, status, data)
    !! `fc` for `set_functions`: the objective (`obj_lin`) and the constraints (`cons0`)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    call obj_lin(x, f, status, data)
    if (status == 0) call cons0(x, c, status, data)
    end subroutine fc_obj_lin_cons0

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

    subroutine gjac_grad_lin_jac0(x, g, jac_val, status, data)
    !! `gjac` for `set_functions`: the gradient (`grad_lin`) and the Jacobian values (`jac0`)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g       !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    call grad_lin(x, g, status, data)
    if (status == 0) call jac0(x, jac_val, status, data)
    end subroutine gjac_grad_lin_jac0


end program test_termination
