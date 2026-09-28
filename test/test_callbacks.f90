program test_callbacks

    !! Tests of the user-function interface (see [[sqpopt_problem_module]]):
    !!
    !! * user data: a `class(*)` object given to `set_functions` is passed
    !!   to every user function and to the `report` callback, and updates to
    !!   it are seen by the caller;
    !! * `status > 0`: the objective reports that it can't be evaluated for
    !!   `x1 > 2.8`; the solver must back off from those points and still
    !!   converge;
    !! * `status < 0`: a user function asks the solver to stop.
    !!
    !!   minimize   (x1-3)^2 + (x2-1)^2
    !!   subject to x1 + x2 <= 3        (solution x* = (2.5, 0.5))

    use sqpopt_module,            only: sqpopt_type
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_types_module,      only: sqpopt_success, sqpopt_user_requested_stop
    use sqpopt_kinds,             only: wp => sqpopt_module_wp

    implicit none

    type :: my_data
        !! an example user data object
        integer  :: n_f = 0, n_g = 0, n_c = 0, n_jac = 0, n_report = 0 !! calls seen by each function
        integer  :: n_refused = 0      !! evaluations refused with `status > 0`
        integer  :: stop_after = -1    !! if > 0, the objective asks to stop on this call
    end type my_data

    type(sqpopt_type)            :: solver
    type(sqpopt_problem_type)    :: problem
    type(sqpopt_options_type)    :: options
    type(sqpopt_linesearch_type) :: linesearch
    type(sqpopt_qp_solver_type)  :: qp_solver
    type(my_data), target        :: data
    real(wp) :: xsol(2), lam(1)
    integer  :: istat

    write(*,*) '----------------------------'
    write(*,*) 'test_callbacks'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-10.0_wp,-10.0_wp], x_ub=[10.0_wp,10.0_wp], c_lb=[-1.0e20_wp], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv, data=data)

    ! (step caps off, so the first step reaches the region where f is refused)
    linesearch%major_step_limit = huge(1.0_wp)
    qp_solver%max_step = huge(1.0_wp)

    ! ---- user data, and status > 0 ----
    call solver%initialize(problem=problem, options=options, linesearch=linesearch, qp_solver=qp_solver, report=report)
    call solver%solve([0.0_wp, 0.0_wp], istat)
    call solver%get_solution(xsol, lam)
    print '(A,2F10.6,A,I0,A,6I4)', 'x=', xsol, ' istat=', istat, ' calls f,g,c,jac,report,refused:', &
        data%n_f, data%n_g, data%n_c, data%n_jac, data%n_report, data%n_refused
    if (istat /= sqpopt_success) error stop 'test_callbacks FAILED: did not converge'
    if (maxval(abs(xsol - [2.5_wp, 0.5_wp])) > 1.0e-5_wp) error stop 'test_callbacks FAILED: wrong solution'
    if (data%n_refused == 0) error stop 'test_callbacks FAILED: the refused region was never reached'
    if (min(data%n_f, data%n_g, data%n_c, data%n_jac, data%n_report) == 0) &
        error stop 'test_callbacks FAILED: user data not passed to every function'

    ! ---- status < 0: the objective asks to stop on its 3rd call ----
    data = my_data(stop_after=3)
    call solver%solve([0.0_wp, 0.0_wp], istat)
    print '(A,I0,A,I0,2A)', 'stop request: istat=', istat, ' f calls=', data%n_f, '  ', solver%status_message()
    if (istat /= sqpopt_user_requested_stop) error stop 'test_callbacks FAILED: stop request ignored'
    if (data%n_f /= 3) error stop 'test_callbacks FAILED: objective called again after the stop request'

    print '(A)', 'test_callbacks PASSED'

    contains

    subroutine obj(x, f, status, data)
    !! the objective
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = (x(1)-3.0_wp)**2 + (x(2)-1.0_wp)**2
    select type (data)
    type is (my_data)
        data%n_f = data%n_f + 1
        if (x(1) > 2.8_wp) then
            status = 1                          ! can't evaluate here
            data%n_refused = data%n_refused + 1
        end if
        if (data%n_f == data%stop_after) status = -1
    end select
    end subroutine obj

    subroutine grad(x, g, status, data)
    !! the objective's gradient
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    g = [2.0_wp*(x(1)-3.0_wp), 2.0_wp*(x(2)-1.0_wp)]
    select type (data)
    type is (my_data)
        data%n_g = data%n_g + 1
    end select
    associate(unused => status); end associate
    end subroutine grad

    subroutine cons(x, c, status, data)
    !! the constraints
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    c(1) = x(1) + x(2)
    select type (data)
    type is (my_data)
        data%n_c = data%n_c + 1
    end select
    associate(unused => status); end associate
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    !! the nonzero values of the constraint Jacobian (in the sparsity pattern's order)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    jac_val = 1.0_wp
    select type (data)
    type is (my_data)
        data%n_jac = data%n_jac + 1
    end select
    associate(unused => x); end associate
    associate(unused => status); end associate
    end subroutine jacv

    subroutine report(iter, x, f, c, lambda, user_stop, data)
    !! the `report` callback: counts its calls in the user data
    integer,                intent(in)    :: iter      !! major iteration number
    real(wp), dimension(:), intent(in)    :: x         !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: c         !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:), intent(in)    :: lambda    !! constraint multipliers `dimension(m)`
    real(wp),               intent(in)    :: f         !! objective value at `x`
    logical,                intent(out)   :: user_stop !! set to `.true.` to stop the solver
    class(*), optional,     intent(inout) :: data      !! the user data passed to `set_functions` (if any)
    user_stop = .false.
    select type (data)
    type is (my_data)
        data%n_report = data%n_report + 1
    end select
    associate(unused => iter); end associate
    associate(unused => x); end associate
    associate(unused => f); end associate
    associate(unused => c); end associate
    associate(unused => lambda); end associate
    end subroutine report

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

    subroutine gjac_grad_jacv(x, g, jac_val, accuracy, status, data)
    !! `gjac` for `set_functions`: the gradient (`grad`) and the Jacobian values (`jacv`)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g       !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    associate(unused => accuracy); end associate
    call grad(x, g, status, data)
    if (status == 0) call jacv(x, jac_val, status, data)
    end subroutine gjac_grad_jacv


end program test_callbacks
