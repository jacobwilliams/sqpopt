program test_nonfinite

    !! Regression test: non-finite (NaN) function values are handled.
    !!
    !! * The objective returns NaN for `x1 > 2.8` (e.g. a domain error in the
    !!   user's model), and the first full QP step from `x0=(0,0)` lands at
    !!   `x1 = 3.5`, in that region. Every line search mode, and the trust
    !!   region, must reject that trial point and backtrack/shrink to a
    !!   finite one, then converge normally:
    !!
    !!     minimize   (x1-3)^2 + (x2-1)^2
    !!     subject to x1 + x2 <= 3
    !!
    !!   known solution: x* = (2.5, 0.5)
    !!
    !!   (the `major_step_limit` and `max_step` step caps are disabled so the
    !!   first trial really is the full step; by default they would already
    !!   shorten it.)
    !!
    !! * A NaN at the starting point must stop the solver with
    !!   `sqpopt_function_error`.
    !!
    !! (NaN comparisons are never made by the solver -- it checks values with
    !! `ieee_is_finite` -- so this also passes with `-ffpe-trap=invalid`.)

    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use sqpopt_module,              only: sqpopt_type
    use sqpopt_problem_module,      only: sqpopt_problem_type
    use sqpopt_options_module,      only: sqpopt_options_type
    use sqpopt_linesearch_module,   only: sqpopt_linesearch_type, sqpopt_linesearch_armijo, sqpopt_linesearch_exact, &
                                          sqpopt_linesearch_watchdog, sqpopt_linesearch_filter, &
                                          sqpopt_linesearch_funnel
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_qp_solver_module,    only: sqpopt_qp_solver_type
    use sqpopt_types_module,        only: sqpopt_success, sqpopt_function_error
    use sqpopt_kinds,               only: wp => sqpopt_module_wp

    implicit none

    integer,  parameter :: ls_modes(5) = [sqpopt_linesearch_armijo, sqpopt_linesearch_exact, &
                                          sqpopt_linesearch_watchdog, sqpopt_linesearch_filter, &
                                          sqpopt_linesearch_funnel]
    real(wp), parameter :: xexpect(2) = [2.5_wp, 0.5_wp]

    type(sqpopt_type)              :: solver
    type(sqpopt_problem_type)      :: problem
    type(sqpopt_options_type)      :: options
    type(sqpopt_trust_region_type) :: trust_region
    type(sqpopt_linesearch_type)   :: linesearch
    type(sqpopt_qp_solver_type)    :: qp_solver
    real(wp) :: xsol(2), lam(1)
    integer  :: istat, i, tr
    integer  :: n_nan  !! number of objective evaluations that returned NaN

    write(*,*) '----------------------------'
    write(*,*) 'test_nonfinite'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-10.0_wp,-10.0_wp], x_ub=[10.0_wp,10.0_wp], c_lb=[-1.0e20_wp], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv)

    do tr = 0, 1
        do i = 1, size(ls_modes)
            options = sqpopt_options_type()
            options%linesearch_mode = ls_modes(i)
            trust_region = sqpopt_trust_region_type()
            trust_region%enabled = tr == 1
            trust_region%radius0 = 10.0_wp  !! so the first trust-region step also reaches the NaN region
            linesearch = sqpopt_linesearch_type()
            linesearch%major_step_limit = huge(1.0_wp)
            qp_solver = sqpopt_qp_solver_type()
            qp_solver%max_step = huge(1.0_wp)
            call solver%initialize(problem=problem, options=options, trust_region=trust_region, &
                                   linesearch=linesearch, qp_solver=qp_solver)
            n_nan = 0
            call solver%solve([0.0_wp, 0.0_wp], istat)
            call solver%get_solution(xsol, lam)
            print '(A,L1,A,I0,A,2F10.6,A,I0,A,I0)', 'trust_region=', tr == 1, ' linesearch_mode=', ls_modes(i), &
                ': x=', xsol, '  istat=', istat, '  NaN evaluations=', n_nan
            if (n_nan == 0) error stop 'test_nonfinite FAILED: the NaN region was never reached (test is not meaningful)'
            if (istat /= sqpopt_success) error stop 'test_nonfinite FAILED: did not reach sqpopt_success'
            if (maxval(abs(xsol-xexpect)) > 1.0e-5_wp) error stop 'test_nonfinite FAILED: wrong solution'
        end do
    end do

    ! a NaN at the starting point:
    call solver%initialize(problem=problem)
    call solver%solve([2.9_wp, 0.0_wp], istat)
    print '(A,I0,2A)', 'NaN at x0: istat=', istat, '  ', solver%status_message()
    if (istat /= sqpopt_function_error) error stop 'test_nonfinite FAILED: NaN at x0 not reported'

    print '(A)', 'test_nonfinite PASSED'

    contains

    subroutine obj(x, f, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),               intent(out) :: f
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    if (x(1) > 2.8_wp) then
        f = ieee_value(f, ieee_quiet_nan)
        n_nan = n_nan + 1
    else
        f = (x(1)-3.0_wp)**2 + (x(2)-1.0_wp)**2
    end if
    end subroutine obj

    subroutine grad(x, g, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    g = [2.0_wp*(x(1)-3.0_wp), 2.0_wp*(x(2)-1.0_wp)]
    end subroutine grad

    subroutine cons(x, c, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    c(1) = x(1) + x(2)
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: jac_val
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    associate(unused => x); end associate
    jac_val = 1.0_wp
    end subroutine jacv

    subroutine fc_obj_cons(x, f, c, status, data)
    !! `fc` for `set_functions`: the objective (`obj`) and the constraints (`cons`)
    real(wp), dimension(:), intent(in)    :: x
    real(wp),               intent(out)   :: f
    real(wp), dimension(:), intent(out)   :: c
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    call obj(x, f, status, data)
    if (status == 0) call cons(x, c, status, data)
    end subroutine fc_obj_cons

    subroutine gjac_grad_jacv(x, g, jac_val, status, data)
    !! `gjac` for `set_functions`: the gradient (`grad`) and the Jacobian values (`jacv`)
    real(wp), dimension(:), intent(in)    :: x
    real(wp), dimension(:), intent(out)   :: g
    real(wp), dimension(:), intent(out)   :: jac_val
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    call grad(x, g, status, data)
    if (status == 0) call jacv(x, jac_val, status, data)
    end subroutine gjac_grad_jacv


end program test_nonfinite
