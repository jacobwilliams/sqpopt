program test_resolve

    !! Regression test: calling `solve` twice on the same solver object must
    !! give identical results (same solution, status, and number of function
    !! evaluations), i.e. no state (penalty parameter, filter, watchdog
    !! best-point, trust-region radius, Hessian) may carry over from one
    !! solve to the next. Checked for every line search mode, with and
    !! without trust-region globalization.
    !!
    !!   minimize   (1-x1)^2 + 100*(x2-x1^2)^2    (Rosenbrock)
    !!   subject to x1^2 + x2^2 <= 1.5

    use sqpopt_module,              only: sqpopt_type
    use sqpopt_problem_module,      only: sqpopt_problem_type
    use sqpopt_options_module,      only: sqpopt_options_type
    use sqpopt_linesearch_module,   only: sqpopt_linesearch_armijo, sqpopt_linesearch_exact, &
                                          sqpopt_linesearch_watchdog, sqpopt_linesearch_filter, &
                                          sqpopt_linesearch_funnel
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_kinds,               only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: ls_modes(5) = [sqpopt_linesearch_armijo, sqpopt_linesearch_exact, &
                                         sqpopt_linesearch_watchdog, sqpopt_linesearch_filter, &
                                         sqpopt_linesearch_funnel]

    type(sqpopt_type)              :: solver
    type(sqpopt_problem_type)      :: problem
    type(sqpopt_options_type)      :: options
    type(sqpopt_trust_region_type) :: trust_region
    real(wp) :: x1(2), x2(2), lam(1)
    integer  :: istat1, istat2, nf1, nf2, i, tr
    integer  :: n_f  !! number of objective evaluations

    write(*,*) '----------------------------'
    write(*,*) 'test_resolve'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-5.0_wp,-5.0_wp], x_ub=[5.0_wp,5.0_wp], c_lb=[-1.0e20_wp], c_ub=[1.5_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv)

    do tr = 0, 1
        do i = 1, size(ls_modes)
            options = sqpopt_options_type()
            options%linesearch_mode = ls_modes(i)
            trust_region = sqpopt_trust_region_type()
            trust_region%enabled = tr == 1
            call solver%initialize(problem=problem, options=options, trust_region=trust_region)

            n_f = 0
            call solver%solve([-1.2_wp, 1.0_wp], istat1)
            call solver%get_solution(x1, lam)
            nf1 = n_f

            n_f = 0
            call solver%solve([-1.2_wp, 1.0_wp], istat2)
            call solver%get_solution(x2, lam)
            nf2 = n_f

            print '(A,L1,A,I0,A,2I4,A,2I6,A,ES10.2)', 'trust_region=', tr == 1, ' linesearch_mode=', ls_modes(i), &
                ': istat=', istat1, istat2, '  n_f=', nf1, nf2, '  |dx|=', maxval(abs(x1-x2))
            if (istat1 /= istat2 .or. nf1 /= nf2 .or. any(x1 /= x2)) &
                error stop 'test_resolve FAILED: second solve differs from the first'
        end do
    end do

    print '(A)', 'test_resolve PASSED'

    contains

    subroutine obj(x, f, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),               intent(out) :: f
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    n_f = n_f + 1
    f = (1.0_wp-x(1))**2 + 100.0_wp*(x(2)-x(1)**2)**2
    end subroutine obj

    subroutine grad(x, g, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    g(1) = -2.0_wp*(1.0_wp-x(1)) - 400.0_wp*x(1)*(x(2)-x(1)**2)
    g(2) = 200.0_wp*(x(2)-x(1)**2)
    end subroutine grad

    subroutine cons(x, c, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    c(1) = x(1)**2 + x(2)**2
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: jac_val
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    jac_val = 2.0_wp*x
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


end program test_resolve
