program test_infeasible

    !! Regression test: an infeasible problem must be reported as
    !! `sqpopt_infeasible` (not run to `max_iter`), for every QP solver mode.
    !!
    !!   minimize   x1^2 + x2^2
    !!   subject to 0 <= c1 = x1 <= 1
    !!              2 <= c2 = x1 <= 3
    !!
    !! The closest the constraints can get to feasibility (in the l2 sense
    !! the infeasibility test uses) is x1 = 1.5.
    !!
    !! Run with the Armijo, filter, and funnel line searches (the latter two's
    !! failure path goes through feasibility restoration), and with the
    !! trust region and the filter (whose failure path goes through the
    !! restoration phase: without it, i.e. with `sqpopt_restoration_gauss_newton`,
    !! the trust region stops at x1 = 1 with `sqpopt_line_search_failed`).
    !! The closest the constraints can get to feasibility in the l1 sense
    !! is any x1 in [1,2], so the restoration phase (which minimizes the l1
    !! violation) must hand over to the Gauss-Newton step to reach x1 = 1.5.

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_auto, sqpopt_qp_dense, sqpopt_qp_reduced_hessian
    use sqpopt_types_module,     only: sqpopt_infeasible
    use sqpopt_linesearch_module, only: sqpopt_linesearch_armijo, sqpopt_linesearch_filter, sqpopt_linesearch_funnel
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: modes(3) = [sqpopt_qp_auto, sqpopt_qp_dense, sqpopt_qp_reduced_hessian]

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_trust_region_type) :: trust_region
    character(len=*), parameter :: labels(4) = [character(len=24) :: 'armijo', 'filter', 'funnel', &
                                                'trust region + filter']
    real(wp) :: xsol(2), lam(2)
    integer  :: istat, i, ls

    write(*,*) '----------------------------'
    write(*,*) 'test_infeasible'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m=2)
    call problem%set_bounds(x_lb=[-10.0_wp,-10.0_wp], x_ub=[10.0_wp,10.0_wp], &
                             c_lb=[0.0_wp,2.0_wp], c_ub=[1.0_wp,3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,2], icol=[1,1])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv)

    do ls = 1, size(labels)
    do i = 1, size(modes)
        options = sqpopt_options_type()
        options%qp_solver_mode = modes(i)
        select case (ls)
        case (1);    options%linesearch_mode = sqpopt_linesearch_armijo
        case (3);    options%linesearch_mode = sqpopt_linesearch_funnel
        case default; options%linesearch_mode = sqpopt_linesearch_filter
        end select
        trust_region = sqpopt_trust_region_type()
        trust_region%enabled = ls >= 4
        call solver%initialize(problem=problem, options=options, trust_region=trust_region)
        call solver%solve([0.5_wp, 0.0_wp], istat)
        call solver%get_solution(xsol, lam)
        print '(A,A,I0,A,2F10.4,A,I0,2A)', labels(ls), &
            ' qp_solver_mode=', modes(i), ': x=', xsol, &
            '  istat=', istat, '  ', solver%status_message()
        if (istat /= sqpopt_infeasible) error stop 'test_infeasible FAILED: infeasibility not detected'
        if (abs(xsol(1)-1.5_wp) > 1.0e-4_wp) error stop 'test_infeasible FAILED: not at the least-infeasible point'
    end do
    end do

    print '(A)', 'test_infeasible PASSED'

    contains

    subroutine obj(x, f, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),               intent(out) :: f
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    f = x(1)**2 + x(2)**2
    end subroutine obj

    subroutine grad(x, g, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    g = 2.0_wp*x
    end subroutine grad

    subroutine cons(x, c, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    c = x(1)
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


end program test_infeasible
