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
    !! Run with both the (default) Armijo line search and the filter line
    !! search (whose failure path goes through feasibility restoration).
    !!
    !! The legacy composite-step heuristic (`sqpopt_qp_composite`) cannot tell
    !! an inconsistent linearization from a consistent one, so it is only
    !! required to stop with a failure status (not success, and not by
    !! running to `max_iter`).

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_auto, sqpopt_qp_composite, sqpopt_qp_dense, sqpopt_qp_reduced_hessian
    use sqpopt_types_module,     only: sqpopt_infeasible, sqpopt_success, sqpopt_max_iter_reached
    use sqpopt_linesearch_module, only: sqpopt_linesearch_armijo, sqpopt_linesearch_filter
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: modes(4) = [sqpopt_qp_auto, sqpopt_qp_composite, sqpopt_qp_dense, sqpopt_qp_reduced_hessian]

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: xsol(2), lam(2)
    integer  :: istat, i, ls

    write(*,*) '----------------------------'
    write(*,*) 'test_infeasible'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=2)
    call problem%set_bounds(x_lb=[-10.0_wp,-10.0_wp], x_ub=[10.0_wp,10.0_wp], &
                             c_lb=[0.0_wp,2.0_wp], c_ub=[1.0_wp,3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,2], icol=[1,1])
    call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)

    do ls = 1, 2
    do i = 1, size(modes)
        options = sqpopt_options_type()
        options%qp_solver_mode = modes(i)
        options%linesearch_mode = merge(sqpopt_linesearch_armijo, sqpopt_linesearch_filter, ls == 1)
        call solver%initialize(problem=problem, options=options)
        call solver%solve([0.5_wp, 0.0_wp], istat)
        call solver%get_solution(xsol, lam)
        print '(A,I0,A,I0,A,2F10.4,A,I0,2A)', 'linesearch_mode=', options%linesearch_mode, &
            ' qp_solver_mode=', modes(i), ': x=', xsol, &
            '  istat=', istat, '  ', solver%status_message()
        if (modes(i) == sqpopt_qp_composite) then
            if (istat == sqpopt_success .or. istat == sqpopt_max_iter_reached) &
                error stop 'test_infeasible FAILED: composite mode did not stop with a failure status'
        else
            if (istat /= sqpopt_infeasible) error stop 'test_infeasible FAILED: infeasibility not detected'
            if (abs(xsol(1)-1.5_wp) > 1.0e-4_wp) error stop 'test_infeasible FAILED: not at the least-infeasible point'
        end if
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

end program test_infeasible
