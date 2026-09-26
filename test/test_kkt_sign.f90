program test_kkt_sign

    !! Regression test: a starting point where the constraints are active
    !! with *wrong-sign* multipliers must not be reported as converged.
    !!
    !!   minimize   -x1 - x2
    !!   subject to  0 <= c1 = x1 <= 10
    !!               0 <= c2 = x2 <= 10
    !!
    !! started on the lower bounds, x0 = (0,0). There `g = J^T lambda` holds
    !! with `lambda = (-1,-1)`, but that is the wrong sign at a lower bound,
    !! so x0 is not a KKT point; the solution is x* = (10,10).
    !! Every QP solver mode must find it.

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_auto, sqpopt_qp_composite, sqpopt_qp_dense, sqpopt_qp_reduced_hessian
    use sqpopt_types_module,     only: sqpopt_success
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: modes(4) = [sqpopt_qp_auto, sqpopt_qp_composite, sqpopt_qp_dense, sqpopt_qp_reduced_hessian]
    real(wp), parameter :: xexpect(2) = [10.0_wp, 10.0_wp]

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: xsol(2), lam(2)
    integer  :: istat, i

    write(*,*) '----------------------------'
    write(*,*) 'test_kkt_sign'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=2)
    call problem%set_bounds(x_lb=[-100.0_wp,-100.0_wp], x_ub=[100.0_wp,100.0_wp], &
                             c_lb=[0.0_wp,0.0_wp], c_ub=[10.0_wp,10.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,2], icol=[1,2])
    call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)

    do i = 1, size(modes)
        options = sqpopt_options_type()
        options%qp_solver_mode = modes(i)
        call solver%initialize(problem=problem, options=options)
        call solver%solve([0.0_wp, 0.0_wp], istat)
        call solver%get_solution(xsol, lam)
        print '(A,I0,A,2F10.4,A,2F10.4,A,I0)', 'qp_solver_mode=', modes(i), ': x=', xsol, '  lambda=', lam, '  istat=', istat
        if (istat /= sqpopt_success) error stop 'test_kkt_sign FAILED: did not reach sqpopt_success'
        if (maxval(abs(xsol-xexpect)) > 1.0e-6_wp) error stop 'test_kkt_sign FAILED: wrong solution'
    end do

    print '(A)', 'test_kkt_sign PASSED'

    contains

    subroutine obj(x, f, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),               intent(out) :: f
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    f = -x(1) - x(2)
    end subroutine obj

    subroutine grad(x, g, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    associate(unused => x); end associate
    g = -1.0_wp
    end subroutine grad

    subroutine cons(x, c, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    c = x
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: jac_val
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    associate(unused => x); end associate
    jac_val = 1.0_wp
    end subroutine jacv

end program test_kkt_sign
