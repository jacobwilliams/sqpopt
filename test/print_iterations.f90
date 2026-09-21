program print_iterations

    !! Demonstrates the optional `report` callback (see [[sqpopt_types_module]]
    !! for the `sqpopt_report_func` interface): a user-supplied subroutine
    !! that `sqpopt` calls once per major iteration, here just printing the
    !! iteration number, point, objective, constraints, and multipliers to
    !! the console. Run with `fpm run --example print_iterations`.
    !!
    !! Problem: minimize (x1-2)^2 + (x2-3)^2  s.t.  x1+x2 = 4
    !! known solution: x* = (1.5, 2.5), f* = 0.5

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    integer :: istat

    write(*,*) '----------------------------'
    write(*,*) 'print_iterations'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m_eq=1, m_ineq=0)
    call problem%set_bounds(x_lb=[-huge(1.0_wp),-huge(1.0_wp)], x_ub=[huge(1.0_wp),huge(1.0_wp)], &
                             c_lb=[4.0_wp], c_ub=[4.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)

    options%max_iter = 100
    x0 = [0.0_wp, 0.0_wp]

    print '(A5,2A12,A12,A12)', 'iter', 'x1', 'x2', 'f', 'c1'
    call solver%initialize(problem=problem, options=options, report=report)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'solution: x = ', xsol
    print '(A,I0)',     'istat       = ', istat

    contains

    subroutine report(iter, x, f, c, lambda, user_stop)
    integer,                intent(in)  :: iter
    real(wp), dimension(:), intent(in)  :: x
    real(wp),               intent(in)  :: f
    real(wp), dimension(:), intent(in)  :: c
    real(wp), dimension(:), intent(in)  :: lambda
    logical,                intent(out) :: user_stop
    print '(I5,2F12.6,F12.6,F12.6)', iter, x, f, c
    user_stop = .false.
    end subroutine report

    subroutine obj(x, f)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),                intent(out) :: f
    f = (x(1)-2.0_wp)**2 + (x(2)-3.0_wp)**2
    end subroutine obj

    subroutine grad(x, g)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    g(1) = 2.0_wp*(x(1)-2.0_wp)
    g(2) = 2.0_wp*(x(2)-3.0_wp)
    end subroutine grad

    subroutine cons(x, c)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    c(1) = x(1) + x(2)
    end subroutine cons

    subroutine jacv(x, jac_val)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: jac_val
    jac_val = [1.0_wp, 1.0_wp]
    end subroutine jacv

end program print_iterations
