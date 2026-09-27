program test_degenerate

    !! Regression test: attraction to a point where the constraint
    !! qualification fails (see `options%elastic_multiplier_limit`).
    !!
    !!     minimize   x1 + x2
    !!     subject to x1^2 + x2^2 = 1,   -1 <= x <= 1
    !!
    !! The minimum is x* = -(1,1)/sqrt(2). From the box corners (-1,1) and
    !! (1,-1), the iterates run along a bound toward (-1,0) or (0,-1), where
    !! the circle is tangent to that bound: there the linearized constraint
    !! lets each step cover only half the remaining distance, the multiplier
    !! doubles every iteration, and the solver used to stop there (stalled,
    !! failed line search, or, with the sparse QP, a false success). The
    !! elastic re-solve lets the iterates leave. Each corner is solved with
    !! the filter, funnel, and Armijo line searches, and both QP solvers.
    !! (From (1,1), the iterates stay on the diagonal and converge to the
    !! maximum, a first-order KKT point; that corner is not tested.)

    use sqpopt_module,            only: sqpopt_type
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_dense, sqpopt_qp_reduced_hessian
    use sqpopt_linesearch_module, only: sqpopt_linesearch_filter, sqpopt_linesearch_funnel, sqpopt_linesearch_armijo
    use sqpopt_types_module,      only: sqpopt_success
    use sqpopt_kinds,             only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: tol = 1.0e-6_wp
    real(wp), dimension(2), parameter :: x_star = -[1.0_wp, 1.0_wp]/sqrt(2.0_wp)
    real(wp), dimension(2,3), parameter :: corners = reshape([-1.0_wp,1.0_wp, 1.0_wp,-1.0_wp, -1.0_wp,-1.0_wp], [2,3])
    integer,  dimension(3),   parameter :: modes = [sqpopt_linesearch_filter, sqpopt_linesearch_funnel, &
                                                    sqpopt_linesearch_armijo]
    character(len=6), dimension(3), parameter :: mode_names = ['filter', 'funnel', 'armijo']
    integer,  dimension(2),   parameter :: qps = [sqpopt_qp_dense, sqpopt_qp_reduced_hessian]
    character(len=6), dimension(2), parameter :: qp_names = ['dense ', 'sparse']

    integer :: ls, qp, k, n_fail

    write(*,*) '----------------------------'
    write(*,*) 'test_degenerate'
    write(*,*) '----------------------------'

    n_fail = 0
    do ls = 1, size(modes)
        do qp = 1, size(qps)
            do k = 1, size(corners, 2)
                call run(ls, qp, k)
            end do
        end do
    end do
    if (n_fail > 0) error stop 'test_degenerate FAILED'
    write(*,*) 'test_degenerate PASSED'

contains

    subroutine run(ls, qp, k)
    integer, intent(in) :: ls, qp, k
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x(2), lam(1)
    integer  :: istat
    logical  :: ok

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-1.0_wp,-1.0_wp], x_ub=[1.0_wp,1.0_wp], c_lb=[1.0_wp], c_ub=[1.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc, gjac=gjac)
    options%linesearch_mode = modes(ls)
    options%qp_solver_mode  = qps(qp)
    call solver%initialize(problem=problem, options=options)
    call solver%solve(corners(:,k), istat)
    call solver%get_solution(x, lam)

    ok = istat == sqpopt_success .and. maxval(abs(x - x_star)) <= tol
    print '(A,1X,A,A,2F6.2,A,2F10.6,A,I0,A,ES9.2,A)', mode_names(ls), qp_names(qp), '  x0=', corners(:,k), &
        '  x=', x, '  istat=', istat, '  lambda=', lam(1), merge('       ', '  <-- X', ok)
    if (.not. ok) n_fail = n_fail + 1
    end subroutine run

    subroutine fc(x, f, c, status, data)
    real(wp), dimension(:), intent(in)    :: x
    real(wp),               intent(out)   :: f
    real(wp), dimension(:), intent(out)   :: c
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    f    = x(1) + x(2)
    c(1) = x(1)**2 + x(2)**2
    end subroutine fc

    subroutine gjac(x, g, jac_val, status, data)
    real(wp), dimension(:), intent(in)    :: x
    real(wp), dimension(:), intent(out)   :: g
    real(wp), dimension(:), intent(out)   :: jac_val
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    g       = [1.0_wp, 1.0_wp]
    jac_val = [2.0_wp*x(1), 2.0_wp*x(2)]
    end subroutine gjac

end program test_degenerate
