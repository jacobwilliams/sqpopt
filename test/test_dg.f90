program test_dg

    !! example from A. Duran and I.E. Grossmann
    !! [taken from Sect 5 of FilterSQP manual]
    !! https://www.mcs.anl.gov/~leyffer/papers/SQP_manual.pdf

    !! this will use an operator overloading approach for automatic differentiation

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type
    use sqpopt_types_module,     only: sqpopt_success
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    integer,parameter :: n = 6
    integer,parameter :: m = 6

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_qp_solver_type) :: qp_solver
    real(wp),dimension(n) :: x0, xsol
    real(wp),dimension(m) :: lam
    real(wp) :: fsol
    integer :: istat

    real(wp),dimension(*),parameter :: x_lb = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    real(wp),dimension(*),parameter :: x_ub = [2.0_wp, 2.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
    real(wp),dimension(*),parameter :: c_lb = [0.0_wp, -2.0_wp, -big, -big, -big, -big]
    real(wp),dimension(*),parameter :: c_ub = [big, big, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]

    ! solution:
    real(wp), dimension(n), parameter :: xexpect = &
        [1.147_wp, 0.547_wp, 1.0_wp, 0.273_wp, 0.300_wp, 0.000_wp] !! [x,y]
    real(wp), parameter :: fexpect = 0.759_wp

    integer,parameter :: perturb_mode = 1  !! perturbation is `dx=dpert`
    integer,parameter :: cache_size = 0    !! `0` indicates not to use cache
    integer,parameter :: sparsity_mode = 4 !!
    real(wp), dimension(n), parameter :: dpert = &
        [1.0e-6_wp, 1.0e-6_wp, 1.0e-6_wp, 1.0e-6_wp, 1.0e-6_wp, 1.0e-6_wp]

    write(*,*) '----------------------------'
    write(*,*) 'test_dg'
    write(*,*) '----------------------------'

    return ! not finished yet


    ! ! equality constraints (1,2) first, then the inequalities (3,4,5):
    ! call problem%set_problem_size(n=n, m=m)
    ! call problem%set_bounds(x_lb=x_lb, x_ub=x_ub, c_lb=c_lb, c_ub=c_ub)

    ! ! call problem%set_jacobian_sparsity(nnz=10, &
    ! !     irow=[1,1,1, 2,2, 3,3, 4,4, 5], &
    ! !     icol=[1,2,3, 4,5, 6,7, 8,9, 10])
    ! call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)

    ! options%max_iter = 300
    ! options%ktol     = 1.0e-4_wp
    ! options%ctol     = 1.0e-6_wp
    ! x0 = [1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp] ! initial guess

    ! qp_solver%max_step = 5.0_wp

    ! call solver%initialize(problem=problem, options=options, qp_solver=qp_solver)
    ! call solver%solve(x0, istat)
    ! call solver%get_solution(xsol, lam)
    ! call obj(xsol, fsol)

    ! print '(A,10F10.5)', 'test_medium: x       = ', xsol
    ! print '(A,10F10.5)', 'test_medium: x_true  = ', xexpect
    ! print '(A,F12.7)',   'test_medium: x_error = ', norm2(xexpect - xsol)
    ! print '(A,F12.7)',   'test_medium: f       = ', fsol
    ! print '(A,F12.7)',   'test_medium: f_true  = ', fexpect
    ! print '(A,I0)',      'test_medium: istat   = ', istat

    ! if (istat /= sqpopt_success) error stop 'test_medium FAILED: did not converge'
    ! if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_medium FAILED: wrong solution'
    ! print *, 'test_medium PASSED'

    ! contains

    ! subroutine obj(x, f)
    ! real(wp), dimension(:), intent(in)  :: x
    ! real(wp),               intent(out) :: f

    ! associate(x1 => x(1), x2 => x(2), x3 => x(3), &
    !           y1 => x(4), y2 => x(5), y3 => x(6))
    !     f = 5*y1 + 6*y2 + 8*y3 + 10*x1 - 7*x3 - 18*log(x2 + 1) &
    !         - 19.2_wp*log(x1-x2+1) + 10.0_wp
    ! end associate
    ! end subroutine obj

    ! subroutine grad(x, g)
    ! real(wp), dimension(:), intent(in)  :: x
    ! real(wp), dimension(:), intent(out) :: g

    ! ! use an operator overloading approach for automatic differentiation

    ! ! use numdiff to compute
    ! ! g(1)  = 2.0_wp*(x(1)-1.0_wp)
    ! ! g(2)  = 2.0_wp*(x(2)-2.0_wp)
    ! ! g(3)  = 2.0_wp*(x(3)-2.0_wp)
    ! ! g(4)  = 2.0_wp*(x(4)-3.0_wp)
    ! ! g(5)  = 2.0_wp*(x(5)-3.0_wp)
    ! ! g(6)  = 2.0_wp*(x(6)-3.0_wp)
    ! ! g(7)  = 2.0_wp*(x(7)-3.0_wp)
    ! ! g(8)  = 2.0_wp*(x(8)-1.0_wp)
    ! ! g(9)  = 2.0_wp*(x(9)-1.0_wp)
    ! ! g(10) = 2.0_wp*(x(10)-5.0_wp)
    ! end subroutine grad

    ! subroutine cons(x, c)
    ! real(wp), dimension(:), intent(in)  :: x
    ! real(wp), dimension(:), intent(out) :: c

    ! associate(x1 => x(1), x2 => x(2), x3 => x(3), &
    !           y1 => x(4), y2 => x(5), y3 => x(6))
    !     c(1) = 0.8_wp * log(x2+1) + 0.96*log(x1-x2+1) - 0.8_wp*x3
    !     c(2) = log(x2+1) + 1.2_wp*log(x1-x2+1) - x3 - 2*y3
    !     c(3) = x2-x1
    !     c(4) = x2-2*y1
    !     c(5) = x1-x2-2*y2
    !     c(6) = y1+y2
    ! end associate

    ! end subroutine cons

    ! subroutine jacv(x, jac_val)
    ! real(wp), dimension(:), intent(in)  :: x
    ! real(wp), dimension(:), intent(out) :: jac_val

    ! ! use numdiff to compute...
    ! ! order matches set_jacobian_sparsity: rows [1,1,1,2,2,3,3,4,4,5]
    ! ! jac_val(1)  = 2.0_wp*x(1)
    ! ! jac_val(2)  = 2.0_wp*x(2)
    ! ! jac_val(3)  = 2.0_wp*x(3)
    ! ! jac_val(4)  = x(5)
    ! ! jac_val(5)  = x(4)
    ! ! jac_val(6)  = 2.0_wp*x(6)
    ! ! jac_val(7)  = 2.0_wp*x(7)
    ! ! jac_val(8)  = x(9)
    ! ! jac_val(9)  = x(8)
    ! ! jac_val(10) = 2.0_wp*x(10)
    ! end subroutine jacv

end program test_dg
