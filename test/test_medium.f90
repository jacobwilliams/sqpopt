program test_medium

    !! A medium-scale nonlinear program (10 variables, 5 nonlinear
    !! constraints, plus bounds on every variable) used to exercise
    !! `sqpopt` at a larger scale than [[test_basic]] while remaining
    !! well-behaved (each constraint couples only a small, disjoint
    !! group of variables, so the KKT system decouples into 5
    !! independent, essentially convex sub-problems with a verifiable
    !! closed-form solution -- unlike `test_hs71`, which is a genuinely
    !! hard corner-point case). This gives a reliable "does the whole
    !! pipeline scale and still converge tightly" regression test.
    !!
    !!   minimize   (x1-1)^2  + (x2-2)^2  + (x3-2)^2            [block 1: sphere]
    !!            + (x4-3)^2  + (x5-3)^2                        [block 2: hyperbola]
    !!            + (x6-3)^2  + (x7-3)^2                        [block 3: disk]
    !!            + (x8-1)^2  + (x9-1)^2                        [block 4: hyperbola]
    !!            + (x10-5)^2                                   [block 5: box-active]
    !!
    !!   subject to c1: x1^2 + x2^2 + x3^2 = 36        (nonlinear equality)
    !!              c2: x4*x5             = 12         (nonlinear equality)
    !!              c3: x6^2 + x7^2      <= 8           (nonlinear inequality)
    !!              c4: x8*x9            >= 4           (nonlinear inequality)
    !!              c5: x10^2            <= 16          (nonlinear inequality)
    !!              0.1 <= x1,...,x9 <= 10,  0.1 <= x10 <= 3
    !!
    !! Each block's target point (the unconstrained minimizer of its
    !! quadratic objective piece) is chosen to lie outside its feasible
    !! region, so every constraint (or, for block 5, the bound) is
    !! forced active and the problem has one verifiable KKT point:
    !!
    !!  * block 1: minimizing ||x-t||^2 s.t. ||x||^2=R^2 (t=(1,2,2), R=6)
    !!    is solved by scaling t onto the sphere: x*=t*R/||t|| = (2,4,4).
    !!  * block 2: minimizing (x4-3)^2+(x5-3)^2 s.t. x4*x5=12; the
    !!    hyperbola branch in the positive quadrant is symmetric in the
    !!    (symmetric) target, giving x4*=x5*=sqrt(12) = 2*sqrt(3).
    !!  * block 3: minimizing ||x-t||^2 s.t. ||x||^2<=8 (t=(3,3)) with t
    !!    infeasible: same projection as block 1, x*=t*sqrt(8)/||t||=(2,2).
    !!  * block 4: minimizing (x8-1)^2+(x9-1)^2 s.t. x8*x9>=4; the region
    !!    {x*y>=4, x,y>0} is convex (it lies above a hyperbola), and by
    !!    the same symmetry argument as block 2, x8*=x9*=2.
    !!  * block 5: minimizing (x10-5)^2 s.t. 0.1<=x10<=3 (the nonlinear
    !!    inequality x10^2<=16 is deliberately slack here, so the upper
    !!    *bound* is what's actually active): x10*=3.
    !!
    !! giving x* = (2, 4, 4, 2*sqrt(3), 2*sqrt(3), 2, 2, 2, 2, 3) and
    !! f* = 9 + 2*(2*sqrt(3)-3)^2 + 2 + 2 + 4 = 17.4307808...

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type
    use sqpopt_types_module,   only: sqpopt_success
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides
    real(wp), parameter :: sqrt3 = 1.7320508075688772_wp

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_qp_solver_type) :: qp_solver
    real(wp) :: x0(10), xsol(10), lam(5)
    real(wp), parameter :: xexpect(10) = &
        [2.0_wp, 4.0_wp, 4.0_wp, 2.0_wp*sqrt3, 2.0_wp*sqrt3, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 3.0_wp]
    real(wp), parameter :: fexpect = 17.4307808_wp
    real(wp) :: fsol
    integer :: istat

    write(*,*) '----------------------------'
    write(*,*) 'test_medium'
    write(*,*) '----------------------------'

    ! equality constraints (1,2) first, then the inequalities (3,4,5):
    call problem%set_problem_size(n=10, m_eq=2, m_ineq=3)
    call problem%set_bounds( &
        x_lb=[0.1_wp,0.1_wp,0.1_wp,0.1_wp,0.1_wp,0.1_wp,0.1_wp,0.1_wp,0.1_wp,0.1_wp], &
        x_ub=[10.0_wp,10.0_wp,10.0_wp,10.0_wp,10.0_wp,10.0_wp,10.0_wp,10.0_wp,10.0_wp,3.0_wp], &
        c_lb=[36.0_wp, 12.0_wp, -big,   4.0_wp, -big], &
        c_ub=[36.0_wp, 12.0_wp,  8.0_wp, big,   16.0_wp])
    call problem%set_jacobian_sparsity(nnz=10, &
        irow=[1,1,1, 2,2, 3,3, 4,4, 5], &
        icol=[1,2,3, 4,5, 6,7, 8,9, 10])
    call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)

    options%max_iter = 300
    options%ktol     = 1.0e-4_wp
    options%ctol     = 1.0e-6_wp
    x0 = [1.0_wp,1.0_wp,1.0_wp, 1.0_wp,1.0_wp, 1.0_wp,1.0_wp, 1.0_wp,1.0_wp, 1.0_wp]

    ! the default max_step=2.0 trust-region cap applies to the *full* n=10
    ! step norm, which artificially couples the 5 otherwise-independent
    ! blocks together early on; widen it so each block can move freely:
    qp_solver%max_step = 5.0_wp

    call solver%initialize(problem=problem, options=options, qp_solver=qp_solver)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)
    call obj(xsol, fsol)

    print '(A,10F10.5)', 'test_medium: x       = ', xsol
    print '(A,10F10.5)', 'test_medium: x_true  = ', xexpect
    print '(A,F12.7)',   'test_medium: x_error = ', norm2(xexpect - xsol)
    print '(A,F12.7)',   'test_medium: f       = ', fsol
    print '(A,F12.7)',   'test_medium: f_true  = ', fexpect
    print '(A,I0)',      'test_medium: istat   = ', istat

    if (istat /= sqpopt_success) error stop 'test_medium FAILED: did not converge'
    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_medium FAILED: wrong solution'
    print *, 'test_medium PASSED'

    contains

    subroutine obj(x, f)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),                intent(out) :: f
    f = (x(1)-1.0_wp)**2 + (x(2)-2.0_wp)**2 + (x(3)-2.0_wp)**2 &
      + (x(4)-3.0_wp)**2 + (x(5)-3.0_wp)**2 &
      + (x(6)-3.0_wp)**2 + (x(7)-3.0_wp)**2 &
      + (x(8)-1.0_wp)**2 + (x(9)-1.0_wp)**2 &
      + (x(10)-5.0_wp)**2
    end subroutine obj

    subroutine grad(x, g)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    g(1)  = 2.0_wp*(x(1)-1.0_wp)
    g(2)  = 2.0_wp*(x(2)-2.0_wp)
    g(3)  = 2.0_wp*(x(3)-2.0_wp)
    g(4)  = 2.0_wp*(x(4)-3.0_wp)
    g(5)  = 2.0_wp*(x(5)-3.0_wp)
    g(6)  = 2.0_wp*(x(6)-3.0_wp)
    g(7)  = 2.0_wp*(x(7)-3.0_wp)
    g(8)  = 2.0_wp*(x(8)-1.0_wp)
    g(9)  = 2.0_wp*(x(9)-1.0_wp)
    g(10) = 2.0_wp*(x(10)-5.0_wp)
    end subroutine grad

    subroutine cons(x, c)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    c(1) = x(1)**2 + x(2)**2 + x(3)**2
    c(2) = x(4)*x(5)
    c(3) = x(6)**2 + x(7)**2
    c(4) = x(8)*x(9)
    c(5) = x(10)**2
    end subroutine cons

    subroutine jacv(x, jac_val)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: jac_val
    ! order matches set_jacobian_sparsity: rows [1,1,1,2,2,3,3,4,4,5]
    jac_val(1)  = 2.0_wp*x(1)
    jac_val(2)  = 2.0_wp*x(2)
    jac_val(3)  = 2.0_wp*x(3)
    jac_val(4)  = x(5)
    jac_val(5)  = x(4)
    jac_val(6)  = 2.0_wp*x(6)
    jac_val(7)  = 2.0_wp*x(7)
    jac_val(8)  = x(9)
    jac_val(9)  = x(8)
    jac_val(10) = 2.0_wp*x(10)
    end subroutine jacv

end program test_medium
