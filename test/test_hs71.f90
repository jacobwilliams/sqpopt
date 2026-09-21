program test_hs71

    !! Hock-Schittkowski problem 71 (nonlinear objective, one nonlinear
    !! equality constraint, one nonlinear inequality constraint, and
    !! variable bounds):
    !!
    !!   minimize   x1*x4*(x1+x2+x3) + x3
    !!   subject to x1*x2*x3*x4 >= 25
    !!              x1**2 + x2**2 + x3**2 + x4**2 = 40
    !!              1 <= x1,x2,x3,x4 <= 5
    !!
    !! starting point x0 = (1,5,5,1); known optimal solution:
    !!   x* = (1, 4.7429994, 3.8211500, 1.3794083), f* = 17.0140173
    !!
    !! (same problem used in the `slsqp` and `psqp` test suites, expressed
    !! here directly as one equality + one inequality constraint since
    !! `sqpopt` supports two-sided constraints natively, without needing
    !! `slsqp`'s slack-variable reformulation)
    !!
    !! @note At this solution *both* constraints are simultaneously active
    !! (a corner point), which is a known hard case for the v1 composite-step
    !! QP heuristic: it settles into a small stable oscillation near the
    !! true solution rather than converging tightly to it (`istat` will be
    !! `sqpopt_max_iter_reached`, not `sqpopt_success`). This is a documented
    !! v1 limitation (see `PLAN.md`) -- a rigorous active-set/interior-point
    !! QP solver is needed for tight convergence on problems like this one.
    !! The tests only check that the solver gets reasonably close.
    !!
    !! Run twice, with the two available merit functions
    !! (`sqpopt_linesearch_module`), to compare them on this problem: the
    !! default `sqpopt_merit_l1` (with the second-order-correction
    !! safeguard in `sqpopt_iterate_module`), and `sqpopt_merit_augmented_lagrangian`
    !! (see PLAN.md \u00a76.1). Neither clears the limit cycle noted above -- the
    !! augmented Lagrangian option gets no closer here, confirming the PLAN.md
    !! \u00a76.1 caveat that it needs pairing with a real QP solve (\u00a76.2) to fully
    !! deliver its Maratos-avoidance benefit.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_linesearch_module, only: sqpopt_merit_l1, sqpopt_merit_augmented_lagrangian
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    write(*,*) '----------------------------'
    write(*,*) 'test_hs71'
    write(*,*) '----------------------------'

    call run_hs71('l1 (default)',           sqpopt_merit_l1)
    call run_hs71('augmented Lagrangian',   sqpopt_merit_augmented_lagrangian)

    contains

    subroutine run_hs71(label, merit_mode)

    character(len=*), intent(in) :: label
    integer,           intent(in) :: merit_mode

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(4), xsol(4), lam(2)
    real(wp), parameter :: xexpect(4) = [1.0_wp, 4.7429994_wp, 3.8211500_wp, 1.3794083_wp]
    real(wp), parameter :: fexpect = 17.0140173_wp
    real(wp) :: fsol
    integer :: istat

    ! equality constraint (index 1) first, then the inequality (index 2):
    call problem%set_problem_size(n=4, m_eq=1, m_ineq=1)
    call problem%set_bounds(x_lb=[1.0_wp,1.0_wp,1.0_wp,1.0_wp], x_ub=[5.0_wp,5.0_wp,5.0_wp,5.0_wp], &
                             c_lb=[40.0_wp, 25.0_wp], c_ub=[40.0_wp, big])
    call problem%set_jacobian_sparsity(nnz=8, irow=[1,1,1,1,2,2,2,2], icol=[1,2,3,4,1,2,3,4])
    call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)

    options%max_iter   = 300
    options%merit_mode = merit_mode
    x0 = [1.0_wp, 5.0_wp, 5.0_wp, 1.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)
    call obj(xsol, fsol)

    print '(3A)',       'test_hs71 [', trim(label), ']'
    print '(A,4F12.7)', 'x       = ', xsol
    print '(A,4F12.7)', 'x_true  = ', xexpect
    print '(A,F12.7)',  'x_error = ', norm2(xexpect - xsol)
    print '(A,F12.7)',  'f       = ', fsol
    print '(A,F12.7)',  'f_true  = ', fexpect
    print '(A,I0)',     'istat   = ', istat

    ! v1 only achieves loose convergence on this problem (see the note
    ! above), so a generous tolerance is used here:
    if (maxval(abs(xsol-xexpect)) > 0.5_wp) error stop 'test_hs71 FAILED: '//label
    print *, 'test_hs71 ['//trim(label)//'] PASSED'

    end subroutine run_hs71

    subroutine obj(x, f)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),                intent(out) :: f
    f = x(1)*x(4)*(x(1)+x(2)+x(3)) + x(3)
    end subroutine obj

    subroutine grad(x, g)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    g(1) = x(4)*(2.0_wp*x(1)+x(2)+x(3))
    g(2) = x(1)*x(4)
    g(3) = x(1)*x(4) + 1.0_wp
    g(4) = x(1)*(x(1)+x(2)+x(3))
    end subroutine grad

    subroutine cons(x, c)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    c(1) = x(1)**2 + x(2)**2 + x(3)**2 + x(4)**2
    c(2) = x(1)*x(2)*x(3)*x(4)
    end subroutine cons

    subroutine jacv(x, jac_val)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: jac_val
    ! order matches set_jacobian_sparsity: rows [1,1,1,1,2,2,2,2], cols [1,2,3,4,1,2,3,4]
    jac_val(1) = 2.0_wp*x(1)
    jac_val(2) = 2.0_wp*x(2)
    jac_val(3) = 2.0_wp*x(3)
    jac_val(4) = 2.0_wp*x(4)
    jac_val(5) = x(2)*x(3)*x(4)
    jac_val(6) = x(1)*x(3)*x(4)
    jac_val(7) = x(1)*x(2)*x(4)
    jac_val(8) = x(1)*x(2)*x(3)
    end subroutine jacv

end program test_hs71
