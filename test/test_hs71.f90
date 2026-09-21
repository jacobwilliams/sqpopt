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
    !!
    !! Run five ways, to compare the available merit functions, line
    !! searches, and QP solvers on this problem: the default
    !! `sqpopt_merit_l1` (with the second-order-correction safeguard in
    !! `sqpopt_iterate_module`), `sqpopt_merit_augmented_lagrangian` (see
    !! PLAN.md section 6.1), and `sqpopt_linesearch_watchdog` (Powell's
    !! VF13 watchdog technique, see PLAN.md section 6.3) all still use the
    !! v1 composite-step QP and only achieve loose convergence (the
    !! watchdog line search gets noticeably closer than the other two, but
    !! none reach `sqpopt_success`). The last two, `sqpopt_qp_dense` (see
    !! `DENSE_QP_PLAN.md`) and `sqpopt_qp_reduced_hessian` (see
    !! `REDUCED_HESSIAN_QP_PLAN.md`), replace the v1 QP heuristic with a
    !! real active-set QP solve (dense and sparse/matrix-free,
    !! respectively) and **both** reach `sqpopt_success`, converging to
    !! the known solution tightly -- confirming the recurring conclusion
    !! (PLAN.md sections 6.1/6.3) that a real QP solve, not another
    !! merit-function/line-search patch, is what was needed here.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_linesearch_module, only: sqpopt_merit_l1, sqpopt_merit_augmented_lagrangian, &
                                        sqpopt_linesearch_armijo, sqpopt_linesearch_watchdog
    use sqpopt_qp_solver_module, only: sqpopt_qp_composite, sqpopt_qp_dense, sqpopt_qp_reduced_hessian
    use sqpopt_types_module,   only: sqpopt_success
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    integer :: i_obj, i_grad, i_cons, i_jac

    write(*,*) '----------------------------'
    write(*,*) 'test_hs71'
    write(*,*) '----------------------------'

    call run_hs71('l1 (default)',           sqpopt_merit_l1,                   sqpopt_linesearch_armijo,   sqpopt_qp_composite)
    call run_hs71('augmented Lagrangian',   sqpopt_merit_augmented_lagrangian, sqpopt_linesearch_armijo,   sqpopt_qp_composite)
    call run_hs71('watchdog',               sqpopt_merit_l1,                   sqpopt_linesearch_watchdog, sqpopt_qp_composite)
    call run_hs71('dense QP',               sqpopt_merit_l1,                   sqpopt_linesearch_armijo,   sqpopt_qp_dense)
    call run_hs71('reduced-Hessian QP',     sqpopt_merit_l1,                   sqpopt_linesearch_armijo,   sqpopt_qp_reduced_hessian)

    contains

    subroutine run_hs71(label, merit_mode, linesearch_mode, qp_mode)

    character(len=*), intent(in) :: label
    integer,           intent(in) :: merit_mode
    integer,           intent(in) :: linesearch_mode
    integer,           intent(in) :: qp_mode

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(4), xsol(4), lam(2)
    real(wp), parameter :: xexpect(4) = [1.0_wp, 4.7429994_wp, 3.8211500_wp, 1.3794083_wp]
    real(wp), parameter :: fexpect = 17.0140173_wp
    real(wp) :: fsol
    integer :: istat

    i_obj  = 0
    i_grad = 0
    i_cons = 0
    i_jac  = 0

    ! equality constraint (index 1) first, then the inequality (index 2):
    call problem%set_problem_size(n=4, m_eq=1, m_ineq=1)
    call problem%set_bounds(x_lb=[1.0_wp,1.0_wp,1.0_wp,1.0_wp], x_ub=[5.0_wp,5.0_wp,5.0_wp,5.0_wp], &
                             c_lb=[40.0_wp, 25.0_wp], c_ub=[40.0_wp, big])
    call problem%set_jacobian_sparsity(nnz=8, irow=[1,1,1,1,2,2,2,2], icol=[1,2,3,4,1,2,3,4])
    call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)

    options%max_iter        = 3000
    options%merit_mode      = merit_mode
    options%linesearch_mode = linesearch_mode
    options%qp_solver_mode  = qp_mode
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
    print '(A,I0)',     'i_obj   = ', i_obj
    print '(A,I0)',     'i_grad  = ', i_grad
    print '(A,I0)',     'i_cons  = ', i_cons
    print '(A,I0)',     'i_jac   = ', i_jac

    if (qp_mode == sqpopt_qp_dense .or. qp_mode == sqpopt_qp_reduced_hessian) then
        ! a real QP solve should converge tightly (see DENSE_QP_PLAN.md/REDUCED_HESSIAN_QP_PLAN.md):
        if (istat /= sqpopt_success) error stop 'test_hs71 FAILED: '//label//' did not reach sqpopt_success'
        if (maxval(abs(xsol-xexpect)) > 1.0e-4_wp) error stop 'test_hs71 FAILED: '//label//' wrong solution'
    else
        ! v1's composite-step QP only achieves loose convergence on this
        ! problem (see the note above), so a generous tolerance is used here:
        if (maxval(abs(xsol-xexpect)) > 0.5_wp) error stop 'test_hs71 FAILED: '//label
    end if
    print '(A)', 'test_hs71 ['//trim(label)//'] PASSED'
    print '(A)', ''

    end subroutine run_hs71

    subroutine obj(x, f)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),                intent(out) :: f
    f = x(1)*x(4)*(x(1)+x(2)+x(3)) + x(3)
    i_obj = i_obj + 1
    end subroutine obj

    subroutine grad(x, g)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    g(1) = x(4)*(2.0_wp*x(1)+x(2)+x(3))
    g(2) = x(1)*x(4)
    g(3) = x(1)*x(4) + 1.0_wp
    g(4) = x(1)*(x(1)+x(2)+x(3))
    i_grad = i_grad + 1
    end subroutine grad

    subroutine cons(x, c)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    c(1) = x(1)**2 + x(2)**2 + x(3)**2 + x(4)**2
    c(2) = x(1)*x(2)*x(3)*x(4)
    i_cons = i_cons + 1
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
    i_jac = i_jac + 1
    end subroutine jacv

end program test_hs71
