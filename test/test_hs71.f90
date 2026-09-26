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
    !! (a corner point); every run below is required to reach
    !! `sqpopt_success` tightly.
    !!
    !! Run nine ways, to compare the available merit functions, penalty
    !! updates, line searches, and QP solvers on this problem: the \(\ell_1\)
    !! merit with the Armijo line search and `sqpopt_qp_auto` (which picks
    !! the dense QP solver for a problem this small), the model-based penalty
    !! updates (`sqpopt_penalty_model`: Byrd-Nocedal for \(\ell_1\), and
    !! Gill-Murray-Saunders-Wright for the augmented Lagrangian),
    !! `sqpopt_merit_augmented_lagrangian` (see PLAN.md section 6.1),
    !! `sqpopt_linesearch_watchdog` (Powell's VF13 watchdog technique, see
    !! PLAN.md section 6.3), `sqpopt_qp_dense` and
    !! `sqpopt_qp_reduced_hessian` explicitly (see `DENSE_QP_PLAN.md` and
    !! `REDUCED_HESSIAN_QP_PLAN.md`), the reduced-Hessian solver with tuned
    !! `LSQR` tolerances (`lsqr_atol`/`lsqr_btol`, exposed as user-settable
    !! fields on `sqpopt_reduced_hessian_qp_type`), and the filter line
    !! search with the dense QP solver. The function-call counts
    !! (`i_obj`/`i_grad`/`i_cons`/`i_jac` below) are printed for comparison.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_linesearch_module, only: sqpopt_merit_l1, sqpopt_merit_augmented_lagrangian, &
                                        sqpopt_linesearch_armijo, sqpopt_linesearch_watchdog, sqpopt_linesearch_filter, &
                                        sqpopt_penalty_multipliers, sqpopt_penalty_model
    use sqpopt_qp_solver_module, only: sqpopt_qp_auto, sqpopt_qp_dense, sqpopt_qp_reduced_hessian, sqpopt_qp_solver_type
    use sqpopt_types_module,   only: sqpopt_success, sqpopt_results_type
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    integer :: i_obj, i_grad, i_cons, i_jac
    logical :: outside_bounds !! set if any function is evaluated outside the variable bounds `1 <= x <= 5`

    write(*,*) '----------------------------'
    write(*,*) 'test_hs71'
    write(*,*) '----------------------------'

    call run_hs71('l1 + Armijo',            sqpopt_merit_l1,                   sqpopt_linesearch_armijo,   sqpopt_qp_auto)
    call run_hs71('augmented Lagrangian',   sqpopt_merit_augmented_lagrangian, sqpopt_linesearch_armijo,   sqpopt_qp_auto)
    call run_hs71('l1, Byrd-Nocedal penalty', sqpopt_merit_l1,                 sqpopt_linesearch_armijo,   sqpopt_qp_auto, &
                  penalty_update=sqpopt_penalty_model)
    call run_hs71('AL, GMSW penalty',       sqpopt_merit_augmented_lagrangian, sqpopt_linesearch_armijo,   sqpopt_qp_auto, &
                  penalty_update=sqpopt_penalty_model)
    call run_hs71('watchdog',               sqpopt_merit_l1,                   sqpopt_linesearch_watchdog, sqpopt_qp_auto)
    call run_hs71('dense QP',               sqpopt_merit_l1,                   sqpopt_linesearch_armijo,   sqpopt_qp_dense)
    call run_hs71('reduced-Hessian QP',     sqpopt_merit_l1,                   sqpopt_linesearch_armijo,   sqpopt_qp_reduced_hessian)
    ! the reduced-Hessian solver's LSQR tolerances are user-tunable (see
    ! PLAN.md section 6.2); kept as an example of setting them (loosening
    ! them trades QP accuracy for speed):
    call run_hs71('rh, tuned LSQR (atol=btol=5e-10)', sqpopt_merit_l1, sqpopt_linesearch_armijo, sqpopt_qp_reduced_hessian, &
                  lsqr_atol=5.0e-10_wp, lsqr_btol=5.0e-10_wp)
    ! filter method (Fletcher & Leyffer, no merit function/penalty parameter
    ! at all) paired with a real active-set QP solve -- also reaches
    ! sqpopt_success tightly, confirming the filter line search is a viable
    ! drop-in alternative globalization strategy alongside the merit-based ones:
    call run_hs71('filter + dense QP',      sqpopt_merit_l1,                   sqpopt_linesearch_filter,  sqpopt_qp_dense)

    contains

    subroutine run_hs71(label, merit_mode, linesearch_mode, qp_mode, lsqr_atol, lsqr_btol, lsqr_itnlim, penalty_update)

    character(len=*), intent(in) :: label
    integer,           intent(in) :: merit_mode
    integer,           intent(in) :: linesearch_mode
    integer,           intent(in) :: qp_mode
    real(wp), intent(in), optional :: lsqr_atol, lsqr_btol !! LSQR tuning (sqpopt_qp_reduced_hessian mode only)
    integer,  intent(in), optional :: lsqr_itnlim          !! LSQR tuning (sqpopt_qp_reduced_hessian mode only)
    integer,  intent(in), optional :: penalty_update       !! penalty update rule (default `sqpopt_penalty_multipliers`)

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_qp_solver_type) :: qp_solver
    real(wp) :: x0(4), xsol(4), lam(2)
    real(wp), parameter :: xexpect(4) = [1.0_wp, 4.7429994_wp, 3.8211500_wp, 1.3794083_wp]
    real(wp), parameter :: fexpect = 17.0140173_wp
    real(wp) :: fsol
    integer :: istat

    i_obj  = 0
    i_grad = 0
    i_cons = 0
    i_jac  = 0
    outside_bounds = .false.

    ! equality constraint (index 1) first, then the inequality (index 2):
    call problem%set_problem_size(n=4, m_eq=1, m_ineq=1)
    call problem%set_bounds(x_lb=[1.0_wp,1.0_wp,1.0_wp,1.0_wp], x_ub=[5.0_wp,5.0_wp,5.0_wp,5.0_wp], &
                             c_lb=[40.0_wp, 25.0_wp], c_ub=[40.0_wp, big])
    call problem%set_jacobian_sparsity(nnz=8, irow=[1,1,1,1,2,2,2,2], icol=[1,2,3,4,1,2,3,4])
    call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv)

    options%max_iter        = 3000
    options%merit_mode      = merit_mode
    options%linesearch_mode = linesearch_mode
    options%penalty_update  = sqpopt_penalty_multipliers
    if (present(penalty_update)) options%penalty_update = penalty_update
    options%qp_solver_mode  = qp_mode
    qp_solver%mode          = qp_mode
    if (present(lsqr_atol))   qp_solver%sparse_qp%lsqr_atol   = lsqr_atol
    if (present(lsqr_btol))   qp_solver%sparse_qp%lsqr_btol   = lsqr_btol
    if (present(lsqr_itnlim)) qp_solver%sparse_qp%lsqr_itnlim = lsqr_itnlim
    x0 = [1.0_wp, 5.0_wp, 5.0_wp, 1.0_wp]

    call solver%initialize(problem=problem, options=options, qp_solver=qp_solver)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)
    block
        type(sqpopt_results_type) :: results
        call solver%get_results(results)
        fsol = results%f
    end block

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

    ! the solver must never evaluate the functions outside the variable bounds:
    if (outside_bounds) error stop 'test_hs71 FAILED: '//label//' evaluated a function outside the variable bounds'

    ! every variant uses a real QP solve, so must converge tightly:
    if (istat /= sqpopt_success) error stop 'test_hs71 FAILED: '//label//' did not reach sqpopt_success'
    if (maxval(abs(xsol-xexpect)) > 1.0e-4_wp) error stop 'test_hs71 FAILED: '//label//' wrong solution'
    print '(A)', 'test_hs71 ['//trim(label)//'] PASSED'
    print '(A)', ''

    end subroutine run_hs71

    subroutine obj(x, f, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),                intent(out) :: f
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    f = x(1)*x(4)*(x(1)+x(2)+x(3)) + x(3)
    i_obj = i_obj + 1
    if (any(x < 1.0_wp) .or. any(x > 5.0_wp)) outside_bounds = .true.
    end subroutine obj

    subroutine grad(x, g, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    g(1) = x(4)*(2.0_wp*x(1)+x(2)+x(3))
    g(2) = x(1)*x(4)
    g(3) = x(1)*x(4) + 1.0_wp
    g(4) = x(1)*(x(1)+x(2)+x(3))
    i_grad = i_grad + 1
    if (any(x < 1.0_wp) .or. any(x > 5.0_wp)) outside_bounds = .true.
    end subroutine grad

    subroutine cons(x, c, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    c(1) = x(1)**2 + x(2)**2 + x(3)**2 + x(4)**2
    c(2) = x(1)*x(2)*x(3)*x(4)
    i_cons = i_cons + 1
    if (any(x < 1.0_wp) .or. any(x > 5.0_wp)) outside_bounds = .true.
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: jac_val
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
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
    if (any(x < 1.0_wp) .or. any(x > 5.0_wp)) outside_bounds = .true.
    end subroutine jacv

end program test_hs71
