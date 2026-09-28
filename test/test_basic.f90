program test_basic

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_hessian_module, only: sqpopt_hessian_sr1
    use sqpopt_linesearch_module, only: sqpopt_linesearch_exact, sqpopt_linesearch_watchdog, sqpopt_merit_augmented_lagrangian, &
                                         sqpopt_linesearch_type, sqpopt_linesearch_filter
    use sqpopt_qp_solver_module, only: sqpopt_qp_dense, sqpopt_qp_reduced_hessian
    use sqpopt_types_module,   only: sqpopt_success, sqpopt_stalled
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    write(*,*) '----------------------------'
    write(*,*) 'test_basic'
    write(*,*) '----------------------------'

    call test_equality_constrained()
    call test_inequality_constrained()
    call test_bounds_only()
    call test_sr1_hessian_mode()
    call test_exact_linesearch_mode()
    call test_watchdog_linesearch_mode()
    call test_augmented_lagrangian_merit()
    call test_dense_qp_mode()
    call test_reduced_hessian_qp_mode()
    call test_major_step_limit()
    call test_print_level_and_stalled_progress()
    call test_filter_linesearch_mode()
    call test_filter_linesearch_mode_equality()
    call test_trust_region_mode()
    call test_trust_region_filter_mode()

    contains

    !> minimize (x1-2)^2 + (x2-3)^2  s.t.  x1+x2 = 4
    !! known solution: x* = (1.5, 2.5), f* = 0.5
    subroutine test_equality_constrained()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.5_wp, 2.5_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-big,-big], x_ub=[big,big], c_lb=[4.0_wp], c_ub=[4.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter = 100
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_equality_constrained: x      = ', xsol
    print '(A,2F12.6)', 'test_equality_constrained: x_true = ', xexpect
    print '(A,I0)',     'test_equality_constrained: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_equality_constrained FAILED'
    print *, 'test_equality_constrained PASSED'

    end subroutine test_equality_constrained

    subroutine obj1(x, f, status, data)
    !! the objective
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),                intent(out)  :: f      !! objective value at `x`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = (x(1)-2.0_wp)**2 + (x(2)-3.0_wp)**2
    end subroutine obj1

    subroutine grad1(x, g, status, data)
    !! the objective's gradient
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    g(1) = 2.0_wp*(x(1)-2.0_wp)
    g(2) = 2.0_wp*(x(2)-3.0_wp)
    end subroutine grad1

    !> no constraints (used for the bounds-only test, where m=0)
    subroutine cons0(x, c, status, data)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    associate(unused => x); end associate
    associate(unused => c); end associate
    end subroutine cons0

    !> no constraints (used for the bounds-only test, where m=0)
    subroutine jacv0(x, jac_val, status, data)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    associate(unused => x); end associate
    associate(unused => jac_val); end associate
    end subroutine jacv0

    subroutine cons1(x, c, status, data)
    !! the constraints
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    c(1) = x(1) + x(2)
    end subroutine cons1

    subroutine jacv1(x, jac_val, status, data)
    !! the nonzero values of the constraint Jacobian (in the sparsity pattern's order)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    jac_val = [1.0_wp, 1.0_wp]
    end subroutine jacv1

    !> minimize (x1-2)^2 + (x2-3)^2  s.t.  x1+x2 <= 3,  x1,x2 >= 0
    !! known solution: x* = (1, 2), f* = 2  (the inequality is active)
    subroutine test_inequality_constrained()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter = 100
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_inequality_constrained: x      = ', xsol
    print '(A,2F12.6)', 'test_inequality_constrained: x_true = ', xexpect
    print '(A,I0)',     'test_inequality_constrained: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_inequality_constrained FAILED'
    print *, 'test_inequality_constrained PASSED'

    end subroutine test_inequality_constrained

    !> minimize (x1-2)^2 + (x2-3)^2  s.t.  x1 <= 1  (no general constraints)
    !! known solution: x* = (1, 3), f* = 1  (only the bound on x1 is active)
    subroutine test_bounds_only()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(0)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 3.0_wp]
    integer :: istat
    integer, dimension(0) :: no_rows

    call problem%set_problem_size(n=2, m=0)
    call problem%set_bounds(x_lb=[-big,-big], x_ub=[1.0_wp,big], c_lb=[real(wp) ::], c_ub=[real(wp) ::])
    call problem%set_jacobian_sparsity(nnz=0, irow=no_rows, icol=no_rows)
    call problem%set_functions(fc=fc_obj1_cons0, gjac=gjac_grad1_jacv0)

    options%max_iter = 100
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_bounds_only: x      = ', xsol
    print '(A,2F12.6)', 'test_bounds_only: x_true = ', xexpect
    print '(A,I0)',     'test_bounds_only: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_bounds_only FAILED'
    print *, 'test_bounds_only PASSED'

    end subroutine test_bounds_only

    !> same problem as `test_equality_constrained`, but using the limited-memory
    !! SR1 Hessian approximation instead of the default BFGS.
    subroutine test_sr1_hessian_mode()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.5_wp, 2.5_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-big,-big], x_ub=[big,big], c_lb=[4.0_wp], c_ub=[4.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter     = 100
    options%hessian_mode = sqpopt_hessian_sr1
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_sr1_hessian_mode: x      = ', xsol
    print '(A,2F12.6)', 'test_sr1_hessian_mode: x_true = ', xexpect
    print '(A,I0)',     'test_sr1_hessian_mode: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_sr1_hessian_mode FAILED'
    print *, 'test_sr1_hessian_mode PASSED'

    end subroutine test_sr1_hessian_mode

    !> same problem as `test_inequality_constrained`, but using the
    !! exact (fmin-based) line search mode instead of the default Armijo.
    subroutine test_exact_linesearch_mode()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter       = 100
    options%linesearch_mode = sqpopt_linesearch_exact
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_exact_linesearch_mode: x      = ', xsol
    print '(A,2F12.6)', 'test_exact_linesearch_mode: x_true = ', xexpect
    print '(A,I0)',     'test_exact_linesearch_mode: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_exact_linesearch_mode FAILED'
    print *, 'test_exact_linesearch_mode PASSED'

    end subroutine test_exact_linesearch_mode

    !> same problem as `test_inequality_constrained`, but using the
    !! watchdog line search mode instead of the default Armijo.
    subroutine test_watchdog_linesearch_mode()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter        = 100
    options%linesearch_mode = sqpopt_linesearch_watchdog
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_watchdog_linesearch_mode: x      = ', xsol
    print '(A,2F12.6)', 'test_watchdog_linesearch_mode: x_true = ', xexpect
    print '(A,I0)',     'test_watchdog_linesearch_mode: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_watchdog_linesearch_mode FAILED'
    print *, 'test_watchdog_linesearch_mode PASSED'

    end subroutine test_watchdog_linesearch_mode

    !> same problem as `test_inequality_constrained`, but using Fletcher &
    !! Leyffer's filter method (no merit function/penalty parameter).
    subroutine test_filter_linesearch_mode()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter        = 100
    options%linesearch_mode = sqpopt_linesearch_filter
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_filter_linesearch_mode: x      = ', xsol
    print '(A,2F12.6)', 'test_filter_linesearch_mode: x_true = ', xexpect
    print '(A,I0)',     'test_filter_linesearch_mode: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_filter_linesearch_mode FAILED'
    print *, 'test_filter_linesearch_mode PASSED'

    end subroutine test_filter_linesearch_mode

    !> same problem as `test_equality_constrained`, with the filter line
    !! search: exercises the "always feasible" fallback (h stays near zero
    !! throughout, since `x0` starts feasible and the problem has only one
    !! equality constraint), which needs the plain-descent-in-f safeguard
    !! described in the module docs (otherwise the filter test alone would
    !! accept any h<=ctol point regardless of f).
    subroutine test_filter_linesearch_mode_equality()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.5_wp, 2.5_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-big,-big], x_ub=[big,big], c_lb=[4.0_wp], c_ub=[4.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter        = 100
    options%linesearch_mode = sqpopt_linesearch_filter
    x0 = [4.0_wp, 0.0_wp] !! feasible starting point (satisfies x1+x2=4)

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_filter_linesearch_mode_equality: x      = ', xsol
    print '(A,2F12.6)', 'test_filter_linesearch_mode_equality: x_true = ', xexpect
    print '(A,I0)',     'test_filter_linesearch_mode_equality: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_filter_linesearch_mode_equality FAILED'
    print *, 'test_filter_linesearch_mode_equality PASSED'

    end subroutine test_filter_linesearch_mode_equality

    !> same problem as `test_inequality_constrained`, with trust-region
    !! globalization enabled instead of a line search (classical merit-
    !! ratio acceptance, since `linesearch%mode` defaults to
    !! `sqpopt_linesearch_armijo`, not `sqpopt_linesearch_filter`).
    subroutine test_trust_region_mode()

    type(sqpopt_type)             :: solver
    type(sqpopt_problem_type)     :: problem
    type(sqpopt_options_type)     :: options
    type(sqpopt_trust_region_type) :: trust_region
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter        = 100
    trust_region%enabled    = .true.
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options, trust_region=trust_region)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_trust_region_mode: x      = ', xsol
    print '(A,2F12.6)', 'test_trust_region_mode: x_true = ', xexpect
    print '(A,I0)',     'test_trust_region_mode: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_trust_region_mode FAILED'
    print *, 'test_trust_region_mode PASSED'

    end subroutine test_trust_region_mode

    !> same problem, with trust-region globalization AND filter-based
    !! acceptance (`linesearch%mode = sqpopt_linesearch_filter`) -- the
    !! *literal* Fletcher & Leyffer filter-SQP combination (see
    !! `plan/TRUST_REGION_PLAN.md` and `sqpopt_trust_region_module`'s docs).
    subroutine test_trust_region_filter_mode()

    type(sqpopt_type)             :: solver
    type(sqpopt_problem_type)     :: problem
    type(sqpopt_options_type)     :: options
    type(sqpopt_trust_region_type) :: trust_region
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter        = 100
    options%linesearch_mode = sqpopt_linesearch_filter
    trust_region%enabled    = .true.
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options, trust_region=trust_region)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_trust_region_filter_mode: x      = ', xsol
    print '(A,2F12.6)', 'test_trust_region_filter_mode: x_true = ', xexpect
    print '(A,I0)',     'test_trust_region_filter_mode: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_trust_region_filter_mode FAILED'
    print *, 'test_trust_region_filter_mode PASSED'

    end subroutine test_trust_region_filter_mode

    !> same problem as `test_inequality_constrained`, but using the smooth
    !! augmented Lagrangian merit function instead of the default l1 one.
    subroutine test_augmented_lagrangian_merit()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter  = 100
    options%merit_mode = sqpopt_merit_augmented_lagrangian
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_augmented_lagrangian_merit: x      = ', xsol
    print '(A,2F12.6)', 'test_augmented_lagrangian_merit: x_true = ', xexpect
    print '(A,I0)',     'test_augmented_lagrangian_merit: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_augmented_lagrangian_merit FAILED'
    print *, 'test_augmented_lagrangian_merit PASSED'

    end subroutine test_augmented_lagrangian_merit

    !> same problem as `test_inequality_constrained`, but using the dense
    !! dense active-set QP solver explicitly (instead of `sqpopt_qp_auto`).
    subroutine test_dense_qp_mode()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter       = 100
    options%qp_solver_mode = sqpopt_qp_dense
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_dense_qp_mode: x      = ', xsol
    print '(A,2F12.6)', 'test_dense_qp_mode: x_true = ', xexpect
    print '(A,I0)',     'test_dense_qp_mode: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_dense_qp_mode FAILED'
    print *, 'test_dense_qp_mode PASSED'

    end subroutine test_dense_qp_mode

    !> same problem as `test_inequality_constrained`, but using the sparse
    !! (projected-CG) reduced-Hessian active-set QP solver.
    subroutine test_reduced_hessian_qp_mode()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter       = 100
    options%qp_solver_mode = sqpopt_qp_reduced_hessian
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_reduced_hessian_qp_mode: x      = ', xsol
    print '(A,2F12.6)', 'test_reduced_hessian_qp_mode: x_true = ', xexpect
    print '(A,I0)',     'test_reduced_hessian_qp_mode: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_reduced_hessian_qp_mode FAILED'
    print *, 'test_reduced_hessian_qp_mode PASSED'

    end subroutine test_reduced_hessian_qp_mode

    !> same problem as `test_inequality_constrained`, but with a tight
    !! `major_step_limit` (SNOPT-inspired option) that forces the very
    !! first major iteration's step to be shrunk well below what the
    !! (uncapped) QP step and line search would otherwise take; confirms
    !! the option doesn't prevent eventual convergence, just slows it down.
    subroutine test_major_step_limit()

    type(sqpopt_type)            :: solver
    type(sqpopt_problem_type)    :: problem
    type(sqpopt_options_type)    :: options
    type(sqpopt_linesearch_type) :: linesearch
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter = 200
    linesearch%major_step_limit = 0.05_wp
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options, linesearch=linesearch)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_major_step_limit: x      = ', xsol
    print '(A,2F12.6)', 'test_major_step_limit: x_true = ', xexpect
    print '(A,I0)',     'test_major_step_limit: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_major_step_limit FAILED'
    print *, 'test_major_step_limit PASSED'

    end subroutine test_major_step_limit

    !> same problem as `test_equality_constrained`, but with `print_level=1`
    !! (confirms the per-iteration diagnostic printing doesn't break anything)
    !! and a very tight `ktol`. On this quadratic
    !! problem the KKT residual can reach exactly zero, so the run may end
    !! either on the KKT test (`sqpopt_success`) or on the `ftol`/`xtol`
    !! stalled-progress test (`sqpopt_stalled`); both are accepted here (the
    !! stalled-progress test itself is unit-tested in `test_convergence`).
    subroutine test_print_level_and_stalled_progress()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.5_wp, 2.5_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-big,-big], x_ub=[big,big], c_lb=[4.0_wp], c_ub=[4.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj1_cons1, gjac=gjac_grad1_jacv1)

    options%max_iter    = 100
    options%print_level = 1
    options%ktol         = 1.0e-15_wp
    options%ftol         = 1.0e-8_wp
    options%xtol         = 1.0e-8_wp
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_print_level_and_stalled_progress: x      = ', xsol
    print '(A,2F12.6)', 'test_print_level_and_stalled_progress: x_true = ', xexpect
    print '(A,I0)',     'test_print_level_and_stalled_progress: istat  = ', istat

    if (istat /= sqpopt_success .and. istat /= sqpopt_stalled) error stop 'test_print_level_and_stalled_progress FAILED: istat'
    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_print_level_and_stalled_progress FAILED: wrong x'
    print *, 'test_print_level_and_stalled_progress PASSED'

    end subroutine test_print_level_and_stalled_progress

    subroutine fc_obj1_cons0(x, f, c, status, data)
    !! `fc` for `set_functions`: the objective (`obj1`) and the constraints (`cons0`)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    call obj1(x, f, status, data)
    if (status == 0) call cons0(x, c, status, data)
    end subroutine fc_obj1_cons0

    subroutine fc_obj1_cons1(x, f, c, status, data)
    !! `fc` for `set_functions`: the objective (`obj1`) and the constraints (`cons1`)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    call obj1(x, f, status, data)
    if (status == 0) call cons1(x, c, status, data)
    end subroutine fc_obj1_cons1

    subroutine gjac_grad1_jacv0(x, g, jac_val, accuracy, status, data)
    !! `gjac` for `set_functions`: the gradient (`grad1`) and the Jacobian values (`jacv0`)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g       !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    associate(unused => accuracy); end associate
    call grad1(x, g, status, data)
    if (status == 0) call jacv0(x, jac_val, status, data)
    end subroutine gjac_grad1_jacv0

    subroutine gjac_grad1_jacv1(x, g, jac_val, accuracy, status, data)
    !! `gjac` for `set_functions`: the gradient (`grad1`) and the Jacobian values (`jacv1`)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g       !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    associate(unused => accuracy); end associate
    call grad1(x, g, status, data)
    if (status == 0) call jacv1(x, jac_val, status, data)
    end subroutine gjac_grad1_jacv1


end program test_basic
