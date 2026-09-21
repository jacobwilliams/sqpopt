program test_basic

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type, sqpopt_hessian_sr1
    use sqpopt_linesearch_module, only: sqpopt_linesearch_exact, sqpopt_linesearch_watchdog, sqpopt_merit_augmented_lagrangian, &
                                         sqpopt_linesearch_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_dense, sqpopt_qp_reduced_hessian, sqpopt_qp_solver_type, sqpopt_bounds_vector
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
    call test_vector_bound_enforcement()
    call test_major_step_limit()
    call test_print_level_and_stalled_progress()

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

    call problem%set_problem_size(n=2, m_eq=1, m_ineq=0)
    call problem%set_bounds(x_lb=[-big,-big], x_ub=[big,big], c_lb=[4.0_wp], c_ub=[4.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

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

    subroutine obj1(x, f)
    real(wp), dimension(:), intent(in)  :: x
    real(wp),                intent(out) :: f
    f = (x(1)-2.0_wp)**2 + (x(2)-3.0_wp)**2
    end subroutine obj1

    subroutine grad1(x, g)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: g
    g(1) = 2.0_wp*(x(1)-2.0_wp)
    g(2) = 2.0_wp*(x(2)-3.0_wp)
    end subroutine grad1

    !> no constraints (used for the bounds-only test, where m=0)
    subroutine cons0(x, c)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    associate(unused => x); end associate
    associate(unused => c); end associate
    end subroutine cons0

    !> no constraints (used for the bounds-only test, where m=0)
    subroutine jacv0(x, jac_val)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: jac_val
    associate(unused => x); end associate
    associate(unused => jac_val); end associate
    end subroutine jacv0

    subroutine cons1(x, c)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: c
    c(1) = x(1) + x(2)
    end subroutine cons1

    subroutine jacv1(x, jac_val)
    real(wp), dimension(:), intent(in)  :: x
    real(wp), dimension(:), intent(out) :: jac_val
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

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

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

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=0)
    call problem%set_bounds(x_lb=[-big,-big], x_ub=[1.0_wp,big], c_lb=[real(wp) ::], c_ub=[real(wp) ::])
    call problem%set_jacobian_sparsity(nnz=0, irow=no_rows, icol=no_rows)
    call problem%set_functions(f=obj1, g=grad1, c=cons0, jac=jacv0)

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

    call problem%set_problem_size(n=2, m_eq=1, m_ineq=0)
    call problem%set_bounds(x_lb=[-big,-big], x_ub=[big,big], c_lb=[4.0_wp], c_ub=[4.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

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

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

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

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

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

    !> same problem as `test_inequality_constrained`, but using the smooth
    !! augmented Lagrangian merit function instead of the default l1 one.
    subroutine test_augmented_lagrangian_merit()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

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
    !! active-set QP solver instead of the default v1 composite step.
    subroutine test_dense_qp_mode()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

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

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

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

    !> same problem as `test_inequality_constrained`, but using the
    !! composite-step QP's `sqpopt_bounds_vector` bound-enforcement mode
    !! (uniformly rescale the whole step, instead of clipping only the
    !! violating components) instead of the default `sqpopt_bounds_scalar`.
    subroutine test_vector_bound_enforcement()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_qp_solver_type) :: qp_solver
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

    options%max_iter = 100
    qp_solver%bound_enforcement = sqpopt_bounds_vector
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options, qp_solver=qp_solver)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_vector_bound_enforcement: x      = ', xsol
    print '(A,2F12.6)', 'test_vector_bound_enforcement: x_true = ', xexpect
    print '(A,I0)',     'test_vector_bound_enforcement: istat  = ', istat

    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_vector_bound_enforcement FAILED'
    print *, 'test_vector_bound_enforcement PASSED'

    end subroutine test_vector_bound_enforcement

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

    call problem%set_problem_size(n=2, m_eq=0, m_ineq=1)
    call problem%set_bounds(x_lb=[0.0_wp,0.0_wp], x_ub=[big,big], c_lb=[-big], c_ub=[3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

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
    !! and a `ktol` tight enough that it's very unlikely to ever be met exactly,
    !! relying instead on the `ftol`/`xtol` stalled-progress criterion to reach
    !! `sqpopt_success` once the iterates stop moving.
    subroutine test_print_level_and_stalled_progress()

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp) :: x0(2), xsol(2), lam(1)
    real(wp), parameter :: xexpect(2) = [1.5_wp, 2.5_wp]
    integer :: istat

    call problem%set_problem_size(n=2, m_eq=1, m_ineq=0)
    call problem%set_bounds(x_lb=[-big,-big], x_ub=[big,big], c_lb=[4.0_wp], c_ub=[4.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(f=obj1, g=grad1, c=cons1, jac=jacv1)

    options%max_iter    = 100
    options%print_level = 1
    options%ktol         = 1.0e-15_wp !! essentially unreachable
    options%ftol         = 1.0e-8_wp
    options%xtol         = 1.0e-8_wp
    x0 = [0.0_wp, 0.0_wp]

    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)

    print '(A,2F12.6)', 'test_print_level_and_stalled_progress: x      = ', xsol
    print '(A,2F12.6)', 'test_print_level_and_stalled_progress: x_true = ', xexpect
    print '(A,I0)',     'test_print_level_and_stalled_progress: istat  = ', istat

    if (istat /= 0) error stop 'test_print_level_and_stalled_progress FAILED: istat'
    if (maxval(abs(xsol-xexpect)) > 1.0e-3_wp) error stop 'test_print_level_and_stalled_progress FAILED: wrong x'
    print *, 'test_print_level_and_stalled_progress PASSED'

    end subroutine test_print_level_and_stalled_progress

end program test_basic
