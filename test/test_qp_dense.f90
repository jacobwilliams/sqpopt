program test_qp_dense

    !! Standalone unit tests for [[sqpopt_qp_dense_module]] (the dense
    !! active-set QP solver, see `DENSE_QP_PLAN.md`), called *directly*
    !! (bypassing the outer NLP loop entirely) on small, hand-verifiable
    !! QPs of the form `min 0.5*p^T*H*p + g^T*p` s.t. linear
    !! constraints/bounds, with `H = I` (a freshly-initialized
    !! [[sqpopt_hessian_type]] gives exactly this, since `gamma=1` and no
    !! `(s,y)` history means `hv_product(v) = v`), so each solution can be
    !! checked directly against the QP's own KKT/Lagrange conditions.

    use sqpopt_qp_dense_module, only: sqpopt_dense_qp_type
    use sqpopt_hessian_module,  only: sqpopt_hessian_type
    use sqpopt_types_module,    only: sqpopt_sparse_matrix, sqpopt_success
    use sqpopt_kinds,           only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp

    write(*,*) '----------------------------'
    write(*,*) 'test_qp_dense'
    write(*,*) '----------------------------'

    call test_bounds_only()
    call test_equality()
    call test_inequality()
    call test_two_active_at_once()

    contains

    !> minimize 0.5*(p1^2+p2^2) - 2*p1 - 3*p2  s.t.  p1 <= 1  (no general constraints)
    !! unconstrained minimizer is (2,3); the bound on p1 is active.
    !! KKT: p1=1 (bound), p2=3 (free, stationary: p2-3=0)
    subroutine test_bounds_only()

    type(sqpopt_dense_qp_type) :: qp
    type(sqpopt_hessian_type)  :: hessian
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: g(2), x(2), c(0), x_lb(2), x_ub(2), c_lb(0), c_ub(0)
    real(wp) :: p(2), lambda(0)
    real(wp), parameter :: pexpect(2) = [1.0_wp, 3.0_wp]
    integer :: istat

    call hessian%initialize(n=2, max_history=5)
    jac%nrows = 0; jac%ncols = 2; jac%nnz = 0
    allocate(jac%irow(0), jac%icol(0), jac%val(0))
    g = [-2.0_wp, -3.0_wp]
    x = [0.0_wp, 0.0_wp]
    x_lb = [-big, -big]
    x_ub = [1.0_wp, big]

    call qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    print '(A,2F12.7)', 'test_bounds_only: p      = ', p
    print '(A,2F12.7)', 'test_bounds_only: p_true = ', pexpect
    print '(A,I0)',     'test_bounds_only: istat  = ', istat
    if (istat /= sqpopt_success) error stop 'test_bounds_only FAILED: istat'
    if (maxval(abs(p-pexpect)) > 1.0e-6_wp) error stop 'test_bounds_only FAILED: wrong p'
    print *, 'test_bounds_only PASSED'

    end subroutine test_bounds_only

    !> minimize 0.5*(p1^2+p2^2) - 2*p1 - 3*p2  s.t.  p1+p2 = 1
    !! KKT: p1=2+lam, p2=3+lam, p1+p2=1 => lam=-2 => p=(0,1)
    subroutine test_equality()

    type(sqpopt_dense_qp_type) :: qp
    type(sqpopt_hessian_type)  :: hessian
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: g(2), x(2), c(1), x_lb(2), x_ub(2), c_lb(1), c_ub(1)
    real(wp) :: p(2), lambda(1)
    real(wp), parameter :: pexpect(2) = [0.0_wp, 1.0_wp]
    integer :: istat

    call hessian%initialize(n=2, max_history=5)
    jac%nrows = 1; jac%ncols = 2; jac%nnz = 2
    jac%irow = [1,1]; jac%icol = [1,2]; jac%val = [1.0_wp, 1.0_wp]
    g = [-2.0_wp, -3.0_wp]
    x = [0.0_wp, 0.0_wp]
    c = [0.0_wp]
    x_lb = [-big, -big]
    x_ub = [big, big]
    c_lb = [1.0_wp]
    c_ub = [1.0_wp]

    call qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    print '(A,2F12.7)', 'test_equality: p      = ', p
    print '(A,2F12.7)', 'test_equality: p_true = ', pexpect
    print '(A,F12.7)',  'test_equality: lambda = ', lambda(1)
    print '(A,I0)',     'test_equality: istat  = ', istat
    if (istat /= sqpopt_success) error stop 'test_equality FAILED: istat'
    if (maxval(abs(p-pexpect)) > 1.0e-6_wp) error stop 'test_equality FAILED: wrong p'
    print *, 'test_equality PASSED'

    end subroutine test_equality

    !> minimize 0.5*(p1^2+p2^2) - 2*p1 - 3*p2  s.t.  p1+p2 <= 3,  p1,p2 >= 0
    !! (the linearization, at x=(0,0), of test_basic's test_inequality_constrained)
    !! KKT: p1=2+lam, p2=3+lam, p1+p2=3 => lam=-1 => p=(1,2)
    subroutine test_inequality()

    type(sqpopt_dense_qp_type) :: qp
    type(sqpopt_hessian_type)  :: hessian
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: g(2), x(2), c(1), x_lb(2), x_ub(2), c_lb(1), c_ub(1)
    real(wp) :: p(2), lambda(1)
    real(wp), parameter :: pexpect(2) = [1.0_wp, 2.0_wp]
    integer :: istat

    call hessian%initialize(n=2, max_history=5)
    jac%nrows = 1; jac%ncols = 2; jac%nnz = 2
    jac%irow = [1,1]; jac%icol = [1,2]; jac%val = [1.0_wp, 1.0_wp]
    g = [-2.0_wp, -3.0_wp]
    x = [0.0_wp, 0.0_wp]
    c = [0.0_wp]
    x_lb = [0.0_wp, 0.0_wp]
    x_ub = [big, big]
    c_lb = [-big]
    c_ub = [3.0_wp]

    call qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    print '(A,2F12.7)', 'test_inequality: p      = ', p
    print '(A,2F12.7)', 'test_inequality: p_true = ', pexpect
    print '(A,F12.7)',  'test_inequality: lambda = ', lambda(1)
    print '(A,I0)',     'test_inequality: istat  = ', istat
    if (istat /= sqpopt_success) error stop 'test_inequality FAILED: istat'
    if (maxval(abs(p-pexpect)) > 1.0e-6_wp) error stop 'test_inequality FAILED: wrong p'
    print *, 'test_inequality PASSED'

    end subroutine test_inequality

    !> minimize 0.5*(p1^2+p2^2) - 2*p1 - 5*p2  s.t.  p1+p2 <= 1,  p2 <= 0.3
    !! both start inactive at p=0 and become active one after another
    !! (exercises two sequential "add to working set" steps, no drops).
    !! KKT (both active): p1=2+lam1, p2=5+lam1+lam2, p1+p2=1, p2=0.3
    !!   => p1=0.7, p2=0.3, lam1=-1.3, lam2=-3.4
    subroutine test_two_active_at_once()

    type(sqpopt_dense_qp_type) :: qp
    type(sqpopt_hessian_type)  :: hessian
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: g(2), x(2), c(1), x_lb(2), x_ub(2), c_lb(1), c_ub(1)
    real(wp) :: p(2), lambda(1)
    real(wp), parameter :: pexpect(2) = [0.7_wp, 0.3_wp]
    integer :: istat

    call hessian%initialize(n=2, max_history=5)
    jac%nrows = 1; jac%ncols = 2; jac%nnz = 2
    jac%irow = [1,1]; jac%icol = [1,2]; jac%val = [1.0_wp, 1.0_wp]
    g = [-2.0_wp, -5.0_wp]
    x = [0.0_wp, 0.0_wp]
    c = [0.0_wp]
    x_lb = [-big, -big]
    x_ub = [big, 0.3_wp]
    c_lb = [-big]
    c_ub = [1.0_wp]

    call qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    print '(A,2F12.7)', 'test_two_active_at_once: p      = ', p
    print '(A,2F12.7)', 'test_two_active_at_once: p_true = ', pexpect
    print '(A,F12.7)',  'test_two_active_at_once: lambda = ', lambda(1)
    print '(A,I0)',     'test_two_active_at_once: istat  = ', istat
    if (istat /= sqpopt_success) error stop 'test_two_active_at_once FAILED: istat'
    if (maxval(abs(p-pexpect)) > 1.0e-6_wp) error stop 'test_two_active_at_once FAILED: wrong p'
    print *, 'test_two_active_at_once PASSED'

    end subroutine test_two_active_at_once

end program test_qp_dense
