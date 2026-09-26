program test_hessian_consistency

    !! Unit test of the limited-memory Hessian approximations in
    !! [[sqpopt_hessian_module]], using exact curvature pairs `y = A s` from
    !! a fixed symmetric positive-definite matrix `A`:
    !!
    !! * BFGS: the forward product `B*v` must satisfy the secant condition
    !!   for the newest pair (`B s_k = y_k`), and must be the exact inverse
    !!   of the two-loop-recursion product `H*v` (`B*(H*v) = v`).
    !! * SR1: the compact L-SR1 product must satisfy the secant condition for
    !!   *every* stored pair, including after the oldest pair has been
    !!   discarded from a full history buffer.
    !! * Powell-damped BFGS: a negative-curvature pair (`s^T y < 0`) is still
    !!   used (after damping, `s^T B s = 0.2 s^T B_old s`) and `B` stays
    !!   positive definite; with damping off, the pair is skipped instead.

    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    integer,  parameter :: n = 5
    integer,  parameter :: n_pairs = 4
    real(wp), parameter :: tol = 1.0e-10_wp

    type(sqpopt_hessian_type) :: h
    real(wp) :: a(n,n), ss(n,n_pairs), yy(n,n_pairs), v(n), d(n), bv(n), bs(n)
    integer  :: i, k

    write(*,*) '----------------------------'
    write(*,*) 'test_hessian_consistency'
    write(*,*) '----------------------------'

    a = 0.0_wp
    do i = 1, n
        a(i,i) = real(i, wp)
    end do
    a(1,2) = 0.5_wp; a(2,1) = 0.5_wp
    a(3,5) = 0.3_wp; a(5,3) = 0.3_wp

    ss(:,1) = [ 1.0_wp,  0.2_wp, -0.3_wp,  0.4_wp,  0.1_wp]
    ss(:,2) = [ 0.1_wp,  1.0_wp,  0.5_wp, -0.2_wp,  0.3_wp]
    ss(:,3) = [-0.4_wp,  0.3_wp,  1.0_wp,  0.2_wp, -0.5_wp]
    ss(:,4) = [ 0.2_wp, -0.1_wp,  0.3_wp,  1.0_wp,  0.6_wp]
    do k = 1, n_pairs
        yy(:,k) = matmul(a, ss(:,k))
    end do
    v = [1.0_wp, -2.0_wp, 0.5_wp, 3.0_wp, -1.0_wp]

    ! ---- BFGS ----
    call h%initialize(n, 10)
    do k = 1, n_pairs
        call h%update_bfgs(ss(:,k), yy(:,k))
    end do
    if (h%n_history /= n_pairs) error stop 'test_hessian_consistency FAILED: BFGS pairs not stored'
    call h%hv_product(ss(:,n_pairs), bs)
    print '(A,ES10.2)', 'BFGS newest-pair secant error   = ', norm2(bs - yy(:,n_pairs))
    if (norm2(bs - yy(:,n_pairs)) > tol) error stop 'test_hessian_consistency FAILED: BFGS secant condition'
    call h%inverse_vector_product(v, d)
    call h%hv_product(d, bv)
    print '(A,ES10.2)', 'BFGS ||B*(H*v) - v||            = ', norm2(bv - v)
    if (norm2(bv - v) > tol*norm2(v)) error stop 'test_hessian_consistency FAILED: BFGS forward/inverse mismatch'

    ! ---- Powell damping: a negative-curvature pair ----
    block
        real(wp) :: s_neg(n), sbs_old, sbs_new, vbv
        call h%initialize(n, 10)
        call h%update_bfgs(ss(:,1), yy(:,1))
        call h%update_bfgs(ss(:,2), yy(:,2))
        s_neg = ss(:,3)
        call h%hv_product(s_neg, bs)
        sbs_old = dot_product(s_neg, bs)
        call h%update_bfgs(s_neg, -s_neg)
        if (h%n_history /= 3) error stop 'test_hessian_consistency FAILED: damped pair not stored'
        call h%hv_product(s_neg, bs)
        sbs_new = dot_product(s_neg, bs)
        print '(A,ES10.2)', 'damped BFGS s^T B s / (0.2 s^T B_old s) - 1 = ', sbs_new/(0.2_wp*sbs_old) - 1.0_wp
        if (abs(sbs_new/(0.2_wp*sbs_old) - 1.0_wp) > tol) error stop 'test_hessian_consistency FAILED: damping'
        do i = 1, n  ! positive definiteness, spot-checked on the unit vectors and v
            d = 0.0_wp; d(i) = 1.0_wp
            call h%hv_product(d, bv)
            vbv = dot_product(d, bv)
            if (vbv <= 0.0_wp) error stop 'test_hessian_consistency FAILED: damped B not positive definite'
        end do
        call h%hv_product(v, bv)
        if (dot_product(v, bv) <= 0.0_wp) error stop 'test_hessian_consistency FAILED: damped B not positive definite'

        h%damping = .false.
        call h%update_bfgs(ss(:,4), -ss(:,4))
        if (h%n_history /= 3) error stop 'test_hessian_consistency FAILED: undamped negative-curvature pair not skipped'
    end block

    ! ---- SR1 (history of 3, so the oldest of the 4 pairs is discarded) ----
    call h%initialize(n, 3, use_sr1=.true.)
    do k = 1, n_pairs
        call h%update_sr1(ss(:,k), yy(:,k))
    end do
    if (h%n_history /= 3) error stop 'test_hessian_consistency FAILED: SR1 pairs not stored'
    do k = 2, n_pairs
        call h%hv_product(ss(:,k), bs)
        print '(A,I0,A,ES10.2)', 'SR1 secant error, pair ', k, '         = ', norm2(bs - yy(:,k))
        if (norm2(bs - yy(:,k)) > tol) error stop 'test_hessian_consistency FAILED: SR1 secant condition'
    end do

    print '(A)', 'test_hessian_consistency PASSED'

end program test_hessian_consistency
