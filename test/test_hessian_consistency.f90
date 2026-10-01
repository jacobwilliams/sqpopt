program test_hessian_consistency

    !! Unit test of the limited-memory Hessian approximations in
    !! [[sqpopt_hessian_module]], using exact curvature pairs `y = A s` from
    !! a fixed symmetric positive-definite matrix `A`:
    !!
    !! * BFGS: the forward product `B*v` must satisfy the secant condition
    !!   for the newest pair (`B s_k = y_k`), and must be the exact inverse
    !!   of the two-loop-recursion product `H*v` (`B*(H*v) = v`), also once
    !!   the circular history buffer has wrapped around.
    !! * SR1: the compact L-SR1 product must satisfy the secant condition for
    !!   *every* stored pair, including after the oldest pair has been
    !!   discarded from a full history buffer.
    !! * The diagonal (BFGS with a wrapped buffer, and SR1) must match
    !!   `e_i^T B e_i`, also when it comes from the cache, and after a
    !!   further update (which must invalidate the cache).
    !! * The dense matrix (`dense`, which the dense QP solver uses) must
    !!   match the products with the unit vectors, in the same cases.
    !! * Powell-damped BFGS: a negative-curvature pair (`s^T y < 0`) is still
    !!   used (after damping, `s^T B s = 0.2 s^T B_old s`) and `B` stays
    !!   positive definite; with damping off, the pair is skipped instead.
    !! * SR1 inverse: `inverse_vector_product` (a matrix-free CG solve, see
    !!   `hessian_cg_solve`) must be the inverse of the forward product.
    !! * Exact mode: with `A` given as a lower-triangle sparsity pattern, the
    !!   product and diagonal must be those of `A`, the inverse product (CG)
    !!   its inverse, and the shift must grow on `reset` (with the product
    !!   then `(A + shift*I) v`) and decay (by 3, then to 0) on `set_values`.

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

    ! ---- BFGS with a full (wrapped-around) circular history buffer ----
    call h%initialize(n, 3)
    do k = 1, n_pairs
        call h%update_bfgs(ss(:,k), yy(:,k))
    end do
    if (h%n_history /= 3) error stop 'test_hessian_consistency FAILED: BFGS pairs not stored (wrapped buffer)'
    call h%hv_product(ss(:,n_pairs), bs)
    call h%inverse_vector_product(v, d)
    call h%hv_product(d, bv)
    print '(A,2ES10.2)', 'BFGS (wrapped buffer) secant, ||B*(H*v) - v|| = ', norm2(bs - yy(:,n_pairs)), norm2(bv - v)
    if (norm2(bs - yy(:,n_pairs)) > tol .or. norm2(bv - v) > tol*norm2(v)) &
        error stop 'test_hessian_consistency FAILED: BFGS with a wrapped buffer'

    ! ---- the diagonal (BFGS, wrapped buffer), and its cache ----
    call check_diagonal('BFGS')
    call check_dense('BFGS')
    call h%update_bfgs(ss(:,1), yy(:,1))  ! (wraps again: the cached diagonal must be recomputed)
    call check_diagonal('BFGS after an update')
    call check_dense('BFGS after an update')

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
        call check_dense('BFGS damped')
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
    call check_diagonal('SR1')
    call check_dense('SR1')

    ! ---- SR1 inverse (the matrix-free CG solve) ----
    call h%inverse_vector_product(v, d)
    call h%hv_product(d, bv)
    print '(A,ES10.2)', 'SR1 ||B*(B^-1*v) - v||          = ', norm2(bv - v)
    if (norm2(bv - v) > 1.0e-8_wp*norm2(v)) error stop 'test_hessian_consistency FAILED: SR1 forward/inverse mismatch'

    ! ---- exact mode ----
    block
        integer,  allocatable :: hrow(:), hcol(:)
        real(wp), allocatable :: hval(:)
        real(wp) :: diag(n), shift
        integer  :: j
        ! the lower triangle of A (its nonzeros only), row by row:
        allocate(hrow(0), hcol(0), hval(0))
        do i = 1, n
            do j = 1, i
                if (a(i,j) /= 0.0_wp) then
                    hrow = [hrow, i]; hcol = [hcol, j]; hval = [hval, a(i,j)]
                end if
            end do
        end do
        call h%initialize(n, 3)
        call h%set_exact(hrow, hcol)
        call h%set_values(hval, decay=.true.)
        call h%update_bfgs(ss(:,1), yy(:,1))     ! (quasi-Newton updates do nothing in this mode)
        call h%update_sr1(ss(:,2), yy(:,2))
        call h%hv_product(v, bv)
        print '(A,ES10.2)', 'exact ||H*v - A*v||             = ', norm2(bv - matmul(a, v))
        if (norm2(bv - matmul(a, v)) > tol) error stop 'test_hessian_consistency FAILED: exact product'
        call h%diagonal(diag)
        if (norm2(diag - [(a(i,i), i=1,n)]) > tol) error stop 'test_hessian_consistency FAILED: exact diagonal'
        call h%inverse_vector_product(v, d)
        print '(A,ES10.2)', 'exact ||A*(H^-1*v) - v||        = ', norm2(matmul(a, d) - v)
        if (norm2(matmul(a, d) - v) > 1.0e-8_wp*norm2(v)) error stop 'test_hessian_consistency FAILED: exact inverse'

        ! the shift: grows on reset (from shift_min*max|A_ij|, x10), and the
        ! product is then that of A + shift*I:
        call h%reset()
        shift = h%shift
        if (abs(shift - h%shift_min*maxval(abs(a))) > tol) error stop 'test_hessian_consistency FAILED: first shift'
        call h%reset()
        if (abs(h%shift - 10.0_wp*shift) > tol) error stop 'test_hessian_consistency FAILED: shift increase'
        shift = h%shift
        call h%hv_product(v, bv)
        print '(A,ES10.2)', 'exact shifted ||H*v - (A+dI)*v|| = ', norm2(bv - matmul(a, v) - shift*v)
        if (norm2(bv - matmul(a, v) - shift*v) > tol) error stop 'test_hessian_consistency FAILED: shifted product'
        call h%diagonal(diag)
        if (norm2(diag - [(a(i,i), i=1,n)] - shift) > tol) error stop 'test_hessian_consistency FAILED: shifted diagonal'
        call check_dense('exact shifted')

        ! ...and decays on set_values (unless decay=.false.), then drops to 0:
        call h%set_values(hval, decay=.false.)
        if (h%shift /= shift) error stop 'test_hessian_consistency FAILED: shift decayed without decay'
        call h%set_values(hval, decay=.true.)
        if (abs(h%shift - shift/3.0_wp) > tol) error stop 'test_hessian_consistency FAILED: shift decay'
        ! (10x the minimum: 10/3 and 10/9 are still above it, 10/27 is below, so 0)
        call h%set_values(hval, decay=.true.)
        if (abs(h%shift - shift/9.0_wp) > tol) error stop 'test_hessian_consistency FAILED: shift decay'
        call h%set_values(hval, decay=.true.)
        if (h%shift /= 0.0_wp) error stop 'test_hessian_consistency FAILED: small shift not dropped to 0'

        ! initialize switches the exact mode off again:
        call h%initialize(n, 3)
        if (h%exact) error stop 'test_hessian_consistency FAILED: exact mode kept by initialize'
    end block

    print '(A)', 'test_hessian_consistency PASSED'

contains

    subroutine check_diagonal(label)
    !! `h%diagonal` must match the diagonal of `B`, from products with the
    !! unit vectors (called twice, to check the cached value too)
    character(len=*), intent(in) :: label !! the case, for the message
    real(wp) :: diag(n), e(n), be(n), ref(n)
    integer :: j, pass
    do j = 1, n
        e = 0.0_wp; e(j) = 1.0_wp
        call h%hv_product(e, be)
        ref(j) = be(j)
    end do
    do pass = 1, 2
        call h%diagonal(diag)
        if (pass == 1) print '(A,ES10.2)', label//' diagonal error = ', norm2(diag - ref)
        if (norm2(diag - ref) > tol*norm2(ref)) error stop 'test_hessian_consistency FAILED: '//label//' diagonal'
    end do
    end subroutine check_diagonal

    subroutine check_dense(label)
    !! `h%dense` must be the matrix of `B`, from products with the unit
    !! vectors (it is built differently, see [[hessian_dense]])
    character(len=*), intent(in) :: label !! the case, for the message
    real(wp) :: hd(n,n), ref(n,n), e(n)
    integer :: j
    do j = 1, n
        e = 0.0_wp; e(j) = 1.0_wp
        call h%hv_product(e, ref(:,j))
    end do
    call h%dense(hd)
    print '(A,ES10.2)', label//' dense matrix error = ', norm2(hd - ref)
    if (norm2(hd - ref) > tol*norm2(ref)) error stop 'test_hessian_consistency FAILED: '//label//' dense matrix'
    if (any(hd /= transpose(hd))) error stop 'test_hessian_consistency FAILED: '//label//' dense matrix is not symmetric'
    end subroutine check_dense

end program test_hessian_consistency
