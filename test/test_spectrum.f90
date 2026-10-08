program test_spectrum

    !! Tests of the symmetric eigensolver ([[sqpopt_eigen_module]]) and of
    !! the spectral decomposition of the limited-memory Hessians
    !! ([[sqpopt_spectral_module]]):
    !!
    !! * the eigensolver (the QL method, Jacobi, and LAPACK's `DSYEV` in a
    !!   build with LAPACK) on random symmetric matrices, and on ones with
    !!   repeated, zero, and widely spread eigenvalues: `A V = V Λ`, `V'V = I`,
    !!   the eigenvalues ascending, and the same eigenvalues without the
    !!   eigenvectors;
    !! * the spectral decomposition of random BFGS and SR1 matrices (more
    !!   pairs than variables included, so the low-rank part is rank
    !!   deficient; one variable, where the SR1 factor `y - θs` is only
    !!   rounding error; and with a shift): its eigenvalues against those of
    !!   the dense matrix, and its eigenvectors (`B v = λ v`, orthonormal);
    !! * the trust-region step: the optimality conditions
    !!   `(B + μI) p = -g`, `μ >= 0`, `||p|| <= Δ`, `μ (Δ - ||p||) = 0`, and
    !!   `B + μI` positive semidefinite, for positive definite and
    !!   indefinite matrices, small and large radii, and the hard case;
    !! * the exact Hessian is reported as not available.

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_eigen_module,   only: symmetric_eigen, sqpopt_eigen_jacobi, sqpopt_eigen_lapack, sqpopt_eigen_ql, &
                                     sqpopt_eigen_has_lapack, sqpopt_eigen_ok, sqpopt_eigen_not_available
    use sqpopt_spectral_module, only: sqpopt_spectrum_type, sqpopt_spectrum_ok, sqpopt_spectrum_not_available

    implicit none

    integer :: n_fail

    write(*,*) '----------------------------'
    write(*,*) 'test_spectrum'
    write(*,*) '----------------------------'

    call seed_rng()
    n_fail = 0
    call test_eigen(sqpopt_eigen_ql, 'QL')
    call test_eigen(sqpopt_eigen_jacobi, 'Jacobi')
    if (sqpopt_eigen_has_lapack) then
        call test_eigen(sqpopt_eigen_lapack, 'LAPACK')
    else
        call check_not_available()
    end if
    call test_decomposition()
    call test_trust_region()
    call test_exact()

    if (n_fail > 0) then
        print '(I0,A)', n_fail, ' check(s) failed'
        error stop 'test_spectrum FAILED'
    end if
    print '(A)', 'test_spectrum PASSED'

    contains

    subroutine seed_rng()
    !! a fixed random sequence (so failures are reproducible)
    integer :: n, i
    integer, dimension(:), allocatable :: seed
    call random_seed(size=n)
    allocate(seed(n))
    seed = [(4242 + 31*i, i=1,n)]
    call random_seed(put=seed)
    end subroutine seed_rng

    real(wp) function urand(a, b)
    !! uniform random number in `[a,b)`
    real(wp), intent(in) :: a !! lower end
    real(wp), intent(in) :: b !! upper end
    real(wp) :: r
    call random_number(r)
    urand = a + (b-a)*r
    end function urand

    subroutine check(ok, msg)
    !! count a failed check, and say what it was
    logical,          intent(in) :: ok  !! the condition
    character(len=*), intent(in) :: msg !! what was checked
    if (.not. ok) then
        n_fail = n_fail + 1
        print '(A)', '  FAILED: '//msg
    end if
    end subroutine check

    subroutine eigen_residuals(a, v, w, res, orth)
    !! `max |A v - λ v|` (relative to `max |A|`) and `max |V'V - I|`
    real(wp), dimension(:,:), intent(in)  :: a    !! the matrix
    real(wp), dimension(:,:), intent(in)  :: v    !! its eigenvectors
    real(wp), dimension(:),   intent(in)  :: w    !! its eigenvalues
    real(wp),                 intent(out) :: res  !! the residual
    real(wp),                 intent(out) :: orth !! the loss of orthogonality
    integer :: j, i
    real(wp), dimension(:), allocatable :: r
    res = 0.0_wp
    orth = 0.0_wp
    allocate(r(size(w)))
    do j = 1, size(w)
        r = matmul(a, v(:,j)) - w(j)*v(:,j)
        res = max(res, maxval(abs(r)))
        do i = 1, size(w)
            orth = max(orth, abs(dot_product(v(:,i), v(:,j)) - merge(1.0_wp, 0.0_wp, i == j)))
        end do
    end do
    res = res/max(maxval(abs(a)), tiny(1.0_wp))
    end subroutine eigen_residuals

    subroutine test_eigen(method, name)
    !! the eigensolver on random and special symmetric matrices
    integer,          intent(in) :: method !! the method
    character(len=*), intent(in) :: name   !! its name, for the output
    integer, parameter :: sizes(*) = [1, 2, 3, 5, 10, 30, 120]
    real(wp), dimension(:,:), allocatable :: a, v
    real(wp), dimension(:), allocatable :: w
    real(wp) :: res, orth, worst_res, worst_orth
    integer :: t, n, i, j, istat, kind
    worst_res = 0.0_wp
    worst_orth = 0.0_wp
    do t = 1, size(sizes)
        n = sizes(t)
        do kind = 1, 4
            if (allocated(a)) deallocate(a, v, w)
            allocate(a(n,n), v(n,n), w(n))
            select case (kind)
            case (1) ! random
                do j = 1, n
                    do i = j, n
                        a(i,j) = urand(-1.0_wp, 1.0_wp); a(j,i) = a(i,j)
                    end do
                end do
            case (2) ! identity plus rank one (n-1 repeated eigenvalues)
                v(:,1) = [(urand(-1.0_wp, 1.0_wp), i=1,n)]
                a = spread(v(:,1), 2, n)*spread(v(:,1), 1, n)
                do i = 1, n
                    a(i,i) = a(i,i) + 1.0_wp
                end do
            case (3) ! zero
                a = 0.0_wp
            case (4) ! widely spread (graded diagonal plus a small symmetric perturbation)
                do j = 1, n
                    do i = j, n
                        a(i,j) = 1.0e-3_wp*urand(-1.0_wp, 1.0_wp); a(j,i) = a(i,j)
                    end do
                    a(j,j) = 10.0_wp**(8.0_wp*real(j-1,wp)/max(real(n-1,wp),1.0_wp) - 4.0_wp)
                end do
            end select
            v = a
            call symmetric_eigen(v, w, istat, method)
            call check(istat == sqpopt_eigen_ok, name//': not converged')
            call eigen_residuals(a, v, w, res, orth)
            worst_res = max(worst_res, res)
            worst_orth = max(worst_orth, orth)
            call check(res <= 1.0e-11_wp*n, name//': A v /= lambda v')
            call check(orth <= 1.0e-11_wp*n, name//': V''V /= I')
            if (n > 1) call check(all(w(2:) >= w(:n-1)), name//': not ascending')
            block
                real(wp), dimension(:,:), allocatable :: b
                real(wp), dimension(:), allocatable :: w2
                b = a
                allocate(w2(n))
                call symmetric_eigen(b, w2, istat, method, vectors=.false.)
                call check(maxval(abs(w2 - w)) <= 1.0e-11_wp*n*max(1.0_wp, maxval(abs(a))), &
                           name//': eigenvalues without the vectors differ')
            end block
        end do
    end do
    print '(A,A,A,ES9.2,A,ES9.2)', 'eigensolver (', name, '): largest residual ', worst_res, &
        ', loss of orthogonality ', worst_orth
    end subroutine test_eigen

    subroutine check_not_available()
    !! asking for LAPACK in a build without it
    real(wp) :: a(2,2), w(2)
    integer :: istat
    a = reshape([2.0_wp, 1.0_wp, 1.0_wp, 2.0_wp], [2,2])
    call symmetric_eigen(a, w, istat, sqpopt_eigen_lapack)
    call check(istat == sqpopt_eigen_not_available, 'LAPACK without LAPACK: not reported')
    end subroutine check_not_available

    subroutine random_hessian(hess, n, k, sr1, shift)
    !! a limited-memory matrix from `k` random pairs `(s, A s)` (A symmetric:
    !! positive definite for BFGS, indefinite for SR1)
    type(sqpopt_hessian_type), intent(out) :: hess  !! the matrix
    integer,                   intent(in)  :: n     !! order
    integer,                   intent(in)  :: k     !! number of pairs offered
    logical,                   intent(in)  :: sr1   !! SR1 (else BFGS)
    real(wp),                  intent(in)  :: shift !! the shift
    real(wp), dimension(:,:), allocatable :: a
    real(wp), dimension(:), allocatable :: s
    integer :: i, j
    allocate(a(n,n), s(n))
    do j = 1, n
        do i = 1, n
            a(i,j) = urand(-1.0_wp, 1.0_wp)
        end do
    end do
    if (sr1) then
        a = 0.5_wp*(a + transpose(a))
    else
        a = matmul(transpose(a), a)
        do i = 1, n
            a(i,i) = a(i,i) + 0.1_wp
        end do
    end if
    call hess%initialize(n, k, use_sr1=sr1)
    do j = 1, k
        s = [(urand(-1.0_wp, 1.0_wp), i=1,n)]
        if (sr1) then
            call hess%update_sr1(s, matmul(a, s))
        else
            call hess%update_bfgs(s, matmul(a, s))
        end if
    end do
    hess%shift = shift
    end subroutine random_hessian

    subroutine test_decomposition()
    !! the spectral decomposition against the dense matrix's eigenvalues
    type(sqpopt_hessian_type)  :: hess
    type(sqpopt_spectrum_type) :: sp
    real(wp), dimension(:,:), allocatable :: h, v, pv
    real(wp), dimension(:), allocatable :: w, all_w, bv
    real(wp) :: err, worst, vres, scale
    integer :: trial, n, k, istat, i, j
    logical :: sr1
    worst = 0.0_wp
    vres = 0.0_wp
    do trial = 1, 60
        n = 1 + mod(trial*7, 13)
        if (trial > 50) n = 40
        k = 1 + mod(trial*5, 11)
        sr1 = mod(trial, 2) == 0
        call random_hessian(hess, n, k, sr1, merge(0.3_wp, 0.0_wp, mod(trial, 5) == 0))
        call sp%compute(hess, istat, vectors=.true.)
        call check(istat == sqpopt_spectrum_ok, 'spectrum: not computed')
        if (istat /= sqpopt_spectrum_ok) cycle
        ! the dense matrix's eigenvalues
        allocate(h(n,n), v(n,n), w(n), all_w(n), bv(n))
        call hess%dense(h)
        v = h
        call symmetric_eigen(v, w, istat, sqpopt_eigen_jacobi)
        ! the spectrum's: lambda, and theta n-q times
        all_w(1:sp%q) = sp%lambda
        all_w(sp%q+1:n) = sp%theta
        call sort(all_w)
        scale = max(1.0_wp, maxval(abs(w)))
        err = maxval(abs(all_w - w))/scale
        worst = max(worst, err)
        call check(err <= 1.0e-8_wp, 'spectrum: eigenvalues differ from the dense matrix''s')
        call check(abs(sp%lambda_min() - w(1)) <= 1.0e-8_wp*scale, 'spectrum: lambda_min')
        call check(abs(sp%lambda_max() - w(n)) <= 1.0e-8_wp*scale, 'spectrum: lambda_max')
        ! the eigenvectors
        allocate(pv(n, sp%q))
        do j = 1, sp%q
            call sp%eigenvector(hess, j, pv(:,j))
        end do
        do j = 1, sp%q
            call hess%hv_product(pv(:,j), bv)
            vres = max(vres, maxval(abs(bv - sp%lambda(j)*pv(:,j)))/scale)
            do i = 1, sp%q
                vres = max(vres, abs(dot_product(pv(:,i), pv(:,j)) - merge(1.0_wp, 0.0_wp, i == j)))
            end do
        end do
        deallocate(h, v, w, all_w, bv, pv)
    end do
    call check(vres <= 1.0e-8_wp, 'spectrum: eigenvectors')
    print '(A,ES9.2,A,ES9.2)', 'spectral decomposition: largest eigenvalue error ', worst, &
        ', eigenvector residual ', vres
    end subroutine test_decomposition

    subroutine sort(x)
    !! sort ascending (insertion sort)
    real(wp), dimension(:), intent(inout) :: x !! the values
    integer :: i, j
    real(wp) :: t
    do i = 2, size(x)
        t = x(i)
        j = i - 1
        do while (j >= 1)
            if (x(j) <= t) exit
            x(j+1) = x(j)
            j = j - 1
        end do
        x(j+1) = t
    end do
    end subroutine sort

    subroutine test_trust_region()
    !! the trust-region step's optimality conditions
    type(sqpopt_hessian_type)  :: hess
    type(sqpopt_spectrum_type) :: sp
    real(wp), dimension(:,:), allocatable :: h, v
    real(wp), dimension(:), allocatable :: w, g, p, r
    real(wp) :: delta, mu, scale, worst
    integer :: trial, n, k, istat, i, j, n_hard, n_boundary
    logical :: sr1, hard_trial
    worst = 0.0_wp
    n_hard = 0
    n_boundary = 0
    do trial = 1, 120
        n = 2 + mod(trial*3, 12)
        k = 1 + mod(trial*7, 9)
        sr1 = mod(trial, 3) /= 0
        hard_trial = sr1 .and. mod(trial, 4) == 0
        call random_hessian(hess, n, k, sr1, 0.0_wp)
        call sp%compute(hess, istat, vectors=.true.)
        if (istat /= sqpopt_spectrum_ok) then
            call check(.false., 'trust region: spectrum not computed')
            cycle
        end if
        allocate(h(n,n), v(n,n), w(n), g(n), p(n), r(n))
        call hess%dense(h)
        v = h
        call symmetric_eigen(v, w, istat, sqpopt_eigen_jacobi)
        g = [(urand(-1.0_wp, 1.0_wp), i=1,n)]
        if (hard_trial .and. w(1) < -1.0e-3_wp) then
            ! the hard case: g orthogonal to the eigenvectors of the smallest eigenvalue
            do j = 1, n
                if (w(j) <= w(1) + 1.0e-10_wp) g = g - dot_product(v(:,j), g)*v(:,j)
            end do
            n_hard = n_hard + 1
        end if
        delta = 10.0_wp**urand(-2.0_wp, 2.0_wp)
        if (hard_trial) delta = 1.0e3_wp
        call sp%trust_region_step(hess, g, delta, p, istat, mu)
        call check(istat == sqpopt_spectrum_ok, 'trust region: no step')
        scale = max(1.0_wp, maxval(abs(w)), mu)*max(1.0_wp, norm2(p)) + norm2(g)
        r = matmul(h, p) + mu*p + g
        worst = max(worst, norm2(r)/scale)
        call check(norm2(r) <= 1.0e-7_wp*scale, 'trust region: (B + mu I) p /= -g')
        call check(mu >= 0.0_wp, 'trust region: mu < 0')
        call check(norm2(p) <= delta*(1.0_wp + 1.0e-8_wp), 'trust region: ||p|| > delta')
        call check(mu*(delta - norm2(p)) <= 1.0e-7_wp*max(1.0_wp, mu)*delta, 'trust region: not complementary')
        call check(w(1) + mu >= -1.0e-8_wp*max(1.0_wp, maxval(abs(w))), 'trust region: B + mu I not PSD')
        if (norm2(p) >= delta*(1.0_wp - 1.0e-8_wp)) n_boundary = n_boundary + 1
        deallocate(h, v, w, g, p, r)
    end do
    print '(A,ES9.2,A,I0,A,I0,A)', 'trust-region step: largest residual ', worst, ' (', n_boundary, &
        ' on the boundary, ', n_hard, ' hard cases)'
    call check(n_hard > 0, 'trust region: no hard case tested')
    end subroutine test_trust_region

    subroutine test_exact()
    !! the exact Hessian has no spectral decomposition here
    type(sqpopt_hessian_type)  :: hess
    type(sqpopt_spectrum_type) :: sp
    integer :: istat
    call hess%initialize(2, 3)
    call hess%set_exact([1, 2], [1, 2])
    call sp%compute(hess, istat)
    call check(istat == sqpopt_spectrum_not_available, 'exact Hessian: not reported as not available')
    end subroutine test_exact

end program test_spectrum
