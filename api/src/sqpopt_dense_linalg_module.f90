!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Small, self-contained dense linear algebra helpers. Most are used only
!  by [[sqpopt_qp_dense_module]] (the dense QP solver mode), the only place
!  where dense `n x n` arrays are formed. The LU factorization and the
!  inertia are for the small matrices (of the order of the quasi-Newton
!  memory) of [[sqpopt_hessian_module]] and [[sqpopt_kkt_module]]. The
!  Bunch-Kaufman \( LDL^T \) factorization is for the dense backend of the
!  sparse solver ([[sqpopt_dense_ldl_module]]). Not linked to any
!  external dependency: classic, textbook Householder QR and modified
!  Cholesky, small enough to validate directly against known small
!  matrices.

    module sqpopt_dense_linalg_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    public :: dense_null_space
    public :: dense_modified_cholesky
    public :: dense_solve_cholesky
    public :: dense_cholesky_curvature
    public :: dense_lu_factor
    public :: dense_lu_solve
    public :: dense_symmetric_inertia
    public :: dense_ldl_factor
    public :: dense_ldl_solve

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  compute an orthonormal basis `Z` (`dimension(n,n_z)`) for the null
!  space of `A` (`dimension(m_a,n)`, `m_a<=n` assumed full row rank), via
!  Householder QR of \( A^T \): \( A^T = QR \), so `null(A)` is spanned by
!  the last `n-m_a` columns of `Q`. `Z` is allocated here (`intent(out),
!  allocatable`) with `dimension(n,n_z)`, `n_z = n-m_a`.

    subroutine dense_null_space(a, m_a, n, z, n_z)

    integer,                    intent(in)  :: m_a  !! number of rows of `a` (may be 0)
    integer,                    intent(in)  :: n    !! number of columns of `a` (and rows of `z`)
    real(wp), dimension(:,:),   intent(in)  :: a    !! `dimension(m_a,n)`
    real(wp), dimension(:,:), allocatable, intent(out) :: z    !! `dimension(n,n_z)`, `n_z = n-m_a`
    integer,                    intent(out) :: n_z  !! number of columns of `z`

    real(wp), dimension(:,:), allocatable :: q
    real(wp), dimension(:,:), allocatable :: r
    real(wp), dimension(:), allocatable :: v
    integer  :: i, j, k, kmax
    real(wp) :: alpha, vnorm, s
    real(wp), parameter :: tol = 1.0e-13_wp

    allocate(q(n,n), r(n,max(m_a,1)), v(n))

    n_z = n - m_a
    allocate(z(n, max(n_z,0)))

    ! Q starts as the identity; accumulate Householder reflectors into it:
    q = 0.0_wp
    do i = 1, n
        q(i,i) = 1.0_wp
    end do

    if (m_a <= 0) then
        if (n_z > 0) z(:,1:n_z) = q(:,1:n_z)
        return
    end if

    r(:,1:m_a) = transpose(a(1:m_a,1:n))

    kmax = min(m_a, n)
    do k = 1, kmax

        vnorm = norm2(r(k:n,k))
        if (vnorm < tol) cycle  !! column already (numerically) zero below the diagonal

        alpha = -sign(vnorm, r(k,k))
        v(k:n) = r(k:n,k)
        v(k)   = v(k) - alpha
        vnorm  = norm2(v(k:n))
        if (vnorm < tol) cycle
        v(k:n) = v(k:n)/vnorm

        ! apply the reflector to R (from the left):
        do j = k, m_a
            s = 2.0_wp*dot_product(v(k:n), r(k:n,j))
            r(k:n,j) = r(k:n,j) - s*v(k:n)
        end do

        ! accumulate the reflector into Q (apply from the right):
        do i = 1, n
            s = 2.0_wp*dot_product(q(i,k:n), v(k:n))
            q(i,k:n) = q(i,k:n) - s*v(k:n)
        end do

    end do

    if (n_z > 0) z(:,1:n_z) = q(:, m_a+1:n)

    end subroutine dense_null_space
!*******************************************************************************

!*******************************************************************************
!>
!  modified Cholesky factorization \( A = LL^T \) of a symmetric
!  `dimension(n,n)` matrix `a` that is not guaranteed positive definite
!  (the limited-memory BFGS/SR1 Hessian restricted to a subspace has no
!  such guarantee). Uses a simple Gill-Murray-Wright-style safeguard: any
!  pivot that would be non-positive is instead replaced by a small
!  positive floor, so the factorization always succeeds and returns a
!  (perturbed) positive-definite `L`.

    subroutine dense_modified_cholesky(a, n, l)

    integer,                  intent(in)  :: n !! order of the matrix `a`
    real(wp), dimension(n,n), intent(in)  :: a !! symmetric matrix to be factorized
    real(wp), dimension(n,n), intent(out) :: l !! lower-triangular Cholesky factor of `a`

    integer  :: i, j
    real(wp) :: piv
    real(wp), parameter :: delta = 1.0e-10_wp  !! floor applied to a non-positive pivot

    l = 0.0_wp
    do j = 1, n
        piv = a(j,j) - dot_product(l(j,1:j-1), l(j,1:j-1))
        if (piv < delta) piv = delta
        l(j,j) = sqrt(piv)
        do i = j+1, n
            l(i,j) = (a(i,j) - dot_product(l(i,1:j-1), l(j,1:j-1)))/l(j,j)
        end do
    end do

    end subroutine dense_modified_cholesky
!*******************************************************************************

!*******************************************************************************
!>
!  Cholesky factorization \( A = LL^T \) of a symmetric `dimension(n,n)`
!  matrix `a` that may not be positive definite, *without* perturbing it.
!  If `a` is (numerically) positive definite, `ok=.true.` and `l` is its
!  Cholesky factor. Otherwise the factorization stops at the first column
!  `j` whose pivot \( a_{jj} - \lVert l_{j,1:j-1} \rVert^2 \) is not
!  sufficiently positive (relative to the largest diagonal element), and
!  returns `ok=.false.` with a direction `d` of nonpositive curvature,
!  \( d^T A d \le \) (roughly) 0:
!  $$ d = \begin{bmatrix} -L_{11}^{-T} l_j \\ 1 \\ 0 \end{bmatrix} $$
!  for which \( d^T A d \) equals that pivot (the Schur complement).

    subroutine dense_cholesky_curvature(a, n, l, ok, d)

    integer,                  intent(in)  :: n  !! order of the matrix `a`
    real(wp), dimension(n,n), intent(in)  :: a  !! symmetric matrix to be factorized
    real(wp), dimension(n,n), intent(out) :: l  !! lower-triangular Cholesky factor (valid only if `ok`)
    logical,                  intent(out) :: ok !! true if `a` is positive definite
    real(wp), dimension(n),   intent(out) :: d  !! direction of nonpositive curvature (only if `.not. ok`)

    real(wp), parameter :: rel_tol = 1.0e-10_wp !! pivots below `rel_tol*max|a_jj|` count as nonpositive

    integer  :: i, j
    real(wp) :: piv, tol

    l  = 0.0_wp
    d  = 0.0_wp
    ok = .true.
    if (n == 0) return

    tol = 0.0_wp
    do j = 1, n
        tol = max(tol, abs(a(j,j)))
    end do
    tol = rel_tol*max(tol, tiny(1.0_wp))

    do j = 1, n
        piv = a(j,j) - dot_product(l(j,1:j-1), l(j,1:j-1))
        if (piv <= tol) then
            ! nonpositive curvature: d = [-L11^{-T} l_j; 1; 0]
            ok   = .false.
            d(j) = 1.0_wp
            do i = j-1, 1, -1
                d(i) = -(l(j,i) + dot_product(l(i+1:j-1,i), d(i+1:j-1)))/l(i,i)
            end do
            return
        end if
        l(j,j) = sqrt(piv)
        do i = j+1, n
            l(i,j) = (a(i,j) - dot_product(l(i,1:j-1), l(j,1:j-1)))/l(j,j)
        end do
    end do

    end subroutine dense_cholesky_curvature
!*******************************************************************************

!*******************************************************************************
!>
!  solve \( Ax=b \) given the Cholesky factor `l` (\( A=ll^T \), from
!  [[dense_modified_cholesky]]) via forward and back substitution.

    subroutine dense_solve_cholesky(l, n, b, x)

    integer,                  intent(in)  :: n !! order of the matrix `l`
    real(wp), dimension(n,n), intent(in)  :: l !! lower-triangular Cholesky factor of `a`
    real(wp), dimension(n),   intent(in)  :: b !! right-hand side vector of the linear system `Ax=b`
    real(wp), dimension(n),   intent(out) :: x !! solution vector of the linear system `Ax=b`

    real(wp), dimension(:), allocatable :: y
    integer :: i

    allocate(y(n))

    ! forward substitution: L y = b
    do i = 1, n
        y(i) = (b(i) - dot_product(l(i,1:i-1), y(1:i-1)))/l(i,i)
    end do

    ! back substitution: L^T x = y
    do i = n, 1, -1
        x(i) = (y(i) - dot_product(l(i+1:n,i), x(i+1:n)))/l(i,i)
    end do

    end subroutine dense_solve_cholesky
!*******************************************************************************

!*******************************************************************************
!>
!  in-place LU factorization with partial pivoting of a small dense
!  matrix `a` (e.g. the middle matrix of [[sqpopt_hessian_module]], whose order is
!  `2*max_history` at most, independent of `n`).
!  `ok` is false if a pivot is negligible relative to the matrix's largest
!  element.

    pure subroutine dense_lu_factor(a, piv, ok)

    real(wp), dimension(:,:), intent(inout) :: a   !! matrix, overwritten by its `L` (unit, below the diagonal) and `U` factors
    integer,  dimension(:),   intent(out)   :: piv !! `piv(p)` is the row swapped with row `p` at step `p`
    logical,                  intent(out)   :: ok  !! false if `a` is (numerically) singular

    integer :: i, p, k2
    real(wp) :: amax, tol

    k2 = size(a,1)
    ok = .true.
    tol = 1.0e-14_wp*max(maxval(abs(a)), tiny(1.0_wp))
    do p = 1, k2
        piv(p) = p - 1 + maxloc(abs(a(p:k2,p)), dim=1)
        amax = abs(a(piv(p),p))
        if (amax <= tol) then
            ok = .false.
            return
        end if
        if (piv(p) /= p) a([p,piv(p)],:) = a([piv(p),p],:)
        do i = p+1, k2
            a(i,p) = a(i,p)/a(p,p)
            a(i,p+1:k2) = a(i,p+1:k2) - a(i,p)*a(p,p+1:k2)
        end do
    end do

    end subroutine dense_lu_factor
!*******************************************************************************

!*******************************************************************************
!>
!  solve \( A x = b \) in place (`b` is overwritten by `x`), given the LU
!  factors from [[dense_lu_factor]].

    pure subroutine dense_lu_solve(a, piv, b)

    real(wp), dimension(:,:), intent(in)    :: a   !! LU factors
    integer,  dimension(:),   intent(in)    :: piv !! row pivots
    real(wp), dimension(:),   intent(inout) :: b   !! right-hand side, overwritten by the solution

    integer :: i, k2

    k2 = size(a,1)
    do i = 1, k2
        if (piv(i) /= i) b([i,piv(i)]) = b([piv(i),i])
    end do
    do i = 2, k2
        b(i) = b(i) - dot_product(a(i,1:i-1), b(1:i-1))
    end do
    do i = k2, 1, -1
        b(i) = (b(i) - dot_product(a(i,i+1:k2), b(i+1:k2)))/a(i,i)
    end do

    end subroutine dense_lu_solve
!*******************************************************************************

!*******************************************************************************
!>
!  in-place Bunch-Kaufman factorization \( P A P^T = L D L^T \) of the
!  symmetric matrix `a` (its lower triangle is used), with `L` unit lower
!  triangular and `D` block diagonal, of 1 by 1 and 2 by 2 blocks (Golub &
!  Van Loan, *Matrix Computations*, 4.4.4). The factors and the pivots are
!  stored exactly as LAPACK's `DSYTRF` stores them with `uplo = 'L'` (it is
!  the unblocked algorithm of LAPACK's `DSYTF2`): `D` on and next to the
!  diagonal, `L` below it, and `piv(k) > 0` for a 1 by 1 block (row and
!  column `k` were swapped with `piv(k)`), or `piv(k) = piv(k+1) < 0` for a
!  2 by 2 block in rows `k` and `k+1` (row `k+1` was swapped with
!  `-piv(k)`). So the inertia of `a` is that of `D` (Sylvester's law), and
!  [[dense_ldl_solve]] solves with the factors. It never fails: a column
!  that is exactly zero gets a zero 1 by 1 pivot, and is not eliminated.
!  The loops run down the columns. Costs \( n^3/3 \) operations.

    pure subroutine dense_ldl_factor(a, piv)

    real(wp), dimension(:,:), intent(inout) :: a   !! the symmetric matrix `dimension(n,n)`, whose lower triangle is
                                                   !! overwritten by `D` and `L`
    integer,  dimension(:),   intent(out)   :: piv !! the pivots `dimension(n)` (see above)

    real(wp), parameter :: alpha = (1.0_wp + sqrt(17.0_wp))/8.0_wp !! (the threshold that bounds the growth of
                                                                   !! the elements, Bunch & Kaufman's)
    integer :: n, k, kk, kp, kstep, imax, i, j
    real(wp) :: absakk, colmax, rowmax, t, d11, d22, d21, wk, wkp1

    n = size(a,1)
    k = 1
    do while (k <= n)
        kstep = 1
        absakk = abs(a(k,k))
        imax = k
        colmax = 0.0_wp
        do i = k+1, n
            if (abs(a(i,k)) > colmax) then
                imax = i
                colmax = abs(a(i,k))
            end if
        end do
        kp = k
        if (max(absakk, colmax) > 0.0_wp .and. absakk < alpha*colmax) then
            ! the largest off-diagonal element of row and column imax:
            rowmax = 0.0_wp
            do j = k, imax-1
                rowmax = max(rowmax, abs(a(imax,j)))
            end do
            do i = imax+1, n
                rowmax = max(rowmax, abs(a(i,imax)))
            end do
            if (absakk*rowmax >= alpha*colmax**2) then
                kp = k                      ! (a 1 by 1 pivot, with no interchange)
            else if (abs(a(imax,imax)) >= alpha*rowmax) then
                kp = imax                   ! (a 1 by 1 pivot: row and column imax)
            else
                kp = imax                   ! (a 2 by 2 pivot: rows k and imax)
                kstep = 2
            end if
        end if

        ! swap rows and columns kk and kp of the trailing matrix (its lower triangle):
        kk = k + kstep - 1
        if (kp /= kk) then
            do i = kp+1, n
                t = a(i,kk); a(i,kk) = a(i,kp); a(i,kp) = t
            end do
            do j = kk+1, kp-1
                t = a(j,kk); a(j,kk) = a(kp,j); a(kp,j) = t
            end do
            t = a(kk,kk); a(kk,kk) = a(kp,kp); a(kp,kp) = t
            if (kstep == 2) then
                t = a(k+1,k); a(k+1,k) = a(kp,k); a(kp,k) = t
            end if
        end if

        ! eliminate the pivot's columns from the trailing matrix:
        if (kstep == 1) then
            piv(k) = kp
            if (a(k,k) /= 0.0_wp) then
                d11 = 1.0_wp/a(k,k)
                do j = k+1, n
                    t = -d11*a(j,k)
                    do i = j, n
                        a(i,j) = a(i,j) + t*a(i,k)
                    end do
                end do
                do i = k+1, n
                    a(i,k) = d11*a(i,k)
                end do
            end if
        else
            piv(k)   = -kp
            piv(k+1) = -kp
            d21 = a(k+1,k)
            d11 = a(k+1,k+1)/d21
            d22 = a(k,k)/d21
            t   = 1.0_wp/(d11*d22 - 1.0_wp)
            d21 = t/d21
            do j = k+2, n
                wk   = d21*(d11*a(j,k) - a(j,k+1))
                wkp1 = d21*(d22*a(j,k+1) - a(j,k))
                do i = j, n
                    a(i,j) = a(i,j) - (a(i,k)*wk + a(i,k+1)*wkp1)
                end do
                a(j,k)   = wk
                a(j,k+1) = wkp1
            end do
        end if
        k = k + kstep
    end do

    end subroutine dense_ldl_factor
!*******************************************************************************

!*******************************************************************************
!>
!  solve \( A x = b \) in place (`b` is overwritten by `x`), given the
!  factors from [[dense_ldl_factor]] (or from LAPACK's `DSYTRF` with
!  `uplo = 'L'`: this is LAPACK's `DSYTRS` for one right-hand side). Every
!  1 by 1 pivot and 2 by 2 block of `D` must be nonsingular.

    pure subroutine dense_ldl_solve(a, piv, b)

    real(wp), dimension(:,:), intent(in)    :: a   !! the factors `dimension(n,n)`
    integer,  dimension(:),   intent(in)    :: piv !! the pivots `dimension(n)`
    real(wp), dimension(:),   intent(inout) :: b   !! the right-hand side `dimension(n)`, overwritten by the
                                                   !! solution

    integer :: n, k, kp, i
    real(wp) :: t, akm1k, akm1, ak, denom, bkm1, bk

    n = size(a,1)

    ! solve L D y = P b:
    k = 1
    do while (k <= n)
        if (piv(k) > 0) then
            kp = piv(k)
            if (kp /= k) then
                t = b(k); b(k) = b(kp); b(kp) = t
            end if
            do i = k+1, n
                b(i) = b(i) - b(k)*a(i,k)
            end do
            b(k) = b(k)/a(k,k)
            k = k + 1
        else
            kp = -piv(k)
            if (kp /= k+1) then
                t = b(k+1); b(k+1) = b(kp); b(kp) = t
            end if
            do i = k+2, n
                b(i) = b(i) - b(k)*a(i,k) - b(k+1)*a(i,k+1)
            end do
            akm1k = a(k+1,k)
            akm1  = a(k,k)/akm1k
            ak    = a(k+1,k+1)/akm1k
            denom = akm1*ak - 1.0_wp
            bkm1  = b(k)/akm1k
            bk    = b(k+1)/akm1k
            b(k)   = (ak*bkm1 - bk)/denom
            b(k+1) = (akm1*bk - bkm1)/denom
            k = k + 2
        end if
    end do

    ! solve L^T P x = y:
    k = n
    do while (k >= 1)
        if (piv(k) > 0) then
            do i = k+1, n
                b(k) = b(k) - a(i,k)*b(i)
            end do
            kp = piv(k)
            if (kp /= k) then
                t = b(k); b(k) = b(kp); b(kp) = t
            end if
            k = k - 1
        else
            do i = k+1, n
                b(k)   = b(k)   - a(i,k)*b(i)
                b(k-1) = b(k-1) - a(i,k-1)*b(i)
            end do
            kp = -piv(k)
            if (kp /= k) then
                t = b(k); b(k) = b(kp); b(kp) = t
            end if
            k = k - 2
        end if
    end do

    end subroutine dense_ldl_solve
!*******************************************************************************

!*******************************************************************************
!>
!  the inertia of the symmetric matrix `a`: its numbers of positive,
!  negative, and zero eigenvalues. The matrix is reduced to tridiagonal form
!  by Householder reflections (a congruence with an orthogonal matrix, so
!  the eigenvalues are unchanged), and the signs are counted from the
!  pivots of that form's \( LDL^T \) factorization (its Sturm sequence at
!  zero; Golub & Van Loan, *Matrix Computations*, 8.3 and 8.4), which has
!  the same inertia as the matrix (Sylvester's law). A pivot counts as zero
!  if it is below `zero_tol` (default \( 10^{-12} \)) times the matrix's
!  largest element. A zero
!  pivot that is coupled to the next row is not a zero eigenvalue (the
!  matrix \( [0, 1; 1, 0] \) has the eigenvalues \( \pm 1 \)): it is taken
!  with that row as a 2 by 2 pivot, as in the Bunch-Kaufman factorization,
!  and the block's two eigenvalues are counted. Costs \( O(n^3) \): meant
!  for small matrices.

    pure subroutine dense_symmetric_inertia(a, n_positive, n_negative, n_zero, zero_tol)

    real(wp), dimension(:,:), intent(in)  :: a          !! the symmetric matrix `dimension(n,n)`
    integer,                  intent(out) :: n_positive !! number of positive eigenvalues
    integer,                  intent(out) :: n_negative !! number of negative eigenvalues
    integer,                  intent(out) :: n_zero     !! number of zero eigenvalues
    real(wp), optional,       intent(in)  :: zero_tol   !! relative tolerance of a zero pivot (default `1e-12`)

    real(wp), parameter :: default_zero_tol = 1.0e-12_wp
    real(wp) :: tol
    real(wp), dimension(:,:), allocatable :: b
    real(wp), dimension(:), allocatable :: v, p, w
    real(wp) :: alpha, vnorm, q, small, e, d, det, mean, radius
    integer :: n, k, i

    allocate(b(size(a,1), size(a,1)), v(size(a,1)), p(size(a,1)), w(size(a,1)))

    n = size(a,1)
    n_positive = 0
    n_negative = 0
    n_zero     = 0
    if (n == 0) return
    b = 0.5_wp*(a + transpose(a))
    tol = default_zero_tol
    if (present(zero_tol)) tol = zero_tol
    small = tol*max(maxval(abs(b)), tiny(1.0_wp))

    ! Householder reduction to tridiagonal form:
    do k = 1, n-2
        alpha = norm2(b(k+1:n,k))
        if (alpha <= 0.0_wp) cycle
        if (b(k+1,k) > 0.0_wp) alpha = -alpha
        v(k+1:n) = b(k+1:n,k)
        v(k+1)   = v(k+1) - alpha
        vnorm = norm2(v(k+1:n))
        if (vnorm <= 0.0_wp) cycle
        v(k+1:n) = v(k+1:n)/vnorm
        ! B <- (I - 2vv^T) B (I - 2vv^T) on the trailing block:
        p(k+1:n) = matmul(b(k+1:n,k+1:n), v(k+1:n))
        w(k+1:n) = p(k+1:n) - dot_product(v(k+1:n), p(k+1:n))*v(k+1:n)
        do i = k+1, n
            b(k+1:n,i) = b(k+1:n,i) - 2.0_wp*(v(k+1:n)*w(i) + w(k+1:n)*v(i))
        end do
        b(k+1,k) = alpha
        b(k,k+1) = alpha
        b(k+2:n,k) = 0.0_wp
        b(k,k+2:n) = 0.0_wp
    end do

    ! the signs of the pivots of the tridiagonal matrix (`q` is the pivot of
    ! row `i`: its diagonal element, less the effect of the rows above):
    q = b(1,1)
    i = 1
    do while (i <= n)
        if (i < n) then
            e = b(i+1,i)
            d = b(i+1,i+1)
        else
            e = 0.0_wp
            d = 0.0_wp
        end if
        if (abs(q) >= small) then
            ! a 1 by 1 pivot:
            call count(q, n_positive, n_negative, n_zero)
            q = d - e**2/q
            i = i + 1
        else if (abs(e) < small) then
            ! a zero pivot that the next row doesn't depend on: a zero eigenvalue
            call count(q, n_positive, n_negative, n_zero)
            q = d
            i = i + 1
        else
            ! a zero pivot, coupled to the next row: a 2 by 2 pivot [q, e; e, d],
            ! with the eigenvalues mean +/- radius
            mean   = 0.5_wp*(q + d)
            radius = hypot(0.5_wp*(q - d), e)
            call count(mean + radius, n_positive, n_negative, n_zero)
            call count(mean - radius, n_positive, n_negative, n_zero)
            if (i + 2 <= n) then
                ! the next pivot: the block's inverse has q/det in its last element
                det = q*d - e**2
                if (det /= 0.0_wp) then
                    q = b(i+2,i+2) - b(i+2,i+1)**2*q/det
                else
                    q = b(i+2,i+2)
                end if
            end if
            i = i + 2
        end if
    end do

    contains

        pure subroutine count(eigenvalue, n_pos, n_neg, n_null)
        !! add an eigenvalue (or a pivot) to the count of its sign
        real(wp), intent(in)    :: eigenvalue !! the eigenvalue
        integer,  intent(inout) :: n_pos      !! number of positive eigenvalues so far
        integer,  intent(inout) :: n_neg      !! number of negative eigenvalues so far
        integer,  intent(inout) :: n_null     !! number of zero eigenvalues so far
        if (abs(eigenvalue) < small) then
            n_null = n_null + 1
        else if (eigenvalue > 0.0_wp) then
            n_pos = n_pos + 1
        else
            n_neg = n_neg + 1
        end if
        end subroutine count

    end subroutine dense_symmetric_inertia
!*******************************************************************************

    end module sqpopt_dense_linalg_module
!*******************************************************************************
