!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Small, self-contained dense linear algebra helpers used only by
!  [[sqpopt_qp_dense_module]] (the opt-in dense QP solver mode). No other
!  part of `sqpopt` uses dense arrays -- these are only ever formed/used
!  when the user explicitly selects the dense QP mode. Not linked to any
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

    real(wp), dimension(n,n)     :: q
    real(wp), dimension(n,max(m_a,1)) :: r
    real(wp), dimension(n) :: v
    integer  :: i, j, k, kmax
    real(wp) :: alpha, vnorm, s
    real(wp), parameter :: tol = 1.0e-13_wp

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

    integer,                  intent(in)  :: n
    real(wp), dimension(n,n), intent(in)  :: a
    real(wp), dimension(n,n), intent(out) :: l

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
!  solve \( Ax=b \) given the Cholesky factor `l` (\( A=ll^T \), from
!  [[dense_modified_cholesky]]) via forward and back substitution.

    subroutine dense_solve_cholesky(l, n, b, x)

    integer,                  intent(in)  :: n
    real(wp), dimension(n,n), intent(in)  :: l
    real(wp), dimension(n),   intent(in)  :: b
    real(wp), dimension(n),   intent(out) :: x

    real(wp), dimension(n) :: y
    integer :: i

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

    end module sqpopt_dense_linalg_module
!*******************************************************************************
