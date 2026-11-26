!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Eigenvalues and eigenvectors of small dense symmetric matrices, for the
!  spectral decomposition of the limited-memory Hessians (see
!  [[sqpopt_spectral_module]], whose matrices are of the order of the
!  number of stored pairs).
!
!  Three methods:
!
!  * `sqpopt_eigen_ql`: Householder reduction to tridiagonal form (with the
!    reflections accumulated), then the implicit QL method with Wilkinson's
!    shift (the algorithm of EISPACK's `tred2`/`tql2`), always available,
!    in any real kind. About \( 4n^3/3 \) operations for the eigenvalues,
!    and \( 6n^3 \) or so with the eigenvectors.
!  * `sqpopt_eigen_jacobi`: the cyclic Jacobi method (plane rotations that
!    zero each off-diagonal element in turn, sweep after sweep, with
!    Rutishauser's threshold strategy), always available. The simplest and
!    the most accurate (small eigenvalues to high relative accuracy), but
!    about \( 20 n^3 \) operations.
!  * `sqpopt_eigen_lapack`: LAPACK's `DSYEV` (the same algorithm as
!    `sqpopt_eigen_ql`, in blocked code), only in a library built with
!    LAPACK (the `HAS_LAPACK` preprocessor directive, in double precision).
!
!  `sqpopt_eigen_auto` picks LAPACK where it is available, and
!  `sqpopt_eigen_ql` otherwise. (On the Hock-Schittkowski problems with the
!  convexified SR1 matrix, whose eigenproblems have orders up to 100,
!  the eigensolver's time was 3.4 s with Jacobi and 0.6 s with LAPACK.)

    module sqpopt_eigen_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

#ifdef HAS_LAPACK
#if defined(REAL32) || defined(REAL128)
#error "HAS_LAPACK needs the default (double precision) real kind: don't define REAL32 or REAL128 with it"
#endif
    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            !! LAPACK: all eigenvalues and, optionally, eigenvectors of a symmetric matrix
            import :: wp
            implicit none
            character, intent(in)    :: jobz      !! `'V'`: eigenvalues and eigenvectors, `'N'`: eigenvalues only
            character, intent(in)    :: uplo      !! which triangle of `a` is used (`'L'`)
            integer,   intent(in)    :: n         !! order
            integer,   intent(in)    :: lda       !! leading dimension of `a`
            real(wp),  intent(inout) :: a(lda,*)  !! the matrix, overwritten by the eigenvectors (with `'V'`)
            real(wp),  intent(out)   :: w(*)      !! the eigenvalues, in ascending order
            real(wp),  intent(inout) :: work(*)   !! workspace
            integer,   intent(in)    :: lwork     !! its size (`-1`: a query)
            integer,   intent(out)   :: info      !! `0` on success
        end subroutine dsyev
    end interface
    logical, parameter, public :: sqpopt_eigen_has_lapack = .true.  !! whether `sqpopt_eigen_lapack` is available
                                                                    !! (a library built with `HAS_LAPACK`)
#else
    logical, parameter, public :: sqpopt_eigen_has_lapack = .false. !! whether `sqpopt_eigen_lapack` is available
                                                                    !! (a library built with `HAS_LAPACK`)
#endif

    ! the methods:
    integer, parameter, public :: sqpopt_eigen_auto   = 0 !! LAPACK if available, else `sqpopt_eigen_ql`
    integer, parameter, public :: sqpopt_eigen_jacobi = 1 !! the cyclic Jacobi method (always available)
    integer, parameter, public :: sqpopt_eigen_lapack = 2 !! LAPACK's `DSYEV` (only with `HAS_LAPACK`)
    integer, parameter, public :: sqpopt_eigen_ql     = 3 !! tridiagonal reduction and the implicit QL method
                                                          !! (always available)

    ! statuses:
    integer, parameter, public :: sqpopt_eigen_ok              = 0 !! converged
    integer, parameter, public :: sqpopt_eigen_not_converged   = 1 !! no convergence (Jacobi: the sweep limit;
                                                                   !! LAPACK: `info > 0`)
    integer, parameter, public :: sqpopt_eigen_not_available   = 2 !! `sqpopt_eigen_lapack` without LAPACK

    public :: symmetric_eigen

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  all the eigenvalues, in ascending order, and the eigenvectors of a
!  symmetric matrix. Only the lower triangle of `a` is read; on return
!  its columns are the orthonormal eigenvectors (column `j` for `w(j)`),
!  unless `vectors` is false (then `a` is overwritten). `istat` is
!  `sqpopt_eigen_ok`, or another status (then `w` and `a` are not
!  meaningful).

    subroutine symmetric_eigen(a, w, istat, method, vectors)

    real(wp), dimension(:,:), intent(inout) :: a       !! the matrix `dimension(n,n)`; on return, the eigenvectors
    real(wp), dimension(:),   intent(out)   :: w       !! the eigenvalues, ascending `dimension(n)`
    integer,                  intent(out)   :: istat   !! status (see above)
    integer,  optional,       intent(in)    :: method  !! `sqpopt_eigen_auto` (default), `sqpopt_eigen_ql`,
                                                       !! `sqpopt_eigen_jacobi`, or `sqpopt_eigen_lapack`
    logical,  optional,       intent(in)    :: vectors !! whether to compute the eigenvectors (default `.true.`)

    integer :: m
    logical :: want

    m = sqpopt_eigen_auto
    if (present(method)) m = method
    if (m == sqpopt_eigen_auto) m = merge(sqpopt_eigen_lapack, sqpopt_eigen_ql, sqpopt_eigen_has_lapack)
    want = .true.
    if (present(vectors)) want = vectors

    select case (m)
    case (sqpopt_eigen_lapack)
        call eigen_lapack(a, w, istat, want)
    case (sqpopt_eigen_jacobi)
        call eigen_jacobi(a, w, istat)
    case default
        call eigen_ql(a, w, istat, want)
    end select

    end subroutine symmetric_eigen
!*******************************************************************************

!*******************************************************************************
!>
!  the cyclic Jacobi method (see the module documentation). Each rotation
!  zeroes `a(p,q)`, with \( t = \tan\phi \) the smaller root of
!  \( t^2 + 2\theta t - 1 = 0 \), \( \theta = (a_{qq} - a_{pp})/(2 a_{pq}) \)
!  (so the rotation angle is at most \( \pi/4 \), for stability). In the
!  first sweeps an element is only rotated away if it is larger than a
!  threshold; after that, an element negligible next to both diagonal
!  elements is set to zero without a rotation. The sweeps stop when the
!  off-diagonal part is negligible next to the matrix's norm.

    subroutine eigen_jacobi(a, w, istat)

    real(wp), dimension(:,:), intent(inout) :: a     !! the matrix `dimension(n,n)`; on return, the eigenvectors
    real(wp), dimension(:),   intent(out)   :: w     !! the eigenvalues, ascending `dimension(n)`
    integer,                  intent(out)   :: istat !! `sqpopt_eigen_ok` or `sqpopt_eigen_not_converged`

    integer, parameter :: max_sweeps = 60
    real(wp), dimension(:,:), allocatable :: v
    integer  :: n, i, j, p, q, sweep
    real(wp) :: off, fro, thresh, theta, t, c, s, tau, apq, g, akp, akq

    n = size(w)
    istat = sqpopt_eigen_ok
    if (n == 0) return
    ! (symmetric: from the lower triangle)
    do j = 1, n
        do i = 1, j-1
            a(i,j) = a(j,i)
        end do
    end do
    allocate(v(n,n))
    v = 0.0_wp
    do i = 1, n
        v(i,i) = 1.0_wp
    end do

    fro = sqrt(sum(a**2))
    do sweep = 1, max_sweeps
        off = 0.0_wp
        do q = 2, n
            do p = 1, q-1
                off = off + a(p,q)**2
            end do
        end do
        off = sqrt(2.0_wp*off)
        if (off <= epsilon(1.0_wp)*fro*1.0e-2_wp .or. off <= tiny(1.0_wp)) exit
        ! (in the first sweeps, only the larger elements: Rutishauser's threshold)
        thresh = 0.0_wp
        if (sweep < 4) thresh = 0.2_wp*off/real(n*n, wp)
        do q = 2, n
            do p = 1, q-1
                apq = a(p,q)
                g = 100.0_wp*abs(apq)
                ! (negligible next to both diagonal elements: set to zero)
                if (sweep > 4 .and. abs(a(p,p)) + g == abs(a(p,p)) .and. abs(a(q,q)) + g == abs(a(q,q))) then
                    a(p,q) = 0.0_wp
                    a(q,p) = 0.0_wp
                    cycle
                end if
                if (abs(apq) <= thresh .or. apq == 0.0_wp) cycle
                theta = (a(q,q) - a(p,p))/(2.0_wp*apq)
                if (abs(theta) > 1.0e30_wp) then
                    t = 0.5_wp/theta                      ! (theta^2 would overflow)
                else
                    t = sign(1.0_wp, theta)/(abs(theta) + sqrt(theta**2 + 1.0_wp))
                end if
                c   = 1.0_wp/sqrt(t**2 + 1.0_wp)
                s   = t*c
                tau = s/(1.0_wp + c)
                a(p,p) = a(p,p) - t*apq
                a(q,q) = a(q,q) + t*apq
                a(p,q) = 0.0_wp
                a(q,p) = 0.0_wp
                do i = 1, n
                    if (i == p .or. i == q) cycle
                    akp = a(i,p)
                    akq = a(i,q)
                    a(i,p) = akp - s*(akq + tau*akp)
                    a(i,q) = akq + s*(akp - tau*akq)
                    a(p,i) = a(i,p)
                    a(q,i) = a(i,q)
                end do
                do i = 1, n
                    akp = v(i,p)
                    akq = v(i,q)
                    v(i,p) = akp - s*(akq + tau*akp)
                    v(i,q) = akq + s*(akp - tau*akq)
                end do
            end do
        end do
    end do
    if (sweep > max_sweeps) istat = sqpopt_eigen_not_converged

    do i = 1, n
        w(i) = a(i,i)
    end do
    a = v
    call sort_ascending(a, w)

    end subroutine eigen_jacobi
!*******************************************************************************

!*******************************************************************************
!>
!  Householder reduction to tridiagonal form, then the implicit QL method
!  (see the module documentation).
!
!  The reduction applies the reflections \( H_k = I - 2 v v^T \),
!  \( k = 1, \dots, n-2 \), that zero column `k` below the subdiagonal, to
!  both sides of the trailing submatrix (\( A \leftarrow A - 2(vw^T + wv^T) \),
!  with \( p = Av \) and \( w = p - (v^Tp)v \)), and, for the eigenvectors,
!  accumulates their product \( Q \). The QL iterations then drive the
!  subdiagonal to zero by plane rotations, each sweep with the shift of the
!  eigenvalue of the leading 2 by 2 block nearer its corner (Wilkinson's),
!  deflating where a subdiagonal element is negligible next to its
!  neighbouring diagonal elements, and applying the rotations to \( Q \).

    subroutine eigen_ql(a, w, istat, vectors)

    real(wp), dimension(:,:), intent(inout) :: a       !! the matrix `dimension(n,n)`; on return, the eigenvectors
    real(wp), dimension(:),   intent(out)   :: w       !! the eigenvalues, ascending `dimension(n)`
    integer,                  intent(out)   :: istat   !! `sqpopt_eigen_ok` or `sqpopt_eigen_not_converged`
    logical,                  intent(in)    :: vectors !! whether to compute the eigenvectors

    integer, parameter :: max_iter = 60 !! QL iterations per eigenvalue
    real(wp), dimension(:,:), allocatable :: q
    real(wp), dimension(:), allocatable :: e, v, pv
    real(wp) :: alpha, vnorm, vp, f, g, r, s, c, b, p, dd
    integer  :: n, i, j, k, l, m, iter
    logical  :: underflow

    n = size(w)
    istat = sqpopt_eigen_ok
    if (n == 0) return
    ! (symmetric: from the lower triangle)
    do j = 1, n
        do i = 1, j-1
            a(i,j) = a(j,i)
        end do
    end do
    allocate(e(n), v(n), pv(n))
    if (vectors) then
        allocate(q(n,n))
        q = 0.0_wp
        do i = 1, n
            q(i,i) = 1.0_wp
        end do
    end if

    ! ---- tridiagonal form ----
    do k = 1, n-2
        ! the reflection that zeroes a(k+2:n,k): v = x - alpha*e1, normalized
        alpha = norm2(a(k+1:n,k))
        if (alpha == 0.0_wp) cycle
        if (a(k+1,k) > 0.0_wp) alpha = -alpha
        v(k+1:n) = a(k+1:n,k)
        v(k+1) = v(k+1) - alpha
        vnorm = norm2(v(k+1:n))
        if (vnorm == 0.0_wp) cycle
        v(k+1:n) = v(k+1:n)/vnorm
        ! the trailing submatrix: A22 <- A22 - 2 (v w' + w v'), w = p - (v'p) v, p = A22 v
        do j = k+1, n
            pv(j) = dot_product(a(k+1:n,j), v(k+1:n))
        end do
        vp = dot_product(v(k+1:n), pv(k+1:n))
        pv(k+1:n) = pv(k+1:n) - vp*v(k+1:n)
        do j = k+1, n
            a(k+1:n,j) = a(k+1:n,j) - 2.0_wp*(v(k+1:n)*pv(j) + pv(k+1:n)*v(j))
        end do
        ! column (and row) k
        a(k+1,k) = alpha
        a(k+2:n,k) = 0.0_wp
        a(k,k+1) = alpha
        a(k,k+2:n) = 0.0_wp
        ! Q <- Q H
        if (vectors) then
            do i = 1, n
                f = 2.0_wp*dot_product(q(i,k+1:n), v(k+1:n))
                q(i,k+1:n) = q(i,k+1:n) - f*v(k+1:n)
            end do
        end if
    end do
    do i = 1, n
        w(i) = a(i,i)
    end do
    do i = 1, n-1
        e(i) = a(i+1,i)
    end do
    e(n) = 0.0_wp

    ! ---- the implicit QL method ----
    do l = 1, n
        iter = 0
        do
            ! (the first negligible subdiagonal element at or after l)
            m = n
            do i = l, n-1
                dd = abs(w(i)) + abs(w(i+1))
                if (abs(e(i)) <= epsilon(1.0_wp)*dd) then
                    m = i
                    exit
                end if
            end do
            if (m == l) exit
            iter = iter + 1
            if (iter > max_iter) then
                istat = sqpopt_eigen_not_converged
                return
            end if
            ! Wilkinson's shift, from the leading 2 by 2 block
            g = (w(l+1) - w(l))/(2.0_wp*e(l))
            r = hypot(g, 1.0_wp)
            g = w(m) - w(l) + e(l)/(g + sign(r, g))
            s = 1.0_wp
            c = 1.0_wp
            p = 0.0_wp
            underflow = .false.
            do i = m-1, l, -1
                f = s*e(i)
                b = c*e(i)
                r = hypot(f, g)
                e(i+1) = r
                if (r == 0.0_wp) then
                    ! (an underflow: the rest of the sweep is skipped)
                    w(i+1) = w(i+1) - p
                    e(m) = 0.0_wp
                    underflow = .true.
                    exit
                end if
                s = f/r
                c = g/r
                g = w(i+1) - p
                r = (w(i) - g)*s + 2.0_wp*c*b
                p = s*r
                w(i+1) = g + p
                g = c*r - b
                if (vectors) then
                    do j = 1, n
                        f = q(j,i+1)
                        q(j,i+1) = s*q(j,i) + c*f
                        q(j,i)   = c*q(j,i) - s*f
                    end do
                end if
            end do
            if (underflow) cycle
            w(l) = w(l) - p
            e(l) = g
            e(m) = 0.0_wp
        end do
    end do

    if (vectors) then
        a = q
        call sort_ascending(a, w)
    else
        call sort_values(w)
    end if

    end subroutine eigen_ql
!*******************************************************************************

!*******************************************************************************
!>
!  sort the eigenvalues into ascending order (without eigenvectors).

    pure subroutine sort_values(w)

    real(wp), dimension(:), intent(inout) :: w !! the eigenvalues `dimension(n)`

    integer  :: i, j
    real(wp) :: t

    do i = 1, size(w) - 1
        j = minloc(w(i:), dim=1) + i - 1
        if (j == i) cycle
        t = w(i); w(i) = w(j); w(j) = t
    end do

    end subroutine sort_values
!*******************************************************************************
!>
!  LAPACK's `DSYEV` (see the module documentation), or
!  `sqpopt_eigen_not_available` in a library built without LAPACK.

    subroutine eigen_lapack(a, w, istat, vectors)

    real(wp), dimension(:,:), intent(inout) :: a       !! the matrix `dimension(n,n)`; on return, the eigenvectors
    real(wp), dimension(:),   intent(out)   :: w       !! the eigenvalues, ascending `dimension(n)`
    integer,                  intent(out)   :: istat   !! status (see [[symmetric_eigen]])
    logical,                  intent(in)    :: vectors !! whether to compute the eigenvectors

#ifdef HAS_LAPACK
    real(wp), dimension(:), allocatable :: work
    real(wp) :: query(1)
    integer :: n, info
    character :: jobz

    n = size(w)
    istat = sqpopt_eigen_ok
    if (n == 0) return
    jobz = merge('V', 'N', vectors)
    call dsyev(jobz, 'L', n, a, size(a,1), w, query, -1, info)
    allocate(work(max(1, int(query(1)))))
    call dsyev(jobz, 'L', n, a, size(a,1), w, work, size(work), info)
    if (info /= 0) istat = sqpopt_eigen_not_converged
#else
    if (vectors) a = 0.0_wp
    w = 0.0_wp
    istat = sqpopt_eigen_not_available
#endif

    end subroutine eigen_lapack
!*******************************************************************************

!*******************************************************************************
!>
!  sort the eigenvalues into ascending order, with their eigenvectors
!  (selection sort: the orders here are small).

    pure subroutine sort_ascending(v, w)

    real(wp), dimension(:,:), intent(inout) :: v !! the eigenvectors (columns) `dimension(n,n)`
    real(wp), dimension(:),   intent(inout) :: w !! the eigenvalues `dimension(n)`

    integer  :: i, j, r
    real(wp) :: t

    do i = 1, size(w) - 1
        j = minloc(w(i:), dim=1) + i - 1
        if (j == i) cycle
        t = w(i); w(i) = w(j); w(j) = t
        do r = 1, size(v,1)
            t = v(r,i); v(r,i) = v(r,j); v(r,j) = t
        end do
    end do

    end subroutine sort_ascending
!*******************************************************************************

    end module sqpopt_eigen_module
!*******************************************************************************
