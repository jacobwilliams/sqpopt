!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The spectral decomposition of the limited-memory Hessian approximations
!  (BFGS and SR1), from their compact representation, without forming an
!  `n x n` matrix (Burdakov, Gratton, Toint, Zikrin 2017; Brust, Burdakov,
!  Erway, Marcia, *Algorithm 1030: SC-SR1*, ACM TOMS 48(4), 2022, section 2.3).
!
!  The matrix is \( B = \theta I + U C U^T \), with \( \theta \) the initial
!  scaling plus the shift, `U` the `n x r` low-rank factor (`r` = `k` for SR1,
!  `2k` for BFGS, for `k` stored pairs) and \( C = \sigma M^{-1} \) the small
!  middle matrix (see [[sqpopt_hessian_module]]). With the eigendecomposition
!  of the Gram matrix \( U^TU = V D V^T \), keeping the `q` eigenvalues that
!  are not negligible, \( U = Q R \) with \( Q = U V D^{-1/2} \) orthonormal
!  (a thin QR factorization, from the stored inner products only) and
!  \( R = D^{1/2} V^T \). Then
!
!  $$ B = \theta I + Q (R C R^T) Q^T, \qquad R C R^T = W \hat\Lambda W^T $$
!
!  so the eigenvalues of `B` are \( \theta + \hat\lambda_i \) (`q` of them),
!  with the eigenvectors \( P = Q W = U (V D^{-1/2} W) \), and \( \theta \)
!  with multiplicity `n - q` (every vector orthogonal to `P`). It costs two
!  symmetric eigenproblems of order at most `r` (see [[sqpopt_eigen_module]])
!  and no operation on vectors of size `n`. The eigenvectors are never
!  formed: only the `r x q` matrix \( V D^{-1/2} W \) is kept, and products
!  with `P` are products with `U` (`O(nr)` each), so a trust-region step
!  costs about `4nr` operations, as many as an L-BFGS product.
!
!  Two uses (see [[sqpopt_iterate_module]] and [[sqpopt_qp_solver_module]]):
!  the smallest eigenvalue, to shift an indefinite SR1 matrix to a positive
!  definite one before a QP (`hessian%convexify`), and the minimizer of the
!  QP's objective within the step cap, an L2 trust-region subproblem solved
!  in closed form in the eigenvector basis (the unconstrained step of the
!  SR1 matrix).

    module sqpopt_spectral_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_eigen_module,   only: symmetric_eigen, sqpopt_eigen_ok

    implicit none

    private

    ! statuses:
    integer, parameter, public :: sqpopt_spectrum_ok            = 0 !! success
    integer, parameter, public :: sqpopt_spectrum_not_available = 1 !! not a limited-memory matrix (exact mode)
    integer, parameter, public :: sqpopt_spectrum_failed        = 2 !! an eigenproblem didn't converge (or LAPACK was
                                                                    !! asked for without it)
    integer, parameter, public :: sqpopt_spectrum_no_vectors    = 3 !! the eigenvectors weren't computed

    real(wp), parameter :: gram_tol = 1.0e-12_wp !! eigenvalues of \( U^TU \) below this times the size of the vectors
                                                 !! that `U` is formed from (see [[hessian_low_rank_scale]]) are
                                                 !! dropped: directions in which `U`'s columns are dependent, or only
                                                 !! rounding error

    type, public :: sqpopt_spectrum_type
        !! the spectral decomposition of a limited-memory Hessian (see the
        !! module documentation).
        integer  :: n = 0           !! order of the matrix
        integer  :: q = 0           !! number of the eigenvalues in `lambda` (the rank of the low-rank part)
        real(wp) :: theta = 0.0_wp  !! the eigenvalue of multiplicity `n - q` (the scaling plus the shift)
        real(wp), dimension(:),   allocatable :: lambda !! the other eigenvalues, ascending `dimension(q)`
        real(wp), dimension(:,:), allocatable :: coef   !! their eigenvectors are \( P = U \) `coef`, `U` the low-rank
                                                        !! factor of the Hessian it was computed for (if computed)
                                                        !! `dimension(r,q)`
        logical  :: has_vectors = .false. !! whether `coef` was computed

        contains

        procedure, public :: compute           => spectrum_compute
        procedure, public :: lambda_min        => spectrum_lambda_min
        procedure, public :: lambda_max        => spectrum_lambda_max
        procedure, public :: eigenvector       => spectrum_eigenvector
        procedure, public :: trust_region_step => spectrum_trust_region_step

    end type sqpopt_spectrum_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  the spectral decomposition of the limited-memory matrix `hessian` (with
!  its current shift), and, with `vectors`, the coefficients of the
!  eigenvectors of the `q` eigenvalues that are not `theta` (see the module
!  documentation). The products with the eigenvectors use `hessian`, which
!  must not change in between.

    subroutine spectrum_compute(me, hessian, istat, vectors)

    class(sqpopt_spectrum_type), intent(inout) :: me
    type(sqpopt_hessian_type),   intent(inout) :: hessian !! the Hessian approximation
    integer,                     intent(out)   :: istat   !! `sqpopt_spectrum_ok`, or another status (see the
                                                          !! module's constants)
    logical, optional,           intent(in)    :: vectors !! also the eigenvectors (default `.false.`)

    real(wp), dimension(:,:), allocatable :: gram, c, t1, rmat, k, rc
    real(wp), dimension(:),   allocatable :: d, khat
    integer, dimension(:), allocatable :: keep
    integer  :: r, j, i, q, eig_stat
    real(wp) :: dtol
    logical  :: want_vectors

    want_vectors = .false.
    if (present(vectors)) want_vectors = vectors
    me%has_vectors = .false.
    me%q = 0
    if (allocated(me%lambda)) deallocate(me%lambda)
    if (allocated(me%coef)) deallocate(me%coef)
    if (hessian%exact) then
        istat = sqpopt_spectrum_not_available
        return
    end if
    istat = sqpopt_spectrum_ok
    me%n = hessian%n
    me%theta = 1.0_wp/hessian%gamma + hessian%shift
    r = hessian%low_rank_size()
    allocate(me%lambda(0))
    if (r == 0) then
        call no_low_rank_part()
        return
    end if

    ! the Gram matrix U'U = V D V' (the thin QR factorization of U):
    allocate(gram(r,r), d(r))
    call hessian%low_rank_gram(gram)
    call symmetric_eigen(gram, d, eig_stat, hessian%eigen_solver)
    if (eig_stat /= sqpopt_eigen_ok) then
        istat = sqpopt_spectrum_failed
        return
    end if
    ! (the eigenvalues are ascending: keep the significant ones, at most n,
    ! the largest if there are more, which only round-off could give)
    allocate(keep(r))
    dtol = gram_tol*max(d(r), hessian%low_rank_scale())
    q = 0
    do j = r, 1, -1
        if (q == me%n) exit
        if (.not. (d(j) > dtol .and. d(j) > 0.0_wp)) exit
        q = q + 1
        keep(q) = j
    end do
    if (q == 0) then
        call no_low_rank_part()
        return
    end if

    ! T1 = V D^{-1/2} (so Q = U T1), and R = D^{1/2} V':
    allocate(t1(r,q), rmat(q,r))
    do j = 1, q
        t1(:,j)   = gram(:,keep(j))/sqrt(d(keep(j)))
        rmat(j,:) = gram(:,keep(j))*sqrt(d(keep(j)))
    end do

    ! K = R C R' = W Lambda W':
    allocate(c(r,r), rc(q,r), k(q,q), khat(q))
    call hessian%low_rank_middle_inverse(c)
    do i = 1, r
        do j = 1, q
            rc(j,i) = dot_product(rmat(j,:), c(:,i))
        end do
    end do
    do i = 1, q
        do j = 1, i
            k(i,j) = dot_product(rc(i,:), rmat(j,:))
            k(j,i) = k(i,j)
        end do
    end do
    call symmetric_eigen(k, khat, eig_stat, hessian%eigen_solver, vectors=want_vectors)
    if (eig_stat /= sqpopt_eigen_ok) then
        istat = sqpopt_spectrum_failed
        return
    end if
    deallocate(me%lambda)
    allocate(me%lambda(q))
    me%lambda = me%theta + khat
    me%q = q

    if (want_vectors) then
        ! P = Q W = U (T1 W):
        allocate(me%coef(r,q))
        do i = 1, q
            do j = 1, r
                me%coef(j,i) = dot_product(t1(j,:), k(:,i))
            end do
        end do
        me%has_vectors = .true.
    end if

    contains

        subroutine no_low_rank_part()
        !! the matrix is \( \theta I \): no eigenvalues besides \( \theta \), and no eigenvectors
        if (want_vectors) then
            allocate(me%coef(0, 0))
            me%has_vectors = .true.
        end if
        end subroutine no_low_rank_part

    end subroutine spectrum_compute
!*******************************************************************************

!*******************************************************************************
!>
!  the smallest eigenvalue (of all `n`).

    pure real(wp) function spectrum_lambda_min(me) result(lmin)

    class(sqpopt_spectrum_type), intent(in) :: me

    lmin = me%theta
    if (me%q > 0) then
        if (me%q < me%n) then
            lmin = min(me%lambda(1), me%theta)
        else
            lmin = me%lambda(1)
        end if
    end if

    end function spectrum_lambda_min
!*******************************************************************************

!*******************************************************************************
!>
!  the largest eigenvalue (of all `n`).

    pure real(wp) function spectrum_lambda_max(me) result(lmax)

    class(sqpopt_spectrum_type), intent(in) :: me

    lmax = me%theta
    if (me%q > 0) then
        if (me%q < me%n) then
            lmax = max(me%lambda(me%q), me%theta)
        else
            lmax = me%lambda(me%q)
        end if
    end if

    end function spectrum_lambda_max
!*******************************************************************************

!*******************************************************************************
!>
!  the eigenvector of `lambda(j)`, \( U \) `coef(:,j)` (`O(nr)`), with the
!  Hessian the decomposition was computed for.

    subroutine spectrum_eigenvector(me, hessian, j, v)

    class(sqpopt_spectrum_type), intent(in)  :: me
    type(sqpopt_hessian_type),   intent(in)  :: hessian !! the Hessian the decomposition was computed for
    integer,                     intent(in)  :: j       !! which (`1..q`)
    real(wp), dimension(:),      intent(out) :: v       !! the eigenvector `dimension(n)`

    call hessian%low_rank_multiply(me%coef(:,j), v)

    end subroutine spectrum_eigenvector
!*******************************************************************************

!*******************************************************************************
!>
!  the minimizer of \( g^Tp + \frac12 p^TBp \) subject to
!  \( \lVert p \rVert_2 \le \Delta \) (the L2 trust-region subproblem), from
!  the spectral decomposition (with its eigenvectors). In the eigenvector
!  basis the problem is separable: with \( a = P^Tg \) and \( g_\perp \) the
!  part of `g` orthogonal to `P`,
!
!  $$ p(\mu) = -P\,\mathrm{diag}\Big(\frac{1}{\lambda_i + \mu}\Big)\,a
!     - \frac{g_\perp}{\theta + \mu} $$
!
!  and the solution is \( p(0) \) if `B` is positive definite and
!  \( \lVert p(0) \rVert \le \Delta \), else \( p(\mu) \) with
!  \( \mu > \max(0, -\lambda_{\min}) \) and \( \lVert p(\mu) \rVert = \Delta \)
!  (a root of the secular equation \( 1/\lVert p(\mu) \rVert = 1/\Delta \),
!  found by Newton's method, safeguarded by bisection), or, in the "hard
!  case" (`g` orthogonal to the eigenspace of \( \lambda_{\min} \le 0 \), and
!  \( \lVert p(-\lambda_{\min}) \rVert < \Delta \)), \( p(-\lambda_{\min}) \)
!  plus a multiple of an eigenvector of \( \lambda_{\min} \) that brings it
!  to the boundary (Moré and Sorensen 1983).

    subroutine spectrum_trust_region_step(me, hessian, g, delta, p, istat, mu)

    class(sqpopt_spectrum_type), intent(in)  :: me
    type(sqpopt_hessian_type),   intent(in)  :: hessian !! the Hessian the decomposition was computed for (for the
                                                       !! products with its eigenvectors)
    real(wp), dimension(:),      intent(in)  :: g     !! the gradient `dimension(n)`
    real(wp),                    intent(in)  :: delta !! the radius \( \Delta > 0 \)
    real(wp), dimension(:),      intent(out) :: p     !! the step `dimension(n)`
    integer,                     intent(out) :: istat !! `sqpopt_spectrum_ok` or `sqpopt_spectrum_no_vectors`
    real(wp), optional,          intent(out) :: mu    !! the multiplier \( \mu \) of the constraint

    real(wp), parameter :: rtol = 1.0e-10_wp !! relative accuracy of \( \lVert p \rVert = \Delta \)
    integer,  parameter :: max_iter = 100
    real(wp), dimension(:), allocatable :: a, gperp, z, w, cw, e
    real(wp) :: gp2, lmin, sig, lo, hi, nrm, dnrm, tol, eps, tau, phi, dphi
    logical  :: has_perp, hard
    integer  :: i, it, n, q, m, r

    n = me%n
    q = me%q
    p = 0.0_wp
    if (present(mu)) mu = 0.0_wp
    if (.not. me%has_vectors) then
        istat = sqpopt_spectrum_no_vectors
        return
    end if
    istat = sqpopt_spectrum_ok
    has_perp = q < n

    ! the gradient in the eigenvector basis, a = P'g = coef' (U'g), and the
    ! rest of it, gperp = g - P a:
    r = size(me%coef, 1)
    allocate(a(q), gperp(n), w(r), cw(r), z(n))
    if (r > 0) call hessian%low_rank_multiply_transpose(g, w)
    do i = 1, q
        a(i) = dot_product(me%coef(:,i), w)
    end do
    gperp = g
    if (q > 0) then
        call coef_times(a, cw)
        call hessian%low_rank_multiply(cw, z)
        gperp = gperp - z
    end if
    gp2 = 0.0_wp
    if (has_perp) gp2 = dot_product(gperp, gperp)
    lmin = me%lambda_min()
    eps = epsilon(1.0_wp)*max(1.0_wp, abs(me%lambda_max()), abs(lmin))
    tol = sqrt(epsilon(1.0_wp))*max(norm2(g), tiny(1.0_wp))

    ! the interior solution:
    if (lmin > eps) then
        if (norm_p(0.0_wp) <= delta) then
            call form_p(0.0_wp)
            return
        end if
    end if

    ! the hard case: g (nearly) orthogonal to the eigenspace of lambda_min
    lo = max(0.0_wp, -lmin)
    hard = lmin <= eps
    if (hard) then
        do i = 1, q
            if (me%lambda(i) <= lmin + eps .and. abs(a(i)) > tol) hard = .false.
        end do
        if (has_perp .and. me%theta <= lmin + eps .and. sqrt(gp2) > tol) hard = .false.
    end if
    if (hard) then
        nrm = norm_p(lo, exclude=.true.)
        if (nrm <= delta) then
            call form_p(lo, exclude=.true.)
            ! (to the boundary along an eigenvector of lambda_min)
            if (q > 0 .and. (.not. has_perp .or. me%lambda(1) <= me%theta)) then
                call hessian%low_rank_multiply(me%coef(:,1), z)
            else
                ! (a unit vector orthogonal to P: a coordinate vector with P's part
                ! taken out; since the squared norms of P's rows add up to q < n,
                ! one of the first q+1 is at least partly outside P's span)
                allocate(e(n))
                do m = 1, min(q+1, n)
                    e = 0.0_wp
                    e(m) = 1.0_wp
                    z = e
                    if (q > 0) then
                        call hessian%low_rank_multiply_transpose(e, w)
                        do i = 1, q
                            a(i) = dot_product(me%coef(:,i), w)
                        end do
                        call coef_times(a, cw)
                        call hessian%low_rank_multiply(cw, e)
                        z = z - e
                    end if
                    if (norm2(z) >= 0.5_wp*sqrt(1.0_wp - real(q, wp)/real(n, wp))) exit
                end do
                z = z/norm2(z)
            end if
            tau = sqrt(max(delta**2 - norm2(p)**2, 0.0_wp))
            if (dot_product(g, z) > 0.0_wp) tau = -tau
            p = p + tau*z
            if (present(mu)) mu = lo
            return
        end if
    end if

    ! the root of 1/||p(mu)|| = 1/delta on (lo, hi]: ||p(mu)|| > delta near lo,
    ! and <= delta at hi = lo + ||g||/delta (where every lambda_i + mu >= ||g||/delta)
    hi = max(lo, 0.0_wp) + max(norm2(g)/delta, tiny(1.0_wp))
    if (lmin < 0.0_wp) hi = -lmin + norm2(g)/delta
    sig = hi
    do it = 1, max_iter
        nrm = norm_p(sig)
        if (abs(nrm - delta) <= rtol*delta) exit
        if (nrm > delta) then
            lo = sig
        else
            hi = sig
        end if
        ! Newton on phi(mu) = 1/||p|| - 1/delta: phi' = (sum c_i^2/(l_i+mu)^3)/||p||^3
        phi  = 1.0_wp/nrm - 1.0_wp/delta
        dnrm = dnorm_p(sig)
        dphi = dnrm/nrm**3
        if (dphi > 0.0_wp) then
            sig = sig - phi/dphi
        else
            sig = 0.5_wp*(lo + hi)
        end if
        if (.not. (sig > lo .and. sig < hi)) sig = 0.5_wp*(lo + hi)
    end do
    call form_p(sig)
    if (present(mu)) mu = sig

    contains

        real(wp) function norm_p(s, exclude)
        !! \( \lVert p(s) \rVert \) (with `exclude`, without the eigenspace of \( \lambda_{\min} \))
        real(wp),          intent(in) :: s       !! the multiplier \( \mu \)
        logical, optional, intent(in) :: exclude !! leave out the components of \( \lambda_{\min} \)
        real(wp) :: t
        integer :: j
        logical :: ex
        ex = .false.
        if (present(exclude)) ex = exclude
        t = 0.0_wp
        do j = 1, q
            if (ex .and. me%lambda(j) <= lmin + eps) cycle
            t = t + (a(j)/(me%lambda(j) + s))**2
        end do
        if (has_perp .and. .not. (ex .and. me%theta <= lmin + eps)) t = t + gp2/(me%theta + s)**2
        norm_p = sqrt(t)
        end function norm_p

        real(wp) function dnorm_p(s)
        !! \( \sum_i c_i^2/(\lambda_i + s)^3 \), so that \( d\lVert p \rVert/ds = -\) this \( /\lVert p \rVert \)
        real(wp), intent(in) :: s !! the multiplier \( \mu \)
        integer :: j
        dnorm_p = 0.0_wp
        do j = 1, q
            dnorm_p = dnorm_p + a(j)**2/(me%lambda(j) + s)**3
        end do
        if (has_perp) dnorm_p = dnorm_p + gp2/(me%theta + s)**3
        end function dnorm_p

        subroutine form_p(s, exclude)
        !! the step \( p(s) \) (with `exclude`, without the eigenspace of \( \lambda_{\min} \))
        real(wp),          intent(in) :: s       !! the multiplier \( \mu \)
        logical, optional, intent(in) :: exclude !! leave out the components of \( \lambda_{\min} \)
        real(wp), dimension(:), allocatable :: b
        integer :: j
        logical :: ex
        ex = .false.
        if (present(exclude)) ex = exclude
        p = 0.0_wp
        if (q > 0) then
            ! p = -P b = -U (coef b), b_j = a_j/(lambda_j + s)
            allocate(b(q))
            do j = 1, q
                b(j) = 0.0_wp
                if (ex .and. me%lambda(j) <= lmin + eps) cycle
                b(j) = a(j)/(me%lambda(j) + s)
            end do
            call coef_times(b, cw)
            call hessian%low_rank_multiply(cw, p)
            p = -p
        end if
        if (has_perp .and. .not. (ex .and. me%theta <= lmin + eps)) p = p - gperp/(me%theta + s)
        end subroutine form_p

        subroutine coef_times(b, cb)
        !! `cb` = `coef` times `b`
        real(wp), dimension(:), intent(in)  :: b  !! the vector `dimension(q)`
        real(wp), dimension(:), intent(out) :: cb !! the product `dimension(r)`
        integer :: j
        cb = 0.0_wp
        do j = 1, q
            cb = cb + me%coef(:,j)*b(j)
        end do
        end subroutine coef_times

    end subroutine spectrum_trust_region_step
!*******************************************************************************

    end module sqpopt_spectral_module
!*******************************************************************************
