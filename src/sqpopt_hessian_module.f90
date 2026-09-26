!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Limited-memory quasi-Newton approximation to the Hessian of the
!  Lagrangian. Rather than storing a dense \( n \times n \) matrix, only
!  the last `max_history` step/gradient-change vector pairs \( (s,y) \)
!  are kept (`max_history` is a small constant, independent of `n`), and
!  Hessian(-inverse)-vector products are formed matrix-free using the
!  standard two-loop recursion. Supports BFGS (with a curvature-condition
!  skip rule) and SR1 updates.

    module sqpopt_hessian_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    type, public :: sqpopt_hessian_type
        !! stores and updates a limited-memory approximation to the
        !! Hessian of the Lagrangian (never forms a dense `n x n` matrix).

        integer :: n           = 0  !! problem size
        integer :: max_history = 0  !! number of `(s,y)` pairs retained (independent of `n`)
        integer :: n_history   = 0  !! number of pairs currently stored (`<= max_history`)
        logical :: use_sr1     = .false. !! if true, use the limited-memory SR1 update instead of BFGS

        real(wp), dimension(:,:), allocatable :: s     !! stored step vectors `dimension(n,max_history)`
        real(wp), dimension(:,:), allocatable :: y     !! stored Lagrangian gradient-change vectors `dimension(n,max_history)`
        real(wp), dimension(:),   allocatable :: rho   !! `1/(y^T s)` for each stored pair `dimension(max_history)` (BFGS)
        real(wp) :: gamma = 1.0_wp  !! scaling of the initial Hessian \( H_0 = \gamma I \)

        contains

        procedure, public :: initialize             => hessian_initialize
        procedure, public :: update_bfgs             => hessian_update_bfgs
        procedure, public :: update_sr1              => hessian_update_sr1
        procedure, public :: hv_product              => hessian_vector_product
        procedure, public :: inverse_vector_product  => hessian_inverse_vector_product
        procedure, public :: reset                   => hessian_reset

    end type sqpopt_hessian_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the limited-memory Hessian approximation (equivalent to
!  \( H_0 = \gamma I \), with no `(s,y)` pairs stored).

    subroutine hessian_initialize(me, n, max_history, use_sr1)

    class(sqpopt_hessian_type), intent(inout) :: me
    integer, intent(in) :: n            !! problem size
    integer, intent(in) :: max_history  !! number of `(s,y)` pairs to retain
    logical, intent(in), optional :: use_sr1  !! if true, use SR1 instead of BFGS (default `.false.`)

    me%n           = n
    me%max_history = max_history
    me%n_history   = 0
    me%gamma       = 1.0_wp
    me%use_sr1     = .false.
    if (present(use_sr1)) me%use_sr1 = use_sr1

    if (allocated(me%s))     deallocate(me%s)
    if (allocated(me%y))     deallocate(me%y)
    if (allocated(me%rho))   deallocate(me%rho)
    allocate(me%s(n,max_history), me%y(n,max_history))
    allocate(me%rho(max_history))

    end subroutine hessian_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  update the limited-memory Hessian approximation using the BFGS
!  update formula, given the step \( s = x_{k+1} - x_k \) and the change
!  in the Lagrangian gradient \( y = \nabla_x \mathcal{L}_{k+1} - \nabla_x \mathcal{L}_k \).
!  The update is skipped if the curvature condition \( s^T y \) is not
!  sufficiently positive (standard "cautious update" safeguard). The
!  oldest pair is discarded once `max_history` pairs are stored.

    subroutine hessian_update_bfgs(me, s, y)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: s  !! step vector `dimension(n)`
    real(wp), dimension(:), intent(in) :: y  !! Lagrangian gradient change `dimension(n)`

    real(wp) :: sty, yty

    ! a near-zero step carries no reliable curvature information and risks
    ! an ill-conditioned (huge `rho`) update, so skip it outright:
    if (norm2(s) <= 1.0e-10_wp) return

    sty = dot_product(s, y)

    ! skip the update if the curvature condition is not sufficiently satisfied:
    if (sty <= 1.0e-10_wp*max(norm2(s)*norm2(y), 1.0_wp)) return

    call hessian_push_pair(me, s, y)
    me%rho(me%n_history) = 1.0_wp/sty

    yty = dot_product(y, y)
    if (yty > 0.0_wp) me%gamma = sty/yty

    end subroutine hessian_update_bfgs
!*******************************************************************************

!*******************************************************************************
!>
!  update the limited-memory Hessian approximation using the symmetric
!  rank-1 (SR1) update formula. The update is skipped if the standard
!  SR1 safeguard \( |(y-Bs)^Ts| \ge \epsilon \lVert y-Bs \rVert \lVert s
!  \rVert \) is not satisfied. The oldest pair is discarded once
!  `max_history` pairs are stored.

    subroutine hessian_update_sr1(me, s, y)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: s  !! step vector `dimension(n)`
    real(wp), dimension(:), intent(in) :: y  !! Lagrangian gradient change `dimension(n)`

    real(wp), dimension(me%n) :: bs, w
    real(wp) :: denom, sty, yty

    ! a near-zero step carries no reliable curvature information and risks
    ! an ill-conditioned update, so skip it outright:
    if (norm2(s) <= 1.0e-10_wp) return

    call hessian_vector_product(me, s, bs)
    w = y - bs
    denom = dot_product(w, s)

    ! standard SR1 safeguard against a near-singular update:
    if (abs(denom) < 1.0e-8_wp*max(norm2(w)*norm2(s), 1.0e-12_wp)) return

    call hessian_push_pair(me, s, y)

    sty = dot_product(s, y)
    yty = dot_product(y, y)
    if (sty > 0.0_wp .and. yty > 0.0_wp) me%gamma = sty/yty

    end subroutine hessian_update_sr1
!*******************************************************************************

!*******************************************************************************
!>
!  push a new `(s,y)` pair into the circular history buffer, discarding
!  the oldest pair if `max_history` pairs are already stored.

    subroutine hessian_push_pair(me, s, y)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: s  !! step vector `dimension(n)`
    real(wp), dimension(:), intent(in) :: y  !! Lagrangian gradient change `dimension(n)`

    if (me%n_history == me%max_history) then
        ! buffer full: discard the oldest pair (column 1), shift the rest down:
        me%s(:,1:me%max_history-1)     = me%s(:,2:me%max_history)
        me%y(:,1:me%max_history-1)     = me%y(:,2:me%max_history)
        me%rho(1:me%max_history-1)     = me%rho(2:me%max_history)
    else
        me%n_history = me%n_history + 1
    end if
    me%s(:,me%n_history) = s
    me%y(:,me%n_history) = y

    end subroutine hessian_push_pair
!*******************************************************************************

!*******************************************************************************
!>
!  compute the matrix-free Hessian-vector product \( h_v = H v \), used
!  by the QP subproblem solver in place of an explicit dense matrix.
!  Uses the compact representation of either the BFGS or the SR1
!  matrix (Byrd, Nocedal & Schnabel, 1994), depending on `use_sr1`. Both
!  are rebuilt from the stored `(s,y)` pairs and the *current* scaling
!  `gamma` on every call, so the product is always consistent with them.

    subroutine hessian_vector_product(me, v, hv)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: v   !! input vector `dimension(n)`
    real(wp), dimension(:), intent(out) :: hv  !! result `dimension(n)`

    integer :: i, k, k2
    real(wp), dimension(:,:), allocatable :: mid   !! the `2k x 2k` middle matrix
    real(wp), dimension(:), allocatable :: rhs, sol
    real(wp) :: theta
    logical :: ok

    k = me%n_history

    if (me%use_sr1) then

        ! compact SR1 representation (Byrd, Nocedal & Schnabel, 1994):
        ! B*v = theta*v + Psi * M^{-1} * Psi^T v, with Psi = Y - theta*S
        ! and M = D + L + L^T - theta*S^T S:
        theta = 1.0_wp/me%gamma
        if (k == 0) then
            hv = theta*v
            return
        end if
        allocate(mid(k,k), rhs(k), sol(k))
        block
            integer :: p, q
            do p = 1, k
                rhs(p) = dot_product(me%y(:,p) - theta*me%s(:,p), v)
                do q = 1, k
                    ! s_max(p,q)^T y_min(p,q) covers D (p==q) and L + L^T (p/=q):
                    mid(p,q) = dot_product(me%s(:,max(p,q)), me%y(:,min(p,q))) &
                               - theta*dot_product(me%s(:,p), me%s(:,q))
                end do
            end do
        end block
        call hessian_solve_small_system(k, mid, rhs, sol, ok)
        hv = theta*v
        if (ok) then
            do i = 1, k
                hv = hv + (me%y(:,i) - theta*me%s(:,i))*sol(i)
            end do
        end if
        deallocate(mid, rhs, sol)

    else

        theta = 1.0_wp/me%gamma
        if (k == 0) then
            hv = theta*v
            return
        end if

        ! compact BFGS representation (Byrd, Nocedal & Schnabel, 1994):
        ! B*v = theta*v - [theta*S Y] * M^{-1} * [theta*S^T v; Y^T v]
        k2 = 2*k
        allocate(mid(k2,k2), rhs(k2), sol(k2))
        mid = 0.0_wp
        do i = 1, k
            rhs(i)   = theta*dot_product(me%s(:,i), v)
            rhs(k+i) = dot_product(me%y(:,i), v)
        end do
        block
            integer :: p, q
            do p = 1, k
                do q = 1, k
                    mid(p,q) = theta*dot_product(me%s(:,p), me%s(:,q))  !! theta*S^T S block
                end do
            end do
            do p = 1, k
                do q = 1, p-1
                    ! L(p,q) = s_p^T y_q for p>q goes in the upper-right block,
                    ! and L^T in the lower-left block:
                    mid(p,k+q) = dot_product(me%s(:,p), me%y(:,q))      !! L (upper-right block)
                    mid(k+q,p) = mid(p,k+q)                             !! L^T (lower-left block)
                end do
                mid(k+p,k+p) = -dot_product(me%s(:,p), me%y(:,p))       !! -D (diagonal)
            end do
        end block

        call hessian_solve_small_system(k2, mid, rhs, sol, ok)
        if (.not. ok) then
            hv = theta*v  !! fall back to the (safe) initial scaling if the small system is singular
        else
            hv = theta*v
            do i = 1, k
                hv = hv - theta*me%s(:,i)*sol(i) - me%y(:,i)*sol(k+i)
            end do
        end if
        deallocate(mid, rhs, sol)

    end if

    end subroutine hessian_vector_product
!*******************************************************************************

!*******************************************************************************
!>
!  solve the small (`k2 x k2`, `k2 = 2*max_history` at most) dense linear
!  system arising in the compact BFGS representation, via Gaussian
!  elimination with partial pivoting. `k2` is a small constant
!  independent of the problem size `n`.

    subroutine hessian_solve_small_system(k2, a, b, x, ok)

    integer, intent(in) :: k2 !! order of the small dense system
    real(wp), dimension(k2,k2), intent(inout) :: a !! coefficient matrix of the small dense system
    real(wp), dimension(k2), intent(inout)    :: b !! right-hand side vector of the small dense system
    real(wp), dimension(k2), intent(out)      :: x !! solution vector of the small dense system
    logical, intent(out) :: ok !! indicates whether the small system was solved successfully

    integer :: i, p, piv
    real(wp) :: fac, amax

    ok = .true.
    do p = 1, k2
        piv = p
        amax = abs(a(p,p))
        do i = p+1, k2
            if (abs(a(i,p)) > amax) then
                amax = abs(a(i,p))
                piv = i
            end if
        end do
        if (amax < 1.0e-13_wp) then
            ok = .false.
            x = 0.0_wp
            return
        end if
        if (piv /= p) then
            a([p,piv],:) = a([piv,p],:)
            b([p,piv])   = b([piv,p])
        end if
        do i = p+1, k2
            fac = a(i,p)/a(p,p)
            a(i,p:k2) = a(i,p:k2) - fac*a(p,p:k2)
            b(i)      = b(i) - fac*b(p)
        end do
    end do

    do i = k2, 1, -1
        x(i) = (b(i) - dot_product(a(i,i+1:k2), x(i+1:k2)))/a(i,i)
    end do

    end subroutine hessian_solve_small_system
!*******************************************************************************

!*******************************************************************************
!>
!  compute the matrix-free inverse-Hessian-vector product \( d = H^{-1} v \).
!  For the BFGS mode this uses the classic two-loop recursion. For the
!  SR1 mode (no simple closed-form inverse exists) this uses a small
!  internal conjugate-gradient solve driven by [[hessian_vector_product]].

    subroutine hessian_inverse_vector_product(me, v, d)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: v  !! input vector `dimension(n)`
    real(wp), dimension(:), intent(out) :: d  !! result `dimension(n)`

    real(wp), dimension(me%n_history) :: alpha !! temporary storage for the two-loop recursion coefficients in L-BFGS
    real(wp), dimension(me%n) :: q !! temporary vector used in the two-loop recursion
    integer :: i !! loop index for the two-loop recursion

    if (me%use_sr1) then
        call hessian_cg_solve(me, v, d)
        return
    end if

    ! standard L-BFGS two-loop recursion:
    q = v
    do i = me%n_history, 1, -1
        alpha(i) = me%rho(i)*dot_product(me%s(:,i), q)
        q = q - alpha(i)*me%y(:,i)
    end do

    d = me%gamma*q

    do i = 1, me%n_history
        block
            real(wp) :: beta
            beta = me%rho(i)*dot_product(me%y(:,i), d)
            d = d + me%s(:,i)*(alpha(i) - beta)
        end block
    end do

    end subroutine hessian_inverse_vector_product
!*******************************************************************************

!*******************************************************************************
!>
!  approximate \( H^{-1} v \) via a small conjugate-gradient solve of
!  \( H d = v \), using only [[hessian_vector_product]] (matrix-free).
!  Used for the SR1 mode, which has no simple closed-form inverse.

    subroutine hessian_cg_solve(me, v, d)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: v
    real(wp), dimension(:), intent(out) :: d

    real(wp), dimension(me%n) :: r, p, hp
    real(wp) :: rs_old, rs_new, alpha, beta, pHp
    integer :: iter, maxit

    d = 0.0_wp
    r = v
    p = r
    rs_old = dot_product(r, r)
    if (rs_old <= 0.0_wp) return

    maxit = min(me%n, max(2*me%max_history, 10))
    do iter = 1, maxit
        call hessian_vector_product(me, p, hp)
        pHp = dot_product(p, hp)
        if (abs(pHp) < 1.0e-14_wp) exit  !! indefinite/singular direction: stop and return current estimate
        alpha = rs_old/pHp
        d = d + alpha*p
        r = r - alpha*hp
        rs_new = dot_product(r, r)
        if (sqrt(rs_new) <= 1.0e-10_wp*max(norm2(v), 1.0_wp)) exit
        beta = rs_new/rs_old
        p = r + beta*p
        rs_old = rs_new
    end do

    end subroutine hessian_cg_solve
!*******************************************************************************

!*******************************************************************************
!>
!  reset the Hessian approximation, discarding all stored `(s,y)` pairs.

    subroutine hessian_reset(me)

    class(sqpopt_hessian_type), intent(inout) :: me

    me%n_history = 0
    me%gamma     = 1.0_wp

    end subroutine hessian_reset
!*******************************************************************************

    end module sqpopt_hessian_module
!*******************************************************************************
