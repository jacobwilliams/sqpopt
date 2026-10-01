!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Limited-memory quasi-Newton approximation to the Hessian of the
!  Lagrangian. Rather than storing a dense \( n \times n \) matrix, only
!  the last `max_history` step/gradient-change vector pairs \( (s,y) \)
!  are kept (`max_history` is a small constant, independent of `n`), and
!  Hessian(-inverse)-vector products are formed matrix-free using the
!  compact representations and the standard two-loop recursion. Supports
!  (Powell-damped) BFGS and SR1 updates.
!
!  The pairs are kept in a circular buffer (no data is moved when the
!  oldest pair is discarded), together with their inner products
!  \( S^TS \) and \( S^TY \), which are updated in `O(nk)` as each pair is
!  added. The small `2k x 2k` (BFGS) or `k x k` (SR1) middle matrix of the
!  compact representation is formed from them and LU-factored only when
!  the pairs or the scaling change, not on every product, so each
!  [[hessian_vector_product]] costs only `O(nk)`. The diagonal (see
!  [[hessian_diagonal]]) is likewise kept until the pairs or the scaling
!  change. The dense QP solver, which needs the whole matrix, gets it from
!  [[hessian_dense]].
!
!  **Exact mode** (see [[hessian_set_exact]]): instead of the quasi-Newton
!  approximation, the user's sparse Hessian of the Lagrangian is used,
!  with the values set at each major iteration by [[hessian_set_values]].
!  The products are formed from its nonzeros, plus a shift
!  \( \delta I \) (\( \delta \ge 0 \)), a simple *primal inertia
!  correction* for an indefinite Hessian (as in IPOPT and Uno, but without
!  factorizing): wherever the solver would restart a quasi-Newton
!  approximation ([[hessian_reset]]: the QP step was not a descent
!  direction, a step failed, or a run of steps was very short), \( \delta \)
!  is increased tenfold, from `shift_min` times the size of the Hessian's
!  largest element; at each new major iteration that follows a good step,
!  it is divided by 3 (and set to 0 once below that minimum). The quasi-Newton updates do nothing
!  in this mode. With `options%inertia_control`, the shift that an
!  iteration needs is instead found from a factorization (see
!  [[sqpopt_inertia_module]]), and only the increases after a failed step
!  or a run of very short ones are carried over to the next iteration.
!
!  **For the sparse factorizations** (see [[sqpopt_kkt_module]]), the
!  quasi-Newton matrix is available in its compact form, a multiple of the
!  identity plus a matrix of low rank `r`:
!
!  $$ B = \theta I + \sigma\, U M^{-1} U^T $$
!
!  with \( U = [\theta S \; Y] \), \( \sigma = -1 \) (BFGS, `r = 2k`) or
!  \( U = Y - \theta S \), \( \sigma = +1 \) (SR1, `r = k`), and the middle
!  matrix \( M \) of [[factor_middle_matrix]]: see [[hessian_low_rank_size]],
!  [[hessian_low_rank_multiply]], [[hessian_low_rank_column]], and
!  [[hessian_low_rank_middle]]. The shift
!  \( \delta I \) also applies to the quasi-Newton matrices (it is zero
!  unless the inertia control sets it, for SR1).

    module sqpopt_hessian_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_dense_linalg_module, only: lu_factor => dense_lu_factor, lu_solve => dense_lu_solve

    implicit none

    private

    ! the Hessian modes (`options%hessian_mode`):
    integer, parameter, public :: sqpopt_hessian_bfgs   = 1  !! limited-memory (damped) BFGS quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_sr1    = 2  !! limited-memory symmetric rank-1 (SR1) quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_exact  = 3  !! user-supplied exact sparse Hessian of the Lagrangian (the
                                                             !! `hess` function of [[set_functions]], with the pattern of
                                                             !! [[set_hessian_sparsity]]; see the module documentation)

    type, public :: sqpopt_hessian_type
        !! stores and updates a limited-memory approximation to the
        !! Hessian of the Lagrangian (never stores a dense `n x n` matrix:
        !! [[hessian_dense]] writes one into the dense QP solver's array).

        integer :: n           = 0  !! problem size
        integer :: max_history = 0  !! number of `(s,y)` pairs retained (independent of `n`)
        integer :: n_history   = 0  !! number of pairs currently stored (`<= max_history`)
        logical :: use_sr1     = .false. !! if true, use the limited-memory SR1 update instead of BFGS
        logical :: damping     = .true.  !! if true, use Powell's damped BFGS update (see [[hessian_update_bfgs]]),
                                         !! else skip any update that fails the curvature condition

        integer :: first       = 1  !! column of `s`/`y` holding the oldest pair (the buffer is circular:
                                    !! the `i`-th oldest pair is in column [[pair_col]]`(i)`)
        real(wp), dimension(:,:), allocatable :: s     !! stored step vectors `dimension(n,max_history)`
        real(wp), dimension(:,:), allocatable :: y     !! stored Lagrangian gradient-change vectors `dimension(n,max_history)`
        real(wp), dimension(:),   allocatable :: rho   !! `1/(y^T s)` for each stored pair `dimension(max_history)` (BFGS)
        real(wp), dimension(:,:), allocatable :: ss    !! `ss(a,b)` \( = s_a^Ts_b \) for the pairs in columns `a`, `b`
                                                       !! `dimension(max_history,max_history)`
        real(wp), dimension(:,:), allocatable :: sy    !! `sy(a,b)` \( = s_a^Ty_b \) for the pairs in columns `a`, `b`
                                                       !! `dimension(max_history,max_history)`
        real(wp) :: gamma  = 1.0_wp !! scaling of the initial Hessian \( H_0 = \gamma I \)
        real(wp) :: gamma0 = 1.0_wp !! `gamma` before any update (and after a [[hessian_reset]])

        ! exact mode (see the module documentation):
        logical  :: exact = .false. !! if true, use the user's sparse Hessian (see [[hessian_set_exact]])
        real(wp) :: shift_min = 1.0e-4_wp !! smallest nonzero shift, relative to \( \max(1,\max|H_{ij}|) \)
        real(wp) :: shift_max = 1.0e10_wp !! largest shift, relative to the same
        real(wp) :: shift = 0.0_wp  !! the current shift \( \delta \) (added to the diagonal in every mode; in the
                                    !! quasi-Newton modes it is zero unless the inertia control sets it)
        integer, dimension(:), allocatable :: h_irow !! sparsity pattern: row indices
        integer, dimension(:), allocatable :: h_icol !! sparsity pattern: column indices
        real(wp), dimension(:), allocatable :: h_val !! current nonzero values

        ! cached LU factorization of the compact representation's middle matrix
        ! (internal; rebuilt by [[hessian_vector_product]] when `mid_valid` is false):
        logical  :: mid_valid = .false.  !! whether `mid_lu`/`mid_piv` match the current pairs and `gamma`
        logical  :: mid_ok    = .false.  !! whether the middle matrix is nonsingular
        real(wp), dimension(:,:), allocatable :: mid_lu  !! LU factors of the middle matrix
        integer,  dimension(:),   allocatable :: mid_piv !! row pivots of the LU factorization
        logical  :: diag_valid = .false. !! whether `diag` matches the current factorization
        real(wp), dimension(:),   allocatable :: diag    !! the diagonal of the approximation (see [[hessian_diagonal]])

        contains

        procedure, public :: initialize             => hessian_initialize
        procedure, public :: update_bfgs             => hessian_update_bfgs
        procedure, public :: update_sr1              => hessian_update_sr1
        procedure, public :: hv_product              => hessian_vector_product
        procedure, public :: diagonal                => hessian_diagonal
        procedure, public :: dense                   => hessian_dense
        procedure, public :: inverse_vector_product  => hessian_inverse_vector_product
        procedure, public :: reset                   => hessian_reset
        procedure, public :: set_exact               => hessian_set_exact
        procedure, public :: set_values              => hessian_set_values
        procedure, public :: magnitude               => hessian_size
        procedure, public :: low_rank_size           => hessian_low_rank_size
        procedure, public :: low_rank_multiply       => hessian_low_rank_multiply
        procedure, public :: low_rank_column         => hessian_low_rank_column
        procedure, public :: low_rank_multiply_transpose => hessian_low_rank_multiply_transpose
        procedure, public :: low_rank_middle         => hessian_low_rank_middle

    end type sqpopt_hessian_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the limited-memory Hessian approximation (equivalent to
!  \( H_0 = \gamma I \), with no `(s,y)` pairs stored).

    subroutine hessian_initialize(me, n, max_history, use_sr1, scale0)

    class(sqpopt_hessian_type), intent(inout) :: me
    integer, intent(in) :: n            !! problem size
    integer, intent(in) :: max_history  !! number of `(s,y)` pairs to retain
    logical, intent(in), optional :: use_sr1  !! if true, use SR1 instead of BFGS (default `.false.`)
    real(wp), intent(in), optional :: scale0  !! initial Hessian approximation \( B_0 = \) `scale0` \( I \) (default 1)

    me%n           = n
    me%max_history = max_history
    me%n_history   = 0
    me%first       = 1
    me%gamma0      = 1.0_wp
    if (present(scale0)) me%gamma0 = 1.0_wp/scale0
    me%gamma       = me%gamma0
    me%mid_valid   = .false.
    me%use_sr1     = .false.
    if (present(use_sr1)) me%use_sr1 = use_sr1
    me%exact       = .false.
    me%shift       = 0.0_wp

    if (allocated(me%s))     deallocate(me%s)
    if (allocated(me%y))     deallocate(me%y)
    if (allocated(me%rho))   deallocate(me%rho)
    if (allocated(me%ss))    deallocate(me%ss)
    if (allocated(me%sy))    deallocate(me%sy)
    allocate(me%s(n,max_history), me%y(n,max_history))
    allocate(me%rho(max_history), me%ss(max_history,max_history), me%sy(max_history,max_history))
    me%diag_valid  = .false.

    end subroutine hessian_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  update the limited-memory Hessian approximation using the BFGS
!  update formula, given the step \( s = x_{k+1} - x_k \) and the change
!  in the Lagrangian gradient \( y = \nabla_x \mathcal{L}_{k+1} - \nabla_x \mathcal{L}_k \).
!
!  The Hessian of the Lagrangian need not be positive definite, so
!  \( s^T y \) can be small or negative. If `damping` is on (the
!  default), Powell's damping (as in `slsqp` and Nocedal & Wright, Proc.
!  18.2) is applied: when \( s^T y < 0.2\, s^T B s \), `y` is replaced by
!  \( \theta y + (1-\theta) B s \) with
!  \( \theta = 0.8\, s^T B s / (s^T B s - s^T y) \), so that
!  \( s^T y = 0.2\, s^T B s > 0 \) and the update keeps `B` positive
!  definite while still using the new curvature information. Otherwise,
!  or if `damping` is off, the update is skipped when \( s^T y \) is not
!  sufficiently positive. The oldest pair is discarded once `max_history`
!  pairs are stored.

    subroutine hessian_update_bfgs(me, s, y)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: s  !! step vector `dimension(n)`
    real(wp), dimension(:), intent(in) :: y  !! Lagrangian gradient change `dimension(n)`

    real(wp), dimension(size(s)) :: y_used, bs
    real(wp) :: sty, yty, sbs, theta

    if (me%exact) return

    ! a near-zero step carries no reliable curvature information and risks
    ! an ill-conditioned (huge `rho`) update, so skip it outright:
    if (norm2(s) <= 1.0e-10_wp) return

    y_used = y
    sty    = dot_product(s, y)

    if (me%damping) then
        ! Powell's damping: blend `y` with `B*s` so that `s^T y = 0.2 s^T B s`
        ! whenever the curvature condition is (nearly) violated:
        call hessian_vector_product(me, s, bs)
        sbs = dot_product(s, bs)
        if (sbs > 0.0_wp .and. sty < 0.2_wp*sbs) then
            theta  = 0.8_wp*sbs/(sbs - sty)
            y_used = theta*y + (1.0_wp - theta)*bs
            sty    = dot_product(s, y_used)
        end if
    end if

    ! skip the update if the curvature condition is not sufficiently satisfied:
    if (sty <= 1.0e-10_wp*max(norm2(s)*norm2(y_used), 1.0_wp)) return

    call hessian_push_pair(me, s, y_used)
    me%rho(pair_col(me, me%n_history)) = 1.0_wp/sty

    yty = dot_product(y_used, y_used)
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

    if (me%exact) return

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
!  the oldest pair if `max_history` pairs are already stored (by
!  overwriting its column and advancing `first`, so no data is moved), and
!  update the inner products `ss` and `sy` with the new pair (`O(nk)`).

    subroutine hessian_push_pair(me, s, y)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: s  !! step vector `dimension(n)`
    real(wp), dimension(:), intent(in) :: y  !! Lagrangian gradient change `dimension(n)`

    integer :: i, j, c

    if (me%n_history == me%max_history) then
        me%first = mod(me%first, me%max_history) + 1   ! the oldest pair's column becomes the newest
    else
        me%n_history = me%n_history + 1
    end if
    c = pair_col(me, me%n_history)
    me%s(:,c) = s
    me%y(:,c) = y
    do i = 1, me%n_history
        j = pair_col(me, i)
        me%ss(c,j) = dot_product(s, me%s(:,j))
        me%ss(j,c) = me%ss(c,j)
        me%sy(c,j) = dot_product(s, me%y(:,j))
        me%sy(j,c) = dot_product(me%s(:,j), y)
    end do
    me%mid_valid = .false.

    end subroutine hessian_push_pair
!*******************************************************************************

!*******************************************************************************
!>
!  the column of `s`/`y`/`rho` holding the `i`-th oldest stored pair
!  (`i=1` is the oldest, `i=n_history` the newest).

    pure integer function pair_col(me, i)

    class(sqpopt_hessian_type), intent(in) :: me
    integer,                    intent(in) :: i !! age rank of the pair (`1` = oldest, `n_history` = newest)

    pair_col = mod(me%first + i - 2, me%max_history) + 1

    end function pair_col
!*******************************************************************************

!*******************************************************************************
!>
!  compute the matrix-free Hessian-vector product \( h_v = H v \), used
!  by the QP subproblem solver in place of an explicit dense matrix.
!  Uses the compact representation of either the BFGS or the SR1
!  matrix (Byrd, Nocedal & Schnabel, 1994), depending on `use_sr1`, with
!  its small middle matrix factored once (see [[factor_middle_matrix]])
!  and reused until the pairs or `gamma` change.

    subroutine hessian_vector_product(me, v, hv)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: v   !! input vector `dimension(n)`
    real(wp), dimension(:), intent(out) :: hv  !! result `dimension(n)`

    integer :: i, k, c
    real(wp), dimension(:), allocatable :: w
    real(wp) :: theta

    if (me%exact) then
        ! H*v from the nonzeros (each off-diagonal one stands for both
        ! (i,j) and (j,i)), plus the shift:
        hv = me%shift*v
        do k = 1, size(me%h_val)
            i = me%h_irow(k)
            c = me%h_icol(k)
            hv(i) = hv(i) + me%h_val(k)*v(c)
            if (i /= c) hv(c) = hv(c) + me%h_val(k)*v(i)
        end do
        return
    end if

    k = me%n_history
    theta = 1.0_wp/me%gamma
    hv = (theta + me%shift)*v
    if (k == 0) return

    if (.not. me%mid_valid) call factor_middle_matrix(me)
    if (.not. me%mid_ok) return  ! singular middle matrix: fall back to the (safe) initial scaling

    if (me%use_sr1) then

        ! B*v = theta*v + Psi * M^{-1} * Psi^T v, with Psi = Y - theta*S
        allocate(w(k))
        do i = 1, k
            c = pair_col(me, i)
            w(i) = dot_product(me%y(:,c), v) - theta*dot_product(me%s(:,c), v)
        end do
        call lu_solve(me%mid_lu, me%mid_piv, w)
        do i = 1, k
            c = pair_col(me, i)
            hv = hv + (me%y(:,c) - theta*me%s(:,c))*w(i)
        end do

    else

        ! B*v = theta*v - [theta*S Y] * M^{-1} * [theta*S^T v; Y^T v]
        allocate(w(2*k))
        do i = 1, k
            c = pair_col(me, i)
            w(i)   = theta*dot_product(me%s(:,c), v)
            w(k+i) = dot_product(me%y(:,c), v)
        end do
        call lu_solve(me%mid_lu, me%mid_piv, w)
        do i = 1, k
            c = pair_col(me, i)
            hv = hv - theta*me%s(:,c)*w(i) - me%y(:,c)*w(k+i)
        end do

    end if

    end subroutine hessian_vector_product
!*******************************************************************************

!*******************************************************************************
!>
!  the diagonal of the Hessian approximation \( H \), from its compact
!  representation (at a cost of \( O(n k^2) \) for `k` pairs, without
!  forming \( H \)): \( H_{ii} = \theta \mp \psi_i^T M^{-1} \psi_i \), with
!  \( \psi_i \) the `i`-th row of \( \Psi = [\theta S \; Y] \) (BFGS, minus)
!  or of \( \Psi = Y - \theta S \) (SR1, plus), i.e. the row sums of
!  \( \Psi \circ (\Psi M^{-1}) \). Used as a (Jacobi) preconditioner. The
!  result is kept until the pairs or the scaling change.

    subroutine hessian_diagonal(me, d)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(out) :: d   !! the diagonal `dimension(n)`

    integer :: i, k, c, nk, j
    real(wp) :: theta
    real(wp), dimension(:,:), allocatable :: psi, minv, pm

    if (me%exact) then
        d = me%shift
        do k = 1, size(me%h_val)
            if (me%h_irow(k) == me%h_icol(k)) d(me%h_irow(k)) = d(me%h_irow(k)) + me%h_val(k)
        end do
        return
    end if

    k = me%n_history
    theta = 1.0_wp/me%gamma
    d = theta + me%shift
    if (k == 0) return
    if (.not. me%mid_valid) call factor_middle_matrix(me)
    if (.not. me%mid_ok) return
    if (me%diag_valid) then
        d = me%diag + me%shift
        return
    end if

    nk = merge(k, 2*k, me%use_sr1)
    allocate(psi(me%n, nk), minv(nk, nk))
    do i = 1, k
        c = pair_col(me, i)
        if (me%use_sr1) then
            psi(:,i) = me%y(:,c) - theta*me%s(:,c)
        else
            psi(:,i)   = theta*me%s(:,c)
            psi(:,k+i) = me%y(:,c)
        end if
    end do
    minv = 0.0_wp
    do j = 1, nk
        minv(j,j) = 1.0_wp
        call lu_solve(me%mid_lu, me%mid_piv, minv(:,j))
    end do
    pm = matmul(psi, minv)
    me%diag = theta + merge(1.0_wp, -1.0_wp, me%use_sr1)*sum(psi*pm, dim=2)
    me%diag_valid = .true.
    d = me%diag + me%shift

    end subroutine hessian_diagonal
!*******************************************************************************

!*******************************************************************************
!>
!  the Hessian approximation \( H \) as a dense matrix, in the caller's
!  array (for the dense QP solver, the only place a dense \( n \times n \)
!  matrix is used).
!
!  For BFGS, \( H \) is built by applying the stored updates, oldest
!  first, to \( H_0 = \theta I \):
!
!  $$ H \leftarrow H - \frac{(Hs)(Hs)^T}{s^THs} + \frac{yy^T}{y^Ts} $$
!
!  which is the matrix of the compact representation that
!  [[hessian_vector_product]] uses (Byrd, Nocedal & Schnabel, 1994), at
!  \( 2n^2 \) multiplications per pair (only the upper triangle is
!  updated). Forming it from `n` products instead costs about
!  \( 4nk + 4k^2 \) each, i.e. \( 4n^2 + 4nk \) per pair, which is 4 times
!  as much with `k = n` pairs. A pair is skipped if \( s^THs \) isn't
!  positive, which can only happen through round-off: every stored pair has
!  \( y^Ts > 0 \), so each update keeps \( H \) positive definite. Where
!  [[hessian_vector_product]] falls back to \( \theta I \) (it takes the
!  middle matrix to be singular), so does this, so that the QP's matrix is
!  the one the products elsewhere in the iteration are made with.
!
!  For SR1 (whose updates needn't all be defined one at a time) \( H \) is
!  formed from the products with the unit vectors (and symmetrized), and
!  in exact mode from the nonzeros and the shift.

    subroutine hessian_dense(me, h)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:,:),   intent(out)   :: h   !! the matrix `dimension(n,n)` (symmetric, both triangles set)

    integer :: i, j, k, c, n
    real(wp) :: sbs, a, b
    real(wp), dimension(:), allocatable :: bs, e

    n = me%n

    if (me%exact) then
        h = 0.0_wp
        do j = 1, n
            h(j,j) = me%shift
        end do
        do k = 1, size(me%h_val)
            i = me%h_irow(k)
            c = me%h_icol(k)
            h(i,c) = h(i,c) + me%h_val(k)
            if (i /= c) h(c,i) = h(c,i) + me%h_val(k)
        end do
        return
    end if

    if (me%use_sr1) then
        allocate(e(n))
        do j = 1, n
            e = 0.0_wp
            e(j) = 1.0_wp
            call hessian_vector_product(me, e, h(:,j))
        end do
        h = 0.5_wp*(h + transpose(h))   ! (symmetrize away the products' round-off asymmetry)
        return
    end if

    h = 0.0_wp
    do j = 1, n
        h(j,j) = 1.0_wp/me%gamma
    end do
    if (me%n_history > 0 .and. .not. me%mid_valid) call factor_middle_matrix(me)
    if (me%n_history == 0 .or. .not. me%mid_ok) then
        ! (as [[hessian_vector_product]]: the initial scaling, and the shift)
        do j = 1, n
            h(j,j) = h(j,j) + me%shift
        end do
        return
    end if
    allocate(bs(n))
    do k = 1, me%n_history
        c = pair_col(me, k)
        ! bs = H*s, from the upper triangle of H:
        bs = 0.0_wp
        do j = 1, n
            bs(1:j-1) = bs(1:j-1) + h(1:j-1,j)*me%s(j,c)
            bs(j)     = bs(j) + dot_product(h(1:j,j), me%s(1:j,c))
        end do
        sbs = dot_product(me%s(:,c), bs)
        if (.not. sbs > 0.0_wp) cycle
        do j = 1, n
            a = bs(j)/sbs
            b = me%rho(c)*me%y(j,c)
            h(1:j,j) = h(1:j,j) - a*bs(1:j) + b*me%y(1:j,c)
        end do
    end do
    do j = 2, n
        h(j,1:j-1) = h(1:j-1,j)
    end do
    do j = 1, n
        h(j,j) = h(j,j) + me%shift
    end do

    end subroutine hessian_dense
!*******************************************************************************

!*******************************************************************************
!>
!  form and LU-factor the middle matrix of the compact representation for
!  the current pairs and `gamma` (\( \theta = 1/\gamma \)):
!
!  * BFGS: \( M = \begin{bmatrix} \theta S^TS & L \\ L^T & -D \end{bmatrix} \)
!  * SR1: \( M = D + L + L^T - \theta S^TS \)
!
!  with \( L_{pq} = s_p^Ty_q \) for \( p>q \) and \( D = \text{diag}(s_p^Ty_p) \),
!  from the stored inner products (`O(k^2)`, plus `O(k^3)` for the LU).

    subroutine factor_middle_matrix(me)

    class(sqpopt_hessian_type), intent(inout) :: me

    integer :: nk

    nk = merge(me%n_history, 2*me%n_history, me%use_sr1)
    if (allocated(me%mid_lu))  deallocate(me%mid_lu)
    if (allocated(me%mid_piv)) deallocate(me%mid_piv)
    allocate(me%mid_lu(nk,nk), me%mid_piv(nk))
    call middle_matrix(me, me%mid_lu)

    call lu_factor(me%mid_lu, me%mid_piv, me%mid_ok)
    me%mid_valid  = .true.
    me%diag_valid = .false.

    end subroutine factor_middle_matrix
!*******************************************************************************

!*******************************************************************************
!>
!  form the middle matrix \( M \) of the compact representation (see
!  [[factor_middle_matrix]]) for the current pairs and `gamma`, in `a`.

    pure subroutine middle_matrix(me, a)

    class(sqpopt_hessian_type), intent(in)  :: me
    real(wp), dimension(:,:),   intent(out) :: a !! the matrix: `dimension(k,k)` (SR1) or `dimension(2k,2k)` (BFGS)

    integer :: k, p, q, cp, cq
    real(wp) :: theta

    k = me%n_history
    theta = 1.0_wp/me%gamma

    if (me%use_sr1) then
        do p = 1, k
            cp = pair_col(me, p)
            do q = 1, k
                cq = pair_col(me, q)
                ! s_max(p,q)^T y_min(p,q) covers D (p==q) and L + L^T (p/=q):
                a(p,q) = me%sy(pair_col(me, max(p,q)), pair_col(me, min(p,q))) - theta*me%ss(cp,cq)
            end do
        end do
    else
        a = 0.0_wp
        do p = 1, k
            cp = pair_col(me, p)
            do q = 1, k
                cq = pair_col(me, q)
                a(p,q) = theta*me%ss(cp,cq)                         ! theta*S^T S
            end do
            do q = 1, p-1
                cq = pair_col(me, q)
                a(p,k+q) = me%sy(cp,cq)                              ! L (upper-right block)
                a(k+q,p) = a(p,k+q)                                  ! L^T (lower-left block)
            end do
            a(k+p,k+p) = -me%sy(cp,cp)                               ! -D
        end do
    end if

    end subroutine middle_matrix
!*******************************************************************************

!*******************************************************************************
!>
!  the rank `r` of the low-rank part of the compact representation (see the
!  module documentation): `2k` (BFGS) or `k` (SR1) for `k` stored pairs. It
!  is `0` when the matrix is just \( \theta I \) (no pairs yet, or a
!  singular middle matrix, where [[hessian_vector_product]] falls back on
!  that too), and in exact mode.

    integer function hessian_low_rank_size(me) result(r)

    class(sqpopt_hessian_type), intent(inout) :: me

    r = 0
    if (me%exact .or. me%n_history == 0) return
    if (.not. me%mid_valid) call factor_middle_matrix(me)
    if (.not. me%mid_ok) return
    r = merge(me%n_history, 2*me%n_history, me%use_sr1)

    end function hessian_low_rank_size
!*******************************************************************************

!*******************************************************************************
!>
!  the product \( v = U w \) with the low-rank factor of the compact
!  representation (see the module documentation).

    pure subroutine hessian_low_rank_multiply(me, w, v)

    class(sqpopt_hessian_type), intent(in)  :: me
    real(wp), dimension(:),     intent(in)  :: w !! the coefficients `dimension(r)`
    real(wp), dimension(:),     intent(out) :: v !! the product `dimension(n)`

    integer :: i, k, c
    real(wp) :: theta

    k = me%n_history
    theta = 1.0_wp/me%gamma
    v = 0.0_wp
    do i = 1, k
        c = pair_col(me, i)
        if (me%use_sr1) then
            v = v + (me%y(:,c) - theta*me%s(:,c))*w(i)
        else
            v = v + theta*me%s(:,c)*w(i) + me%y(:,c)*w(k+i)
        end if
    end do

    end subroutine hessian_low_rank_multiply
!*******************************************************************************

!*******************************************************************************
!>
!  column `j` of the low-rank factor \( U \) of the compact representation
!  (see the module documentation).

    pure subroutine hessian_low_rank_column(me, j, u)

    class(sqpopt_hessian_type), intent(in)  :: me
    integer,                    intent(in)  :: j !! the column (`1..r`)
    real(wp), dimension(:),     intent(out) :: u !! the column `dimension(n)`

    integer :: k
    real(wp) :: theta

    k = me%n_history
    theta = 1.0_wp/me%gamma
    if (me%use_sr1) then
        u = me%y(:,pair_col(me, j)) - theta*me%s(:,pair_col(me, j))
    else if (j <= k) then
        u = theta*me%s(:,pair_col(me, j))
    else
        u = me%y(:,pair_col(me, j-k))
    end if

    end subroutine hessian_low_rank_column
!*******************************************************************************

!*******************************************************************************
!>
!  the product \( w = U^T v \) with the low-rank factor of the compact
!  representation (see the module documentation).

    pure subroutine hessian_low_rank_multiply_transpose(me, v, w)

    class(sqpopt_hessian_type), intent(in)  :: me
    real(wp), dimension(:),     intent(in)  :: v !! the vector `dimension(n)`
    real(wp), dimension(:),     intent(out) :: w !! the product `dimension(r)`

    integer :: i, k, c
    real(wp) :: theta

    k = me%n_history
    theta = 1.0_wp/me%gamma
    do i = 1, k
        c = pair_col(me, i)
        if (me%use_sr1) then
            w(i) = dot_product(me%y(:,c), v) - theta*dot_product(me%s(:,c), v)
        else
            w(i)   = theta*dot_product(me%s(:,c), v)
            w(k+i) = dot_product(me%y(:,c), v)
        end if
    end do

    end subroutine hessian_low_rank_multiply_transpose
!*******************************************************************************

!*******************************************************************************
!>
!  the middle matrix \( M \) of the compact representation, and the sign
!  \( \sigma \) of its term (see the module documentation).

    pure subroutine hessian_low_rank_middle(me, a, sigma)

    class(sqpopt_hessian_type), intent(in)  :: me
    real(wp), dimension(:,:),   intent(out) :: a     !! the matrix \( M \) `dimension(r,r)`
    real(wp),                   intent(out) :: sigma !! `-1` (BFGS) or `+1` (SR1)

    call middle_matrix(me, a)
    sigma = merge(1.0_wp, -1.0_wp, me%use_sr1)

    end subroutine hessian_low_rank_middle
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
    integer :: c !! column of the `i`-th pair
    real(wp) :: beta

    if (me%use_sr1 .or. me%exact) then
        call hessian_cg_solve(me, v, d)
        return
    end if

    ! standard L-BFGS two-loop recursion:
    q = v
    do i = me%n_history, 1, -1
        c = pair_col(me, i)
        alpha(i) = me%rho(c)*dot_product(me%s(:,c), q)
        q = q - alpha(i)*me%y(:,c)
    end do

    d = me%gamma*q

    do i = 1, me%n_history
        c = pair_col(me, i)
        beta = me%rho(c)*dot_product(me%y(:,c), d)
        d = d + me%s(:,c)*(alpha(i) - beta)
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
    real(wp), dimension(:), intent(in)  :: v !! right-hand side `dimension(n)`
    real(wp), dimension(:), intent(out) :: d !! solution of \( B d = v \) `dimension(n)`

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
!  reset the Hessian approximation, discarding all stored `(s,y)` pairs
!  (and any shift). In exact mode, increase the shift instead (see the
!  module documentation).

    subroutine hessian_reset(me)

    class(sqpopt_hessian_type), intent(inout) :: me

    if (me%exact) then
        ! increase the shift (primal inertia correction, see the module documentation):
        me%shift = min(max(me%shift_min*hessian_size(me), 10.0_wp*me%shift), me%shift_max*hessian_size(me))
        return
    end if

    me%n_history = 0
    me%first     = 1
    me%gamma     = me%gamma0
    me%mid_valid = .false.
    me%shift     = 0.0_wp

    end subroutine hessian_reset
!*******************************************************************************

!*******************************************************************************
!>
!  switch to the exact mode (see the module documentation), with the
!  sparsity pattern `irow`/`icol` of the Hessian of the Lagrangian (each
!  off-diagonal element given once, see [[set_hessian_sparsity]]). The
!  values are set by [[hessian_set_values]]; until then they are zero.

    subroutine hessian_set_exact(me, irow, icol)

    class(sqpopt_hessian_type), intent(inout) :: me
    integer, dimension(:), intent(in) :: irow !! row indices `dimension(nnz)`
    integer, dimension(:), intent(in) :: icol !! column indices `dimension(nnz)`

    me%exact  = .true.
    me%h_irow = irow
    me%h_icol = icol
    if (allocated(me%h_val)) deallocate(me%h_val)
    allocate(me%h_val(size(irow)))
    me%h_val = 0.0_wp
    me%shift = 0.0_wp

    end subroutine hessian_set_exact
!*******************************************************************************

!*******************************************************************************
!>
!  set the values of the exact Hessian for a new major iteration, and, if
!  `decay`, decrease the shift (divided by 3, and set to 0 once below its
!  minimum).

    subroutine hessian_set_values(me, val, decay)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: val   !! the nonzero values `dimension(nnz)`
    logical,                intent(in) :: decay !! whether to decrease the shift (after a good step)

    me%h_val = val
    if (decay) me%shift = me%shift/3.0_wp
    if (me%shift < me%shift_min*hessian_size(me)) me%shift = 0.0_wp

    end subroutine hessian_set_values
!*******************************************************************************

!*******************************************************************************
!>
!  the size of the Hessian, which the shift is relative to:
!  \( \max(1, \max_{ij} |H_{ij}|) \) for the exact Hessian, and
!  \( \max(1, \theta) \) for a quasi-Newton one.

    pure function hessian_size(me) result(h)

    class(sqpopt_hessian_type), intent(in) :: me
    real(wp) :: h

    h = 1.0_wp
    if (me%exact) then
        if (size(me%h_val) > 0) h = max(1.0_wp, maxval(abs(me%h_val)))
    else
        h = max(1.0_wp, 1.0_wp/me%gamma)
    end if

    end function hessian_size
!*******************************************************************************

    end module sqpopt_hessian_module
!*******************************************************************************
