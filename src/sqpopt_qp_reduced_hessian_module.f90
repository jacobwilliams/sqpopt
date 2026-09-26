!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  **Sparse** active-set QP solver for the linearized SQP subproblem (the
!  default for larger problems, see `sqpopt_qp_auto`). Like
!  [[sqpopt_qp_dense_module]], this enforces the linearized constraints
!  and bounds exactly, but stays fully sparse/matrix-free: instead of
!  forming a dense Hessian and an orthonormal null-space basis `Z`, it
!  gets any null-space projection it needs by solving a least-squares
!  problem with `LSQR`, and solves the reduced-space Newton system with
!  **projected conjugate gradients** (Gould, Hribar & Nocedal 1998; Nocedal
!  & Wright, *Numerical Optimization*, Ch. 16) instead of a direct
!  factorization.
!
!  The formulation is the same as the dense solver's: general constraints
!  and bounds are uniform two-sided "rows" (stored in compressed-row form,
!  so each row's product costs only its own nonzeros), and every general
!  constraint violated at `p=0` gets an elastic slack with an \( \ell_1 \)
!  penalty, so `p=0` plus those slacks is a feasible starting point and
!  inconsistent linearized constraints are detected (the penalty weight is
!  raised up to `elastic_weight_max`, then `istat=sqpopt_infeasible`); see
!  [[sqpopt_qp_dense_module]] for the details (here, the slacks also get a
!  small proximal curvature, so CG steps along them stay bounded). Projected
!  CG stops on
!  (relative) convergence, and follows any direction of nonpositive
!  curvature (e.g. along an elastic slack, or an indefinite SR1 Hessian) to
!  the nearest blocking row. Rows are only added to the initial working set
!  if they are (numerically) linearly independent of it.

    module sqpopt_qp_reduced_hessian_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_qp_solve_failed, &
                                     sqpopt_infeasible, sqpopt_infinity
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use lsqr_module,           only: lsqr_solver_ez

    implicit none

    private

    type, public :: sqpopt_reduced_hessian_qp_type
        !! options for the sparse (projected-CG) active-set QP solver.

        integer  :: max_iter     = 100       !! minimum limit on the number of active-set iterations per QP solve
                                             !! (the actual limit is `max(max_iter, 10*(number of rows+1))`)
        integer  :: max_pcg_iter = 0         !! maximum projected-CG iterations per active-set face
                                             !! (`<=0` => twice the number of unknowns)
        real(wp) :: active_tol   = 1.0e-8_wp !! relative tolerance for a row being at a bound, and for the sign
                                             !! of a multiplier (relative to the largest multiplier)
        real(wp) :: opt_tol      = 1.0e-10_wp !! relative tolerance on the projected-gradient stationarity test
                                              !! (relative to \( 1+\lVert Hp+g \rVert_\infty \))
        real(wp) :: pcg_rtol     = 1.0e-10_wp !! projected CG stops once the projected residual has been reduced
                                              !! by this factor (or meets `opt_tol`)
        real(wp) :: feas_tol     = 1.0e-6_wp  !! an elastic slack larger than `feas_tol*max(1,|bound|)` at the
                                              !! solution counts as a violated linearized constraint
        real(wp) :: elastic_weight     = 1.0e4_wp  !! initial elastic penalty weight, relative to
                                                   !! \( \max(1,\lVert g \rVert_\infty) \)
        real(wp) :: elastic_weight_max = 1.0e8_wp  !! largest elastic penalty weight tried (same scaling) before the
                                                   !! linearized constraints are declared inconsistent (lower than
                                                   !! the dense solver's, since the iterative projections' accuracy
                                                   !! is relative to the penalty weight)

        ! `LSQR` settings, used for every null-space projection and multiplier
        ! solve in this module (see [[lsqr_module]] for the precise meaning of
        ! each). `0` for `lsqr_atol`/`lsqr_btol`/`lsqr_conlim` means "let LSQR
        ! use its own machine-precision-based default"; loosening these trades
        ! QP-solve accuracy for speed:
        real(wp) :: lsqr_atol   = 0.0_wp !! `LSQR` relative error tolerance in `A` (0 => `LSQR` default)
        real(wp) :: lsqr_btol   = 0.0_wp !! `LSQR` relative error tolerance in `b` (0 => `LSQR` default)
        real(wp) :: lsqr_conlim = 0.0_wp !! `LSQR` upper limit on `cond(Abar)` (0 => `LSQR` default)
        integer  :: lsqr_itnlim = 0      !! `LSQR` maximum iterations per solve (`<=0` => `2*(rows+columns)+10`)

        contains

        procedure, public :: solve => solve_reduced_hessian_qp

    end type sqpopt_reduced_hessian_qp_type

    type :: csr_rows
        !! the combined constraint rows, in compressed-row form
        integer :: nrows = 0
        integer :: ncols = 0
        integer,  dimension(:), allocatable :: ptr  !! row `k` is entries `ptr(k):ptr(k+1)-1`
        integer,  dimension(:), allocatable :: col
        real(wp), dimension(:), allocatable :: val
    end type csr_rows

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda` (see the module-level
!  documentation). `istat` is as for [[sqpopt_qp_dense_module]]'s solver.

    subroutine solve_reduced_hessian_qp(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    class(sqpopt_reduced_hessian_qp_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! matrix-free Hessian approximation
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! sparse constraint Jacobian, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x       !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g       !! objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c       !! constraint values `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: x_lb, x_ub !! variable bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c_lb, c_ub !! constraint bounds `dimension(m)`
    real(wp), dimension(:),     intent(out)   :: p       !! search direction `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: lambda  !! Lagrange multiplier estimate `dimension(m)`
    integer,                    intent(out)   :: istat   !! status code (see [[sqpopt_types_module]])

    integer :: n, m, nv, nt, mtot, k, i, it, maxit, max_pcg, itnlim
    type(csr_rows) :: rows
    type(sqpopt_sparse_matrix) :: ja   !! the working set's general rows, restricted to the free unknowns
    real(wp), dimension(:), allocatable :: row_lb, row_ub, u, hu_g, gproj, d_total, d_extra, coeff, s_sign, s0
    integer,  dimension(:), allocatable :: status, orig_idx, coeff_idx, slack_row
    logical,  dimension(:), allocatable :: is_equality, fixed
    real(wp) :: rho, rho_max, gscale, alpha, alpha_cap, scale
    integer  :: n_active, blocking, blocking_side
    logical  :: at_face_optimum, truncated

    n = size(g)
    m = size(c)

    ! ---- elastic slacks: one for each general row violated at p=0 ----
    allocate(row_lb(m), row_ub(m), s_sign(m), slack_row(m))
    row_lb = c_lb - c
    row_ub = c_ub - c
    nv = 0
    do i = 1, m
        s_sign(i) = 0.0_wp
        if (row_lb(i) > me%feas_tol*max(1.0_wp, abs(row_lb(i)))) then
            s_sign(i) = 1.0_wp
        else if (row_ub(i) < -me%feas_tol*max(1.0_wp, abs(row_ub(i)))) then
            s_sign(i) = -1.0_wp
        end if
        if (s_sign(i) /= 0.0_wp) then
            nv = nv + 1
            slack_row(nv) = i
        end if
    end do
    nt   = n + nv
    mtot = m + nt
    max_pcg = merge(me%max_pcg_iter, 2*nt, me%max_pcg_iter > 0)
    itnlim  = merge(me%lsqr_itnlim, 2*(mtot+nt)+10, me%lsqr_itnlim > 0)

    ! ---- the combined constraint rows, in compressed-row form ----
    call build_rows()
    row_lb = [row_lb, x_lb - x, spread(0.0_wp, 1, nv)]
    row_ub = [row_ub, x_ub - x, spread(sqpopt_infinity, 1, nv)]
    allocate(is_equality(mtot))
    do k = 1, mtot
        is_equality(k) = row_ub(k)-row_lb(k) <= me%active_tol*max(1.0_wp, abs(row_lb(k)))
    end do

    gscale  = 1.0_wp
    if (n > 0) gscale = max(1.0_wp, maxval(abs(g)))
    rho     = me%elastic_weight*gscale
    rho_max = me%elastic_weight_max*gscale

    ! ---- feasible starting point: p=0, slacks just large enough ----
    allocate(u(nt)); u = 0.0_wp
    u(1:n) = min(max(0.0_wp, row_lb(m+1:m+n)), row_ub(m+1:m+n))
    do k = 1, nv
        i = slack_row(k)
        u(n+k) = merge(row_lb(i), -row_ub(i), s_sign(i) > 0.0_wp)
    end do
    s0 = u(n+1:nt)

    ! ---- initial working set ----
    allocate(status(mtot)); status = 0
    call initial_working_set()

    ! ---- active-set iterations ----
    istat = sqpopt_qp_solve_failed
    allocate(coeff(0), coeff_idx(0), gproj(nt), d_total(nt), d_extra(nt))
    maxit = max(me%max_iter, 10*(mtot+1))

    do it = 1, maxit

        call build_working_set(rows, status, m, ja, fixed, orig_idx, n_active)

        hu_g = gradient(u)
        ! (relative to the gradient on every free unknown, including the
        ! elastic slacks, whose large penalty weight limits how accurately
        ! the iterative projections can be computed):
        scale = 1.0_wp + maxval(abs(merge(0.0_wp, hu_g, fixed)))
        call project_null(ja, fixed, hu_g, gproj)
        at_face_optimum = .false.

        if (norm2(gproj) <= me%opt_tol*scale) then

            at_face_optimum = .true.

        else

            call projected_cg(ja, fixed, hu_g, gproj, me%opt_tol*scale, d_total, d_extra, truncated)

            ! first, the step accumulated by CG (which may itself be blocked):
            alpha_cap = 1.0_wp
            call ratio_test(u, d_total, alpha_cap, alpha, blocking, blocking_side)
            u = u + alpha*d_total
            if (blocking /= 0 .and. alpha < alpha_cap - 1.0e-12_wp) then
                status(blocking) = blocking_side
                cycle
            end if

            if (truncated) then
                ! then, if CG found nonpositive curvature, move along that
                ! (downhill) direction, not bounded by 1, to the nearest blocking row:
                alpha_cap = huge(1.0_wp)
                call ratio_test(u, d_extra, alpha_cap, alpha, blocking, blocking_side)
                if (blocking == 0) exit  ! unbounded QP
                u = u + alpha*d_extra
                status(blocking) = blocking_side
                cycle
            end if
            at_face_optimum = .true.

        end if

        if (at_face_optimum) then

            ! multipliers for the working set, from H*u+g = G^T*lambda_G + (bound
            ! multipliers on the fixed unknowns): lambda_G is the least-squares
            ! solution on the free unknowns, and each bound multiplier is then
            ! the remaining residual in its fixed coordinate:
            hu_g = gradient(u)
            coeff     = [real(wp) ::]
            coeff_idx = [integer ::]
            if (n_active > 0) then
                block
                    real(wp), dimension(ja%nrows) :: lam_g
                    real(wp), dimension(nt) :: resid
                    type(lsqr_solver_ez) :: lsqr
                    integer :: istop, idx, j
                    lam_g = 0.0_wp
                    if (ja%nrows > 0) then
                        call lsqr%initialize(nt, ja%nrows, ja%val, ja%icol, ja%irow, &
                                              atol=me%lsqr_atol, btol=me%lsqr_btol, conlim=me%lsqr_conlim, itnlim=itnlim)
                        call lsqr%solve(merge(0.0_wp, hu_g, fixed), 0.0_wp, lam_g, istop)
                    end if
                    resid = hu_g
                    do idx = 1, ja%nrows  ! (the first ja%nrows entries of orig_idx are the general rows)
                        k = orig_idx(idx)
                        do j = rows%ptr(k), rows%ptr(k+1)-1
                            resid(rows%col(j)) = resid(rows%col(j)) - rows%val(j)*lam_g(idx)
                        end do
                    end do
                    deallocate(coeff)
                    allocate(coeff(n_active))
                    coeff(1:ja%nrows) = lam_g
                    do idx = ja%nrows+1, n_active
                        coeff(idx) = resid(orig_idx(idx) - m)
                    end do
                    coeff_idx = orig_idx
                end block
            end if

            ! drop the inequality row with the most wrongly-signed multiplier, if any:
            block
                integer  :: worst_idx, idx
                real(wp) :: worst, tol_mult
                ! (relative to the largest multiplier, but not counting the elastic
                ! slacks' bounds, whose multipliers are the large penalty weight):
                tol_mult = me%active_tol
                if (size(coeff) > 0) tol_mult = me%active_tol*max(1.0_wp, &
                    maxval(abs(coeff), mask=coeff_idx <= m+n))
                worst     = tol_mult
                worst_idx = 0
                do idx = 1, n_active
                    k = orig_idx(idx)
                    if (is_equality(k)) cycle
                    if (status(k) == -1 .and. -coeff(idx) > worst) then
                        worst = -coeff(idx); worst_idx = k
                    else if (status(k) == 1 .and. coeff(idx) > worst) then
                        worst = coeff(idx); worst_idx = k
                    end if
                end do
                if (worst_idx /= 0) then
                    status(worst_idx) = 0
                    cycle
                end if
            end block

            ! optimal for the current elastic weight. Any slack still positive?
            if (any(u(n+1:nt) > me%feas_tol*max(1.0_wp, s0))) then
                if (rho < rho_max) then
                    rho = min(100.0_wp*rho, rho_max)
                    cycle
                end if
                istat = sqpopt_infeasible
            else
                istat = sqpopt_success
            end if
            exit

        end if

    end do

    p = u(1:n)
    lambda = 0.0_wp
    do k = 1, size(coeff_idx)
        if (coeff_idx(k) <= m) lambda(coeff_idx(k)) = coeff(k)
    end do

    contains

        subroutine build_rows()
        !! the combined rows: `J` (plus the slack columns), then a unit row per unknown
        integer, dimension(:), allocatable :: cnt, pos
        integer :: kk, r
        rows%nrows = mtot
        rows%ncols = nt
        allocate(cnt(mtot)); cnt = 1           ! (every general row has room for its slack; unused if none)
        cnt(1:m) = 0
        do kk = 1, jac%nnz
            cnt(jac%irow(kk)) = cnt(jac%irow(kk)) + 1
        end do
        do kk = 1, nv
            cnt(slack_row(kk)) = cnt(slack_row(kk)) + 1
        end do
        allocate(rows%ptr(mtot+1))
        rows%ptr(1) = 1
        do r = 1, mtot
            rows%ptr(r+1) = rows%ptr(r) + cnt(r)
        end do
        allocate(rows%col(rows%ptr(mtot+1)-1), rows%val(rows%ptr(mtot+1)-1))
        pos = rows%ptr(1:mtot)
        do kk = 1, jac%nnz
            r = jac%irow(kk)
            rows%col(pos(r)) = jac%icol(kk); rows%val(pos(r)) = jac%val(kk); pos(r) = pos(r) + 1
        end do
        do kk = 1, nv
            r = slack_row(kk)
            rows%col(pos(r)) = n+kk; rows%val(pos(r)) = s_sign(r); pos(r) = pos(r) + 1
        end do
        do kk = 1, nt
            r = m + kk
            rows%col(pos(r)) = kk; rows%val(pos(r)) = 1.0_wp; pos(r) = pos(r) + 1
        end do
        end subroutine build_rows

        function gradient(v) result(gr)
        !! the gradient of the (elastic) QP objective at `v`: `H*v_p + g`, then `rho + delta*s` for each slack
        real(wp), dimension(:), intent(in) :: v
        real(wp), dimension(size(v)) :: gr
        call hext_product(v, gr)
        gr(1:n)    = gr(1:n) + g
        gr(n+1:nt) = gr(n+1:nt) + rho
        end function gradient

        subroutine hext_product(v, hv)
        !! the (elastic) QP Hessian times `v`: `H` on `p`, and a small proximal
        !! curvature `delta` on the slacks. (With zero curvature there, CG
        !! steps along a slack would be unboundedly long; the slacks are zero
        !! at any feasible solution, so `delta` doesn't change the solution
        !! then, and only perturbs the size, not the positivity, of the
        !! slacks for inconsistent constraints.)
        real(wp), dimension(:), intent(in)  :: v
        real(wp), dimension(:), intent(out) :: hv
        call hessian%hv_product(v(1:n), hv(1:n))
        hv(n+1:nt) = gscale*v(n+1:nt)
        end subroutine hext_product

        subroutine project_null(ja, fixed, v, out)
        !! project `v` onto the null space of the working set: zero in the
        !! `fixed` coordinates (exactly), then `out = v - ja^T z` on the rest,
        !! with `z` the minimum-norm least-squares solution of `ja^T z ~ v`
        !! (`ja` has no entries in the fixed columns)
        type(sqpopt_sparse_matrix), intent(in)  :: ja
        logical,  dimension(:),     intent(in)  :: fixed
        real(wp), dimension(:),     intent(in)  :: v
        real(wp), dimension(:),     intent(out) :: out
        type(lsqr_solver_ez) :: lsqr
        real(wp), dimension(ja%nrows) :: z
        integer :: istop, kk
        out = merge(0.0_wp, v, fixed)
        if (ja%nrows == 0) return
        call lsqr%initialize(nt, ja%nrows, ja%val, ja%icol, ja%irow, &
                              atol=me%lsqr_atol, btol=me%lsqr_btol, conlim=me%lsqr_conlim, itnlim=itnlim)
        call lsqr%solve(out, 0.0_wp, z, istop)
        do kk = 1, ja%nnz
            out(ja%icol(kk)) = out(ja%icol(kk)) - ja%val(kk)*z(ja%irow(kk))
        end do
        end subroutine project_null

        subroutine projected_cg(ja, fixed, hu_g0, gproj0, abs_tol, d_total, d_extra, truncated)
        !! projected conjugate gradients on the current face, from `u`: returns
        !! the accumulated step `d_total`; if a direction of nonpositive
        !! curvature is found, it is returned (oriented downhill) in `d_extra`
        !! with `truncated=.true.`
        type(sqpopt_sparse_matrix), intent(in)  :: ja
        logical,  dimension(:),     intent(in)  :: fixed   !! the unknowns fixed at a bound by the working set
        real(wp), dimension(:),     intent(in)  :: hu_g0   !! `H*u+g` at `u`
        real(wp), dimension(:),     intent(in)  :: gproj0  !! its projection onto the face
        real(wp),                   intent(in)  :: abs_tol !! absolute stopping tolerance on the projected residual
        real(wp), dimension(:),     intent(out) :: d_total, d_extra
        logical,                    intent(out) :: truncated
        real(wp), dimension(nt) :: r, gp, dvec, hd, tmp
        real(wp) :: rg_old, rg_new, kappa, alpha, beta, tol
        integer :: j, max_it
        d_total = 0.0_wp
        d_extra = 0.0_wp
        truncated = .false.
        r    = hu_g0
        gp   = gproj0
        dvec = -gp
        ! (`r^T P r = |P r|^2` for the orthogonal projector `P`; the latter is
        ! computed much more accurately when `r` is large and `P r` is small)
        rg_old = dot_product(gp, gp)
        tol = max(abs_tol, me%pcg_rtol*norm2(gproj0))
        ! in exact arithmetic CG converges within the dimension of the null
        ! space; iterating further only accumulates rounding error:
        max_it = min(max_pcg, max(1, count(.not. fixed) - ja%nrows))
        do j = 1, max_it
            if (norm2(gp) <= tol) exit
            call hext_product(dvec, hd)
            kappa = dot_product(dvec, hd)
            if (kappa <= 1.0e-10_wp*norm2(dvec)*norm2(hd)) then
                ! (numerically) zero or negative curvature along `dvec`, which is
                ! a descent direction for the model at `u+d_total` (re-projected,
                ! so that following it can't drift off the working set):
                call project_null(ja, fixed, dvec, d_extra)
                truncated = norm2(d_extra) > 1.0e-8_wp*norm2(dvec)
                return
            end if
            alpha   = rg_old/kappa
            d_total = d_total + alpha*dvec
            r       = r + alpha*hd
            call project_null(ja, fixed, r, gp)
            rg_new  = dot_product(gp, gp)
            if (abs(rg_old) <= tiny(1.0_wp)) exit
            beta   = rg_new/rg_old
            dvec   = -gp + beta*dvec
            rg_old = rg_new
        end do
        ! make sure the step stays exactly on the working set:
        call project_null(ja, fixed, d_total, tmp)
        d_total = tmp
        end subroutine projected_cg

        subroutine ratio_test(base, d, alpha_cap, alpha, blocking, blocking_side)
        !! the largest `alpha <= alpha_cap` for which `base+alpha*d` satisfies every
        !! row not in the working set, and the row (and side) that blocks first
        real(wp), dimension(:), intent(in)  :: base, d
        real(wp),               intent(in)  :: alpha_cap
        real(wp),               intent(out) :: alpha
        integer,                intent(out) :: blocking, blocking_side
        real(wp) :: rate, alpha_k, val, dnorm
        integer :: kk
        alpha = alpha_cap
        blocking = 0
        blocking_side = 0
        dnorm = norm2(d)
        do kk = 1, mtot
            if (status(kk) /= 0) cycle
            rate = row_dot(rows, kk, d)
            if (abs(rate) <= 1.0e-12_wp*row_norm(rows, kk)*dnorm) cycle
            val = row_dot(rows, kk, base)
            if (rate > 0.0_wp) then
                if (row_ub(kk) >= sqpopt_infinity) cycle
                alpha_k = max((row_ub(kk) - val)/rate, 0.0_wp)
                if (alpha_k < alpha) then
                    alpha = alpha_k; blocking = kk; blocking_side = 1
                end if
            else
                if (row_lb(kk) <= -sqpopt_infinity) cycle
                alpha_k = max((row_lb(kk) - val)/rate, 0.0_wp)
                if (alpha_k < alpha) then
                    alpha = alpha_k; blocking = kk; blocking_side = -1
                end if
            end if
        end do
        end subroutine ratio_test

        subroutine initial_working_set()
        !! add the rows that are at a bound at `u` (equality rows first) to the
        !! working set, skipping any whose component outside the span of the
        !! rows already added is negligible (i.e., that are linearly dependent)
        type(sqpopt_sparse_matrix) :: ja_cur
        integer, dimension(:), allocatable :: idx_cur
        logical, dimension(:), allocatable :: fixed_cur
        real(wp), dimension(nt) :: a, r
        real(wp) :: val
        integer :: kk, pass, side, na, j
        do pass = 1, 2
            do kk = 1, mtot
                if (status(kk) /= 0) cycle
                if ((pass == 1) .neqv. is_equality(kk)) cycle
                val = row_dot(rows, kk, u)
                if (abs(val-row_lb(kk)) <= me%active_tol*max(1.0_wp, abs(row_lb(kk)))) then
                    side = -1
                else if (abs(val-row_ub(kk)) <= me%active_tol*max(1.0_wp, abs(row_ub(kk)))) then
                    side = 1
                else
                    cycle
                end if
                a = 0.0_wp
                do j = rows%ptr(kk), rows%ptr(kk+1)-1
                    a(rows%col(j)) = rows%val(j)
                end do
                call build_working_set(rows, status, m, ja_cur, fixed_cur, idx_cur, na)
                if (na >= nt) return
                call project_null(ja_cur, fixed_cur, a, r)
                if (norm2(r) > 1.0e-8_wp*norm2(a)) status(kk) = side
            end do
        end do
        end subroutine initial_working_set

    end subroutine solve_reduced_hessian_qp
!*******************************************************************************

!*******************************************************************************
!>
!  dot product of row `k` of `rows` with `v`.

    pure function row_dot(rows, k, v) result(s)

    type(csr_rows),         intent(in) :: rows
    integer,                intent(in) :: k
    real(wp), dimension(:), intent(in) :: v
    real(wp) :: s

    integer :: j

    s = 0.0_wp
    do j = rows%ptr(k), rows%ptr(k+1)-1
        s = s + rows%val(j)*v(rows%col(j))
    end do

    end function row_dot
!*******************************************************************************

!*******************************************************************************
!>
!  2-norm of row `k` of `rows`.

    pure function row_norm(rows, k) result(s)

    type(csr_rows), intent(in) :: rows
    integer,        intent(in) :: k
    real(wp) :: s

    s = norm2(rows%val(rows%ptr(k):rows%ptr(k+1)-1))

    end function row_norm
!*******************************************************************************

!*******************************************************************************
!>
!  the working set, split into its active variable-bound rows -- which
!  simply fix those unknowns (`fixed`) -- and its active general rows,
!  gathered (renumbered `1..ja%nrows`) into the COO matrix `ja` restricted
!  to the *free* unknowns (entries in fixed columns are dropped). Handling
!  the bounds this way, rather than as rows in `ja`, keeps the least-squares
!  problems small and makes the fixed coordinates of every projection
!  exactly zero. `orig_idx` gives the index in `rows` of each active row:
!  the general rows first (in the order of `ja`'s rows), then the bounds.

    subroutine build_working_set(rows, status, m, ja, fixed, orig_idx, n_active)

    type(csr_rows),             intent(in)  :: rows
    integer,  dimension(:),     intent(in)  :: status   !! 0 = inactive, -1/+1 = active at the lower/upper bound
    integer,                    intent(in)  :: m        !! number of general rows (the rest are the bounds on each unknown)
    type(sqpopt_sparse_matrix), intent(out) :: ja       !! the active general rows, on the free unknowns
    logical,  dimension(:), allocatable, intent(out) :: fixed    !! `dimension(ncols)`: unknowns fixed at a bound
    integer,  dimension(:), allocatable, intent(out) :: orig_idx !! index in `rows` of each active row
    integer,                    intent(out) :: n_active !! number of active rows

    integer :: k, j, idx, nnz_a, n_gen

    allocate(fixed(rows%ncols))
    fixed = status(m+1:m+rows%ncols) /= 0

    n_active = count(status /= 0)
    n_gen    = count(status(1:m) /= 0)
    allocate(orig_idx(n_active))
    idx = 0
    do k = 1, m
        if (status(k) /= 0) then
            idx = idx + 1
            orig_idx(idx) = k
        end if
    end do
    do k = m+1, rows%nrows
        if (status(k) /= 0) then
            idx = idx + 1
            orig_idx(idx) = k
        end if
    end do

    nnz_a = 0
    do idx = 1, n_gen
        k = orig_idx(idx)
        do j = rows%ptr(k), rows%ptr(k+1)-1
            if (.not. fixed(rows%col(j))) nnz_a = nnz_a + 1
        end do
    end do
    ja%nrows = n_gen
    ja%ncols = rows%ncols
    ja%nnz   = nnz_a
    allocate(ja%irow(nnz_a), ja%icol(nnz_a), ja%val(nnz_a))
    nnz_a = 0
    do idx = 1, n_gen
        k = orig_idx(idx)
        do j = rows%ptr(k), rows%ptr(k+1)-1
            if (fixed(rows%col(j))) cycle
            nnz_a = nnz_a + 1
            ja%irow(nnz_a) = idx
            ja%icol(nnz_a) = rows%col(j)
            ja%val(nnz_a)  = rows%val(j)
        end do
    end do

    end subroutine build_working_set
!*******************************************************************************

    end module sqpopt_qp_reduced_hessian_module
!*******************************************************************************
