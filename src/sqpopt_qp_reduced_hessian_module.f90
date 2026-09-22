!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Opt-in **sparse** active-set QP solver for the linearized SQP
!  subproblem. Like
!  [[sqpopt_qp_dense_module]], this enforces the linearized constraints
!  and bounds *exactly* (unlike `sqpopt_qp_solver_module`'s default composite
!  step), but stays fully sparse/matrix-free: instead of forming a dense
!  `n x n` Hessian and an `n x (n-m_a)` orthonormal null-space basis
!  (`Z`), it gets any null-space projection it needs by re-solving a
!  small least-squares problem with `LSQR` (the same technique the
!  composite step already uses for its tangential step), and solves the
!  reduced-space Newton system with **projected conjugate gradients**
!  (Gould, Hribar & Nocedal 1998; Nocedal & Wright, *Numerical
!  Optimization*, Ch. 16) instead of a direct Cholesky factorization.
!
!  As in [[sqpopt_qp_dense_module]], general constraints and variable
!  bounds are treated uniformly as `m+n` candidate two-sided "rows" on
!  `p` directly (no `w=(p,s)` slack padding -- mathematically equivalent,
!  simpler): the `m` rows of the sparse Jacobian, followed by `n` unit
!  rows for the variable bounds. Equality rows (`row_lb==row_ub`) are
!  permanently active and never leave the working set. The outer
!  active-set control logic (ratio test to add a row, Lagrange-multiplier
!  sign check to drop one) is otherwise identical to the dense solver --
!  only the linear algebra inside each iteration differs.

    module sqpopt_qp_reduced_hessian_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_qp_solve_failed
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use lsqr_module,           only: lsqr_solver_ez

    implicit none

    private

    type, public :: sqpopt_reduced_hessian_qp_type
        !! options and workspace for the sparse (projected-CG) active-set QP solver.

        integer  :: max_iter     = 100      !! maximum number of active-set changes allowed per QP solve
        integer  :: max_pcg_iter = 0        !! maximum projected-CG iterations per active-set face (`<=0` => use `n`)
        real(wp) :: active_tol   = 1.0e-8_wp !! tolerance used to detect an (in)active/equality row
        real(wp) :: opt_tol      = 1.0e-8_wp !! tolerance on the projected-residual stationarity test

        ! `LSQR` settings, used for every `project_null`/`project_onto_active`/multiplier
        ! solve in this module (see [[lsqr_module]] for the precise meaning of each --
        ! `0` for `lsqr_atol`/`lsqr_btol`/`lsqr_conlim` means "let LSQR use its own
        ! machine-precision-based default", which is tighter than usually necessary and
        ! can mean more internal LSQR iterations per call; loosening these (and/or
        ! raising `lsqr_itnlim`) is the main lever for trading QP-solve accuracy for
        ! speed in this QP mode:
        real(wp) :: lsqr_atol   = 0.0_wp !! `LSQR` relative error tolerance in `A` (0 => `LSQR` default)
        real(wp) :: lsqr_btol   = 0.0_wp !! `LSQR` relative error tolerance in `b` (0 => `LSQR` default)
        real(wp) :: lsqr_conlim = 0.0_wp !! `LSQR` upper limit on `cond(Abar)` (0 => `LSQR` default)
        integer  :: lsqr_itnlim = 100    !! `LSQR` maximum iterations per solve

        contains

        procedure, public :: solve => solve_reduced_hessian_qp

    end type sqpopt_reduced_hessian_qp_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda`, using the sparse
!  projected-CG active-set method (see the module-level documentation).

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

    integer :: n, m, mtot, k, it, n_active, max_pcg
    type(sqpopt_sparse_matrix) :: arows  !! combined m+n rows: J's rows, then n identity (bound) rows
    real(wp), dimension(:), allocatable :: row_lb, row_ub, u, hu_g, gproj, d_total, d_extra, coeff
    logical,  dimension(:), allocatable :: is_equality
    integer,  dimension(:), allocatable :: status
    type(sqpopt_sparse_matrix) :: ja
    integer,  dimension(:), allocatable :: orig_idx
    real(wp), dimension(:), allocatable :: rhs_active
    real(wp) :: alpha, rate, alpha_k, worst, val, alpha_cap
    integer  :: blocking, blocking_side, worst_idx, idx
    logical  :: at_face_optimum, truncated

    n = size(g)
    m = size(c)
    mtot = m + n
    max_pcg = merge(me%max_pcg_iter, n, me%max_pcg_iter > 0)

    ! ---- combine general constraints (rows 1..m) and variable bounds (rows m+1..m+n)
    !      into one sparse set of two-sided "rows" on p ----
    allocate(arows%irow(jac%nnz+n), arows%icol(jac%nnz+n), arows%val(jac%nnz+n))
    arows%nrows = mtot
    arows%ncols = n
    arows%nnz   = jac%nnz + n
    if (jac%nnz > 0) then
        arows%irow(1:jac%nnz) = jac%irow(1:jac%nnz)
        arows%icol(1:jac%nnz) = jac%icol(1:jac%nnz)
        arows%val(1:jac%nnz)  = jac%val(1:jac%nnz)
    end if
    do k = 1, n
        arows%irow(jac%nnz+k) = m+k
        arows%icol(jac%nnz+k) = k
        arows%val(jac%nnz+k)  = 1.0_wp
    end do

    allocate(row_lb(mtot), row_ub(mtot))
    row_lb(1:m) = c_lb - c
    row_ub(1:m) = c_ub - c
    row_lb(m+1:mtot) = x_lb - x
    row_ub(m+1:mtot) = x_ub - x

    allocate(is_equality(mtot))
    is_equality(1:m)      = (c_ub-c_lb) <= me%active_tol
    is_equality(m+1:mtot) = (x_ub-x_lb) <= me%active_tol

    allocate(status(mtot))
    status = 0
    where (is_equality) status = -1  !! permanently active

    allocate(u(n)); u = 0.0_wp

    ! ---- phase 1: bootstrap a feasible-for-the-initial-working-set starting point ----
    call project_onto_active(arows, row_lb, row_ub, status, mtot, n, u, &
                              me%lsqr_atol, me%lsqr_btol, me%lsqr_conlim, me%lsqr_itnlim)
    do k = 1, n
        if (is_equality(m+k)) cycle
        if (u(k) < row_lb(m+k) - me%active_tol) then
            status(m+k) = -1
        else if (u(k) > row_ub(m+k) + me%active_tol) then
            status(m+k) = 1
        end if
    end do
    call sparse_row_value(arows, u, m, is_equality, status, row_lb, row_ub, me%active_tol)
    call project_onto_active(arows, row_lb, row_ub, status, mtot, n, u, &
                              me%lsqr_atol, me%lsqr_btol, me%lsqr_conlim, me%lsqr_itnlim)

    ! ---- phase 2: active-set iterations ----
    istat = sqpopt_qp_solve_failed
    allocate(coeff(0))

    do it = 1, max(me%max_iter, 10*(mtot+1))

        call build_active_set(arows, row_lb, row_ub, status, mtot, n, ja, rhs_active, orig_idx, n_active)

        if (allocated(hu_g)) deallocate(hu_g)
        allocate(hu_g(n))
        call hessian%hv_product(u, hu_g)
        hu_g = hu_g + g

        if (allocated(gproj)) deallocate(gproj)
        allocate(gproj(n))
        call project_null(ja, n_active, n, hu_g, gproj, me%lsqr_atol, me%lsqr_btol, me%lsqr_conlim, me%lsqr_itnlim)

        at_face_optimum = .false.

        if (norm2(gproj) <= me%opt_tol) then

            at_face_optimum = .true.

        else

            if (allocated(d_total)) deallocate(d_total)
            if (allocated(d_extra)) deallocate(d_extra)
            allocate(d_total(n), d_extra(n))
            call projected_cg(hessian, ja, n_active, n, hu_g, gproj, max_pcg, me%opt_tol, me%active_tol, &
                               me%lsqr_atol, me%lsqr_btol, me%lsqr_conlim, me%lsqr_itnlim, &
                               d_total, d_extra, truncated)

            ! ---- ratio test against every currently-inactive row ----
            if (truncated) then
                alpha_cap = huge(1.0_wp)
            else
                alpha_cap = 1.0_wp
            end if
            alpha = alpha_cap
            blocking = 0
            blocking_side = 0
            do k = 1, mtot
                if (status(k) /= 0) cycle
                if (truncated) then
                    rate = sparse_dot_row(arows, k, d_extra)
                else
                    rate = sparse_dot_row(arows, k, d_total)
                end if
                if (rate > me%active_tol) then
                    alpha_k = ratio_alpha(arows, k, u, d_total, d_extra, truncated, row_ub(k), .true.)
                    if (alpha_k < alpha) then
                        alpha = max(alpha_k, 0.0_wp); blocking = k; blocking_side = 1
                    end if
                else if (rate < -me%active_tol) then
                    alpha_k = ratio_alpha(arows, k, u, d_total, d_extra, truncated, row_lb(k), .false.)
                    if (alpha_k < alpha) then
                        alpha = max(alpha_k, 0.0_wp); blocking = k; blocking_side = -1
                    end if
                end if
            end do

            if (truncated) then
                u = u + d_total + alpha*d_extra
            else
                u = u + alpha*d_total
            end if

            if (blocking /= 0 .and. alpha < alpha_cap - 1.0e-10_wp) then
                status(blocking) = blocking_side
                cycle
            else if (truncated) then
                ! unbounded direction with no blocking row: not a well-posed
                ! bounded QP face -- bail out defensively:
                istat = sqpopt_qp_solve_failed
                exit
            else
                at_face_optimum = .true.
            end if

        end if

        if (at_face_optimum) then

            if (n_active == 0) then
                istat = sqpopt_success
                exit
            end if

            if (allocated(hu_g)) deallocate(hu_g)
            allocate(hu_g(n))
            call hessian%hv_product(u, hu_g)
            hu_g = hu_g + g

            block
                real(wp), dimension(n_active) :: coeff_local
                type(lsqr_solver_ez) :: lsqr
                integer :: istop
                call lsqr%initialize(n, n_active, ja%val, ja%icol, ja%irow, &
                                      atol=me%lsqr_atol, btol=me%lsqr_btol, conlim=me%lsqr_conlim, itnlim=me%lsqr_itnlim)
                call lsqr%solve(hu_g, 0.0_wp, coeff_local, istop)
                if (allocated(coeff)) deallocate(coeff)
                allocate(coeff(n_active))
                coeff = coeff_local
            end block

            worst = me%active_tol
            worst_idx = 0
            do idx = 1, n_active
                k = orig_idx(idx)
                if (is_equality(k)) cycle
                if (status(k) == -1) then
                    if (-coeff(idx) > worst) then
                        worst = -coeff(idx); worst_idx = k
                    end if
                else
                    if (coeff(idx) > worst) then
                        worst = coeff(idx); worst_idx = k
                    end if
                end if
            end do

            if (worst_idx == 0) then
                istat = sqpopt_success
                exit
            else
                status(worst_idx) = 0
                cycle
            end if

        end if

    end do

    p = u
    lambda = 0.0_wp
    do idx = 1, size(orig_idx)
        k = orig_idx(idx)
        if (k <= m) lambda(k) = coeff(idx)
    end do

    contains

    !> ratio-test step length for row `k` along the candidate move, expressed
    !! generically for both the "scale the whole PCG step" case
    !! (`truncated=.false.`, tests `u+alpha*d_total`, `alpha` in `[0,1]`) and
    !! the "scale only the truncation direction" case (`truncated=.true.`,
    !! tests `u+d_total+alpha*d_extra`, `alpha>=0`).
    pure function ratio_alpha(arows, k, u, d_total, d_extra, truncated, bound_val, is_upper) result(a)
    type(sqpopt_sparse_matrix), intent(in) :: arows !! sparse combined constraint matrix
    integer,                    intent(in) :: k !! row index in the sparse combined constraint matrix
    real(wp), dimension(:),     intent(in) :: u !! current solution vector
    real(wp), dimension(:),     intent(in) :: d_total !! candidate move direction
    real(wp), dimension(:),     intent(in) :: d_extra !! extra move direction for truncated step
    logical,                    intent(in) :: truncated !! whether the step is truncated
    logical,                    intent(in) :: is_upper !! whether the bound is an upper bound
    real(wp),                   intent(in) :: bound_val !! value of the bound
    real(wp) :: a
    real(wp) :: base_val, rate
    if (truncated) then
        base_val = sparse_dot_row(arows, k, u + d_total)
        rate     = sparse_dot_row(arows, k, d_extra)
    else
        base_val = sparse_dot_row(arows, k, u)
        rate     = sparse_dot_row(arows, k, d_total)
    end if
    a = (bound_val - base_val)/rate
    end function ratio_alpha

    end subroutine solve_reduced_hessian_qp
!*******************************************************************************

!*******************************************************************************
!>
!  dot product of row `k` of the sparse combined row set `arows` with a
!  dense vector `v` (`dimension(n)`), i.e. `arows(k,:) . v`.

    pure function sparse_dot_row(arows, k, v) result(s)

    type(sqpopt_sparse_matrix), intent(in) :: arows !! sparse combined constraint matrix
    integer,                    intent(in) :: k !! row index in the sparse combined constraint matrix
    real(wp), dimension(:),     intent(in) :: v !! dense vector to be dotted with row `k` of `arows`
    real(wp) :: s !! result of the dot product of row `k` of `arows` with vector `v`

    integer :: j

    s = 0.0_wp
    do j = 1, arows%nnz
        if (arows%irow(j) == k) s = s + arows%val(j)*v(arows%icol(j))
    end do

    end function sparse_dot_row
!*******************************************************************************

!*******************************************************************************
!>
!  project a vector `v` onto the null space of the active-row matrix
!  `ja` (`n_active x n`): `out = v - ja^T*z`, `z` the minimum-norm
!  least-squares solution of `ja^T*z ~ v`, solved with `LSQR` using the
!  same transpose-orientation trick as `sqpopt_qp_solver_module`'s
!  composite step (swap `irow`/`icol` so `LSQR` sees `ja^T` directly).

    subroutine project_null(ja, n_active, n, v, out, atol, btol, conlim, itnlim)

    type(sqpopt_sparse_matrix), intent(in)  :: ja !! active-row matrix (`n_active x n`)
    integer,                    intent(in)  :: n_active !! number of active rows in `ja`
    integer,                    intent(in)  :: n !! number of columns in `ja`
    real(wp), dimension(:),     intent(in)  :: v !! vector to be projected onto the null space of `ja`
    real(wp), dimension(:),     intent(out) :: out !! projected vector onto the null space of `ja`
    real(wp),                   intent(in)  :: atol     !! `LSQR` tolerances (see [[sqpopt_reduced_hessian_qp_type]])
    real(wp),                   intent(in)  :: btol     !! `LSQR` tolerances (see [[sqpopt_reduced_hessian_qp_type]])
    real(wp),                   intent(in)  :: conlim   !! `LSQR` tolerances (see [[sqpopt_reduced_hessian_qp_type]])
    integer,                    intent(in)  :: itnlim   !! `LSQR` max iterations

    type(lsqr_solver_ez) :: lsqr
    real(wp), dimension(:), allocatable :: z !! minimum-norm least-squares solution of `ja^T*z ~ v`
    real(wp), dimension(n) :: jtz !! `ja^T*z`
    integer :: istop, j

    if (n_active == 0) then
        out = v
        return
    end if

    allocate(z(n_active))
    call lsqr%initialize(n, n_active, ja%val, ja%icol, ja%irow, &
                          atol=atol, btol=btol, conlim=conlim, itnlim=itnlim)
    call lsqr%solve(v, 0.0_wp, z, istop)

    jtz = 0.0_wp
    do j = 1, ja%nnz
        jtz(ja%icol(j)) = jtz(ja%icol(j)) + ja%val(j)*z(ja%irow(j))
    end do

    out = v - jtz

    end subroutine project_null
!*******************************************************************************

!*******************************************************************************
!>
!  projected conjugate gradients (Gould, Hribar & Nocedal 1998; Nocedal &
!  Wright Ch. 16): (approximately) solve the equality-constrained
!  subproblem `min 0.5*u^T*H*u + g^T*u` s.t. `ja*u = rhs_active` (implicit
!  in `ja`/the current working set) starting from a point already on that
!  face, returning the accumulated step `d_total`. If negative/zero
!  curvature is detected (`H` is not guaranteed positive-definite on the
!  null space), CG is truncated (Steihaug-Toint style) and the direction
!  at truncation is returned separately in `d_extra` with `truncated=.true.`
!  (the caller ratio-tests an *unbounded* move along `d_extra` in addition
!  to the already-accumulated `d_total`).

    subroutine projected_cg(hessian, ja, n_active, n, hu_g0, gproj0, max_pcg, opt_tol, curv_tol, &
                             lsqr_atol, lsqr_btol, lsqr_conlim, lsqr_itnlim, &
                             d_total, d_extra, truncated)

    type(sqpopt_hessian_type),  intent(inout) :: hessian
    type(sqpopt_sparse_matrix), intent(in)    :: ja
    integer,                    intent(in)    :: n_active !! number of active rows in `ja`
    integer,                    intent(in)    :: n !! number of columns in `ja`
    integer,                    intent(in)    :: max_pcg !!
    real(wp), dimension(n),     intent(in)    :: hu_g0   !! H*u0+g at the starting point
    real(wp), dimension(n),     intent(in)    :: gproj0  !! project_null(ja,hu_g0) at the starting point
    real(wp),                   intent(in)    :: opt_tol !! optimality tolerance
    real(wp),                   intent(in)    :: curv_tol !! curvature tolerance
    real(wp),                   intent(in)    :: lsqr_atol    !! `LSQR` tolerances
    real(wp),                   intent(in)    :: lsqr_btol    !! `LSQR` tolerances
    real(wp),                   intent(in)    :: lsqr_conlim  !! `LSQR` tolerances
    integer,                    intent(in)    :: lsqr_itnlim  !! `LSQR` max iterations
    real(wp), dimension(n),     intent(out)   :: d_total !! accumulated step
    real(wp), dimension(n),     intent(out)   :: d_extra !! truncation direction (only meaningful if truncated)
    logical,                    intent(out)   :: truncated !! whether CG was truncated due to negative/zero curvature

    real(wp), dimension(n) :: r, gproj, dvec, hd
    real(wp) :: rg_old, rg_new, kappa, alpha, beta
    integer :: j

    d_total = 0.0_wp
    d_extra = 0.0_wp
    truncated = .false.

    r     = hu_g0
    gproj = gproj0
    dvec  = -gproj
    rg_old = dot_product(r, gproj)

    do j = 1, max_pcg

        if (norm2(gproj) <= opt_tol) exit

        call hessian%hv_product(dvec, hd)
        kappa = dot_product(dvec, hd)

        if (kappa <= curv_tol) then
            d_extra   = dvec
            truncated = .true.
            return
        end if

        alpha  = rg_old/kappa
        d_total = d_total + alpha*dvec
        r       = r + alpha*hd
        call project_null(ja, n_active, n, r, gproj, lsqr_atol, lsqr_btol, lsqr_conlim, lsqr_itnlim)
        rg_new = dot_product(r, gproj)
        if (abs(rg_old) <= tiny(1.0_wp)) exit
        beta   = rg_new/rg_old
        dvec   = -gproj + beta*dvec
        rg_old = rg_new

    end do

    end subroutine projected_cg
!*******************************************************************************

!*******************************************************************************
!>
!  gather the currently-active rows (`status/=0`) of the combined
!  `mtot x n` sparse row set into a `n_active x n` sparse sub-matrix `ja`
!  and their target right-hand-side values `rhs_active`, along with
!  `orig_idx`, mapping each active row back to its index in `1..mtot`.

    subroutine build_active_set(arows, row_lb, row_ub, status, mtot, n, ja, rhs_active, orig_idx, n_active)

    type(sqpopt_sparse_matrix), intent(in)  :: arows !! the combined `mtot x n` sparse row set
    integer,                    intent(in)  :: mtot !! total number of rows in the combined sparse row set
    integer,                    intent(in)  :: n    !! number of columns in the combined sparse row set
    real(wp), dimension(mtot),  intent(in)  :: row_lb !! lower bounds for each row
    real(wp), dimension(mtot),  intent(in)  :: row_ub !! upper bounds for each row
    integer,  dimension(mtot),  intent(in)  :: status !! status of each row (0 = inactive, -1 = active at lower bound, 1 = active at upper bound)
    type(sqpopt_sparse_matrix), intent(out) :: ja !! the `n_active x n` sparse sub-matrix of active rows
    real(wp), dimension(:), allocatable, intent(out) :: rhs_active !! right-hand-side values for the active rows
    integer,  dimension(:), allocatable, intent(out) :: orig_idx !! mapping of each active row back to its index in `1..mtot`
    integer,                    intent(out) :: n_active !! number of active rows

    integer, dimension(:), allocatable :: row_map !! mapping of each row in the combined sparse row set to its index in the active set (0 if inactive)
    integer :: k, idx, nnz_a, j

    n_active = count(status /= 0)
    allocate(rhs_active(n_active), orig_idx(n_active))
    allocate(row_map(mtot)); row_map = 0

    idx = 0
    do k = 1, mtot
        if (status(k) /= 0) then
            idx = idx + 1
            row_map(k) = idx
            orig_idx(idx) = k
            rhs_active(idx) = merge(row_lb(k), row_ub(k), status(k) == -1)
        end if
    end do

    nnz_a = count(row_map(arows%irow(1:arows%nnz)) > 0)
    ja%nrows = n_active
    ja%ncols = n
    ja%nnz   = nnz_a
    allocate(ja%irow(nnz_a), ja%icol(nnz_a), ja%val(nnz_a))

    idx = 0
    do j = 1, arows%nnz
        if (row_map(arows%irow(j)) > 0) then
            idx = idx + 1
            ja%irow(idx) = row_map(arows%irow(j))
            ja%icol(idx) = arows%icol(j)
            ja%val(idx)  = arows%val(j)
        end if
    end do

    end subroutine build_active_set
!*******************************************************************************

!*******************************************************************************
!>
!  adjust `u` (in place) by the minimum-norm correction needed so that
!  every currently-active row exactly satisfies its target bound value
!  (bootstraps a feasible-for-the-working-set starting point). Uses `LSQR`
!  in its normal orientation (minimum-norm solution of the underdetermined
!  system `ja*correction = resid`), the same way the composite step
!  already uses it for its own normal step.

    subroutine project_onto_active(arows, row_lb, row_ub, status, mtot, n, u, atol, btol, conlim, itnlim)

    type(sqpopt_sparse_matrix), intent(in)    :: arows ! ! sparse matrix of active rows
    integer,                    intent(in)    :: mtot   !! total number of rows in the original constraint matrix
    integer,                    intent(in)    :: n      !! total number of columns in the original constraint matrix
    real(wp), dimension(mtot),  intent(in)    :: row_lb !! lower bounds for the rows of the original constraint matrix
    real(wp), dimension(mtot),  intent(in)    :: row_ub !! upper bounds for the rows of the original constraint matrix
    integer,  dimension(mtot),  intent(in)    :: status !! status of the rows of the original constraint matrix
    real(wp), dimension(n),     intent(inout) :: u      !! current solution vector
    real(wp),                   intent(in)    :: atol   !! `LSQR` tolerances (see [[sqpopt_reduced_hessian_qp_type]])
    real(wp),                   intent(in)    :: btol   !! `LSQR` tolerances (see [[sqpopt_reduced_hessian_qp_type]])
    real(wp),                   intent(in)    :: conlim !! `LSQR` tolerances (see [[sqpopt_reduced_hessian_qp_type]])
    integer,                    intent(in)    :: itnlim !! `LSQR` max iterations

    type(sqpopt_sparse_matrix) :: ja
    type(lsqr_solver_ez) :: lsqr
    real(wp), dimension(:), allocatable :: rhs_active, resid, correction
    integer,  dimension(:), allocatable :: orig_idx
    integer :: n_active, j, istop

    call build_active_set(arows, row_lb, row_ub, status, mtot, n, ja, rhs_active, orig_idx, n_active)
    if (n_active == 0) return

    allocate(resid(n_active), correction(n))
    resid = rhs_active
    do j = 1, ja%nnz
        resid(ja%irow(j)) = resid(ja%irow(j)) - ja%val(j)*u(ja%icol(j))
    end do

    call lsqr%initialize(n_active, n, ja%val, ja%irow, ja%icol, &
                          atol=atol, btol=btol, conlim=conlim, itnlim=itnlim)
    call lsqr%solve(resid, 0.0_wp, correction, istop)
    u = u + correction

    end subroutine project_onto_active
!*******************************************************************************

!*******************************************************************************
!>
!  update `row_lb`/`row_ub`-relative `status` for any general-constraint
!  row (`1..m`) currently violated at `u` (used during phase 1
!  bootstrapping, mirroring the equivalent bound-row check already done
!  inline in [[solve_reduced_hessian_qp]]).

    subroutine sparse_row_value(arows, u, m, is_equality, status, row_lb, row_ub, tol)

    type(sqpopt_sparse_matrix), intent(in)    :: arows !! sparse matrix of the general constraint rows
    real(wp), dimension(:),     intent(in)    :: u !! current solution vector
    integer,                    intent(in)    :: m !! number of general constraint rows
    logical,  dimension(:),     intent(in)    :: is_equality !! indicates whether each general constraint row is an equality constraint
    integer,  dimension(:),     intent(inout) :: status !! current status of each general constraint row
    real(wp), dimension(:),     intent(in)    :: row_lb !! lower bounds for the general constraint rows
    real(wp), dimension(:),     intent(in)    :: row_ub !! upper bounds for the general constraint rows
    real(wp),                   intent(in)    :: tol !! tolerance for checking constraint violations

    integer :: k !! loop index for the general constraint rows
    real(wp) :: val !! value of the current general constraint row at `u`

    do k = 1, m
        if (is_equality(k)) cycle
        val = sparse_dot_row(arows, k, u)
        if (val < row_lb(k) - tol) then
            status(k) = -1
        else if (val > row_ub(k) + tol) then
            status(k) = 1
        end if
    end do

    end subroutine sparse_row_value
!*******************************************************************************

    end module sqpopt_qp_reduced_hessian_module
!*******************************************************************************
