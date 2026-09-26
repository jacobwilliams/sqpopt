!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  **Sparse** active-set QP solver for the linearized SQP subproblem (the
!  default for larger problems, see `sqpopt_qp_auto`). Like
!  [[sqpopt_qp_dense_module]], this enforces the linearized constraints
!  and bounds exactly, but never forms a dense matrix: the Hessian is only
!  used through products, and active variable bounds simply fix their
!  variables. The null space of the working set is handled in one of two
!  ways (`null_space`):
!
!  * `sqpopt_null_space_lu` (the default, SQOPT-style): every general row
!    gets a slack variable (\( J p - s = 0 \)), so all the constraints
!    are simple bounds on the unknowns, and the working set is the set of
!    unknowns fixed at a bound. The free unknowns are split into a square
!    nonsingular **basis** `B` (sparse LU factors, from `LUSOL`) and the
!    **superbasic** rest `S`, so the null space is spanned by
!    \( Z = [-B^{-1} S; I] \). Every product with `Z` or \( Z^T \) is
!    one solve with `B`'s factors, the reduced-space Newton system
!    \( Z^T H Z \, d_S = -Z^T (Hv+g) \) is solved by conjugate gradients,
!    and the multipliers come from one solve with \( B^T \). All of these
!    are direct (not iterative) solves. Fixing or freeing a superbasic
!    unknown leaves `B` unchanged; fixing a basic one swaps in a superbasic
!    column with a Bartels-Golub update of the factors (`lu8rpc`), so `B` is
!    only refactorized occasionally.
!  * `sqpopt_null_space_lsqr`: orthogonal projections onto the null space,
!    each a least-squares solve with `LSQR`, with **projected conjugate
!    gradients** (Gould, Hribar & Nocedal 1998; Nocedal & Wright,
!    *Numerical Optimization*, Ch. 16).
!
!  The formulation is the same as the dense solver's: general constraints
!  and bounds are uniform two-sided "rows" (stored in compressed-row form,
!  so each row's product costs only its own nonzeros), and every general
!  constraint violated at `p=0` gets an elastic slack with an \( \ell_1 \)
!  penalty, so the starting step plus those slacks is feasible, and
!  inconsistent linearized constraints are detected (the penalty weight is
!  raised up to `elastic_weight_max`, then `istat=sqpopt_infeasible`); see
!  [[sqpopt_qp_dense_module]] for the details (here, the slacks also get a
!  small proximal curvature, so CG steps along them stay bounded). Projected
!  CG stops on
!  (relative) convergence, and follows any direction of nonpositive
!  curvature (e.g. along an elastic slack, or an indefinite SR1 Hessian) to
!  the nearest blocking row. The initial working set is a linearly
!  independent subset of the rows at a bound, picked by one rank-revealing
!  sparse LU factorization (`LUSOL`); every row added after that is
!  independent by construction (a blocking row is never in the span of
!  the working set).

    module sqpopt_qp_reduced_hessian_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_qp_solve_failed, &
                                     sqpopt_infeasible, sqpopt_infinity
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_linalg_module,  only: independent_columns, sqpopt_lu_type
    use lsqr_module,           only: lsqr_solver_ez

    implicit none

    private

    integer, parameter, public :: sqpopt_null_space_lu   = 1 !! (default) basis partition of the working set, with sparse
                                                              !! LU factors of the basis (SQOPT-style, see the module docs)
    integer, parameter, public :: sqpopt_null_space_lsqr = 2 !! orthogonal projections onto the null space, by `LSQR`

    type, public :: sqpopt_reduced_hessian_qp_type
        !! options for the sparse (projected-CG) active-set QP solver.

        integer  :: null_space   = sqpopt_null_space_lu !! how the null space of the working set is handled
                                                        !! (see the `sqpopt_null_space_*` constants)
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

        ! `LSQR` settings, used for the minimum-norm starting step, and (with
        ! `null_space=sqpopt_null_space_lsqr`) for every null-space projection
        ! and multiplier solve (see [[lsqr_module]] for the precise meaning of
        ! each). `0` for `lsqr_atol`/`lsqr_btol`/`lsqr_conlim` means "let LSQR
        ! use its own machine-precision-based default"; loosening these trades
        ! QP-solve accuracy for speed:
        real(wp) :: lsqr_atol   = 0.0_wp !! `LSQR` relative error tolerance in `A` (0 => `LSQR` default)
        real(wp) :: lsqr_btol   = 0.0_wp !! `LSQR` relative error tolerance in `b` (0 => `LSQR` default)
        real(wp) :: lsqr_conlim = 0.0_wp !! `LSQR` upper limit on `cond(Abar)` (0 => `LSQR` default)
        integer  :: lsqr_itnlim = 0      !! `LSQR` maximum iterations per solve (`<=0` => `2*(rows+columns)+10`)

        logical  :: warm_start = .true.  !! start from the previous solve's final working set, if the problem
                                         !! size is unchanged (see [[sqpopt_qp_dense_module]]: the crash and
                                         !! warm starts work the same way here, with `LSQR` computing the
                                         !! minimum-norm starting step)

        integer :: n_iter = 0 !! number of active-set iterations taken by the last solve (output)

        ! internal state (the working set at the end of the previous solve, for warm starts):
        integer, dimension(:), allocatable :: warm_status !! side (-1/0/+1) of each general row and variable bound

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
    real(wp), dimension(:), allocatable :: row_lb, row_ub, u, hu_g, gproj, d_total, d_extra, coeff, s_sign, s0, jp0
    real(wp), dimension(size(g)) :: p0
    integer,  dimension(:), allocatable :: status, orig_idx, coeff_idx, slack_row
    logical,  dimension(:), allocatable :: is_equality, fixed
    real(wp) :: rho, rho_max, gscale, alpha, alpha_cap, scale
    integer  :: n_active, blocking, blocking_side
    logical  :: at_face_optimum, truncated
    logical  :: lu_ok
    ! the basis method's state (see `basis_active_set`):
    integer, parameter :: max_updates = 100  !! refactorize `B` after this many column replacements
    integer  :: nn     !! number of unknowns `v = (p, e, s)`
    integer  :: nupd   !! column replacements since `B` was last factorized
    integer,  dimension(:), allocatable :: cptr, crow  !! the constraint matrix `[J E -I]`, by columns
    real(wp), dimension(:), allocatable :: cval
    real(wp), dimension(:), allocatable :: v, vlb, vub !! the unknowns and their bounds
    integer,  dimension(:), allocatable :: state       !! 0 = free, -1/+1 = fixed at the lower/upper bound
    integer,  dimension(:), allocatable :: bpos        !! position of each unknown in the basis (0 if not basic)
    integer,  dimension(:), allocatable :: bvar        !! the basic unknowns
    integer,  dimension(:), allocatable :: sup         !! the superbasic unknowns
    logical,  dimension(:), allocatable :: is_eqv      !! unknowns with equal bounds
    type(sqpopt_lu_type) :: blu                        !! LU factors of `B`

    n = size(g)
    m = size(c)

    ! ---- starting step: crash or warm start (see [[sqpopt_qp_dense_module]]) ----
    call starting_step(p0)

    ! ---- elastic slacks: one for each general row violated at p0 ----
    allocate(row_lb(m), row_ub(m), s_sign(m), slack_row(m), s0(m), jp0(m))
    row_lb = c_lb - c
    row_ub = c_ub - c
    jp0 = 0.0_wp
    do k = 1, jac%nnz
        jp0(jac%irow(k)) = jp0(jac%irow(k)) + jac%val(k)*p0(jac%icol(k))
    end do
    nv = 0
    do i = 1, m
        s_sign(i) = 0.0_wp
        if (jp0(i) < row_lb(i) - me%feas_tol*max(1.0_wp, abs(row_lb(i)))) then
            s_sign(i) = 1.0_wp
        else if (jp0(i) > row_ub(i) + me%feas_tol*max(1.0_wp, abs(row_ub(i)))) then
            s_sign(i) = -1.0_wp
        end if
        if (s_sign(i) /= 0.0_wp) then
            nv = nv + 1
            slack_row(nv) = i
            s0(nv) = merge(row_lb(i) - jp0(i), jp0(i) - row_ub(i), s_sign(i) > 0.0_wp)
        end if
    end do
    s0 = s0(1:nv)
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

    ! ---- feasible starting point: p0, with the slacks just large enough ----
    allocate(u(nt))
    u(1:n)    = p0
    u(n+1:nt) = s0

    ! ---- the basis method (the default), falling back to the LSQR method if it fails ----
    if (me%null_space == sqpopt_null_space_lu) then
        call basis_active_set(lu_ok)
        if (lu_ok) return
        u(1:n)    = p0
        u(n+1:nt) = s0
        rho       = me%elastic_weight*gscale
    end if

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

    me%n_iter = min(it, maxit)
    p = u(1:n)
    lambda = 0.0_wp
    do k = 1, size(coeff_idx)
        if (coeff_idx(k) <= m) lambda(coeff_idx(k)) = coeff(k)
    end do

    ! remember the final working set (general rows and variable bounds) for a warm start:
    me%warm_status = status(1:m+n)

    contains

        subroutine starting_step(p0)
        !! the minimum-norm step satisfying the initial working-set guess (the
        !! previous solve's final working set, or the equality constraints and
        !! fixed variables): guessed bounds fix their variables, and `LSQR`
        !! gives the minimum-norm solution for the free ones of the guessed
        !! general rows. Any violated variable bounds are added to the guess
        !! (up to 4 rounds), then the step is clipped to the bounds.
        real(wp), dimension(n), intent(out) :: p0
        real(wp), dimension(n) :: blb, bub
        integer,  dimension(m+n) :: guess
        logical,  dimension(n) :: fix
        integer,  dimension(m) :: gmap
        real(wp), dimension(:), allocatable :: rhs, pf
        integer,  dimension(:), allocatable :: ir, ic
        real(wp), dimension(:), allocatable :: vv
        type(lsqr_solver_ez) :: lsqr
        integer :: kk, round, ng, nnz_g, istop
        logical :: added

        blb = x_lb - x
        bub = x_ub - x
        guess = 0
        if (me%warm_start .and. allocated(me%warm_status)) then
            if (size(me%warm_status) == m+n) guess = me%warm_status
        end if
        if (all(guess == 0)) then
            do kk = 1, m
                if (c_ub(kk) - c_lb(kk) <= 0.0_wp) guess(kk) = -1   ! equality constraint
            end do
            do kk = 1, n
                if (bub(kk) - blb(kk) <= 0.0_wp) guess(m+kk) = -1   ! fixed variable
            end do
        end if
        ! (a guessed side whose bound is infinite can't be targeted)
        do kk = 1, m
            if (guess(kk) == -1 .and. c_lb(kk)-c(kk) <= -sqpopt_infinity) guess(kk) = 0
            if (guess(kk) ==  1 .and. c_ub(kk)-c(kk) >=  sqpopt_infinity) guess(kk) = 0
        end do
        do kk = 1, n
            if (guess(m+kk) == -1 .and. blb(kk) <= -sqpopt_infinity) guess(m+kk) = 0
            if (guess(m+kk) ==  1 .and. bub(kk) >=  sqpopt_infinity) guess(m+kk) = 0
        end do

        do round = 1, 4
            ! the fixed variables, at their guessed bounds:
            p0 = 0.0_wp
            fix = guess(m+1:m+n) /= 0
            where (guess(m+1:m+n) == -1) p0 = blb
            where (guess(m+1:m+n) ==  1) p0 = bub
            ! minimum-norm solution for the free variables of the guessed general rows:
            ng = 0
            gmap = 0
            do kk = 1, m
                if (guess(kk) /= 0) then
                    ng = ng + 1
                    gmap(kk) = ng
                end if
            end do
            if (ng > 0) then
                allocate(rhs(ng))
                do kk = 1, m
                    if (gmap(kk) > 0) rhs(gmap(kk)) = merge(c_lb(kk)-c(kk), c_ub(kk)-c(kk), guess(kk) == -1)
                end do
                nnz_g = 0
                do kk = 1, jac%nnz
                    if (gmap(jac%irow(kk)) == 0) cycle
                    if (fix(jac%icol(kk))) then
                        rhs(gmap(jac%irow(kk))) = rhs(gmap(jac%irow(kk))) - jac%val(kk)*p0(jac%icol(kk))
                    else
                        nnz_g = nnz_g + 1
                    end if
                end do
                if (nnz_g > 0) then
                    allocate(ir(nnz_g), ic(nnz_g), vv(nnz_g), pf(n))
                    nnz_g = 0
                    do kk = 1, jac%nnz
                        if (gmap(jac%irow(kk)) == 0 .or. fix(jac%icol(kk))) cycle
                        nnz_g = nnz_g + 1
                        ir(nnz_g) = gmap(jac%irow(kk)); ic(nnz_g) = jac%icol(kk); vv(nnz_g) = jac%val(kk)
                    end do
                    call lsqr%initialize(ng, n, vv, ir, ic, atol=me%lsqr_atol, btol=me%lsqr_btol, &
                                          conlim=me%lsqr_conlim, itnlim=2*(ng+n)+10)
                    call lsqr%solve(rhs, 0.0_wp, pf, istop)
                    where (.not. fix) p0 = pf
                    deallocate(ir, ic, vv, pf)
                end if
                deallocate(rhs)
            end if
            ! add any violated variable bounds to the guess, and try again:
            added = .false.
            do kk = 1, n
                if (guess(m+kk) /= 0) cycle
                if (p0(kk) < blb(kk)) then
                    guess(m+kk) = -1; added = .true.
                else if (p0(kk) > bub(kk)) then
                    guess(m+kk) = 1; added = .true.
                end if
            end do
            if (.not. added) exit
        end do
        p0 = min(max(p0, blb), bub)

        end subroutine starting_step

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

        subroutine basis_active_set(ok)
        !! the active-set iterations of the basis method (`null_space =
        !! sqpopt_null_space_lu`), SQOPT-style. Each general row gets a slack
        !! variable `s` (`J p + E e - s = 0`, with the elastic slacks `e`), so
        !! the unknowns are `v = (p, e, s)` with simple bounds only, and the
        !! working set is just the set of *fixed* unknowns (each at one of its
        !! bounds). The free unknowns are split into `m` **basic** ones, whose
        !! columns form a nonsingular basis `B` (with sparse LU factors), and
        !! the **superbasic** rest; a step `d_S` in the superbasics moves the
        !! basics by `d_B = -B^{-1} S d_S`, which stays on the constraints.
        !!
        !! Fixing a superbasic unknown or freeing a fixed one doesn't change
        !! `B`; fixing a basic one replaces its column of `B` by that of a
        !! superbasic (the one with the largest pivot), with a Bartels-Golub
        !! update of the factors (`lu8rpc`) rather than a refactorization.
        !!
        !! Sets `p`, `lambda`, `istat`, `me%n_iter`, and `me%warm_status`, and
        !! `ok=.true.`; or `ok=.false.` if the factorization failed.
        logical, intent(out) :: ok

        integer  :: j, k, l, r, iter, maxit_b, blk, side, st
        real(wp), dimension(:), allocatable :: gv, y, rs, dtot, dext
        logical,  dimension(:), allocatable :: chosen
        real(wp) :: alpha_b, sc, worst, tol_mult, dj
        logical  :: trunc, cg_done

        ok = .false.
        nn = nt + m

        ! ---- the constraint matrix `[J E -I]` (`m x nn`), by columns ----
        allocate(cptr(nn+1)); cptr = 0
        do r = 1, m
            do l = rows%ptr(r), rows%ptr(r+1)-1
                cptr(rows%col(l)+1) = cptr(rows%col(l)+1) + 1
            end do
            cptr(nt+r+1) = 1
        end do
        cptr(1) = 1
        do j = 1, nn
            cptr(j+1) = cptr(j+1) + cptr(j)
        end do
        allocate(crow(cptr(nn+1)-1), cval(cptr(nn+1)-1))
        block
            integer, dimension(nn) :: pos
            pos = cptr(1:nn)
            do r = 1, m
                do l = rows%ptr(r), rows%ptr(r+1)-1
                    j = rows%col(l)
                    crow(pos(j)) = r; cval(pos(j)) = rows%val(l); pos(j) = pos(j) + 1
                end do
                crow(pos(nt+r)) = r; cval(pos(nt+r)) = -1.0_wp
            end do
        end block

        ! ---- the unknowns and their bounds ----
        allocate(v(nn), vlb(nn), vub(nn), is_eqv(nn))
        v(1:nt) = u
        do r = 1, m
            v(nt+r) = row_dot(rows, r, u)
        end do
        vlb(1:nt) = row_lb(m+1:m+nt);  vub(1:nt) = row_ub(m+1:m+nt)
        vlb(nt+1:nn) = row_lb(1:m);    vub(nt+1:nn) = row_ub(1:m)
        do j = 1, nn
            is_eqv(j) = vub(j)-vlb(j) <= me%active_tol*max(1.0_wp, abs(vlb(j)))
        end do

        ! ---- initial basis and working set ----
        ! The unknowns at a bound are the candidates for the working set. The
        ! basis is picked from all the columns by one rank-revealing LU, on the
        ! row-normalized matrix, with each column scaled by how much it should
        ! be *free*: the free row slacks most (as exact unit columns, which
        ! LUSOL takes first: an inactive row's slack is its natural basic
        ! variable), then the other free unknowns, then inequality rows'
        ! slacks, then variable bounds, then equality constraints. The
        ! candidates not picked are fixed (the working set); those picked stay
        ! free, which drops them from the working set only if needed for a
        ! nonsingular basis (i.e. if the candidates are linearly dependent).
        ! (Normalizing the columns instead can make a poor basis: e.g. a
        ! variable with a single small entry would become a unit column.)
        allocate(state(nn), bpos(nn), bvar(m), chosen(nn)); state = 0; bpos = 0
        do j = 1, nn
            state(j) = at_bound_v(j)
            if (is_eqv(j)) state(j) = -1
        end do
        if (m > 0) then
            block
                real(wp), dimension(size(cval)) :: wv
                integer,  dimension(size(cval)) :: wc
                real(wp), dimension(m) :: rn
                real(wp) :: wt
                rn = 0.0_wp
                do j = 1, nt
                    do l = cptr(j), cptr(j+1)-1
                        rn(crow(l)) = rn(crow(l)) + cval(l)**2
                    end do
                end do
                rn = sqrt(rn)
                where (rn <= 0.0_wp) rn = 1.0_wp
                do j = 1, nn
                    if (state(j) == 0) then
                        wt = 1.0_wp
                    else if (is_eqv(j)) then
                        wt = 1.0e-6_wp
                    else if (j <= nt) then
                        wt = 1.0e-4_wp
                    else
                        wt = 1.0e-2_wp
                    end if
                    do l = cptr(j), cptr(j+1)-1
                        wv(l) = wt*cval(l)/rn(crow(l))
                        wc(l) = j
                    end do
                    if (j > nt .and. state(j) == 0) wv(cptr(j)) = -1.0_wp
                end do
                call independent_columns(m, nn, crow, wc, wv, 1.0e-8_wp, 1.0e-6_wp*epsilon(1.0_wp)**0.67_wp, &
                                         chosen, st)
            end block
            if (st /= 0 .or. count(chosen) /= m) return
            k = 0
            do j = 1, nn
                if (.not. chosen(j)) cycle
                k = k + 1
                bvar(k) = j
                bpos(j) = k
                state(j) = 0
            end do
            if (.not. factorize_basis()) return
        end if
        call update_basics()

        ! ---- active-set iterations ----
        istat   = sqpopt_qp_solve_failed
        maxit_b = max(me%max_iter, 10*(mtot+1))
        allocate(gv(nn), y(m), dtot(nn), dext(nn))
        y = 0.0_wp
        nupd = 0
        cg_done = .false.
        do iter = 1, maxit_b

            sup = pack([(j, j=1,nn)], state == 0 .and. bpos == 0)
            gv(1:nt) = gradient(v(1:nt))
            gv(nt+1:nn) = 0.0_wp
            if (allocated(rs)) deallocate(rs)
            allocate(rs(size(sup)))
            call zt_times(gv, rs, y)
            sc = 1.0_wp + maxval(abs(merge(gv, 0.0_wp, state == 0)))

            ! (after an unblocked CG step, the face is solved: the recomputed
            ! reduced gradient can't always meet the tolerance again, because of
            ! cancellation with the large elastic multipliers)
            if (norm2(rs) > me%opt_tol*sc .and. .not. cg_done) then

                call reduced_cg(rs, me%opt_tol*sc, dtot, dext, trunc)

                ! first, the step accumulated by CG (which may itself be blocked):
                call ratio_test_v(dtot, 1.0_wp, alpha_b, blk, side)
                v = v + alpha_b*dtot
                if (blk /= 0 .and. alpha_b < 1.0_wp - 1.0e-12_wp) then
                    if (.not. fix(blk, side)) return
                    cycle
                end if
                if (.not. trunc) then
                    cg_done = .true.     ! (then at the face optimum)
                    cycle
                end if
                ! then, if CG found nonpositive curvature, move along that
                ! (downhill) direction, not bounded by 1, to the nearest blocking bound:
                call ratio_test_v(dext, huge(1.0_wp), alpha_b, blk, side)
                if (blk == 0) exit  ! unbounded QP
                v = v + alpha_b*dext
                if (.not. fix(blk, side)) return
                cycle

            end if

            cg_done = .false.

            ! optimal on this face. Free the fixed unknown with the most wrongly
            ! signed reduced cost `g_j - a_j^T y` (its multiplier), if any
            ! (relative to the largest one, not counting the elastic slacks'
            ! bounds, whose multipliers are the large penalty weight):
            tol_mult = 0.0_wp
            do j = 1, nn
                if (state(j) == 0 .or. (j > n .and. j <= nt)) cycle
                tol_mult = max(tol_mult, abs(gv(j) - col_dot(j, y)))
            end do
            tol_mult = me%active_tol*max(1.0_wp, tol_mult)
            worst = tol_mult
            k = 0
            do j = 1, nn
                if (state(j) == 0 .or. is_eqv(j)) cycle
                dj = gv(j) - col_dot(j, y)
                if (state(j) == -1 .and. -dj > worst) then
                    worst = -dj; k = j
                else if (state(j) == 1 .and. dj > worst) then
                    worst = dj; k = j
                end if
            end do
            if (k /= 0) then
                state(k) = 0
                cycle
            end if

            ! optimal for the current elastic weight. Any slack still positive?
            if (any(v(n+1:nt) > me%feas_tol*max(1.0_wp, s0))) then
                if (rho < rho_max) then
                    rho = min(100.0_wp*rho, rho_max)
                    cycle
                end if
                istat = sqpopt_infeasible
            else
                istat = sqpopt_success
            end if
            exit

        end do

        me%n_iter = min(iter, maxit_b)
        u = v(1:nt)
        p = v(1:n)
        lambda = 0.0_wp
        do r = 1, m
            if (state(nt+r) /= 0) lambda(r) = y(r)
        end do
        me%warm_status = [state(nt+1:nn), state(1:n)]
        ok = .true.

        end subroutine basis_active_set

        integer function at_bound_v(jj)
        !! -1 or +1 if unknown `jj` is at its lower or upper bound, else 0
        integer, intent(in) :: jj
        at_bound_v = 0
        if (vlb(jj) > -sqpopt_infinity) then
            if (abs(v(jj)-vlb(jj)) <= me%active_tol*max(1.0_wp, abs(vlb(jj)))) at_bound_v = -1
        end if
        if (at_bound_v == 0 .and. vub(jj) < sqpopt_infinity) then
            if (abs(v(jj)-vub(jj)) <= me%active_tol*max(1.0_wp, abs(vub(jj)))) at_bound_v = 1
        end if
        end function at_bound_v

        pure real(wp) function col_dot(jj, w)
        !! `a_jj^T w`
        integer,                intent(in) :: jj
        real(wp), dimension(:), intent(in) :: w
        col_dot = dot_product(cval(cptr(jj):cptr(jj+1)-1), w(crow(cptr(jj):cptr(jj+1)-1)))
        end function col_dot

        logical function factorize_basis()
        !! factorize `B` (the columns `bvar`) from scratch
        integer,  dimension(:), allocatable :: bi, bj
        real(wp), dimension(:), allocatable :: bv
        integer :: kk, stat
        bi = [(crow(cptr(bvar(kk)):cptr(bvar(kk)+1)-1), kk=1,m)]
        bj = [(spread(kk, 1, cptr(bvar(kk)+1)-cptr(bvar(kk))), kk=1,m)]
        bv = [(cval(cptr(bvar(kk)):cptr(bvar(kk)+1)-1), kk=1,m)]
        call blu%factorize(m, bi, bj, bv, 1.0e-12_wp, stat)
        factorize_basis = stat == 0
        nupd = 0
        end function factorize_basis

        subroutine update_basics()
        !! the basic unknowns from the constraints, given the others:
        !! `B v_B = -(sum of a_j v_j over the nonbasic j)` (this also removes
        !! any drift from the constraints)
        real(wp), dimension(m) :: rhs, vb
        integer :: jj
        if (m == 0) return
        rhs = 0.0_wp
        do jj = 1, nn
            if (bpos(jj) /= 0 .or. v(jj) == 0.0_wp) cycle
            rhs(crow(cptr(jj):cptr(jj+1)-1)) = rhs(crow(cptr(jj):cptr(jj+1)-1)) - cval(cptr(jj):cptr(jj+1)-1)*v(jj)
        end do
        call blu%solve(rhs, vb, transpose=.false.)
        v(bvar) = vb
        end subroutine update_basics

        subroutine z_times(vs, d)
        !! `d = Z vs`: `vs` on the superbasics, then the basics from `B d_B = -S vs`
        real(wp), dimension(:), intent(in)  :: vs
        real(wp), dimension(:), intent(out) :: d
        real(wp), dimension(m) :: rhs, db
        integer :: kk, jj
        d = 0.0_wp
        if (size(sup) > 0) d(sup) = vs
        if (m == 0) return
        rhs = 0.0_wp
        do kk = 1, size(sup)
            jj = sup(kk)
            rhs(crow(cptr(jj):cptr(jj+1)-1)) = rhs(crow(cptr(jj):cptr(jj+1)-1)) - cval(cptr(jj):cptr(jj+1)-1)*vs(kk)
        end do
        call blu%solve(rhs, db, transpose=.false.)
        d(bvar) = db
        end subroutine z_times

        subroutine zt_times(w, rr, yy)
        !! `rr = Z^T w = w_S - S^T yy`, with `B^T yy = w_B` (at a face
        !! optimum, `yy` are the general rows' multipliers)
        real(wp), dimension(:), intent(in)  :: w
        real(wp), dimension(:), intent(out) :: rr
        real(wp), dimension(:), intent(out) :: yy
        integer :: kk
        if (m > 0) then
            call blu%solve(w(bvar), yy, transpose=.true.)
        end if
        do kk = 1, size(sup)
            rr(kk) = w(sup(kk))
            if (m > 0) rr(kk) = rr(kk) - col_dot(sup(kk), yy)
        end do
        end subroutine zt_times

        subroutine hv_product(d, hd)
        !! the QP Hessian (in `v`) times `d`: zero on the row slacks
        real(wp), dimension(:), intent(in)  :: d
        real(wp), dimension(:), intent(out) :: hd
        call hext_product(d(1:nt), hd(1:nt))
        hd(nt+1:nn) = 0.0_wp
        end subroutine hv_product

        subroutine reduced_cg(rg0, abs_tol, d_total, d_extra, truncated)
        !! conjugate gradients on the reduced Hessian `Z^T H Z` (in the
        !! superbasics): returns the accumulated step `d_total = Z d_S`; if a
        !! direction of nonpositive curvature is found, it is returned
        !! (downhill) in `d_extra` with `truncated=.true.`. Every step stays
        !! on the constraints by construction.
        real(wp), dimension(:), intent(in)  :: rg0     !! reduced gradient
        real(wp),               intent(in)  :: abs_tol !! absolute stopping tolerance on it
        real(wp), dimension(:), intent(out) :: d_total, d_extra
        logical,                intent(out) :: truncated
        real(wp), dimension(size(rg0)) :: rr, ds, hd
        real(wp), dimension(m)  :: yy
        real(wp), dimension(nn) :: zd, hzd
        real(wp) :: rr_old, rr_new, kappa, alpha, tol
        integer :: jj, max_it
        d_total = 0.0_wp
        d_extra = 0.0_wp
        truncated = .false.
        rr = rg0
        ds = -rr
        rr_old = dot_product(rr, rr)
        tol = max(abs_tol, me%pcg_rtol*norm2(rg0))
        max_it = min(max_pcg, max(1, size(rg0)))
        do jj = 1, max_it
            if (sqrt(rr_old) <= tol) exit
            call z_times(ds, zd)
            call hv_product(zd, hzd)
            call zt_times(hzd, hd, yy)
            kappa = dot_product(ds, hd)
            if (kappa <= 1.0e-10_wp*norm2(ds)*norm2(hd)) then
                d_extra = zd
                truncated = norm2(zd) > 0.0_wp
                return
            end if
            alpha   = rr_old/kappa
            d_total = d_total + alpha*zd
            rr      = rr + alpha*hd
            rr_new  = dot_product(rr, rr)
            ds      = -rr + (rr_new/rr_old)*ds
            rr_old  = rr_new
        end do
        end subroutine reduced_cg

        subroutine ratio_test_v(d, alpha_cap, alpha, blocking, blocking_side)
        !! the largest `alpha <= alpha_cap` for which `v+alpha*d` satisfies the
        !! bounds of every free unknown, and the unknown (and bound) that blocks first
        real(wp), dimension(:), intent(in)  :: d
        real(wp),               intent(in)  :: alpha_cap
        real(wp),               intent(out) :: alpha
        integer,                intent(out) :: blocking, blocking_side
        real(wp) :: alpha_k, dtol
        integer :: jj
        alpha = alpha_cap
        blocking = 0
        blocking_side = 0
        dtol = 1.0e-12_wp*maxval(abs(d))
        do jj = 1, nn
            if (state(jj) /= 0 .or. abs(d(jj)) <= dtol) cycle
            if (d(jj) > 0.0_wp) then
                if (vub(jj) >= sqpopt_infinity) cycle
                alpha_k = max((vub(jj) - v(jj))/d(jj), 0.0_wp)
                if (alpha_k < alpha) then
                    alpha = alpha_k; blocking = jj; blocking_side = 1
                end if
            else
                if (vlb(jj) <= -sqpopt_infinity) cycle
                alpha_k = max((vlb(jj) - v(jj))/d(jj), 0.0_wp)
                if (alpha_k < alpha) then
                    alpha = alpha_k; blocking = jj; blocking_side = -1
                end if
            end if
        end do
        end subroutine ratio_test_v

        logical function fix(jj, sd)
        !! fix unknown `jj` at its bound `sd` (add it to the working set). If
        !! it is basic, the superbasic with the largest pivot in its row of
        !! `B^{-1} S` replaces it in the basis. `.false.` if the factorization failed.
        integer, intent(in) :: jj, sd
        real(wp), dimension(m) :: e_p, w
        real(wp) :: piv, best
        integer  :: pp, kk, q, stat
        fix = .true.
        state(jj) = sd
        v(jj) = merge(vlb(jj), vub(jj), sd == -1)
        pp = bpos(jj)
        if (pp /= 0) then
            e_p = 0.0_wp
            e_p(pp) = 1.0_wp
            call blu%solve(e_p, w, transpose=.true.)
            q = 0
            best = 0.0_wp
            do kk = 1, size(sup)
                if (sup(kk) == jj) cycle
                piv = abs(col_dot(sup(kk), w))
                if (piv > best) then
                    best = piv; q = sup(kk)
                end if
            end do
            if (q == 0 .or. best <= 1.0e-12_wp*max(1.0_wp, maxval(abs(w)))) then
                fix = .false.   ! (can't happen in exact arithmetic: the basic unknown couldn't have moved)
                return
            end if
            bpos(jj) = 0
            bpos(q)  = pp
            bvar(pp) = q
            nupd = nupd + 1
            stat = 1
            if (nupd <= max_updates) &
                call blu%replace_column(pp, crow(cptr(q):cptr(q+1)-1), cval(cptr(q):cptr(q+1)-1), stat)
            if (stat /= 0) then
                fix = factorize_basis()
                if (.not. fix) return
            end if
        end if
        call update_basics()
        end function fix


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
        !! add a linearly independent subset of the rows that are at a bound
        !! at `u` to the working set. The subset is picked by one rank-revealing
        !! LU factorization of the candidate rows (see [[independent_columns]]),
        !! each normalized and then scaled by its priority: equality rows (kept
        !! first), then variable bounds, then inequality rows. If the
        !! factorization fails, the rows are checked one at a time instead
        !! (with an `LSQR` projection each).
        real(wp), parameter :: w_bound = 1.0e-2_wp  !! priority weight of an (inequality) variable bound
        real(wp), parameter :: w_ineq  = 1.0e-4_wp  !! priority weight of an inequality general row
        integer,  dimension(mtot) :: cand, cside
        integer,  dimension(:), allocatable :: ir, ic
        real(wp), dimension(:), allocatable :: vv
        logical,  dimension(:), allocatable :: indep
        real(wp) :: wt, rn
        integer  :: kk, pass, side, nc, nz, j, l, lu_stat

        ! the candidates: rows at a bound (equality rows and bounds first, as
        ! for the fallback below):
        nc = 0
        do pass = 1, 2
            do kk = 1, mtot
                if (status(kk) /= 0) cycle
                if ((pass == 1) .neqv. (is_equality(kk) .or. kk > m)) cycle
                side = at_bound(kk)
                if (side == 0) cycle
                nc = nc + 1
                cand(nc)  = kk
                cside(nc) = side
            end do
        end do
        if (nc == 0) return

        ! the candidates as the columns of an `nt x nc` matrix:
        nz = 0
        do j = 1, nc
            nz = nz + rows%ptr(cand(j)+1) - rows%ptr(cand(j))
        end do
        allocate(ir(nz), ic(nz), vv(nz), indep(nc))
        nz = 0
        do j = 1, nc
            kk = cand(j)
            if (is_equality(kk)) then
                wt = 1.0_wp
            else if (kk > m) then
                wt = w_bound
            else
                wt = w_ineq
            end if
            rn = row_norm(rows, kk)
            if (rn > 0.0_wp) wt = wt/rn
            do l = rows%ptr(kk), rows%ptr(kk+1)-1
                nz = nz + 1
                ir(nz) = rows%col(l)
                ic(nz) = j
                vv(nz) = wt*rows%val(l)
            end do
        end do
        call independent_columns(nt, nc, ir, ic, vv, 1.0e-8_wp, epsilon(1.0_wp)**0.67_wp*w_ineq, indep, lu_stat)

        if (lu_stat == 0) then
            do j = 1, nc
                if (indep(j)) status(cand(j)) = cside(j)
            end do
        else
            call add_independent_rows_lsqr(cand(1:nc), cside(1:nc))
        end if

        end subroutine initial_working_set

        subroutine add_independent_rows_lsqr(cand, cside)
        !! the fallback for [[initial_working_set]]: add the candidate rows `cand`
        !! (at bound side `cside`) in order, skipping any whose component outside
        !! the span of the rows already added is negligible (i.e., that is
        !! linearly dependent). Each check costs an `LSQR` solve.
        integer, dimension(:), intent(in) :: cand, cside
        type(sqpopt_sparse_matrix) :: ja_cur
        integer, dimension(:), allocatable :: idx_cur
        logical, dimension(:), allocatable :: fixed_cur
        real(wp), dimension(nt) :: a, r
        integer :: jj, kk, na, j
        do jj = 1, size(cand)
            kk = cand(jj)
            a = 0.0_wp
            do j = rows%ptr(kk), rows%ptr(kk+1)-1
                a(rows%col(j)) = rows%val(j)
            end do
            call build_working_set(rows, status, m, ja_cur, fixed_cur, idx_cur, na)
            if (na >= nt) return
            call project_null(ja_cur, fixed_cur, a, r)
            if (norm2(r) > 1.0e-8_wp*norm2(a)) status(kk) = cside(jj)
        end do
        end subroutine add_independent_rows_lsqr

        integer function at_bound(kk)
        !! -1 or +1 if row `kk` is at its lower or upper bound at `u`, else 0
        integer, intent(in) :: kk
        real(wp) :: val
        val = row_dot(rows, kk, u)
        if (abs(val-row_lb(kk)) <= me%active_tol*max(1.0_wp, abs(row_lb(kk)))) then
            at_bound = -1
        else if (abs(val-row_ub(kk)) <= me%active_tol*max(1.0_wp, abs(row_ub(kk)))) then
            at_bound = 1
        else
            at_bound = 0
        end if
        end function at_bound

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
