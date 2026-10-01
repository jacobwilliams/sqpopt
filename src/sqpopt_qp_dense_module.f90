!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  **Dense** active-set QP solver for the linearized SQP subproblem
!  (the default for small problems, see `sqpopt_qp_auto`). Forms a dense
!  Jacobian and a dense Hessian from `sqpopt`'s usual sparse/matrix-free
!  representations each time it is called, then solves the QP with
!  classical primal active-set machinery (Nocedal & Wright, *Numerical
!  Optimization*, Ch. 16.5): a null-space/reduced-Hessian step within the
!  current working set (Householder QR + Cholesky, see
!  [[sqpopt_dense_linalg_module]]), a ratio test to add a newly-binding
!  bound/constraint to the working set, and a Lagrange-multiplier sign
!  check to drop one when the current working set's minimizer has been
!  reached. Since it forms `O(n^2)`/`O(mn)` dense arrays every call, it is
!  only suitable for small-to-moderate `n`/`m`.
!
!  General constraints and variable bounds are treated uniformly as
!  "rows" of a single constraint matrix, each two-sided
!  (`row_lb<=row^T u<=row_ub`); equality rows (`row_lb==row_ub`) never
!  leave the working set once in it.
!
!  **Elastic mode.** A primal active-set method needs a feasible starting
!  point. Rather than a separate phase-1 problem, every general constraint
!  that is violated at `p=0` is made *elastic*: it gets a nonnegative slack
!  `s_i` (entering the row with the sign that can restore feasibility)
!  with an \( \ell_1 \) penalty \( \rho s_i \) in the objective (SNOPT's
!  elastic mode, Fletcher's S\( \ell_1 \)QP). `p=0` plus those slacks is
!  then feasible, and all other rows are kept hard. For a large enough
!  weight \( \rho \) (more than the largest multiplier), the slacks are zero
!  at the solution whenever the linearized constraints are consistent,
!  and the solution is that of the original QP. If slacks are still
!  positive at the solution, \( \rho \) is increased (warm-starting from
!  the current working set) up to `elastic_weight_max`; if they are still
!  positive then, the linearized constraints are inconsistent, and
!  `istat=sqpopt_infeasible` is returned, with `p` the elastic solution
!  (which minimizes the \( \ell_1 \) violation of the linearized
!  constraints, then the QP objective).
!
!  **Forced elastic mode** (`force_sign`, `force_weight`, set for one solve
!  by [[sqpopt_qp_solver_module]] when the SQP iteration re-solves a QP
!  with diverging-multiplier constraints elastic, see
!  `options%elastic_multiplier_limit`): the chosen rows get an elastic
!  slack in the given direction even if they aren't violated at the
!  starting step, every slack has the fixed weight `force_weight` (not
!  raised), and positive slacks at the solution are accepted
!  (`istat=sqpopt_success`): the solution is that of the \( \ell_1 \)
!  penalty QP, whose multipliers on those rows are at most the weight.
!
!  **Crash and warm start.** Rather than from `p=0`, the iterations start
!  from the minimum-norm step that satisfies an initial guess of the
!  working set: the previous solve's final working set (a *warm start*,
!  when `warm_start` is on and the problem size is unchanged -- near an SQP
!  solution the active set settles, and the QP then finishes in one or two
!  iterations), or else the equality constraints and fixed variables (a
!  *crash* start, so equality constraints don't each need an elastic slack
!  and an iteration to remove it). Any variable bounds this step violates
!  are added to the guess and the step recomputed (a few rounds), then it
!  is clipped to the bounds; only the rows it still violates get elastic
!  slacks. The guess only sets the starting point: optimality is still
!  established by the active-set iterations.
!
!  **Robustness.** Rows are only ever added to the working set if they are
!  linearly independent of it (the initial working set is built with a
!  Gram-Schmidt independence check; a row blocking a step along a
!  null-space direction is independent by construction), so the null-space
!  basis is always well defined. If the reduced Hessian \( Z^THZ \) is not
!  positive definite (e.g. an indefinite SR1 approximation, or zero
!  curvature along an elastic slack), a direction of nonpositive curvature
!  is followed to the nearest blocking row instead of taking a (huge or
!  undefined) Newton step. Tolerances are relative to the problem scale.

    module sqpopt_qp_dense_module

    use sqpopt_kinds,              only: wp => sqpopt_module_wp
    use sqpopt_types_module,       only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_qp_solve_failed, &
                                         sqpopt_infeasible, sqpopt_infinity
    use sqpopt_hessian_module,     only: sqpopt_hessian_type
    use sqpopt_dense_linalg_module, only: dense_null_space, dense_modified_cholesky, dense_solve_cholesky, &
                                          dense_cholesky_curvature

    implicit none

    private

    type, public :: sqpopt_dense_qp_type
        !! options for the dense active-set QP solver.

        integer  :: max_iter   = 100        !! minimum limit on the number of active-set iterations per QP solve
                                            !! (the actual limit is `max(max_iter, 10*(number of rows+1))`)
        real(wp) :: active_tol = 1.0e-8_wp  !! relative tolerance for a row being at a bound, and for the sign of
                                            !! a multiplier (relative to the largest multiplier)
        real(wp) :: opt_tol    = 1.0e-10_wp !! relative tolerance on the reduced-gradient stationarity test
                                            !! (relative to \( 1+\lVert Hp+g \rVert_\infty \))
        real(wp) :: feas_tol   = 1.0e-6_wp  !! an elastic slack larger than `feas_tol*max(1,|bound|)` at the
                                            !! solution counts as a violated linearized constraint
        real(wp) :: elastic_weight     = 1.0e4_wp  !! initial elastic penalty weight \( \rho \), relative to
                                                   !! \( \max(1,\lVert g \rVert_\infty) \)
        real(wp) :: elastic_weight_max = 1.0e10_wp !! largest elastic penalty weight tried (same scaling) before the
                                                   !! linearized constraints are declared inconsistent
        logical  :: warm_start = .true.            !! start from the previous solve's final working set, if the
                                                   !! problem size is unchanged (see the module-level documentation)

        ! forced elastic mode (internal inputs, set for one solve by [[sqpopt_qp_solver_module]], see
        ! `solve_qp_subproblem`): the rows with a nonzero `force_sign` are made elastic even if not violated
        ! at the starting step, in that direction (`+1`: the slack relaxes the row's lower bound, `-1` its
        ! upper bound), and every slack has the fixed weight `force_weight` (not raised); positive slacks at
        ! the solution are then accepted (`istat=sqpopt_success`), since they are the point of the solve
        integer,  dimension(:), allocatable :: force_sign
        real(wp) :: force_weight = 0.0_wp

        integer :: n_iter = 0 !! number of active-set iterations taken by the last solve (output)
        integer :: n_working = 0 !! number of general rows and variable bounds in the final working set of the
                                 !! last solve (output)
        integer :: n_slacks  = 0 !! number of elastic slacks in the last solve (output; see the module docs)
        logical :: negative_curvature = .false. !! whether the last solve found a direction of negative curvature
                                                !! of the Hessian in the variables (output; the QP was nonconvex)

        ! internal state (the working set at the end of the previous solve, for warm starts):
        integer, dimension(:), allocatable :: warm_status !! side (-1/0/+1) of each general row and variable bound

        contains

        procedure, public :: solve => solve_dense_qp

    end type sqpopt_dense_qp_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda` (see the module-level
!  documentation). `istat` is `sqpopt_success`, `sqpopt_infeasible` (the
!  linearized constraints are inconsistent; `p` is the elastic solution),
!  or `sqpopt_qp_solve_failed` (iteration limit, or unbounded along a
!  direction of nonpositive curvature; `p` is the last iterate).

    subroutine solve_dense_qp(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    class(sqpopt_dense_qp_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian    !! matrix-free Hessian approximation (densified here)
    type(sqpopt_sparse_matrix), intent(in)    :: jac        !! sparse constraint Jacobian, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x          !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g          !! objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c          !! constraint values `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: x_lb, x_ub !! variable bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c_lb, c_ub !! constraint bounds `dimension(m)`
    real(wp), dimension(:),     intent(out)   :: p          !! search direction `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: lambda     !! Lagrange multiplier estimate `dimension(m)`
    integer,                    intent(out)   :: istat      !! status code (see [[sqpopt_types_module]])

    integer :: n, m, nv, nt, mtot, k, i, it, n_z, n_active, maxit
    real(wp), dimension(:,:), allocatable :: h, arows, ja, z, jd
    real(wp), dimension(:),   allocatable :: row_lb, row_ub, u, gext, hu_g, rg, dvec, coeff, s_sign, rhs_active, s0
    real(wp), dimension(size(g)) :: p0
    real(wp), dimension(:),   allocatable :: jp0
    integer,  dimension(:),   allocatable :: status, orig_idx, coeff_idx, slack_row
    logical,  dimension(:),   allocatable :: is_equality
    real(wp) :: rho, rho_max, gscale, alpha, alpha_cap, scale
    integer  :: blocking, blocking_side
    logical  :: pd, at_face_optimum

    n = size(g)
    m = size(c)

    me%negative_curvature = .false.
    me%n_working = 0
    me%n_slacks  = 0

    ! ---- the dense Jacobian ----
    allocate(jd(m,n)); jd = 0.0_wp
    do k = 1, jac%nnz
        jd(jac%irow(k), jac%icol(k)) = jd(jac%irow(k), jac%icol(k)) + jac%val(k)
    end do

    ! ---- starting step: crash or warm start (see the module docs) ----
    call starting_step(p0)

    ! ---- elastic slacks: one for each general row violated at p0 ----
    allocate(row_lb(m), row_ub(m))
    row_lb = c_lb - c
    row_ub = c_ub - c
    allocate(s_sign(m), slack_row(m), s0(m))
    jp0 = matmul(jd, p0)
    nv = 0
    do i = 1, m
        s_sign(i) = 0.0_wp
        if (jp0(i) < row_lb(i) - me%feas_tol*max(1.0_wp, abs(row_lb(i)))) then
            s_sign(i) = 1.0_wp       ! J_i p + s >= row_lb, i.e. s makes up the shortfall
        else if (jp0(i) > row_ub(i) + me%feas_tol*max(1.0_wp, abs(row_ub(i)))) then
            s_sign(i) = -1.0_wp      ! J_i p - s <= row_ub
        end if
        if (s_sign(i) == 0.0_wp .and. forced()) then
            if (me%force_sign(i) /= 0) s_sign(i) = real(sign(1, me%force_sign(i)), wp)
        end if
        if (s_sign(i) /= 0.0_wp) then
            nv = nv + 1
            slack_row(nv) = i
            s0(nv) = max(0.0_wp, merge(row_lb(i) - jp0(i), jp0(i) - row_ub(i), s_sign(i) > 0.0_wp))
        end if
    end do
    s0 = s0(1:nv)
    nt   = n + nv       ! unknowns: p, then the slacks
    mtot = m + nt       ! rows: general constraints, then bounds on every unknown

    ! ---- the combined constraint rows ----
    allocate(arows(mtot,nt)); arows = 0.0_wp
    arows(1:m,1:n) = jd
    do k = 1, nv
        arows(slack_row(k), n+k) = s_sign(slack_row(k))
    end do
    do k = 1, nt
        arows(m+k,k) = 1.0_wp
    end do
    row_lb = [row_lb, x_lb - x, spread(0.0_wp, 1, nv)]
    row_ub = [row_ub, x_ub - x, spread(sqpopt_infinity, 1, nv)]

    allocate(is_equality(mtot))
    do k = 1, mtot
        is_equality(k) = row_ub(k)-row_lb(k) <= me%active_tol*max(1.0_wp, abs(row_lb(k)))
    end do

    ! ---- the dense Hessian (zero in the slack directions) and gradient ----
    allocate(h(n,n))
    call hessian%dense(h)

    gscale  = 1.0_wp
    if (n > 0) gscale = max(1.0_wp, maxval(abs(g)))
    rho     = me%elastic_weight*gscale
    rho_max = me%elastic_weight_max*gscale
    if (forced()) then
        rho     = me%force_weight
        rho_max = me%force_weight
    end if
    allocate(gext(nt))
    gext(1:n)    = g
    gext(n+1:nt) = rho

    ! ---- feasible starting point: p0, with the slacks just large enough ----
    allocate(u(nt))
    u(1:n)    = p0
    u(n+1:nt) = s0

    ! ---- initial working set: rows at a bound at u, if linearly independent ----
    allocate(status(mtot)); status = 0
    call initial_working_set()

    ! ---- active-set iterations ----
    istat = sqpopt_qp_solve_failed
    allocate(coeff(0), coeff_idx(0))
    maxit = max(me%max_iter, 10*(mtot+1))

    do it = 1, maxit

        call build_active_set(arows, row_lb, row_ub, status, mtot, nt, ja, rhs_active, orig_idx, n_active)
        call dense_null_space(ja, n_active, nt, z, n_z)

        hu_g = gradient(u)
        scale = 1.0_wp
        if (n > 0) scale = 1.0_wp + maxval(abs(hu_g(1:n)))
        at_face_optimum = .false.

        if (n_z <= 0) then

            at_face_optimum = .true.

        else

            rg = matmul(transpose(z), hu_g)
            block
                real(wp), dimension(n_z,n_z) :: zthz, l_fac
                real(wp), dimension(n_z)     :: dz, dcurv
                zthz = matmul(transpose(z(1:n,:)), matmul(h, z(1:n,:)))
                call dense_cholesky_curvature(zthz, n_z, l_fac, pd, dcurv)
                if (pd) then
                    if (norm2(rg) <= me%opt_tol*scale) then
                        at_face_optimum = .true.
                    else
                        ! Newton step to the minimizer on this face:
                        call dense_solve_cholesky(l_fac, n_z, -rg, dz)
                        dvec = matmul(z, dz)
                        alpha_cap = 1.0_wp
                    end if
                else
                    ! nonpositive curvature on this face: move along it, downhill,
                    ! to the nearest blocking row (the step is not bounded by 1).
                    ! (Negative curvature in the variables -- not just zero
                    ! curvature, e.g. along an elastic slack -- is reported.)
                    block
                        real(wp), dimension(n) :: dp
                        dp = matmul(z(1:n,:), dcurv)
                        if (dot_product(dp, matmul(h, dp)) < &
                            -1.0e-8_wp*max(1.0_wp, maxval(abs(h)))*dot_product(dp, dp)) me%negative_curvature = .true.
                    end block
                    if (dot_product(rg, dcurv) > 0.0_wp) dcurv = -dcurv
                    dvec = matmul(z, dcurv)
                    alpha_cap = huge(1.0_wp)
                end if
            end block

            if (.not. at_face_optimum) then
                call ratio_test(dvec, alpha_cap, alpha, blocking, blocking_side)
                if (blocking == 0 .and. .not. pd) then
                    ! no row blocks the nonpositive-curvature direction; if it is
                    ! flat (no descent), try the other way before giving up:
                    if (abs(dot_product(hu_g, dvec)) <= me%opt_tol*scale*norm2(dvec)) then
                        dvec = -dvec
                        call ratio_test(dvec, alpha_cap, alpha, blocking, blocking_side)
                    end if
                    if (blocking == 0) exit  ! unbounded QP
                end if
                u = u + alpha*dvec
                if (blocking /= 0 .and. alpha < alpha_cap - 1.0e-12_wp) then
                    status(blocking) = blocking_side
                    cycle  ! enlarge the working set and re-derive the null space
                end if
                at_face_optimum = .true.  ! took the full Newton step
            end if

        end if

        if (at_face_optimum) then

            hu_g = gradient(u)

            ! multipliers for the working set: least-squares solve of ja^T*coeff = H*u+g
            if (n_active > 0) then
                block
                    real(wp), dimension(n_active,n_active) :: gram, l_fac
                    real(wp), dimension(n_active) :: coeff_local, correction
                    gram = matmul(ja, transpose(ja))
                    call dense_modified_cholesky(gram, n_active, l_fac)
                    call dense_solve_cholesky(l_fac, n_active, matmul(ja, hu_g), coeff_local)
                    ! one step of iterative refinement (the normal equations square
                    ! the condition number of `ja`, which costs accuracy when rows
                    ! are nearly dependent):
                    call dense_solve_cholesky(l_fac, n_active, &
                                              matmul(ja, hu_g - matmul(transpose(ja), coeff_local)), correction)
                    coeff_local = coeff_local + correction
                    coeff     = coeff_local
                    coeff_idx = orig_idx
                end block
            else
                coeff     = [real(wp) ::]
                coeff_idx = [integer ::]
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
            if (any(u(n+1:nt) > me%feas_tol*max(1.0_wp, s0)) .and. .not. forced()) then
                if (rho < rho_max) then
                    ! not yet known to be inconsistent: raise the weight and continue
                    ! from the current point and working set:
                    rho = min(100.0_wp*rho, rho_max)
                    gext(n+1:nt) = rho
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
    me%n_working = count(status(1:m+n) /= 0)
    me%n_slacks  = nv
    p = u(1:n)
    lambda = 0.0_wp
    do k = 1, size(coeff_idx)
        if (coeff_idx(k) <= m) lambda(coeff_idx(k)) = coeff(k)
    end do

    ! remember the final working set (general rows and variable bounds) for a warm start:
    me%warm_status = status(1:m+n)

    contains

        pure logical function forced()
        !! whether this solve is in forced elastic mode (see `force_sign`)
        forced = me%force_weight > 0.0_wp .and. allocated(me%force_sign)
        if (forced) forced = size(me%force_sign) == m
        end function forced

        subroutine starting_step(p0)
        !! the minimum-norm step satisfying the initial working-set guess (the
        !! previous solve's final working set, or the equality constraints and
        !! fixed variables), with any violated variable bounds added to the
        !! guess (up to 4 rounds), then clipped to the bounds
        real(wp), dimension(n), intent(out) :: p0 !! the starting step `dimension(n)`
        real(wp), dimension(n) :: blb, bub
        integer,  dimension(m+n) :: guess   ! side (-1/+1, 0 = not in the guess) of each general row / bound
        real(wp), dimension(n, n) :: qb     ! orthonormal basis of the selected rows
        real(wp), dimension(n, n) :: t      ! lower-triangular coefficients: row_j = sum_i t(j,i) qb(:,i)
        real(wp), dimension(n) :: a, r, y, bsel
        real(wp) :: target
        integer :: kk, nb, round, j, jj
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

        do round = 1, 4
            ! select a linearly independent subset of the guessed rows (modified
            ! Gram-Schmidt, with reorthogonalization), recording the coefficients:
            nb = 0
            do kk = 1, m+n
                if (guess(kk) == 0 .or. nb >= n) cycle
                if (kk <= m) then
                    a = jd(kk,:)
                    target = merge(c_lb(kk)-c(kk), c_ub(kk)-c(kk), guess(kk) == -1)
                else
                    a = 0.0_wp; a(kk-m) = 1.0_wp
                    target = merge(blb(kk-m), bub(kk-m), guess(kk) == -1)
                end if
                if (abs(target) >= sqpopt_infinity) cycle
                r = a
                t(nb+1,:) = 0.0_wp
                do jj = 1, 2
                    do j = 1, nb
                        y(1) = dot_product(qb(:,j), r)
                        t(nb+1,j) = t(nb+1,j) + y(1)
                        r = r - y(1)*qb(:,j)
                    end do
                end do
                if (norm2(r) > 1.0e-10_wp*norm2(a)) then
                    nb = nb + 1
                    t(nb,nb) = norm2(r)
                    qb(:,nb) = r/t(nb,nb)
                    bsel(nb) = target
                end if
            end do
            ! minimum-norm solution of (selected rows)*p0 = targets: p0 = qb*y,
            ! with t*y = targets (forward substitution):
            do j = 1, nb
                y(j) = (bsel(j) - dot_product(t(j,1:j-1), y(1:j-1)))/t(j,j)
            end do
            p0 = 0.0_wp
            if (nb > 0) p0 = matmul(qb(:,1:nb), y(1:nb))
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

        function gradient(v) result(gr)
        !! the gradient of the (elastic) QP objective at `v`: `H*v_p + g`, then `rho` for each slack
        real(wp), dimension(:), intent(in) :: v !! the unknowns: the step, then the elastic slacks `dimension(nt)`
        real(wp), dimension(size(v)) :: gr
        gr(1:n)    = matmul(h, v(1:n)) + g
        gr(n+1:nt) = rho
        end function gradient

        subroutine ratio_test(d, alpha_cap, alpha, blocking, blocking_side)
        !! the largest `alpha <= alpha_cap` for which `u+alpha*d` satisfies every
        !! row not in the working set, and the row (and side) that blocks first
        real(wp), dimension(:), intent(in)  :: d             !! search direction in the unknowns `dimension(nt)`
        real(wp),               intent(in)  :: alpha_cap     !! largest step length to consider
        real(wp),               intent(out) :: alpha         !! step length
        integer,                intent(out) :: blocking      !! the blocking row (`0` if none blocks before `alpha_cap`)
        integer,                intent(out) :: blocking_side !! its bound: `-1` lower, `+1` upper
        real(wp) :: rate, alpha_k, val, dnorm
        integer :: kk
        alpha = alpha_cap
        blocking = 0
        blocking_side = 0
        dnorm = norm2(d)
        do kk = 1, mtot
            if (status(kk) /= 0) cycle
            rate = dot_product(arows(kk,:), d)
            if (abs(rate) <= 1.0e-12_wp*norm2(arows(kk,:))*dnorm) cycle  ! (numerically) parallel to the row
            val = dot_product(arows(kk,:), u)
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
        !! working set, skipping any that are linearly dependent on those
        !! already added (modified Gram-Schmidt, with reorthogonalization)
        real(wp), dimension(nt, min(mtot,nt)) :: qb
        real(wp), dimension(nt) :: r
        real(wp) :: val
        integer :: kk, pass, nb, side
        nb = 0
        do pass = 1, 2
            do kk = 1, mtot
                if (nb >= nt) return
                if (status(kk) /= 0) cycle
                if ((pass == 1) .neqv. is_equality(kk)) cycle
                val = dot_product(arows(kk,:), u)
                if (abs(val-row_lb(kk)) <= me%active_tol*max(1.0_wp, abs(row_lb(kk)))) then
                    side = -1
                else if (abs(val-row_ub(kk)) <= me%active_tol*max(1.0_wp, abs(row_ub(kk)))) then
                    side = 1
                else
                    cycle
                end if
                r = arows(kk,:)
                if (nb > 0) then
                    r = r - matmul(qb(:,1:nb), matmul(r, qb(:,1:nb)))
                    r = r - matmul(qb(:,1:nb), matmul(r, qb(:,1:nb)))
                end if
                if (norm2(r) > 1.0e-10_wp*norm2(arows(kk,:))) then
                    nb = nb + 1
                    qb(:,nb) = r/norm2(r)
                    status(kk) = side
                end if
            end do
        end do
        end subroutine initial_working_set

    end subroutine solve_dense_qp
!*******************************************************************************

!*******************************************************************************
!>
!  gather the currently-active rows (`status/=0`) of the combined
!  `mtot x n` row set into a dense `n_active x n` matrix `ja` and their
!  target right-hand-side values `rhs_active` (`row_lb` or `row_ub`,
!  whichever side is active), along with `orig_idx`, mapping each active
!  row back to its index in `1..mtot`.

    subroutine build_active_set(arows, row_lb, row_ub, status, mtot, n, ja, rhs_active, orig_idx, n_active)

    integer,                            intent(in)  :: mtot !! total number of rows in the combined constraint matrix
    integer,                            intent(in)  :: n !! number of columns in the constraint matrix
    real(wp), dimension(mtot,n),        intent(in)  :: arows !! combined constraint matrix (mtot rows, n columns)
    real(wp), dimension(mtot),          intent(in)  :: row_lb !! lower bounds for each row in the combined constraint matrix
    real(wp), dimension(mtot),          intent(in)  :: row_ub !! upper bounds for each row in the combined constraint matrix
    integer,  dimension(mtot),          intent(in)  :: status !! status of each row (0 if inactive, -1 if active at lower bound, 1 if active at upper bound)
    real(wp), dimension(:,:), allocatable, intent(out) :: ja !! dense matrix of currently-active rows
    real(wp), dimension(:),   allocatable, intent(out) :: rhs_active !! right-hand-side values for the active rows
    integer,  dimension(:),   allocatable, intent(out) :: orig_idx !! mapping from active rows to their original indices in 1..mtot
    integer,                            intent(out) :: n_active !! number of currently-active rows

    integer :: k, idx

    n_active = count(status /= 0)
    allocate(ja(n_active,n), rhs_active(n_active), orig_idx(n_active))

    idx = 0
    do k = 1, mtot
        if (status(k) /= 0) then
            idx = idx + 1
            ja(idx,:)      = arows(k,:)
            rhs_active(idx) = merge(row_lb(k), row_ub(k), status(k) == -1)
            orig_idx(idx)  = k
        end if
    end do

    end subroutine build_active_set
!*******************************************************************************

    end module sqpopt_qp_dense_module
!*******************************************************************************
