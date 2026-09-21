!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Opt-in **dense** active-set QP solver for the linearized SQP
!  subproblem. Forms a dense Jacobian and a
!  dense Hessian from `sqpopt`'s usual sparse/matrix-free representations
!  each time it is called, then solves the QP with classical dense
!  active-set machinery (Nocedal & Wright, *Numerical Optimization*, Ch.
!  16.3): a null-space/reduced-Hessian step within the current working
!  set (Householder QR + modified Cholesky, see
!  [[sqpopt_dense_linalg_module]]), a ratio test to add a newly-binding
!  bound/constraint to the working set, and a Lagrange-multiplier sign
!  check to drop one when the current working set's exact minimizer has
!  been reached.
!
!  Unlike `sqpopt_qp_solver_module`'s default composite-step heuristic,
!  this
!  enforces the linearized constraints and bounds **exactly** within the
!  QP itself. Unlike a general large-scale sparse QP, it is only suitable
!  for small-to-moderate `n`/`m`, since it forms `O(n^2)`/`O(mn)` dense
!  arrays every call -- this is why it is a separate, explicitly opt-in
!  mode (`sqpopt_qp_dense`) rather than the default.
!
!  General constraints and variable bounds are treated uniformly as
!  `m+n` candidate "rows" (the `m` rows of the Jacobian, then `n` unit
!  rows for the variable bounds on `p`), each two-sided
!  (`row_lb<=row^T p<=row_ub`); equality rows (`row_lb==row_ub`) are
!  permanently active and never leave the working set.
!
!  @note Bounds and general constraints are both just two-sided rows
!  here, so there is no need to introduce a separate slack vector `s` for
!  the general constraints -- working directly in `p`-space (`n`
!  unknowns) is simpler and mathematically equivalent.

    module sqpopt_qp_dense_module

    use sqpopt_kinds,              only: wp => sqpopt_module_wp
    use sqpopt_types_module,       only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_qp_solve_failed
    use sqpopt_hessian_module,     only: sqpopt_hessian_type
    use sqpopt_dense_linalg_module, only: dense_null_space, dense_modified_cholesky, dense_solve_cholesky

    implicit none

    private

    type, public :: sqpopt_dense_qp_type
        !! options and workspace for the dense active-set QP solver.

        integer  :: max_iter   = 100       !! maximum number of active-set changes allowed per QP solve
        real(wp) :: active_tol = 1.0e-8_wp !! tolerance used to detect an (in)active/equality row
        real(wp) :: opt_tol    = 1.0e-8_wp !! tolerance on the reduced-gradient stationarity test

        contains

        procedure, public :: solve => solve_dense_qp

    end type sqpopt_dense_qp_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda`, using a dense active-set
!  method (see the module-level documentation).

    subroutine solve_dense_qp(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    class(sqpopt_dense_qp_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! matrix-free Hessian approximation (densified here)
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! sparse constraint Jacobian, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x       !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g       !! objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c       !! constraint values `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: x_lb, x_ub !! variable bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c_lb, c_ub !! constraint bounds `dimension(m)`
    real(wp), dimension(:),     intent(out)   :: p       !! search direction `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: lambda  !! Lagrange multiplier estimate `dimension(m)`
    integer,                    intent(out)   :: istat   !! status code (see [[sqpopt_types_module]])

    integer :: n, m, mtot, k, it, n_z, n_active
    real(wp), dimension(:,:), allocatable :: h, arows, jdense, z
    real(wp), dimension(:),   allocatable :: row_lb, row_ub, u, hu_g, rg, d, dz, coeff
    real(wp), dimension(:,:), allocatable :: ja, gram, l_fac
    real(wp), dimension(:),   allocatable :: rhs_active
    integer,  dimension(:),   allocatable :: orig_idx
    logical,  dimension(:),   allocatable :: is_equality
    integer,  dimension(:),   allocatable :: status
    real(wp) :: alpha, rate, alpha_k, worst, val
    integer  :: blocking, blocking_side, worst_idx, idx
    logical  :: at_face_optimum

    n = size(g)
    m = size(c)
    mtot = m + n

    ! ---- form the dense Jacobian from the sparse COO representation ----
    allocate(jdense(m,n)); jdense = 0.0_wp
    do k = 1, jac%nnz
        jdense(jac%irow(k), jac%icol(k)) = jdense(jac%irow(k), jac%icol(k)) + jac%val(k)
    end do

    ! ---- form the dense Hessian by densifying the matrix-free operator column by column ----
    allocate(h(n,n))
    block
        real(wp), dimension(n) :: e, he
        do k = 1, n
            e = 0.0_wp
            e(k) = 1.0_wp
            call hessian%hv_product(e, he)
            h(:,k) = he
        end do
    end block
    h = 0.5_wp*(h + transpose(h))  !! symmetrize away any tiny roundoff asymmetry

    ! ---- combine general constraints (rows 1..m) and variable bounds (rows m+1..m+n)
    !      into one set of two-sided "rows" on p, row_lb(k) <= row(k,:).p <= row_ub(k) ----
    allocate(arows(mtot,n)); arows = 0.0_wp
    arows(1:m,:) = jdense
    do k = 1, n
        arows(m+k,k) = 1.0_wp
    end do
    allocate(row_lb(mtot), row_ub(mtot))
    row_lb(1:m) = c_lb - c
    row_ub(1:m) = c_ub - c
    row_lb(m+1:mtot) = x_lb - x
    row_ub(m+1:mtot) = x_ub - x

    allocate(is_equality(mtot))
    is_equality(1:m)        = (c_ub-c_lb) <= me%active_tol
    is_equality(m+1:mtot)   = (x_ub-x_lb) <= me%active_tol

    allocate(status(mtot))
    status = 0
    where (is_equality) status = -1  !! permanently active; -1 vs +1 is arbitrary since row_lb==row_ub there

    allocate(u(n)); u = 0.0_wp

    ! ---- phase 1: bootstrap a feasible-for-the-initial-working-set starting point ----
    call project_onto_active(arows, row_lb, row_ub, status, mtot, n, u)
    do k = 1, n
        if (is_equality(m+k)) cycle
        if (u(k) < row_lb(m+k) - me%active_tol) then
            status(m+k) = -1
        else if (u(k) > row_ub(m+k) + me%active_tol) then
            status(m+k) = 1
        end if
    end do
    do k = 1, m
        if (is_equality(k)) cycle
        val = dot_product(arows(k,:), u)
        if (val < row_lb(k) - me%active_tol) then
            status(k) = -1
        else if (val > row_ub(k) + me%active_tol) then
            status(k) = 1
        end if
    end do
    call project_onto_active(arows, row_lb, row_ub, status, mtot, n, u)

    ! ---- phase 2: active-set iterations ----
    istat = sqpopt_qp_solve_failed
    allocate(coeff(0))  !! defined once optimality is reached; harmless placeholder until then

    do it = 1, max(me%max_iter, 10*(mtot+1))

        call build_active_set(arows, row_lb, row_ub, status, mtot, n, ja, rhs_active, orig_idx, n_active)
        call dense_null_space(ja, n_active, n, z, n_z)

        at_face_optimum = .false.

        if (n_z == 0) then

            at_face_optimum = .true.

        else

            if (allocated(hu_g)) deallocate(hu_g)
            allocate(hu_g(n)); hu_g = matmul(h,u) + g
            if (allocated(rg)) deallocate(rg)
            allocate(rg(n_z)); rg = matmul(transpose(z), hu_g)

            if (norm2(rg) <= me%opt_tol) then

                at_face_optimum = .true.

            else

                block
                    real(wp), dimension(n_z,n_z) :: zthz, l_fac_z
                    real(wp), dimension(n_z)     :: dz_local
                    zthz = matmul(transpose(z), matmul(h, z))
                    call dense_modified_cholesky(zthz, n_z, l_fac_z)
                    call dense_solve_cholesky(l_fac_z, n_z, -rg, dz_local)
                    if (allocated(d)) deallocate(d)
                    allocate(d(n))
                    d = matmul(z, dz_local)
                end block

                ! ---- ratio test: how far can we move along d before an inactive row binds? ----
                alpha = 1.0_wp
                blocking = 0
                blocking_side = 0
                do k = 1, mtot
                    if (status(k) /= 0) cycle
                    rate = dot_product(arows(k,:), d)
                    if (rate > me%active_tol) then
                        alpha_k = (row_ub(k) - dot_product(arows(k,:), u))/rate
                        if (alpha_k < alpha) then
                            alpha = max(alpha_k, 0.0_wp)
                            blocking = k
                            blocking_side = 1
                        end if
                    else if (rate < -me%active_tol) then
                        alpha_k = (row_lb(k) - dot_product(arows(k,:), u))/rate
                        if (alpha_k < alpha) then
                            alpha = max(alpha_k, 0.0_wp)
                            blocking = k
                            blocking_side = -1
                        end if
                    end if
                end do

                u = u + alpha*d

                if (blocking /= 0 .and. alpha < 1.0_wp - 1.0e-10_wp) then
                    status(blocking) = blocking_side
                    cycle  !! enlarge the working set and re-derive the null space
                else
                    at_face_optimum = .true.  !! took the full reduced-Newton step
                end if

            end if

        end if

        if (at_face_optimum) then

            ! ---- multiplier check: are all active inequality rows correctly signed? ----
            if (allocated(hu_g)) deallocate(hu_g)
            allocate(hu_g(n)); hu_g = matmul(h,u) + g

            if (n_active == 0) then
                istat = sqpopt_success
                exit
            end if

            block
                real(wp), dimension(n_active,n_active) :: gram_local, l_fac_a
                real(wp), dimension(n_active) :: coeff_local
                gram_local = matmul(ja, transpose(ja))
                call dense_modified_cholesky(gram_local, n_active, l_fac_a)
                call dense_solve_cholesky(l_fac_a, n_active, matmul(ja, hu_g), coeff_local)
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
                        worst = -coeff(idx)
                        worst_idx = k
                    end if
                else
                    if (coeff(idx) > worst) then
                        worst = coeff(idx)
                        worst_idx = k
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

    end subroutine solve_dense_qp
!*******************************************************************************

!*******************************************************************************
!>
!  gather the currently-active rows (`status/=0`) of the combined
!  `mtot x n` row set into a dense `n_active x n` matrix `ja` and their
!  target right-hand-side values `rhs_active` (`row_lb` or `row_ub`,
!  whichever side is active), along with `orig_idx`, mapping each active
!  row back to its index in `1..mtot` (`<=m` for a general constraint,
!  `>m` for a variable bound).

    subroutine build_active_set(arows, row_lb, row_ub, status, mtot, n, ja, rhs_active, orig_idx, n_active)

    integer,                            intent(in)  :: mtot, n
    real(wp), dimension(mtot,n),        intent(in)  :: arows
    real(wp), dimension(mtot),          intent(in)  :: row_lb, row_ub
    integer,  dimension(mtot),          intent(in)  :: status
    real(wp), dimension(:,:), allocatable, intent(out) :: ja
    real(wp), dimension(:),   allocatable, intent(out) :: rhs_active
    integer,  dimension(:),   allocatable, intent(out) :: orig_idx
    integer,                            intent(out) :: n_active

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

!*******************************************************************************
!>
!  adjust `u` (in place) by the minimum-norm correction needed so that
!  every currently-active row (`status/=0`) exactly satisfies its target
!  bound value -- used to bootstrap a feasible-for-the-working-set
!  starting point.

    subroutine project_onto_active(arows, row_lb, row_ub, status, mtot, n, u)

    integer,                     intent(in)    :: mtot, n
    real(wp), dimension(mtot,n), intent(in)    :: arows
    real(wp), dimension(mtot),   intent(in)    :: row_lb, row_ub
    integer,  dimension(mtot),   intent(in)    :: status
    real(wp), dimension(n),      intent(inout) :: u

    real(wp), dimension(:,:), allocatable :: ja
    real(wp), dimension(:),   allocatable :: rhs_active
    integer,  dimension(:),   allocatable :: orig_idx
    integer :: n_active

    call build_active_set(arows, row_lb, row_ub, status, mtot, n, ja, rhs_active, orig_idx, n_active)
    if (n_active == 0) return

    block
        real(wp), dimension(n_active,n_active) :: gram, l_fac
        real(wp), dimension(n_active) :: resid, y
        resid = rhs_active - matmul(ja, u)
        gram  = matmul(ja, transpose(ja))
        call dense_modified_cholesky(gram, n_active, l_fac)
        call dense_solve_cholesky(l_fac, n_active, resid, y)
        u = u + matmul(transpose(ja), y)
    end block

    end subroutine project_onto_active
!*******************************************************************************

    end module sqpopt_qp_dense_module
!*******************************************************************************
