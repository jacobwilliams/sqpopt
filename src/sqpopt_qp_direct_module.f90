!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  A direct solution of the QP subproblem (`options%direct_qp`), by sparse
!  factorizations of the KKT matrix of its working set (see
!  [[sqpopt_kkt_module]]), so it is only available in a build with MUMPS.
!
!  For a given working set (the general rows and variable bounds held at
!  a bound), the QP is an equality-constrained one, and its solution is one
!  linear solve:
!
!  $$ \begin{bmatrix} H & J_a^T \\ J_a & 0 \end{bmatrix}
!     \begin{bmatrix} p \\ -\lambda \end{bmatrix} =
!     \begin{bmatrix} -g \\ b_a - c_a \end{bmatrix} $$
!
!  with the variables at a bound of the working set fixed there. That step
!  solves the whole QP if it is *feasible* (it satisfies the rows and bounds
!  outside the working set) and *optimal* (the multipliers of the working
!  set's inequalities have the right signs), and the face is convex (no
!  negative curvature, which the factorization tells).
!
!  [[direct_qp_step]] starts from the working set that the QP solver would
!  start from: near a solution of the nonlinear problem, that is already
!  the right one, and the QP costs one factorization and one solve. If the
!  step isn't feasible and optimal, the working set is changed as in a
!  *primal-dual active-set method* (Hintermüller, Ito & Kunisch, 2002):
!  every row and bound the step violates is added, every one whose
!  multiplier has the wrong sign is dropped, and the new face is factored
!  and solved. After a few changes without success, it gives up, and the
!  caller solves the QP with its active-set solver instead (see
!  [[solve_qp_subproblem]]). So the direct method only ever replaces a QP
!  solve by a solution that satisfies the QP's optimality conditions.
!
!  Adding every violated bound at once can leave a row of the working set
!  with no free variable, or otherwise dependent on the others: the KKT
!  matrix is then singular. Such a face is factored again with a small
!  regularization (\( -\epsilon I \) in the rows' block, see
!  [[sqpopt_kkt_module]]), which satisfies the dependent rows only in a
!  heavily weighted least-squares sense. That gives them large multipliers,
!  so the bounds that over-determine them are dropped at the next change.
!  A step is only accepted if it satisfies the working set's rows, so a
!  regularized face is accepted only if its rows were consistent. A
!  regularized solve leaves each row short by \( \epsilon \) times its
!  multiplier, so if that is all that is wrong with the step, the face is
!  solved once more with a much smaller \( \epsilon \). (That matters
!  when the solver finds a matrix singular that isn't: near the solution of
!  a chain of 100,000 circle constraints, the first regularization left the
!  rows `7e-8` short, and the method gave up on a QP it had solved.) A face
!  that is still singular (a direction of zero curvature) ends the method.
!
!  A face with negative curvature has no minimizer. With inertia control
!  (see [[sqpopt_inertia_module]]), the Hessian's shift is raised until the
!  face is convex, and the method continues from there with the shifted
!  Hessian (the shift only ever increases, so this can't cycle). Without
!  it, the method gives up at such a face.
!
!  The factorization can't be updated when the working set changes: every
!  change is a new factorization. That is cheap next to the active-set
!  solvers on large problems (their iterations each solve with the reduced
!  Hessian), but it is why this isn't used for a QP that needs many
!  changes, such as the first one from an infeasible starting point.

    module sqpopt_qp_direct_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_all_finite
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_kkt_module,     only: sqpopt_kkt_type
    use sqpopt_inertia_module, only: sqpopt_inertia_type
    use sqpopt_linalg_module,  only: sparse_matvec, sparse_matvec_transpose

    implicit none

    private

    public :: direct_qp_step, direct_outcome_text

    real(wp), parameter :: stationarity_tol = 1.0e-8_wp !! tolerance on the stationarity residual of a solve, in
                                                        !! the free variables, relative to the terms that it is
                                                        !! the difference of (the gradient, and the Hessian times
                                                        !! the step)
    real(wp), parameter :: reg_factor = 1.0e-8_wp  !! the regularization of a singular face, relative to the square
                                                   !! of the largest element of the face's part of the Jacobian (the
                                                   !! working set's rows, in the variables that are not fixed)
    real(wp), parameter :: reg_factor_fine = 1.0e-11_wp !! the smaller one, tried if the step of a regularized face
                                                   !! only fails to satisfy the working set's rows: a regularized
                                                   !! solve leaves each row short by the regularization times its
                                                   !! multiplier

    ! how [[direct_qp_step]] ended:
    integer, parameter, public :: sqpopt_direct_solved     = 0 !! it found the QP's solution
    integer, parameter, public :: sqpopt_direct_nonconvex  = 1 !! a face had negative curvature
    integer, parameter, public :: sqpopt_direct_singular   = 2 !! a face had no unique minimizer, or
                                                               !! inconsistent rows
    integer, parameter, public :: sqpopt_direct_max_changes = 3 !! too many changes of the working set
    integer, parameter, public :: sqpopt_direct_failed     = 4 !! a factorization or a solve failed
    integer, parameter, public :: sqpopt_direct_inaccurate = 5 !! a solve was too inaccurate

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  try to solve the QP subproblem directly (see the module documentation),
!  starting from the working set `status`. If it succeeds
!  (`outcome = sqpopt_direct_solved`), `p` and `lambda` are the QP's
!  solution and `status` its working set. Otherwise they are not to be
!  used, `status` is unchanged, and `outcome` says why it gave up (see the
!  `sqpopt_direct_*` constants).

    subroutine direct_qp_step(kkt, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, max_changes, tol, &
                              status, p, lambda, n_changes, outcome, inertia)

    type(sqpopt_kkt_type),      intent(inout) :: kkt     !! the KKT matrix (see [[sqpopt_kkt_module]])
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! the Hessian (exact, with its shift, or quasi-Newton)
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! sparse constraint Jacobian, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x       !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g       !! objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c       !! current constraint values `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: x_lb    !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: x_ub    !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c_lb    !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: c_ub    !! constraint upper bounds `dimension(m)`
    integer,                    intent(in)    :: max_changes !! maximum number of changes of the working set
    real(wp),                   intent(in)    :: tol     !! relative tolerance for a violated row or bound, and for
                                                         !! the sign of a multiplier
    integer, dimension(:),      intent(inout) :: status  !! the working set: the side (`-1` lower, `+1` upper, `0`
                                                         !! not in it) of each general row, then of each variable
                                                         !! bound `dimension(m+n)`
    real(wp), dimension(:),     intent(out)   :: p       !! the step `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: lambda  !! the constraint multipliers `dimension(m)`
    integer,                    intent(out)   :: n_changes !! number of changes of the working set that were made
    integer,                    intent(out)   :: outcome   !! how it ended (a `sqpopt_direct_*` constant)
    type(sqpopt_inertia_type), optional, intent(inout) :: inertia !! the inertia control, to raise the Hessian's shift
                                                                  !! at a nonconvex face (used if `inertia%enabled`)

    integer :: n, m, i, j, k
    integer, dimension(:), allocatable :: st, st_new
    real(wp), dimension(:), allocatable :: v
    real(wp), dimension(:), allocatable :: p_fixed, hp, z
    real(wp), dimension(:), allocatable :: jp
    real(wp) :: dual_tol, lin
    logical  :: ok, convex, shifted, consistent
    logical  :: regularized !! whether the current face's KKT matrix is singular, and so regularized
    logical  :: fine        !! whether the current face is being solved again with the smaller regularization
    real(wp) :: reg, target, jmax

    allocate(st(size(status)), st_new(size(status)))
    allocate(v(size(x) + size(c)), p_fixed(size(x)), hp(size(x)), z(size(x)), jp(size(c)))

    n = size(x)
    m = size(c)
    outcome   = sqpopt_direct_failed
    n_changes = 0
    p      = 0.0_wp
    lambda = 0.0_wp
    st = status
    fine = .false.
    regularized = .false.

    do

        if (.not. fine) then
            call kkt%factor(hessian, jac, st, ok)
            if (.not. ok) return
            if (kkt%n_negative > 0) then
                ! a nonconvex face: raise the shift until it is convex, if we may
                convex = .false.
                if (present(inertia)) call inertia%correct(kkt, hessian, jac, st, shifted, convex)
                if (.not. convex) then
                    outcome = sqpopt_direct_nonconvex
                    return
                end if
            end if
            regularized = kkt%singular
        end if
        if (regularized) then
            ! dependent rows: regularize them (see the module documentation),
            ! relative to the elements of the Jacobian that are in the face's matrix
            reg = merge(reg_factor_fine, reg_factor, fine)
            jmax = 0.0_wp
            do k = 1, jac%nnz
                if (st(jac%irow(k)) /= 0 .and. st(m + jac%icol(k)) == 0) jmax = max(jmax, abs(jac%val(k)))
            end do
            if (jmax > 0.0_wp) reg = max(reg*jmax**2, tiny(1.0_wp))   ! (not zero, if the square underflows)
            call kkt%factor(hessian, jac, st, ok, reg=reg)
            if (.not. ok) return
            ! (a face without a unique minimizer is left to the active-set solver)
            if (kkt%singular .or. kkt%n_negative > 0) then
                outcome = sqpopt_direct_singular
                return
            end if
        end if

        ! the right-hand side: the variables at a bound of the working set are
        ! fixed there, and their part of the step moves to the right-hand side
        do j = 1, n
            select case (st(m+j))
            case (:-1);   p_fixed(j) = x_lb(j) - x(j)
            case (1:);    p_fixed(j) = x_ub(j) - x(j)
            case default; p_fixed(j) = 0.0_wp
            end select
        end do
        call hessian%hv_product(p_fixed, hp)
        call sparse_matvec(jac, p_fixed, jp)
        do j = 1, n
            v(j) = merge(-g(j) - hp(j), p_fixed(j), st(m+j) == 0)
        end do
        do i = 1, m
            select case (st(i))
            case (:-1);   v(n+i) = c_lb(i) - c(i) - jp(i)
            case (1:);    v(n+i) = c_ub(i) - c(i) - jp(i)
            case default; v(n+i) = 0.0_wp
            end select
        end do

        call kkt%solve(v, ok, hessian=hessian)
        if (.not. ok) return
        p = v(1:n)
        do i = 1, m
            lambda(i) = merge(-v(n+i), 0.0_wp, st(i) /= 0)
        end do

        ! the bound multipliers: z = g + H p - J^T lambda
        call hessian%hv_product(p, hp)
        call sparse_matvec_transpose(jac, lambda, z)
        z = g + hp - z
        call sparse_matvec(jac, p, jp)

        ! (the step must satisfy the stationarity condition in the free
        ! variables: this checks the accuracy of the solve. The residual is
        ! relative to the terms that it is the difference of: for a long step,
        ! H p is the sum of products much larger than itself, and its roundoff
        ! is relative to those)
        if (maxval(abs(z), mask=st(m+1:) == 0) > stationarity_tol*(1.0_wp + maxval(abs(g)) + maxval(abs(hp)) + &
                (hessian%magnitude() + hessian%shift)*maxval(abs(p)))) then
            outcome = sqpopt_direct_inaccurate
            return
        end if

        dual_tol = max(1.0_wp, maxval(abs(g)))
        if (m > 0) dual_tol = max(dual_tol, maxval(abs(lambda)))
        dual_tol = tol*dual_tol

        ! the new working set: add what the step violates, and drop what has a
        ! multiplier of the wrong sign (an equality, or a fixed variable, stays)
        st_new = st
        consistent = .true.
        do i = 1, m
            lin = c(i) + jp(i)
            if (st(i) /= 0) then
                ! (a regularized face may not satisfy its rows)
                target = merge(c_lb(i), c_ub(i), st(i) < 0)
                if (abs(lin - target) > tol*max(1.0_wp, abs(target))) consistent = .false.
            end if
            if (st(i) == 0) then
                if (lin < c_lb(i) - tol*max(1.0_wp, abs(c_lb(i)))) then
                    st_new(i) = -1
                else if (lin > c_ub(i) + tol*max(1.0_wp, abs(c_ub(i)))) then
                    st_new(i) = 1
                end if
            else if (c_ub(i) - c_lb(i) > 0.0_wp) then
                if ((st(i) < 0 .and. lambda(i) < -dual_tol) .or. (st(i) > 0 .and. lambda(i) > dual_tol)) st_new(i) = 0
            end if
        end do
        do j = 1, n
            if (st(m+j) == 0) then
                if (x(j) + p(j) < x_lb(j) - tol*max(1.0_wp, abs(x_lb(j)))) then
                    st_new(m+j) = -1
                else if (x(j) + p(j) > x_ub(j) + tol*max(1.0_wp, abs(x_ub(j)))) then
                    st_new(m+j) = 1
                end if
            else if (x_ub(j) - x_lb(j) > 0.0_wp) then
                if ((st(m+j) < 0 .and. z(j) < -dual_tol) .or. (st(m+j) > 0 .and. z(j) > dual_tol)) st_new(m+j) = 0
            end if
        end do

        if (all(st_new == st)) then
            if (.not. consistent .and. regularized .and. .not. fine) then
                ! (only the rows are not satisfied: perhaps because of the regularization)
                fine = .true.
                cycle
            else if (.not. consistent) then
                ! (the working set's rows are inconsistent, and nothing can be dropped)
                outcome = sqpopt_direct_singular
            else if (sqpopt_all_finite(p) .and. sqpopt_all_finite(lambda)) then
                ! feasible and optimal: the QP's solution
                outcome = sqpopt_direct_solved
                status  = st
            end if
            return
        end if
        if (n_changes >= max_changes) then
            outcome = sqpopt_direct_max_changes
            return
        end if
        n_changes = n_changes + 1
        st = st_new
        fine = .false.

    end do

    end subroutine direct_qp_step
!*******************************************************************************

!*******************************************************************************
!>
!  a short description of how [[direct_qp_step]] ended, for the detailed log.

    pure function direct_outcome_text(outcome) result(text)

    integer, intent(in) :: outcome !! a `sqpopt_direct_*` constant
    character(len=:), allocatable :: text

    select case (outcome)
    case (sqpopt_direct_solved);      text = 'solved'
    case (sqpopt_direct_nonconvex);   text = 'a face is nonconvex'
    case (sqpopt_direct_singular);    text = 'a face is singular'
    case (sqpopt_direct_max_changes); text = 'too many changes of the working set'
    case (sqpopt_direct_inaccurate);  text = 'a solve was inaccurate'
    case default;                     text = 'the factorization failed'
    end select

    end function direct_outcome_text
!*******************************************************************************

    end module sqpopt_qp_direct_module
!*******************************************************************************
