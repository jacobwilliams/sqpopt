!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The QP subproblem solved by DAQP (`options%qp_solver_mode = sqpopt_qp_daqp`):
!  the dual active-set solver for dense convex QPs of D. Arnström, A.
!  Bemporad, and D. Axehill (*A dual active-set solver for embedded quadratic
!  programming using recursive LDLᵀ updates*, IEEE Transactions on Automatic
!  Control 67(8), 2022), through its Fortran port, the `daqp-fortran` fpm
!  package (`daqp_module`).
!
!  Like the dense QP solver ([[sqpopt_qp_dense_module]]), it forms the dense
!  Hessian and Jacobian, and solves
!
!  $$ \min_p \; \tfrac12 p^T H p + g^T p \quad \text{s.t.} \quad
!     x_l - x \le p \le x_u - x, \quad c_l - c \le J p \le c_u - c $$
!
!  with the variable bounds as DAQP's simple bounds (which it handles without
!  forming their rows), and the linearized constraints as its general rows.
!  DAQP factors \( H \) once per QP (Cholesky), and then updates the
!  \( LDL^T \) factors of the working set's Gram matrix as rows enter and leave
!  it, where the dense solver refactors at every iteration. Being a dual
!  method, it needs no feasible starting point.
!
!  It only solves *convex* QPs with *consistent* constraints. When DAQP
!  doesn't solve the QP (an indefinite \( H \), inconsistent linearized
!  constraints, cycling, the iteration limit, ...), [[solve_daqp_qp]] says so,
!  and [[sqpopt_qp_solver_module]] solves it with the dense active-set solver
!  instead, which has the elastic mode and follows negative curvature; the
!  forced elastic re-solves go to the dense solver directly. A singular
!  \( H \) goes to the dense solver too: DAQP's proximal-point loop for
!  semidefinite Hessians is turned off (`eps_prox = 0`), because on the
!  Hock-Schittkowski problems it treated indefinite SR1 Hessians as
!  semidefinite and ran into the iteration limit (17 times the time, and
!  worse results), and did worse with the exact Hessian too.
!
!  **Warm start.** Each QP starts from the working set given (the previous
!  QP's final one, which this solver shares with the dense one, see
!  [[working_set]]), mapped to DAQP's numbering; DAQP drops the rows of it
!  that are linearly dependent on the others. Otherwise (after `setup`), it
!  starts from the equality constraints and fixed variables.
!
!  This module is the only code that uses `daqp_module`.

    module sqpopt_qp_daqp_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_qp_solve_failed, &
                                     sqpopt_infinity, sqpopt_out_of_memory
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_log_module,     only: fmt_i
    use daqp_module,           only: daqp_type, daqp_wp, daqp_ip, daqp_inf, daqp_optimal, daqp_soft_optimal, &
                                     daqp_no_freedom, daqp_optimal_inexact, daqp_infeasible, daqp_cycling, &
                                     daqp_unbounded, daqp_iteration_limit, daqp_nonconvex, daqp_overdetermined, &
                                     daqp_time_limit, daqp_unsupported, daqp_invalid_input, daqp_out_of_memory, &
                                     daqp_not_setup

    implicit none

    private

    type, public :: sqpopt_daqp_qp_type
        !! options of the DAQP QP solver, and the outcome of its last solve.

        integer  :: max_iter   = 1000      !! DAQP's limit on the active-set iterations of a QP solve (its
                                           !! `iter_limit`; at the limit, the dense QP solver takes over)
        real(wp) :: primal_tol = 1.0e-12_wp !! DAQP's tolerance on the violation of a constraint (its
                                           !! `primal_tol`, on the rows normalized by DAQP)
        real(wp) :: dual_tol   = 1.0e-12_wp !! DAQP's tolerance on the sign of a multiplier (its `dual_tol`)
        logical  :: warm_start = .true.    !! start from the previous QP's final working set, if the problem size
                                           !! is unchanged (see the module documentation)

        integer :: n_iter    = 0 !! number of DAQP's iterations in the last solve (output)
        integer :: n_working = 0 !! number of general rows and variable bounds in DAQP's final working set of the
                                 !! last solve (output)
        integer :: status    = 0 !! DAQP's status of the last solve (output; `daqp_optimal` if it solved the QP,
                                 !! see [[daqp_status_text]])

        contains

        procedure, public :: solve => solve_daqp_qp

    end type sqpopt_daqp_qp_type

    public :: daqp_status_text

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem with DAQP (see the module
!  documentation). `istat` is `sqpopt_success` if DAQP solved it, with `p`
!  the step, `lambda` the multipliers of the linearized constraints
!  (\( \lambda_i \ge 0 \) at a lower bound, as the other QP solvers return
!  them), and `warm_status` the final working set. Otherwise it is
!  `sqpopt_qp_solve_failed` (DAQP's status is in `me%status`), or
!  `sqpopt_out_of_memory` if a dense matrix (here, or in DAQP) can't be
!  allocated; `p` and
!  `lambda` are then zero, `warm_status` is unchanged, and the caller solves
!  the QP another way.

    subroutine solve_daqp_qp(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, warm_status, p, lambda, istat)

    class(sqpopt_daqp_qp_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! matrix-free Hessian approximation (densified here)
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! sparse constraint Jacobian, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x       !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g       !! objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c       !! constraint values `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: x_lb    !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: x_ub    !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c_lb    !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: c_ub    !! constraint upper bounds `dimension(m)`
    integer,  dimension(:), allocatable, intent(inout) :: warm_status !! the working set to start from (if
                                                         !! allocated with size `m+n`, and `me%warm_start`), and,
                                                         !! if solved, the final one: the side (`-1` lower, `+1`
                                                         !! upper, `0` none) of each general row, then of each
                                                         !! variable bound
    real(wp), dimension(:),     intent(out)   :: p       !! the step `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: lambda  !! multipliers of the linearized constraints `dimension(m)`
    integer,                    intent(out)   :: istat   !! status code (see above)

    type(daqp_type) :: qp
    real(wp),      dimension(:,:), allocatable :: h
    real(daqp_wp), dimension(:,:), allocatable :: hd, ad
    real(daqp_wp), dimension(:),   allocatable :: f, bu, bl, xd, lam
    integer(daqp_ip), dimension(:), allocatable :: active
    logical,          dimension(:), allocatable :: at_lower
    integer(daqp_ip) :: dstat
    integer :: n, m, i, k, nw, row, side, alloc_stat
    logical :: lower

    n = size(g)
    m = size(c)
    p      = 0.0_wp
    lambda = 0.0_wp
    me%n_iter    = 0
    me%n_working = 0
    me%status    = daqp_not_setup
    istat        = sqpopt_out_of_memory

    ! ---- the dense Hessian and Jacobian, in DAQP's real kind ----
    allocate(h(n,n), stat=alloc_stat)
    if (alloc_stat /= 0) return
    call hessian%dense(h)
    allocate(hd(n,n), stat=alloc_stat)
    if (alloc_stat /= 0) return
    hd = real(h, daqp_wp)
    deallocate(h)
    allocate(ad(m,n), stat=alloc_stat)
    if (alloc_stat /= 0) return
    ad = 0.0_daqp_wp
    do k = 1, jac%nnz
        ad(jac%irow(k), jac%icol(k)) = ad(jac%irow(k), jac%icol(k)) + real(jac%val(k), daqp_wp)
    end do

    ! ---- the bounds: the variables' (DAQP's simple bounds), then the rows' ----
    allocate(f(n), bu(n+m), bl(n+m), xd(n), lam(n+m))
    f = real(g, daqp_wp)
    do i = 1, n
        call set_bounds(i, x_lb(i), x_ub(i), x(i))
    end do
    do i = 1, m
        call set_bounds(n+i, c_lb(i), c_ub(i), c(i))
    end do

    qp%iter_limit = int(me%max_iter, daqp_ip)
    qp%primal_tol = real(me%primal_tol, daqp_wp)
    qp%dual_tol   = real(me%dual_tol, daqp_wp)
    qp%eps_prox   = 0.0_daqp_wp   ! (no proximal-point loop: see the module documentation)
    call qp%setup(hd, f, ad, bu, bl, dstat)
    deallocate(hd, ad)
    if (dstat < 0) then
        call not_solved(dstat)
        return
    end if

    ! ---- the warm start, in DAQP's numbering (bounds first) ----
    nw = 0
    if (me%warm_start .and. allocated(warm_status)) then
        if (size(warm_status) == m+n) nw = count(warm_status /= 0)
    end if
    allocate(active(nw), at_lower(nw))
    k = 0
    if (nw > 0) then
        do i = 1, m+n
            if (warm_status(i) == 0) cycle
            lower = warm_status(i) < 0
            ! (general row `i`, or the bound of variable `i-m`)
            row = merge(n+i, i-m, i <= m)
            ! (a side without a bound can't be active)
            if (lower .and. bl(row) <= -daqp_inf) cycle
            if (.not. lower .and. bu(row) >= daqp_inf) cycle
            k = k + 1
            active(k)   = int(row, daqp_ip)
            at_lower(k) = lower
        end do
    end if

    call qp%solve(xd, lam, dstat, active=active(1:k), at_lower=at_lower(1:k))
    if (dstat /= daqp_optimal) then
        call not_solved(dstat)
        return
    end if

    ! ---- the solution, multipliers, and working set ----
    p = real(xd, wp)
    do i = 1, m
        ! (DAQP's multipliers are negative at a lower bound)
        lambda(i) = -real(lam(n+i), wp)
    end do
    call qp%get_working_set(active, at_lower)
    if (allocated(warm_status)) deallocate(warm_status)
    allocate(warm_status(m+n), source=0)
    do k = 1, size(active)
        side = merge(-1, 1, at_lower(k))
        if (active(k) <= n) then
            warm_status(m+active(k)) = side
        else
            warm_status(active(k)-n) = side
        end if
    end do
    me%status    = dstat
    me%n_iter    = int(qp%iter)
    me%n_working = size(active)
    istat        = sqpopt_success

    contains

        subroutine set_bounds(row, lb, ub, v)
        !! the bounds of DAQP's row `row` on the step: `lb-v` and `ub-v`, or
        !! DAQP's infinity where `lb` or `ub` is infinite
        integer,  intent(in) :: row !! DAQP's row
        real(wp), intent(in) :: lb  !! lower bound on the variable or constraint
        real(wp), intent(in) :: ub  !! upper bound on the variable or constraint
        real(wp), intent(in) :: v   !! its value at `x`
        if (lb <= -sqpopt_infinity) then
            bl(row) = -daqp_inf
        else
            bl(row) = real(lb - v, daqp_wp)
        end if
        if (ub >= sqpopt_infinity) then
            bu(row) = daqp_inf
        else
            bu(row) = real(ub - v, daqp_wp)
        end if
        end subroutine set_bounds

        subroutine not_solved(daqp_istat)
        !! the outputs when DAQP didn't solve the QP
        integer(daqp_ip), intent(in) :: daqp_istat !! DAQP's status
        me%status = int(daqp_istat)
        me%n_iter = int(qp%iter)
        if (daqp_istat == daqp_out_of_memory) then
            istat = sqpopt_out_of_memory
        else
            istat = sqpopt_qp_solve_failed
        end if
        p      = 0.0_wp
        lambda = 0.0_wp
        end subroutine not_solved

    end subroutine solve_daqp_qp
!*******************************************************************************

!*******************************************************************************
!>
!  a short description of a DAQP status (`sqpopt_daqp_qp_type%status`), for
!  the printed output.

    pure function daqp_status_text(status) result(text)

    integer, intent(in) :: status !! DAQP's status
    character(len=:), allocatable :: text

    select case (status)
    case (daqp_optimal);         text = 'optimal'
    case (daqp_soft_optimal);    text = 'optimal with soft constraints violated'
    case (daqp_no_freedom);      text = 'no freedom'
    case (daqp_optimal_inexact); text = 'optimal after cycling, inexact'
    case (daqp_infeasible);      text = 'infeasible'
    case (daqp_cycling);         text = 'cycling'
    case (daqp_unbounded);       text = 'unbounded'
    case (daqp_iteration_limit); text = 'iteration limit'
    case (daqp_nonconvex);       text = 'nonconvex'
    case (daqp_overdetermined);  text = 'overdetermined initial working set'
    case (daqp_time_limit);      text = 'time limit'
    case (daqp_unsupported);     text = 'unsupported'
    case (daqp_invalid_input);   text = 'invalid input'
    case (daqp_out_of_memory);   text = 'out of memory'
    case (daqp_not_setup);       text = 'not set up'
    case default;                text = 'status '//fmt_i(status)
    end select

    end function daqp_status_text
!*******************************************************************************

    end module sqpopt_qp_daqp_module
!*******************************************************************************
