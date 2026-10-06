!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Solves the quadratic programming (QP) subproblem generated at each
!  major SQP iteration:
!
!  $$ \min_{p} \; \frac{1}{2} p^T H p + g^T p $$
!
!  subject to the linearized constraints and variable bounds:
!
!  $$ c_l - c(x) \le J(x) p \le c_u - c(x) $$
!  $$ x_l - x \le p \le x_u - x $$
!
!  The subproblem is solved by one of the active-set QP solvers, selected
!  by `mode`:
!
!  * `sqpopt_qp_dense`: a dense active-set QP (see [[sqpopt_qp_dense_module]]);
!  * `sqpopt_qp_reduced_hessian`: a sparse/matrix-free active-set QP (see
!    [[sqpopt_qp_reduced_hessian_module]]);
!  * `sqpopt_qp_daqp`: DAQP's dual active-set solver for dense convex QPs
!    (see [[sqpopt_qp_daqp_module]]), with `sqpopt_qp_dense` as its
!    fallback: a QP that DAQP doesn't solve (a nonconvex one, inconsistent
!    linearized constraints, ...), and a forced elastic re-solve, are solved
!    by the dense solver. The two share the working set that each QP starts
!    from;
!  * `sqpopt_qp_auto` (the default): `sqpopt_qp_dense` when
!    `n <= auto_dense_max_n`, else `sqpopt_qp_reduced_hessian`.
!
!  Before either is run, the *unconstrained step* is tried (unless
!  `unconstrained_step` is false): with the limited-memory BFGS Hessian,
!  which is positive definite, the minimizer of the QP's objective alone is
!  \( p = -H^{-1} g \), which the two-loop recursion gives in `O(nk)`
!  operations for `k` stored pairs (see [[hessian_inverse_vector_product]]).
!  If that step satisfies the bounds and the linearized constraints, it is
!  the solution of the QP, with zero multipliers, and no QP solver is run.
!  That is every QP of a problem without constraints whose bounds aren't
!  active, where the active-set solvers would otherwise work with every
!  variable free (the sparse one by conjugate gradients on a reduced
!  Hessian of order `n`). When the step isn't feasible, the attempt costs
!  one such product.
!
!  With `direct` (`options%direct_qp`), the
!  subproblem is then tried directly, by factoring the KKT matrix of the
!  working set that the active-set solver would start from (see
!  [[sqpopt_qp_direct_module]]); the active-set solver is only run if that
!  doesn't give the solution within `direct_max_changes` changes of the
!  working set.
!
!  Both enforce the variable bounds and the linearized constraints exactly
!  as part of the QP solve. The step length is then capped at
!  `max_step*step_scale` (a trust-region-style safeguard), where the major
!  iterations adapt `step_scale` like a trust radius: it doubles after a
!  capped step that the line search accepted in full, and halves (down to
!  1) after a shortened one. Each solve starts with
!  `step_scale = max(1, ||x0||_inf/max_step)`, so the cap starts at
!  \( \max(\) `max_step` \(, \lVert x_0 \rVert_\infty) \): relative to the
!  size of the variables, and a solution far away is still reached in a
!  few iterations.

    module sqpopt_qp_solver_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use, intrinsic :: iso_fortran_env, only: int64
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_out_of_memory, sqpopt_success, sqpopt_all_finite
    use sqpopt_linalg_module,  only: sparse_matvec
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_qp_dense_module, only: sqpopt_dense_qp_type
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_reduced_hessian_qp_type
    use sqpopt_qp_daqp_module,   only: sqpopt_daqp_qp_type
    use sqpopt_kkt_module,       only: sqpopt_kkt_type
    use sqpopt_inertia_module,   only: sqpopt_inertia_type
    use sqpopt_qp_direct_module, only: direct_qp_step, sqpopt_direct_solved

    implicit none

    private

    integer, parameter, public :: sqpopt_qp_auto            = 0  !! (default) `sqpopt_qp_dense` if `n <= auto_dense_max_n`,
                                                                  !! else `sqpopt_qp_reduced_hessian`
    integer, parameter, public :: sqpopt_qp_dense           = 2  !! dense active-set QP solver (see [[sqpopt_qp_dense_module]])
    integer, parameter, public :: sqpopt_qp_reduced_hessian = 3  !! sparse (projected-CG) active-set QP solver (see [[sqpopt_qp_reduced_hessian_module]])
                                                                  !! (value 1 was the removed composite-step heuristic)
    integer, parameter, public :: sqpopt_qp_daqp            = 4  !! DAQP's dual active-set solver for dense convex QPs, with
                                                                  !! `sqpopt_qp_dense` as its fallback (see [[sqpopt_qp_daqp_module]])

    type, public :: sqpopt_qp_solver_type
        !! workspace and options for the QP subproblem solver.

        integer  :: mode                = sqpopt_qp_auto       !! which QP algorithm to use (see the `sqpopt_qp_*` constants)
        integer  :: auto_dense_max_n    = 200                  !! `mode==sqpopt_qp_auto` uses the dense QP solver for problems
                                                                 !! with at most this many variables, and the sparse
                                                                 !! reduced-Hessian QP solver for larger ones
        real(wp) :: max_step           = 2.0_wp                 !! trust-region-style cap on \( \lVert p \rVert_2 \):
                                                                 !! `p` is rescaled if it is longer than `max_step*step_scale`
                                                                 !! (the cap starts at \( \max(\) `max_step`
                                                                 !! \(, \lVert x_0 \rVert_\infty) \), see the module docs)
        real(wp) :: step_scale         = 1.0_wp                 !! the adaptive factor on `max_step` (internal state,
                                                                 !! see the module docs; reset on each `solve`)
        logical  :: capped             = .false.                !! whether the last step was capped (output)
        integer  :: n_elastic          = 0                      !! elastic re-solves in this solve (internal state, see
                                                                 !! [[sqpopt_iterate_module]])
        integer  :: n_short            = 0                      !! consecutive very short line-search steps (internal
                                                                 !! state, see [[sqpopt_iterate_module]])
        integer :: n_iter = 0 !! number of active-set iterations taken by the last QP solve (output)
        integer :: n_working = 0 !! size of the final working set of the last QP solve (output)
        integer :: n_slacks  = 0 !! number of elastic slacks in the last QP solve (output)
        real(wp) :: time = 0.0_wp !! wall-clock time spent in QP solves (output, seconds; reset on each `solve`)
        logical :: negative_curvature = .false. !! whether the last QP solve found negative curvature of the Hessian
                                                !! in the variables (output; see [[sqpopt_iterate_module]])
        logical :: out_of_memory = .false. !! whether a QP solve of this `solve` returned `sqpopt_out_of_memory`
                                           !! (output; it stays set, so the solver stops whichever step asked
                                           !! for that QP; reset on each `solve`)
        logical :: unconstrained_step = .true. !! whether to try the unconstrained step first, with the
                                           !! limited-memory BFGS Hessian: if \( -H^{-1} g \) satisfies the bounds
                                           !! and the linearized constraints, it is the QP's solution, and no QP
                                           !! solver is run (see the module docs)
        logical :: unconstrained_used = .false. !! whether the last QP was solved by the unconstrained step (output)
        integer :: n_unconstrained = 0     !! number of QPs of this `solve` solved by the unconstrained step (output)
        logical :: direct = .false.        !! whether to try the direct method first (overwritten from
                                           !! `options%direct_qp`; see [[sqpopt_qp_direct_module]])
        integer :: direct_max_changes = 10 !! the direct method gives up, and the active-set solver is run,
                                           !! after this many changes of the working set (each one is a
                                           !! factorization) without reaching the solution
        real(wp) :: direct_tol = 1.0e-8_wp !! the direct method's relative tolerance for a violated row or
                                           !! bound, and for the sign of a multiplier
        logical :: direct_used = .false.   !! whether the last QP was solved by the direct method (output)
        integer :: direct_outcome = -1     !! how the direct method ended in the last QP solve: a `sqpopt_direct_*`
                                           !! constant (see [[sqpopt_qp_direct_module]]), or `-1` if it wasn't
                                           !! tried (output)
        integer :: direct_changes = 0      !! the changes of the working set it made there (output)
        integer :: n_solves = 0            !! number of QP solves in this `solve` (output)
        integer :: n_direct = 0            !! of which, by the direct method (output)
        logical :: daqp_fallback = .false. !! whether DAQP didn't solve the last QP, and the dense QP solver did
                                           !! (output; with `mode==sqpopt_qp_daqp`)
        logical :: daqp_used = .false.     !! whether DAQP solved the last QP (output; with `mode==sqpopt_qp_daqp`,
                                           !! the dense QP solver solved it otherwise: a fallback, or a forced
                                           !! elastic re-solve)
        integer :: n_daqp_fallbacks = 0    !! number of QPs of this `solve` that DAQP didn't solve, solved by the
                                           !! dense QP solver instead (output)
        type(sqpopt_dense_qp_type)           :: dense_qp    !! the dense QP solver (used only when `mode` is
                                                            !! `sqpopt_qp_dense`, or `sqpopt_qp_daqp` as its fallback)
        type(sqpopt_reduced_hessian_qp_type) :: sparse_qp   !! the sparse QP solver (used only when `mode==sqpopt_qp_reduced_hessian`)
        type(sqpopt_daqp_qp_type)            :: daqp_qp     !! the DAQP QP solver (used only when `mode==sqpopt_qp_daqp`)

        contains

        procedure, public :: solve => solve_qp_subproblem
        procedure, public :: mode_name
        procedure, public :: solver_name
        procedure, public :: working_set
        procedure, public :: starting_working_set
        procedure, public :: elastic_slacks

    end type sqpopt_qp_solver_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda`, dispatching to the
!  active-set solver selected by `me%mode` (see [[sqpopt_qp_dense_module]],
!  [[sqpopt_qp_reduced_hessian_module]], [[sqpopt_qp_daqp_module]]; with
!  `sqpopt_qp_daqp`, a QP that DAQP doesn't solve, and a forced elastic
!  solve, go to the dense solver, and `me%daqp_fallback` is set in the
!  first case). If the solver returns
!  `istat=sqpopt_out_of_memory`, `me%out_of_memory` is also set (and stays
!  set), for [[sqpopt_iterate]] to stop the solve.
!
!  With the limited-memory BFGS Hessian (and `me%unconstrained_step`, and
!  not in a forced elastic solve), the unconstrained step is tried first
!  (see the module documentation): if it is feasible, it is the QP's
!  solution (`me%unconstrained_used`), no solver is run, and the next QP
!  starts from an empty working set.
!
!  If `kkt` is given (and `me%direct` is set, and this is not a forced
!  elastic solve), the direct method is tried
!  next (see [[direct_qp_step]]), from the working set that the active-set
!  solver would start from. If it finds the QP's solution, the active-set
!  solver is not run (`me%direct_used`), and the working set it ended with
!  is kept as that solver's warm start. With `inertia`, the direct method
!  may increase the exact Hessian's shift (`hessian%shift`).
!
!  The time spent here is added to `me%time`, without the time that `kkt`
!  spent factoring and solving (which its solver counts itself).

    subroutine solve_qp_subproblem(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat, &
                                   elastic_sign, elastic_weight, kkt, inertia)

    class(sqpopt_qp_solver_type), intent(inout) :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! matrix-free Hessian approximation (never a dense `n x n` matrix)
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! sparse constraint Jacobian, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x       !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: g       !! objective gradient `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c       !! current constraint values `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: x_lb    !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: x_ub    !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c_lb    !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)    :: c_ub    !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:),     intent(out)   :: p       !! computed search direction `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: lambda  !! Lagrange multipliers for the linearized constraints `dimension(m)`
    integer,                    intent(out)   :: istat   !! status code (see [[sqpopt_types_module]])
    integer,  dimension(:), optional, intent(in) :: elastic_sign   !! forced elastic mode: the rows to make
                                                                  !! elastic, and in which direction (`+1`
                                                                  !! relaxes a row's lower bound, `-1` its upper
                                                                  !! bound, `0` none) `dimension(m)`
    real(wp),               optional, intent(in) :: elastic_weight !! forced elastic mode: the fixed \( \ell_1 \)
                                                                  !! weight of every elastic slack
    type(sqpopt_kkt_type),  optional, intent(inout) :: kkt         !! the KKT matrix, for the direct method (see
                                                                  !! [[sqpopt_kkt_module]])
    type(sqpopt_inertia_type), optional, intent(inout) :: inertia  !! the inertia control, for the direct method
                                                                  !! at a nonconvex face (see
                                                                  !! [[sqpopt_inertia_module]])

    logical :: forced
    integer(int64) :: t0, t1, rate
    real(wp) :: t_kkt !! the time `kkt` spent in its solver during this call

    call system_clock(t0, rate)
    forced = present(elastic_sign) .and. present(elastic_weight)
    me%n_solves    = me%n_solves + 1
    me%unconstrained_used = .false.
    me%daqp_fallback  = .false.
    me%daqp_used      = .false.
    me%direct_used    = .false.
    me%direct_outcome = -1
    me%direct_changes = 0
    t_kkt = 0.0_wp

    if (me%unconstrained_step .and. .not. forced) call solve_unconstrained()
    if (.not. me%unconstrained_used) then
        if (present(kkt) .and. me%direct .and. .not. forced) call solve_direct()
        if (.not. me%direct_used) call solve_active_set()
    end if

    if (istat == sqpopt_out_of_memory) me%out_of_memory = .true.

    ! trust-region-style safeguard on the step length:
    me%capped = norm2(p) > me%max_step*me%step_scale
    if (me%capped) p = p*(me%max_step*me%step_scale/norm2(p))

    ! the solvers only satisfy the bounds to within their own
    ! tolerances, so make sure `x+p` (and hence every `x+alpha*p`,
    ! `0<=alpha<=1`) is exactly within them -- the user functions must
    ! never be evaluated outside the variable bounds:
    p = min(max(x + p, x_lb), x_ub) - x

    call system_clock(t1)
    me%time = me%time + max(real(t1 - t0, wp)/real(rate, wp) - t_kkt, 0.0_wp)

    contains

        subroutine solve_unconstrained()
        !! try the unconstrained step \( -H^{-1} g \) of the limited-memory
        !! BFGS Hessian (see the module documentation). Sets
        !! `me%unconstrained_used` if it is feasible, and so solves the QP.
        real(wp), dimension(:), allocatable :: jp, neg_g
        integer :: n, m
        ! (only for a positive definite matrix with a cheap inverse)
        if (hessian%exact .or. hessian%use_sr1 .or. hessian%shift /= 0.0_wp) return
        n = size(g)
        m = size(c)
        allocate(jp(m), neg_g(n))
        neg_g = -g
        call hessian%inverse_vector_product(neg_g, p)
        if (.not. sqpopt_all_finite(p)) return
        if (.not. dot_product(g, p) <= 0.0_wp) return
        if (any(x + p < x_lb) .or. any(x + p > x_ub)) return
        if (m > 0) then
            call sparse_matvec(jac, p, jp)
            if (any(c + jp < c_lb) .or. any(c + jp > c_ub)) return
        end if
        me%unconstrained_used = .true.
        me%n_unconstrained = me%n_unconstrained + 1
        lambda = 0.0_wp
        istat  = sqpopt_success
        me%n_iter    = 0
        me%n_working = 0
        me%n_slacks  = 0
        me%negative_curvature = .false.
        ! (the next QP starts from the empty working set)
        if (uses_dense_working_set(me, n)) then
            if (allocated(me%dense_qp%warm_status)) deallocate(me%dense_qp%warm_status)
            allocate(me%dense_qp%warm_status(m + n), source=0)
        else
            if (allocated(me%sparse_qp%warm_status)) deallocate(me%sparse_qp%warm_status)
            allocate(me%sparse_qp%warm_status(m + n), source=0)
        end if
        end subroutine solve_unconstrained

        subroutine solve_direct()
        !! try the direct method (see [[sqpopt_qp_direct_module]]), from the
        !! working set that the active-set solver would start from. Sets
        !! `me%direct_used` if it solved the QP.
        integer, dimension(:), allocatable :: status
        if (.not. kkt%enabled) return
        t_kkt = kkt%solver%time
        call me%starting_working_set(x_lb, x_ub, c_lb, c_ub, status)
        call direct_qp_step(kkt, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, me%direct_max_changes, &
                            me%direct_tol, status, p, lambda, me%direct_changes, me%direct_outcome, inertia=inertia)
        t_kkt = kkt%solver%time - t_kkt
        me%direct_used = me%direct_outcome == sqpopt_direct_solved
        if (.not. me%direct_used) return
        istat = sqpopt_success
        me%n_direct  = me%n_direct + 1
        me%n_iter    = me%direct_changes
        me%n_working = count(status /= 0)
        me%n_slacks  = 0
        me%negative_curvature = .false.
        ! (the next QP starts from this working set)
        if (uses_dense_working_set(me, size(g))) then
            me%dense_qp%warm_status = status
        else
            me%sparse_qp%warm_status = status
        end if
        end subroutine solve_direct

        subroutine solve_active_set()
        !! solve the QP with the active-set solver selected by `me%mode`
        select case (resolved_mode(me, size(g)))
        case (sqpopt_qp_daqp)
            if (.not. forced) then
                call me%daqp_qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, me%dense_qp%warm_status, &
                                      p, lambda, istat)
                if (istat == sqpopt_success) then
                    me%daqp_used = .true.
                    me%n_iter    = me%daqp_qp%n_iter
                    me%n_working = me%daqp_qp%n_working
                    me%n_slacks  = 0
                    me%negative_curvature = .false.
                    return
                end if
                me%daqp_fallback    = .true.
                me%n_daqp_fallbacks = me%n_daqp_fallbacks + 1
            end if
            call solve_dense()
        case (sqpopt_qp_dense)
            call solve_dense()
        case default ! sqpopt_qp_reduced_hessian
            if (forced) then
                me%sparse_qp%force_sign   = elastic_sign
                me%sparse_qp%force_weight = elastic_weight
            end if
            call me%sparse_qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)
            if (forced) then
                deallocate(me%sparse_qp%force_sign)
                me%sparse_qp%force_weight = 0.0_wp
            end if
            me%n_iter = me%sparse_qp%n_iter
            me%n_working = me%sparse_qp%n_working
            me%n_slacks  = me%sparse_qp%n_slacks
            me%negative_curvature = me%sparse_qp%negative_curvature
        end select
        end subroutine solve_active_set

        subroutine solve_dense()
        !! solve the QP with the dense active-set solver
        if (forced) then
            me%dense_qp%force_sign   = elastic_sign
            me%dense_qp%force_weight = elastic_weight
        end if
        call me%dense_qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)
        if (forced) then
            deallocate(me%dense_qp%force_sign)
            me%dense_qp%force_weight = 0.0_wp
        end if
        me%n_iter = me%dense_qp%n_iter
        me%n_working = me%dense_qp%n_working
        me%n_slacks  = me%dense_qp%n_slacks
        me%negative_curvature = me%dense_qp%negative_curvature
        end subroutine solve_dense

    end subroutine solve_qp_subproblem
!*******************************************************************************

!*******************************************************************************
!>
!  the QP algorithm actually used for a problem with `n` variables:
!  `me%mode`, with `sqpopt_qp_auto` resolved to a specific solver.

    pure integer function resolved_mode(me, n)

    class(sqpopt_qp_solver_type), intent(in) :: me
    integer,                      intent(in) :: n !! number of variables

    if (me%mode == sqpopt_qp_auto) then
        resolved_mode = merge(sqpopt_qp_dense, sqpopt_qp_reduced_hessian, n <= me%auto_dense_max_n)
    else
        resolved_mode = me%mode
    end if

    end function resolved_mode
!*******************************************************************************

!*******************************************************************************
!>
!  whether the working set of the QPs with `n` variables is kept by the
!  dense QP solver: with `sqpopt_qp_dense`, and with `sqpopt_qp_daqp`,
!  which shares it with the dense solver, its fallback.

    pure logical function uses_dense_working_set(me, n)

    class(sqpopt_qp_solver_type), intent(in) :: me
    integer,                      intent(in) :: n !! number of variables

    uses_dense_working_set = any(resolved_mode(me, n) == [sqpopt_qp_dense, sqpopt_qp_daqp])

    end function uses_dense_working_set
!*******************************************************************************

!*******************************************************************************
!>
!  the final working set of the last QP solve with `n` variables and `m`
!  constraints: the side (`-1` lower, `+1` upper, `0` not in the working
!  set) of each general row, then of each variable bound. `status` is
!  returned unallocated if there is none (no QP of that size has been
!  solved yet).

    subroutine working_set(me, n, m, status)

    class(sqpopt_qp_solver_type),       intent(in)  :: me
    integer,                            intent(in)  :: n      !! number of variables
    integer,                            intent(in)  :: m      !! number of constraints
    integer, dimension(:), allocatable, intent(out) :: status !! the working set `dimension(m+n)` (see above)

    if (uses_dense_working_set(me, n)) then
        if (allocated(me%dense_qp%warm_status)) status = me%dense_qp%warm_status
    else
        if (allocated(me%sparse_qp%warm_status)) status = me%sparse_qp%warm_status
    end if
    if (allocated(status)) then
        if (size(status) /= m + n) deallocate(status)
    end if

    end subroutine working_set
!*******************************************************************************

!*******************************************************************************
!>
!  the elastic slacks of the last QP solve with `n` variables: the general
!  row of each, and its value at the QP's solution (by how much the step
!  violates that linearized constraint). They are only kept if the
!  active-set solvers' `keep_slacks` is set (the solver sets it for its
!  diagnostics); `row` and `slack` have size zero if there are none.

    subroutine elastic_slacks(me, n, row, slack)

    class(sqpopt_qp_solver_type),        intent(in)  :: me
    integer,                             intent(in)  :: n     !! number of variables
    integer,  dimension(:), allocatable, intent(out) :: row   !! the row of each elastic slack
    real(wp), dimension(:), allocatable, intent(out) :: slack !! its value

    if (me%n_slacks > 0) then
        ! (with `sqpopt_qp_daqp`, slacks only come from its fallback, the dense solver)
        if (uses_dense_working_set(me, n)) then
            if (allocated(me%dense_qp%slack_row)) then
                allocate(row(size(me%dense_qp%slack_row)), slack(size(me%dense_qp%slack_row)))
                row   = me%dense_qp%slack_row
                slack = me%dense_qp%slack_value
            end if
        else
            if (allocated(me%sparse_qp%slack_row)) then
                allocate(row(size(me%sparse_qp%slack_row)), slack(size(me%sparse_qp%slack_row)))
                row   = me%sparse_qp%slack_row
                slack = me%sparse_qp%slack_value
            end if
        end if
    end if
    if (.not. allocated(row)) allocate(row(0), slack(0))

    end subroutine elastic_slacks
!*******************************************************************************

!*******************************************************************************
!>
!  the working set that the next QP solve will start from, for a problem
!  with the given bounds: the final working set of the last one (see
!  [[working_set]]), or, if there is none, the equality constraints and
!  the fixed variables (the active-set solvers' crash start).

    subroutine starting_working_set(me, x_lb, x_ub, c_lb, c_ub, status)

    class(sqpopt_qp_solver_type),       intent(in)  :: me
    real(wp), dimension(:),             intent(in)  :: x_lb   !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),             intent(in)  :: x_ub   !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),             intent(in)  :: c_lb   !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),             intent(in)  :: c_ub   !! constraint upper bounds `dimension(m)`
    integer, dimension(:), allocatable, intent(out) :: status !! the working set `dimension(m+n)` (see [[working_set]])

    integer :: n, m

    n = size(x_lb)
    m = size(c_lb)
    call me%working_set(n, m, status)
    if (allocated(status)) return
    allocate(status(m+n))
    status = 0
    where (c_ub - c_lb <= 0.0_wp) status(1:m) = -1
    where (x_ub - x_lb <= 0.0_wp) status(m+1:m+n) = -1

    end subroutine starting_working_set
!*******************************************************************************

!*******************************************************************************
!>
!  the name of the QP algorithm used for a problem with `n` variables
!  (`mode` with `sqpopt_qp_auto` resolved), for the printed output.

    function mode_name(me, n) result(name)

    class(sqpopt_qp_solver_type), intent(in) :: me
    integer,                      intent(in) :: n !! number of variables (which `sqpopt_qp_auto` depends on)
    character(len=:), allocatable :: name

    select case (resolved_mode(me, n))
    case (sqpopt_qp_dense); name = 'dense QP'
    case (sqpopt_qp_daqp);  name = 'DAQP QP'
    case default;           name = 'sparse QP'
    end select
    if (me%mode == sqpopt_qp_auto) name = name//' (auto)'

    end function mode_name
!*******************************************************************************

!*******************************************************************************
!>
!  the name of the QP solver that solved the last QP with `n` variables, for
!  the printed output: [[mode_name]], except that with `sqpopt_qp_daqp`, a
!  QP that DAQP didn't solve (a fallback, or a forced elastic re-solve) was
!  solved by the dense QP solver.

    function solver_name(me, n) result(name)

    class(sqpopt_qp_solver_type), intent(in) :: me
    integer,                      intent(in) :: n !! number of variables
    character(len=:), allocatable :: name

    if (resolved_mode(me, n) == sqpopt_qp_daqp .and. .not. me%daqp_used) then
        name = 'dense QP'
    else
        name = me%mode_name(n)
    end if

    end function solver_name
!*******************************************************************************

    end module sqpopt_qp_solver_module
!*******************************************************************************
