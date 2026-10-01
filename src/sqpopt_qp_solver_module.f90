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
!  The subproblem is solved by one of two active-set QP solvers, selected
!  by `mode`:
!
!  * `sqpopt_qp_dense`: a dense active-set QP (see [[sqpopt_qp_dense_module]]);
!  * `sqpopt_qp_reduced_hessian`: a sparse/matrix-free active-set QP (see
!    [[sqpopt_qp_reduced_hessian_module]]);
!  * `sqpopt_qp_auto` (the default): `sqpopt_qp_dense` when
!    `n <= auto_dense_max_n`, else `sqpopt_qp_reduced_hessian`.
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
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_out_of_memory
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_qp_dense_module, only: sqpopt_dense_qp_type
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_reduced_hessian_qp_type

    implicit none

    private

    integer, parameter, public :: sqpopt_qp_auto            = 0  !! (default) `sqpopt_qp_dense` if `n <= auto_dense_max_n`,
                                                                  !! else `sqpopt_qp_reduced_hessian`
    integer, parameter, public :: sqpopt_qp_dense           = 2  !! dense active-set QP solver (see [[sqpopt_qp_dense_module]])
    integer, parameter, public :: sqpopt_qp_reduced_hessian = 3  !! sparse (projected-CG) active-set QP solver (see [[sqpopt_qp_reduced_hessian_module]])
                                                                  !! (value 1 was the removed composite-step heuristic)

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
        type(sqpopt_dense_qp_type)           :: dense_qp    !! the dense QP solver (used only when `mode==sqpopt_qp_dense`)
        type(sqpopt_reduced_hessian_qp_type) :: sparse_qp   !! the sparse QP solver (used only when `mode==sqpopt_qp_reduced_hessian`)

        contains

        procedure, public :: solve => solve_qp_subproblem
        procedure, public :: mode_name
        procedure, public :: working_set

    end type sqpopt_qp_solver_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda`, dispatching to the
!  active-set solver selected by `me%mode` (see [[sqpopt_qp_dense_module]],
!  [[sqpopt_qp_reduced_hessian_module]]). If the solver returns
!  `istat=sqpopt_out_of_memory`, `me%out_of_memory` is also set (and stays
!  set), for [[sqpopt_iterate]] to stop the solve.

    subroutine solve_qp_subproblem(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat, &
                                   elastic_sign, elastic_weight)

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

    logical :: forced
    integer(int64) :: t0, t1, rate

    call system_clock(t0, rate)
    forced = present(elastic_sign) .and. present(elastic_weight)

    select case (resolved_mode(me, size(g)))
    case (sqpopt_qp_dense)
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

    if (istat == sqpopt_out_of_memory) me%out_of_memory = .true.

    ! trust-region-style safeguard on the step length:
    me%capped = norm2(p) > me%max_step*me%step_scale
    if (me%capped) p = p*(me%max_step*me%step_scale/norm2(p))

    ! the active-set solvers only satisfy the bounds to within their own
    ! tolerances, so make sure `x+p` (and hence every `x+alpha*p`,
    ! `0<=alpha<=1`) is exactly within them -- the user functions must
    ! never be evaluated outside the variable bounds:
    p = min(max(x + p, x_lb), x_ub) - x

    call system_clock(t1)
    me%time = me%time + real(t1 - t0, wp)/real(rate, wp)

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

    select case (resolved_mode(me, n))
    case (sqpopt_qp_dense)
        if (allocated(me%dense_qp%warm_status)) status = me%dense_qp%warm_status
    case default
        if (allocated(me%sparse_qp%warm_status)) status = me%sparse_qp%warm_status
    end select
    if (allocated(status)) then
        if (size(status) /= m + n) deallocate(status)
    end if

    end subroutine working_set
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
    case default;           name = 'sparse QP'
    end select
    if (me%mode == sqpopt_qp_auto) name = name//' (auto)'

    end function mode_name
!*******************************************************************************

    end module sqpopt_qp_solver_module
!*******************************************************************************
