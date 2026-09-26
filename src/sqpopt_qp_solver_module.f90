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
!  1) after a shortened one. So the cap starts at `max_step`, but a
!  solution far away is still reached in a few iterations.

    module sqpopt_qp_solver_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix
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
        real(wp) :: max_step           = 2.0_wp                 !! initial trust-region-style cap on \( \lVert p \rVert_2 \):
                                                                 !! `p` is rescaled if it is longer than `max_step*step_scale`
        real(wp) :: step_scale         = 1.0_wp                 !! the adaptive factor on `max_step` (internal state,
                                                                 !! see the module docs; reset on each `solve`)
        logical  :: capped             = .false.                !! whether the last step was capped (output)
        integer  :: n_short            = 0                      !! consecutive very short line-search steps (internal
                                                                 !! state, see [[sqpopt_iterate_module]])
        integer :: n_iter = 0 !! number of active-set iterations taken by the last QP solve (output)
        type(sqpopt_dense_qp_type)           :: dense_qp    !! the dense QP solver (used only when `mode==sqpopt_qp_dense`)
        type(sqpopt_reduced_hessian_qp_type) :: sparse_qp   !! the sparse QP solver (used only when `mode==sqpopt_qp_reduced_hessian`)

        contains

        procedure, public :: solve => solve_qp_subproblem

    end type sqpopt_qp_solver_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  solve the linearized QP subproblem for the search direction `p` and
!  the associated Lagrange multipliers `lambda`, dispatching to the
!  active-set solver selected by `me%mode` (see [[sqpopt_qp_dense_module]],
!  [[sqpopt_qp_reduced_hessian_module]]).

    subroutine solve_qp_subproblem(me, hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

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

    select case (resolved_mode(me, size(g)))
    case (sqpopt_qp_dense)
        call me%dense_qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)
        me%n_iter = me%dense_qp%n_iter
    case default ! sqpopt_qp_reduced_hessian
        call me%sparse_qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)
        me%n_iter = me%sparse_qp%n_iter
    end select

    ! trust-region-style safeguard on the step length:
    me%capped = norm2(p) > me%max_step*me%step_scale
    if (me%capped) p = p*(me%max_step*me%step_scale/norm2(p))

    ! the active-set solvers only satisfy the bounds to within their own
    ! tolerances, so make sure `x+p` (and hence every `x+alpha*p`,
    ! `0<=alpha<=1`) is exactly within them -- the user functions must
    ! never be evaluated outside the variable bounds:
    p = min(max(x + p, x_lb), x_ub) - x

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

    end module sqpopt_qp_solver_module
!*******************************************************************************
