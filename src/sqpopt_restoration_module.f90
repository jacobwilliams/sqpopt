!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  A simple feasibility-restoration step, used by [[sqpopt_iterate_module]]
!  when the QP subproblem reports that the linearized constraints are
!  inconsistent (`sqpopt_infeasible`), so that no QP step can be trusted.
!  Instead of a QP step, a Gauss-Newton step on the constraint violation
!  is taken:
!
!  $$ \min_x \; \tfrac12 \lVert r_c(x) \rVert_2^2 \quad \text{s.t.} \quad x_l \le x \le x_u $$
!
!  where \( r_c \) is the signed violation of each constraint's bounds.
!  Repeated restoration steps converge to either a feasible point (after
!  which the normal SQP iterations resume) or a point that is stationary
!  for the violation, which [[check_convergence]] then reports as
!  `sqpopt_infeasible`.
!
!  @note This is a lightweight stand-in for a full elastic-mode QP /
!  feasibility-restoration phase (see `plan/ROADMAP.md`, F2).

    module sqpopt_restoration_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_line_search_failed
    use sqpopt_problem_module, only: sqpopt_problem_type
    use lsqr_module,           only: lsqr_solver_ez

    implicit none

    private

    public :: restoration_step

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  take one Gauss-Newton step toward feasibility of the constraints (see
!  the module-level documentation), with a backtracking Armijo search on
!  \( \tfrac12 \lVert r_c \rVert_2^2 \). If no decrease is found, `x_new=x`
!  and `istat=sqpopt_line_search_failed`.

    subroutine restoration_step(problem, jac, x, c, max_step, x_new, alpha, istat)

    type(sqpopt_problem_type),  intent(inout) :: problem  !! problem definition
    type(sqpopt_sparse_matrix), intent(in)    :: jac      !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)    :: x        !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)    :: c        !! constraint values at `x` `dimension(m)`
    real(wp),                   intent(in)    :: max_step !! cap on \( \lVert p \rVert_2 \)
    real(wp), dimension(:),     intent(out)   :: x_new    !! new point `dimension(n)`
    real(wp),                   intent(out)   :: alpha    !! accepted step length (`0` if none)
    integer,                    intent(out)   :: istat    !! status code (see [[sqpopt_types_module]])

    real(wp), parameter :: sigma     = 1.0e-4_wp !! Armijo sufficient-decrease parameter
    real(wp), parameter :: backtrack = 0.5_wp    !! step-length reduction factor
    integer,  parameter :: max_ls    = 30        !! maximum number of backtracking steps

    real(wp), dimension(size(c)) :: rc, c_trial
    real(wp), dimension(size(x)) :: p, x_trial
    real(wp) :: h0, h_trial, dh0
    type(lsqr_solver_ez) :: lsqr
    integer :: istop, it

    rc = violation(c, problem%c_lb, problem%c_ub)
    h0 = 0.5_wp*dot_product(rc, rc)

    ! Gauss-Newton step: minimum-norm solution of J*p = -r_c, then made to
    ! respect the variable bounds and the step-length cap:
    call lsqr%initialize(problem%m, problem%n, jac%val, jac%irow, jac%icol)
    call lsqr%solve(-rc, 0.0_wp, p, istop)
    p = min(max(x+p, problem%x_lb), problem%x_ub) - x
    if (norm2(p) > max_step) p = p*(max_step/norm2(p))

    ! directional derivative of 0.5*|r_c|^2 along p: r_c^T J p
    block
        real(wp), dimension(size(c)) :: jp
        integer :: k
        jp = 0.0_wp
        do k = 1, jac%nnz
            jp(jac%irow(k)) = jp(jac%irow(k)) + jac%val(k)*p(jac%icol(k))
        end do
        dh0 = dot_product(rc, jp)
    end block

    alpha = 1.0_wp
    if (dh0 < 0.0_wp) then
        do it = 1, max_ls
            x_trial = x + alpha*p
            call problem%eval_c(x_trial, c_trial)
            rc = violation(c_trial, problem%c_lb, problem%c_ub)
            h_trial = 0.5_wp*dot_product(rc, rc)
            if (h_trial <= h0 + sigma*alpha*dh0) then
                x_new = x_trial
                istat = sqpopt_success
                return
            end if
            alpha = backtrack*alpha
        end do
    end if

    ! no decrease in the violation along the Gauss-Newton direction:
    alpha = 0.0_wp
    x_new = x
    istat = sqpopt_line_search_failed

    end subroutine restoration_step
!*******************************************************************************

!*******************************************************************************
!>
!  the signed violation of the constraint bounds: `c-c_lb` below the lower
!  bound, `c-c_ub` above the upper bound, and `0` in between.

    pure function violation(c, c_lb, c_ub) result(rc)

    real(wp), dimension(:), intent(in) :: c    !! constraint values `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_lb !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_ub !! constraint upper bounds `dimension(m)`
    real(wp), dimension(size(c)) :: rc

    rc = min(c-c_lb, 0.0_wp) + max(c-c_ub, 0.0_wp)

    end function violation
!*******************************************************************************

    end module sqpopt_restoration_module
!*******************************************************************************
