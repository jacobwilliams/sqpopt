!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Second-order correction (SOC) step, used by the line searches in
!  [[sqpopt_linesearch_module]] and by the trust-region step in
!  [[sqpopt_trust_region_module]] when a full step has been rejected.
!
!  For strongly nonlinear constraints, a genuinely good step `p` can be
!  rejected by the merit function or filter because the *linearized*
!  constraint prediction \( c + Jp \) differs from the true value
!  \( c(x+p) \) (the Maratos effect). The SOC step adds a small
!  correction `d` that accounts for the true constraint values at `x+p`,
!  reusing the Jacobian at `x` (no new derivative evaluations):
!
!  $$ J_S\, d = -\left( c_S(x+p) - b_S \right) $$
!
!  solved in the minimum-norm sense with `LSQR` (or directly, with
!  `options%direct_least_squares`: see [[sqpopt_least_squares_module]]).
!  Here `S` is the set of
!  constraints that are active in the linearization (\( c+Jp \) at a bound
!  `b`), or violated at `x+p` (with `b` the violated bound). The corrected
!  step `p+d` is then projected onto the variable bounds. The correction is
!  not used if it is larger than the step itself (it is meant to be a small
!  correction; a large one comes from a nearly singular \( J_S \)), or if
!  `LSQR` stops at its iteration limit.
!
!@note Variables at a bound at `x+p` are not held fixed in the solve for
!      `d`, so the final projection can undo part of the correction. (Holding
!      them fixed, or re-solving with the variables that the correction takes
!      outside their bounds held fixed, was tried: it saves function
!      evaluations on some problems, but makes the correction large where
!      the free variables' columns of \( J_S \) are small, e.g. TP13, and
!      sends TP116 to a different local solution.)

    module sqpopt_soc_module

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_sparse_matrix, sqpopt_all_finite
    use sqpopt_linalg_module, only: sparse_matvec
    use lsqr_module,          only: lsqr_solver_ez
    use sqpopt_least_squares_module, only: sqpopt_least_squares_type

    implicit none

    private

    public :: soc_step

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  compute the second-order-corrected step `p_soc` for the trial step `p`
!  (see the module-level documentation). `ok` is false (and `p_soc=p`) if
!  there is nothing to correct or the correction is not usable. With
!  `least_squares` (if it is enabled), the correction is computed by a
!  direct solve, and by `LSQR` only if that fails.

    subroutine soc_step(jac, x, p, c, c_trial, c_lb, c_ub, x_lb, x_ub, p_soc, ok, least_squares)

    type(sqpopt_sparse_matrix), intent(in)  :: jac     !! constraint Jacobian at `x`, `dimension(m,n)`
    real(wp), dimension(:),     intent(in)  :: x       !! current point `dimension(n)`
    real(wp), dimension(:),     intent(in)  :: p       !! the rejected trial step `dimension(n)`
    real(wp), dimension(:),     intent(in)  :: c       !! constraint values at `x` `dimension(m)`
    real(wp), dimension(:),     intent(in)  :: c_trial !! constraint values at `x+p` `dimension(m)`
    real(wp), dimension(:),     intent(in)  :: c_lb    !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)  :: c_ub    !! constraint upper bounds `dimension(m)`
    real(wp), dimension(:),     intent(in)  :: x_lb    !! variable lower bounds `dimension(n)`
    real(wp), dimension(:),     intent(in)  :: x_ub    !! variable upper bounds `dimension(n)`
    real(wp), dimension(:),     intent(out) :: p_soc   !! corrected step `dimension(n)`
    logical,                    intent(out) :: ok      !! true if `p_soc` is a usable corrected step
    type(sqpopt_least_squares_type), optional, intent(inout) :: least_squares !! the direct least-squares solver
                                                                              !! (see [[sqpopt_least_squares_module]])

    real(wp), parameter :: act_tol = 1.0e-6_wp !! relative tolerance for "at a bound" in the linearization
    real(wp), parameter :: soc_max_ratio = 1.0_wp !! the correction is rejected if \( \lVert d \rVert > \)
                                                  !! `soc_max_ratio` \( \lVert p \rVert \)
    integer,  parameter :: lsqr_itnlim_stop = 5 !! `LSQR`'s `istop` for "iteration limit reached"

    real(wp), dimension(size(c)) :: c_lin, resid
    integer,  dimension(size(c)) :: row_map
    real(wp), dimension(size(x)) :: d
    real(wp), dimension(:), allocatable :: rhs
    integer,  dimension(:), allocatable :: irow, icol
    real(wp), dimension(:), allocatable :: val
    type(lsqr_solver_ez) :: lsqr
    integer :: i, k, m_s, nnz_s, istop
    logical :: solved

    p_soc = p
    ok    = .false.
    if (size(c) == 0) return

    call sparse_matvec(jac, p, c_lin)
    c_lin = c + c_lin

    ! select the rows to correct, and the residual of each from its bound:
    row_map = 0
    m_s = 0
    do i = 1, size(c)
        if (abs(c_lin(i)-c_lb(i)) <= act_tol*max(1.0_wp, abs(c_lb(i)))) then
            resid(i) = c_trial(i) - c_lb(i)      ! active at the lower bound (or an equality)
        else if (abs(c_lin(i)-c_ub(i)) <= act_tol*max(1.0_wp, abs(c_ub(i)))) then
            resid(i) = c_trial(i) - c_ub(i)      ! active at the upper bound
        else if (c_trial(i) < c_lb(i)) then
            resid(i) = c_trial(i) - c_lb(i)      ! inactive in the linearization, but violated at x+p
        else if (c_trial(i) > c_ub(i)) then
            resid(i) = c_trial(i) - c_ub(i)
        else
            cycle
        end if
        m_s = m_s + 1
        row_map(i) = m_s
    end do
    if (m_s == 0) return

    ! the minimum-norm correction, by a direct solve if there is one:
    solved = .false.
    if (present(least_squares)) call least_squares%min_norm(jac, row_map > 0, -resid, d, solved)

    if (.not. solved) then
        ! the sub-Jacobian of the selected rows, and the right-hand side:
        nnz_s = count(row_map(jac%irow(1:jac%nnz)) > 0)
        allocate(irow(nnz_s), icol(nnz_s), val(nnz_s), rhs(m_s))
        k = 0
        do i = 1, jac%nnz
            if (row_map(jac%irow(i)) > 0) then
                k = k + 1
                irow(k) = row_map(jac%irow(i))
                icol(k) = jac%icol(i)
                val(k)  = jac%val(i)
            end if
        end do
        do i = 1, size(c)
            if (row_map(i) > 0) rhs(row_map(i)) = -resid(i)
        end do

        call lsqr%initialize(m_s, size(x), val, irow, icol, itnlim=2*(m_s+size(x))+10)
        call lsqr%solve(rhs, 0.0_wp, d, istop)
        if (istop == lsqr_itnlim_stop .or. .not. sqpopt_all_finite(d)) return
    end if
    if (norm2(d) > soc_max_ratio*norm2(p)) return  ! (the correction is meant to be small relative to the step)

    p_soc = min(max(x + p + d, x_lb), x_ub) - x
    ok    = .true.

    end subroutine soc_step
!*******************************************************************************

    end module sqpopt_soc_module
!*******************************************************************************
