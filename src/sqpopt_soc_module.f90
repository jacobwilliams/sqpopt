!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Second-order correction (SOC) step, shared by the line-search-based
!  major iteration ([[sqpopt_iterate_module]]) and the trust-region-based
!  one ([[sqpopt_trust_region_module]]).

    module sqpopt_soc_module

    use sqpopt_kinds,              only: wp => sqpopt_module_wp
    use sqpopt_types_module,       only: sqpopt_sparse_matrix
    use sqpopt_problem_module,     only: sqpopt_problem_type
    use sqpopt_linesearch_module,  only: sqpopt_linesearch_type
    use sqpopt_linalg_module,      only: sparse_matvec
    use lsqr_module,                only: lsqr_solver_ez

    implicit none

    private

    public :: second_order_correction

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  second-order correction: refine `p` using the true (nonlinear)
!  constraint residual at `x+p` rather than its linear prediction, and
!  replace `p` with the corrected step if it improves the merit function
!  (a standard remedy for the Maratos effect near nonlinear constraints).

    subroutine second_order_correction(problem, linesearch, jac, x, c, lambda, p)

    type(sqpopt_problem_type),    intent(inout) :: problem
    type(sqpopt_linesearch_type), intent(inout) :: linesearch
    type(sqpopt_sparse_matrix),   intent(in)    :: jac
    real(wp), dimension(:),       intent(in)    :: x
    real(wp), dimension(:),       intent(in)    :: c
    real(wp), dimension(:),       intent(in)    :: lambda
    real(wp), dimension(:),       intent(inout) :: p

    real(wp), dimension(size(c)) :: jp, c_p, resid, c_soc
    real(wp), dimension(size(p)) :: p_corr, p_soc
    real(wp) :: f_p, f_soc, phi_p, phi_soc
    type(lsqr_solver_ez) :: lsqr
    integer :: istop

    call sparse_matvec(jac, p, jp)
    call problem%eval_c(x+p, c_p)
    resid = c_p - (c+jp)  !! nonlinear residual left uncorrected by the linear model

    call lsqr%initialize(problem%m, problem%n, jac%val, jac%irow, jac%icol)
    call lsqr%solve(-resid, 0.0_wp, p_corr, istop)
    p_soc = p + p_corr

    call problem%eval_f(x+p, f_p)
    call linesearch%eval_merit(f_p, c_p, problem%c_lb, problem%c_ub, lambda, phi_p)

    call problem%eval_f(x+p_soc, f_soc)
    call problem%eval_c(x+p_soc, c_soc)
    call linesearch%eval_merit(f_soc, c_soc, problem%c_lb, problem%c_ub, lambda, phi_soc)

    if (phi_soc < phi_p) p = p_soc

    end subroutine second_order_correction
!*******************************************************************************

    end module sqpopt_soc_module
!*******************************************************************************
