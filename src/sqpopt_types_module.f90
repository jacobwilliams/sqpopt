!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Common types and status codes used throughout `sqpopt`. Sparse matrices
!  are always stored in coordinate (COO) format with 1-based indices --
!  the same convention used by the `lusol`, `LSQR`, and `LSMR` dependencies
!  -- so that they may be passed to those solvers without conversion.

    module sqpopt_types_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    ! solver status/exit codes:
    integer, parameter, public :: sqpopt_success              = 0  !! converged successfully
    integer, parameter, public :: sqpopt_max_iter_reached     = 1  !! stopped: maximum number of iterations reached
    integer, parameter, public :: sqpopt_infeasible           = 2  !! stopped: problem appears to be infeasible
    integer, parameter, public :: sqpopt_line_search_failed   = 3  !! stopped: line search failed to find an acceptable step
    integer, parameter, public :: sqpopt_qp_solve_failed      = 4  !! stopped: QP subproblem solver failed
    integer, parameter, public :: sqpopt_user_requested_stop  = 5  !! stopped: user requested stop
    integer, parameter, public :: sqpopt_error                = -1 !! stopped: an unspecified error occurred

    type, public :: sqpopt_sparse_matrix
        !! a sparse matrix stored in coordinate (COO) format
        integer :: nrows = 0  !! number of rows
        integer :: ncols = 0  !! number of columns
        integer :: nnz   = 0  !! number of nonzero elements
        integer,  dimension(:), allocatable :: irow  !! row indices `dimension(nnz)`
        integer,  dimension(:), allocatable :: icol  !! column indices `dimension(nnz)`
        real(wp), dimension(:), allocatable :: val   !! nonzero values `dimension(nnz)`
    end type sqpopt_sparse_matrix

    abstract interface
        subroutine sqpopt_report_func(iter, x, f, c, lambda, user_stop)
            !! user-supplied callback invoked once per major SQP iteration
            !! for progress monitoring (see `sqpopt_type%initialize`'s
            !! `report` argument). Set `user_stop=.true.` to have the
            !! solver stop after the current iteration
            !! (`istat=sqpopt_user_requested_stop`).
            import :: wp
            implicit none
            integer,                intent(in)  :: iter      !! major iteration number (starts at 1)
            real(wp), dimension(:), intent(in)  :: x         !! current optimization variables `dimension(n)`
            real(wp),               intent(in)  :: f         !! current objective function value
            real(wp), dimension(:), intent(in)  :: c         !! current constraint values `dimension(m)`
            real(wp), dimension(:), intent(in)  :: lambda    !! current Lagrange multiplier estimate `dimension(m)`
            logical,                intent(out) :: user_stop !! set `.true.` to request the solver stop
        end subroutine sqpopt_report_func
    end interface

    public :: sqpopt_report_func

    end module sqpopt_types_module
!*******************************************************************************
