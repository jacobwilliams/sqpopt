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
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    private

    ! solver status/exit codes:
    integer, parameter, public :: sqpopt_success              = 0  !! converged successfully
    integer, parameter, public :: sqpopt_max_iter_reached     = 1  !! stopped: maximum number of iterations reached
    integer, parameter, public :: sqpopt_infeasible           = 2  !! stopped: problem appears to be infeasible
    integer, parameter, public :: sqpopt_line_search_failed   = 3  !! stopped: line search failed to find an acceptable step
    integer, parameter, public :: sqpopt_qp_solve_failed      = 4  !! stopped: QP subproblem solver failed
    integer, parameter, public :: sqpopt_user_requested_stop  = 5  !! stopped: user requested stop
    integer, parameter, public :: sqpopt_invalid_input        = 6  !! stopped: the problem definition or options are invalid
    integer, parameter, public :: sqpopt_stalled              = 7  !! stopped: the point is feasible but the objective and variables
                                                                    !! have stopped changing (see `ftol`/`xtol`) before the KKT test
                                                                    !! was satisfied -- usually an acceptable, if imprecise, solution
    integer, parameter, public :: sqpopt_function_error       = 8  !! stopped: the problem functions returned a non-finite value
                                                                    !! (NaN or Inf) at the current point (trial points with
                                                                    !! non-finite values are rejected instead)
    integer, parameter, public :: sqpopt_error                = -1 !! stopped: an unspecified error occurred

    real(wp), parameter, public :: sqpopt_infinity = 1.0e20_wp !! any bound with magnitude `>= sqpopt_infinity` is treated
                                                                !! as absent (bounds are clamped to `[-sqpopt_infinity,
                                                                !! sqpopt_infinity]` so that e.g. `-huge(1.0_wp)` is safe to use)

    public :: sqpopt_status_message
    public :: sqpopt_all_finite

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

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  a human-readable description of a solver status code.

    pure function sqpopt_status_message(istat) result(msg)

    integer, intent(in) :: istat !! status code (one of the `sqpopt_*` status parameters)
    character(len=:), allocatable :: msg

    select case (istat)
    case (sqpopt_success);             msg = 'converged successfully'
    case (sqpopt_max_iter_reached);    msg = 'maximum number of iterations reached'
    case (sqpopt_infeasible);          msg = 'problem appears to be (locally) infeasible'
    case (sqpopt_line_search_failed);  msg = 'line search failed to find an acceptable step'
    case (sqpopt_qp_solve_failed);     msg = 'QP subproblem solver failed'
    case (sqpopt_user_requested_stop); msg = 'user requested stop'
    case (sqpopt_invalid_input);       msg = 'invalid problem definition or options'
    case (sqpopt_stalled);             msg = 'feasible, but no further progress is being made'
    case (sqpopt_function_error);      msg = 'the problem functions returned a non-finite value (NaN or Inf)'
    case (sqpopt_error);               msg = 'an unspecified error occurred'
    case default;                      msg = 'unknown status code'
    end select

    end function sqpopt_status_message
!*******************************************************************************

!*******************************************************************************
!>
!  true if every element of `v` is finite (not NaN or +/-Inf). Uses
!  `ieee_is_finite`, which (unlike an ordinary comparison) does not
!  raise an IEEE invalid exception on a NaN.

    pure logical function sqpopt_all_finite(v)

    real(wp), dimension(:), intent(in) :: v

    sqpopt_all_finite = all(ieee_is_finite(v))

    end function sqpopt_all_finite
!*******************************************************************************

    end module sqpopt_types_module
!*******************************************************************************
