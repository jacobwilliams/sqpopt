!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Common types and status codes used throughout `sqpopt`.

    module sqpopt_types_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    ! solver status/exit codes:
    integer, parameter, public :: sqpopt_success             = 0   !! converged successfully
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

    end module sqpopt_types_module
!*******************************************************************************
