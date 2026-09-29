!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The Fortran side of the Python bindings (see `python/sqpopt`), built into
!  a Python extension by [PRIK](https://github.com/PyNumLab/prik).
!
!  PRIK can't wrap the library's API directly: `set_functions` keeps
!  procedure pointers to the user functions after it returns, and passes
!  `class(*)` user data, neither of which it supports. So this module
!  exposes a single procedure, [[solve]], whose arguments are all
!  intrinsic-type scalars and arrays plus four callbacks (Python callables),
!  that are only called during the solve. The callbacks are forwarded to
!  the library through module (private) procedure pointers, which are saved
!  and restored around each solve, so that solves may be nested (e.g. a
!  callback that itself calls `solve`).
!
!  The implementation is in `sqpopt_python_core` (this module only
!  forwards to it, so that PRIK doesn't need to read it), and the
!  callbacks' interfaces are in `sqpopt_python_interfaces`.

    module sqpopt_python

    use, intrinsic :: iso_fortran_env, only: real64, int32
    use sqpopt_python_interfaces, only: py_fc_func, py_gjac_func, py_hess_func, py_report_func

    implicit none

    private

    public :: solve, n_info

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  the number of values [[solve]] returns in `iinfo` and `rinfo`, and the
!  signature of the options numbering, which the Python side checks
!  against its own

    subroutine n_info(n_iinfo, n_rinfo, n_options, options_signature)

    use sqpopt_python_core, only: sqpopt_python_n_info

    integer(int32), intent(out) :: n_iinfo           !! size of `iinfo`
    integer(int32), intent(out) :: n_rinfo           !! size of `rinfo`
    integer(int32), intent(out) :: n_options         !! number of options
    integer(int32), intent(out) :: options_signature !! checksum of the options' names and kinds

    call sqpopt_python_n_info(n_iinfo, n_rinfo, n_options, options_signature)

    end subroutine n_info
!*******************************************************************************

!*******************************************************************************
!>
!  solve a problem (see `sqpopt_python_core`'s `sqpopt_python_solve` for the
!  outputs `iinfo` and `rinfo`). The outputs `x`, `lambda`, `z`, `c`,
!  `iinfo`, and `rinfo` must be allocated by the caller.

    subroutine solve(fc, gjac, hess, report, use_hess, use_report, x0, x_lb, x_ub, c_lb, c_ub, &
                     jac_irow, jac_icol, hess_irow, hess_icol, opt_id, opt_val, &
                     lambda0, use_lambda0, output_file, x, lambda, z, c, iinfo, rinfo, message)

    use sqpopt_python_core, only: sqpopt_python_solve

    procedure(py_fc_func)     :: fc     !! objective and constraints
    procedure(py_gjac_func)   :: gjac   !! gradient and Jacobian values
    procedure(py_hess_func)   :: hess   !! Hessian of the Lagrangian values (only called if `use_hess /= 0`)
    procedure(py_report_func) :: report !! progress callback (only called if `use_report /= 0`)
    integer(int32),                 intent(in)    :: use_hess   !! whether to use `hess` (and the exact-Hessian pattern)
    integer(int32),                 intent(in)    :: use_report !! whether to call `report`
    real(real64),   dimension(:),   intent(in)    :: x0         !! starting point `dimension(n)`
    real(real64),   dimension(:),   intent(in)    :: x_lb       !! variable lower bounds `dimension(n)`
    real(real64),   dimension(:),   intent(in)    :: x_ub       !! variable upper bounds `dimension(n)`
    real(real64),   dimension(:),   intent(in)    :: c_lb       !! constraint lower bounds `dimension(m)`
    real(real64),   dimension(:),   intent(in)    :: c_ub       !! constraint upper bounds `dimension(m)`
    integer(int32), dimension(:),   intent(in)    :: jac_irow   !! Jacobian pattern: row indices (1-based)
    integer(int32), dimension(:),   intent(in)    :: jac_icol   !! Jacobian pattern: column indices (1-based)
    integer(int32), dimension(:),   intent(in)    :: hess_irow  !! Hessian pattern: row indices (1-based; each
                                                                !! off-diagonal element once)
    integer(int32), dimension(:),   intent(in)    :: hess_icol  !! Hessian pattern: column indices (1-based)
    integer(int32), dimension(:),   intent(in)    :: opt_id     !! options to set, by number
    real(real64),   dimension(:),   intent(in)    :: opt_val    !! their values (integers and logicals as reals)
    real(real64),   dimension(:),   intent(in)    :: lambda0    !! starting constraint multipliers `dimension(m)`
                                                                !! (only used if `use_lambda0 /= 0`)
    integer(int32),                 intent(in)    :: use_lambda0 !! whether to start from `lambda0`
    character(len=*),               intent(in)    :: output_file !! file for the printed output (replaced; blank:
                                                                 !! `options%output_unit`)
    real(real64),   dimension(:),   intent(inout) :: x          !! solution `dimension(n)`
    real(real64),   dimension(:),   intent(inout) :: lambda     !! constraint multipliers `dimension(m)`
    real(real64),   dimension(:),   intent(inout) :: z          !! variable-bound multipliers `dimension(n)`
    real(real64),   dimension(:),   intent(inout) :: c          !! constraint values at `x` `dimension(m)`
    integer(int32), dimension(:),   intent(inout) :: iinfo      !! integer results
    real(real64),   dimension(:),   intent(inout) :: rinfo      !! real results
    character(len=256),             intent(out)   :: message    !! description of the status

    call sqpopt_python_solve(fc, gjac, hess, report, use_hess, use_report, x0, x_lb, x_ub, c_lb, c_ub, &
                             jac_irow, jac_icol, hess_irow, hess_icol, opt_id, opt_val, &
                             lambda0, use_lambda0, output_file, x, lambda, z, c, iinfo, rinfo, message)

    end subroutine solve
!*******************************************************************************

    end module sqpopt_python
!*******************************************************************************
