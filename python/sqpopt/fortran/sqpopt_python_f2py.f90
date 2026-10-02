!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The routines f2py wraps for the Python bindings (see `python/sqpopt`),
!  as described to it by the signature file `_sqpopt.pyf` (so that f2py
!  never reads the Fortran sources): external procedures with explicit-shape
!  arrays, which forward to the `sqpopt_python` module.

!*******************************************************************************
!>
!  the sizes of `iinfo` and `rinfo`, and the number and signature of the
!  options (see `sqpopt_python_n_info`)

    subroutine sqpopt_py_info(n_iinfo, n_rinfo, n_options, options_signature)

    use, intrinsic :: iso_fortran_env, only: int32
    use sqpopt_python, only: sqpopt_python_n_info

    implicit none

    integer(int32), intent(out) :: n_iinfo           !! size of `iinfo`
    integer(int32), intent(out) :: n_rinfo           !! size of `rinfo`
    integer(int32), intent(out) :: n_options         !! number of options
    integer(int32), intent(out) :: options_signature !! checksum of the options' names and kinds

    call sqpopt_python_n_info(n_iinfo, n_rinfo, n_options, options_signature)

    end subroutine sqpopt_py_info
!*******************************************************************************

!*******************************************************************************
!>
!  solve a problem (see `sqpopt_python_solve`, which this forwards to). The
!  sizes are given explicitly, and the Python caller passes every array with
!  at least one element (see `_sqpopt.pyf`); only the first `n`, `m`, ...
!  elements are used.

    subroutine sqpopt_py_solve(fc, gjac, hess, report, use_hess, use_report, n, m, nnz, nnzh, nopt, &
                               x0, x_lb, x_ub, c_lb, c_ub, jac_irow, jac_icol, hess_irow, hess_icol, &
                               opt_id, opt_val, lambda0, use_lambda0, max_step, use_max_step, output_file, &
                               diagnostics_file, x, lambda, z, c, iinfo, rinfo, message, diag_index, diag_value, &
                               report_text, problem_report)

    use, intrinsic :: iso_fortran_env, only: real64, int32
    use sqpopt_python, only: sqpopt_python_solve, py_fc_func, py_gjac_func, py_hess_func, py_report_func, &
                             sqpopt_python_n_iinfo, sqpopt_python_n_rinfo, sqpopt_python_n_diag, &
                             sqpopt_python_len_report

    implicit none

    procedure(py_fc_func)     :: fc     !! objective and constraints
    procedure(py_gjac_func)   :: gjac   !! gradient and Jacobian values
    procedure(py_hess_func)   :: hess   !! Hessian of the Lagrangian values (only called if `use_hess /= 0`)
    procedure(py_report_func) :: report !! progress callback (only called if `use_report /= 0`)
    integer(int32),                  intent(in)  :: use_hess    !! whether to use `hess` (and the Hessian pattern)
    integer(int32),                  intent(in)  :: use_report  !! whether to call `report`
    integer(int32),                  intent(in)  :: n           !! number of variables
    integer(int32),                  intent(in)  :: m           !! number of constraints
    integer(int32),                  intent(in)  :: nnz         !! number of Jacobian nonzeros
    integer(int32),                  intent(in)  :: nnzh        !! number of Hessian nonzeros
    integer(int32),                  intent(in)  :: nopt        !! number of options given
    real(real64),   dimension(n),    intent(in)  :: x0          !! starting point
    real(real64),   dimension(n),    intent(in)  :: x_lb        !! variable lower bounds
    real(real64),   dimension(n),    intent(in)  :: x_ub        !! variable upper bounds
    real(real64),   dimension(m),    intent(in)  :: c_lb        !! constraint lower bounds
    real(real64),   dimension(m),    intent(in)  :: c_ub        !! constraint upper bounds
    integer(int32), dimension(nnz),  intent(in)  :: jac_irow    !! Jacobian pattern: row indices (1-based)
    integer(int32), dimension(nnz),  intent(in)  :: jac_icol    !! Jacobian pattern: column indices (1-based)
    integer(int32), dimension(nnzh), intent(in)  :: hess_irow   !! Hessian pattern: row indices (1-based; each
                                                                !! off-diagonal element once)
    integer(int32), dimension(nnzh), intent(in)  :: hess_icol   !! Hessian pattern: column indices (1-based)
    integer(int32), dimension(nopt), intent(in)  :: opt_id      !! options to set, by number
    real(real64),   dimension(nopt), intent(in)  :: opt_val     !! their values (integers and logicals as reals)
    real(real64),   dimension(m),    intent(in)  :: lambda0     !! starting constraint multipliers (only used if
                                                                !! `use_lambda0 /= 0`)
    integer(int32),                  intent(in)  :: use_lambda0 !! whether to start from `lambda0`
    real(real64),   dimension(n),    intent(in)  :: max_step    !! the largest change of each variable in one major
                                                                !! iteration (only used if `use_max_step /= 0`)
    integer(int32),                  intent(in)  :: use_max_step !! whether to limit the steps by `max_step`
    character(len=*),                intent(in)  :: output_file !! file for the printed output (replaced; blank:
                                                                !! `options%output_unit`)
    character(len=*),                intent(in)  :: diagnostics_file !! file for the history of the iterations
                                                                !! (replaced; blank: none)
    real(real64),   dimension(n),    intent(out) :: x           !! solution
    real(real64),   dimension(m),    intent(out) :: lambda      !! constraint multipliers
    real(real64),   dimension(n),    intent(out) :: z           !! variable-bound multipliers
    real(real64),   dimension(m),    intent(out) :: c           !! constraint values at `x`
    integer(int32), dimension(sqpopt_python_n_iinfo), intent(out) :: iinfo !! integer results (see `sqpopt_python_solve`)
    real(real64),   dimension(sqpopt_python_n_rinfo), intent(out) :: rinfo !! real results (see `sqpopt_python_solve`)
    character(len=256),              intent(out) :: message     !! description of the status
    integer(int32), dimension(sqpopt_python_n_diag), intent(out) :: diag_index !! the lists of the diagnosis: indices
                                                                               !! (see `sqpopt_python_solve`)
    real(real64),   dimension(sqpopt_python_n_diag), intent(out) :: diag_value !! the lists of the diagnosis: values
    character(len=sqpopt_python_len_report), intent(out) :: report_text    !! the diagnosis of the solve
    character(len=sqpopt_python_len_report), intent(out) :: problem_report !! the diagnostics' report on the
                                                                           !! starting point

    x = x0
    lambda = 0.0_real64
    z = 0.0_real64
    c = 0.0_real64
    call sqpopt_python_solve(fc, gjac, hess, report, use_hess, use_report, x0, x_lb, x_ub, c_lb, c_ub, &
                             jac_irow, jac_icol, hess_irow, hess_icol, opt_id, opt_val, &
                             lambda0, use_lambda0, max_step, use_max_step, output_file, diagnostics_file, &
                             x, lambda, z, c, iinfo, rinfo, message, diag_index, diag_value, report_text, problem_report)

    end subroutine sqpopt_py_solve
!*******************************************************************************
