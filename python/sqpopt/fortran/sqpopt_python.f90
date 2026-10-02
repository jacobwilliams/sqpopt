!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The Fortran side of the Python bindings (see `python/sqpopt`), whose
!  extension is built by f2py (the wrapped routines are in
!  `sqpopt_python_f2py.f90`, described to f2py by `_sqpopt.pyf`).
!
!  f2py can't wrap the library's API directly: `set_functions` keeps
!  procedure pointers to the user functions after it returns, and passes
!  `class(*)` user data. So this module provides a single solve
!  ([[sqpopt_python_solve]]), whose arguments are all intrinsic-type
!  scalars and arrays plus four callbacks (Python callables), which are only
!  called during the solve. The callbacks are forwarded to the library
!  through module (private) procedure pointers, which are saved and
!  restored around each solve, so that solves may be nested (e.g. a
!  callback that itself calls `solve`).
!
!  The callbacks have f2py's calling convention: explicit-shape arrays with
!  their sizes, the `status` (or `stop`) flag first (the only value the
!  Python function returns), and the outputs written in place (a scalar
!  output is a 1-element array), which also works for zero-size arrays.
!
!  The options are set by number (see `sqpopt_python_options`, which is
!  generated from the Python options schema by the build).

    module sqpopt_python

    use, intrinsic :: iso_fortran_env, only: real64, int32, output_unit

    implicit none

    private

    integer(int32), parameter, public :: sqpopt_python_n_iinfo = 17 !! size of `iinfo` (see [[sqpopt_python_solve]])
    integer(int32), parameter, public :: sqpopt_python_n_rinfo = 8  !! size of `rinfo`

    abstract interface
        subroutine py_fc_func(n, m, status, x, f, c)
            !! the objective and the constraints
            import :: real64, int32
            implicit none
            integer(int32),               intent(in)    :: n      !! number of variables
            integer(int32),               intent(in)    :: m      !! number of constraints
            integer(int32),               intent(inout) :: status !! `0` on entry; `> 0`: can't evaluate at `x`;
                                                                  !! `< 0`: stop
            real(real64), dimension(n),   intent(in)    :: x      !! point
            real(real64), dimension(1),   intent(inout) :: f      !! objective value at `x`
            real(real64), dimension(m),   intent(inout) :: c      !! constraint values at `x`
        end subroutine py_fc_func
        subroutine py_gjac_func(n, nnz, status, accuracy, x, g, jac_val)
            !! the objective gradient and the nonzero values of the constraint Jacobian
            import :: real64, int32
            implicit none
            integer(int32),               intent(in)    :: n        !! number of variables
            integer(int32),               intent(in)    :: nnz      !! number of Jacobian nonzeros
            integer(int32),               intent(inout) :: status   !! `0` on entry; `> 0`: can't evaluate at `x`;
                                                                    !! `< 0`: stop
            integer(int32),               intent(in)    :: accuracy !! requested accuracy (`1` fast, `2` accurate)
            real(real64), dimension(n),   intent(in)    :: x        !! point
            real(real64), dimension(n),   intent(inout) :: g        !! objective gradient at `x`
            real(real64), dimension(nnz), intent(inout) :: jac_val  !! Jacobian values at `x`, in the pattern's order
        end subroutine py_gjac_func
        subroutine py_hess_func(n, m, nnz, status, x, lambda, hess_val)
            !! the nonzero values of the Hessian of the Lagrangian
            import :: real64, int32
            implicit none
            integer(int32),               intent(in)    :: n        !! number of variables
            integer(int32),               intent(in)    :: m        !! number of constraints
            integer(int32),               intent(in)    :: nnz      !! number of Hessian nonzeros
            integer(int32),               intent(inout) :: status   !! `0` on entry; `> 0`: can't evaluate at `x`;
                                                                    !! `< 0`: stop
            real(real64), dimension(n),   intent(in)    :: x        !! point
            real(real64), dimension(m),   intent(in)    :: lambda   !! constraint multipliers
            real(real64), dimension(nnz), intent(inout) :: hess_val !! Hessian values at `x`, in the pattern's order
        end subroutine py_hess_func
        subroutine py_report_func(n, m, stop, iter, x, f, c, lambda)
            !! called once per major iteration
            import :: real64, int32
            implicit none
            integer(int32),               intent(in)    :: n      !! number of variables
            integer(int32),               intent(in)    :: m      !! number of constraints
            integer(int32),               intent(inout) :: stop   !! `0` on entry; set nonzero to stop the solver
            integer(int32),               intent(in)    :: iter   !! major iteration number
            real(real64), dimension(n),   intent(in)    :: x      !! point
            real(real64),                 intent(in)    :: f      !! objective value at `x`
            real(real64), dimension(m),   intent(in)    :: c      !! constraint values at `x`
            real(real64), dimension(m),   intent(in)    :: lambda !! constraint multipliers
        end subroutine py_report_func
    end interface

    procedure(py_fc_func),     pointer :: py_fc     => null() !! the callbacks of the current solve
    procedure(py_gjac_func),   pointer :: py_gjac   => null()
    procedure(py_hess_func),   pointer :: py_hess   => null()
    procedure(py_report_func), pointer :: py_report => null()

    public :: py_fc_func, py_gjac_func, py_hess_func, py_report_func
    public :: sqpopt_python_solve, sqpopt_python_n_info

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  the number of values [[sqpopt_python_solve]] returns in `iinfo` and `rinfo`, and the
!  signature of the options numbering (see `sqpopt_python_options`), which
!  the Python side checks against its own

    subroutine sqpopt_python_n_info(n_iinfo, n_rinfo, n_options, options_signature)

    use sqpopt_python_options, only: sqpopt_python_n_options, sqpopt_python_options_signature

    integer(int32), intent(out) :: n_iinfo           !! size of `iinfo`
    integer(int32), intent(out) :: n_rinfo           !! size of `rinfo`
    integer(int32), intent(out) :: n_options         !! number of options
    integer(int32), intent(out) :: options_signature !! checksum of the options' names and kinds

    n_iinfo = sqpopt_python_n_iinfo
    n_rinfo = sqpopt_python_n_rinfo
    n_options = sqpopt_python_n_options
    options_signature = sqpopt_python_options_signature

    end subroutine sqpopt_python_n_info
!*******************************************************************************

!*******************************************************************************
!>
!  solve a problem (see the library's `sqpopt_type%solve`). The outputs
!  `x`, `lambda`, `z`, `c`, `iinfo`, and `rinfo` must be allocated by the
!  caller (`iinfo` and `rinfo` with `sqpopt_python_n_iinfo` and
!  `sqpopt_python_n_rinfo` elements).
!
!  `iinfo`: `istat`, `iterations`, `n_eval_fc`, `n_eval_gjac`,
!  `n_eval_hess`, `n_qp_iterations`, `derivative_switch_iteration`,
!  `n_soc`, `n_hessian_resets`, `n_restoration_steps`,
!  `n_restoration_phases`, `n_elastic`, `n_escape`, `n_factorizations`,
!  `n_qp_solves`, `n_direct_qp`, `n_unconstrained_qp`.
!
!  `rinfo`: `f`, `kkt_error`, `feasibility_error`, `stationarity_error`,
!  `time`, `time_functions`, `time_qp`, `time_factorization`.

    subroutine sqpopt_python_solve(fc, gjac, hess, report, use_hess, use_report, x0, x_lb, x_ub, c_lb, c_ub, &
                     jac_irow, jac_icol, hess_irow, hess_icol, opt_id, opt_val, &
                     lambda0, use_lambda0, max_step, use_max_step, output_file, &
                     x, lambda, z, c, iinfo, rinfo, message)

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_types_module,     only: sqpopt_results_type, sqpopt_invalid_input
    use sqpopt_python_options,   only: sqpopt_python_set_option

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
    integer(int32), dimension(:),   intent(in)    :: opt_id     !! options to set, by number (see `sqpopt_python_options`)
    real(real64),   dimension(:),   intent(in)    :: opt_val    !! their values (integers and logicals as reals)
    real(real64),   dimension(:),   intent(in)    :: lambda0    !! starting constraint multipliers `dimension(m)`
                                                                !! (only used if `use_lambda0 /= 0`)
    integer(int32),                 intent(in)    :: use_lambda0 !! whether to start from `lambda0`
    real(real64),   dimension(:),   intent(in)    :: max_step   !! the largest change of each variable in one major
                                                                !! iteration `dimension(n)` (only used if
                                                                !! `use_max_step /= 0`)
    integer(int32),                 intent(in)    :: use_max_step !! whether to limit the steps by `max_step`
    character(len=*),               intent(in)    :: output_file !! file for the printed output (replaced; blank:
                                                                 !! `options%output_unit`)
    real(real64),   dimension(:),   intent(inout) :: x          !! solution `dimension(n)`
    real(real64),   dimension(:),   intent(inout) :: lambda     !! constraint multipliers `dimension(m)`
    real(real64),   dimension(:),   intent(inout) :: z          !! variable-bound multipliers `dimension(n)`
    real(real64),   dimension(:),   intent(inout) :: c          !! constraint values at `x` `dimension(m)`
    integer(int32), dimension(:),   intent(inout) :: iinfo      !! integer results (see above)
    real(real64),   dimension(:),   intent(inout) :: rinfo      !! real results (see above)
    character(len=256),             intent(out)   :: message    !! description of the status

    type(sqpopt_type)              :: solver
    type(sqpopt_problem_type)      :: problem
    type(sqpopt_options_type)      :: options
    type(sqpopt_hessian_type)      :: hessian
    type(sqpopt_qp_solver_type)    :: qp_solver
    type(sqpopt_linesearch_type)   :: linesearch
    type(sqpopt_trust_region_type) :: trust_region
    type(sqpopt_results_type)      :: r
    procedure(py_fc_func),     pointer :: fc0
    procedure(py_gjac_func),   pointer :: gjac0
    procedure(py_hess_func),   pointer :: hess0
    procedure(py_report_func), pointer :: report0
    integer :: istat, i, n, m
    integer :: unit !! the unit of `output_file` (`-1` if none)
    integer :: ios
    logical :: ok

    n = size(x0)
    m = size(c_lb)
    unit = -1
    iinfo = 0
    rinfo = 0.0_real64

    ! (the callbacks of an enclosing solve, restored on exit)
    fc0 => py_fc
    gjac0 => py_gjac
    hess0 => py_hess
    report0 => py_report
    py_fc => fc
    py_gjac => gjac
    py_hess => hess
    py_report => report

    do i = 1, min(size(opt_id), size(opt_val))
        call sqpopt_python_set_option(opt_id(i), opt_val(i), options, hessian, qp_solver, linesearch, &
                                      trust_region, ok)
        if (.not. ok) then
            iinfo(1) = sqpopt_invalid_input
            message = 'invalid option number'
            call restore()
            return
        end if
    end do

    if (len_trim(output_file) > 0) then
        open(newunit=unit, file=trim(output_file), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            unit = -1
            iinfo(1) = sqpopt_invalid_input
            message = 'output_file can''t be opened: '//trim(output_file)
            call restore()
            return
        end if
        options%output_unit = unit
    end if

    call problem%set_problem_size(n=n, m=m)
    call problem%set_bounds(x_lb, x_ub, c_lb, c_ub)
    if (use_max_step /= 0) call problem%set_max_step(max_step)
    call problem%set_jacobian_sparsity(size(jac_irow), jac_irow, jac_icol)
    if (use_hess /= 0) then
        call problem%set_functions(fc=fc_bridge, gjac=gjac_bridge, hess=hess_bridge)
        call problem%set_hessian_sparsity(size(hess_irow), hess_irow, hess_icol)
    else
        call problem%set_functions(fc=fc_bridge, gjac=gjac_bridge)
    end if

    if (use_report /= 0) then
        call solver%initialize(problem=problem, options=options, hessian=hessian, qp_solver=qp_solver, &
                               linesearch=linesearch, trust_region=trust_region, report=report_bridge)
    else
        call solver%initialize(problem=problem, options=options, hessian=hessian, qp_solver=qp_solver, &
                               linesearch=linesearch, trust_region=trust_region)
    end if
    if (use_lambda0 /= 0) then
        call solver%solve(x0, istat, lambda0)
    else
        call solver%solve(x0, istat)
    end if
    call solver%get_results(r)
    flush(output_unit, iostat=ios)   ! (so the printed log isn't interleaved with Python's output)

    if (allocated(r%x) .and. size(x) == n) x = r%x
    if (allocated(r%lambda) .and. size(lambda) == m) lambda = r%lambda
    if (allocated(r%z) .and. size(z) == n) z = r%z
    if (allocated(r%c) .and. size(c) == m) c = r%c
    iinfo(1:sqpopt_python_n_iinfo) = [r%istat, r%iterations, r%n_eval_fc, r%n_eval_gjac, r%n_eval_hess, r%n_qp_iterations, &
                   r%derivative_switch_iteration, r%n_soc, r%n_hessian_resets, r%n_restoration_steps, &
                   r%n_restoration_phases, r%n_elastic, r%n_escape, r%n_factorizations, r%n_qp_solves, r%n_direct_qp, &
                   r%n_unconstrained_qp]
    rinfo(1:sqpopt_python_n_rinfo) = [r%f, r%kkt_error, r%feasibility_error, r%stationarity_error, r%time, &
                                      r%time_functions, r%time_qp, r%time_factorization]
    message = ''
    if (allocated(r%message)) message = r%message

    call restore()

    contains

        subroutine restore()
        !! restore the callbacks of an enclosing solve, and close `output_file`
        if (unit /= -1) close(unit, iostat=ios)
        py_fc => fc0
        py_gjac => gjac0
        py_hess => hess0
        py_report => report0
        end subroutine restore

    end subroutine sqpopt_python_solve
!*******************************************************************************

!*******************************************************************************
!>
!  `fc` for the library: calls the current solve's Python `fc`

    subroutine fc_bridge(x, f, c, status, data)
    real(real64), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(real64),               intent(out)   :: f      !! objective value at `x`
    real(real64), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                    intent(inout) :: status !! `0` on entry; `> 0`: can't evaluate at `x`; `< 0`: stop
    class(*), optional,         intent(inout) :: data   !! (not used)
    real(real64), dimension(1) :: f1
    f1 = 0.0_real64
    c = 0.0_real64
    call py_fc(size(x), size(c), status, x, f1, c)
    f = f1(1)
    end subroutine fc_bridge

!*******************************************************************************
!>
!  `gjac` for the library: calls the current solve's Python `gjac`

    subroutine gjac_bridge(x, g, jac_val, accuracy, status, data)
    real(real64), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(real64), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(real64), dimension(:), intent(out)   :: jac_val  !! Jacobian values at `x`, in the pattern's order
    integer,                    intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                    intent(inout) :: status   !! `0` on entry; `> 0`: can't evaluate at `x`; `< 0`: stop
    class(*), optional,         intent(inout) :: data     !! (not used)
    g = 0.0_real64
    jac_val = 0.0_real64
    call py_gjac(size(x), size(jac_val), status, accuracy, x, g, jac_val)
    end subroutine gjac_bridge

!*******************************************************************************
!>
!  `hess` for the library: calls the current solve's Python `hess`

    subroutine hess_bridge(x, lambda, hess_val, status, data)
    real(real64), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(real64), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(real64), dimension(:), intent(out)   :: hess_val !! Hessian values at `x`, in the pattern's order
    integer,                    intent(inout) :: status   !! `0` on entry; `> 0`: can't evaluate at `x`; `< 0`: stop
    class(*), optional,         intent(inout) :: data     !! (not used)
    hess_val = 0.0_real64
    call py_hess(size(x), size(lambda), size(hess_val), status, x, lambda, hess_val)
    end subroutine hess_bridge

!*******************************************************************************
!>
!  `report` for the library: calls the current solve's Python `report`

    subroutine report_bridge(iter, x, f, c, lambda, user_stop, data)
    integer,                    intent(in)    :: iter      !! major iteration number
    real(real64), dimension(:), intent(in)    :: x         !! point `dimension(n)`
    real(real64),               intent(in)    :: f         !! objective value at `x`
    real(real64), dimension(:), intent(in)    :: c         !! constraint values at `x` `dimension(m)`
    real(real64), dimension(:), intent(in)    :: lambda    !! constraint multipliers `dimension(m)`
    logical,                    intent(out)   :: user_stop !! set to `.true.` to stop the solver
    class(*), optional,         intent(inout) :: data      !! (not used)
    integer(int32) :: stop
    stop = 0
    call py_report(size(x), size(c), stop, iter, x, f, c, lambda)
    user_stop = stop /= 0
    end subroutine report_bridge
!*******************************************************************************

    end module sqpopt_python
!*******************************************************************************
