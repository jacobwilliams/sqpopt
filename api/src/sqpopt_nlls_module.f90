!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Nonlinear least-squares problems, solved with the SQP solver
!  ([[sqpopt_nlls_type]]): an optional layer on top of [[sqpopt_type]],
!  which the rest of the library doesn't know about.
!
!  The problem is
!
!  $$ \min_x \; \tfrac12 \sum_{i=1}^{l} r_i(x)^2 \quad \text{s.t.} \quad
!     c_l \le c(x) \le c_u, \quad x_l \le x \le x_u $$
!
!  for residuals \( r_i \) (e.g. the differences between a model and
!  data) and, optionally, constraints \( c \). It is solved in the form
!  proposed by Schittkowski (the code DFNLP): with a new variable
!  \( z_i \) for each residual,
!
!  $$ \min_{x,z} \; \tfrac12 z^T z \quad \text{s.t.} \quad r(x) - z = 0,
!     \quad c_l \le c(x) \le c_u, \quad x_l \le x \le x_u $$
!
!  which is an ordinary nonlinear program for [[sqpopt_type]], in
!  \( n + l \) variables and \( l + m \) constraints. Its QP subproblems
!  are always consistent as far as the residuals go (any step in \( x \)
!  has a \( z \) to go with it), and their step is a Gauss-Newton step
!  with a correction from the Hessian approximation, as in the special
!  codes for least-squares problems.
!
!  **Use it** for data fitting, and for a system of equations that can't
!  all hold (more equations than unknowns): given as equality constraints
!  to [[sqpopt_type]] they make every QP subproblem inconsistent, and the
!  solve ends as `sqpopt_infeasible` at best, after many restoration
!  steps.
!
!  **The Hessian.** With the quasi-Newton modes of `options%hessian_mode`
!  (the default), the solver approximates the Hessian of the transformed
!  problem as it does any other. With `sqpopt_hessian_exact`, this module
!  supplies the *Gauss-Newton* Hessian: the identity on \( z \), and zero
!  on \( x \) (the second derivatives of the residuals and of the
!  constraints are neglected, as in the Gauss-Newton method). The solver's
!  Hessian shift then acts as the damping of the Levenberg-Marquardt
!  method. It usually takes fewer iterations, mostly so when the residuals
!  are small at the solution, but it fails more often.
!
!  **What to expect.** On 94 overdetermined systems of equations from
!  CUTEst (at most 250 iterations), the default L-BFGS reached the best
!  sum of squares found on 59, with 6,842 residual evaluations in all, and
!  the Gauss-Newton Hessian on 68, with 13,154 (one or the other: 79).
!  Given to [[sqpopt_type]] as equality constraints: 52, with 52,958. A
!  code written for least squares does better on such problems (scipy's
!  `least_squares`: 87, with 10,326), so prefer one (e.g. MINPACK) for an
!  unconstrained fit that is hard. What this interface adds is the
!  solver's constraints and sparsity. It is a layer on the solver, not a
!  least-squares method of its own: the solver is unchanged by it.

    module sqpopt_nlls_module

    use sqpopt_kinds,               only: wp => sqpopt_module_wp
    use sqpopt_types_module,        only: sqpopt_results_type, sqpopt_infinity, sqpopt_all_finite
    use sqpopt_module,              only: sqpopt_type
    use sqpopt_problem_module,      only: sqpopt_problem_type
    use sqpopt_options_module,      only: sqpopt_options_type
    use sqpopt_hessian_module,      only: sqpopt_hessian_type, sqpopt_hessian_exact
    use sqpopt_qp_solver_module,    only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module,   only: sqpopt_linesearch_type
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type

    implicit none

    private

    abstract interface

        subroutine sqpopt_residual_func(x, r, c, status, data)
            !! evaluates the residuals \( r(x) \) and the constraints \( c(x) \)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)    :: x      !! variables `dimension(n)`
            real(wp), dimension(:), intent(out)   :: r      !! residuals `dimension(l)`
            real(wp), dimension(:), intent(out)   :: c      !! constraints `dimension(m)` (`m` may be 0)
            integer,                intent(inout) :: status !! `0` on entry; `>0`: can't evaluate here, `<0`: stop
            class(*), optional,     intent(inout) :: data   !! user data (see `sqpopt_nlls_type%initialize`)
        end subroutine sqpopt_residual_func

        subroutine sqpopt_residual_jac_func(x, rjac_val, cjac_val, accuracy, status, data)
            !! evaluates the nonzero values of the Jacobians of the residuals
            !! and of the constraints, each in the order of its sparsity pattern
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)    :: x        !! variables `dimension(n)`
            real(wp), dimension(:), intent(out)   :: rjac_val !! nonzeros of the residuals' Jacobian
            real(wp), dimension(:), intent(out)   :: cjac_val !! nonzeros of the constraints' Jacobian (may be none)
            integer,                intent(in)    :: accuracy !! `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
            integer,                intent(inout) :: status   !! `0` on entry; `>0`: can't evaluate here, `<0`: stop
            class(*), optional,     intent(inout) :: data     !! user data (see `sqpopt_nlls_type%initialize`)
        end subroutine sqpopt_residual_jac_func

    end interface

    public :: sqpopt_residual_func, sqpopt_residual_jac_func

    type :: nlls_state
        !! what the solver's callbacks need: the user's functions and the
        !! sizes. The solver is given a pointer to it as its user data.
        integer :: n = 0     !! number of variables `x`
        integer :: l = 0     !! number of residuals
        integer :: m = 0     !! number of constraints
        integer :: nnz_r = 0 !! nonzeros of the residuals' Jacobian
        procedure(sqpopt_residual_func),     pointer, nopass :: residuals => null() !! the user's residuals and constraints
        procedure(sqpopt_residual_jac_func), pointer, nopass :: jacobian  => null() !! the user's Jacobians
        class(*), pointer :: user_data => null() !! the user's data (if any)
        ! the residuals and constraints at the starting point, which `solve`
        ! evaluates to start `z` from: kept for the solver's first evaluation there
        logical :: have_start = .false.                  !! whether they are kept
        real(wp), dimension(:), allocatable :: x_start   !! the point `dimension(n)`
        real(wp), dimension(:), allocatable :: rc_start  !! `r` then `c` there `dimension(l+m)`
    end type nlls_state

    type, public :: sqpopt_nlls_type
        !! a nonlinear least-squares problem and its solver (see the module
        !! documentation). `initialize` defines the problem, `solve` solves it
        !! from a starting point, and `get_solution` returns the variables
        !! and the residuals there. Call `destroy` when done with it.

        private

        type(sqpopt_type) :: solver                   !! the SQP solver, for the transformed problem
        type(nlls_state), pointer :: state => null()  !! what its callbacks need (see `nlls_state`)
        real(wp), dimension(:), allocatable :: x_lb   !! lower bounds of the variables `dimension(n)`
        real(wp), dimension(:), allocatable :: x_ub   !! upper bounds of the variables `dimension(n)`
        real(wp), dimension(:), allocatable :: c_lb   !! lower bounds of the constraints `dimension(m)`
        real(wp), dimension(:), allocatable :: c_ub   !! upper bounds of the constraints `dimension(m)`
        real(wp) :: ctol = 0.0_wp                     !! the feasibility tolerance of the options

        contains

        procedure, public :: initialize     => nlls_initialize
        procedure, public :: solve          => nlls_solve
        procedure, public :: get_solution   => nlls_get_solution
        procedure, public :: get_results    => nlls_get_results
        procedure, public :: status_message => nlls_status_message
        procedure, public :: destroy        => nlls_destroy

    end type sqpopt_nlls_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  define the least-squares problem: `n` variables, `n_residuals` residuals
!  with the sparsity pattern `rjac_irow`/`rjac_icol` of their Jacobian
!  (1-based row and column of each nonzero, as for
!  [[set_jacobian_sparsity]]), and optionally bounds on the variables and
!  constraints `c_lb <= c(x) <= c_ub` with the pattern
!  `cjac_irow`/`cjac_icol` of their Jacobian (give all four, or none).
!  The solver's components are configured as for [[sqpopt_initialize]].
!  With `options%hessian_mode = sqpopt_hessian_exact`, the Gauss-Newton
!  Hessian is used (see the module documentation).
!
!  Invalid input (sizes that don't match, indices out of range) is
!  reported by `solve`, as `sqpopt_invalid_input`.

    subroutine nlls_initialize(me, n, n_residuals, residuals, jacobian, rjac_irow, rjac_icol, x_lb, x_ub, &
                               c_lb, c_ub, cjac_irow, cjac_icol, options, hessian, qp_solver, linesearch, &
                               trust_region, data)

    class(sqpopt_nlls_type), intent(inout) :: me
    integer,                 intent(in)    :: n           !! number of variables
    integer,                 intent(in)    :: n_residuals !! number of residuals \( l \)
    procedure(sqpopt_residual_func)        :: residuals   !! evaluates the residuals and the constraints
    procedure(sqpopt_residual_jac_func)    :: jacobian    !! evaluates their Jacobians' nonzeros
    integer,  dimension(:),  intent(in)    :: rjac_irow   !! residuals' Jacobian pattern: row (residual) indices
    integer,  dimension(:),  intent(in)    :: rjac_icol   !! residuals' Jacobian pattern: column (variable) indices
    real(wp), dimension(:),  intent(in), optional :: x_lb !! lower bounds of the variables `dimension(n)` (default: none)
    real(wp), dimension(:),  intent(in), optional :: x_ub !! upper bounds of the variables `dimension(n)` (default: none)
    real(wp), dimension(:),  intent(in), optional :: c_lb !! lower bounds of the constraints `dimension(m)`
    real(wp), dimension(:),  intent(in), optional :: c_ub !! upper bounds of the constraints `dimension(m)`
    integer,  dimension(:),  intent(in), optional :: cjac_irow !! constraints' Jacobian pattern: row indices
    integer,  dimension(:),  intent(in), optional :: cjac_icol !! constraints' Jacobian pattern: column indices
    type(sqpopt_options_type),      intent(in), optional :: options      !! solver options
    type(sqpopt_hessian_type),      intent(in), optional :: hessian      !! Hessian approximation settings
    type(sqpopt_qp_solver_type),    intent(in), optional :: qp_solver    !! QP solver settings
    type(sqpopt_linesearch_type),   intent(in), optional :: linesearch   !! line search settings
    type(sqpopt_trust_region_type), intent(in), optional :: trust_region !! trust-region settings
    class(*), target, intent(inout), optional :: data !! user data, passed to `residuals` and `jacobian`

    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: opts
    integer,  dimension(:), allocatable :: irow, icol, hrow
    real(wp), dimension(:), allocatable :: lb, ub, clb, cub
    integer :: l, m, nnz_r, nnz_c, k

    call me%destroy()

    l = max(n_residuals, 0)
    m = 0
    nnz_c = 0
    if (present(c_lb)) m = size(c_lb)
    if (present(cjac_irow)) nnz_c = size(cjac_irow)
    nnz_r = size(rjac_irow)

    allocate(me%state)
    me%state%n = n
    me%state%l = l
    me%state%m = m
    me%state%nnz_r = nnz_r
    me%state%residuals => residuals
    me%state%jacobian  => jacobian
    if (present(data)) me%state%user_data => data

    ! the bounds: those of `x`, and none on `z`; the residual rows are equalities
    allocate(me%x_lb(max(n, 0)), me%x_ub(max(n, 0)), me%c_lb(m), me%c_ub(m))
    me%x_lb = -sqpopt_infinity
    me%x_ub =  sqpopt_infinity
    if (present(x_lb)) then
        if (size(x_lb) == n) me%x_lb = x_lb
    end if
    if (present(x_ub)) then
        if (size(x_ub) == n) me%x_ub = x_ub
    end if
    if (present(c_lb)) me%c_lb = c_lb
    me%c_ub = sqpopt_infinity
    if (present(c_ub)) then
        if (size(c_ub) == m) me%c_ub = c_ub
    end if
    allocate(lb(max(n, 0) + l), ub(max(n, 0) + l), clb(l + m), cub(l + m))
    lb = -sqpopt_infinity
    ub =  sqpopt_infinity
    lb(1:max(n, 0)) = me%x_lb
    ub(1:max(n, 0)) = me%x_ub
    clb = 0.0_wp
    cub = 0.0_wp
    clb(l+1:) = me%c_lb
    cub(l+1:) = me%c_ub

    ! the Jacobian of the transformed constraints: the residuals' Jacobian,
    ! `-1` for each `z_i` in its row, then the constraints' Jacobian
    allocate(irow(nnz_r + l + nnz_c), icol(nnz_r + l + nnz_c))
    irow(1:nnz_r) = rjac_irow
    icol(1:nnz_r) = rjac_icol(1:min(nnz_r, size(rjac_icol)))
    do k = 1, l
        irow(nnz_r + k) = k
        icol(nnz_r + k) = n + k
    end do
    if (nnz_c > 0) then
        irow(nnz_r+l+1:) = l + cjac_irow
        icol(nnz_r+l+1:) = cjac_icol(1:min(nnz_c, size(cjac_icol)))
    end if

    if (present(options)) then
        opts = options
    else
        opts = sqpopt_options_type()
    end if
    me%ctol = opts%ctol

    call problem%set_problem_size(n=max(n, 0) + l, m=l + m)
    call problem%set_bounds(lb, ub, clb, cub)
    call problem%set_jacobian_sparsity(nnz_r + l + nnz_c, irow, icol)
    if (opts%hessian_mode == sqpopt_hessian_exact) then
        ! (the Gauss-Newton Hessian: the identity on `z`)
        allocate(hrow(l))
        do k = 1, l
            hrow(k) = n + k
        end do
        call problem%set_functions(fc=nlls_fc, gjac=nlls_gjac, hess=nlls_hess, data=me%state)
        call problem%set_hessian_sparsity(l, hrow, hrow)
    else
        call problem%set_functions(fc=nlls_fc, gjac=nlls_gjac, data=me%state)
    end if

    call me%solver%initialize(problem=problem, options=opts, hessian=hessian, qp_solver=qp_solver, &
                              linesearch=linesearch, trust_region=trust_region)

    end subroutine nlls_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  solve the problem from the starting point `x0` (moved inside the bounds,
!  as [[sqpopt_solve]] does). `istat` is the solver's status code (see
!  [[sqpopt_types_module]]): `sqpopt_success` at a least-squares solution.
!
!  The residuals are evaluated at the starting point, to start \( z \)
!  from them if the constraints hold there (and from zero otherwise, as
!  Schittkowski recommends); the solver's first evaluation reuses them.
!  A second `solve` starts from the configuration given to `initialize`, as
!  for [[sqpopt_type]].

    subroutine nlls_solve(me, x0, istat)

    class(sqpopt_nlls_type), intent(inout) :: me
    real(wp), dimension(:),  intent(in)    :: x0    !! starting point `dimension(n)`
    integer,                 intent(out)   :: istat !! status code (see [[sqpopt_types_module]])

    real(wp), dimension(:), allocatable :: y0
    integer :: n, l, m, status
    logical :: feasible

    if (.not. associated(me%state)) then
        ! (not initialized: let the solver report the invalid problem)
        call me%solver%solve(x0, istat)
        return
    end if
    n = me%state%n
    l = me%state%l
    m = me%state%m

    allocate(y0(max(n, 0) + l))
    y0 = 0.0_wp
    me%state%have_start = .false.
    if (size(x0) == n .and. n > 0) then
        y0(1:n) = min(max(x0, me%x_lb), me%x_ub)
        if (sqpopt_all_finite(x0)) then
            if (allocated(me%state%x_start)) deallocate(me%state%x_start, me%state%rc_start)
            allocate(me%state%x_start(n), me%state%rc_start(l + m))
            me%state%x_start = y0(1:n)
            status = 0
            if (associated(me%state%user_data)) then
                call me%state%residuals(me%state%x_start, me%state%rc_start(1:l), me%state%rc_start(l+1:), status, &
                                        me%state%user_data)
            else
                call me%state%residuals(me%state%x_start, me%state%rc_start(1:l), me%state%rc_start(l+1:), status)
            end if
            if (status == 0 .and. sqpopt_all_finite(me%state%rc_start)) then
                me%state%have_start = .true.
                feasible = .true.
                if (m > 0) feasible = all(me%state%rc_start(l+1:) >= me%c_lb - me%ctol) .and. &
                                      all(me%state%rc_start(l+1:) <= me%c_ub + me%ctol)
                if (feasible) y0(n+1:) = me%state%rc_start(1:l)
            end if
        end if
        call me%solver%solve(y0, istat)
    else
        ! (a starting point of the wrong size: reported by the solver)
        call me%solver%solve(x0, istat)
    end if

    end subroutine nlls_solve
!*******************************************************************************

!*******************************************************************************
!>
!  the solution of the last `solve`: the variables, and optionally the
!  residuals there, the sum of their squares, and the multipliers of the
!  constraints (with the sign convention of [[sqpopt_results_type]]).

    subroutine nlls_get_solution(me, x, residuals, sum_of_squares, lambda)

    class(sqpopt_nlls_type), intent(in)  :: me
    real(wp), dimension(:),  intent(out) :: x                        !! variables `dimension(n)`
    real(wp), dimension(:),  intent(out), optional :: residuals      !! residuals \( r(x) \) `dimension(l)`
    real(wp),                intent(out), optional :: sum_of_squares !! \( \sum r_i(x)^2 \)
    real(wp), dimension(:),  intent(out), optional :: lambda         !! constraint multipliers `dimension(m)`

    type(sqpopt_results_type) :: r
    integer :: n, l, i
    real(wp) :: ri

    x = 0.0_wp
    if (present(residuals)) residuals = 0.0_wp
    if (present(sum_of_squares)) sum_of_squares = 0.0_wp
    if (present(lambda)) lambda = 0.0_wp
    if (.not. associated(me%state)) return
    call me%solver%get_results(r)
    n = me%state%n
    l = me%state%l
    if (.not. allocated(r%x)) return
    if (size(r%x) /= n + l .or. size(r%c) /= l + me%state%m) return

    x = r%x(1:n)
    ! (the residual is `z` plus the value of its constraint `r - z`, which is zero at a feasible point)
    do i = 1, l
        ri = r%x(n+i) + r%c(i)
        if (present(residuals)) residuals(i) = ri
        if (present(sum_of_squares)) sum_of_squares = sum_of_squares + ri**2
    end do
    if (present(lambda)) lambda = r%lambda(l+1:)

    end subroutine nlls_get_solution
!*******************************************************************************

!*******************************************************************************
!>
!  the results of the last `solve`, as the SQP solver returns them for the
!  *transformed* problem (see the module documentation): `x` is the
!  variables followed by \( z \), `f` is \( \tfrac12 z^T z \), and `c` and
!  `lambda` are those of the residual equations followed by the
!  constraints. The status, the counts, and the times are as for any solve
!  (`n_eval_fc` and `n_eval_gjac` count the calls of `residuals` and
!  `jacobian`).

    subroutine nlls_get_results(me, results)

    class(sqpopt_nlls_type),   intent(in)  :: me
    type(sqpopt_results_type), intent(out) :: results !! the results of the last `solve`

    call me%solver%get_results(results)

    end subroutine nlls_get_results
!*******************************************************************************

!*******************************************************************************
!>
!  a description of the status of the last `solve`.

    function nlls_status_message(me) result(msg)

    class(sqpopt_nlls_type), intent(in) :: me
    character(len=:), allocatable :: msg

    msg = me%solver%status_message()

    end function nlls_status_message
!*******************************************************************************

!*******************************************************************************
!>
!  free everything.

    subroutine nlls_destroy(me)

    class(sqpopt_nlls_type), intent(inout) :: me

    call me%solver%destroy()
    if (associated(me%state)) deallocate(me%state)
    me%state => null()
    if (allocated(me%x_lb)) deallocate(me%x_lb)
    if (allocated(me%x_ub)) deallocate(me%x_ub)
    if (allocated(me%c_lb)) deallocate(me%c_lb)
    if (allocated(me%c_ub)) deallocate(me%c_ub)

    end subroutine nlls_destroy
!*******************************************************************************

!*******************************************************************************
!>
!  the objective \( \tfrac12 z^T z \) and the constraints \( r(x) - z \)
!  and \( c(x) \) of the transformed problem, for the solver. `data` is the
!  wrapper's state (`nlls_state`).

    subroutine nlls_fc(x, f, c, status, data)

    real(wp), dimension(:), intent(in)    :: x      !! the variables, then `z` `dimension(n+l)`
    real(wp),               intent(out)   :: f      !! \( \tfrac12 z^T z \)
    real(wp), dimension(:), intent(out)   :: c      !! \( r(x) - z \), then \( c(x) \) `dimension(l+m)`
    integer,                intent(inout) :: status !! `0` on entry; the user function's status on exit
    class(*), optional,     intent(inout) :: data   !! the wrapper's state (`nlls_state`)

    integer :: n, l, i

    f = 0.0_wp
    c = 0.0_wp
    if (.not. present(data)) return
    select type (data)
    type is (nlls_state)
        n = data%n
        l = data%l
        if (data%have_start) then
            ! (the starting point, which `solve` has evaluated already)
            data%have_start = .false.
            if (all(x(1:n) == data%x_start)) then
                c = data%rc_start
            else
                call evaluate()
            end if
        else
            call evaluate()
        end if
        do i = 1, l
            c(i) = c(i) - x(n+i)
            f = f + x(n+i)**2
        end do
        f = 0.5_wp*f
    end select

    contains

        subroutine evaluate()
        !! call the user's function: the residuals and the constraints, straight into `c`
        select type (data)
        type is (nlls_state)
            if (associated(data%user_data)) then
                call data%residuals(x(1:n), c(1:l), c(l+1:), status, data%user_data)
            else
                call data%residuals(x(1:n), c(1:l), c(l+1:), status)
            end if
        end select
        end subroutine evaluate

    end subroutine nlls_fc
!*******************************************************************************

!*******************************************************************************
!>
!  the objective gradient \( (0, z) \) and the Jacobian's nonzeros (the
!  residuals' Jacobian, `-1` for each \( z_i \), then the constraints'
!  Jacobian) of the transformed problem, for the solver.

    subroutine nlls_gjac(x, g, jac_val, accuracy, status, data)

    real(wp), dimension(:), intent(in)    :: x        !! the variables, then `z` `dimension(n+l)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient `dimension(n+l)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzeros of the Jacobian (see above)
    integer,                intent(in)    :: accuracy !! requested accuracy, passed on to the user's function
    integer,                intent(inout) :: status   !! `0` on entry; the user function's status on exit
    class(*), optional,     intent(inout) :: data     !! the wrapper's state (`nlls_state`)

    integer :: n, l, nnz_r

    g = 0.0_wp
    jac_val = 0.0_wp
    if (.not. present(data)) return
    select type (data)
    type is (nlls_state)
        n = data%n
        l = data%l
        nnz_r = data%nnz_r
        g(n+1:n+l) = x(n+1:n+l)
        if (associated(data%user_data)) then
            call data%jacobian(x(1:n), jac_val(1:nnz_r), jac_val(nnz_r+l+1:), accuracy, status, data%user_data)
        else
            call data%jacobian(x(1:n), jac_val(1:nnz_r), jac_val(nnz_r+l+1:), accuracy, status)
        end if
        jac_val(nnz_r+1:nnz_r+l) = -1.0_wp
    end select

    end subroutine nlls_gjac
!*******************************************************************************

!*******************************************************************************
!>
!  the Gauss-Newton Hessian of the transformed problem, for the solver (with
!  `options%hessian_mode = sqpopt_hessian_exact`): the identity on \( z \).
!  Its pattern is the diagonal elements of \( z \).

    subroutine nlls_hess(x, lambda, hess_val, status, data)

    real(wp), dimension(:), intent(in)    :: x        !! the variables, then `z` `dimension(n+l)` (unused)
    real(wp), dimension(:), intent(in)    :: lambda   !! multipliers `dimension(l+m)` (unused: the residuals'
                                                      !! second derivatives are neglected)
    real(wp), dimension(:), intent(out)   :: hess_val !! the Hessian's nonzeros `dimension(l)`
    integer,                intent(inout) :: status   !! `0` on entry (left unchanged)
    class(*), optional,     intent(inout) :: data     !! the wrapper's state (`nlls_state`) (unused)

    hess_val = 1.0_wp

    end subroutine nlls_hess
!*******************************************************************************

    end module sqpopt_nlls_module
!*******************************************************************************
