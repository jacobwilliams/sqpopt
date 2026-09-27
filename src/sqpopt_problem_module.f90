!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Defines the nonlinear program (NLP) that is to be solved:
!
!  $$ \min_{x} f(x) $$
!
!  subject to:
!
!  $$ c_l \le c(x) \le c_u $$
!  $$ x_l \le x \le x_u $$
!
!  where \( c(x) \) is a vector of nonlinear constraints. An equality
!  constraint is specified by setting the corresponding elements of
!  \( c_l \) and \( c_u \) to the same value, and an inequality
!  constraint by setting them to different values (using \( \pm\infty \)
!  for one-sided inequalities). Variable bounds are specified similarly
!  via \( x_l \) and \( x_u \).
!
!  The Jacobian of \( c(x) \) and the Hessian of the Lagrangian are
!  never treated as dense \( m \times n \) or \( n \times n \) arrays:
!  their (fixed) sparsity patterns are supplied once as 1-based COO
!  `irow`/`icol` triplets (the same convention used by the `lusol` and
!  `LSQR` dependencies), and the user-supplied evaluation
!  routines only need to fill in the corresponding nonzero *values* on
!  each call.
!
!  **User functions.** The problem is defined by two routines (see
!  [[set_functions]]): `fc`, which evaluates the objective \( f(x) \) and
!  the constraints \( c(x) \) together, and `gjac`, which evaluates the
!  objective gradient \( \nabla f(x) \) and the nonzero values of the
!  constraint Jacobian together (and, optionally, `hess` for an exact
!  Hessian of the Lagrangian). Every user function has two trailing
!  arguments:
!
!  * `status` (`integer, intent(inout)`): `0` on entry. Leave it `0` on
!    success; set it `> 0` if the function cannot be evaluated at `x` (e.g.
!    a domain error) -- the solver then treats the point like one where
!    the function returned NaN, and backs off from it -- or `< 0` to ask the
!    solver to stop (`istat=sqpopt_user_requested_stop`).
!  * `data` (`class(*), intent(inout), optional`): the user data object
!    given to [[set_functions]] (absent if none was given), for passing any
!    context to the functions without module variables; use `select type`
!    to access it.
!
!  **Evaluation layer.** The solver evaluates everything through the
!  type-bound [[eval_f_cached|f]], [[eval_c_cached|c]], [[eval_g_cached|g]],
!  and [[eval_jac_cached|jac]] methods, which handle `status`/`data`, count
!  the user calls, keep a small cache of recent evaluations (so e.g. the
!  point accepted by a line search is not evaluated again, and asking for
!  `c` at a point where `f` was just evaluated doesn't call `fc` again),
!  and apply the objective and constraint scaling (see
!  [[compute_scaling]]): the solver works with \( s_f f \) and
!  \( s_{c,i} c_i \).

    module sqpopt_problem_module

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_success, sqpopt_invalid_input, sqpopt_infinity, sqpopt_all_finite
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    implicit none

    private

    integer, parameter :: cache_size = 4 !! number of recent evaluations of `fc` kept

    public :: sqpopt_fc_func, sqpopt_gjac_func, sqpopt_hessian_func

    type, public :: sqpopt_problem_type
        !! defines the problem to be solved: the problem size, the
        !! variable and constraint bounds, the sparsity patterns of the
        !! constraint Jacobian and Lagrangian Hessian, and the
        !! user-supplied procedures used to evaluate the objective
        !! function, constraints, and their derivatives.

        integer :: n      = 0  !! number of optimization variables (\( n>0 \))
        integer :: m      = 0  !! total number of nonlinear constraints \( c(x) \) (\( m \ge 0 \))
        integer :: m_eq   = 0  !! number of nonlinear equality constraints (\( c_l = c_u \))
        integer :: m_ineq = 0  !! number of nonlinear inequality constraints (\( c_l \ne c_u \))

        real(wp), dimension(:), allocatable :: x_lb  !! lower bounds on the optimization variables `dimension(n)`
        real(wp), dimension(:), allocatable :: x_ub  !! upper bounds on the optimization variables `dimension(n)`

        real(wp), dimension(:), allocatable :: c_lb  !! lower bounds on the constraints `dimension(m)`
        real(wp), dimension(:), allocatable :: c_ub  !! upper bounds on the constraints `dimension(m)`

        integer :: jac_nnz = 0   !! number of nonzero elements in the constraint Jacobian
        integer, dimension(:), allocatable :: jac_irow  !! Jacobian sparsity pattern: row indices `dimension(jac_nnz)`
        integer, dimension(:), allocatable :: jac_icol  !! Jacobian sparsity pattern: column indices `dimension(jac_nnz)`

        integer :: hess_nnz = 0  !! number of nonzero elements in the (symmetric) Hessian of the Lagrangian
        integer, dimension(:), allocatable :: hess_irow !! Hessian sparsity pattern: row indices `dimension(hess_nnz)`
        integer, dimension(:), allocatable :: hess_icol !! Hessian sparsity pattern: column indices `dimension(hess_nnz)`

        procedure(sqpopt_fc_func),        pointer, nopass :: eval_fc   => null() !! evaluates \( f(x) \) and \( c(x) \)
        procedure(sqpopt_gjac_func),      pointer, nopass :: eval_gjac => null() !! evaluates \( \nabla f(x) \) and the
                                                                                !! nonzero values of the Jacobian of \( c(x) \)
        procedure(sqpopt_hessian_func),   pointer, nopass :: eval_hess => null() !! evaluates the nonzero values of the
                                                                                !! Hessian of the Lagrangian (only used
                                                                                !! when an exact Hessian is requested)
        class(*), pointer :: user_data => null() !! user data passed to every user function (see [[set_functions]])

        ! ---- internal evaluation state (see the module-level documentation) ----
        logical  :: stop_requested = .false. !! set when a user function returns `status < 0`
        integer  :: n_eval_fc   = 0 !! number of calls of the user's `fc`
        integer  :: n_eval_gjac = 0 !! number of calls of the user's `gjac`
        real(wp) :: f_scale = 1.0_wp !! objective scale factor \( s_f \)
        real(wp), dimension(:), allocatable :: c_scale !! constraint scale factors \( s_{c,i} \) `dimension(m)`
        integer :: cache_n = 0, cache_next = 1 !! entries used, and the slot for the next one, in the `fc` cache
        real(wp), dimension(:,:), allocatable :: cache_x !! points at which `fc` was evaluated `dimension(n,cache_size)`
        real(wp), dimension(:),   allocatable :: cache_f !! (unscaled) `f` at those points `dimension(cache_size)`
        real(wp), dimension(:,:), allocatable :: cache_c !! (unscaled) `c` at those points `dimension(m,cache_size)`
        logical :: have_gjac = .false. !! whether the one-entry `gjac` cache is filled
        real(wp), dimension(:), allocatable :: cache_xg  !! point at which `gjac` was evaluated `dimension(n)`
        real(wp), dimension(:), allocatable :: cache_g   !! (unscaled) `g` there `dimension(n)`
        real(wp), dimension(:), allocatable :: cache_jac !! (unscaled) Jacobian values there `dimension(jac_nnz)`

        contains

        procedure, public :: set_problem_size      !! set the problem dimensions and allocate the bound arrays
        procedure, public :: set_bounds            !! set the variable and constraint bounds
        procedure, public :: set_jacobian_sparsity !! set the (fixed) sparsity pattern of the constraint Jacobian
        procedure, public :: set_hessian_sparsity  !! set the (fixed) sparsity pattern of the Lagrangian Hessian
        procedure, public :: set_functions         !! attach the user-supplied evaluation procedures (and user data)
        procedure, public :: validate              !! check the problem definition and normalize infinite bounds
        procedure, public :: f   => eval_f_cached   !! evaluate the (scaled) objective
        procedure, public :: c   => eval_c_cached   !! evaluate the (scaled) constraints
        procedure, public :: g   => eval_g_cached   !! evaluate the (scaled) objective gradient
        procedure, public :: jac => eval_jac_cached !! evaluate the (scaled) Jacobian values
        procedure, public :: reset_evaluations     !! empty the caches, zero the counters, and remove any scaling
        procedure, public :: compute_scaling       !! set gradient-based objective/constraint scale factors

    end type sqpopt_problem_type

    abstract interface

        subroutine sqpopt_fc_func(x, f, c, status, data)
            !! evaluates the objective function \( f(x) \) and the nonlinear
            !! constraint vector \( c(x) \)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)    :: x      !! optimization variable vector `dimension(n)`
            real(wp),               intent(out)   :: f      !! value of the objective function
            real(wp), dimension(:), intent(out)   :: c      !! constraint vector `dimension(m)` (`m` may be 0)
            integer,                intent(inout) :: status !! `0` on entry; `>0`: can't evaluate here, `<0`: stop
            class(*), optional,     intent(inout) :: data   !! user data (see [[set_functions]])
        end subroutine sqpopt_fc_func

        subroutine sqpopt_gjac_func(x, g, jac_val, status, data)
            !! evaluates the gradient of the objective function \( \nabla f(x) \)
            !! and the nonzero values of the Jacobian of the constraint vector,
            !! \( J_{ij} = \partial c_i / \partial x_j \), ordered to match the
            !! sparsity pattern set by `set_jacobian_sparsity`
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)    :: x       !! optimization variable vector `dimension(n)`
            real(wp), dimension(:), intent(out)   :: g       !! gradient vector `dimension(n)`
            real(wp), dimension(:), intent(out)   :: jac_val !! nonzero Jacobian values `dimension(jac_nnz)` (may be 0)
            integer,                intent(inout) :: status  !! `0` on entry; `>0`: can't evaluate here, `<0`: stop
            class(*), optional,     intent(inout) :: data    !! user data (see [[set_functions]])
        end subroutine sqpopt_gjac_func

        subroutine sqpopt_hessian_func(x, lambda, hess_val, status, data)
            !! evaluates the nonzero values of the Hessian of the Lagrangian:
            !! \( H = \nabla^2 f(x) - \sum_i \lambda_i \nabla^2 c_i(x) \),
            !! ordered to match the sparsity pattern set by `set_hessian_sparsity`.
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)    :: x        !! optimization variable vector `dimension(n)`
            real(wp), dimension(:), intent(in)    :: lambda   !! Lagrange multipliers `dimension(m)`
            real(wp), dimension(:), intent(out)   :: hess_val !! nonzero Hessian values `dimension(hess_nnz)`
            integer,                intent(inout) :: status   !! `0` on entry; `>0`: can't evaluate here, `<0`: stop
            class(*), optional,     intent(inout) :: data     !! user data (see [[set_functions]])
        end subroutine sqpopt_hessian_func

    end interface

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  set the problem dimensions and allocate the variable/constraint bound arrays.

    subroutine set_problem_size(me, n, m_eq, m_ineq)

    class(sqpopt_problem_type), intent(inout) :: me
    integer, intent(in) :: n       !! number of optimization variables
    integer, intent(in) :: m_eq    !! number of nonlinear equality constraints
    integer, intent(in) :: m_ineq  !! number of nonlinear inequality constraints

    me%n      = n
    me%m_eq   = m_eq
    me%m_ineq = m_ineq
    me%m      = m_eq + m_ineq

    if (allocated(me%x_lb)) deallocate(me%x_lb)
    if (allocated(me%x_ub)) deallocate(me%x_ub)
    if (allocated(me%c_lb)) deallocate(me%c_lb)
    if (allocated(me%c_ub)) deallocate(me%c_ub)
    allocate(me%x_lb(n), me%x_ub(n))
    allocate(me%c_lb(me%m), me%c_ub(me%m))

    end subroutine set_problem_size
!*******************************************************************************

!*******************************************************************************
!>
!  set the variable bounds \( x_l \le x \le x_u \) and the constraint
!  bounds \( c_l \le c(x) \le c_u \) (equality constraints are given by
!  `c_lb(i) == c_ub(i)`). `set_problem_size` must be called first.

    subroutine set_bounds(me, x_lb, x_ub, c_lb, c_ub)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: x_lb !! variable lower bounds `dimension(n)`
    real(wp), dimension(:), intent(in) :: x_ub !! variable upper bounds `dimension(n)`
    real(wp), dimension(:), intent(in) :: c_lb !! constraint lower bounds `dimension(m)`
    real(wp), dimension(:), intent(in) :: c_ub !! constraint upper bounds `dimension(m)`

    me%x_lb = x_lb
    me%x_ub = x_ub
    me%c_lb = c_lb
    me%c_ub = c_ub

    end subroutine set_bounds
!*******************************************************************************

!*******************************************************************************
!>
!  check that the problem is completely and consistently defined, and
!  normalize the bounds: any bound with magnitude `>= sqpopt_infinity` is
!  clamped to `+/-sqpopt_infinity`, so that e.g. `-huge(1.0_wp)` may be
!  used for an absent bound without overflowing expressions like
!  `x_ub-x_lb`. Returns `istat=sqpopt_invalid_input` (and a description in
!  `msg`) on the first problem found, else `sqpopt_success`.

    subroutine validate(me, istat, msg)

    class(sqpopt_problem_type), intent(inout) :: me
    integer,                       intent(out) :: istat !! `sqpopt_success` or `sqpopt_invalid_input`
    character(len=:), allocatable, intent(out) :: msg   !! description of the problem found (empty if none)

    istat = sqpopt_invalid_input
    msg   = ''

    if (me%n <= 0) then
        msg = 'the number of variables n must be > 0 (call set_problem_size)'; return
    end if
    if (me%m < 0 .or. me%m_eq < 0 .or. me%m_ineq < 0) then
        msg = 'the number of constraints must be >= 0'; return
    end if
    if (.not. (allocated(me%x_lb) .and. allocated(me%x_ub) .and. allocated(me%c_lb) .and. allocated(me%c_ub))) then
        msg = 'the bounds have not been set (call set_bounds)'; return
    end if
    if (size(me%x_lb) /= me%n .or. size(me%x_ub) /= me%n) then
        msg = 'the variable bounds x_lb/x_ub must have size n'; return
    end if
    if (size(me%c_lb) /= me%m .or. size(me%c_ub) /= me%m) then
        msg = 'the constraint bounds c_lb/c_ub must have size m'; return
    end if

    me%x_lb = max(me%x_lb, -sqpopt_infinity)
    me%x_ub = min(me%x_ub,  sqpopt_infinity)
    me%c_lb = max(me%c_lb, -sqpopt_infinity)
    me%c_ub = min(me%c_ub,  sqpopt_infinity)

    if (.not. (sqpopt_all_finite(me%x_lb) .and. sqpopt_all_finite(me%x_ub) .and. &
               sqpopt_all_finite(me%c_lb) .and. sqpopt_all_finite(me%c_ub))) then
        msg = 'a bound is NaN'; return
    end if
    if (any(me%x_lb > me%x_ub)) then
        msg = 'a variable lower bound exceeds its upper bound'; return
    end if
    if (any(me%c_lb > me%c_ub)) then
        msg = 'a constraint lower bound exceeds its upper bound'; return
    end if
    if (any(me%x_lb >= sqpopt_infinity) .or. any(me%x_ub <= -sqpopt_infinity) .or. &
        any(me%c_lb >= sqpopt_infinity) .or. any(me%c_ub <= -sqpopt_infinity)) then
        msg = 'a lower bound is +infinity or an upper bound is -infinity'; return
    end if

    if (me%m > 0) then
        if (.not. (allocated(me%jac_irow) .and. allocated(me%jac_icol))) then
            msg = 'the Jacobian sparsity pattern has not been set (call set_jacobian_sparsity)'; return
        end if
    end if
    if (me%jac_nnz < 0) then
        msg = 'the number of Jacobian nonzeros must be >= 0'; return
    end if
    if (me%jac_nnz > 0) then
        if (.not. (allocated(me%jac_irow) .and. allocated(me%jac_icol))) then
            msg = 'the Jacobian sparsity pattern has not been set (call set_jacobian_sparsity)'; return
        end if
        if (size(me%jac_irow) /= me%jac_nnz .or. size(me%jac_icol) /= me%jac_nnz) then
            msg = 'the Jacobian sparsity pattern arrays must have size jac_nnz'; return
        end if
        if (any(me%jac_irow < 1) .or. any(me%jac_irow > me%m)) then
            msg = 'a Jacobian row index is outside 1..m'; return
        end if
        if (any(me%jac_icol < 1) .or. any(me%jac_icol > me%n)) then
            msg = 'a Jacobian column index is outside 1..n'; return
        end if
    else
        ! make sure the (possibly empty) pattern arrays exist, so they can be copied safely:
        if (.not. allocated(me%jac_irow)) allocate(me%jac_irow(0))
        if (.not. allocated(me%jac_icol)) allocate(me%jac_icol(0))
    end if

    if (.not. (associated(me%eval_fc) .and. associated(me%eval_gjac))) then
        msg = 'the problem functions have not been set (call set_functions)'; return
    end if

    istat = sqpopt_success

    end subroutine validate
!*******************************************************************************

!*******************************************************************************
!>
!  set the (fixed) sparsity pattern of the constraint Jacobian, given as
!  1-based COO `irow`/`icol` triplets.

    subroutine set_jacobian_sparsity(me, nnz, irow, icol)

    class(sqpopt_problem_type), intent(inout) :: me
    integer, intent(in) :: nnz  !! number of nonzero Jacobian elements
    integer, dimension(:), intent(in) :: irow  !! row indices `dimension(nnz)`
    integer, dimension(:), intent(in) :: icol  !! column indices `dimension(nnz)`

    me%jac_nnz = nnz
    me%jac_irow = irow(1:nnz)
    me%jac_icol = icol(1:nnz)

    end subroutine set_jacobian_sparsity
!*******************************************************************************

!*******************************************************************************
!>
!  set the (fixed) sparsity pattern of the Hessian of the Lagrangian,
!  given as 1-based COO `irow`/`icol` triplets (only the lower triangle
!  need be supplied, since the Hessian is symmetric).

    subroutine set_hessian_sparsity(me, nnz, irow, icol)

    class(sqpopt_problem_type), intent(inout) :: me
    integer, intent(in) :: nnz  !! number of nonzero Hessian elements
    integer, dimension(:), intent(in) :: irow  !! row indices `dimension(nnz)`
    integer, dimension(:), intent(in) :: icol  !! column indices `dimension(nnz)`

    me%hess_nnz = nnz
    me%hess_irow = irow(1:nnz)
    me%hess_icol = icol(1:nnz)

    end subroutine set_hessian_sparsity
!*******************************************************************************

!*******************************************************************************
!>
!  attach the user-supplied procedures used to evaluate the objective
!  function and constraints (`fc`), their derivatives (`gjac`), and
!  optionally the exact Hessian of the Lagrangian (`hess`), and
!  (optionally) a user data object that is passed to each of them (see the
!  module-level documentation). `fc` returns \( f \) and \( c \) together,
!  and `gjac` returns \( \nabla f \) and the Jacobian values together,
!  since they usually share intermediate results.
!
!  `data` is *pointed to*, not copied, so the functions see (and may
!  update) the caller's object: it must have the `target` (or `pointer`)
!  attribute and exist for as long as the problem is being solved.

    subroutine set_functions(me, fc, gjac, hess, data)

    class(sqpopt_problem_type), intent(inout) :: me
    procedure(sqpopt_fc_func)                 :: fc   !! objective function and constraint vector
    procedure(sqpopt_gjac_func)               :: gjac !! objective gradient and sparse constraint Jacobian values
    procedure(sqpopt_hessian_func), optional  :: hess !! sparse exact Hessian of the Lagrangian values
    class(*), target, optional, intent(inout) :: data !! user data passed to each function

    me%eval_fc   => fc
    me%eval_gjac => gjac
    me%eval_hess => null()
    me%user_data => null()
    if (present(hess)) me%eval_hess => hess
    if (present(data)) me%user_data => data

    end subroutine set_functions
!*******************************************************************************

!*******************************************************************************
!>
!  empty the evaluation caches (sizing them for the current `n`, `m`, and
!  Jacobian pattern), zero the evaluation counters, clear any stop
!  request, and reset the scale factors to 1.

    subroutine reset_evaluations(me)

    class(sqpopt_problem_type), intent(inout) :: me

    if (allocated(me%cache_x))   deallocate(me%cache_x)
    if (allocated(me%cache_f))   deallocate(me%cache_f)
    if (allocated(me%cache_c))   deallocate(me%cache_c)
    if (allocated(me%cache_xg))  deallocate(me%cache_xg)
    if (allocated(me%cache_g))   deallocate(me%cache_g)
    if (allocated(me%cache_jac)) deallocate(me%cache_jac)
    if (allocated(me%c_scale))   deallocate(me%c_scale)
    allocate(me%cache_x(me%n,cache_size), me%cache_f(cache_size), me%cache_c(me%m,cache_size))
    allocate(me%cache_xg(me%n), me%cache_g(me%n), me%cache_jac(max(me%jac_nnz,0)))
    allocate(me%c_scale(me%m))
    me%cache_n = 0; me%cache_next = 1
    me%have_gjac = .false.
    me%n_eval_fc = 0; me%n_eval_gjac = 0
    me%stop_requested = .false.
    me%f_scale = 1.0_wp
    me%c_scale = 1.0_wp

    end subroutine reset_evaluations
!*******************************************************************************

!*******************************************************************************
!>
!  gradient-based scaling (as in IPOPT): at the point `x`, the objective
!  is scaled by \( s_f = \min(1, g_{max}/\lVert \nabla f \rVert_\infty) \)
!  and each constraint by \( s_{c,i} = \min(1, g_{max}/\lVert \nabla c_i
!  \rVert_\infty) \), so that no scaled gradient is larger than `max_gradient`
!  (problems whose gradients are all at most `max_gradient` are not
!  scaled). The (finite) constraint bounds are scaled to match. Call
!  [[reset_evaluations]] first; the evaluations used here are cached, so
!  they are not repeated by the first major iteration.

    subroutine compute_scaling(me, x, max_gradient)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: x            !! point at which to measure the gradients `dimension(n)`
    real(wp),               intent(in) :: max_gradient !! \( g_{max} \)

    real(wp), dimension(me%n) :: g
    real(wp), dimension(max(me%jac_nnz,0)) :: jval
    real(wp), dimension(me%m) :: row_max
    integer :: k

    call raw_gjac(me, x, g, jval)
    if (sqpopt_all_finite(g) .and. me%n > 0) then
        if (maxval(abs(g)) > max_gradient) me%f_scale = max_gradient/maxval(abs(g))
    end if

    if (me%m > 0) then
        if (sqpopt_all_finite(jval)) then
            row_max = 0.0_wp
            do k = 1, me%jac_nnz
                row_max(me%jac_irow(k)) = max(row_max(me%jac_irow(k)), abs(jval(k)))
            end do
            where (row_max > max_gradient) me%c_scale = max_gradient/row_max
        end if
        where (abs(me%c_lb) < sqpopt_infinity) me%c_lb = me%c_scale*me%c_lb
        where (abs(me%c_ub) < sqpopt_infinity) me%c_ub = me%c_scale*me%c_ub
    end if

    end subroutine compute_scaling
!*******************************************************************************

!*******************************************************************************
!>
!  the (scaled) objective \( s_f f(x) \). The underlying value is reused
!  from one of the last few evaluations if `x` is exactly (bitwise) the
!  same point -- the solver evaluates `f` at each trial point of a line
!  search and then again at the start of the next major iteration, at the
!  accepted point, so this saves one call per major iteration. NaN if the
!  user function failed (`status /= 0`).

    subroutine eval_f_cached(me, x, f)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x  !! point `dimension(n)`
    real(wp),                   intent(out)   :: f  !! scaled objective function value at `x`

    real(wp), dimension(me%m) :: c

    call raw_fc(me, x, f, c)
    f = me%f_scale*f

    end subroutine eval_f_cached
!*******************************************************************************

!*******************************************************************************
!>
!  the (scaled) constraints \( s_{c,i} c_i(x) \), with the same caching as
!  [[eval_f_cached]].

    subroutine eval_c_cached(me, x, c)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x  !! point `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: c  !! scaled constraint values at `x` `dimension(m)`

    real(wp) :: f

    call raw_fc(me, x, f, c)
    c = me%c_scale*c

    end subroutine eval_c_cached
!*******************************************************************************

!*******************************************************************************
!>
!  the (scaled) objective gradient \( s_f \nabla f(x) \), reusing the last
!  evaluation if `x` is the same point.

    subroutine eval_g_cached(me, x, g)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x  !! point `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: g  !! scaled gradient at `x` `dimension(n)`

    real(wp), dimension(max(me%jac_nnz,0)) :: jac_val

    call raw_gjac(me, x, g, jac_val)
    g = me%f_scale*g

    end subroutine eval_g_cached
!*******************************************************************************

!*******************************************************************************
!>
!  the (scaled) Jacobian values \( s_{c,i} \partial c_i/\partial x_j \),
!  reusing the last evaluation if `x` is the same point.

    subroutine eval_jac_cached(me, x, jac_val)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:),     intent(out)   :: jac_val !! scaled Jacobian values `dimension(jac_nnz)`

    integer :: k
    real(wp), dimension(me%n) :: g

    call raw_gjac(me, x, g, jac_val)
    do k = 1, me%jac_nnz
        jac_val(k) = me%c_scale(me%jac_irow(k))*jac_val(k)
    end do

    end subroutine eval_jac_cached
!*******************************************************************************

!*******************************************************************************
!>
!  the unscaled objective and constraints, from the cache or the user's `fc`
!  (whose call is counted, and whose `status` is handled by [[check_status]]).

    subroutine raw_fc(me, x, f, c)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x
    real(wp),                   intent(out)   :: f
    real(wp), dimension(:),     intent(out)   :: c

    integer :: k, status
    real(wp), dimension(1+size(c)) :: v

    if (.not. allocated(me%cache_x)) call reset_evaluations(me)
    do k = 1, me%cache_n
        if (all(me%cache_x(:,k) == x)) then
            f = me%cache_f(k)
            c = me%cache_c(:,k)
            return
        end if
    end do

    if (me%stop_requested) then   ! (no more user calls once a stop has been requested)
        f = ieee_value(1.0_wp, ieee_quiet_nan)
        c = f
        return
    end if
    status = 0
    if (associated(me%user_data)) then
        call me%eval_fc(x, f, c, status, me%user_data)
    else
        call me%eval_fc(x, f, c, status)
    end if
    me%n_eval_fc = me%n_eval_fc + 1
    v = [f, c]
    call check_status(me, status, v)
    f = v(1)
    c = v(2:)

    me%cache_x(:,me%cache_next) = x
    me%cache_f(me%cache_next)   = f
    me%cache_c(:,me%cache_next) = c
    me%cache_n    = min(me%cache_n + 1, cache_size)
    me%cache_next = mod(me%cache_next, cache_size) + 1

    end subroutine raw_fc
!*******************************************************************************

!*******************************************************************************
!>
!  the unscaled gradient and Jacobian values, from the one-entry cache or
!  the user's `gjac` (whose call is counted, and whose `status` is handled
!  by [[check_status]]).

    subroutine raw_gjac(me, x, g, jac_val)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x
    real(wp), dimension(:),     intent(out)   :: g
    real(wp), dimension(:),     intent(out)   :: jac_val

    integer :: status
    real(wp), dimension(size(g)+size(jac_val)) :: v

    if (.not. allocated(me%cache_xg)) call reset_evaluations(me)
    if (me%have_gjac) then
        if (all(me%cache_xg == x)) then
            g       = me%cache_g
            jac_val = me%cache_jac
            return
        end if
    end if

    if (me%stop_requested) then
        g       = ieee_value(1.0_wp, ieee_quiet_nan)
        jac_val = ieee_value(1.0_wp, ieee_quiet_nan)
        return
    end if
    status = 0
    if (associated(me%user_data)) then
        call me%eval_gjac(x, g, jac_val, status, me%user_data)
    else
        call me%eval_gjac(x, g, jac_val, status)
    end if
    me%n_eval_gjac = me%n_eval_gjac + 1
    v = [g, jac_val]
    call check_status(me, status, v)
    g       = v(1:size(g))
    jac_val = v(size(g)+1:)

    me%cache_xg  = x
    me%cache_g   = g
    me%cache_jac = jac_val
    me%have_gjac = .true.

    end subroutine raw_gjac
!*******************************************************************************

!*******************************************************************************
!>
!  handle a user function's returned `status`: if it is nonzero the
!  outputs `v` are replaced by NaN (so the solver rejects the point, or
!  reports `sqpopt_function_error` if it is the current point), and if it
!  is negative a stop is requested (after which no user function is called
!  again: uncached evaluations just return NaN).

    subroutine check_status(me, status, v)

    class(sqpopt_problem_type), intent(inout) :: me
    integer,                    intent(in)    :: status
    real(wp), dimension(:),     intent(inout) :: v

    if (status == 0) return
    if (status < 0) me%stop_requested = .true.
    v = ieee_value(1.0_wp, ieee_quiet_nan)

    end subroutine check_status
!*******************************************************************************

    end module sqpopt_problem_module
!*******************************************************************************
