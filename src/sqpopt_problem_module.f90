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
!  `irow`/`icol` triplets (the same convention used by the `lusol`,
!  `LSQR`, and `LSMR` dependencies), and the user-supplied evaluation
!  routines only need to fill in the corresponding nonzero *values* on
!  each call.
!
!  **User functions.** Every user function has two trailing arguments:
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
!  point accepted by a line search is not evaluated again), and apply the
!  objective and constraint scaling (see [[compute_scaling]]): the solver
!  works with \( s_f f \) and \( s_{c,i} c_i \).

    module sqpopt_problem_module

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_success, sqpopt_invalid_input, sqpopt_infinity, sqpopt_all_finite
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    implicit none

    private

    integer, parameter :: cache_size = 4 !! number of recent evaluations of `f` (and of `c`) kept

    public :: sqpopt_objective_func, sqpopt_gradient_func, sqpopt_constraint_func
    public :: sqpopt_jacobian_func, sqpopt_hessian_func

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

        procedure(sqpopt_objective_func), pointer, nopass :: eval_f    => null() !! evaluates \( f(x) \)
        procedure(sqpopt_gradient_func),  pointer, nopass :: eval_g    => null() !! evaluates \( \nabla f(x) \)
        procedure(sqpopt_constraint_func),pointer, nopass :: eval_c    => null() !! evaluates \( c(x) \)
        procedure(sqpopt_jacobian_func),  pointer, nopass :: eval_jac  => null() !! evaluates the nonzero values of the
                                                                                !! Jacobian of \( c(x) \)
        procedure(sqpopt_hessian_func),   pointer, nopass :: eval_hess => null() !! evaluates the nonzero values of the
                                                                                !! Hessian of the Lagrangian (only used
                                                                                !! when an exact Hessian is requested)
        class(*), pointer :: user_data => null() !! user data passed to every user function (see [[set_functions]])

        ! ---- internal evaluation state (see the module-level documentation) ----
        logical  :: stop_requested = .false. !! set when a user function returns `status < 0`
        integer  :: n_eval_f   = 0 !! number of calls of the user's `f`
        integer  :: n_eval_g   = 0 !! number of calls of the user's `g`
        integer  :: n_eval_c   = 0 !! number of calls of the user's `c`
        integer  :: n_eval_jac = 0 !! number of calls of the user's `jac`
        real(wp) :: f_scale = 1.0_wp !! objective scale factor \( s_f \)
        real(wp), dimension(:), allocatable :: c_scale !! constraint scale factors \( s_{c,i} \) `dimension(m)`
        integer :: cache_nf = 0, cache_next_f = 1 !! entries used, and the slot for the next one, in the `f` cache
        integer :: cache_nc = 0, cache_next_c = 1 !! the same, for the `c` cache
        real(wp), dimension(:,:), allocatable :: cache_xf !! points at which `f` was evaluated `dimension(n,cache_size)`
        real(wp), dimension(:),   allocatable :: cache_f  !! (unscaled) `f` at those points `dimension(cache_size)`
        real(wp), dimension(:,:), allocatable :: cache_xc !! points at which `c` was evaluated `dimension(n,cache_size)`
        real(wp), dimension(:,:), allocatable :: cache_c  !! (unscaled) `c` at those points `dimension(m,cache_size)`
        logical :: have_g = .false., have_jac = .false. !! whether the one-entry `g`/`jac` caches are filled
        real(wp), dimension(:), allocatable :: cache_xg, cache_g   !! point and (unscaled) `g` there `dimension(n)`
        real(wp), dimension(:), allocatable :: cache_xj, cache_jac !! point and (unscaled) Jacobian values there

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

        subroutine sqpopt_objective_func(x, f, status, data)
            !! evaluates the objective function \( f(x) \)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)    :: x      !! optimization variable vector `dimension(n)`
            real(wp),               intent(out)   :: f      !! value of the objective function
            integer,                intent(inout) :: status !! `0` on entry; `>0`: can't evaluate here, `<0`: stop
            class(*), optional,     intent(inout) :: data   !! user data (see [[set_functions]])
        end subroutine sqpopt_objective_func

        subroutine sqpopt_gradient_func(x, g, status, data)
            !! evaluates the gradient of the objective function \( \nabla f(x) \)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)    :: x      !! optimization variable vector `dimension(n)`
            real(wp), dimension(:), intent(out)   :: g      !! gradient vector `dimension(n)`
            integer,                intent(inout) :: status !! `0` on entry; `>0`: can't evaluate here, `<0`: stop
            class(*), optional,     intent(inout) :: data   !! user data (see [[set_functions]])
        end subroutine sqpopt_gradient_func

        subroutine sqpopt_constraint_func(x, c, status, data)
            !! evaluates the nonlinear constraint vector \( c(x) \)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)    :: x      !! optimization variable vector `dimension(n)`
            real(wp), dimension(:), intent(out)   :: c      !! constraint vector `dimension(m)`
            integer,                intent(inout) :: status !! `0` on entry; `>0`: can't evaluate here, `<0`: stop
            class(*), optional,     intent(inout) :: data   !! user data (see [[set_functions]])
        end subroutine sqpopt_constraint_func

        subroutine sqpopt_jacobian_func(x, jac_val, status, data)
            !! evaluates the nonzero values of the Jacobian of the constraint
            !! vector: \( J_{ij} = \partial c_i / \partial x_j \), ordered to
            !! match the sparsity pattern set by `set_jacobian_sparsity`.
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)    :: x       !! optimization variable vector `dimension(n)`
            real(wp), dimension(:), intent(out)   :: jac_val !! nonzero Jacobian values `dimension(jac_nnz)`
            integer,                intent(inout) :: status  !! `0` on entry; `>0`: can't evaluate here, `<0`: stop
            class(*), optional,     intent(inout) :: data    !! user data (see [[set_functions]])
        end subroutine sqpopt_jacobian_func

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

    if (.not. (associated(me%eval_f) .and. associated(me%eval_g) .and. &
               associated(me%eval_c) .and. associated(me%eval_jac))) then
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
!  function, constraints, and their derivatives, and (optionally) a user
!  data object that is passed to each of them (see the module-level
!  documentation).
!
!  `data` is *pointed to*, not copied, so the functions see (and may
!  update) the caller's object: it must have the `target` (or `pointer`)
!  attribute and exist for as long as the problem is being solved.

    subroutine set_functions(me, f, g, c, jac, hess, data)

    class(sqpopt_problem_type), intent(inout) :: me
    procedure(sqpopt_objective_func)          :: f    !! objective function
    procedure(sqpopt_gradient_func)           :: g    !! objective function gradient
    procedure(sqpopt_constraint_func)         :: c    !! constraint vector
    procedure(sqpopt_jacobian_func)           :: jac  !! sparse constraint Jacobian values
    procedure(sqpopt_hessian_func), optional  :: hess !! sparse exact Hessian of the Lagrangian values
    class(*), target, optional, intent(inout) :: data !! user data passed to each function

    me%eval_f    => f
    me%eval_g    => g
    me%eval_c    => c
    me%eval_jac  => jac
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

    if (allocated(me%cache_xf))  deallocate(me%cache_xf)
    if (allocated(me%cache_f))   deallocate(me%cache_f)
    if (allocated(me%cache_xc))  deallocate(me%cache_xc)
    if (allocated(me%cache_c))   deallocate(me%cache_c)
    if (allocated(me%cache_xg))  deallocate(me%cache_xg)
    if (allocated(me%cache_g))   deallocate(me%cache_g)
    if (allocated(me%cache_xj))  deallocate(me%cache_xj)
    if (allocated(me%cache_jac)) deallocate(me%cache_jac)
    if (allocated(me%c_scale))   deallocate(me%c_scale)
    allocate(me%cache_xf(me%n,cache_size), me%cache_f(cache_size))
    allocate(me%cache_xc(me%n,cache_size), me%cache_c(me%m,cache_size))
    allocate(me%cache_xg(me%n), me%cache_g(me%n), me%cache_xj(me%n), me%cache_jac(max(me%jac_nnz,0)))
    allocate(me%c_scale(me%m))
    me%cache_nf = 0; me%cache_next_f = 1
    me%cache_nc = 0; me%cache_next_c = 1
    me%have_g   = .false.
    me%have_jac = .false.
    me%n_eval_f = 0; me%n_eval_g = 0; me%n_eval_c = 0; me%n_eval_jac = 0
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

    call raw_g(me, x, g)
    if (sqpopt_all_finite(g) .and. me%n > 0) then
        if (maxval(abs(g)) > max_gradient) me%f_scale = max_gradient/maxval(abs(g))
    end if

    if (me%m > 0) then
        call raw_jac(me, x, jval)
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

    call raw_f(me, x, f)
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

    call raw_c(me, x, c)
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

    call raw_g(me, x, g)
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

    call raw_jac(me, x, jac_val)
    do k = 1, me%jac_nnz
        jac_val(k) = me%c_scale(me%jac_irow(k))*jac_val(k)
    end do

    end subroutine eval_jac_cached
!*******************************************************************************

!*******************************************************************************
!>
!  the unscaled objective, from the cache or the user function.

    subroutine raw_f(me, x, f)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x
    real(wp),                   intent(out)   :: f

    integer :: k

    if (.not. allocated(me%cache_xf)) call reset_evaluations(me)
    do k = 1, me%cache_nf
        if (all(me%cache_xf(:,k) == x)) then
            f = me%cache_f(k)
            return
        end if
    end do

    if (me%stop_requested) then   ! (no more user calls once a stop has been requested)
        f = ieee_value(1.0_wp, ieee_quiet_nan)
        return
    end if
    call call_f(me, x, f)
    me%cache_xf(:,me%cache_next_f) = x
    me%cache_f(me%cache_next_f)    = f
    me%cache_nf     = min(me%cache_nf + 1, cache_size)
    me%cache_next_f = mod(me%cache_next_f, cache_size) + 1

    end subroutine raw_f
!*******************************************************************************

!*******************************************************************************
!>
!  the unscaled constraints, from the cache or the user function.

    subroutine raw_c(me, x, c)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x
    real(wp), dimension(:),     intent(out)   :: c

    integer :: k

    if (.not. allocated(me%cache_xc)) call reset_evaluations(me)
    if (me%m == 0) return
    do k = 1, me%cache_nc
        if (all(me%cache_xc(:,k) == x)) then
            c = me%cache_c(:,k)
            return
        end if
    end do

    if (me%stop_requested) then
        c = ieee_value(1.0_wp, ieee_quiet_nan)
        return
    end if
    call call_c(me, x, c)
    me%cache_xc(:,me%cache_next_c) = x
    me%cache_c(:,me%cache_next_c)  = c
    me%cache_nc     = min(me%cache_nc + 1, cache_size)
    me%cache_next_c = mod(me%cache_next_c, cache_size) + 1

    end subroutine raw_c
!*******************************************************************************

!*******************************************************************************
!>
!  the unscaled gradient, from the one-entry cache or the user function.

    subroutine raw_g(me, x, g)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x
    real(wp), dimension(:),     intent(out)   :: g

    integer :: status

    if (.not. allocated(me%cache_xg)) call reset_evaluations(me)
    if (me%have_g) then
        if (all(me%cache_xg == x)) then
            g = me%cache_g
            return
        end if
    end if

    if (me%stop_requested) then
        g = ieee_value(1.0_wp, ieee_quiet_nan)
        return
    end if
    status = 0
    if (associated(me%user_data)) then
        call me%eval_g(x, g, status, me%user_data)
    else
        call me%eval_g(x, g, status)
    end if
    me%n_eval_g = me%n_eval_g + 1
    call check_status(me, status, g)

    me%cache_xg = x
    me%cache_g  = g
    me%have_g   = .true.

    end subroutine raw_g
!*******************************************************************************

!*******************************************************************************
!>
!  the unscaled Jacobian values, from the one-entry cache or the user function.

    subroutine raw_jac(me, x, jac_val)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x
    real(wp), dimension(:),     intent(out)   :: jac_val

    integer :: status

    if (.not. allocated(me%cache_xj)) call reset_evaluations(me)
    if (me%have_jac) then
        if (all(me%cache_xj == x)) then
            jac_val = me%cache_jac
            return
        end if
    end if

    if (me%stop_requested) then
        jac_val = ieee_value(1.0_wp, ieee_quiet_nan)
        return
    end if
    status = 0
    if (associated(me%user_data)) then
        call me%eval_jac(x, jac_val, status, me%user_data)
    else
        call me%eval_jac(x, jac_val, status)
    end if
    me%n_eval_jac = me%n_eval_jac + 1
    call check_status(me, status, jac_val)

    me%cache_xj  = x
    me%cache_jac = jac_val
    me%have_jac  = .true.

    end subroutine raw_jac
!*******************************************************************************

!*******************************************************************************
!>
!  call the user's objective function (passing `status` and the user data),
!  count the call, and handle the returned status.

    subroutine call_f(me, x, f)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x
    real(wp),                   intent(out)   :: f

    integer :: status
    real(wp), dimension(1) :: fv

    status = 0
    if (associated(me%user_data)) then
        call me%eval_f(x, f, status, me%user_data)
    else
        call me%eval_f(x, f, status)
    end if
    me%n_eval_f = me%n_eval_f + 1
    fv = f
    call check_status(me, status, fv)
    f = fv(1)

    end subroutine call_f
!*******************************************************************************

!*******************************************************************************
!>
!  call the user's constraint function (passing `status` and the user data),
!  count the call, and handle the returned status.

    subroutine call_c(me, x, c)

    class(sqpopt_problem_type), intent(inout) :: me
    real(wp), dimension(:),     intent(in)    :: x
    real(wp), dimension(:),     intent(out)   :: c

    integer :: status

    status = 0
    if (associated(me%user_data)) then
        call me%eval_c(x, c, status, me%user_data)
    else
        call me%eval_c(x, c, status)
    end if
    me%n_eval_c = me%n_eval_c + 1
    call check_status(me, status, c)

    end subroutine call_c
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
