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

    module sqpopt_problem_module

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_success, sqpopt_invalid_input, sqpopt_infinity

    implicit none

    private

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

        integer :: jac_nnz = 0  !! number of nonzero elements in the constraint Jacobian
        integer, dimension(:), allocatable :: jac_irow  !! Jacobian sparsity pattern: row indices `dimension(jac_nnz)`
        integer, dimension(:), allocatable :: jac_icol  !! Jacobian sparsity pattern: column indices `dimension(jac_nnz)`

        integer :: hess_nnz = 0  !! number of nonzero elements in the (symmetric) Hessian of the Lagrangian
        integer, dimension(:), allocatable :: hess_irow !! Hessian sparsity pattern: row indices `dimension(hess_nnz)`
        integer, dimension(:), allocatable :: hess_icol !! Hessian sparsity pattern: column indices `dimension(hess_nnz)`

        procedure(sqpopt_objective_func), pointer, nopass :: eval_f    => null() !! evaluates \( f(x) \)
        procedure(sqpopt_gradient_func),  pointer, nopass :: eval_g    => null() !! evaluates \( \nabla f(x) \)
        procedure(sqpopt_constraint_func),pointer, nopass :: eval_c    => null() !! evaluates \( c(x) \)
        procedure(sqpopt_jacobian_func),  pointer, nopass :: eval_jac  => null() !! evaluates the nonzero values of the Jacobian of \( c(x) \)
        procedure(sqpopt_hessian_func),   pointer, nopass :: eval_hess => null() !! evaluates the nonzero values of the Hessian of the Lagrangian (only used when an exact Hessian is requested)

        contains

        procedure, public :: set_problem_size      !! set the problem dimensions and allocate the bound arrays
        procedure, public :: set_bounds            !! set the variable and constraint bounds
        procedure, public :: set_jacobian_sparsity !! set the (fixed) sparsity pattern of the constraint Jacobian
        procedure, public :: set_hessian_sparsity  !! set the (fixed) sparsity pattern of the Lagrangian Hessian
        procedure, public :: set_functions         !! attach the user-supplied evaluation procedures
        procedure, public :: validate              !! check the problem definition and normalize infinite bounds

    end type sqpopt_problem_type

    abstract interface

        subroutine sqpopt_objective_func(x, f)
            !! evaluates the objective function \( f(x) \)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in) :: x  !! optimization variable vector `dimension(n)`
            real(wp), intent(out)              :: f  !! value of the objective function
        end subroutine sqpopt_objective_func

        subroutine sqpopt_gradient_func(x, g)
            !! evaluates the gradient of the objective function \( \nabla f(x) \)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)  :: x  !! optimization variable vector `dimension(n)`
            real(wp), dimension(:), intent(out) :: g  !! gradient vector `dimension(n)`
        end subroutine sqpopt_gradient_func

        subroutine sqpopt_constraint_func(x, c)
            !! evaluates the nonlinear constraint vector \( c(x) \)
            !! (equality constraints first, followed by inequality constraints)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)  :: x  !! optimization variable vector `dimension(n)`
            real(wp), dimension(:), intent(out) :: c  !! constraint vector `dimension(m)`
        end subroutine sqpopt_constraint_func

        subroutine sqpopt_jacobian_func(x, jac_val)
            !! evaluates the nonzero values of the Jacobian of the constraint
            !! vector: \( J_{ij} = \partial c_i / \partial x_j \), ordered to
            !! match the sparsity pattern set by `set_jacobian_sparsity`.
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)  :: x       !! optimization variable vector `dimension(n)`
            real(wp), dimension(:), intent(out) :: jac_val !! nonzero Jacobian values `dimension(jac_nnz)`
        end subroutine sqpopt_jacobian_func

        subroutine sqpopt_hessian_func(x, lambda, hess_val)
            !! evaluates the nonzero values of the Hessian of the Lagrangian:
            !! \( H = \nabla^2 f(x) - \sum_i \lambda_i \nabla^2 c_i(x) \),
            !! ordered to match the sparsity pattern set by `set_hessian_sparsity`.
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in)  :: x        !! optimization variable vector `dimension(n)`
            real(wp), dimension(:), intent(in)  :: lambda   !! Lagrange multipliers `dimension(m)`
            real(wp), dimension(:), intent(out) :: hess_val !! nonzero Hessian values `dimension(hess_nnz)`
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

    if (any(me%x_lb /= me%x_lb) .or. any(me%x_ub /= me%x_ub) .or. &
        any(me%c_lb /= me%c_lb) .or. any(me%c_ub /= me%c_ub)) then
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
!  function, constraints, and (optionally) their derivatives.

    subroutine set_functions(me, f, g, c, jac, hess)

    class(sqpopt_problem_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: f     !! objective function
    procedure(sqpopt_gradient_func)   :: g     !! objective function gradient
    procedure(sqpopt_constraint_func) :: c     !! constraint vector
    procedure(sqpopt_jacobian_func)   :: jac   !! sparse constraint Jacobian values
    procedure(sqpopt_hessian_func), optional :: hess !! sparse exact Hessian of the Lagrangian values (optional)

    me%eval_f   => f
    me%eval_g   => g
    me%eval_c   => c
    me%eval_jac => jac
    if (present(hess)) me%eval_hess => hess

    end subroutine set_functions
!*******************************************************************************

    end module sqpopt_problem_module
!*******************************************************************************
