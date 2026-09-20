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

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

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
        procedure, public :: set_jacobian_sparsity !! set the (fixed) sparsity pattern of the constraint Jacobian
        procedure, public :: set_hessian_sparsity  !! set the (fixed) sparsity pattern of the Lagrangian Hessian
        procedure, public :: set_functions         !! attach the user-supplied evaluation procedures

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

    ! TODO: implement

    end subroutine set_problem_size
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

    ! TODO: implement

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

    ! TODO: implement

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

    ! TODO: implement

    end subroutine set_functions
!*******************************************************************************

    end module sqpopt_problem_module
!*******************************************************************************
