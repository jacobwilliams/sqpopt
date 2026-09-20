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

    module sqpopt_problem_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    type, public :: sqpopt_problem_type
        !! defines the problem to be solved: the problem size, the
        !! variable and constraint bounds, and the user-supplied
        !! procedures used to evaluate the objective function,
        !! constraints, and their derivatives.

        integer :: n      = 0  !! number of optimization variables (\( n>0 \))
        integer :: m      = 0  !! total number of nonlinear constraints \( c(x) \) (\( m \ge 0 \))
        integer :: m_eq   = 0  !! number of nonlinear equality constraints (\( c_l = c_u \))
        integer :: m_ineq = 0  !! number of nonlinear inequality constraints (\( c_l \ne c_u \))

        real(wp), dimension(:), allocatable :: x_lb  !! lower bounds on the optimization variables `dimension(n)`
        real(wp), dimension(:), allocatable :: x_ub  !! upper bounds on the optimization variables `dimension(n)`

        real(wp), dimension(:), allocatable :: c_lb  !! lower bounds on the constraints `dimension(m)`
        real(wp), dimension(:), allocatable :: c_ub  !! upper bounds on the constraints `dimension(m)`

        procedure(sqpopt_objective_func), pointer, nopass :: eval_f    => null() !! evaluates \( f(x) \)
        procedure(sqpopt_gradient_func),  pointer, nopass :: eval_g    => null() !! evaluates \( \nabla f(x) \)
        procedure(sqpopt_constraint_func),pointer, nopass :: eval_c    => null() !! evaluates \( c(x) \)
        procedure(sqpopt_jacobian_func),  pointer, nopass :: eval_jac  => null() !! evaluates the Jacobian of \( c(x) \)
        procedure(sqpopt_hessian_func),   pointer, nopass :: eval_hess => null() !! evaluates the Hessian of the Lagrangian (only used when an exact Hessian is requested)

        contains

        procedure, public :: set_problem_size  !! set the problem dimensions and allocate the bound arrays
        procedure, public :: set_functions     !! attach the user-supplied evaluation procedures

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

        subroutine sqpopt_jacobian_func(x, jac)
            !! evaluates the Jacobian of the constraint vector: \( J_{ij} = \partial c_i / \partial x_j \)
            import :: wp
            implicit none
            real(wp), dimension(:),   intent(in)  :: x    !! optimization variable vector `dimension(n)`
            real(wp), dimension(:,:), intent(out) :: jac  !! constraint Jacobian `dimension(m,n)`
        end subroutine sqpopt_jacobian_func

        subroutine sqpopt_hessian_func(x, lambda, h)
            !! evaluates the Hessian of the Lagrangian:
            !! \( H = \nabla^2 f(x) - \sum_i \lambda_i \nabla^2 c_i(x) \)
            import :: wp
            implicit none
            real(wp), dimension(:),   intent(in)  :: x      !! optimization variable vector `dimension(n)`
            real(wp), dimension(:),   intent(in)  :: lambda !! Lagrange multipliers `dimension(m)`
            real(wp), dimension(:,:), intent(out) :: h      !! Hessian of the Lagrangian `dimension(n,n)`
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
!  attach the user-supplied procedures used to evaluate the objective
!  function, constraints, and (optionally) their derivatives.

    subroutine set_functions(me, f, g, c, jac, hess)

    class(sqpopt_problem_type), intent(inout) :: me
    procedure(sqpopt_objective_func)  :: f     !! objective function
    procedure(sqpopt_gradient_func)   :: g     !! objective function gradient
    procedure(sqpopt_constraint_func) :: c     !! constraint vector
    procedure(sqpopt_jacobian_func)   :: jac   !! constraint Jacobian
    procedure(sqpopt_hessian_func), optional :: hess !! exact Hessian of the Lagrangian (optional)

    ! TODO: implement

    end subroutine set_functions
!*******************************************************************************

    end module sqpopt_problem_module
!*******************************************************************************
