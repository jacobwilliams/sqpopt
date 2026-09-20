!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Module containing the main object-oriented interface to the `sqpopt`
!  Sequential Quadratic Programming solver. It is used via the
!  [[sqpopt_type]] class, which is the only public entity in this module.
!  The other components of the algorithm (problem definition, options,
!  Hessian approximation, QP subproblem solver, line search, and
!  convergence checking) are each implemented in their own module so
!  that they may be developed, tested, and swapped out independently.

    module sqpopt_module

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_types_module,      only: sqpopt_success, sqpopt_error
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type
    use sqpopt_hessian_module,    only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type
    use sqpopt_iterate_module,    only: sqpopt_iterate

    implicit none

    private

    type, public :: sqpopt_type
        !! main class for the SQP optimizer.

        private

        type(sqpopt_problem_type)    :: problem      !! the nonlinear program to be solved
        type(sqpopt_options_type)    :: options      !! solver options
        type(sqpopt_hessian_type)    :: hessian      !! Hessian of the Lagrangian approximation
        type(sqpopt_qp_solver_type)  :: qp_solver    !! QP subproblem solver
        type(sqpopt_linesearch_type) :: linesearch   !! merit function / line search

        real(wp), dimension(:), allocatable :: x       !! current/final optimization variables
        real(wp), dimension(:), allocatable :: lambda  !! current/final Lagrange multipliers

        integer :: iter  = 0  !! number of major iterations performed
        integer :: istat = 0  !! solver status code (see [[sqpopt_types_module]])

        contains

        private

        procedure, public :: initialize   => sqpopt_initialize
        procedure, public :: set_problem  => sqpopt_set_problem
        procedure, public :: set_options  => sqpopt_set_options
        procedure, public :: solve        => sqpopt_solve
        procedure, public :: get_solution => sqpopt_get_solution
        procedure, public :: destroy      => sqpopt_destroy

    end type sqpopt_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize (or reinitialize) an [[sqpopt_type]] solver instance.

    subroutine sqpopt_initialize(me)

    class(sqpopt_type), intent(inout) :: me

    ! TODO: implement

    end subroutine sqpopt_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  define the nonlinear program to be solved.

    subroutine sqpopt_set_problem(me, problem)

    class(sqpopt_type),        intent(inout) :: me
    type(sqpopt_problem_type), intent(in)    :: problem  !! the problem definition

    ! TODO: implement

    end subroutine sqpopt_set_problem
!*******************************************************************************

!*******************************************************************************
!>
!  set the solver options.

    subroutine sqpopt_set_options(me, options)

    class(sqpopt_type),        intent(inout) :: me
    type(sqpopt_options_type), intent(in)    :: options  !! the solver options

    ! TODO: implement

    end subroutine sqpopt_set_options
!*******************************************************************************

!*******************************************************************************
!>
!  solve the nonlinear program starting from the initial guess `x0`,
!  running major SQP iterations until convergence or a stopping
!  criterion is met.

    subroutine sqpopt_solve(me, x0, istat)

    class(sqpopt_type),     intent(inout) :: me
    real(wp), dimension(:), intent(in)    :: x0     !! initial guess for the optimization variables `dimension(n)`
    integer,                intent(out)   :: istat  !! status code (see [[sqpopt_types_module]])

    ! TODO: implement

    end subroutine sqpopt_solve
!*******************************************************************************

!*******************************************************************************
!>
!  return the current (or final) solution and associated Lagrange multipliers.

    subroutine sqpopt_get_solution(me, x, lambda)

    class(sqpopt_type),     intent(in)  :: me
    real(wp), dimension(:), intent(out) :: x       !! optimization variables `dimension(n)`
    real(wp), dimension(:), intent(out) :: lambda  !! Lagrange multipliers `dimension(m)`

    ! TODO: implement

    end subroutine sqpopt_get_solution
!*******************************************************************************

!*******************************************************************************
!>
!  destroy the solver instance, deallocating all internal arrays.

    subroutine sqpopt_destroy(me)

    class(sqpopt_type), intent(inout) :: me

    ! TODO: implement

    end subroutine sqpopt_destroy
!*******************************************************************************

    end module sqpopt_module
!*******************************************************************************