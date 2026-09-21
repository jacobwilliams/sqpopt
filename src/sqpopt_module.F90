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
!  Internally, sparse (COO) storage is used by default for the
!  constraint Jacobian and the Lagrangian Hessian is a matrix-free
!  limited-memory operator -- dense `n x n`/`m x n` arrays are never formed.

    module sqpopt_module

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_types_module,      only: sqpopt_success, sqpopt_error, sqpopt_max_iter_reached
    use sqpopt_problem_module,    only: sqpopt_problem_type
    use sqpopt_options_module,    only: sqpopt_options_type, sqpopt_hessian_sr1
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
        procedure, public :: solve        => sqpopt_solve
        procedure, public :: get_solution => sqpopt_get_solution
        procedure, public :: destroy      => sqpopt_destroy

    end type sqpopt_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize (or reinitialize) an [[sqpopt_type]] solver instance.

    subroutine sqpopt_initialize(me, problem, options, hessian, qp_solver, linesearch)

    class(sqpopt_type), intent(inout) :: me
    type(sqpopt_problem_type),optional,intent(in)    :: problem      !! the nonlinear program to be solved
    type(sqpopt_options_type),optional,intent(in)    :: options      !! solver options
    type(sqpopt_hessian_type),optional,intent(in)    :: hessian      !! Hessian of the Lagrangian approximation
    type(sqpopt_qp_solver_type),optional,intent(in)  :: qp_solver    !! QP subproblem solver
    type(sqpopt_linesearch_type),optional,intent(in) :: linesearch   !! merit function / line search

    ! use inputs if present, else use defaults:
    if (present(problem))    then; me%problem = problem; else; me%problem = sqpopt_problem_type(); end if
    if (present(options))    then; me%options = options; else; me%options = sqpopt_options_type(); end if
    if (present(hessian))    then; me%hessian = hessian; else; me%hessian = sqpopt_hessian_type(); end if
    if (present(qp_solver))  then; me%qp_solver = qp_solver; else; me%qp_solver = sqpopt_qp_solver_type(); end if
    if (present(linesearch)) then; me%linesearch = linesearch; else; me%linesearch = sqpopt_linesearch_type(); end if

    if (allocated(me%x))      deallocate(me%x)
    if (allocated(me%lambda)) deallocate(me%lambda)
    me%iter  = 0
    me%istat = 0

    end subroutine sqpopt_initialize
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

    real(wp), dimension(:), allocatable :: x_prev, gl_prev  !! quasi-Newton state (unallocated until the 2nd iteration)
    logical :: converged
    integer :: iter_istat, iter

    me%x = x0
    if (allocated(me%lambda)) deallocate(me%lambda)
    allocate(me%lambda(me%problem%m))
    me%lambda = 0.0_wp

    call me%hessian%initialize(me%problem%n, me%options%lbfgs_memory, &
                                use_sr1=(me%options%hessian_mode == sqpopt_hessian_sr1))
    me%qp_solver%linear_solver_mode = me%options%linear_solver_mode
    me%linesearch%mode              = me%options%linesearch_mode
    me%linesearch%merit_mode        = me%options%merit_mode

    do iter = 1, me%options%max_iter
        me%iter = iter
        call sqpopt_iterate(me%problem, me%options, me%hessian, me%qp_solver, me%linesearch, &
                             me%x, me%lambda, x_prev, gl_prev, converged, iter_istat)
        if (converged) then
            istat = sqpopt_success
            me%istat = istat
            return
        end if
    end do

    istat    = sqpopt_max_iter_reached
    me%istat = istat

    end subroutine sqpopt_solve
!*******************************************************************************

!*******************************************************************************
!>
!  return the current (or final) solution and associated Lagrange multipliers.

    subroutine sqpopt_get_solution(me, x, lambda)

    class(sqpopt_type),     intent(in)  :: me
    real(wp), dimension(:), intent(out) :: x       !! optimization variables `dimension(n)`
    real(wp), dimension(:), intent(out) :: lambda  !! Lagrange multipliers `dimension(m)`

    x      = me%x
    lambda = me%lambda

    end subroutine sqpopt_get_solution
!*******************************************************************************

!*******************************************************************************
!>
!  destroy the solver instance, deallocating all internal arrays.

    subroutine sqpopt_destroy(me)

    class(sqpopt_type), intent(inout) :: me

    call me%initialize()

    end subroutine sqpopt_destroy
!*******************************************************************************

    end module sqpopt_module
!*******************************************************************************