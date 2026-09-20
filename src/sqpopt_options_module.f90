!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  User-settable options that control the behavior of the SQP solver.

    module sqpopt_options_module

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_linalg_module, only: sqpopt_linsolve_lusol

    implicit none

    private

    integer, parameter, public :: sqpopt_hessian_bfgs   = 1  !! limited-memory (damped) BFGS quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_sr1    = 2  !! limited-memory symmetric rank-1 (SR1) quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_exact  = 3  !! user-supplied exact sparse Hessian of the Lagrangian

    type, public :: sqpopt_options_type
        !! options that control the SQP algorithm.

        integer  :: max_iter          = 100      !! maximum number of major SQP iterations
        integer  :: print_level       = 0        !! amount of diagnostic printing (0 = silent)
        integer  :: hessian_mode      = sqpopt_hessian_bfgs !! Hessian approximation strategy to use
        integer  :: lbfgs_memory      = 10        !! number of `(s,y)` vector pairs retained by the limited-memory
                                                   !! Hessian approximation (independent of the problem size `n`)
        integer  :: linear_solver_mode = sqpopt_linsolve_lusol !! sparse linear solver used for the QP subproblem
                                                                !! (see [[sqpopt_linalg_module]])

        real(wp) :: ftol  = 1.0e-8_wp   !! convergence tolerance on relative change in the objective function
        real(wp) :: xtol  = 1.0e-8_wp   !! convergence tolerance on relative change in the optimization variables
        real(wp) :: ctol  = 1.0e-8_wp   !! feasibility tolerance on the constraint violation
        real(wp) :: ktol  = 1.0e-6_wp   !! tolerance on the KKT optimality conditions

    end type sqpopt_options_type

    end module sqpopt_options_module
!*******************************************************************************
