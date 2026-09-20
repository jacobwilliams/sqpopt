!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  User-settable options that control the behavior of the SQP solver.

    module sqpopt_options_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    integer, parameter, public :: sqpopt_hessian_bfgs   = 1  !! damped BFGS quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_sr1    = 2  !! symmetric rank-1 (SR1) quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_exact  = 3  !! user-supplied exact Hessian of the Lagrangian

    type, public :: sqpopt_options_type
        !! options that control the SQP algorithm.

        integer  :: max_iter          = 100      !! maximum number of major SQP iterations
        integer  :: print_level       = 0        !! amount of diagnostic printing (0 = silent)
        integer  :: hessian_mode      = sqpopt_hessian_bfgs !! Hessian approximation strategy to use

        real(wp) :: ftol  = 1.0e-8_wp   !! convergence tolerance on relative change in the objective function
        real(wp) :: xtol  = 1.0e-8_wp   !! convergence tolerance on relative change in the optimization variables
        real(wp) :: ctol  = 1.0e-8_wp   !! feasibility tolerance on the constraint violation
        real(wp) :: ktol  = 1.0e-6_wp   !! tolerance on the KKT optimality conditions

        logical  :: use_sparse = .false. !! if true, use sparse linear algebra for the Jacobian/Hessian

    end type sqpopt_options_type

    end module sqpopt_options_module
!*******************************************************************************
