!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  User-settable options that control the behavior of the SQP solver.

    module sqpopt_options_module

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_linesearch_module, only: sqpopt_linesearch_armijo, sqpopt_merit_l1
    use sqpopt_qp_solver_module,  only: sqpopt_qp_composite

    implicit none

    private

    integer, parameter, public :: sqpopt_hessian_bfgs   = 1  !! limited-memory (damped) BFGS quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_sr1    = 2  !! limited-memory symmetric rank-1 (SR1) quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_exact  = 3  !! user-supplied exact sparse Hessian of the Lagrangian

    type, public :: sqpopt_options_type
        !! options that control the SQP algorithm.

        integer  :: max_iter          = 100      !! maximum number of major SQP iterations
        integer  :: print_level       = 0        !! amount of diagnostic printing to `stdout` per major iteration
                                                   !! (`0` = silent, `>=1` = one summary line per iteration)
        integer  :: hessian_mode      = sqpopt_hessian_bfgs !! Hessian approximation strategy to use
        integer  :: lbfgs_memory      = 10        !! number of `(s,y)` vector pairs retained by the limited-memory
                                                   !! Hessian approximation (independent of the problem size `n`)
        integer  :: qp_solver_mode = sqpopt_qp_composite !! QP subproblem algorithm to use
                                                          !! (see [[sqpopt_qp_solver_module]])
        integer  :: linesearch_mode = sqpopt_linesearch_armijo !! line search strategy to use
                                                                !! (see [[sqpopt_linesearch_module]])
        integer  :: merit_mode = sqpopt_merit_l1 !! merit function to use
                                                  !! (see [[sqpopt_linesearch_module]])

        real(wp) :: ftol  = 1.0e-8_wp   !! secondary "stalled progress" stopping criterion: once feasible,
                                        !! also stop if the objective's relative change from the previous
                                        !! iterate is below `ftol` *and* the variables' relative change is
                                        !! below `xtol` (a safeguard against looping to `max_iter` on
                                        !! marginal steps when the KKT test never quite reaches `ktol`)
        real(wp) :: xtol  = 1.0e-8_wp   !! see `ftol`
        real(wp) :: ctol  = 1.0e-8_wp   !! feasibility tolerance on the constraint violation
        real(wp) :: ktol  = 1.0e-6_wp   !! tolerance on the KKT optimality conditions

    end type sqpopt_options_type

    end module sqpopt_options_module
!*******************************************************************************
