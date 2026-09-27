!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  User-settable options that control the behavior of the SQP solver.

    module sqpopt_options_module

    use sqpopt_kinds,             only: wp => sqpopt_module_wp
    use sqpopt_linesearch_module, only: sqpopt_linesearch_filter, sqpopt_merit_l1, sqpopt_penalty_multipliers
    use sqpopt_qp_solver_module,  only: sqpopt_qp_auto
    use sqpopt_restoration_module, only: sqpopt_restoration_phase
    use sqpopt_types_module,      only: sqpopt_infinity
    use, intrinsic :: iso_fortran_env, only: output_unit

    implicit none

    private

    integer, parameter, public :: sqpopt_hessian_bfgs   = 1  !! limited-memory (damped) BFGS quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_sr1    = 2  !! limited-memory symmetric rank-1 (SR1) quasi-Newton Hessian approximation
    integer, parameter, public :: sqpopt_hessian_exact  = 3  !! user-supplied exact sparse Hessian of the Lagrangian (the
                                                             !! `hess` function of [[set_functions]], with the pattern of
                                                             !! [[set_hessian_sparsity]]; see [[sqpopt_hessian_module]])

    type, public :: sqpopt_options_type
        !! options that control the SQP algorithm.

        integer  :: max_iter          = 100      !! maximum number of major SQP iterations
        integer  :: print_level       = 0        !! amount of printing to `output_unit`: `0` = none, `1` = one line
                                                   !! per major iteration plus a final summary, `2` = also the
                                                   !! penalty parameter, step norm, and QP iterations
        integer  :: output_unit       = output_unit !! Fortran unit for the printed output (default: standard output)
        integer  :: hessian_mode      = sqpopt_hessian_bfgs !! Hessian approximation strategy to use
        integer  :: lbfgs_memory      = 10        !! number of `(s,y)` vector pairs retained by the limited-memory
                                                   !! Hessian approximation (independent of the problem size `n`)
        integer  :: qp_solver_mode = sqpopt_qp_auto !! QP subproblem algorithm to use
                                                     !! (see [[sqpopt_qp_solver_module]])
        integer  :: linesearch_mode = sqpopt_linesearch_filter !! line search strategy to use
                                                                !! (see [[sqpopt_linesearch_module]])
        integer  :: merit_mode = sqpopt_merit_l1 !! merit function to use
                                                  !! (see [[sqpopt_linesearch_module]])
        integer  :: penalty_update = sqpopt_penalty_multipliers !! how the merit function's penalty parameter is
                                                                 !! updated (see [[update_penalty_parameter]])
        integer  :: restoration_mode = sqpopt_restoration_phase !! feasibility restoration strategy (see
                                                                 !! [[sqpopt_restoration_module]]): a restoration phase
                                                                 !! (`sqpopt_restoration_phase`), or a single Gauss-Newton
                                                                 !! step each time (`sqpopt_restoration_gauss_newton`)
        real(wp) :: restoration_exit_factor = 0.9_wp !! a restoration phase ends once the \( \ell_1 \) violation is below this
                                                     !! factor times its value where the phase started, and the point is
                                                     !! acceptable to the filter or funnel (`0 < factor < 1`)
        integer  :: restoration_max_iter = 50        !! maximum number of iterations of a restoration phase
        integer  :: max_consecutive_failures = 5 !! stop (with the failing component's status code, e.g.
                                                  !! `sqpopt_line_search_failed` or `sqpopt_qp_solve_failed`)
                                                  !! after this many consecutive major iterations in which the
                                                  !! QP solve or the line search/trust-region step failed

        real(wp) :: ftol  = 1.0e-8_wp   !! secondary "stalled progress" stopping criterion: once feasible,
                                        !! also stop if the objective's relative change from the previous
                                        !! iterate is below `ftol` *and* the variables' relative change is
                                        !! below `xtol` (a safeguard against looping to `max_iter` on
                                        !! marginal steps when the KKT test never quite reaches `ktol`)
        real(wp) :: xtol  = 1.0e-8_wp   !! see `ftol`
        real(wp) :: ctol  = 1.0e-8_wp   !! feasibility tolerance on the constraint violation
        real(wp) :: ktol  = 1.0e-6_wp   !! tolerance on the KKT optimality conditions

        real(wp) :: acceptable_ktol = 1.0e-4_wp !! looser "acceptable" KKT tolerance (as in IPOPT): if the KKT test
                                                !! with `acceptable_ktol`/`acceptable_ctol` holds for `acceptable_iter`
                                                !! consecutive iterations (but the normal one doesn't), stop with
                                                !! `istat=sqpopt_acceptable`
        real(wp) :: acceptable_ctol = 1.0e-6_wp !! looser "acceptable" feasibility tolerance (see `acceptable_ktol`)
        integer  :: acceptable_iter = 15        !! consecutive acceptable iterations needed (`0` disables the test)
        integer  :: stall_iter = 3              !! the stalled-progress test (`ftol`/`xtol`) must hold for this many
                                                !! consecutive iterations before the solver stops with
                                                !! `sqpopt_stalled` (a single negligible step isn't a stall)

        integer  :: max_evals       = 0         !! stop (`istat=sqpopt_max_evals_reached`) after this many calls of the
                                                !! objective and constraint function `fc` (`0` = no limit)
        real(wp) :: max_time        = 0.0_wp    !! stop (`istat=sqpopt_time_limit_reached`) after this much wall-clock
                                                !! time, in seconds (`0` = no limit)
        real(wp) :: obj_lower_limit = -sqpopt_infinity !! stop (`istat=sqpopt_unbounded`) if the objective falls below
                                                       !! this at a feasible point (`-sqpopt_infinity` = no limit)

        logical  :: scaling = .true.            !! gradient-based scaling of the objective and constraints (as in
                                                !! IPOPT): at the starting point, each of them whose gradient has an
                                                !! element larger than `scaling_max_gradient` is scaled down so that
                                                !! its largest element is `scaling_max_gradient`. The tolerances
                                                !! `ktol`/`ctol` apply to the scaled problem; all results are
                                                !! reported for the original one
        real(wp) :: scaling_max_gradient = 100.0_wp !! see `scaling`

        real(wp) :: hessian_scale0 = 1.0_wp     !! initial Hessian approximation \( B_0 = \) `hessian_scale0` \( \times I \)
                                                !! (replaced by the usual quasi-Newton scaling after the first update)

    end type sqpopt_options_type

    end module sqpopt_options_module
!*******************************************************************************
