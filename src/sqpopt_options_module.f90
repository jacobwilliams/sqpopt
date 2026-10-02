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
    use sqpopt_hessian_module,    only: sqpopt_hessian_bfgs
    use sqpopt_problem_module,    only: sqpopt_derivatives_accurate
    use, intrinsic :: iso_fortran_env, only: output_unit

    implicit none

    private

    type, public :: sqpopt_options_type
        !! options that control the SQP algorithm.

        integer  :: max_iter          = 100      !! maximum number of major SQP iterations
        integer  :: print_level       = 0        !! amount of printing to `output_unit`: `0` = none; `1` = the
                                                   !! problem and method, one line per major iteration (objective,
                                                   !! infeasibility, KKT error, step length, evaluations, and flags
                                                   !! for the events of the iteration), and a final summary; `2` =
                                                   !! also the step norm, QP iterations, evaluations per iteration,
                                                   !! largest multiplier, unscaled stationarity, the globalization's
                                                   !! state, and the Hessian's; `3` = also the details of each
                                                   !! iteration (QP solves, every line-search or trust-region trial,
                                                   !! restoration phases, Hessian resets), the constraint scale
                                                   !! factors, and the solution (variables and constraints, with
                                                   !! their bounds, multipliers, and which are active)
        integer  :: output_unit       = output_unit !! Fortran unit for the printed output (default: standard output)
        integer  :: diagnostic_level  = 0        !! diagnostics of the solve, to help when it doesn't converge (see
                                                   !! [[sqpopt_diagnostics_module]]): `0` = none (no cost); `1` = a
                                                   !! diagnosis of the final point (the constraints and variables
                                                   !! with the largest errors and multipliers, the active set and
                                                   !! its degeneracy, the last QP, and what the final status most
                                                   !! likely means); `2` = also a report on the starting point (the
                                                   !! problem's scaling and structure), and what the iterations
                                                   !! showed (the rate of convergence, constraints whose changes
                                                   !! disagree with their derivatives, a working set that keeps
                                                   !! flipping, where the time went), which costs work of the order
                                                   !! of `n + m` and the Jacobian's nonzeros per iteration; `3` =
                                                   !! also the diagnostics that call the user's functions (two
                                                   !! calls of `fc` next to the starting point). Levels `0` to `2`
                                                   !! call the user's functions exactly as often as level `0`, and
                                                   !! no level changes the iterates. The diagnostics are returned
                                                   !! in `results%diagnosis`, and printed to `output_unit` if
                                                   !! `print_level >= 1` (with `print_level = 0`, nothing is
                                                   !! printed at any diagnostic level)
        integer  :: diagnostics_unit  = -1       !! with `diagnostic_level >= 2`: a Fortran unit (open for
                                                   !! formatted writing) for the history of the iterations, one
                                                   !! line per major iteration with comma-separated values and a
                                                   !! header line (`-1`: none). It is data for plotting and
                                                   !! post-processing, not part of the log, so it is written
                                                   !! whatever `print_level` is
        integer  :: hessian_mode      = sqpopt_hessian_bfgs !! Hessian approximation strategy to use
        logical  :: inertia_control   = .false.   !! with `hessian_mode = sqpopt_hessian_exact` or
                                                   !! `sqpopt_hessian_sr1` (the Hessians that can be indefinite):
                                                   !! find the shift \( \delta \) that makes \( H + \delta I \)
                                                   !! convex on the QP's working set from the inertia of the KKT
                                                   !! matrix, by a sparse factorization (see
                                                   !! [[sqpopt_inertia_module]]). Without it, the exact Hessian's
                                                   !! shift comes from the QP solver's tests alone, and SR1 is not
                                                   !! corrected. It needs a library built with MUMPS (the
                                                   !! `HAS_MUMPS` preprocessor directive, see `sqpopt_has_mumps`),
                                                   !! and is invalid without it. It is not used with BFGS, which is
                                                   !! positive definite
        logical  :: direct_qp         = .false.   !! first try to solve each QP subproblem directly, by sparse
                                                   !! factorizations of the KKT matrix of its working set, starting
                                                   !! from the working set of the previous QP (see
                                                   !! [[sqpopt_qp_direct_module]]); the active-set QP solver is
                                                   !! only run if that doesn't give the solution after a few
                                                   !! changes of the working set (`qp_solver%direct_max_changes`).
                                                   !! Meant for large problems, where the active-set solvers take
                                                   !! most of the time. With the exact Hessian, use it with
                                                   !! `inertia_control`. With a quasi-Newton Hessian, each
                                                   !! factorization costs two more solves per stored pair, so the
                                                   !! automatic memory is short with it (see `lbfgs_memory`). It
                                                   !! needs a library built with MUMPS, and is invalid without it
        logical  :: direct_least_squares = .false. !! compute the Gauss-Newton restoration steps and the
                                                   !! second-order corrections by a sparse factorization instead
                                                   !! of the iterative `LSQR` (see
                                                   !! [[sqpopt_least_squares_module]]), with any Hessian mode. It
                                                   !! pays on large problems that take such steps and whose
                                                   !! constraints are coupled, where `LSQR` needs many iterations
                                                   !! (a chain of 100,000 circle constraints: 88.7 s with `LSQR`,
                                                   !! 1.1 s with this); on small problems it changes little. It
                                                   !! needs a library built with MUMPS, and is invalid without it
        integer  :: factorization_threads = 1     !! number of OpenMP threads the sparse factorizations use
                                                   !! (`inertia_control`, `direct_qp`, and `direct_least_squares`):
                                                   !! `1` (the default) for none, a larger number for that many,
                                                   !! or `0` to leave it to the OpenMP environment
                                                   !! (`OMP_NUM_THREADS`, or every core). It needs MUMPS and its
                                                   !! BLAS to be built with OpenMP (conda-forge's are). It only
                                                   !! pays on large problems whose factors are dense enough: a
                                                   !! 3-D grid matrix of order 216,000 was factored 2.3 times
                                                   !! faster on 4 threads, but banded problems and a 2-D grid
                                                   !! gained nothing, and on small problems threads cost a lot
                                                   !! (see [[sqpopt_symmetric_solver_module]])
        integer  :: lbfgs_memory      = 0         !! number of `(s,y)` vector pairs retained by the limited-memory
                                                   !! Hessian approximation. `0` (the default) picks it from the problem
                                                   !! size `n`: \( \max(10, \min(n, 100)) \) -- more pairs help the
                                                   !! larger problems (on a quadratic, BFGS with `n` pairs converges in
                                                   !! about `n` iterations), but more pairs than variables keep stale
                                                   !! curvature and slow the small ones down. With `direct_qp`, `0`
                                                   !! picks 10 instead: each factorization then costs two solves per
                                                   !! pair, and work that grows with the square of their number (on
                                                   !! a control problem with 10,001 variables, 4.9 s with 10 pairs
                                                   !! and 18 s with 100, in about the same number of iterations)
        integer  :: qp_solver_mode = sqpopt_qp_auto !! QP subproblem algorithm to use
                                                     !! (see [[sqpopt_qp_solver_module]])
        integer  :: linesearch_mode = sqpopt_linesearch_filter !! line search strategy to use
                                                                !! (see [[sqpopt_linesearch_module]])
        integer  :: merit_mode = sqpopt_merit_l1 !! merit function to use
                                                  !! (see [[sqpopt_merit_module]])
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
        integer  :: derivative_accuracy = sqpopt_derivatives_accurate !! the accuracy of the derivatives the solver
                                                   !! asks the user's `gjac` for (its `accuracy` argument): with
                                                   !! `sqpopt_derivatives_fast`, it asks for fast (e.g.
                                                   !! forward-difference) derivatives until it is near a
                                                   !! solution, then for accurate ones (see
                                                   !! [[sqpopt_iterate_module]] and `derivative_switch_tol`); with
                                                   !! `sqpopt_derivatives_accurate` (the default), always for
                                                   !! accurate ones
        real(wp) :: derivative_switch_tol = 1.0e-5_wp !! with `derivative_accuracy = sqpopt_derivatives_fast`: the
                                                   !! KKT and feasibility errors (of the scaled problem) below which
                                                   !! the solver is near enough a solution to switch to accurate
                                                   !! derivatives. It also switches when a step fails, when progress
                                                   !! stalls, and before stopping (at a solution, or as infeasible),
                                                   !! so this is only a shortcut: switching earlier costs accurate
                                                   !! derivatives on iterations that don't need them (on the HS
                                                   !! suite, anything from about `10*ktol` down did equally well)
        real(wp) :: elastic_multiplier_limit = 30.0_wp !! when a constraint's multiplier diverges (its push
                                                     !! \( |\lambda_i| \lVert \nabla c_i \rVert_\infty \) exceeds this
                                                     !! \( \times \max(1,\lVert g \rVert_\infty) \) and keeps growing,
                                                     !! a sign that the constraint qualification fails nearby), the QP
                                                     !! is re-solved with that constraint elastic, so the iterates can
                                                     !! leave (see [[sqpopt_iterate_module]]); `0` disables this
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

        real(wp) :: dual_inf_tol = 1.0_wp  !! tolerance on the stationarity residual of the *unscaled* problem
                                           !! (as IPOPT's `dual_inf_tol`): with automatic scaling, `ktol` applies
                                           !! to the scaled problem, and an objective scaled far down (e.g. at a
                                           !! poor starting point) would make it very loose in the original
                                           !! units; convergence also requires this
        real(wp) :: acceptable_ktol = 1.0e-4_wp !! looser "acceptable" KKT tolerance (as in IPOPT): if the KKT test
                                                !! with `acceptable_ktol`/`acceptable_ctol` holds for `acceptable_iter`
                                                !! consecutive iterations (but the normal one doesn't), stop with
                                                !! `istat=sqpopt_acceptable`
        real(wp) :: acceptable_ctol = 1.0e-6_wp !! looser "acceptable" feasibility tolerance (see `acceptable_ktol`)
        integer  :: acceptable_iter = 15        !! consecutive acceptable iterations needed (`0` disables the test)
        real(wp) :: acceptable_obj_change_tol = 1.0e-3_wp !! an iteration only counts as acceptable if the objective
                                                !! changed by at most this from the previous one, relative to
                                                !! `max(1, |f|)` (of the original problem; as IPOPT's
                                                !! `acceptable_obj_change_tol`, which is off by default there): the
                                                !! solver doesn't stop at the acceptable level while the objective
                                                !! is still changing fast, which happens when the problem's scaling
                                                !! makes `acceptable_ktol` loose. The Hock-Schittkowski results are
                                                !! the same for any value from `1e-4` to `0.1`; a large value
                                                !! (`1e20`) turns the test off
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
        real(wp) :: scaling_min_value = 1.0e-8_wp   !! the smallest scale factor of `scaling` (as IPOPT's
                                                    !! `nlp_scaling_min_value`): a function whose gradient at the
                                                    !! starting point is larger than `scaling_max_gradient` divided
                                                    !! by this is scaled by this factor, and no further (`0` = no
                                                    !! limit)

        real(wp) :: hessian_scale0 = 1.0_wp     !! initial Hessian approximation \( B_0 = \) `hessian_scale0` \( \times I \)
                                                !! (replaced by the usual quasi-Newton scaling after the first update)

    end type sqpopt_options_type

    end module sqpopt_options_module
!*******************************************************************************
