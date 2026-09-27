program test_hs_suite

    !! Benchmark and regression test on the Schittkowski/Hock-Schittkowski
    !! test problem collection (305 problems with validated optimal
    !! solutions; see [[hs_problems_module]]).
    !!
    !! Each problem is solved with the default options (and `max_iter=1000`)
    !! from its standard starting point, and classified as:
    !!
    !! * **solved**: feasible (violation <= `feas_tol`), with an objective
    !!   value within `rel_tol` (relative) of the validated optimum, or better;
    !! * **local**: feasible and converged (`sqpopt_success`, `_stalled`, or
    !!   `_acceptable`), but at a worse objective value -- a different local
    !!   solution;
    !! * **failed**: anything else.
    !!
    !! Each problem's derivatives are first checked against finite
    !! differences (see [[check_derivatives]]); where they don't match, this
    !! harness supplies finite-difference derivatives instead (marked `fd` in
    !! the table). That is now only the 16 problems that have no analytic
    !! derivatives at all (TP332, 348, 356, 357, 362, 364, 365, 366, 369,
    !! 370, 371, 377, 378, 391, 392, 393): the collection intends them to be
    !! solved with numerical derivatives. The analytic derivatives that were
    !! wrong in the original collection (18 problems) are fixed in
    !! [[schittkowski_problems_module]].
    !!
    !! Evaluation counts are shown next to those of Schittkowski's NLPQLP on
    !! the same problems (from the collection's `TEST.DAT`; its gradients
    !! were forward-difference approximations, each costing `n` extra
    !! function evaluations, which are not included in its counts here).
    !!
    !! **Regression test:** the test fails if any problem outside
    !! `known_unsolved` is not solved. (If a problem in that list becomes
    !! solved, it is reported, so that the list can be updated.)
    !!
    !! **Report:** the results are also written as a Markdown report (summary,
    !! efficiency vs. NLPQLP, termination statuses, and a per-problem table)
    !! to `test/hs_suite_results.md`, or to the file given as a (positional)
    !! command-line argument:
    !!
    !!    fpm test test_hs_suite --profile release -- results.md
    !!
    !! **Configuration:** the globalization can be changed from the command
    !! line, for comparisons (see `tools/hs_performance_table.sh`):
    !!
    !! * `--linesearch=filter|funnel|armijo|watchdog|exact` (`options%linesearch_mode`)
    !! * `--restoration=phase|gauss-newton` (`options%restoration_mode`)
    !! * `--hessian=bfgs|sr1|exact` (`options%hessian_mode`; the collection
    !!   has no second derivatives, so `exact` uses the Hessian of the
    !!   Lagrangian computed by central differences of the analytic gradient
    !!   and Jacobian, see [[hess_fd]], except on the problems that use
    !!   finite-difference first derivatives, which keep BFGS)
    !! * `--trust-region` (`trust_region%enabled = .true.`: the trust-region
    !!   globalization, with the filter or funnel test in those modes, else
    !!   the merit-function ratio test)
    !! * `--merit=l1|al` (`options%merit_mode`)
    !! * `--penalty=multipliers|model` (`options%penalty_update`)
    !! * `--no-interpolate` (`linesearch%interpolate = .false.`)
    !! * `--nonmonotone=N` (`linesearch%nonmonotone_len = N`)
    !! * `--qp=auto|dense|sparse|sparse-lsqr` (`options%qp_solver_mode`; the
    !!   HS problems are small, so `auto` picks the dense QP for all of them:
    !!   `sparse` runs them through the sparse QP instead, and `sparse-lsqr`
    !!   through its `LSQR` null-space method)
    !!
    !! **Web data:** `--web-data=FILE` also writes the results as a JavaScript
    !! data file for the interactive results page of the user guide
    !! (`web/hs_results.html`, which reads `web/js/hs_results_data.js`; see
    !! [[write_web_data]]). It doesn't change the configuration.
    !!
    !! For debugging, `--problem=N` solves only problem `TPN`, and
    !! `--print=L` sets `options%print_level = L`.
    !!
    !! If any of these options is given, the regression test is skipped (the
    !! baseline is for the defaults). The last line printed before the
    !! regression test is a machine-readable summary:
    !! `summary: solved=... local=... failed=... nf=... ng=...`.

    use, intrinsic :: ieee_arithmetic, only: ieee_set_halting_mode, ieee_all, ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, compiler_version
    use hs_problems_module
    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type, sqpopt_hessian_bfgs, sqpopt_hessian_sr1, sqpopt_hessian_exact
    use sqpopt_qp_solver_module,  only: sqpopt_qp_solver_type, sqpopt_qp_auto, sqpopt_qp_dense, &
                                        sqpopt_qp_reduced_hessian
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_null_space_lu, sqpopt_null_space_lsqr
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type, sqpopt_linesearch_armijo, sqpopt_linesearch_exact, &
                                        sqpopt_linesearch_watchdog, sqpopt_linesearch_filter, sqpopt_merit_l1, &
                                        sqpopt_linesearch_funnel, &
                                        sqpopt_merit_augmented_lagrangian, sqpopt_penalty_multipliers, &
                                        sqpopt_penalty_model
    use sqpopt_types_module,   only: sqpopt_results_type, sqpopt_success, sqpopt_stalled, sqpopt_acceptable, &
                                     sqpopt_status_message
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_restoration_module,  only: sqpopt_restoration_phase, sqpopt_restoration_gauss_newton
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    real(dp), parameter :: rel_tol  = 1.0e-4_dp  !! relative objective tolerance for "solved" (4 significant digits)
    real(dp), parameter :: feas_tol = 1.0e-6_dp  !! constraint/bound violation tolerance for "solved"/"local"

    !> problems not (yet) solved by `sqpopt` with the default options -- the
    !! regression baseline (see the program documentation). As of 2026-09-26:
    !! 276 of the 305 problems solved, 29 local solutions, 0 failures.
    integer, dimension(*), parameter :: known_unsolved = [ &
          2,  16,  25,  33,  38,  54,  55,  57,  59,  87,  97,  98, 105, 109, 202, &
        213, 236, 239, 265, 272, 283, 287, 304, 305, 312, 327, 338, 340, 362 ]

    type :: problem_context
        !! the user data passed to the problem functions
        integer :: id = 0
        integer :: n = 0, m = 0
        logical :: fd_g = .false., fd_jac = .false. !! use central differences for the gradient/Jacobian
    end type problem_context

    type :: run_record
        !! the outcome of one problem, for the report
        integer  :: id = 0, n = 0, m = 0, me = 0
        integer  :: istat = 0, iterations = 0
        integer  :: nf = 0, ng = 0                    !! evaluation counts: calls of `fc` and of `gjac`
        integer  :: q_nf = 0, q_ndf = 0               !! NLPQLP's evaluation counts
        real(dp) :: f = 0.0_dp, f_star = 0.0_dp, rel = 0.0_dp, viol = 0.0_dp, kkt = 0.0_dp
        character(len=6) :: outcome = ''
        logical  :: fd = .false.          !! central-difference derivatives were used
        logical  :: regression = .false. !! not solved, but not in `known_unsolved`
    end type run_record

    character(len=*), parameter :: default_report_file = 'test/hs_suite_results.md'

    type(problem_context), target :: ctx
    type(run_record), dimension(hs_n_problems) :: rec
    type(hs_problem) :: p
    integer :: k, n_solved, n_local, n_failed, n_fd, n_regressions, n_improved
    integer :: sum_nf, sum_ng, sum_nlpqlp_nf, sum_nlpqlp_ndf
    character(len=6) :: outcome
    integer(int64) :: t0, t1, rate
    character(len=:), allocatable :: report_file
    character(len=:), allocatable :: web_data_file !! `--web-data=FILE` (empty: none)
    character(len=:), allocatable :: cfg_string    !! the configuration options given, for the web data

    ! the configuration (see the program documentation):
    integer :: cfg_linesearch  = sqpopt_linesearch_filter
    integer :: cfg_merit       = sqpopt_merit_l1
    integer :: cfg_penalty     = sqpopt_penalty_multipliers
    logical :: cfg_interpolate = .true.
    integer :: cfg_nonmonotone = 0
    logical :: cfg_trust_region = .false. !! `--trust-region`
    integer :: cfg_restoration  = sqpopt_restoration_phase !! `--restoration=`
    integer :: cfg_hessian      = sqpopt_hessian_bfgs      !! `--hessian=`
    integer :: cfg_problem     = 0  !! `--problem=N`: solve only this problem (`0` = all)
    integer :: cfg_print       = 0  !! `--print=L`: `options%print_level`
    integer :: cfg_qp          = sqpopt_qp_auto
    integer :: cfg_null_space  = sqpopt_null_space_lu
    logical :: cfg_default     = .true.   !! whether every setting is the default (then the regression test runs)

    call ieee_set_halting_mode(ieee_all, .false.)  ! (trial points may produce NaN/Inf, which the solver handles)
    call parse_arguments()

    write(*,*) '----------------------------'
    write(*,*) 'test_hs_suite'
    write(*,*) '----------------------------'
    write(*,'(A5,3A4,A4,A7,2A17,2A10,2A6,2A7,A4)') 'TP', 'n', 'm', 'me', 'st', 'result', 'f', 'f*', &
        'rel.err', 'viol', 'nf', 'ng', 'Q:nf', 'Q:ndf', 'fd'

    n_solved = 0; n_local = 0; n_failed = 0; n_fd = 0; n_regressions = 0; n_improved = 0
    sum_nf = 0; sum_ng = 0; sum_nlpqlp_nf = 0; sum_nlpqlp_ndf = 0
    call system_clock(t0, rate)

    do k = 1, hs_n_problems
        if (cfg_problem /= 0 .and. hs_problem_ids(k) /= cfg_problem) cycle
        call run_problem(hs_problem_ids(k), k)
    end do

    call system_clock(t1)
    write(*,'(A)') ''
    write(*,'(A,I0)') 'problems:                 ', hs_n_problems
    write(*,'(A,I0)') 'solved:                   ', n_solved
    write(*,'(A,I0)') 'local solutions:          ', n_local
    write(*,'(A,I0)') 'failed:                   ', n_failed
    write(*,'(A,I0)') 'with FD derivatives:      ', n_fd
    write(*,'(A,2(I0,A))') 'evaluations (solved problems): sqpopt ', sum_nf, ' f, ', sum_ng, ' g'
    write(*,'(A,2(I0,A))') '                               NLPQLP ', sum_nlpqlp_nf, ' f, ', sum_nlpqlp_ndf, ' g'
    write(*,'(A,F0.2,A)') 'time: ', real(t1-t0, dp)/real(rate, dp), ' s'

    if (cfg_problem == 0) then
        call write_report(report_file, real(t1-t0, dp)/real(rate, dp))
        write(*,'(A)') 'report written to: '//report_file
        if (len(web_data_file) > 0) then
            call write_web_data(web_data_file)
            write(*,'(A)') 'web data written to: '//web_data_file
        end if
    end if
    write(*,'(5(A,I0))') 'summary: solved=', n_solved, ' local=', n_local, ' failed=', n_failed, &
                         ' nf=', sum_nf, ' ng=', sum_ng

    if (.not. cfg_default) then
        print '(A)', 'test_hs_suite: non-default configuration, regression test skipped'
        stop
    end if

    if (n_improved > 0) then
        write(*,'(I0,A)') n_improved, ' problem(s) in known_unsolved are now solved: update the list'
        write(*,'(A,*(1X,I0))') '   now solved: TP', pack(rec%id, rec%outcome == 'solved' .and. &
                                                       [(any(known_unsolved == rec(k)%id), k=1,hs_n_problems)])
    end if
    if (n_regressions > 0) then
        write(*,'(I0,A)') n_regressions, ' problem(s) not in known_unsolved were not solved'
        do k = 1, hs_n_problems
            if (rec(k)%regression) write(*,'(A,I0,A,A,A,I0,A,A)') '   TP', rec(k)%id, ': ', trim(rec(k)%outcome), &
                ' (istat=', rec(k)%istat, ') ', sqpopt_status_message(rec(k)%istat)
        end do
        error stop 'test_hs_suite FAILED'
    end if
    print '(A)', 'test_hs_suite PASSED'

    contains

    subroutine parse_arguments()
    !! the report file (positional) and the configuration options (see the
    !! program documentation)
    integer :: i, n, ios
    character(len=256) :: arg
    report_file = default_report_file
    web_data_file = ''
    cfg_string = ''
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        if (arg(1:2) /= '--') then
            report_file = trim(arg)
            cycle
        end if
        if (arg(1:11) == '--web-data=') then
            web_data_file = trim(arg(12:))
            if (len(web_data_file) == 0) error stop 'test_hs_suite: bad --web-data value'
            cycle
        end if
        cfg_default = .false.
        cfg_string = trim(cfg_string//' '//trim(arg))
        select case (trim(arg))
        case ('--linesearch=filter');   cfg_linesearch = sqpopt_linesearch_filter
        case ('--linesearch=funnel');   cfg_linesearch = sqpopt_linesearch_funnel
        case ('--linesearch=armijo');   cfg_linesearch = sqpopt_linesearch_armijo
        case ('--linesearch=watchdog'); cfg_linesearch = sqpopt_linesearch_watchdog
        case ('--linesearch=exact');    cfg_linesearch = sqpopt_linesearch_exact
        case ('--merit=l1');            cfg_merit = sqpopt_merit_l1
        case ('--merit=al');            cfg_merit = sqpopt_merit_augmented_lagrangian
        case ('--penalty=multipliers'); cfg_penalty = sqpopt_penalty_multipliers
        case ('--penalty=model');       cfg_penalty = sqpopt_penalty_model
        case ('--no-interpolate');      cfg_interpolate = .false.
        case ('--trust-region');        cfg_trust_region = .true.
        case ('--restoration=phase');   cfg_restoration = sqpopt_restoration_phase
        case ('--restoration=gauss-newton'); cfg_restoration = sqpopt_restoration_gauss_newton
        case ('--hessian=bfgs');        cfg_hessian = sqpopt_hessian_bfgs
        case ('--hessian=sr1');         cfg_hessian = sqpopt_hessian_sr1
        case ('--hessian=exact');       cfg_hessian = sqpopt_hessian_exact
        case ('--qp=auto');             cfg_qp = sqpopt_qp_auto
        case ('--qp=dense');            cfg_qp = sqpopt_qp_dense
        case ('--qp=sparse');           cfg_qp = sqpopt_qp_reduced_hessian
        case ('--qp=sparse-lsqr');      cfg_qp = sqpopt_qp_reduced_hessian; cfg_null_space = sqpopt_null_space_lsqr
        case default
            if (arg(1:14) == '--nonmonotone=') then
                read(arg(15:), *, iostat=ios) n
                if (ios /= 0 .or. n < 0) error stop 'test_hs_suite: bad --nonmonotone value'
                cfg_nonmonotone = n
            else if (arg(1:10) == '--problem=') then
                read(arg(11:), *, iostat=ios) n
                if (ios /= 0 .or. .not. any(hs_problem_ids == n)) error stop 'test_hs_suite: bad --problem value'
                cfg_problem = n
            else if (arg(1:8) == '--print=') then
                read(arg(9:), *, iostat=ios) n
                if (ios /= 0) error stop 'test_hs_suite: bad --print value'
                cfg_print = n
            else
                write(*,'(A)') 'test_hs_suite: unknown option: '//trim(arg)
                error stop 1
            end if
        end select
    end do
    end subroutine parse_arguments

    subroutine run_problem(id, k)
    !! solve problem `id` (the `k`-th in the collection) and report the outcome
    integer, intent(in) :: id, k

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_linesearch_type) :: linesearch
    type(sqpopt_qp_solver_type) :: qp_solver
    type(sqpopt_trust_region_type) :: trust_region
    type(sqpopt_results_type) :: r
    integer, dimension(:), allocatable :: irow, icol, hrow, hcol
    integer  :: i, j, nnz, istat
    real(dp) :: rel, viol
    logical  :: feasible, converged, expected_unsolved

    call hs_setup(id, p)
    ctx = problem_context(id=id, n=p%n, m=p%m)
    call check_derivatives(ctx)
    if (ctx%fd_g .or. ctx%fd_jac) n_fd = n_fd + 1

    ! dense Jacobian pattern (row by row):
    nnz = p%m*p%n
    allocate(irow(nnz), icol(nnz))
    do i = 1, p%m
        do j = 1, p%n
            irow((i-1)*p%n+j) = i
            icol((i-1)*p%n+j) = j
        end do
    end do

    call problem%set_problem_size(n=p%n, m_eq=p%me, m_ineq=p%m-p%me)
    call problem%set_bounds(real(p%x_lb, wp), real(p%x_ub, wp), real(p%c_lb, wp), real(p%c_ub, wp))
    call problem%set_jacobian_sparsity(nnz, irow, icol)
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv, hess=hess_fd, data=ctx)
    ! dense lower-triangle Hessian pattern (row by row):
    block
        integer :: k
        k = 0
        allocate(hrow(p%n*(p%n+1)/2), hcol(p%n*(p%n+1)/2))
        do i = 1, p%n
            do j = 1, i
                k = k + 1
                hrow(k) = i
                hcol(k) = j
            end do
        end do
    end block
    call problem%set_hessian_sparsity(size(hrow), hrow, hcol)
    options%max_iter        = 1000
    options%linesearch_mode = cfg_linesearch
    options%merit_mode      = cfg_merit
    options%penalty_update  = cfg_penalty
    options%qp_solver_mode  = cfg_qp
    options%print_level     = cfg_print
    options%restoration_mode = cfg_restoration
    options%hessian_mode    = cfg_hessian
    if (cfg_hessian == sqpopt_hessian_exact .and. (ctx%fd_g .or. ctx%fd_jac)) options%hessian_mode = sqpopt_hessian_bfgs
    qp_solver%sparse_qp%null_space = cfg_null_space
    linesearch%interpolate     = cfg_interpolate
    linesearch%nonmonotone_len = cfg_nonmonotone

    trust_region%enabled       = cfg_trust_region

    call solver%initialize(problem=problem, options=options, linesearch=linesearch, qp_solver=qp_solver, &
                           trust_region=trust_region)
    call solver%solve(real(p%x0, wp), istat)
    call solver%get_results(r)

    viol = max(r%feasibility_error, maxval(max(p%x_lb - r%x, 0.0_dp) + max(r%x - p%x_ub, 0.0_dp)))
    rel  = (r%f - p%f_star)/max(1.0_dp, abs(p%f_star))
    feasible  = viol <= feas_tol .and. ieee_is_finite(r%f)
    converged = istat == sqpopt_success .or. istat == sqpopt_stalled .or. istat == sqpopt_acceptable
    if (feasible .and. rel <= rel_tol) then
        outcome = 'solved'
        n_solved = n_solved + 1
        sum_nf = sum_nf + r%n_eval_fc
        sum_ng = sum_ng + r%n_eval_gjac
        sum_nlpqlp_nf  = sum_nlpqlp_nf  + hs_nlpqlp_nf(k)
        sum_nlpqlp_ndf = sum_nlpqlp_ndf + hs_nlpqlp_ndf(k)
    else if (feasible .and. converged) then
        outcome = 'local'
        n_local = n_local + 1
    else
        outcome = 'FAILED'
        n_failed = n_failed + 1
    end if

    expected_unsolved = any(known_unsolved == id)
    if (outcome /= 'solved' .and. .not. expected_unsolved) n_regressions = n_regressions + 1
    if (outcome == 'solved' .and. expected_unsolved) n_improved = n_improved + 1

    rec(k) = run_record(id=id, n=p%n, m=p%m, me=p%me, istat=istat, iterations=r%iterations, &
                        nf=r%n_eval_fc, ng=r%n_eval_gjac, &
                        q_nf=hs_nlpqlp_nf(k), q_ndf=hs_nlpqlp_ndf(k), &
                        f=r%f, f_star=p%f_star, rel=rel, viol=viol, kkt=r%kkt_error, outcome=outcome, &
                        fd=ctx%fd_g .or. ctx%fd_jac, regression=outcome /= 'solved' .and. .not. expected_unsolved)

    write(*,'(I5,3I4,I4,A7,2ES17.8,2ES10.2,2I6,2I7,A4,A)') id, p%n, p%m, p%me, istat, outcome, r%f, p%f_star, &
        rel, viol, r%n_eval_fc, r%n_eval_gjac, hs_nlpqlp_nf(k), hs_nlpqlp_ndf(k), &
        merge(' fd', '   ', ctx%fd_g .or. ctx%fd_jac), &
        merge('   <-- regression', '                 ', outcome /= 'solved' .and. .not. expected_unsolved)

    end subroutine run_problem

    subroutine write_report(file, total_time)
    !! write the results of all the problems (`rec`) as a Markdown report
    character(len=*), intent(in) :: file
    real(dp),         intent(in) :: total_time !! wall-clock time of the whole suite (s)

    integer :: u, i, istat, cnt(3), n_fewer, n_ratio
    integer, dimension(:), allocatable :: iters
    real(dp), dimension(:), allocatable :: ratios
    character(len=8)  :: date
    character(len=10) :: time
    character(len=*), parameter :: outcomes(3) = ['solved', 'local ', 'FAILED']

    call date_and_time(date, time)
    open(newunit=u, file=file, status='replace', action='write')

    write(u,'(A)') '# Hock-Schittkowski test suite results'
    write(u,'(A)') ''
    write(u,'(A)') 'Generated by `test/test_hs_suite.f90` on '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)// &
                   ' '//time(1:2)//':'//time(3:4)//', compiled with '//compiler_version()//'.'
    write(u,'(A)') 'Regenerate with `fpm test test_hs_suite --profile release`.'
    write(u,'(A)') ''
    write(u,'(A)') 'Each of the '//s_i(hs_n_problems)//' problems of Schittkowski''s collection is solved with the '// &
                   'default options (and `max_iter = 1000`) from its standard starting point, and classified as:'
    write(u,'(A)') ''
    write(u,'(A)') '* **solved**: feasible (violation <= '//s_r(feas_tol,'(ES8.1)')//') with an objective within '// &
                   s_r(rel_tol,'(ES8.1)')//' (relative) of the validated optimum `f*`, or better;'
    write(u,'(A)') '* **local**: feasible and converged, but at a worse objective value (a different local solution);'
    write(u,'(A)') '* **FAILED**: anything else.'
    write(u,'(A)') ''
    write(u,'(A)') 'Problems whose supplied derivatives don''t match central differences are solved with '// &
                   'central-difference derivatives instead (marked *fd*).'

    ! ---- summary ----
    write(u,'(A)') ''
    write(u,'(A)') '## Summary'
    write(u,'(A)') ''
    write(u,'(A)') '| result | problems | % |'
    write(u,'(A)') '|:--|--:|--:|'
    do i = 1, 3
        cnt(i) = count(rec%outcome == outcomes(i))
        write(u,'(A)') '| '//trim(outcomes(i))//' | '//s_i(cnt(i))//' | '// &
                       s_r(100.0_dp*cnt(i)/hs_n_problems,'(F5.1)')//' |'
    end do
    write(u,'(A)') '| **total** | **'//s_i(hs_n_problems)//'** | |'
    write(u,'(A)') ''
    write(u,'(A)') '* regressions (not solved, but not in `known_unsolved`): '//s_i(count(rec%regression))
    write(u,'(A)') '* problems using central-difference derivatives: '//s_i(count(rec%fd))
    write(u,'(A)') '* total time: '//s_r(total_time,'(F8.2)')//' s'

    ! ---- efficiency ----
    write(u,'(A)') ''
    write(u,'(A)') '## Efficiency on the solved problems'
    write(u,'(A)') ''
    write(u,'(A)') 'NLPQLP''s counts are from the collection''s `TEST.DAT`. Its gradients were forward-difference '// &
                   'approximations, whose extra function evaluations are not included in its counts.'
    write(u,'(A)') ''
    ratios = pack(real(rec%nf, dp)/max(1, rec%q_nf), rec%outcome == 'solved' .and. rec%q_nf > 0)
    n_ratio = size(ratios)
    n_fewer = count(rec%outcome == 'solved' .and. rec%q_nf > 0 .and. rec%nf <= rec%q_nf)
    iters = pack(rec%iterations, rec%outcome == 'solved')
    write(u,'(A)') '| | sqpopt | NLPQLP |'
    write(u,'(A)') '|:--|--:|--:|'
    write(u,'(A)') '| objective evaluations (total) | '//s_i(sum(rec%nf, rec%outcome == 'solved'))//' | '// &
                   s_i(sum(rec%q_nf, rec%outcome == 'solved'))//' |'
    write(u,'(A)') '| gradient evaluations (total) | '//s_i(sum(rec%ng, rec%outcome == 'solved'))//' | '// &
                   s_i(sum(rec%q_ndf, rec%outcome == 'solved'))//' |'
    write(u,'(A)') '| iterations (total / median) | '//s_i(sum(iters))//' / '// &
                   s_r(median(real(iters, dp)),'(F8.1)')//' | |'
    write(u,'(A)') ''
    write(u,'(A)') '* median ratio of objective evaluations (sqpopt / NLPQLP): '//s_r(median(ratios),'(F8.2)')
    write(u,'(A)') '* problems where sqpopt used no more objective evaluations than NLPQLP: '// &
                   s_i(n_fewer)//' of '//s_i(n_ratio)

    ! ---- termination statuses ----
    write(u,'(A)') ''
    write(u,'(A)') '## Termination status'
    write(u,'(A)') ''
    write(u,'(A)') '| istat | status | solved | local | FAILED | total |'
    write(u,'(A)') '|--:|:--|--:|--:|--:|--:|'
    do istat = minval(rec%istat), maxval(rec%istat)
        if (.not. any(rec%istat == istat)) cycle
        do i = 1, 3
            cnt(i) = count(rec%istat == istat .and. rec%outcome == outcomes(i))
        end do
        write(u,'(A)') '| '//s_i(istat)//' | '//sqpopt_status_message(istat)//' | '//s_i(cnt(1))//' | '// &
                       s_i(cnt(2))//' | '//s_i(cnt(3))//' | '//s_i(sum(cnt))//' |'
    end do

    ! ---- unsolved problems ----
    write(u,'(A)') ''
    write(u,'(A)') '## Unsolved problems'
    write(u,'(A)') ''
    write(u,'(A)') '| TP | n | m | me | result | istat | status | iter | f | f* | rel. err | viol | notes |'
    write(u,'(A)') '|--:|--:|--:|--:|:--|--:|:--|--:|--:|--:|--:|--:|:--|'
    do i = 1, hs_n_problems
        if (rec(i)%outcome == 'solved') cycle
        write(u,'(A)') '| '//s_i(rec(i)%id)//' | '//s_i(rec(i)%n)//' | '//s_i(rec(i)%m)//' | '//s_i(rec(i)%me)// &
            ' | '//trim(rec(i)%outcome)//' | '//s_i(rec(i)%istat)//' | '//sqpopt_status_message(rec(i)%istat)// &
            ' | '//s_i(rec(i)%iterations)//' | '//s_r(rec(i)%f,'(ES12.4)')//' | '//s_r(rec(i)%f_star,'(ES12.4)')// &
            ' | '//s_r(rec(i)%rel,'(ES9.2)')//' | '//s_r(rec(i)%viol,'(ES9.2)')//' | '//notes(rec(i))//' |'
    end do

    ! ---- all problems ----
    write(u,'(A)') ''
    write(u,'(A)') '## All problems'
    write(u,'(A)') ''
    write(u,'(A)') '`nf`, `ng`: calls of the objective-and-constraints function `fc` and of the '// &
                   'gradient-and-Jacobian function `gjac`; '// &
                   '`Q:nf`, `Q:ng`: NLPQLP''s objective and gradient evaluations; `rel. err`: `(f - f*)/max(1,|f*|)`; '// &
                   '`viol`: largest constraint or bound violation; `KKT`: final (scaled) KKT error.'
    write(u,'(A)') ''
    write(u,'(A)') '| TP | n | m | me | result | istat | iter | nf | ng | Q:nf | Q:ng | f | f* '// &
                   '| rel. err | viol | KKT | notes |'
    write(u,'(A)') '|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|'
    do i = 1, hs_n_problems
        write(u,'(A)') '| '//s_i(rec(i)%id)//' | '//s_i(rec(i)%n)//' | '//s_i(rec(i)%m)//' | '//s_i(rec(i)%me)// &
            ' | '//trim(rec(i)%outcome)//' | '//s_i(rec(i)%istat)//' | '//s_i(rec(i)%iterations)// &
            ' | '//s_i(rec(i)%nf)//' | '//s_i(rec(i)%ng)// &
            ' | '//s_i(rec(i)%q_nf)//' | '//s_i(rec(i)%q_ndf)// &
            ' | '//s_r(rec(i)%f,'(ES12.4)')//' | '//s_r(rec(i)%f_star,'(ES12.4)')//' | '//s_r(rec(i)%rel,'(ES9.2)')// &
            ' | '//s_r(rec(i)%viol,'(ES9.2)')//' | '//s_r(rec(i)%kkt,'(ES9.2)')//' | '//notes(rec(i))//' |'
    end do

    close(u)

    end subroutine write_report

    subroutine write_web_data(file)
    !! write the results of all the problems (`rec`) as a JavaScript data file,
    !! `window.SQPOPT_HS_RESULTS = {...}` (a script rather than JSON, so that
    !! the page also works when opened as a local file): the configuration,
    !! and one object per problem with the columns of the report (non-finite
    !! values as `null`)
    character(len=*), intent(in) :: file
    integer :: u, i
    character(len=8)  :: date
    character(len=10) :: time
    call date_and_time(date, time)
    open(newunit=u, file=file, status='replace', action='write')
    write(u,'(A)') '// Generated by test/test_hs_suite.f90 (--web-data): the Hock-Schittkowski test suite results'
    write(u,'(A)') '// for the interactive results page (web/hs_results.html). Do not edit.'
    write(u,'(A)') 'window.SQPOPT_HS_RESULTS = {'
    write(u,'(A)') '  "generated": "'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '//time(1:2)//':'//time(3:4)//'",'
    write(u,'(A)') '  "compiler": "'//json_escape(compiler_version())//'",'
    write(u,'(A)') '  "options": "'//json_escape(merge('(defaults)', cfg_string, len(cfg_string) == 0))//'",'
    write(u,'(A)') '  "feas_tol": '//s_js(feas_tol)//', "rel_tol": '//s_js(rel_tol)//','
    write(u,'(A)') '  "problems": ['
    do i = 1, hs_n_problems
        associate (r => rec(i))
        write(u,'(A)') '    {"id": '//s_i(r%id)//', "n": '//s_i(r%n)//', "m": '//s_i(r%m)//', "me": '//s_i(r%me)// &
            ', "outcome": "'//trim(merge('failed', r%outcome, r%outcome == 'FAILED'))//'", "istat": '//s_i(r%istat)// &
            ', "status": "'//json_escape(sqpopt_status_message(r%istat))//'"'// &
            ', "iter": '//s_i(r%iterations)//', "nf": '//s_i(r%nf)//', "ng": '//s_i(r%ng)// &
            ', "q_nf": '//s_i(r%q_nf)//', "q_ng": '//s_i(r%q_ndf)// &
            ', "f": '//s_js(r%f)//', "f_star": '//s_js(r%f_star)//', "rel": '//s_js(r%rel)// &
            ', "viol": '//s_js(r%viol)//', "kkt": '//s_js(r%kkt)// &
            ', "fd": '//merge('true ', 'false', r%fd)//'}'//merge(',', ' ', i < hs_n_problems)
        end associate
    end do
    write(u,'(A)') '  ]'
    write(u,'(A)') '};'
    close(u)
    end subroutine write_web_data

    function s_js(x) result(s)
    !! a real as a JavaScript number literal (`null` if not finite)
    real(dp), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=32) :: buf
    if (.not. ieee_is_finite(x)) then
        s = 'null'
    else
        write(buf,'(ES24.16E3)') x
        s = trim(adjustl(buf))
    end if
    end function s_js

    function json_escape(t) result(s)
    !! `t` with backslashes and double quotes escaped, for a JavaScript string literal
    character(len=*), intent(in) :: t
    character(len=:), allocatable :: s
    integer :: i
    s = ''
    do i = 1, len_trim(t)
        if (t(i:i) == '\' .or. t(i:i) == '"') s = s//'\'
        s = s//t(i:i)
    end do
    end function json_escape

    function notes(r) result(s)
    !! the "notes" column of the report for record `r`
    type(run_record), intent(in) :: r
    character(len=:), allocatable :: s
    s = ''
    if (r%fd) s = '*fd*'
    if (r%regression) s = trim(s//' **regression**')
    s = adjustl(s)
    end function notes

    function s_i(i) result(s)
    !! integer to string
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=20) :: buf
    write(buf,'(I0)') i
    s = trim(buf)
    end function s_i

    function s_r(x, fmt) result(s)
    !! real to string, with the edit descriptor `fmt`
    real(dp),         intent(in) :: x
    character(len=*), intent(in) :: fmt
    character(len=:), allocatable :: s
    character(len=40) :: buf
    write(buf, fmt) x
    s = trim(adjustl(buf))
    end function s_r

    function median(a) result(med)
    !! median of `a` (0 if empty)
    real(dp), dimension(:), intent(in) :: a
    real(dp) :: med
    real(dp), dimension(size(a)) :: b
    real(dp) :: t
    integer :: i, j, n
    n = size(a)
    med = 0.0_dp
    if (n == 0) return
    b = a
    do i = 2, n   ! (insertion sort)
        t = b(i)
        j = i - 1
        do while (j >= 1)
            if (b(j) <= t) exit
            b(j+1) = b(j)
            j = j - 1
        end do
        b(j+1) = t
    end do
    med = merge(b((n+1)/2), 0.5_dp*(b(n/2) + b(n/2+1)), mod(n,2) == 1)
    end function median

    subroutine check_derivatives(ctx)
    !! compare the problem's gradient and Jacobian with finite differences at
    !! the starting point and at a nearby point, and switch to finite
    !! differences for whichever doesn't match. Each is compared at three
    !! difference steps (`h`, `h/100`, `10h`), and only counts as a mismatch
    !! if even the closest of the three disagrees (by more than 1e-3,
    !! relatively): finite differences of badly scaled functions can be
    !! inaccurate at any one step (e.g. TP376, 383), but a wrong derivative
    !! disagrees at all of them.
    type(problem_context), intent(inout) :: ctx
    real(dp), dimension(ctx%n) :: x, g
    real(dp), dimension(ctx%n,3) :: gd
    real(dp), dimension(ctx%m,ctx%n) :: jac
    real(dp), dimension(ctx%m,ctx%n,3) :: jd
    real(dp), dimension(3), parameter :: hf = [1.0_dp, 1.0e-2_dp, 1.0e1_dp]
    integer :: t, j_, k
    do t = 1, 2
        x = p%x0
        if (t == 2) x = x + 0.01_dp*max(1.0_dp, abs(x))*[(sin(real(3*t+j_, dp)), j_=1,ctx%n)]
        x = min(max(x, p%x_lb), p%x_ub)
        call hs_g(ctx%id, x, g)
        do k = 1, 3
            call fd_gradient(ctx%id, x, gd(:,k), hf(k))
        end do
        if (mismatch(g, gd)) ctx%fd_g = .true.
        if (ctx%m > 0) then
            call hs_jac(ctx%id, x, jac)
            do k = 1, 3
                call fd_jacobian(ctx%id, x, jd(:,:,k), hf(k))
            end do
            if (mismatch(reshape(jac, [size(jac)]), reshape(jd, [size(jac), 3]))) ctx%fd_jac = .true.
        end if
    end do
    end subroutine check_derivatives

    logical function mismatch(a, b)
    !! whether analytic values `a` differ from the finite-difference values
    !! `b(:,k)` (at three steps `k`) by more than 1e-3 (relatively) at every
    !! step. A non-finite analytic value where a finite difference is finite
    !! is a mismatch (e.g. the original TP25's gradient took a fractional
    !! power of a negative number); non-finite finite differences are skipped.
    real(dp), dimension(:),   intent(in) :: a
    real(dp), dimension(:,:), intent(in) :: b
    integer :: q, k
    real(dp) :: best
    mismatch = .false.
    do q = 1, size(a)
        best = huge(1.0_dp)
        do k = 1, size(b,2)
            if (.not. ieee_is_finite(b(q,k))) cycle
            if (.not. ieee_is_finite(a(q))) then
                best = huge(1.0_dp)
                exit
            end if
            best = min(best, abs(a(q) - b(q,k))/max(1.0_dp, abs(a(q)), abs(b(q,k))))
        end do
        if (all(.not. ieee_is_finite(b(q,:)))) cycle
        if (best > 1.0e-3_dp) mismatch = .true.
    end do
    end function mismatch

    subroutine fd_gradient(id, x, g, hfac)
    !! finite-difference gradient of problem `id`'s objective (see [[fd_step]])
    integer,                intent(in)  :: id
    real(dp), dimension(:), intent(in)  :: x
    real(dp), dimension(:), intent(out) :: g
    real(dp), optional,     intent(in)  :: hfac !! factor on the default difference step
    real(dp), dimension(size(x)) :: xp
    real(dp) :: f0, f1, f2, h
    integer :: j, side
    call hs_f(id, x, f0)
    do j = 1, size(x)
        call fd_step(x, j, h, side, hfac)
        xp = x; xp(j) = x(j) + h;          call hs_f(id, xp, f1)
        xp(j) = x(j) + merge(-h, 2*h, side == 0); call hs_f(id, xp, f2)
        if (side == 0) then
            g(j) = (f1 - f2)/(2.0_dp*h)
        else
            g(j) = (-3.0_dp*f0 + 4.0_dp*f1 - f2)/(2.0_dp*h)
        end if
    end do
    end subroutine fd_gradient

    subroutine fd_jacobian(id, x, jac, hfac)
    !! finite-difference Jacobian of problem `id`'s constraints (see [[fd_step]])
    integer,                  intent(in)  :: id
    real(dp), dimension(:),   intent(in)  :: x
    real(dp), dimension(:,:), intent(out) :: jac
    real(dp), optional,       intent(in)  :: hfac !! factor on the default difference step
    real(dp), dimension(size(x)) :: xp
    real(dp), dimension(size(jac,1)) :: c0, c1, c2
    real(dp) :: h
    integer :: j, side
    call hs_c(id, x, c0)
    do j = 1, size(x)
        call fd_step(x, j, h, side, hfac)
        xp = x; xp(j) = x(j) + h;          call hs_c(id, xp, c1)
        xp(j) = x(j) + merge(-h, 2*h, side == 0); call hs_c(id, xp, c2)
        if (side == 0) then
            jac(:,j) = (c1 - c2)/(2.0_dp*h)
        else
            jac(:,j) = (-3.0_dp*c0 + 4.0_dp*c1 - c2)/(2.0_dp*h)
        end if
    end do
    end subroutine fd_jacobian

    subroutine fd_step(x, j, h, side, hfac)
    !! the step for differencing along `x(j)`: central (`side=0`, points
    !! `x(j)+h` and `x(j)-h`), or, if that would cross one of the problem's
    !! variable bounds, a second-order one-sided difference into the
    !! interior (`side=1`, points `x(j)+h` and `x(j)+2h`, with `h<0` at an
    !! upper bound). Some problems' functions have a kink at a bound (e.g.
    !! TP358 clips `x` to its bounds), where a central difference would
    !! average across it (TP331, 358, 376, 383 start on such a bound).
    real(dp), dimension(:), intent(in)  :: x
    integer,                intent(in)  :: j
    real(dp),               intent(out) :: h
    integer,                intent(out) :: side
    real(dp), optional,     intent(in)  :: hfac
    h = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp, abs(x(j)))
    if (present(hfac)) h = hfac*h
    side = 0
    if (x(j) - h < p%x_lb(j)) then
        side = 1
    else if (x(j) + h > p%x_ub(j)) then
        side = 1
        h = -h
    end if
    end subroutine fd_step

    ! ---- the problem functions for `sqpopt` (the problem is identified by the user data) ----

    subroutine obj(x, f, status, data)
    real(wp), dimension(:), intent(in)    :: x
    real(wp),               intent(out)   :: f
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    real(dp) :: fd
    select type (data)
    type is (problem_context)
        call hs_f(data%id, real(x, dp), fd)
        f = real(fd, wp)
    end select
    associate(unused => status); end associate
    end subroutine obj

    subroutine grad(x, g, status, data)
    real(wp), dimension(:), intent(in)    :: x
    real(wp), dimension(:), intent(out)   :: g
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    real(dp), dimension(size(x)) :: gd
    select type (data)
    type is (problem_context)
        if (data%fd_g) then
            call fd_gradient(data%id, real(x, dp), gd)
        else
            call hs_g(data%id, real(x, dp), gd)
        end if
        g = real(gd, wp)
    end select
    associate(unused => status); end associate
    end subroutine grad

    subroutine cons(x, c, status, data)
    real(wp), dimension(:), intent(in)    :: x
    real(wp), dimension(:), intent(out)   :: c
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    real(dp), dimension(size(c)) :: cd
    select type (data)
    type is (problem_context)
        call hs_c(data%id, real(x, dp), cd)
        c = real(cd, wp)
    end select
    associate(unused => status); end associate
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    real(wp), dimension(:), intent(in)    :: x
    real(wp), dimension(:), intent(out)   :: jac_val
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    real(dp), dimension(:,:), allocatable :: jd
    select type (data)
    type is (problem_context)
        allocate(jd(data%m, data%n))
        if (data%fd_jac) then
            call fd_jacobian(data%id, real(x, dp), jd)
        else
            call hs_jac(data%id, real(x, dp), jd)
        end if
        jac_val = real(reshape(transpose(jd), [size(jd)]), wp)   ! (row by row, as the pattern)
    end select
    associate(unused => status); end associate
    end subroutine jacv

    subroutine hess_fd(x, lambda, hess_val, status, data)
    !! `hess` for `set_functions`: the Hessian of the Lagrangian
    !! \( \nabla^2 f - \sum_i \lambda_i \nabla^2 c_i \), by central differences
    !! (one-sided at a bound, see [[fd_step]]) of its analytic gradient
    !! \( \nabla f - J^T \lambda \), symmetrized; the lower triangle, row by row
    real(wp), dimension(:), intent(in)    :: x
    real(wp), dimension(:), intent(in)    :: lambda
    real(wp), dimension(:), intent(out)   :: hess_val
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    real(dp), dimension(size(x),size(x)) :: h
    real(dp), dimension(size(x)) :: xd, xp, gl0, gl1, gl2
    real(dp) :: step
    integer :: i, j, k, side
    select type (data)
    type is (problem_context)
        xd = real(x, dp)
        gl0 = lagrangian_gradient(data, xd, lambda)
        do j = 1, size(x)
            call fd_step(xd, j, step, side)
            xp = xd; xp(j) = xd(j) + step;                        gl1 = lagrangian_gradient(data, xp, lambda)
            xp(j) = xd(j) + merge(-step, 2*step, side == 0); gl2 = lagrangian_gradient(data, xp, lambda)
            if (side == 0) then
                h(:,j) = (gl1 - gl2)/(2.0_dp*step)
            else
                h(:,j) = (-3.0_dp*gl0 + 4.0_dp*gl1 - gl2)/(2.0_dp*step)
            end if
        end do
        h = 0.5_dp*(h + transpose(h))
        k = 0
        do i = 1, size(x)
            do j = 1, i
                k = k + 1
                hess_val(k) = real(h(i,j), wp)
            end do
        end do
    end select
    associate(unused => status); end associate
    end subroutine hess_fd

    function lagrangian_gradient(ctx, xx, lambda) result(gl)
    !! the analytic gradient of the Lagrangian, \( \nabla f - J^T \lambda \), for [[hess_fd]]
    type(problem_context),  intent(in) :: ctx
    real(dp), dimension(:), intent(in) :: xx
    real(wp), dimension(:), intent(in) :: lambda
    real(dp), dimension(size(xx)) :: gl
    real(dp), dimension(:,:), allocatable :: jd
    call hs_g(ctx%id, xx, gl)
    if (ctx%m > 0) then
        allocate(jd(ctx%m, ctx%n))
        call hs_jac(ctx%id, xx, jd)
        gl = gl - matmul(real(lambda, dp), jd)
    end if
    end function lagrangian_gradient

    subroutine fc_obj_cons(x, f, c, status, data)
    !! `fc` for `set_functions`: the objective (`obj`) and the constraints (`cons`)
    real(wp), dimension(:), intent(in)    :: x
    real(wp),               intent(out)   :: f
    real(wp), dimension(:), intent(out)   :: c
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    call obj(x, f, status, data)
    if (status == 0) call cons(x, c, status, data)
    end subroutine fc_obj_cons

    subroutine gjac_grad_jacv(x, g, jac_val, status, data)
    !! `gjac` for `set_functions`: the gradient (`grad`) and the Jacobian values (`jacv`)
    real(wp), dimension(:), intent(in)    :: x
    real(wp), dimension(:), intent(out)   :: g
    real(wp), dimension(:), intent(out)   :: jac_val
    integer,                intent(inout) :: status
    class(*), optional,     intent(inout) :: data
    call grad(x, g, status, data)
    if (status == 0) call jacv(x, jac_val, status, data)
    end subroutine gjac_grad_jacv


end program test_hs_suite
