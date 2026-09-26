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
    !! to `test/hs_suite_results.md`, or to the file given as the first
    !! command-line argument:
    !!
    !!    fpm test test_hs_suite --profile release -- results.md

    use, intrinsic :: ieee_arithmetic, only: ieee_set_halting_mode, ieee_all, ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, compiler_version
    use hs_problems_module
    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_linesearch_module, only: sqpopt_linesearch_type
    use sqpopt_types_module,   only: sqpopt_results_type, sqpopt_success, sqpopt_stalled, sqpopt_acceptable, &
                                     sqpopt_status_message
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
        integer  :: nf = 0, ng = 0, nc = 0, njac = 0  !! evaluation counts
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
    integer :: arg_len

    call ieee_set_halting_mode(ieee_all, .false.)  ! (trial points may produce NaN/Inf, which the solver handles)

    write(*,*) '----------------------------'
    write(*,*) 'test_hs_suite'
    write(*,*) '----------------------------'
    write(*,'(A5,3A4,A4,A7,2A17,2A10,2A6,2A7,A4)') 'TP', 'n', 'm', 'me', 'st', 'result', 'f', 'f*', &
        'rel.err', 'viol', 'nf', 'ng', 'Q:nf', 'Q:ndf', 'fd'

    n_solved = 0; n_local = 0; n_failed = 0; n_fd = 0; n_regressions = 0; n_improved = 0
    sum_nf = 0; sum_ng = 0; sum_nlpqlp_nf = 0; sum_nlpqlp_ndf = 0
    call system_clock(t0, rate)

    do k = 1, hs_n_problems
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

    call get_command_argument(1, length=arg_len)
    if (arg_len > 0) then
        allocate(character(len=arg_len) :: report_file)
        call get_command_argument(1, report_file)
    else
        report_file = default_report_file
    end if
    call write_report(report_file, real(t1-t0, dp)/real(rate, dp))
    write(*,'(A)') 'report written to: '//report_file

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

    subroutine run_problem(id, k)
    !! solve problem `id` (the `k`-th in the collection) and report the outcome
    integer, intent(in) :: id, k

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_linesearch_type) :: linesearch
    type(sqpopt_results_type) :: r
    integer, dimension(:), allocatable :: irow, icol
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
    call problem%set_functions(f=obj, g=grad, c=cons, jac=jacv, data=ctx)
    options%max_iter = 1000

    call solver%initialize(problem=problem, options=options, linesearch=linesearch)
    call solver%solve(real(p%x0, wp), istat)
    call solver%get_results(r)

    viol = max(r%feasibility_error, maxval(max(p%x_lb - r%x, 0.0_dp) + max(r%x - p%x_ub, 0.0_dp)))
    rel  = (r%f - p%f_star)/max(1.0_dp, abs(p%f_star))
    feasible  = viol <= feas_tol .and. ieee_is_finite(r%f)
    converged = istat == sqpopt_success .or. istat == sqpopt_stalled .or. istat == sqpopt_acceptable
    if (feasible .and. rel <= rel_tol) then
        outcome = 'solved'
        n_solved = n_solved + 1
        sum_nf = sum_nf + r%n_eval_f
        sum_ng = sum_ng + r%n_eval_g
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
                        nf=r%n_eval_f, ng=r%n_eval_g, nc=r%n_eval_c, njac=r%n_eval_jac, &
                        q_nf=hs_nlpqlp_nf(k), q_ndf=hs_nlpqlp_ndf(k), &
                        f=r%f, f_star=p%f_star, rel=rel, viol=viol, kkt=r%kkt_error, outcome=outcome, &
                        fd=ctx%fd_g .or. ctx%fd_jac, regression=outcome /= 'solved' .and. .not. expected_unsolved)

    write(*,'(I5,3I4,I4,A7,2ES17.8,2ES10.2,2I6,2I7,A4,A)') id, p%n, p%m, p%me, istat, outcome, r%f, p%f_star, &
        rel, viol, r%n_eval_f, r%n_eval_g, hs_nlpqlp_nf(k), hs_nlpqlp_ndf(k), &
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
    write(u,'(A)') '`nf`, `ng`, `nc`, `nJ`: calls of the objective, gradient, constraint, and Jacobian functions; '// &
                   '`Q:nf`, `Q:ng`: NLPQLP''s objective and gradient evaluations; `rel. err`: `(f - f*)/max(1,|f*|)`; '// &
                   '`viol`: largest constraint or bound violation; `KKT`: final (scaled) KKT error.'
    write(u,'(A)') ''
    write(u,'(A)') '| TP | n | m | me | result | istat | iter | nf | ng | nc | nJ | Q:nf | Q:ng | f | f* '// &
                   '| rel. err | viol | KKT | notes |'
    write(u,'(A)') '|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|'
    do i = 1, hs_n_problems
        write(u,'(A)') '| '//s_i(rec(i)%id)//' | '//s_i(rec(i)%n)//' | '//s_i(rec(i)%m)//' | '//s_i(rec(i)%me)// &
            ' | '//trim(rec(i)%outcome)//' | '//s_i(rec(i)%istat)//' | '//s_i(rec(i)%iterations)// &
            ' | '//s_i(rec(i)%nf)//' | '//s_i(rec(i)%ng)//' | '//s_i(rec(i)%nc)//' | '//s_i(rec(i)%njac)// &
            ' | '//s_i(rec(i)%q_nf)//' | '//s_i(rec(i)%q_ndf)// &
            ' | '//s_r(rec(i)%f,'(ES12.4)')//' | '//s_r(rec(i)%f_star,'(ES12.4)')//' | '//s_r(rec(i)%rel,'(ES9.2)')// &
            ' | '//s_r(rec(i)%viol,'(ES9.2)')//' | '//s_r(rec(i)%kkt,'(ES9.2)')//' | '//notes(rec(i))//' |'
    end do

    close(u)

    end subroutine write_report

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

end program test_hs_suite
