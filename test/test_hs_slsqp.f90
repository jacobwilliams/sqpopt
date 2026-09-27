program test_hs_slsqp

    !! The Hock-Schittkowski test collection (305 problems, see
    !! [[hs_problems_module]]) solved with [SLSQP](https://github.com/jacobwilliams/slsqp)
    !! (a dev-dependency), for comparison with SQPOPT (`test_hs_suite.f90`)
    !! on the interactive results page of the user guide (`web/hs_results.html`).
    !!
    !! Each problem is solved from its standard starting point (projected onto
    !! the variable bounds, as SQPOPT does: a few problems start outside them,
    !! e.g. TP45, where SLSQP would otherwise stop at once), with the same
    !! derivatives as `test_hs_suite` (central differences where the analytic
    !! ones don't match, see [[hs_derivatives_module]]), at most 1000
    !! iterations, `acc = 1e-8`, and SLSQP's inexact line search, and
    !! classified the same way:
    !!
    !! * **solved**: feasible (violation <= `feas_tol`), with an objective
    !!   value within `rel_tol` (relative) of the validated optimum, or better;
    !! * **local**: feasible and converged (`istat = 0`), but at a worse
    !!   objective value -- a different local solution;
    !! * **failed**: anything else.
    !!
    !! SLSQP's constraints are equalities `c = 0` (first) and inequalities
    !! `c >= 0`, so each problem's constraints `c_l <= c(x) <= c_u` are mapped
    !! to them: `c - c_l = 0` for an equality, and `c - c_l >= 0` and/or
    !! `c_u - c >= 0` for the finite sides of an inequality. The evaluation
    !! counts are the calls of the function (`f` and `c` together, like
    !! SQPOPT's `fc`) and gradient (`g` and the Jacobian, like `gjac`)
    !! routines.
    !!
    !! There is no regression test (SLSQP is not the code under test). With
    !! `--web-data=FILE`, the results are written as a JavaScript data file
    !! (`window.SQPOPT_HS_SLSQP = {...}`, the same format as
    !! `test_hs_suite`'s), for the results page:
    !!
    !!    fpm test test_hs_slsqp --profile release -- --web-data=web/js/hs_slsqp_data.js

    use slsqp_module,          only: slsqp_solver
    use hs_problems_module
    use hs_derivatives_module, only: p => hs_current, check_derivatives, fd_gradient, fd_jacobian
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, compiler_version

    implicit none

    real(dp), parameter :: rel_tol  = 1.0e-4_dp  !! relative objective tolerance for "solved" (as `test_hs_suite`)
    real(dp), parameter :: feas_tol = 1.0e-6_dp  !! constraint/bound violation tolerance for "solved"/"local"
    real(dp), parameter :: acc      = 1.0e-8_dp  !! SLSQP's accuracy
    integer,  parameter :: max_iter = 1000       !! SLSQP's iteration limit (as `test_hs_suite`)

    type :: run_record
        !! the outcome of one problem
        integer  :: id = 0, n = 0, m = 0, me = 0, istat = 0, iterations = 0, nf = 0, ng = 0
        real(dp) :: f = 0.0_dp, f_star = 0.0_dp, rel = 0.0_dp, viol = 0.0_dp
        character(len=6) :: outcome = ''
        character(len=:), allocatable :: status
        logical  :: fd = .false.
    end type run_record

    type(run_record), dimension(hs_n_problems) :: rec

    ! the current problem (SLSQP's callbacks take no user data):
    integer :: cur_id = 0
    logical :: fd_g = .false., fd_jac = .false.
    integer :: ms = 0, meqs = 0                   !! SLSQP's number of constraints, and of equalities
    integer,  dimension(:), allocatable :: row    !! the problem's constraint behind each SLSQP constraint
    real(dp), dimension(:), allocatable :: sgn    !! its sign (`+1`: `c - b`, `-1`: `b - c`)
    real(dp), dimension(:), allocatable :: rhs    !! its bound `b`
    integer :: n_fc = 0, n_gjac = 0

    character(len=:), allocatable :: web_data_file
    integer :: k, n_solved, n_local, n_failed, sum_nf, sum_ng

    call parse_arguments()

    write(*,*) '----------------------------'
    write(*,*) 'test_hs_slsqp'
    write(*,*) '----------------------------'
    write(*,'(A)') '   TP   n   m  me   istat result                f               f*   rel.err      viol    nf    ng  fd'

    n_solved = 0; n_local = 0; n_failed = 0; sum_nf = 0; sum_ng = 0
    do k = 1, hs_n_problems
        call run_problem(hs_problem_ids(k), rec(k))
        select case (rec(k)%outcome)
        case ('solved')
            n_solved = n_solved + 1
            sum_nf = sum_nf + rec(k)%nf
            sum_ng = sum_ng + rec(k)%ng
        case ('local')
            n_local = n_local + 1
        case default
            n_failed = n_failed + 1
        end select
    end do

    write(*,'(A)') ''
    write(*,'(5(A,I0))') 'summary: solved=', n_solved, ' local=', n_local, ' failed=', n_failed, &
                         ' nf=', sum_nf, ' ng=', sum_ng
    if (len(web_data_file) > 0) then
        call write_web_data(web_data_file)
        write(*,'(A)') 'web data written to: '//web_data_file
    end if
    print '(A)', 'test_hs_slsqp PASSED'

    contains

    subroutine parse_arguments()
    !! `--web-data=FILE` (the only option)
    integer :: i
    character(len=256) :: arg
    web_data_file = ''
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        if (arg(1:11) == '--web-data=') then
            web_data_file = trim(arg(12:))
        else
            write(*,'(A)') 'test_hs_slsqp: unknown option: '//trim(arg)
            error stop 1
        end if
    end do
    end subroutine parse_arguments

    subroutine run_problem(id, r)
    !! solve problem `id` with SLSQP and classify the outcome
    integer,          intent(in)  :: id
    type(run_record), intent(out) :: r

    type(slsqp_solver) :: solver
    real(dp), dimension(:), allocatable :: x, c
    real(dp) :: f, viol
    integer  :: i, istat, iterations
    logical  :: status_ok
    character(len=:), allocatable :: message

    call hs_setup(id, p)
    cur_id = id
    fd_g = .false.; fd_jac = .false.
    call check_derivatives(id, p%n, p%m, fd_g, fd_jac)
    call map_constraints()

    n_fc = 0; n_gjac = 0
    x = min(max(p%x0, p%x_lb), p%x_ub)
    call solver%initialize(p%n, ms, meqs, max_iter, acc, slsqp_func, slsqp_grad, p%x_lb, p%x_ub, &
                           status_ok, linesearch_mode=1, infinite_bound=hs_infinity, iprint=0)
    if (.not. status_ok) error stop 'test_hs_slsqp: SLSQP initialization failed'
    call solver%optimize(x, istat, iterations, message)

    ! the outcome, for the original problem at the final point:
    call hs_f(id, x, f)
    allocate(c(p%m))
    viol = maxval([0.0_dp, max(p%x_lb - x, 0.0_dp) + max(x - p%x_ub, 0.0_dp)])
    if (p%m > 0) then
        call hs_c(id, x, c)
        viol = max(viol, maxval(max(p%c_lb - c, 0.0_dp) + max(c - p%c_ub, 0.0_dp)))
    end if
    if (.not. (ieee_is_finite(viol))) viol = huge(1.0_dp)

    r%id = id; r%n = p%n; r%m = p%m; r%me = p%me
    r%istat = istat; r%iterations = iterations; r%nf = n_fc; r%ng = n_gjac
    r%f = f; r%f_star = p%f_star; r%rel = (f - p%f_star)/max(1.0_dp, abs(p%f_star)); r%viol = viol
    r%status = message
    r%fd = fd_g .or. fd_jac
    if (viol <= feas_tol .and. ieee_is_finite(f) .and. r%rel <= rel_tol) then
        r%outcome = 'solved'
    else if (viol <= feas_tol .and. ieee_is_finite(f) .and. istat == 0) then
        r%outcome = 'local'
    else
        r%outcome = 'failed'
    end if

    write(*,'(I5,3I4,I8,1X,A7,2ES17.8,2ES10.2,2I6,A4)') id, p%n, p%m, p%me, istat, r%outcome, f, p%f_star, &
        r%rel, min(viol, 9.99e99_dp), n_fc, n_gjac, merge(' fd', '   ', r%fd)
    associate (unused => i); end associate

    end subroutine run_problem

    subroutine map_constraints()
    !! SLSQP's constraints for the current problem (see the program documentation)
    integer :: i, j
    ms = 0
    meqs = count(p%c_lb == p%c_ub)
    do i = 1, p%m
        if (p%c_lb(i) == p%c_ub(i)) then
            ms = ms + 1
        else
            if (p%c_lb(i) > -hs_infinity) ms = ms + 1
            if (p%c_ub(i) <  hs_infinity) ms = ms + 1
        end if
    end do
    if (allocated(row)) deallocate(row, sgn, rhs)
    allocate(row(ms), sgn(ms), rhs(ms))
    j = 0
    do i = 1, p%m                      ! the equalities first
        if (p%c_lb(i) == p%c_ub(i)) then
            j = j + 1; row(j) = i; sgn(j) = 1.0_dp; rhs(j) = p%c_lb(i)
        end if
    end do
    do i = 1, p%m                      ! then the finite sides of the inequalities
        if (p%c_lb(i) == p%c_ub(i)) cycle
        if (p%c_lb(i) > -hs_infinity) then
            j = j + 1; row(j) = i; sgn(j) = 1.0_dp; rhs(j) = p%c_lb(i)
        end if
        if (p%c_ub(i) < hs_infinity) then
            j = j + 1; row(j) = i; sgn(j) = -1.0_dp; rhs(j) = p%c_ub(i)
        end if
    end do
    end subroutine map_constraints

    subroutine slsqp_func(me, x, f, c)
    !! SLSQP's problem function: the objective and the mapped constraints
    class(slsqp_solver),    intent(inout) :: me
    real(dp), dimension(:), intent(in)    :: x
    real(dp),               intent(out)   :: f
    real(dp), dimension(:), intent(out)   :: c
    real(dp), dimension(p%m) :: cp
    n_fc = n_fc + 1
    call hs_f(cur_id, x, f)
    if (ms > 0) then
        call hs_c(cur_id, x, cp)
        c(1:ms) = sgn*(cp(row) - rhs)
    end if
    associate (unused => me); end associate
    end subroutine slsqp_func

    subroutine slsqp_grad(me, x, g, a)
    !! SLSQP's gradient function: the objective gradient and the mapped
    !! constraints' Jacobian (central differences where `test_hs_suite` uses them)
    class(slsqp_solver),      intent(inout) :: me
    real(dp), dimension(:),   intent(in)    :: x
    real(dp), dimension(:),   intent(out)   :: g
    real(dp), dimension(:,:), intent(out)   :: a
    real(dp), dimension(p%m, p%n) :: jac
    integer :: j
    n_gjac = n_gjac + 1
    if (fd_g) then
        call fd_gradient(cur_id, x, g)
    else
        call hs_g(cur_id, x, g)
    end if
    a = 0.0_dp
    if (ms > 0) then
        if (fd_jac) then
            call fd_jacobian(cur_id, x, jac)
        else
            call hs_jac(cur_id, x, jac)
        end if
        do j = 1, ms
            a(j, 1:p%n) = sgn(j)*jac(row(j), :)
        end do
    end if
    associate (unused => me); end associate
    end subroutine slsqp_grad

    subroutine write_web_data(file)
    !! the results as a JavaScript data file (the format of `test_hs_suite`'s
    !! `--web-data`, without NLPQLP's counts and the KKT error)
    character(len=*), intent(in) :: file
    integer :: u, i
    character(len=8)  :: date
    character(len=10) :: time
    call date_and_time(date, time)
    open(newunit=u, file=file, status='replace', action='write')
    write(u,'(A)') '// Generated by test/test_hs_slsqp.f90 (--web-data): the Hock-Schittkowski test suite results'
    write(u,'(A)') '// of SLSQP, for the interactive results page (web/hs_results.html). Do not edit.'
    write(u,'(A)') 'window.SQPOPT_HS_SLSQP = {'
    write(u,'(A)') '  "generated": "'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '//time(1:2)//':'//time(3:4)//'",'
    write(u,'(A)') '  "compiler": "'//json_escape(compiler_version())//'",'
    write(u,'(A)') '  "options": "acc = 1e-8, max_iter = '//s_i(max_iter)//', inexact line search, x0 projected onto the bounds",'
    write(u,'(A)') '  "feas_tol": '//s_js(feas_tol)//', "rel_tol": '//s_js(rel_tol)//','
    write(u,'(A)') '  "problems": ['
    do i = 1, hs_n_problems
        associate (r => rec(i))
        write(u,'(A)') '    {"id": '//s_i(r%id)//', "n": '//s_i(r%n)//', "m": '//s_i(r%m)//', "me": '//s_i(r%me)// &
            ', "outcome": "'//trim(r%outcome)//'", "istat": '//s_i(r%istat)// &
            ', "status": "'//json_escape(r%status)//'"'// &
            ', "iter": '//s_i(r%iterations)//', "nf": '//s_i(r%nf)//', "ng": '//s_i(r%ng)// &
            ', "f": '//s_js(r%f)//', "f_star": '//s_js(r%f_star)//', "rel": '//s_js(r%rel)// &
            ', "viol": '//s_js(r%viol)//', "fd": '//merge('true ', 'false', r%fd)//'}'//merge(',', ' ', i < hs_n_problems)
        end associate
    end do
    write(u,'(A)') '  ]'
    write(u,'(A)') '};'
    close(u)
    end subroutine write_web_data

    function s_i(i) result(s)
    !! integer to string
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=16) :: buf
    write(buf,'(I0)') i
    s = trim(buf)
    end function s_i

    function s_js(x) result(s)
    !! a real as a JavaScript number literal (`null` if not finite or huge)
    real(dp), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=32) :: buf
    if (.not. ieee_is_finite(x) .or. abs(x) >= huge(1.0_dp)) then
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

end program test_hs_slsqp
