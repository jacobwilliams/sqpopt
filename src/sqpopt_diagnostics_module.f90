!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Diagnostics of a solve, to help when it doesn't converge (see
!  `options%diagnostic_level`). The iteration log says what happened at
!  each iteration; the diagnostics say which constraints and variables are
!  responsible, and what the final status most likely means. They are
!  opt-in, by level:
!
!  * `1`: a diagnosis of the final point ([[diagnostics_finish]]): the
!    constraints with the largest violations and multipliers, the variables
!    with the largest stationarity residuals, the active set and its
!    degeneracy (multipliers that are zero or have the wrong sign, active
!    gradients that are linearly dependent), the last QP subproblem (its
!    elastic constraints), and a verdict for the final status. It uses the
!    values the solver has at the final point anyway.
!  * `2`: also a report on the starting point ([[diagnostics_start]]: the
!    sizes of the gradients, which show how the problem is scaled, and the
!    structure of the constraints), and what the iterations showed: the
!    rate of convergence, the step lengths, a working set that flips between
!    two states, the iteration that took the most time
!    ([[diagnostics_record]]), and constraints whose changes along the steps
!    disagree with their derivatives ([[diagnostics_point]]). With
!    `options%diagnostics_unit`, the history of the iterations is also
!    written to a file, as comma-separated values.
!  * `3`: also the diagnostics that evaluate the user's functions: `f` and
!    `c` at two points next to the starting point, to see whether it is at
!    the edge of their domain.
!
!  Levels `0` to `2` call the user's functions exactly as often as a solve
!  without diagnostics: everything is computed from values the solver has
!  evaluated (through the problem's evaluation caches), and no level
!  changes the iterates. The results are returned in
!  `results%diagnosis` (see [[sqpopt_diagnosis_type]]), whose `report` and
!  `problem_report` are the texts the solver prints (if
!  `options%print_level >= 1`).
!
!  **The derivative check of level 2** needs no evaluations. Between two
!  consecutive iterates \( x \) and \( x^+ = x + s \), the change of a
!  function is \( \tfrac12 (\nabla c(x) + \nabla c(x^+))^T s \) up to a term
!  of the third order in \( s \) (the trapezoid rule, exact for a quadratic),
!  and the solver has the values and the derivatives at both. A constraint
!  (or the objective) whose actual changes differ from this by more than
!  20% in most of the steps has derivatives that don't match its values.
!  Only steps that change no variable by more than 10% of its size are
!  used, so that the third-order term is small.

    module sqpopt_diagnostics_module

    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use sqpopt_kinds,            only: wp => sqpopt_module_wp
    use sqpopt_types_module,     only: sqpopt_sparse_matrix, sqpopt_results_type, sqpopt_diagnosis_type, &
                                       sqpopt_iter_info, sqpopt_iteration_flags, sqpopt_all_finite, sqpopt_infinity, &
                                       sqpopt_success, sqpopt_acceptable, sqpopt_stalled, sqpopt_max_iter_reached, &
                                       sqpopt_user_requested_stop, sqpopt_max_evals_reached, &
                                       sqpopt_time_limit_reached, sqpopt_infeasible, sqpopt_line_search_failed, &
                                       sqpopt_qp_solve_failed, sqpopt_function_error, sqpopt_unbounded, &
                                       sqpopt_out_of_memory
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type
    use sqpopt_qp_daqp_module,   only: daqp_status_text
    use sqpopt_linalg_module,    only: independent_columns, sparse_matvec_transpose
    use sqpopt_log_module,       only: fmt_e, fmt_i, plural, qp_status_text

    implicit none

    private

    integer,  parameter :: n_listed = 5   !! the most constraints or variables listed for any one finding
    integer,  parameter :: n_window = 20  !! the most iterations (the last ones) the convergence rate is fitted over
    real(wp), parameter :: badly_scaled = 1.0e6_wp !! a ratio of gradient sizes above this is reported as poor scaling
    character(len=1), parameter :: nl = new_line('a') !! ends each line of a report
    integer,  parameter :: line_width = 110 !! the reports' lines are broken to this width

    type, public :: sqpopt_diagnostics_type
        !! the state of the diagnostics of one solve (a local variable of
        !! the solve: nothing is kept from one solve to the next).

        private

        integer, public :: level = 0  !! the diagnostic level (`options%diagnostic_level`; set by [[diagnostics_start]])
        integer :: unit = -1          !! the unit of the history file (`-1`: none)
        type(sqpopt_diagnosis_type) :: found !! what the starting point showed (the rest is set by [[diagnostics_finish]])

        ! the previous iterate, for the derivative check (see the module documentation):
        logical  :: have_point = .false.                 !! whether there is one
        real(wp) :: f_prev = 0.0_wp                      !! the objective there
        real(wp), dimension(:), allocatable :: x_prev    !! the point `dimension(n)`
        real(wp), dimension(:), allocatable :: g_prev    !! the objective gradient there `dimension(n)`
        real(wp), dimension(:), allocatable :: c_prev    !! the constraints there `dimension(m)`
        real(wp), dimension(:), allocatable :: jac_prev  !! the Jacobian's values there `dimension(jac_nnz)`
        real(wp), dimension(:), allocatable :: predicted !! work: the predicted changes of the constraints `dimension(m)`
        integer :: n_checked_f = 0                       !! steps in which the objective's change was checked
        integer :: n_off_f     = 0                       !! of which, it disagreed with the gradients
        integer, dimension(:), allocatable :: n_checked  !! the same for each constraint `dimension(m)`
        integer, dimension(:), allocatable :: n_off

        ! the history of the iterations:
        integer  :: n_err = 0                        !! iterations whose error was recorded
        real(wp), dimension(n_window) :: err = 0.0_wp !! the error (relative to the tolerances) of the last ones
                                                      !! (iteration `k` is element `mod(k-1, n_window) + 1`)
        integer  :: n_alpha = 0                      !! steps taken
        real(wp), dimension(n_window) :: alpha = 1.0_wp !! the step lengths of the last ones (stored like `err`)
        integer, dimension(:), allocatable :: set1   !! the QP's working set after the last step `dimension(m+n)`
        integer, dimension(:), allocatable :: set2   !! and after the one before it
        integer  :: n_flips = 0                      !! consecutive iterations in which it went back to `set2`

        ! the time:
        integer(int64) :: clock = 0     !! the system clock at the end of the last iteration
        real(wp) :: t_user = 0.0_wp     !! the time in the user's functions up to then
        real(wp) :: t_qp   = 0.0_wp     !! the time in the QP solver up to then
        real(wp) :: t_fact = 0.0_wp     !! the time in the factorizations up to then
        real(wp) :: t_total = 0.0_wp    !! the time of the iterations so far
        integer  :: slowest = 0         !! the iteration that took the most time
        real(wp) :: t_slowest = 0.0_wp  !! its time
        real(wp), dimension(3) :: part_slowest = 0.0_wp !! of which: in the user's functions, the QP solver, and the
                                                        !! factorizations

        contains

        procedure, public :: start  => diagnostics_start
        procedure, public :: point  => diagnostics_point
        procedure, public :: record => diagnostics_record
        procedure, public :: finish => diagnostics_finish
        procedure, public :: problem_report => diagnostics_problem_report

    end type sqpopt_diagnostics_type

    public :: sqpopt_diagnostics_write

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  start the diagnostics of a solve: set the level, and from level `2`
!  report on the starting point `x` (the sizes of the gradients there,
!  which show how the problem is scaled, and the structure of the
!  constraints), set up the history, and write the header line of the
!  history file. At level `3`, also evaluate `f` and `c` at two points next
!  to `x`. Call it after the problem's scaling has been computed: the
!  gradient and the Jacobian it asks for are then already in the evaluation
!  cache (and without scaling, they stay there for the first iteration), so
!  it costs no evaluation.

    subroutine diagnostics_start(me, problem, options, x)

    class(sqpopt_diagnostics_type), intent(inout) :: me
    type(sqpopt_problem_type),      intent(inout) :: problem !! the problem (scaled)
    type(sqpopt_options_type),      intent(in)    :: options !! the solver options
    real(wp), dimension(:),         intent(in)    :: x       !! the starting point `dimension(n)`

    real(wp), dimension(:), allocatable :: g, jval, rowmax, colmax, vv
    integer,  dimension(:), allocatable :: rowcnt, colcnt, eqmap, ir, ic
    logical,  dimension(:), allocatable :: indep
    character(len=:), allocatable :: text, line
    integer,  dimension(n_listed) :: idx
    integer  :: n, m, nnz, i, j, k, cnt, i_min, i_max, j_min, j_max, n_fixed, m_eq, nz, stat, ios
    real(wp) :: v, ratio

    me%level = options%diagnostic_level
    me%found%level = me%level
    if (me%level < 2) return

    n   = problem%n
    m   = problem%m
    nnz = max(problem%jac_nnz, 0)
    me%unit = options%diagnostics_unit
    call system_clock(me%clock)
    me%t_user = problem%time_user

    allocate(me%x_prev(n), me%g_prev(n), me%c_prev(m), me%jac_prev(nnz), me%predicted(m), me%n_checked(m), me%n_off(m))
    me%n_checked = 0
    me%n_off     = 0

    if (me%unit /= -1) then
        write(me%unit, '(A)', iostat=ios) 'iteration,objective,infeasibility,kkt_error,stationarity,alpha,step_norm,'// &
            'fc,fc_iteration,qp_iterations,qp_status,largest_multiplier,globalization,hessian,working_set_changes,'// &
            'time,time_functions,time_qp,time_factorization,flags'
    end if

    text = ''
    call add(text, '   diagnostics of the starting point:')

    allocate(g(n), jval(nnz))
    call problem%g(x, g)
    call problem%jac(x, jval)
    if (.not. (sqpopt_all_finite(g) .and. sqpopt_all_finite(jval))) then
        call add(text, '     the gradient or the Jacobian is not finite at the starting point')
        me%found%problem_report = text
        return
    end if

    ! ---- the sizes of the gradients, of the original problem ----
    allocate(rowmax(m), colmax(n), rowcnt(m), colcnt(n))
    rowmax = 0.0_wp
    rowcnt = 0
    colcnt = 0
    do j = 1, n
        colmax(j) = abs(g(j))/problem%f_scale
    end do
    do k = 1, nnz
        i = problem%jac_irow(k)
        j = problem%jac_icol(k)
        v = abs(jval(k))/problem%c_scale(i)
        rowmax(i) = max(rowmax(i), v)
        colmax(j) = max(colmax(j), v)
        rowcnt(i) = rowcnt(i) + 1
        colcnt(j) = colcnt(j) + 1
    end do

    if (n > 0) call add(text, '     objective gradient: largest element '//fmt_e(maxval(abs(g))/problem%f_scale))
    call extremes(rowmax, i_min, i_max)
    if (i_min > 0) then
        ratio = rowmax(i_max)/rowmax(i_min)
        me%found%constraint_gradient_ratio = ratio
        if (i_min == i_max) then
            call add(text, '     constraint gradients: largest element '//fmt_e(rowmax(i_max)))
        else
            call add(text, '     constraint gradients: largest element from '//fmt_e(rowmax(i_min))//' (c('// &
                     fmt_i(i_min)//')) to '//fmt_e(rowmax(i_max))//' (c('//fmt_i(i_max)//')), a ratio of '//fmt_e(ratio))
        end if
        if (ratio > badly_scaled) call add(text, '       the constraints differ widely in scale: the automatic '// &
            'scaling only scales down those above scaling_max_gradient, so consider scaling the small ones up')
    end if
    call extremes(colmax, j_min, j_max)
    if (j_min > 0) then
        ratio = colmax(j_max)/colmax(j_min)
        me%found%variable_gradient_ratio = ratio
        if (j_min /= j_max) call add(text, '     variables: largest gradient or Jacobian element from '// &
            fmt_e(colmax(j_min))//' (x('//fmt_i(j_min)//')) to '//fmt_e(colmax(j_max))//' (x('//fmt_i(j_max)// &
            ')), a ratio of '//fmt_e(ratio))
        if (ratio > badly_scaled) call add(text, '       the variables differ widely in scale: only the objective '// &
            'and the constraints are scaled automatically, so consider rescaling the variables')
    end if
    if (options%scaling .and. options%scaling_min_value > 0.0_wp) then
        cnt = 0
        do i = 1, m
            if (problem%c_scale(i) <= options%scaling_min_value) cnt = cnt + 1
        end do
        if (problem%f_scale <= options%scaling_min_value .or. cnt > 0) then
            line = '     scaling: '
            if (problem%f_scale <= options%scaling_min_value) line = line//'the objective'
            if (problem%f_scale <= options%scaling_min_value .and. cnt > 0) line = line//' and '
            if (cnt > 0) line = line//plural(cnt, 'constraint', 'constraints')
            call add(text, line//' at the smallest scale factor (scaling_min_value = '// &
                     fmt_e(options%scaling_min_value)//'): the gradient at the starting point is extremely large')
        end if
    end if

    ! ---- the structure of the constraints ----
    n_fixed = 0
    do j = 1, n
        if (problem%x_ub(j) - problem%x_lb(j) <= 0.0_wp) n_fixed = n_fixed + 1
    end do
    if (n_fixed > 0) call add(text, '     structure: '//plural(n_fixed, 'variable is', 'variables are')//' fixed')

    cnt = 0
    do i = 1, m
        if (rowcnt(i) == 0) call note(idx, cnt, i)
    end do
    me%found%n_constant_constraints = cnt
    if (cnt > 0) call add(text, '     structure: '//plural(cnt, 'constraint has', 'constraints have')// &
        ' no entry in the Jacobian''s pattern (constant): '//listed('c', idx, cnt))

    cnt = 0
    do i = 1, m
        if (rowcnt(i) == 1) call note(idx, cnt, i)
    end do
    me%found%n_single_variable_constraints = cnt
    if (cnt > 0) call add(text, '     structure: '//plural(cnt, 'constraint depends', 'constraints depend')// &
        ' on a single variable (if linear, a bound on the variable is cheaper): '//listed('c', idx, cnt))

    cnt = 0
    do i = 1, m
        if (rowcnt(i) > 0 .and. rowmax(i) == 0.0_wp) call note(idx, cnt, i)
    end do
    if (cnt > 0) call add(text, '     structure: '//plural(cnt, 'constraint has', 'constraints have')// &
        ' a zero gradient at the starting point: '//listed('c', idx, cnt))

    cnt = 0
    do j = 1, n
        if (colmax(j) == 0.0_wp) call note(idx, cnt, j)
    end do
    if (cnt > 0) call add(text, '     structure: '//plural(cnt, 'variable has', 'variables have')// &
        ' no effect at the starting point (zero gradient and Jacobian column): '//listed('x', idx, cnt))

    ! ---- the equality constraints: too many, or with dependent gradients ----
    allocate(eqmap(m))
    m_eq = 0
    do i = 1, m
        eqmap(i) = 0
        if (problem%c_ub(i) - problem%c_lb(i) <= 0.0_wp) then
            m_eq = m_eq + 1
            eqmap(i) = m_eq
        end if
    end do
    if (m_eq > n - n_fixed) call add(text, '     equalities: '//plural(m_eq, 'equality constraint', &
        'equality constraints')//' for '//plural(n - n_fixed, 'free variable', 'free variables')// &
        ': more than can hold in general')
    if (m_eq > 0) then
        ! (each gradient normalized by its largest element, in the variables that aren't fixed)
        allocate(ir(nnz), ic(nnz), vv(nnz), indep(m_eq))
        nz = 0
        do k = 1, nnz
            i = problem%jac_irow(k)
            j = problem%jac_icol(k)
            if (eqmap(i) == 0 .or. jval(k) == 0.0_wp .or. problem%x_ub(j) - problem%x_lb(j) <= 0.0_wp) cycle
            nz = nz + 1
            ir(nz) = j
            ic(nz) = eqmap(i)
            vv(nz) = jval(k)/(rowmax(i)*problem%c_scale(i))
        end do
        call independent_columns(n, m_eq, ir(1:nz), ic(1:nz), vv(1:nz), 1.0e-8_wp, epsilon(1.0_wp)**0.67_wp, indep, stat)
        if (stat == 0) then
            cnt = 0
            do i = 1, m
                if (eqmap(i) == 0) cycle
                if (.not. indep(eqmap(i))) call note(idx, cnt, i)
            end do
            me%found%n_dependent_equalities = cnt
            if (cnt == 1) then
                call add(text, '     equalities: the gradient of '//listed('c', idx, cnt)//' is linearly dependent '// &
                    'on those of the other equality constraints at the starting point (a redundant or '// &
                    'inconsistent constraint)')
            else if (cnt > 1) then
                call add(text, '     equalities: the gradients of '//listed('c', idx, cnt)//' are linearly '// &
                    'dependent on those of the other equality constraints at the starting point (redundant or '// &
                    'inconsistent constraints)')
            end if
        end if
    end if

    if (me%level >= 3) call probe()

    me%found%problem_report = text

    contains

        subroutine probe()
        !! level 3: evaluate `f` and `c` at two points next to the starting point (each variable moved by
        !! `1e-6` of its size, in alternating directions, within its bounds)
        real(wp), parameter :: delta = 1.0e-6_wp
        real(wp), dimension(:), allocatable :: xp, cp
        real(wp) :: fp
        integer  :: side, jj, n_bad, first_bad
        logical  :: f_bad
        allocate(xp(n), cp(m))
        f_bad     = .false.
        n_bad     = 0
        first_bad = 0
        do side = 1, -1, -2
            do jj = 1, n
                xp(jj) = x(jj) + real(side*(1 - 2*mod(jj, 2)), wp)*delta*max(1.0_wp, abs(x(jj)))
                xp(jj) = min(max(xp(jj), problem%x_lb(jj)), problem%x_ub(jj))
            end do
            call problem%f(xp, fp)
            call problem%c(xp, cp)
            if (problem%stop_requested) exit
            if (.not. ieee_is_finite(fp)) f_bad = .true.
            do jj = 1, m
                if (ieee_is_finite(cp(jj))) cycle
                n_bad = n_bad + 1
                if (first_bad == 0) first_bad = jj
            end do
        end do
        me%found%probe_not_finite = f_bad .or. n_bad > 0
        if (me%found%probe_not_finite) then
            line = '     probe: '
            if (f_bad) line = line//'the objective'
            if (f_bad .and. n_bad > 0) line = line//' and '
            if (n_bad > 0) line = line//'c('//fmt_i(first_bad)//')'
            call add(text, line//' could not be evaluated at a point within '//fmt_e(delta)// &
                     ' (relative) of the starting point: it is at the edge of the functions'' domain; '// &
                     'add bounds that keep the variables inside it')
        else
            call add(text, '     probe: the objective and the constraints are finite at two points next to the '// &
                     'starting point')
        end if
        end subroutine probe

    end subroutine diagnostics_start
!*******************************************************************************

!*******************************************************************************
!>
!  the report on the starting point made by [[diagnostics_start]] (empty
!  below level `2`), for the header of the printed output.

    function diagnostics_problem_report(me) result(text)

    class(sqpopt_diagnostics_type), intent(in) :: me
    character(len=:), allocatable :: text

    if (allocated(me%found%problem_report)) then
        text = me%found%problem_report
    else
        text = ''
    end if

    end function diagnostics_problem_report
!*******************************************************************************

!*******************************************************************************
!>
!  note an iterate, with the values and derivatives the solver evaluated
!  there (level `>= 2`): check the changes of the objective and of each
!  constraint since the previous iterate against their gradients at both
!  (see the module documentation), and keep this one for the next check.
!  Calling it again at the same point (after the derivatives were
!  re-evaluated more accurately) just replaces what is kept.

    subroutine diagnostics_point(me, x, f, g, c, jac)

    class(sqpopt_diagnostics_type), intent(inout) :: me
    real(wp), dimension(:),         intent(in)    :: x   !! the iterate `dimension(n)`
    real(wp),                       intent(in)    :: f   !! the objective there
    real(wp), dimension(:),         intent(in)    :: g   !! its gradient `dimension(n)`
    real(wp), dimension(:),         intent(in)    :: c   !! the constraints there `dimension(m)`
    type(sqpopt_sparse_matrix),     intent(in)    :: jac !! their Jacobian

    real(wp), parameter :: max_move = 0.1_wp !! the largest relative change of a variable in a step that is checked
    real(wp) :: pred, move
    integer  :: i, j, k
    logical  :: moved

    if (me%level < 2) return

    if (me%have_point) then
        moved = .false.
        move  = 0.0_wp
        pred  = 0.0_wp
        do j = 1, size(x)
            if (x(j) /= me%x_prev(j)) moved = .true.
            move = max(move, abs(x(j) - me%x_prev(j))/max(1.0_wp, abs(x(j))))
            pred = pred + 0.5_wp*(g(j) + me%g_prev(j))*(x(j) - me%x_prev(j))
        end do
        if (moved .and. move <= max_move) then
            call compare(f - me%f_prev, pred, max(abs(f), abs(me%f_prev)), me%n_checked_f, me%n_off_f)
            me%predicted = 0.0_wp
            do k = 1, jac%nnz
                i = jac%irow(k)
                j = jac%icol(k)
                me%predicted(i) = me%predicted(i) + 0.5_wp*(jac%val(k) + me%jac_prev(k))*(x(j) - me%x_prev(j))
            end do
            do i = 1, size(c)
                call compare(c(i) - me%c_prev(i), me%predicted(i), max(abs(c(i)), abs(me%c_prev(i))), &
                             me%n_checked(i), me%n_off(i))
            end do
        end if
    end if

    me%x_prev   = x
    me%f_prev   = f
    me%g_prev   = g
    me%c_prev   = c
    me%jac_prev = jac%val
    me%have_point = .true.

    contains

        pure subroutine compare(actual, predicted, size_of, n_checked, n_off)
        !! count a change that is large enough to tell, and whether it disagrees with its prediction
        real(wp), intent(in)    :: actual    !! the change of the function
        real(wp), intent(in)    :: predicted !! the change its gradients predict
        real(wp), intent(in)    :: size_of   !! the size of the function's values
        integer,  intent(inout) :: n_checked !! changes checked so far
        integer,  intent(inout) :: n_off     !! of which, off by more than 20%
        real(wp) :: big
        big = max(abs(actual), abs(predicted))
        if (.not. big > 1.0e-9_wp*max(1.0_wp, size_of)) return   ! (lost in the rounding of the values)
        n_checked = n_checked + 1
        if (abs(actual - predicted) > 0.2_wp*big) n_off = n_off + 1
        end subroutine compare

    end subroutine diagnostics_point
!*******************************************************************************

!*******************************************************************************
!>
!  record an iteration (level `>= 2`): its error and step length, for the
!  rate of convergence; the changes of the QP's working set; and its time.
!  With a history file, also write its line.

    subroutine diagnostics_record(me, iter, info, iter_istat, problem, qp_solver, options, time_factorization)

    class(sqpopt_diagnostics_type), intent(inout) :: me
    integer,                        intent(in)    :: iter       !! the major iteration
    type(sqpopt_iter_info),         intent(in)    :: info       !! what happened in it
    integer,                        intent(in)    :: iter_istat !! its status code
    type(sqpopt_problem_type),      intent(in)    :: problem    !! the problem (for its counts and scale factors)
    type(sqpopt_qp_solver_type),    intent(in)    :: qp_solver  !! the QP solver (for its working set and time)
    type(sqpopt_options_type),      intent(in)    :: options    !! the solver options
    real(wp),                       intent(in)    :: time_factorization !! the time in the sparse factorizations
                                                                        !! so far (seconds)

    integer, dimension(:), allocatable :: set
    integer(int64) :: now, rate
    real(wp) :: e, dt
    real(wp), dimension(3) :: part
    integer  :: changes, ios
    character(len=:), allocatable :: line

    if (me%level < 2) return

    ! the error, relative to the tolerances:
    e = max(info%kkt/options%ktol, info%feas/options%ctol)
    if (ieee_is_finite(e)) then
        me%n_err = me%n_err + 1
        me%err(mod(me%n_err - 1, n_window) + 1) = e
    end if

    changes = -1
    if (info%stepped) then
        me%n_alpha = me%n_alpha + 1
        me%alpha(mod(me%n_alpha - 1, n_window) + 1) = info%alpha

        ! the QP's working set: how much of it changed, and whether it is back
        ! where it was two iterations ago
        call qp_solver%working_set(problem%n, problem%m, set)
        if (allocated(set)) then
            if (allocated(me%set1)) then
                changes = count(set /= me%set1)
                if (changes > 0 .and. allocated(me%set2)) then
                    if (all(set == me%set2)) then
                        me%n_flips = me%n_flips + 1
                    else
                        me%n_flips = 0
                    end if
                else
                    me%n_flips = 0
                end if
            end if
            call move_alloc(me%set1, me%set2)
            call move_alloc(set, me%set1)
        end if
    end if

    ! the time of the iteration, and where it went:
    call system_clock(now, rate)
    dt   = real(now - me%clock, wp)/real(rate, wp)
    part(1) = problem%time_user - me%t_user
    part(2) = qp_solver%time - me%t_qp
    part(3) = time_factorization - me%t_fact
    me%clock  = now
    me%t_user = problem%time_user
    me%t_qp   = qp_solver%time
    me%t_fact = time_factorization
    me%t_total = me%t_total + dt
    if (dt > me%t_slowest) then
        me%t_slowest    = dt
        me%slowest      = iter
        me%part_slowest = part
    end if

    if (me%unit == -1) return
    line = fmt_i(iter)//','//csv(info%f/problem%f_scale)//','//csv(info%feas)//','//csv(info%kkt)//','// &
           csv(info%stat_unscaled)//','
    if (info%stepped) then
        line = line//csv(info%alpha)//','//csv(info%step_norm)//','//fmt_i(problem%n_eval_fc)//','// &
               fmt_i(info%n_fc)//','//fmt_i(info%qp_iter)//','//fmt_i(info%qp_istat)//','//csv(info%lam_max)//','// &
               csv(info%glob)//','//csv(info%hess_measure)//','
        if (changes >= 0) line = line//fmt_i(changes)
    else
        ! (the final point: no step was taken from it)
        line = line//',,'//fmt_i(problem%n_eval_fc)//','//fmt_i(info%n_fc)//',,,,,,'
    end if
    line = line//','//csv(dt)//','//csv(part(1))//','//csv(part(2))//','//csv(part(3))//','// &
           trim(sqpopt_iteration_flags(info, iter_istat))
    write(me%unit, '(A)', iostat=ios) line

    end subroutine diagnostics_record
!*******************************************************************************

!*******************************************************************************
!>
!  diagnose the solve at its final point (level `>= 1`), and set
!  `results%diagnosis`, including the text of its report. `results` must
!  hold the outcome of the solve (the status, the counts, and `x`, `c`,
!  `lambda`, and `z` of the original problem). The objective gradient is
!  taken from the problem's evaluation cache (the solver has just evaluated
!  it at `x`), so nothing is evaluated.

    subroutine diagnostics_finish(me, problem, options, results, x, lambda, c, jac, qp_solver, hessian_shift, &
                                  kkt_singular)

    class(sqpopt_diagnostics_type), intent(inout) :: me
    type(sqpopt_problem_type),      intent(inout) :: problem   !! the problem (scaled)
    type(sqpopt_options_type),      intent(in)    :: options   !! the solver options
    type(sqpopt_results_type),      intent(inout) :: results   !! the results of the solve (its `diagnosis` is set)
    real(wp), dimension(:),         intent(in)    :: x         !! the final point `dimension(n)`
    real(wp), dimension(:),         intent(in)    :: lambda    !! the multipliers of the scaled problem `dimension(m)`
    real(wp), dimension(:),         intent(in)    :: c         !! the constraints of the scaled problem at `x`
                                                               !! `dimension(m)`
    type(sqpopt_sparse_matrix),     intent(in)    :: jac       !! the Jacobian of the scaled problem at `x`
    type(sqpopt_qp_solver_type),    intent(in)    :: qp_solver !! the QP solver (after its last solve)
    real(wp),                       intent(in)    :: hessian_shift !! the shift of the Hessian at the end
    logical,                        intent(in)    :: kkt_singular  !! whether the last factorization of a KKT
                                                                   !! matrix (if any) found it singular

    real(wp), dimension(:), allocatable :: g, r, wn, wm, vv, slack
    integer,  dimension(:), allocatable :: cside, xside, ir, ic, colmap, slack_row
    logical,  dimension(:), allocatable :: indep
    character(len=:), allocatable :: text, line
    real(wp) :: tol, lam_big, tol_mult, lb, ub, rate, e_last, slope
    integer  :: n, m, i, j, k, nz, n_act, n_free, stat, n_fit
    logical  :: failed, have_g

    associate (d => results%diagnosis)

    d = me%found
    d%level = me%level
    if (me%level < 1) return

    n = problem%n
    m = problem%m
    failed = results%istat > sqpopt_stalled
    allocate(g(n), r(n), wn(n), wm(m), cside(m), xside(n))

    ! ---- the stationarity residual of each variable: the gradient of the
    ! Lagrangian, projected onto the bounds as in the convergence test ----
    call problem%g(x, g)
    have_g = sqpopt_all_finite(g) .and. sqpopt_all_finite(lambda) .and. sqpopt_all_finite(jac%val)
    wn = 0.0_wp
    if (have_g) then
        call sparse_matvec_transpose(jac, lambda, r)
        do j = 1, n
            r(j) = g(j) - r(j)
            if (problem%x_ub(j) - problem%x_lb(j) <= options%ctol) then
                r(j) = 0.0_wp
            else if (x(j) - problem%x_lb(j) <= options%ctol) then
                r(j) = min(r(j), 0.0_wp)
            else if (problem%x_ub(j) - x(j) <= options%ctol) then
                r(j) = max(r(j), 0.0_wp)
            end if
            wn(j) = abs(r(j))/problem%f_scale
        end do
    end if
    call largest(wn, d%stationarity_variables, d%stationarity_residuals)

    ! ---- the violation of each constraint, and the multipliers ----
    do i = 1, m
        wm(i) = (max(problem%c_lb(i) - c(i), 0.0_wp) + max(c(i) - problem%c_ub(i), 0.0_wp))/problem%c_scale(i)
    end do
    call largest(wm, d%violated_constraints, d%violations)
    do i = 1, m
        wm(i) = abs(results%lambda(i))
    end do
    call largest(wm, d%multiplier_constraints, d%multipliers)
    do i = 1, size(d%multipliers)
        d%multipliers(i) = results%lambda(d%multiplier_constraints(i))
    end do

    ! ---- the active set: which constraints and bounds are active, with what multipliers ----
    tol = max(options%ctol, 1.0e-8_wp)
    lam_big = 1.0_wp
    do i = 1, m
        if (ieee_is_finite(results%lambda(i))) lam_big = max(lam_big, abs(results%lambda(i)))
    end do
    do j = 1, n
        if (ieee_is_finite(results%z(j))) lam_big = max(lam_big, abs(results%z(j)))
    end do
    tol_mult = 1.0e-8_wp*lam_big
    do i = 1, m
        lb = problem%c_lb(i)/problem%c_scale(i)
        ub = problem%c_ub(i)/problem%c_scale(i)
        cside(i) = 0
        if (ub - lb <= 0.0_wp) then
            cside(i) = 2
        else if (abs(results%c(i) - lb) <= tol*max(1.0_wp, abs(lb))) then
            cside(i) = -1
        else if (abs(results%c(i) - ub) <= tol*max(1.0_wp, abs(ub))) then
            cside(i) = 1
        end if
        if (cside(i) /= 0) d%n_active_constraints = d%n_active_constraints + 1
        if (abs(cside(i)) == 1) call classify(cside(i), results%lambda(i))
    end do
    n_free = 0
    do j = 1, n
        xside(j) = 0
        if (problem%x_ub(j) - problem%x_lb(j) <= 0.0_wp) then
            xside(j) = 2
        else if (abs(results%x(j) - problem%x_lb(j)) <= tol*max(1.0_wp, abs(problem%x_lb(j)))) then
            xside(j) = -1
        else if (abs(results%x(j) - problem%x_ub(j)) <= tol*max(1.0_wp, abs(problem%x_ub(j)))) then
            xside(j) = 1
        end if
        if (xside(j) == 0) then
            n_free = n_free + 1
        else
            d%n_active_bounds = d%n_active_bounds + 1
        end if
        if (abs(xside(j)) == 1) call classify(xside(j), results%z(j))
    end do

    ! ---- the active constraints whose gradients are linearly dependent, in the
    ! free variables (each gradient normalized by its largest element there) ----
    n_act = d%n_active_constraints
    if (n_act > 0 .and. sqpopt_all_finite(jac%val)) then
        allocate(colmap(m), ir(jac%nnz), ic(jac%nnz), vv(jac%nnz), indep(n_act))
        k = 0
        do i = 1, m
            colmap(i) = 0
            wm(i) = 0.0_wp
            if (cside(i) /= 0) then
                k = k + 1
                colmap(i) = k
            end if
        end do
        do k = 1, jac%nnz
            if (xside(jac%icol(k)) == 0) wm(jac%irow(k)) = max(wm(jac%irow(k)), abs(jac%val(k)))
        end do
        nz = 0
        do k = 1, jac%nnz
            i = jac%irow(k)
            j = jac%icol(k)
            if (colmap(i) == 0 .or. xside(j) /= 0 .or. jac%val(k) == 0.0_wp) cycle
            nz = nz + 1
            ir(nz) = j
            ic(nz) = colmap(i)
            vv(nz) = jac%val(k)/wm(i)
        end do
        call independent_columns(n, n_act, ir(1:nz), ic(1:nz), vv(1:nz), 1.0e-8_wp, epsilon(1.0_wp)**0.67_wp, indep, stat)
        if (stat == 0) d%n_dependent = n_act - count(indep)
    else if (n_act == 0) then
        d%n_dependent = 0
    end if

    ! ---- the elastic constraints of the last QP ----
    call qp_solver%elastic_slacks(n, slack_row, slack)
    do k = 1, size(slack)
        slack(k) = slack(k)/problem%c_scale(slack_row(k))
    end do
    call largest(slack, d%elastic_constraints, d%elastic_slacks)
    do k = 1, size(d%elastic_constraints)
        d%elastic_constraints(k) = slack_row(d%elastic_constraints(k))
    end do

    ! ---- the history of the iterations (level >= 2) ----
    rate   = 0.0_wp
    e_last = 0.0_wp
    if (me%level >= 2) then
        call fit_rate(n_fit, slope, e_last)
        if (n_fit >= 5) then
            rate = exp(slope)
            d%convergence_rate = rate
            if (rate < 1.0_wp .and. e_last > 1.0_wp) then
                if (log(e_last)/(-slope) < 1.0e6_wp) d%iterations_needed = ceiling(log(e_last)/(-slope))
            end if
        end if
        d%n_active_set_flips = me%n_flips
        d%objective_derivative_suspect = me%n_off_f >= 3 .and. 2*me%n_off_f >= me%n_checked_f
        do i = 1, m
            wm(i) = 0.0_wp
            if (me%n_off(i) >= 3 .and. 2*me%n_off(i) >= me%n_checked(i)) wm(i) = real(me%n_off(i), wp)/me%n_checked(i)
        end do
        call largest(wm, d%derivative_suspects, vv)
        d%slowest_iteration      = me%slowest
        d%slowest_iteration_time = me%t_slowest
    end if

    ! ---- the report ----
    text = ''
    call add(text, ' sqpopt diagnosis (diagnostic_level = '//fmt_i(me%level)//'):')
    call verdict()

    if (size(d%violations) > 0) then
        if (d%violations(1) > options%ctol) then
            line = '   largest violations: '
            do k = 1, size(d%violations)
                i = d%violated_constraints(k)
                if (k > 1) line = line//'; '
                line = line//'c('//fmt_i(i)//') by '//fmt_e(d%violations(k))//' (value '//fmt_e(results%c(i))// &
                       ', bounds ['//bound(problem%c_lb(i)/problem%c_scale(i))//', '// &
                       bound(problem%c_ub(i)/problem%c_scale(i))//'])'
            end do
            call add(text, line)
        end if
    end if
    if (size(d%stationarity_residuals) > 0 .and. results%istat /= sqpopt_success .and. &
        results%istat /= sqpopt_infeasible .and. results%istat /= sqpopt_function_error) then
        line = '   largest stationarity residuals: '
        do k = 1, size(d%stationarity_residuals)
            j = d%stationarity_variables(k)
            if (k > 1) line = line//'; '
            line = line//'x('//fmt_i(j)//') '//fmt_e(d%stationarity_residuals(k))
            if (xside(j) == -1) line = line//' (at its lower bound)'
            if (xside(j) == 1)  line = line//' (at its upper bound)'
        end do
        call add(text, line)
    end if
    if (size(d%multipliers) > 0 .and. (failed .or. lam_big > badly_scaled)) then
        line = '   largest multipliers: '
        do k = 1, size(d%multipliers)
            if (k > 1) line = line//'; '
            line = line//'c('//fmt_i(d%multiplier_constraints(k))//') '//fmt_e(d%multipliers(k))
        end do
        call add(text, line)
    end if

    line = '   active set: '//fmt_i(d%n_active_constraints)//' of '//plural(m, 'constraint', 'constraints')//' and '// &
           fmt_i(d%n_active_bounds)//' of '//plural(n, 'variable bound', 'variable bounds')//' active, '// &
           plural(n_free, 'free variable', 'free variables')
    if (d%n_weakly_active > 0) line = line//'; '//fmt_i(d%n_weakly_active)//' weakly active (multiplier near zero)'
    if (d%n_wrong_sign > 0) line = line//'; '//fmt_i(d%n_wrong_sign)//' with a multiplier of the wrong sign'
    if (d%n_dependent > 0) line = line//'; '//plural(d%n_dependent, 'active gradient', 'active gradients')// &
                                  ' linearly dependent on the others'
    call add(text, line)
    if (d%n_dependent > 0 .or. d%n_weakly_active > 0) call add(text, '     (a degenerate point: the multipliers '// &
        'may not be unique, and convergence near it can be slow)')

    if (results%n_qp_solves > 0 .and. (failed .or. qp_solver%n_slacks > 0)) then
        line = '   last QP: '//qp_solver%mode_name(n)
        if (qp_solver%unconstrained_used) then
            line = line//', solved by the unconstrained step'
        else if (qp_solver%direct_used) then
            line = line//', solved directly'
        else
            if (qp_solver%daqp_fallback) then
                line = line//', solved by the dense QP solver (DAQP: '//daqp_status_text(qp_solver%daqp_qp%status)//')'
            else if (qp_solver%solver_name(n) /= qp_solver%mode_name(n)) then
                line = line//', solved by the '//qp_solver%solver_name(n)//' solver (an elastic re-solve)'
            end if
            line = line//', '//plural(qp_solver%n_iter, 'iteration', 'iterations')
        end if
        line = line//', working set '//fmt_i(qp_solver%n_working)//' (of '//fmt_i(n)//' variables)'
        if (qp_solver%negative_curvature) line = line//', negative curvature'
        if (hessian_shift > 0.0_wp) line = line//', Hessian shift '//fmt_e(hessian_shift)
        if (kkt_singular) line = line//', singular KKT matrix'
        call add(text, line)
        if (qp_solver%n_slacks > 0) then
            line = '     its linearized constraints were relaxed ('// &
                   plural(qp_solver%n_slacks, 'elastic constraint', 'elastic constraints')//')'
            do k = 1, size(d%elastic_slacks)
                line = line//merge(':', ';', k == 1)//' c('//fmt_i(d%elastic_constraints(k))//') by '// &
                       fmt_e(d%elastic_slacks(k))
            end do
            call add(text, line)
        end if
    end if

    if (me%level >= 2) then
        if (rate > 0.0_wp .and. results%istat /= sqpopt_success) then
            line = '   convergence: the error is '//fmt_e(e_last)//' times its tolerance, and changed by a factor '// &
                   fmt_e(rate)//' per iteration over the last '//fmt_i(n_fit)
            if (d%iterations_needed >= 0) line = line//': about '// &
                plural(d%iterations_needed, 'more iteration', 'more iterations')//' at this rate'
            call add(text, line)
        end if
        call steps_line()
        if (me%n_flips >= 3) call add(text, '   working set: the QP''s working set went back and forth between the '// &
            'same two states in the last '//fmt_i(me%n_flips)//' iterations (constraints entering and leaving in turn)')
        if (d%objective_derivative_suspect) call add(text, '   derivatives: the objective''s changes disagreed with '// &
            'its gradient in '//fmt_i(me%n_off_f)//' of '//plural(me%n_checked_f, 'step', 'steps')// &
            ': check the gradient (or the objective isn''t smooth along the way)')
        if (size(d%derivative_suspects) > 0) then
            line = '   derivatives: changes that disagreed with the Jacobian row: '
            do k = 1, size(d%derivative_suspects)
                i = d%derivative_suspects(k)
                if (k > 1) line = line//'; '
                line = line//'c('//fmt_i(i)//') in '//fmt_i(me%n_off(i))//' of '//plural(me%n_checked(i), 'step', 'steps')
            end do
            call add(text, line//': check these derivatives (or the constraints aren''t smooth along the way)')
        end if
        if (me%slowest > 0 .and. me%t_total > 0.0_wp .and. results%iterations >= 5) then
            if (me%t_slowest > 0.3_wp*me%t_total) then
                call add(text, '   time: iteration '//fmt_i(me%slowest)//' took '//fmt_e(me%t_slowest)//' s, '// &
                    fmt_i(nint(100.0_wp*me%t_slowest/me%t_total))//'% of the iterations'' time (user functions '// &
                    fmt_e(me%part_slowest(1))//' s, QP '//fmt_e(me%part_slowest(2))//' s, factorizations '// &
                    fmt_e(me%part_slowest(3))//' s)')
            end if
        end if
    end if

    d%report = text

    end associate

    contains

        subroutine classify(side, mult)
        !! count an active inequality constraint or bound whose multiplier is near zero, or has the wrong sign
        integer,  intent(in) :: side !! `-1`: at its lower bound, `+1`: at its upper bound
        real(wp), intent(in) :: mult !! its multiplier (`>= 0` at a lower bound, `<= 0` at an upper bound)
        if (.not. ieee_is_finite(mult)) return
        if (abs(mult) <= tol_mult) then
            results%diagnosis%n_weakly_active = results%diagnosis%n_weakly_active + 1
        else if (real(side, wp)*mult > 0.0_wp) then
            results%diagnosis%n_wrong_sign = results%diagnosis%n_wrong_sign + 1
        end if
        end subroutine classify

        subroutine fit_rate(n_fit, slope, e_last)
        !! the least-squares slope of the logarithm of the error over the last iterations recorded
        integer,  intent(out) :: n_fit  !! the iterations used (at most `n_window`)
        real(wp), intent(out) :: slope  !! the slope (the logarithm of the factor per iteration)
        real(wp), intent(out) :: e_last !! the error of the last one
        real(wp) :: sx, sy, sxx, sxy, y, xk
        integer  :: kk, first
        n_fit  = 0
        slope  = 0.0_wp
        e_last = 0.0_wp
        if (me%n_err == 0) return
        e_last = me%err(mod(me%n_err - 1, n_window) + 1)
        first = max(1, me%n_err - n_window + 1)
        sx = 0.0_wp; sy = 0.0_wp; sxx = 0.0_wp; sxy = 0.0_wp
        do kk = first, me%n_err
            y = me%err(mod(kk - 1, n_window) + 1)
            if (.not. y > 0.0_wp) cycle
            n_fit = n_fit + 1
            xk  = real(kk - first, wp)
            y   = log(y)
            sx  = sx + xk
            sy  = sy + y
            sxx = sxx + xk*xk
            sxy = sxy + xk*y
        end do
        if (n_fit >= 2) then
            if (n_fit*sxx - sx*sx > 0.0_wp) slope = (n_fit*sxy - sx*sy)/(n_fit*sxx - sx*sx)
        end if
        end subroutine fit_rate

        subroutine steps_line()
        !! the report's line on the step lengths of the last iterations, if many of them were shortened
        integer  :: kk, first, n_steps, n_short
        real(wp) :: a_min
        if (me%n_alpha == 0) return
        first   = max(1, me%n_alpha - n_window + 1)
        n_steps = me%n_alpha - first + 1
        n_short = 0
        a_min   = 1.0_wp
        do kk = first, me%n_alpha
            if (me%alpha(mod(kk - 1, n_window) + 1) < 1.0_wp) n_short = n_short + 1
            a_min = min(a_min, me%alpha(mod(kk - 1, n_window) + 1))
        end do
        if (2*n_short >= n_steps .and. n_short >= 3) call add(text, '   steps: '//fmt_i(n_short)//' of the last '// &
            fmt_i(n_steps)//' steps were shortened by the line search (the shortest to '//fmt_e(a_min)// &
            ' of the QP''s step): the QP''s model predicts the functions poorly here')
        end subroutine steps_line

        subroutine verdict()
        !! the report's lines on what the final status most likely means
        integer :: kk, n_at_bound, iv
        logical :: in_violated
        associate (d => results%diagnosis)
        select case (results%istat)
        case (sqpopt_success)
            call add(text, '   converged: the KKT conditions hold to the tolerances')
        case (sqpopt_acceptable)
            call add(text, '   stopped at the acceptable tolerances: the KKT error ('//fmt_e(results%kkt_error)// &
                     ') stayed above ktol ('//fmt_e(options%ktol)//') for acceptable_iter iterations')
            call add(text, '     what usually prevents a tighter convergence: inaccurate derivatives, or poor scaling')
        case (sqpopt_stalled)
            call add(text, '   stalled: the point is feasible and the steps have become negligible, but the KKT '// &
                     'error is '//fmt_e(results%kkt_error)//' (ktol '//fmt_e(options%ktol)//')')
            call add(text, '     the usual causes: inaccurate derivatives (the line search then finds no better '// &
                     'point), a degenerate solution, or poor scaling')
        case (sqpopt_max_iter_reached, sqpopt_max_evals_reached, sqpopt_time_limit_reached)
            if (rate > 0.0_wp .and. rate < 1.0_wp) then
                call add(text, '   the limit was reached while the solver was still converging: raise it')
            else if (rate >= 1.0_wp) then
                call add(text, '   the limit was reached, and the error was no longer decreasing: raising the limit '// &
                         'alone is unlikely to help')
            else if (me%level < 2) then
                call add(text, '   the limit was reached (diagnostic_level = 2 tells whether the solver was still '// &
                         'converging)')
            else
                call add(text, '   the limit was reached after too few iterations to tell the rate of convergence')
            end if
            if (4*results%n_hessian_resets >= results%iterations .and. results%n_hessian_resets >= 3) then
                call add(text, '     the Hessian was reset (or its shift raised) in '// &
                         fmt_i(results%n_hessian_resets)//' of '//fmt_i(results%iterations)// &
                         ' iterations: the problem is nonconvex or poorly scaled along the way'// &
                         ' (with the exact Hessian, options%inertia_control finds the shift directly)')
            end if
            if (2*results%n_restoration_steps >= results%iterations .and. results%n_restoration_steps >= 3) then
                call add(text, '     '//fmt_i(results%n_restoration_steps)//' of '//fmt_i(results%iterations)// &
                         ' iterations were restoration steps: the solver spent them getting feasible'// &
                         ' (a starting point nearer the constraints would help)')
            end if
        case (sqpopt_infeasible)
            call add(text, '   infeasible here: the point is stationary for the constraint violation, so the '// &
                     'violated constraints below can''t be satisfied together near it')
            ! (how many of the variables of the violated constraints are at a bound)
            n_at_bound = 0
            wn = 0.0_wp   ! (marks the variables counted)
            do kk = 1, jac%nnz
                if (xside(jac%icol(kk)) == 0 .or. jac%val(kk) == 0.0_wp .or. wn(jac%icol(kk)) /= 0.0_wp) cycle
                in_violated = .false.
                do iv = 1, size(d%violated_constraints)
                    if (d%violated_constraints(iv) == jac%irow(kk)) in_violated = .true.
                end do
                if (in_violated) then
                    n_at_bound = n_at_bound + 1
                    wn(jac%icol(kk)) = 1.0_wp
                end if
            end do
            if (n_at_bound > 0) call add(text, '     '//plural(n_at_bound, 'variable of these constraints is', &
                'variables of these constraints are')//' at a bound, which may be what prevents it')
            call add(text, '     the problem may have no feasible point, or this is a local minimum of the '// &
                     'violation: check the constraints'' bounds, or try another starting point')
        case (sqpopt_line_search_failed)
            call add(text, '   no acceptable step was found along the QP''s direction in '// &
                     fmt_i(options%max_consecutive_failures)//' consecutive iterations')
            if (me%level >= 2 .and. .not. (d%objective_derivative_suspect .or. size(d%derivative_suspects) > 0)) then
                call add(text, '     the derivatives agreed with the functions'' changes along the steps taken, so '// &
                         'the likely causes are a function that isn''t smooth here, or poor scaling')
            else if (me%level < 2) then
                call add(text, '     the usual cause is a wrong or inaccurate derivative (diagnostic_level = 2 '// &
                         'checks them against the steps taken), or a function that isn''t smooth here')
            end if
        case (sqpopt_qp_solve_failed)
            call add(text, '   the QP subproblem solver failed in '//fmt_i(options%max_consecutive_failures)// &
                     ' consecutive iterations (its iteration limit, or a QP that is unbounded below)')
            call add(text, '     see the last QP below; a QP with negative curvature needs a larger Hessian shift'// &
                     ' (options%inertia_control) or a positive definite Hessian (options%hessian_mode)')
        case (sqpopt_function_error)
            line = '   a function value is not finite at the current point:'
            if (.not. ieee_is_finite(results%f)) line = line//' the objective'
            do kk = 1, m
                if (ieee_is_finite(results%c(kk))) cycle
                line = line//' c('//fmt_i(kk)//')'
                exit
            end do
            if (.not. have_g) line = line//' (or a derivative)'
            call add(text, line)
            call add(text, '     add bounds that keep the variables inside the functions'' domain, or return '// &
                     'status > 0 from fc there, so that the solver backs off (gjac and hess are only called at '// &
                     'accepted points, so their status > 0 ends the solve)')
        case (sqpopt_unbounded)
            call add(text, '   the objective fell below options%obj_lower_limit at a feasible point: the problem '// &
                     'is unbounded below, or needs bounds or constraints that are missing')
        case (sqpopt_out_of_memory)
            call add(text, '   an array could not be allocated: for a large problem use the sparse QP solver '// &
                     '(options%qp_solver_mode)')
        case (sqpopt_user_requested_stop)
            call add(text, '   stopped at the user''s request')
        case default
            call add(text, '   status '//fmt_i(results%istat))
        end select

        ! a very large multiplier, at a point that isn't a solution
        if (failed .and. size(d%multipliers) > 0) then
            if (abs(d%multipliers(1)) > badly_scaled*max(1.0_wp, maxval(abs(g))/problem%f_scale)) then
                call add(text, '     the multiplier of c('//fmt_i(d%multiplier_constraints(1))//') is very large ('// &
                         fmt_e(d%multipliers(1))//'): the constraint qualification probably fails near this point'// &
                         ' (dependent or redundant constraints)')
            end if
        end if
        end associate
        end subroutine verdict

    end subroutine diagnostics_finish
!*******************************************************************************

!*******************************************************************************
!>
!  write a report (the lines of `text`, each ended by a newline character)
!  to `unit`. As all the printed output, it never stops the solver.

    subroutine sqpopt_diagnostics_write(unit, text)

    integer,          intent(in) :: unit !! the unit to write to
    character(len=*), intent(in) :: text !! the report

    integer :: first, last, ios

    first = 1
    do while (first <= len(text))
        last = index(text(first:), nl)
        if (last == 0) then
            write(unit, '(A)', iostat=ios) text(first:)
            exit
        end if
        write(unit, '(A)', iostat=ios) text(first:first+last-2)
        first = first + last
    end do

    end subroutine sqpopt_diagnostics_write
!*******************************************************************************

!*******************************************************************************
!>
!  add a line to a report. A line longer than `line_width` is continued on
!  further lines, broken at spaces, and indented two more than the first.

    subroutine add(text, line)

    character(len=:), allocatable, intent(inout) :: text !! the report so far
    character(len=*),              intent(in)    :: line !! the line to add

    integer :: first, last, indent, width, lead

    indent = verify(line, ' ') - 1
    if (indent < 0) indent = 0
    first = 1
    width = line_width
    lead  = indent + 1   ! (a space of the indentation is not a place to break the line)
    do
        if (len(line) - first + 1 <= width) then
            text = text//line(first:)//nl
            exit
        end if
        ! (the last space that keeps this piece within the width; if there is none, the first one after it)
        last = index(line(first:first+width), ' ', back=.true.)
        if (last <= lead) then
            last = index(line(first+width:), ' ')
            if (last == 0) then
                text = text//line(first:)//nl
                exit
            end if
            last = last + width
        end if
        text = text//line(first:first+last-2)//nl//repeat(' ', indent + 2)
        first = first + last
        width = line_width - indent - 2
        lead  = 1
    end do

    end subroutine add
!*******************************************************************************

!*******************************************************************************
!>
!  add an index to a list of at most `n_listed`, and count it.

    pure subroutine note(idx, cnt, i)

    integer, dimension(n_listed), intent(inout) :: idx !! the first indices noted
    integer,                      intent(inout) :: cnt !! how many were noted in all
    integer,                      intent(in)    :: i   !! the index to note

    cnt = cnt + 1
    if (cnt <= n_listed) idx(cnt) = i

    end subroutine note
!*******************************************************************************

!*******************************************************************************
!>
!  the indices of a list made by [[note]], as text: `c(3), c(7) and 4 more`.

    pure function listed(name, idx, cnt) result(str)

    character(len=*),             intent(in) :: name !! what they index: `c` or `x`
    integer, dimension(n_listed), intent(in) :: idx  !! the first indices
    integer,                      intent(in) :: cnt  !! how many there are in all
    character(len=:), allocatable :: str

    integer :: k

    str = ''
    do k = 1, min(cnt, n_listed)
        if (k > 1) str = str//', '
        str = str//name//'('//fmt_i(idx(k))//')'
    end do
    if (cnt > n_listed) str = str//' and '//fmt_i(cnt - n_listed)//' more'

    end function listed
!*******************************************************************************

!*******************************************************************************
!>
!  the indices of the smallest and of the largest positive element of `v`
!  (`0` if it has none).

    pure subroutine extremes(v, i_min, i_max)

    real(wp), dimension(:), intent(in)  :: v     !! the values (none negative)
    integer,                intent(out) :: i_min !! index of the smallest positive one
    integer,                intent(out) :: i_max !! index of the largest one

    integer :: i

    i_min = 0
    i_max = 0
    do i = 1, size(v)
        if (.not. v(i) > 0.0_wp) cycle
        if (i_min == 0) then
            i_min = i
            i_max = i
        else
            if (v(i) < v(i_min)) i_min = i
            if (v(i) > v(i_max)) i_max = i
        end if
    end do

    end subroutine extremes
!*******************************************************************************

!*******************************************************************************
!>
!  the largest positive (and finite) elements of `v`, at most `n_listed` of
!  them, largest first: their indices and their values.

    pure subroutine largest(v, idx, val)

    real(wp), dimension(:),              intent(in)  :: v   !! the values
    integer,  dimension(:), allocatable, intent(out) :: idx !! the indices of the largest ones
    real(wp), dimension(:), allocatable, intent(out) :: val !! their values

    integer,  dimension(n_listed) :: best
    real(wp), dimension(n_listed) :: best_v
    integer :: i, j, k

    k = 0
    do i = 1, size(v)
        if (.not. ieee_is_finite(v(i))) cycle
        if (.not. v(i) > 0.0_wp) cycle
        if (k < n_listed) then
            k = k + 1
        else if (.not. v(i) > best_v(k)) then
            cycle
        end if
        j = k
        do while (j > 1)
            if (best_v(j-1) >= v(i)) exit
            best_v(j) = best_v(j-1)
            best(j)   = best(j-1)
            j = j - 1
        end do
        best_v(j) = v(i)
        best(j)   = i
    end do
    allocate(idx(k), val(k))
    idx = best(1:k)
    val = best_v(1:k)

    end subroutine largest
!*******************************************************************************

!*******************************************************************************
!>
!  a bound, for a report (`-inf` or `inf` if there is none).

    pure function bound(v) result(str)

    real(wp), intent(in) :: v !! the bound
    character(len=:), allocatable :: str

    if (v <= -sqpopt_infinity) then
        str = '-inf'
    else if (v >= sqpopt_infinity) then
        str = 'inf'
    else
        str = fmt_e(v)
    end if

    end function bound
!*******************************************************************************

!*******************************************************************************
!>
!  a real number for the history file (blank if it can't be written).

    pure function csv(v) result(str)

    real(wp), intent(in) :: v !! the number
    character(len=:), allocatable :: str

    character(len=32) :: buf
    integer :: ios

    write(buf, '(ES16.8E3)', iostat=ios) v
    if (ios /= 0) then
        str = ''
    else
        str = trim(adjustl(buf))
    end if

    end function csv
!*******************************************************************************

    end module sqpopt_diagnostics_module
!*******************************************************************************
