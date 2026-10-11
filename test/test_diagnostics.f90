program test_diagnostics

    !! The diagnostics of a solve (`options%diagnostic_level`, see
    !! `sqpopt_diagnostics_module`):
    !!
    !! * they are opt-in, and cost no evaluations up to level 2: with levels
    !!   `0`, `1`, and `2` the solver makes the same calls of `fc` and `gjac`
    !!   and returns the same point (with and without scaling); level `3`
    !!   returns the same point too, with two more calls of `fc`.
    !! * `print_level` overrides them: with `print_level = 0` nothing is
    !!   printed at any diagnostic level, but `results%diagnosis` is set.
    !!   With `print_level = 1`, the report on the starting point and the
    !!   diagnosis are in the output.
    !! * the history file (`diagnostics_unit`) has a header and one line per
    !!   iteration, whatever `print_level` is.
    !! * what they find, on problems made for it: a wrong Jacobian row
    !!   (from the steps taken), inconsistent constraints, an iteration limit
    !!   reached while converging, dependent and single-variable constraints
    !!   and poor scaling at the starting point, a weakly active constraint
    !!   at the solution, a starting point at the edge of the objective's
    !!   domain, and one outside it (where the functions can't be evaluated:
    !!   the diagnosis must cope with the values that are not finite).
    !! * an invalid `diagnostic_level` or `diagnostics_unit` is rejected.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_types_module,   only: sqpopt_results_type, sqpopt_success, sqpopt_infeasible, sqpopt_max_iter_reached, &
                                     sqpopt_invalid_input, sqpopt_function_error
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none
    type(sqpopt_options_type) :: default_options !! default values (default-initialized)

    ! the problems (see `fc`):
    integer, parameter :: p_disk      = 1 !! the Rosenbrock function on the unit disk
    integer, parameter :: p_wrong     = 2 !! the same, with a second constraint (never active) whose Jacobian row
                                          !! is wrong
    integer, parameter :: p_infeas    = 3 !! two constraints that can't both hold
    integer, parameter :: p_structure = 4 !! dependent equalities, a single-variable constraint, and a tiny variable
    integer, parameter :: p_weak      = 5 !! a constraint that is active with a zero multiplier
    integer, parameter :: p_edge      = 6 !! an objective that isn't defined on one side of the starting point
    integer, parameter :: p_outside   = 7 !! the same objective, from a starting point where it isn't defined

    integer :: which !! the problem being solved

    write(*,*) '----------------------------'
    write(*,*) 'test_diagnostics'
    write(*,*) '----------------------------'

    call test_no_cost(scaling=.true.)
    call test_no_cost(scaling=.false.)
    call test_printing()
    call test_history_file()
    call test_wrong_derivative()
    call test_infeasible()
    call test_limit()
    call test_starting_point()
    call test_weakly_active()
    call test_probe()
    call test_function_error()
    call test_invalid()

    print '(A)', 'test_diagnostics PASSED'

    contains

    subroutine setup(problem, x0)
    !! define problem `which`, and its starting point
    type(sqpopt_problem_type),           intent(out) :: problem !! the problem
    real(wp), dimension(:), allocatable, intent(out) :: x0      !! the starting point
    real(wp), parameter :: inf = 1.0e20_wp
    select case (which)
    case (p_disk)
        call problem%set_problem_size(n=2, m=1)
        call problem%set_bounds([-inf, -inf], [inf, inf], [-inf], [1.0_wp])
        call problem%set_jacobian_sparsity(2, [1, 1], [1, 2])
        x0 = [-1.2_wp, 1.0_wp]
    case (p_wrong)
        call problem%set_problem_size(n=2, m=2)
        call problem%set_bounds([-inf, -inf], [inf, inf], [-inf, -inf], [1.0_wp, 100.0_wp])
        call problem%set_jacobian_sparsity(4, [1, 1, 2, 2], [1, 2, 1, 2])
        x0 = [-1.2_wp, 1.0_wp]
    case (p_infeas)
        call problem%set_problem_size(n=2, m=2)
        call problem%set_bounds([-inf, -inf], [inf, inf], [3.0_wp, -inf], [inf, 1.0_wp])
        call problem%set_jacobian_sparsity(4, [1, 1, 2, 2], [1, 2, 1, 2])
        x0 = [0.0_wp, 0.0_wp]
    case (p_structure)
        call problem%set_problem_size(n=3, m=3)
        call problem%set_bounds([-inf, -inf, -inf], [inf, inf, inf], [1.0_wp, 2.0_wp, -inf], [1.0_wp, 2.0_wp, 5.0_wp])
        call problem%set_jacobian_sparsity(5, [1, 1, 2, 2, 3], [1, 2, 1, 2, 3])
        x0 = [0.0_wp, 0.0_wp, 0.0_wp]
    case (p_weak)
        call problem%set_problem_size(n=2, m=1)
        call problem%set_bounds([-inf, -inf], [inf, inf], [-inf], [3.0_wp])
        call problem%set_jacobian_sparsity(2, [1, 1], [1, 2])
        x0 = [0.0_wp, 0.0_wp]
    case (p_edge, p_outside)
        call problem%set_problem_size(n=1, m=0)
        call problem%set_bounds([-inf], [inf], [real(wp) ::], [real(wp) ::])
        call problem%set_jacobian_sparsity(0, [integer ::], [integer ::])
        x0 = [0.0_wp]
        if (which == p_outside) x0 = [-1.0_wp]
    end select
    call problem%set_functions(fc=fc, gjac=gjac)
    end subroutine setup

    subroutine solve(level, print_level, options, r, istat)
    !! solve problem `which` with the given diagnostic and print levels
    integer,                   intent(in)    :: level       !! `options%diagnostic_level`
    integer,                   intent(in)    :: print_level !! `options%print_level`
    type(sqpopt_options_type), intent(inout) :: options     !! the other options (the two levels are set here)
    type(sqpopt_results_type), intent(out)   :: r           !! the results
    integer,                   intent(out)   :: istat       !! the status
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    real(wp), dimension(:), allocatable :: x0
    call setup(problem, x0)
    options%diagnostic_level = level
    options%print_level      = print_level
    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_results(r)
    end subroutine solve

    subroutine test_no_cost(scaling)
    !! levels 0 to 2 make the same evaluations and give the same point; level 3 gives the same point with
    !! two more calls of `fc`
    logical, intent(in) :: scaling !! whether the automatic scaling is on
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r0, r
    integer :: istat, level, p
    integer, dimension(2), parameter :: problems = [p_disk, p_structure]
    do p = 1, size(problems)
        which = problems(p)
        options = default_options
        options%scaling = scaling
        call solve(0, 0, options, r0, istat)
        if (r0%diagnosis%level /= 0 .or. allocated(r0%diagnosis%report)) then
            error stop 'test_diagnostics FAILED: level 0 made a diagnosis'
        end if
        do level = 1, 3
            call solve(level, 0, options, r, istat)
            print '(A,I0,A,L1,A,I0,3(A,I0))', 'problem ', which, ', scaling ', scaling, ', level ', level, &
                ': iterations ', r%iterations, ', fc ', r%n_eval_fc, ', gjac ', r%n_eval_gjac
            if (r%istat /= r0%istat .or. r%iterations /= r0%iterations .or. any(r%x /= r0%x) .or. r%f /= r0%f) then
                error stop 'test_diagnostics FAILED: the diagnostics changed the iterates'
            end if
            if (r%n_eval_gjac /= r0%n_eval_gjac) error stop 'test_diagnostics FAILED: the diagnostics called gjac'
            if (level <= 2 .and. r%n_eval_fc /= r0%n_eval_fc) then
                error stop 'test_diagnostics FAILED: levels 1 and 2 must not call fc'
            end if
            if (level == 3 .and. r%n_eval_fc /= r0%n_eval_fc + 2) then
                error stop 'test_diagnostics FAILED: level 3 must call fc twice more'
            end if
            if (r%diagnosis%level /= level .or. .not. allocated(r%diagnosis%report)) then
                error stop 'test_diagnostics FAILED: no diagnosis in the results'
            end if
            if (allocated(r%diagnosis%problem_report) .neqv. level >= 2) then
                error stop 'test_diagnostics FAILED: the report on the starting point is for levels 2 and 3'
            end if
        end do
    end do
    end subroutine test_no_cost

    subroutine test_printing()
    !! `print_level = 0` prints nothing at any diagnostic level; `print_level = 1` prints both reports
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat, u, ios, n_lines
    logical :: found_start, found_diagnosis
    character(len=1024) :: line
    which = p_disk
    open(newunit=u, status='scratch', action='readwrite', form='formatted')
    options%output_unit = u
    call solve(3, 0, options, r, istat)
    rewind(u)
    read(u, '(A)', iostat=ios) line
    if (ios == 0) error stop 'test_diagnostics FAILED: something was printed with print_level = 0'
    close(u)

    open(newunit=u, status='scratch', action='readwrite', form='formatted')
    options%output_unit = u
    call solve(2, 1, options, r, istat)
    rewind(u)
    found_start     = .false.
    found_diagnosis = .false.
    n_lines = 0
    do
        read(u, '(A)', iostat=ios) line
        if (ios /= 0) exit
        n_lines = n_lines + 1
        if (index(line, 'diagnostics of the starting point') > 0) found_start = .true.
        if (index(line, 'sqpopt diagnosis') > 0) found_diagnosis = .true.
    end do
    close(u)
    print '(A,I0,A,2L2)', 'print_level 1: ', n_lines, ' lines, reports found', found_start, found_diagnosis
    if (.not. (found_start .and. found_diagnosis)) error stop 'test_diagnostics FAILED: a report was not printed'
    print '(A)', r%diagnosis%problem_report//r%diagnosis%report
    end subroutine test_printing

    subroutine test_history_file()
    !! the history file has a header and a line per iteration, also with `print_level = 0`, and only from level 2
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat, u, ios, n_lines, level
    character(len=1024) :: line, header
    which = p_disk
    do level = 1, 2
        open(newunit=u, status='scratch', action='readwrite', form='formatted')
        options%diagnostics_unit = u
        call solve(level, 0, options, r, istat)
        rewind(u)
        header  = ''
        n_lines = 0
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            n_lines = n_lines + 1
            if (n_lines == 1) header = line
        end do
        close(u)
        print '(A,I0,A,I0,A,I0,A)', 'history file, level ', level, ': ', n_lines, ' lines for ', r%iterations, ' iterations'
        if (level == 1 .and. n_lines /= 0) error stop 'test_diagnostics FAILED: a history file at level 1'
        if (level == 2) then
            if (n_lines /= r%iterations + 1) error stop 'test_diagnostics FAILED: lines of the history file'
            if (header(1:20) /= 'iteration,objective,') error stop 'test_diagnostics FAILED: header of the history file'
        end if
    end do
    end subroutine test_history_file

    subroutine test_wrong_derivative()
    !! a Jacobian row that doesn't match its constraint is found from the steps taken (the constraint is
    !! never active, so the solve is the same as without it)
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat
    which = p_wrong
    call solve(2, 0, options, r, istat)
    print '(A)', r%diagnosis%report
    if (size(r%diagnosis%derivative_suspects) /= 1) error stop 'test_diagnostics FAILED: one constraint is suspect'
    if (r%diagnosis%derivative_suspects(1) /= 2) error stop 'test_diagnostics FAILED: the suspect constraint is c(2)'
    if (r%diagnosis%objective_derivative_suspect) error stop 'test_diagnostics FAILED: the objective gradient is right'
    ! (the other problems have the right derivatives)
    which = p_disk
    options = default_options
    call solve(2, 0, options, r, istat)
    if (size(r%diagnosis%derivative_suspects) /= 0 .or. r%diagnosis%objective_derivative_suspect) then
        error stop 'test_diagnostics FAILED: right derivatives were suspected'
    end if
    end subroutine test_wrong_derivative

    subroutine test_infeasible()
    !! inconsistent constraints: the violated ones are named
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat
    which = p_infeas
    call solve(1, 0, options, r, istat)
    print '(A)', r%diagnosis%report
    if (istat /= sqpopt_infeasible) error stop 'test_diagnostics FAILED: the problem should be infeasible'
    if (size(r%diagnosis%violated_constraints) /= 2) error stop 'test_diagnostics FAILED: both constraints are violated'
    if (any(abs(r%diagnosis%violations - 1.0_wp) > 1.0e-6_wp)) error stop 'test_diagnostics FAILED: the violations are 1'
    if (index(r%diagnosis%report, 'infeasible here') == 0) error stop 'test_diagnostics FAILED: no infeasibility verdict'
    end subroutine test_infeasible

    subroutine test_limit()
    !! an iteration limit reached while converging: the rate is below 1
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat
    which = p_disk
    options%max_iter = 9
    call solve(2, 0, options, r, istat)
    print '(A)', r%diagnosis%report
    if (istat /= sqpopt_max_iter_reached) error stop 'test_diagnostics FAILED: the iteration limit should be reached'
    if (.not. (r%diagnosis%convergence_rate > 0.0_wp .and. r%diagnosis%convergence_rate < 1.0_wp)) then
        error stop 'test_diagnostics FAILED: the solver was converging'
    end if
    if (r%diagnosis%iterations_needed < 1) error stop 'test_diagnostics FAILED: no estimate of the iterations needed'
    if (r%diagnosis%slowest_iteration < 1) error stop 'test_diagnostics FAILED: no slowest iteration'
    end subroutine test_limit

    subroutine test_starting_point()
    !! the report on the starting point: dependent equalities, a single-variable constraint, poor scaling
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat
    which = p_structure
    call solve(2, 0, options, r, istat)
    print '(A)', r%diagnosis%problem_report
    if (r%diagnosis%n_dependent_equalities /= 1) error stop 'test_diagnostics FAILED: one equality is dependent'
    if (r%diagnosis%n_single_variable_constraints /= 1) error stop 'test_diagnostics FAILED: one single-variable constraint'
    if (r%diagnosis%n_constant_constraints /= 0) error stop 'test_diagnostics FAILED: no constant constraint'
    if (.not. r%diagnosis%variable_gradient_ratio > 1.0e6_wp) error stop 'test_diagnostics FAILED: poorly scaled variables'
    if (abs(r%diagnosis%constraint_gradient_ratio - 2.0e7_wp) > 1.0_wp) then
        error stop 'test_diagnostics FAILED: the ratio of the constraint gradients'
    end if
    end subroutine test_starting_point

    subroutine test_weakly_active()
    !! a constraint that is active at the solution with a zero multiplier
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat
    which = p_weak
    call solve(1, 0, options, r, istat)
    print '(A)', r%diagnosis%report
    if (istat /= sqpopt_success) error stop 'test_diagnostics FAILED: the weakly active problem should converge'
    if (r%diagnosis%n_active_constraints /= 1 .or. r%diagnosis%n_weakly_active /= 1) then
        error stop 'test_diagnostics FAILED: one weakly active constraint'
    end if
    if (r%diagnosis%n_dependent /= 0 .or. r%diagnosis%n_wrong_sign /= 0) then
        error stop 'test_diagnostics FAILED: no dependent gradient or wrong sign'
    end if
    end subroutine test_weakly_active

    subroutine test_probe()
    !! level 3 finds a starting point at the edge of the objective's domain (and not at level 2)
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat
    which = p_edge
    call solve(3, 0, options, r, istat)
    print '(A)', r%diagnosis%problem_report
    if (.not. r%diagnosis%probe_not_finite) error stop 'test_diagnostics FAILED: the probe should find the edge'
    call solve(2, 0, options, r, istat)
    if (r%diagnosis%probe_not_finite) error stop 'test_diagnostics FAILED: no probe at level 2'
    which = p_disk
    call solve(3, 0, options, r, istat)
    if (r%diagnosis%probe_not_finite) error stop 'test_diagnostics FAILED: the disk problem is defined everywhere'
    end subroutine test_probe

    subroutine test_function_error()
    !! a starting point where the objective can't be evaluated: the diagnosis says so, at every level
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat, level
    which = p_outside
    do level = 1, 3
        call solve(level, 0, options, r, istat)
        if (level == 1) print '(A)', r%diagnosis%report
        if (istat /= sqpopt_function_error) error stop 'test_diagnostics FAILED: the objective is not defined at x0'
        if (index(r%diagnosis%report, 'is not finite at the current point: the objective') == 0) then
            error stop 'test_diagnostics FAILED: no verdict for the function error'
        end if
    end do
    end subroutine test_function_error

    subroutine test_invalid()
    !! an invalid level, and a history unit that isn't open
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat
    which = p_disk
    call solve(4, 0, options, r, istat)
    print '(A,I0,2A)', 'level 4: istat = ', istat, '  ', r%message
    if (istat /= sqpopt_invalid_input) error stop 'test_diagnostics FAILED: diagnostic_level = 4 must be invalid'
    call solve(-1, 0, options, r, istat)
    if (istat /= sqpopt_invalid_input) error stop 'test_diagnostics FAILED: diagnostic_level = -1 must be invalid'
    options%diagnostics_unit = 987
    call solve(2, 0, options, r, istat)
    print '(A,I0,2A)', 'unit not open: istat = ', istat, '  ', r%message
    if (istat /= sqpopt_invalid_input) error stop 'test_diagnostics FAILED: a closed diagnostics_unit must be invalid'
    call solve(1, 0, options, r, istat)   ! (the unit is only used from level 2)
    if (istat /= sqpopt_success) error stop 'test_diagnostics FAILED: diagnostics_unit is not used at level 1'
    end subroutine test_invalid

    subroutine fc(x, f, c, status, data)
    !! the objective and the constraints of problem `which`
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (none)
    select case (which)
    case (p_disk)
        f = 100.0_wp*(x(2) - x(1)**2)**2 + (1.0_wp - x(1))**2
        c(1) = x(1)**2 + x(2)**2
    case (p_wrong)
        f = 100.0_wp*(x(2) - x(1)**2)**2 + (1.0_wp - x(1))**2
        c(1) = x(1)**2 + x(2)**2
        c(2) = x(1) + x(2)**2
    case (p_infeas)
        f = x(1)**2 + x(2)**2
        c(1) = x(1) + x(2)
        c(2) = x(1) + x(2)
    case (p_structure)
        f = (x(1) - 1.0_wp)**2 + (x(2) - 1.0_wp)**2 + (1.0e-7_wp*x(3) - 1.0_wp)**2
        c(1) = x(1) + x(2)
        c(2) = 2.0_wp*x(1) + 2.0_wp*x(2)
        c(3) = 1.0e-7_wp*x(3)
    case (p_weak)
        f = (x(1) - 1.0_wp)**2 + (x(2) - 2.0_wp)**2
        c(1) = x(1) + x(2)
    case (p_edge, p_outside)
        if (x(1) < 0.0_wp) then
            status = 1
            f = 0.0_wp
        else
            f = (sqrt(x(1)) - 1.0_wp)**2 + x(1)
        end if
    end select
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient and the constraint Jacobian of problem `which`
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    select case (which)
    case (p_disk)
        g = [-400.0_wp*x(1)*(x(2) - x(1)**2) - 2.0_wp*(1.0_wp - x(1)), 200.0_wp*(x(2) - x(1)**2)]
        jac_val = 2.0_wp*x
    case (p_wrong)
        g = [-400.0_wp*x(1)*(x(2) - x(1)**2) - 2.0_wp*(1.0_wp - x(1)), 200.0_wp*(x(2) - x(1)**2)]
        ! (the second row is wrong: its last entry should be `2*x(2)`)
        jac_val = [2.0_wp*x(1), 2.0_wp*x(2), 1.0_wp, 0.5_wp*x(2)]
    case (p_infeas)
        g = 2.0_wp*x
        jac_val = 1.0_wp
    case (p_structure)
        g = [2.0_wp*(x(1) - 1.0_wp), 2.0_wp*(x(2) - 1.0_wp), 2.0e-7_wp*(1.0e-7_wp*x(3) - 1.0_wp)]
        jac_val = [1.0_wp, 1.0_wp, 2.0_wp, 2.0_wp, 1.0e-7_wp]
    case (p_weak)
        g = 2.0_wp*(x - [1.0_wp, 2.0_wp])
        jac_val = 1.0_wp
    case (p_edge, p_outside)
        g(1) = 1.0_wp   ! (the derivative of `sqrt` is infinite at 0: only the other term's is given there)
        if (x(1) > 0.0_wp) g(1) = (sqrt(x(1)) - 1.0_wp)/sqrt(x(1)) + 1.0_wp
    end select
    end subroutine gjac

end program test_diagnostics
