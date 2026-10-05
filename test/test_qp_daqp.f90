program test_qp_daqp

    !! Tests of the DAQP QP solver in the SQP iterations
    !! (`options%qp_solver_mode = sqpopt_qp_daqp`, see
    !! [[sqpopt_qp_daqp_module]]), and of its fallback to the dense QP solver:
    !!
    !! * a convex problem: DAQP solves every QP (no fallback), and the solve
    !!   ends where the dense QP solver's does, with the same multiplier;
    !! * inconsistent linearized constraints: DAQP reports them, the dense
    !!   solver's elastic mode takes over (counted in
    !!   `results%n_daqp_fallbacks`), and the solve ends as
    !!   `sqpopt_infeasible`, as with the dense solver;
    !! * a nonconvex problem with the exact Hessian (negative definite):
    !!   DAQP can't factor it, so those QPs fall back to the dense solver
    !!   (which finds the negative curvature, and the solver shifts the
    !!   Hessian), until the shift makes it positive definite, and DAQP
    !!   solves them; the solve is the dense solver's, step for step;
    !! * the detailed log (`print_level = 3`), which reports the fallbacks,
    !!   doesn't change the results.

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_exact
    use sqpopt_qp_solver_module, only: sqpopt_qp_dense, sqpopt_qp_daqp
    use sqpopt_types_module,     only: sqpopt_success, sqpopt_infeasible, sqpopt_results_type
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: tol = 1.0e-6_wp

    write(*,*) '----------------------------'
    write(*,*) 'test_qp_daqp'
    write(*,*) '----------------------------'

    call test_convex()
    call test_infeasible()
    call test_nonconvex()

    print '(A)', 'test_qp_daqp PASSED'

    contains

    subroutine check(ok, msg)
    !! stop with a failure message unless `ok`
    logical,          intent(in) :: ok  !! the condition
    character(len=*), intent(in) :: msg !! what failed
    if (.not. ok) error stop 'test_qp_daqp FAILED: '//msg
    end subroutine check

    subroutine run(problem, options, x0, r)
    !! solve the problem, and return its results
    type(sqpopt_problem_type), intent(in)  :: problem !! the problem
    type(sqpopt_options_type), intent(in)  :: options !! the options
    real(wp), dimension(:),    intent(in)  :: x0      !! the starting point
    type(sqpopt_results_type), intent(out) :: r       !! the results
    type(sqpopt_type) :: solver
    integer :: istat
    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_results(r)
    call check(istat == r%istat, 'istat and results%istat differ')
    end subroutine run

    subroutine test_convex()
    !! min (x1-1)^2 + (x2-2)^2  s.t.  x1 + x2 <= 2,  -5 <= x <= 5: the solution is (0.5, 1.5), with
    !! the multiplier -1 (at the constraint's upper bound)
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r_dense, r_daqp
    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-5.0_wp,-5.0_wp], x_ub=[5.0_wp,5.0_wp], c_lb=[-1.0e20_wp], c_ub=[2.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_convex, gjac=gjac_convex)
    options%qp_solver_mode = sqpopt_qp_dense
    call run(problem, options, [-3.0_wp, 4.0_wp], r_dense)
    options%qp_solver_mode = sqpopt_qp_daqp
    call run(problem, options, [-3.0_wp, 4.0_wp], r_daqp)
    print '(A,2F10.6,A,F10.6,A,I0,A,I0,A,I0)', 'convex:     x=', r_daqp%x, '  lambda=', r_daqp%lambda, &
        '  istat=', r_daqp%istat, '  QPs=', r_daqp%n_qp_solves, '  fallbacks=', r_daqp%n_daqp_fallbacks
    call check(r_daqp%istat == sqpopt_success, 'convex: not solved')
    call check(maxval(abs(r_daqp%x - [0.5_wp, 1.5_wp])) <= tol, 'convex: wrong solution')
    call check(abs(r_daqp%lambda(1) + 1.0_wp) <= tol, 'convex: wrong multiplier')
    call check(r_daqp%n_qp_solves > 0, 'convex: no QP solved')
    call check(r_daqp%n_daqp_fallbacks == 0, 'convex: DAQP fell back to the dense QP solver')
    call check(r_dense%n_daqp_fallbacks == 0, 'convex: fallbacks counted without DAQP')
    call check(maxval(abs(r_daqp%x - r_dense%x)) <= tol, 'convex: not the dense QP solver''s solution')
    end subroutine test_convex

    subroutine test_infeasible()
    !! x1 in [0,1] and x1 in [2,3] (the problem of `test_infeasible`): the linearized constraints are
    !! inconsistent, so the dense solver's elastic mode takes over, and the least-infeasible point is
    !! x1 = 1.5
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    call problem%set_problem_size(n=2, m=2)
    call problem%set_bounds(x_lb=[-10.0_wp,-10.0_wp], x_ub=[10.0_wp,10.0_wp], &
                            c_lb=[0.0_wp,2.0_wp], c_ub=[1.0_wp,3.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,2], icol=[1,1])
    call problem%set_functions(fc=fc_infeasible, gjac=gjac_infeasible)
    options%qp_solver_mode = sqpopt_qp_daqp
    call run(problem, options, [0.5_wp, 0.0_wp], r)
    print '(A,2F10.6,A,I0,A,I0,A,I0)', 'infeasible: x=', r%x, '  istat=', r%istat, &
        '  QPs=', r%n_qp_solves, '  fallbacks=', r%n_daqp_fallbacks
    call check(r%istat == sqpopt_infeasible, 'infeasible: not reported as sqpopt_infeasible')
    call check(abs(r%x(1) - 1.5_wp) <= 1.0e-4_wp, 'infeasible: not at the least-infeasible point')
    call check(r%n_daqp_fallbacks > 0, 'infeasible: no fallback counted')
    end subroutine test_infeasible

    subroutine test_nonconvex()
    !! min -(x1^2 + x2^2) + 0.1 x1  s.t.  x1 + x2 = 0.5,  -1 <= x <= 1, with the exact Hessian (-2I):
    !! the local minima are the ends of the segment, (1,-0.5) and (-0.5,1)
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r_dense, r_daqp, r_log
    integer :: unit, ios
    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-1.0_wp,-1.0_wp], x_ub=[1.0_wp,1.0_wp], c_lb=[0.5_wp], c_ub=[0.5_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_hessian_sparsity(nnz=2, irow=[1,2], icol=[1,2])
    call problem%set_functions(fc=fc_nonconvex, gjac=gjac_nonconvex, hess=hess_nonconvex)
    options%hessian_mode   = sqpopt_hessian_exact
    options%qp_solver_mode = sqpopt_qp_dense
    call run(problem, options, [0.3_wp, 0.2_wp], r_dense)
    options%qp_solver_mode = sqpopt_qp_daqp
    call run(problem, options, [0.3_wp, 0.2_wp], r_daqp)
    print '(A,2F10.6,A,I0,A,I0,A,I0)', 'nonconvex:  x=', r_daqp%x, '  istat=', r_daqp%istat, &
        '  QPs=', r_daqp%n_qp_solves, '  fallbacks=', r_daqp%n_daqp_fallbacks
    call check(r_daqp%istat == sqpopt_success, 'nonconvex: not solved')
    call check(min(maxval(abs(r_daqp%x - [1.0_wp, -0.5_wp])), maxval(abs(r_daqp%x - [-0.5_wp, 1.0_wp]))) <= tol, &
               'nonconvex: not at a local minimum')
    call check(r_daqp%n_daqp_fallbacks > 0, 'nonconvex: no fallback counted')
    call check(all(r_daqp%x == r_dense%x) .and. r_daqp%n_eval_fc == r_dense%n_eval_fc .and. &
               r_daqp%iterations == r_dense%iterations, 'nonconvex: not the dense QP solver''s solve')

    ! the detailed log reports the fallbacks, and doesn't change the results
    open(newunit=unit, status='scratch', action='readwrite', iostat=ios)
    call check(ios == 0, 'could not open a scratch file')
    options%print_level = 3
    options%output_unit = unit
    call run(problem, options, [0.3_wp, 0.2_wp], r_log)
    call check(all(r_log%x == r_daqp%x) .and. r_log%n_eval_fc == r_daqp%n_eval_fc .and. &
               r_log%n_daqp_fallbacks == r_daqp%n_daqp_fallbacks, 'the detailed log changed the results')
    call check(log_has(unit, 'DAQP: not solved (nonconvex)'), 'the log doesn''t report the fallbacks')
    call check(log_has(unit, 'DAQP fallbacks'), 'the summary doesn''t report the fallbacks')
    close(unit)
    end subroutine test_nonconvex

    logical function log_has(unit, text)
    !! whether a line of the file open on `unit` contains `text`
    integer,          intent(in) :: unit !! the unit
    character(len=*), intent(in) :: text !! the text to look for
    character(len=512) :: line
    integer :: ios
    log_has = .false.
    rewind(unit)
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line, text) > 0) then
            log_has = .true.
            exit
        end if
    end do
    end function log_has

    subroutine fc_convex(x, f, c, status, data)
    !! objective and constraint of the convex problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = (x(1)-1.0_wp)**2 + (x(2)-2.0_wp)**2
    c(1) = x(1) + x(2)
    end subroutine fc_convex

    subroutine gjac_convex(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian of the convex problem
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac      !! nonzero values of the constraint Jacobian at `x`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are exact)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    g = [2.0_wp*(x(1)-1.0_wp), 2.0_wp*(x(2)-2.0_wp)]
    jac = [1.0_wp, 1.0_wp]
    end subroutine gjac_convex

    subroutine fc_infeasible(x, f, c, status, data)
    !! objective and constraints of the infeasible problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = x(1)**2 + x(2)**2
    c = [x(1), x(1)]
    end subroutine fc_infeasible

    subroutine gjac_infeasible(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian of the infeasible problem
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac      !! nonzero values of the constraint Jacobian at `x`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are exact)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    g = 2.0_wp*x
    jac = [1.0_wp, 1.0_wp]
    end subroutine gjac_infeasible

    subroutine fc_nonconvex(x, f, c, status, data)
    !! objective and constraint of the nonconvex problem
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f = -(x(1)**2 + x(2)**2) + 0.1_wp*x(1)
    c(1) = x(1) + x(2)
    end subroutine fc_nonconvex

    subroutine gjac_nonconvex(x, g, jac, accuracy, status, data)
    !! objective gradient and constraint Jacobian of the nonconvex problem
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac      !! nonzero values of the constraint Jacobian at `x`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are exact)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    g = [-2.0_wp*x(1) + 0.1_wp, -2.0_wp*x(2)]
    jac = [1.0_wp, 1.0_wp]
    end subroutine gjac_nonconvex

    subroutine hess_nonconvex(x, lambda, hess_val, status, data)
    !! the Hessian of the Lagrangian of the nonconvex problem (its constraint is linear): `-2I`
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    hess_val = [-2.0_wp, -2.0_wp]
    end subroutine hess_nonconvex

end program test_qp_daqp
