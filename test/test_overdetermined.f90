program test_overdetermined

    !! Problems with more equality constraints than variables, which the
    !! solver handles by Gauss-Newton restoration steps (every QP subproblem
    !! is inconsistent; see `restoration_step`):
    !!
    !! * a consistent one whose solution is far away, Brown's badly scaled
    !!   problem as equations (`x1 = 1e6`, `x2 = 2e-6`, `x1*x2 = 2`, from
    !!   `(1, 1)`): the cap on the step's length must grow, or the solver
    !!   moves by `max_step` per iteration and never arrives.
    !! * an inconsistent one with a Jacobian that becomes nearly rank
    !!   deficient, Jennrich and Sampson's problem as equations
    !!   (`2 + 2i = exp(i*x1) + exp(i*x2)`, `i = 1..10`, from `(0.3, 0.4)`):
    !!   the Gauss-Newton direction is then useless, and the steps need
    !!   Levenberg-Marquardt damping. The solver must stop as infeasible at
    !!   the least-squares solution (a sum of squares of 124.362 at
    !!   `x1 = x2 = 0.2578`, without scaling), not at its iteration limit.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_types_module,   only: sqpopt_results_type, sqpopt_success, sqpopt_infeasible
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: p_brown   = 1 !! Brown's badly scaled problem
    integer, parameter :: p_jensamp = 2 !! Jennrich and Sampson's problem
    integer :: which !! the problem being solved

    type(sqpopt_results_type) :: r
    integer :: istat

    write(*,*) '----------------------------'
    write(*,*) 'test_overdetermined'
    write(*,*) '----------------------------'

    call solve(p_brown, .true., r, istat)
    if (istat /= sqpopt_success) error stop 'test_overdetermined FAILED: Brown''s problem should converge'
    if (abs(r%x(1) - 1.0e6_wp) > 1.0e-2_wp .or. abs(r%x(2) - 2.0e-6_wp) > 1.0e-10_wp) then
        error stop 'test_overdetermined FAILED: Brown''s problem, wrong solution'
    end if
    if (r%iterations > 60) error stop 'test_overdetermined FAILED: Brown''s problem took too many iterations'

    call solve(p_jensamp, .false., r, istat)
    if (istat /= sqpopt_infeasible) error stop 'test_overdetermined FAILED: Jennrich-Sampson should end as infeasible'
    if (abs(sum(r%c**2) - 124.362_wp) > 1.0e-2_wp) then
        error stop 'test_overdetermined FAILED: Jennrich-Sampson, not at the least-squares solution'
    end if
    if (r%iterations > 100) error stop 'test_overdetermined FAILED: Jennrich-Sampson took too many iterations'

    ! (with the automatic scaling the constraints are weighted, so the least-squares solution is another one)
    call solve(p_jensamp, .true., r, istat)
    if (istat /= sqpopt_infeasible) error stop 'test_overdetermined FAILED: Jennrich-Sampson (scaled) should end as infeasible'
    if (r%iterations > 100) error stop 'test_overdetermined FAILED: Jennrich-Sampson (scaled) took too many iterations'

    print '(A)', 'test_overdetermined PASSED'

    contains

    subroutine solve(p, scaling, r, istat)
    !! solve problem `p`
    integer,                   intent(in)  :: p       !! the problem
    logical,                   intent(in)  :: scaling !! whether the automatic scaling is on
    type(sqpopt_results_type), intent(out) :: r       !! the results
    integer,                   intent(out) :: istat   !! the status
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    real(wp), parameter :: inf = 1.0e20_wp
    integer :: i, m
    real(wp), dimension(:), allocatable :: zero
    integer,  dimension(:), allocatable :: irow, icol

    which = p
    m = merge(3, 10, p == p_brown)
    allocate(zero(m), irow(2*m), icol(2*m))
    zero = 0.0_wp
    do i = 1, m
        irow(2*i-1:2*i) = i
        icol(2*i-1) = 1
        icol(2*i)   = 2
    end do
    call problem%set_problem_size(n=2, m=m)
    call problem%set_bounds([-inf, -inf], [inf, inf], zero, zero)
    call problem%set_jacobian_sparsity(2*m, irow, icol)
    call problem%set_functions(fc=fc, gjac=gjac)
    options%scaling  = scaling
    options%max_iter = 250
    call solver%initialize(problem=problem, options=options)
    if (p == p_brown) then
        call solver%solve([1.0_wp, 1.0_wp], istat)
    else
        call solver%solve([0.3_wp, 0.4_wp], istat)
    end if
    call solver%get_results(r)
    print '(A,I0,A,L1,A,I0,A,I0,A,I0,A,2ES16.8,A,ES12.5)', 'problem ', p, ', scaling ', scaling, ': istat = ', istat, &
        ', iterations = ', r%iterations, ', fc = ', r%n_eval_fc, ', x = ', r%x, ', sum of squares = ', sum(r%c**2)
    end subroutine solve

    subroutine fc(x, f, c, status, data)
    !! the (constant) objective and the equations
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(2)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! the equations' residuals at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (none)
    integer :: i
    f = 0.0_wp
    if (which == p_brown) then
        c = [x(1) - 1.0e6_wp, x(2) - 2.0e-6_wp, x(1)*x(2) - 2.0_wp]
    else
        do i = 1, size(c)
            c(i) = 2.0_wp + 2.0_wp*i - (exp(i*x(1)) + exp(i*x(2)))
        end do
    end if
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient (zero) and the equations' Jacobian (by rows)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(2)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` `dimension(2*m)`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    integer :: i
    g = 0.0_wp
    if (which == p_brown) then
        jac_val = [1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, x(2), x(1)]
    else
        do i = 1, size(jac_val)/2
            jac_val(2*i-1) = -i*exp(i*x(1))
            jac_val(2*i)   = -i*exp(i*x(2))
        end do
    end if
    end subroutine gjac

end program test_overdetermined
