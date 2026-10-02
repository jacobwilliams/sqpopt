program test_overdetermined

    !! Problems with more equality constraints than variables, which the
    !! solver handles by Gauss-Newton restoration steps (every QP subproblem
    !! is inconsistent; see `restoration_step`):
    !!
    !! * a consistent one whose solution is far away, Brown's badly scaled
    !!   problem as equations (`x1 = 1e6`, `x2 = 2e-6`, `x1*x2 = 2`, from
    !!   `(1, 1)`): the cap on the step's length must grow, or the solver
    !!   moves by `max_step` per iteration and never arrives.
    !! * an inconsistent one, Bard's problem as equations (15 equations
    !!   `y_i = x1 + u_i/(v_i*x2 + w_i*x3)` in 3 variables, from `(1, 1, 1)`):
    !!   the solver must stop as infeasible at the least-squares solution (a
    !!   sum of squares of 8.21487e-3, without scaling), not at its iteration
    !!   limit. (`test_nlls` solves the same problem with the least-squares
    !!   interface, which is the better way to pose it.)

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_types_module,   only: sqpopt_results_type, sqpopt_success, sqpopt_infeasible
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: p_brown   = 1 !! Brown's badly scaled problem
    integer, parameter :: p_bard    = 2 !! Bard's problem
    integer :: which !! the problem being solved
    real(wp), dimension(15), parameter :: bard_y = [0.14_wp, 0.18_wp, 0.22_wp, 0.25_wp, 0.29_wp, 0.32_wp, 0.35_wp, 0.39_wp, &
                                                    0.37_wp, 0.58_wp, 0.73_wp, 0.96_wp, 1.34_wp, 2.10_wp, 4.39_wp]
        !! the data of Bard's problem

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

    call solve(p_bard, .false., r, istat)
    if (istat /= sqpopt_infeasible) error stop 'test_overdetermined FAILED: Bard''s problem should end as infeasible'
    if (abs(sum(r%c**2) - 8.21487e-3_wp) > 1.0e-7_wp) then
        error stop 'test_overdetermined FAILED: Bard''s problem, not at the least-squares solution'
    end if
    if (r%iterations > 50) error stop 'test_overdetermined FAILED: Bard''s problem took too many iterations'

    ! (with the automatic scaling, which weights the constraints, so the least-squares solution may be another one)
    call solve(p_bard, .true., r, istat)
    if (istat /= sqpopt_infeasible) error stop 'test_overdetermined FAILED: Bard''s problem (scaled) should end as infeasible'
    if (r%iterations > 50) error stop 'test_overdetermined FAILED: Bard''s problem (scaled) took too many iterations'

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
    integer :: i, j, m, n
    real(wp), dimension(:), allocatable :: zero, lb, ub, x0
    integer,  dimension(:), allocatable :: irow, icol

    which = p
    n = merge(2, 3, p == p_brown)
    m = merge(3, 15, p == p_brown)
    allocate(zero(m), irow(n*m), icol(n*m), lb(n), ub(n), x0(n))
    zero = 0.0_wp
    lb = -inf
    ub = inf
    x0 = 1.0_wp
    ! (a dense Jacobian, by rows)
    do i = 1, m
        do j = 1, n
            irow(n*(i-1)+j) = i
            icol(n*(i-1)+j) = j
        end do
    end do
    call problem%set_problem_size(n=n, m=m)
    call problem%set_bounds(lb, ub, zero, zero)
    call problem%set_jacobian_sparsity(n*m, irow, icol)
    call problem%set_functions(fc=fc, gjac=gjac)
    options%scaling  = scaling
    options%max_iter = 250
    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_results(r)
    print '(A,I0,A,L1,A,I0,A,I0,A,I0,A,ES12.5,A,3ES16.8)', 'problem ', p, ', scaling ', scaling, ': istat = ', istat, &
        ', iterations = ', r%iterations, ', fc = ', r%n_eval_fc, ', sum of squares = ', sum(r%c**2), ', x = ', r%x
    end subroutine solve

    subroutine fc(x, f, c, status, data)
    !! the (constant) objective and the equations
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! the equations' residuals at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (none)
    integer :: i
    f = 0.0_wp
    if (which == p_brown) then
        c = [x(1) - 1.0e6_wp, x(2) - 2.0e-6_wp, x(1)*x(2) - 2.0_wp]
    else
        do i = 1, 15
            c(i) = bard_y(i) - (x(1) + i/((16 - i)*x(2) + min(i, 16 - i)*x(3)))
        end do
    end if
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient (zero) and the equations' Jacobian (by rows)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` `dimension(n*m)`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    integer :: i
    g = 0.0_wp
    if (which == p_brown) then
        jac_val = [1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, x(2), x(1)]
    else
        do i = 1, 15
            jac_val(3*i-2) = -1.0_wp
            jac_val(3*i-1) = i*(16 - i)/((16 - i)*x(2) + min(i, 16 - i)*x(3))**2
            jac_val(3*i)   = i*min(i, 16 - i)/((16 - i)*x(2) + min(i, 16 - i)*x(3))**2
        end do
    end if
    end subroutine gjac

end program test_overdetermined
