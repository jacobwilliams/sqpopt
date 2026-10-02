program test_nlls

    !! The least-squares interface ([[sqpopt_nlls_type]]):
    !!
    !! * Bard's problem (15 residuals `y_i - (x1 + u_i/(v_i*x2 + w_i*x3))`,
    !!   from `(1, 1, 1)`): must end with `sqpopt_success` at the sum of
    !!   squares 8.21487e-3, with the default L-BFGS Hessian and with the
    !!   Gauss-Newton Hessian. (Given as equality constraints to the solver
    !!   itself, such a problem ends as infeasible at best: see
    !!   `test_overdetermined`.)
    !! * Brown's badly scaled problem (a zero-residual one, in 2 variables).
    !! * a fit with a bound and a constraint: a line `x1 + x2*t` through
    !!   four points, with `x2 <= 0.5` and `x1 + x2 = 1`, and user data.
    !! * a second `solve` gives the same result as the first, and invalid
    !!   input is reported.

    use sqpopt_nlls_module,    only: sqpopt_nlls_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_hessian_module, only: sqpopt_hessian_exact
    use sqpopt_types_module,   only: sqpopt_results_type, sqpopt_success, sqpopt_invalid_input
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: p_brown   = 1 !! Brown's badly scaled problem
    integer, parameter :: p_bard    = 2 !! Bard's problem
    integer, parameter :: p_line    = 3 !! the constrained line fit

    type :: line_data
        !! the points of the line fit (its user data)
        real(wp), dimension(4) :: t = [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp] !! abscissas
        real(wp), dimension(4) :: y = [0.1_wp, 0.9_wp, 2.1_wp, 2.9_wp] !! ordinates
        integer :: n_calls = 0 !! number of residual evaluations
    end type line_data

    integer :: which !! the problem being solved
    real(wp), dimension(15), parameter :: bard_y = [0.14_wp, 0.18_wp, 0.22_wp, 0.25_wp, 0.29_wp, 0.32_wp, 0.35_wp, 0.39_wp, &
                                                    0.37_wp, 0.58_wp, 0.73_wp, 0.96_wp, 1.34_wp, 2.10_wp, 4.39_wp]
        !! the data of Bard's problem

    real(wp), dimension(:), allocatable :: x
    real(wp) :: ss
    integer :: istat, iter, n_fc

    write(*,*) '----------------------------'
    write(*,*) 'test_nlls'
    write(*,*) '----------------------------'

    call solve(p_bard, .false., x, ss, istat, iter, n_fc)
    if (istat /= sqpopt_success) error stop 'test_nlls FAILED: Bard''s problem should converge'
    if (abs(ss - 8.21487e-3_wp) > 1.0e-7_wp) error stop 'test_nlls FAILED: Bard''s problem, wrong sum of squares'
    if (maxval(abs(x - [0.0824106_wp, 1.13304_wp, 2.34370_wp])) > 1.0e-4_wp) then
        error stop 'test_nlls FAILED: Bard''s problem, wrong solution'
    end if

    call solve(p_bard, .true., x, ss, istat, iter, n_fc)
    if (istat /= sqpopt_success) error stop 'test_nlls FAILED: Bard''s problem (Gauss-Newton) should converge'
    if (abs(ss - 8.21487e-3_wp) > 1.0e-7_wp) error stop 'test_nlls FAILED: Bard''s problem (Gauss-Newton), wrong sum of squares'
    if (iter > 15) error stop 'test_nlls FAILED: Bard''s problem (Gauss-Newton) took too many iterations'

    call solve(p_brown, .false., x, ss, istat, iter, n_fc)
    if (istat /= sqpopt_success) error stop 'test_nlls FAILED: Brown''s problem should converge'
    if (ss > 1.0e-10_wp) error stop 'test_nlls FAILED: Brown''s problem, the residuals should be zero'

    call solve(p_brown, .true., x, ss, istat, iter, n_fc)
    if (istat /= sqpopt_success) error stop 'test_nlls FAILED: Brown''s problem (Gauss-Newton) should converge'
    if (ss > 1.0e-10_wp) error stop 'test_nlls FAILED: Brown''s problem (Gauss-Newton), the residuals should be zero'

    call test_line()
    call test_invalid()

    print '(A)', 'test_nlls PASSED'

    contains

    subroutine solve(p, gauss_newton, x, ss, istat, iter, n_fc)
    !! solve problem `p` (without constraints), twice, and check that the two solves agree
    integer,                intent(in)  :: p            !! the problem
    logical,                intent(in)  :: gauss_newton !! whether to use the Gauss-Newton Hessian
    real(wp), dimension(:), allocatable, intent(out) :: x !! the solution
    real(wp),               intent(out) :: ss           !! the sum of squares there
    integer,                intent(out) :: istat        !! the status
    integer,                intent(out) :: iter         !! the number of iterations
    integer,                intent(out) :: n_fc         !! the number of residual evaluations
    type(sqpopt_nlls_type)    :: nlls
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: i, j, l, n
    integer,  dimension(:), allocatable :: irow, icol
    real(wp), dimension(:), allocatable :: res, x0, x_again
    real(wp) :: ss_again

    which = p
    n = merge(2, 3, p == p_brown)
    l = merge(3, 15, p == p_brown)
    allocate(irow(n*l), icol(n*l), res(l), x(n), x0(n), x_again(n))
    ! (dense Jacobians, by rows)
    do i = 1, l
        do j = 1, n
            irow(n*(i-1)+j) = i
            icol(n*(i-1)+j) = j
        end do
    end do
    if (gauss_newton) options%hessian_mode = sqpopt_hessian_exact
    options%max_iter = 500
    call nlls%initialize(n=n, n_residuals=l, residuals=residuals, jacobian=jacobian, rjac_irow=irow, rjac_icol=icol, &
                         options=options)
    x0 = 1.0_wp
    call nlls%solve(x0, istat)
    call nlls%get_solution(x, residuals=res, sum_of_squares=ss)
    call nlls%get_results(r)
    iter = r%iterations
    n_fc = r%n_eval_fc
    print '(A,I0,A,L1,A,I0,A,I0,A,I0,A,ES12.5,A,3ES16.8)', 'problem ', p, ', Gauss-Newton ', gauss_newton, &
        ': istat = ', istat, ', iterations = ', iter, ', evaluations = ', n_fc, ', sum of squares = ', ss, ', x = ', x
    if (abs(sum(res**2) - ss) > 1.0e-12_wp*max(1.0_wp, ss)) error stop 'test_nlls FAILED: residuals and sum of squares disagree'

    call nlls%solve(x0, istat)
    call nlls%get_solution(x_again, sum_of_squares=ss_again)
    call nlls%get_results(r)
    if (any(x_again /= x) .or. ss_again /= ss .or. r%iterations /= iter .or. r%n_eval_fc /= n_fc) then
        error stop 'test_nlls FAILED: a second solve differs from the first'
    end if
    call nlls%destroy()
    end subroutine solve

    subroutine test_line()
    !! the line fit with a bound, a constraint, and user data
    type(sqpopt_nlls_type) :: nlls
    type(line_data) :: points
    real(wp), dimension(2) :: x
    real(wp), dimension(1) :: lambda
    real(wp), dimension(4) :: res
    real(wp) :: ss
    integer :: istat, i
    real(wp), parameter :: inf = 1.0e20_wp

    which = p_line
    call nlls%initialize(n=2, n_residuals=4, residuals=residuals, jacobian=jacobian, &
                         rjac_irow=[(i, i, i = 1, 4)], rjac_icol=[(1, 2, i = 1, 4)], &
                         x_lb=[-inf, -inf], x_ub=[inf, 0.5_wp], c_lb=[1.0_wp], c_ub=[1.0_wp], &
                         cjac_irow=[1, 1], cjac_icol=[1, 2], data=points)
    call nlls%solve([0.0_wp, 0.0_wp], istat)
    call nlls%get_solution(x, residuals=res, sum_of_squares=ss, lambda=lambda)
    print '(A,I0,A,2ES16.8,A,ES12.5,A,ES12.5,A,I0)', 'line fit: istat = ', istat, ', x = ', x, ', sum of squares = ', ss, &
        ', multiplier = ', lambda(1), ', evaluations = ', points%n_calls
    if (istat /= sqpopt_success) error stop 'test_nlls FAILED: the line fit should converge'
    ! (on `x1 + x2 = 1` the sum of squares is least at `x2 = 5.8/6`, so the bound `x2 <= 0.5` is active)
    if (abs(x(1) - 0.5_wp) > 1.0e-6_wp .or. abs(x(2) - 0.5_wp) > 1.0e-6_wp) error stop 'test_nlls FAILED: line fit, wrong solution'
    if (abs(ss - sum((x(1) + x(2)*points%t - points%y)**2)) > 1.0e-8_wp) error stop 'test_nlls FAILED: line fit, wrong sum of squares'
    if (maxval(abs(res - (x(1) + x(2)*points%t - points%y))) > 1.0e-8_wp) error stop 'test_nlls FAILED: line fit, wrong residuals'
    if (points%n_calls == 0) error stop 'test_nlls FAILED: the user data was not passed'
    call nlls%destroy()
    end subroutine test_line

    subroutine test_invalid()
    !! invalid input must be reported by `solve`, not crash
    type(sqpopt_nlls_type) :: nlls
    real(wp), dimension(2) :: x
    integer :: istat

    which = p_brown
    ! (not initialized)
    call nlls%solve([1.0_wp, 1.0_wp], istat)
    if (istat /= sqpopt_invalid_input) error stop 'test_nlls FAILED: a solve without initialize should be invalid input'
    call nlls%get_solution(x)
    ! (a column index out of range)
    call nlls%initialize(n=2, n_residuals=3, residuals=residuals, jacobian=jacobian, &
                         rjac_irow=[1, 1, 2, 2, 3, 3], rjac_icol=[1, 2, 1, 2, 1, 7])
    call nlls%solve([1.0_wp, 1.0_wp], istat)
    if (istat /= sqpopt_invalid_input) error stop 'test_nlls FAILED: a bad Jacobian pattern should be invalid input'
    ! (a starting point of the wrong size)
    call nlls%initialize(n=2, n_residuals=3, residuals=residuals, jacobian=jacobian, &
                         rjac_irow=[1, 1, 2, 2, 3, 3], rjac_icol=[1, 2, 1, 2, 1, 2])
    call nlls%solve([1.0_wp, 1.0_wp, 1.0_wp], istat)
    if (istat /= sqpopt_invalid_input) error stop 'test_nlls FAILED: a starting point of the wrong size should be invalid input'
    call nlls%destroy()
    end subroutine test_invalid

    subroutine residuals(x, r, c, status, data)
    !! the residuals (and the constraint of the line fit)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: r      !! residuals at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraints at `x` (one for the line fit, none otherwise)
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data (a `line_data` for the line fit)
    integer :: i
    select case (which)
    case (p_brown)
        r = [x(1) - 1.0e6_wp, x(2) - 2.0e-6_wp, x(1)*x(2) - 2.0_wp]
    case (p_bard)
        do i = 1, 15
            r(i) = bard_y(i) - (x(1) + i/((16 - i)*x(2) + min(i, 16 - i)*x(3)))
        end do
    case default
        r = 0.0_wp
        if (present(data)) then
            select type (data)
            type is (line_data)
                r = x(1) + x(2)*data%t - data%y
                data%n_calls = data%n_calls + 1
            end select
        end if
        c(1) = x(1) + x(2)
    end select
    end subroutine residuals

    subroutine jacobian(x, rjac_val, cjac_val, accuracy, status, data)
    !! the Jacobians of the residuals (by rows) and of the constraint
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: rjac_val !! nonzeros of the residuals' Jacobian at `x`
    real(wp), dimension(:), intent(out)   :: cjac_val !! nonzeros of the constraint's Jacobian (line fit only)
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data (a `line_data` for the line fit)
    integer :: i
    select case (which)
    case (p_brown)
        rjac_val = [1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, x(2), x(1)]
    case (p_bard)
        do i = 1, 15
            rjac_val(3*i-2) = -1.0_wp
            rjac_val(3*i-1) = i*(16 - i)/((16 - i)*x(2) + min(i, 16 - i)*x(3))**2
            rjac_val(3*i)   = i*min(i, 16 - i)/((16 - i)*x(2) + min(i, 16 - i)*x(3))**2
        end do
    case default
        rjac_val = 0.0_wp
        if (present(data)) then
            select type (data)
            type is (line_data)
                do i = 1, 4
                    rjac_val(2*i-1) = 1.0_wp
                    rjac_val(2*i)   = data%t(i)
                end do
            end select
        end if
        cjac_val = 1.0_wp
    end select
    end subroutine jacobian

end program test_nlls
