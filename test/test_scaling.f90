program test_scaling

    !! Two safeguards against a starting point with a very large gradient.
    !!
    !! The limit on the scale factors of the gradient-based scaling
    !! (`options%scaling_min_value`). The Zakharov function (see
    !! [[scalable_functions_module]]) with 200 variables has a gradient of
    !! about `4e14` at its starting point, and of order 1 near its minimizer.
    !! Scaled by `100/4e14`, its convergence test means nothing in the
    !! original units, and the solve stops far from the minimum (0), at an
    !! objective of about 40 (without the second safeguard, below). With the default limit of `1e-8` on the factor it
    !! must come close to the minimum (to within `f_tol`: the factor is still
    !! small, so the result is less accurate than without scaling, which must
    !! reach the minimum to roundoff).
    !!
    !! And the test on the objective's change in the acceptable-level stop
    !! (`options%acceptable_obj_change_tol`). Schwefel's function 1.2 with 200
    !! variables has a starting gradient of `2e6`, within that limit. Its
    !! scaled problem passes the acceptable test while the objective is still
    !! falling by a large fraction at every iteration: without the test on
    !! the change the solve stops there, at an objective of about 0.8, and
    !! with it, it must go on to an objective a hundred times smaller.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_types_module,   only: sqpopt_results_type, sqpopt_success, sqpopt_acceptable
    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use scalable_functions_module

    implicit none

    integer, parameter :: n = 200
    real(wp), parameter :: f_tol = 1.0e-3_wp !! tolerance on the minimum with the default limit
    type(scalable_function_type) :: fun
    type(sqpopt_results_type) :: r_limit, r_none, r_off, r_change, r_no_change
    integer :: id

    write(*,*) '----------------------------'
    write(*,*) 'test_scaling'
    write(*,*) '----------------------------'

    do id = 1, n_scalable_functions
        if (scalable_function_name(id) == 'zakharov') call scalable_function_setup(id, n, fun)
    end do
    if (fun%n /= n) error stop 'test_scaling FAILED: no zakharov function'

    call run(1.0e-8_wp, .true., 1.0e-3_wp, r_limit)   ! (the defaults)
    call run(0.0_wp, .true., 1.0e20_wp, r_none)        ! (without either safeguard)
    call run(0.0_wp, .false., 1.0e-3_wp, r_off)
    print '(A,ES10.2,A,I0,A,I0)', 'scaling_min_value = 1e-8: f = ', r_limit%f, ', istat = ', r_limit%istat, &
                                  ', iterations = ', r_limit%iterations
    print '(A,ES10.2,A,I0,A,I0)', 'scaling_min_value = 0:    f = ', r_none%f, ', istat = ', r_none%istat, &
                                  ', iterations = ', r_none%iterations
    print '(A,ES10.2,A,I0,A,I0)', 'scaling off:              f = ', r_off%f, ', istat = ', r_off%istat, &
                                  ', iterations = ', r_off%iterations

    if (r_limit%istat /= sqpopt_success .or. r_limit%f > f_tol) then
        error stop 'test_scaling FAILED: the minimum was not reached with the default limit'
    end if
    if (r_off%istat /= sqpopt_success .or. r_off%f > 1.0e-10_wp) then
        error stop 'test_scaling FAILED: the minimum was not reached without scaling'
    end if
    ! (what the limit is for: without it, the solve stops far from the minimum)
    if (.not. r_none%f > 1.0e3_wp*f_tol) then
        error stop 'test_scaling FAILED: the limit made no difference (is the test problem still badly scaled?)'
    end if

    ! the acceptable-level stop, on Schwefel's function 1.2:
    do id = 1, n_scalable_functions
        if (scalable_function_name(id) == 'schwefel12') call scalable_function_setup(id, n, fun)
    end do
    if (fun%name /= 'schwefel12') error stop 'test_scaling FAILED: no schwefel12 function'
    call run(1.0e-8_wp, .true., 1.0e-3_wp, r_change)   ! (the defaults)
    call run(1.0e-8_wp, .true., 1.0e20_wp, r_no_change)
    print '(A,ES10.2,A,I0,A,I0)', 'acceptable_obj_change_tol = 1e-3: f = ', r_change%f, ', istat = ', r_change%istat, &
                                  ', iterations = ', r_change%iterations
    print '(A,ES10.2,A,I0,A,I0)', 'acceptable_obj_change_tol = 1e20: f = ', r_no_change%f, ', istat = ', &
                                  r_no_change%istat, ', iterations = ', r_no_change%iterations
    if (r_no_change%istat /= sqpopt_acceptable) then
        error stop 'test_scaling FAILED: without the test on the change, the solve did not stop as acceptable'
    end if
    if (r_change%istat /= sqpopt_success .and. r_change%istat /= sqpopt_acceptable) then
        error stop 'test_scaling FAILED: the solve did not converge with the test on the change'
    end if
    if (.not. r_change%f < 1.0e-2_wp*r_no_change%f) then
        error stop 'test_scaling FAILED: the test on the change did not make the solve go on'
    end if

    print '(A)', 'test_scaling PASSED'

    contains

    subroutine run(min_value, scaling, obj_change, r)
    !! solve the problem of `fun`
    real(wp),                  intent(in)  :: min_value  !! `options%scaling_min_value`
    logical,                   intent(in)  :: scaling    !! `options%scaling`
    real(wp),                  intent(in)  :: obj_change !! `options%acceptable_obj_change_tol`
    type(sqpopt_results_type), intent(out) :: r         !! the results
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    integer, dimension(0) :: no_rows
    integer :: istat
    call problem%set_problem_size(n=n, m=0)
    call problem%set_bounds(x_lb=fun%x_lb, x_ub=fun%x_ub, c_lb=[real(wp) ::], c_ub=[real(wp) ::])
    call problem%set_jacobian_sparsity(nnz=0, irow=no_rows, icol=no_rows)
    call problem%set_functions(fc=fc, gjac=gjac)
    options%max_iter          = 1000
    options%scaling           = scaling
    options%scaling_min_value = min_value
    options%acceptable_obj_change_tol = obj_change
    call solver%initialize(problem=problem, options=options)
    call solver%solve(fun%x0, istat)
    call solver%get_results(r)
    end subroutine run

    subroutine fc(x, f, c, status, data)
    !! the objective (there are no constraints)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(0)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (none)
    f = fun%f(x)
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient (the Jacobian is empty)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` `dimension(0)`
    integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (none)
    call fun%g(x, g)
    end subroutine gjac

end program test_scaling
