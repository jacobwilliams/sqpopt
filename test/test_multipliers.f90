program test_multipliers

    !! Test of the least-squares multiplier estimate
    !! ([[multiplier_estimate]]), and of its use with the exact Hessian (see
    !! [[sqpopt_iterate]]):
    !!
    !! * the estimate for a small Jacobian, against the solution of its normal
    !!   equations: for all the rows and variables, for a subset of the rows,
    !!   and with a variable at a bound; by `LSQR`, and, in a build with MUMPS,
    !!   by the direct solver too;
    !! * a hanging chain of 200 links of length 0.01 between (0,0) and (1,0):
    !!   minimize the sum of the joints' heights, with the links' lengths as
    !!   equality constraints. The objective is linear, so the exact Hessian
    !!   of the Lagrangian is all multipliers, and the first QPs need large
    !!   shifts. Without the estimate, the multipliers of those QPs made the
    !!   next Hessians more indefinite still, until the solver stopped as
    !!   stalled far from the solution (objective -65.2). It must now converge
    !!   to the catenary (objective -91.12), with multipliers of the right
    !!   size; in a build with MUMPS, also with inertia control and the direct
    !!   QP method.

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_exact
    use sqpopt_least_squares_module, only: sqpopt_least_squares_type, multiplier_estimate
    use sqpopt_symmetric_solver_module, only: sqpopt_has_mumps
    use sqpopt_types_module,     only: sqpopt_success, sqpopt_results_type, sqpopt_sparse_matrix
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    write(*,*) '----------------------------'
    write(*,*) 'test_multipliers'
    write(*,*) '----------------------------'

    call test_estimate()
    call test_chain('exact Hessian', .false., .false.)
    if (sqpopt_has_mumps) then
        call test_chain('exact, inertia control', .true., .false.)
        call test_chain('exact, inertia control, direct', .true., .true.)
    end if

    print '(A)', 'test_multipliers PASSED'

    contains

    subroutine test_estimate()
    !! the estimate for `J = [1 2 0; 0 1 1]` and `g = (1, 2, 3)`
    type(sqpopt_sparse_matrix) :: jac
    type(sqpopt_least_squares_type) :: ls
    real(wp) :: lambda(2), g(3)
    logical :: ok
    integer :: pass

    jac%nrows = 2
    jac%ncols = 3
    jac%nnz   = 4
    jac%irow  = [1, 1, 2, 2]
    jac%icol  = [1, 2, 2, 3]
    jac%val   = [1.0_wp, 2.0_wp, 1.0_wp, 1.0_wp]
    g = [1.0_wp, 2.0_wp, 3.0_wp]

    ! (pass 1: LSQR, since the direct solver `ls` isn't started; pass 2, with MUMPS: the direct solver)
    do pass = 1, merge(2, 1, sqpopt_has_mumps)
        if (pass == 2) then
            call ls%initialize(3, 2, jac%irow, jac%icol, ok)
            if (.not. ok) error stop 'test_multipliers FAILED: the least-squares solver could not be started'
        end if

        ! both rows, every variable free: (J J^T) lambda = J g, with J J^T = [5 2; 2 2], J g = (5, 5)
        lambda = 0.0_wp
        call multiplier_estimate(jac, [.true., .true.], [.true., .true., .true.], g, lambda, ok, least_squares=ls)
        if (.not. ok .or. maxval(abs(lambda - [0.0_wp, 2.5_wp])) > 1.0e-6_wp) error stop 'test_multipliers FAILED: both rows'

        ! the first row only: lambda_1 = (J_1 g)/(J_1 J_1^T) = 5/5; the other multiplier is untouched
        lambda = [0.0_wp, 7.0_wp]
        call multiplier_estimate(jac, [.true., .false.], [.true., .true., .true.], g, lambda, ok, least_squares=ls)
        if (.not. ok .or. maxval(abs(lambda - [1.0_wp, 7.0_wp])) > 1.0e-6_wp) error stop 'test_multipliers FAILED: one row'

        ! the second variable at a bound: J = [1 0 0; 0 0 1] in the free ones, so lambda = (g_1, g_3)
        lambda = 0.0_wp
        call multiplier_estimate(jac, [.true., .true.], [.true., .false., .true.], g, lambda, ok, least_squares=ls)
        if (.not. ok .or. maxval(abs(lambda - [1.0_wp, 3.0_wp])) > 1.0e-6_wp) then
            error stop 'test_multipliers FAILED: a variable at a bound'
        end if
        ! (the direct solver must have been used in pass 2, and only then)
        if ((ls%kkt%solver%n_factor > 0) .neqv. (pass == 2)) error stop 'test_multipliers FAILED: the wrong solver was used'

        ! no rows: nothing to estimate
        if (pass == 1) then
            call multiplier_estimate(jac, [.false., .false.], [.true., .true., .true.], g, lambda, ok)
            if (ok) error stop 'test_multipliers FAILED: an estimate without rows'
        end if
    end do
    call ls%destroy()

    print '(A)', 'test_multipliers [the estimate] PASSED'

    end subroutine test_estimate

    subroutine test_chain(label, inertia, direct)
    !! solve the hanging chain with the exact Hessian (see the program's documentation)
    character(len=*), intent(in) :: label   !! the configuration, for the output
    logical,          intent(in) :: inertia !! `options%inertia_control`
    logical,          intent(in) :: direct  !! `options%direct_qp` and `options%direct_least_squares`

    integer,  parameter :: nn = 200
    real(wp), parameter :: f_star = -91.1200341_wp
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    real(wp) :: x_lb(2*nn+2), x_ub(2*nn+2), x0(2*nn+2), t
    integer :: istat, i

    ! variables: x_0..x_N are 1..N+1, y_0..y_N are N+2..2N+2; the two ends are fixed
    x_lb = -big
    x_ub = big
    x_lb([1, nn+2, 2*nn+2]) = 0.0_wp
    x_ub([1, nn+2, 2*nn+2]) = 0.0_wp
    x_lb(nn+1) = 1.0_wp
    x_ub(nn+1) = 1.0_wp
    call problem%set_problem_size(n=2*nn+2, m=nn)
    call problem%set_bounds(x_lb, x_ub, spread((2.0_wp/nn)**2, 1, nn), spread((2.0_wp/nn)**2, 1, nn))
    ! row i = link i, in (x_{i-1}, x_i, y_{i-1}, y_i)
    call problem%set_jacobian_sparsity(4*nn, [(i, i, i, i, i=1,nn)], [(i, i+1, nn+1+i, nn+2+i, i=1,nn)])
    ! the diagonal, then the subdiagonals of the x block and of the y block
    call problem%set_hessian_sparsity(2*nn+2 + 2*nn, [(i, i=1,2*nn+2), (i+1, i=1,nn), (nn+2+i, i=1,nn)], &
                                                     [(i, i=1,2*nn+2), (i, i=1,nn),   (nn+1+i, i=1,nn)])
    call problem%set_functions(fc=fc_chain, gjac=gjac_chain, hess=hess_chain)
    ! a shallow parabola between the ends
    do i = 0, nn
        t = real(i, wp)/real(nn, wp)
        x0(i+1)    = t
        x0(nn+2+i) = -0.4_wp*t*(1.0_wp - t)
    end do

    options%max_iter             = 1000
    options%hessian_mode         = sqpopt_hessian_exact
    options%inertia_control      = inertia
    options%direct_qp            = direct
    options%direct_least_squares = direct
    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_results(r)
    print '(3A,I0,A,I0,A,F14.8,A,ES9.2,A,ES9.2)', 'hanging chain [', label, ']: istat = ', istat, ', iterations = ', &
        r%iterations, ', f = ', r%f, ', KKT error = ', r%kkt_error, ', largest multiplier = ', maxval(abs(r%lambda))

    if (istat /= sqpopt_success) error stop 'test_multipliers FAILED: '//label//': the chain did not converge'
    if (abs(r%f - f_star) > 1.0e-5_wp*abs(f_star)) error stop 'test_multipliers FAILED: '//label//': wrong objective'
    if (maxval(abs(r%lambda)) > 1.0e5_wp) error stop 'test_multipliers FAILED: '//label//': the multipliers are too large'
    print '(A)', 'test_multipliers [hanging chain, '//label//'] PASSED'

    end subroutine test_chain

    subroutine fc_chain(x, f, c, status, data)
    !! the objective and the constraints
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    integer :: nn
    nn = size(c)
    f = sum(x(nn+2:))
    c = (x(2:nn+1) - x(1:nn))**2 + (x(nn+3:2*nn+2) - x(nn+2:2*nn+1))**2
    end subroutine fc_chain

    subroutine gjac_chain(x, g, jac_val, accuracy, status, data)
    !! the objective's gradient and the Jacobian's nonzeros
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` (in its sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! the accuracy asked for (the derivatives here are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: nn, i
    real(wp) :: dx, dy
    nn = size(jac_val)/4
    g(1:nn+1) = 0.0_wp
    g(nn+2:)  = 1.0_wp
    do i = 1, nn
        dx = x(i+1) - x(i)
        dy = x(nn+2+i) - x(nn+1+i)
        jac_val(4*i-3:4*i) = [-2.0_wp*dx, 2.0_wp*dx, -2.0_wp*dy, 2.0_wp*dy]
    end do
    end subroutine gjac_chain

    subroutine hess_chain(x, lambda, hess_val, status, data)
    !! the Hessian of the Lagrangian, `-sum_i lambda_i * hess(c_i)`: each link adds `-2*lambda_i` to
    !! the diagonal elements of its two joints, and `+2*lambda_i` between them, in the x and y blocks
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    integer :: nn, n, i
    nn = size(lambda)
    n  = 2*nn + 2
    hess_val = 0.0_wp
    do i = 1, nn
        hess_val([i, i+1, nn+1+i, nn+2+i]) = hess_val([i, i+1, nn+1+i, nn+2+i]) - 2.0_wp*lambda(i)
        hess_val(n+i)    = 2.0_wp*lambda(i)
        hess_val(n+nn+i) = 2.0_wp*lambda(i)
    end do
    end subroutine hess_chain

end program test_multipliers
