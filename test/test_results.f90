program test_results

    !! Tests of the results reported by a solve (`get_results`/`get_solution`):
    !!
    !!   minimize   s*((x1-2)^2 + (x2-3)^2)
    !!   subject to s*(x1 + x2) <= 3*s,   x1 <= 0.5
    !!
    !! The solution is x* = (0.5, 2.5), with both the constraint and the
    !! bound on x1 active. From g = J^T lambda + z: lambda = -1 (at the
    !! constraint's upper bound, so <= 0) and z = (-2*s, 0) (x1 at its upper
    !! bound, so <= 0) -- independent of `s` for lambda, since the
    !! constraint is scaled by the same factor as the objective.
    !!
    !! Checked, for s=1 and for s=1e4 (where the automatic gradient-based
    !! scaling is active, so the solver works with a rescaled problem and
    !! must unscale what it reports), with scaling on and off:
    !!
    !! * the solution, `lambda`, and `z`;
    !! * the evaluation counts against the user functions' own counters;
    !! * a warm start from `x*` and `lambda0 = lambda*` converges at once.

    use sqpopt_module,         only: sqpopt_type
    use sqpopt_problem_module, only: sqpopt_problem_type
    use sqpopt_options_module, only: sqpopt_options_type
    use sqpopt_types_module,   only: sqpopt_success, sqpopt_results_type
    use sqpopt_kinds,          only: wp => sqpopt_module_wp

    implicit none

    real(wp) :: s   !! problem scale
    integer  :: n_f, n_g, n_c, n_jac
    integer  :: k, sc

    write(*,*) '----------------------------'
    write(*,*) 'test_results'
    write(*,*) '----------------------------'

    do k = 1, 2
        s = merge(1.0_wp, 1.0e4_wp, k == 1)
        do sc = 0, 1
            call run(s, sc == 1)
        end do
    end do

    print '(A)', 'test_results PASSED'

    contains

    subroutine run(s_in, scaling)
    !! solve the problem with scale `s_in`, and check the reported results
    real(wp), intent(in) :: s_in    !! the problem's scale `s` (see the program docs)
    logical,  intent(in) :: scaling !! `options%scaling`
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    real(wp) :: x(2), lam(1), z(2)
    integer :: istat
    real(wp), parameter :: tol = 1.0e-6_wp

    s = s_in
    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-10.0_wp,-10.0_wp], x_ub=[0.5_wp,10.0_wp], c_lb=[-1.0e20_wp], c_ub=[3.0_wp*s])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv)
    options%scaling = scaling

    n_f = 0; n_g = 0; n_c = 0; n_jac = 0
    call solver%initialize(problem=problem, options=options)
    call solver%solve([0.0_wp, 0.0_wp], istat)
    call solver%get_solution(x, lam, z)
    call solver%get_results(r)
    print '(A,ES8.1,A,L1,A,2F10.6,A,F10.6,A,2ES12.4,A,I0,A,4I4)', 's=', s, ' scaling=', scaling, ' x=', x, &
        ' lambda=', lam, ' z=', z, ' iter=', r%iterations, ' evals=', r%n_eval_fc, r%n_eval_gjac
    if (istat /= sqpopt_success) error stop 'test_results FAILED: did not converge'
    if (maxval(abs(x - [0.5_wp, 2.5_wp])) > tol) error stop 'test_results FAILED: wrong x'
    if (abs(lam(1) + 1.0_wp) > tol) error stop 'test_results FAILED: wrong lambda'
    if (abs(z(1) + 2.0_wp*s) > tol*s .or. abs(z(2)) > tol*s) error stop 'test_results FAILED: wrong z'
    if (abs(r%f - s*2.5_wp) > tol*s) error stop 'test_results FAILED: wrong f'
    if (abs(r%c(1) - 3.0_wp*s) > tol*s) error stop 'test_results FAILED: wrong c'
    if (r%feasibility_error > tol*s) error stop 'test_results FAILED: feasibility error'
    if (any(r%x /= x) .or. any(r%lambda /= lam) .or. any(r%z /= z)) error stop 'test_results FAILED: inconsistent results'
    ! (each `fc` call calls `obj` and `cons` once, and each `gjac` call `grad` and `jacv`)
    if (r%n_eval_fc /= n_f .or. r%n_eval_fc /= n_c .or. r%n_eval_gjac /= n_g .or. r%n_eval_gjac /= n_jac) then
        error stop 'test_results FAILED: evaluation counts'
    end if

    ! warm start from the solution and its multipliers: converged at the first iterate
    options%hessian_scale0 = 2.0_wp
    call solver%initialize(problem=problem, options=options)
    call solver%solve(x, istat, lambda0=lam)
    call solver%get_results(r)
    print '(A,I0,A,I0)', '   warm start: istat=', istat, ' iterations=', r%iterations
    if (istat /= sqpopt_success .or. r%iterations /= 1) error stop 'test_results FAILED: warm start'
    end subroutine run

    subroutine obj(x, f, status, data)
    !! the objective
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    n_f = n_f + 1
    f = s*((x(1)-2.0_wp)**2 + (x(2)-3.0_wp)**2)
    end subroutine obj

    subroutine grad(x, g, status, data)
    !! the objective's gradient
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    n_g = n_g + 1
    g = s*[2.0_wp*(x(1)-2.0_wp), 2.0_wp*(x(2)-3.0_wp)]
    end subroutine grad

    subroutine cons(x, c, status, data)
    !! the constraints
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    n_c = n_c + 1
    c(1) = s*(x(1) + x(2))
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    !! the nonzero values of the constraint Jacobian (in the sparsity pattern's order)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    n_jac = n_jac + 1
    jac_val = s
    end subroutine jacv

    subroutine fc_obj_cons(x, f, c, status, data)
    !! `fc` for `set_functions`: the objective (`obj`) and the constraints (`cons`)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    call obj(x, f, status, data)
    if (status == 0) call cons(x, c, status, data)
    end subroutine fc_obj_cons

    subroutine gjac_grad_jacv(x, g, jac_val, accuracy, status, data)
    !! `gjac` for `set_functions`: the gradient (`grad`) and the Jacobian values (`jacv`)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g       !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    call grad(x, g, status, data)
    if (status == 0) call jacv(x, jac_val, status, data)
    end subroutine gjac_grad_jacv


end program test_results
