program test_maratos

    !! The classic Maratos-effect example (Nocedal & Wright, *Numerical
    !! Optimization*, 2nd ed., Example 15.4):
    !!
    !!   minimize   2*(x1^2 + x2^2 - 1) - x1
    !!   subject to x1^2 + x2^2 = 1
    !!
    !! solution x* = (1, 0), f* = -1. Started on the constraint close to the
    !! solution, at x0 = (cos(0.3), sin(0.3)): full SQP steps there increase
    !! both the objective and the constraint violation, so a merit function
    !! or filter rejects them (and plain backtracking converges slowly) unless
    !! the second-order correction is applied. Every line search mode, and the
    !! trust region, must converge, within a modest evaluation budget.

    use sqpopt_module,              only: sqpopt_type
    use sqpopt_problem_module,      only: sqpopt_problem_type
    use sqpopt_options_module,      only: sqpopt_options_type
    use sqpopt_linesearch_module,   only: sqpopt_linesearch_armijo, sqpopt_linesearch_watchdog, &
                                          sqpopt_linesearch_filter, sqpopt_merit_augmented_lagrangian
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_types_module,        only: sqpopt_success
    use sqpopt_kinds,               only: wp => sqpopt_module_wp

    implicit none

    integer,  parameter :: ls_modes(3) = [sqpopt_linesearch_armijo, sqpopt_linesearch_watchdog, sqpopt_linesearch_filter]
    integer,  parameter :: max_evals = 100  !! evaluation budget for each run
    real(wp), parameter :: xexpect(2) = [1.0_wp, 0.0_wp]
    real(wp), parameter :: theta = 0.3_wp

    type(sqpopt_type)              :: solver
    type(sqpopt_problem_type)      :: problem
    type(sqpopt_options_type)      :: options
    type(sqpopt_trust_region_type) :: trust_region
    real(wp) :: xsol(2), lam(1)
    integer  :: istat, i, tr, merit
    integer  :: n_f  !! number of objective evaluations

    write(*,*) '----------------------------'
    write(*,*) 'test_maratos'
    write(*,*) '----------------------------'

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[-10.0_wp,-10.0_wp], x_ub=[10.0_wp,10.0_wp], c_lb=[0.0_wp], c_ub=[0.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_functions(fc=fc_obj_cons, gjac=gjac_grad_jacv)

    do tr = 0, 1
        do merit = 1, 2
            do i = 1, size(ls_modes)
                if (merit == 2 .and. ls_modes(i) == sqpopt_linesearch_filter) cycle ! the filter has no merit function
                options = sqpopt_options_type()
                options%linesearch_mode = ls_modes(i)
                if (merit == 2) options%merit_mode = sqpopt_merit_augmented_lagrangian
                trust_region = sqpopt_trust_region_type()
                trust_region%enabled = tr == 1
                call solver%initialize(problem=problem, options=options, trust_region=trust_region)
                n_f = 0
                call solver%solve([cos(theta), sin(theta)], istat)
                call solver%get_solution(xsol, lam)
                print '(A,L1,A,I0,A,I0,A,2F12.8,A,I0,A,I0)', 'trust_region=', tr == 1, ' merit_mode=', options%merit_mode, &
                    ' linesearch_mode=', ls_modes(i), ': x=', xsol, '  istat=', istat, '  n_f=', n_f
                if (istat /= sqpopt_success) error stop 'test_maratos FAILED: did not reach sqpopt_success'
                if (maxval(abs(xsol-xexpect)) > 1.0e-6_wp) error stop 'test_maratos FAILED: wrong solution'
                if (n_f > max_evals) error stop 'test_maratos FAILED: too many function evaluations'
            end do
        end do
    end do

    print '(A)', 'test_maratos PASSED'

    contains

    subroutine obj(x, f, status, data)
    !! the objective
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    n_f = n_f + 1
    f = 2.0_wp*(x(1)**2 + x(2)**2 - 1.0_wp) - x(1)
    end subroutine obj

    subroutine grad(x, g, status, data)
    !! the objective's gradient
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g      !! objective gradient at `x` `dimension(n)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    g = [4.0_wp*x(1) - 1.0_wp, 4.0_wp*x(2)]
    end subroutine grad

    subroutine cons(x, c, status, data)
    !! the constraints
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    c(1) = x(1)**2 + x(2)**2 - 1.0_wp
    end subroutine cons

    subroutine jacv(x, jac_val, status, data)
    !! the nonzero values of the constraint Jacobian (in the sparsity pattern's order)
    real(wp), dimension(:), intent(in)    :: x       !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero values of the constraint Jacobian at `x` (in the sparsity pattern's order)
    integer,                intent(inout) :: status  !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data    !! the user data passed to `set_functions` (if any)
    jac_val = 2.0_wp*x
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


end program test_maratos
