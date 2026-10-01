program test_inertia

    !! Test of the inertia control of the exact Hessian
    !! (`options%inertia_control`, see [[sqpopt_inertia_module]]), which
    !! needs a library built with MUMPS (the `HAS_MUMPS` preprocessor
    !! directive):
    !!
    !!    fpm test test_inertia --flag "-DHAS_MUMPS -I$CONDA_PREFIX/include" --link-flag "-ldmumps_seq"
    !!
    !! With MUMPS, it checks the shift that [[inertia_correct]] finds for small
    !! matrices whose inertia is known, and solves a nonconvex problem with
    !! inertia control (with both QP solvers, and with the trust region).
    !! Without MUMPS, it checks that the option is rejected as invalid input.

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_type, sqpopt_hessian_exact
    use sqpopt_inertia_module,   only: sqpopt_inertia_type, sqpopt_has_mumps
    use sqpopt_kkt_module,       only: sqpopt_kkt_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_dense, sqpopt_qp_reduced_hessian
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_types_module,     only: sqpopt_success, sqpopt_invalid_input, sqpopt_results_type, sqpopt_sparse_matrix
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    write(*,*) '----------------------------'
    write(*,*) 'test_inertia'
    write(*,*) '----------------------------'

    if (sqpopt_has_mumps) then
        call test_correction()
        call test_solve('dense QP',     sqpopt_qp_dense,           .false.)
        call test_solve('sparse QP',    sqpopt_qp_reduced_hessian, .false.)
        call test_solve('trust region', sqpopt_qp_dense,           .true.)
    else
        call test_unavailable()
    end if

    print '(A)', 'test_inertia PASSED'

    contains

    subroutine test_correction()
    !! the shift found for `H = diag(1,-1)` (and `diag(1,0)`), with the
    !! constraint row `J = [1 0]`, for several working sets

    type(sqpopt_inertia_type)  :: inertia
    type(sqpopt_kkt_type)      :: kkt
    type(sqpopt_hessian_type)  :: h
    type(sqpopt_sparse_matrix) :: jac
    logical :: changed, ok
    integer :: n_negative, n_factor0
    real(wp) :: shift

    call h%initialize(2, 1)
    call h%set_exact([1, 2], [1, 2])
    call h%set_values([1.0_wp, -1.0_wp], decay=.false.)
    jac%nrows = 1
    jac%ncols = 2
    jac%nnz   = 1
    jac%irow  = [1]
    jac%icol  = [1]
    jac%val   = [1.0_wp]

    call kkt%initialize(2, 1, jac%irow, jac%icol, ok, hess_irow=h%h_irow, hess_icol=h%h_icol)
    if (.not. (ok .and. kkt%enabled)) error stop 'test_inertia FAILED: MUMPS could not be started'
    inertia%enabled = .true.

    ! the constraint is in the working set: its null space is the second
    ! variable, where the curvature is -1, so the shift must exceed 1 (and
    ! the shift tried before it, an 8th of it, must not):
    call inertia%correct(kkt, h, jac, [1, 0, 0], changed, ok, n_negative)
    print '(A,ES10.2,A,I0,A)', 'shift = ', h%shift, ' (', kkt%solver%n_factor, ' factorizations)'
    if (.not. (ok .and. changed)) error stop 'test_inertia FAILED: negative curvature not corrected'
    if (n_negative /= 1) error stop 'test_inertia FAILED: wrong number of negative eigenvalues'
    if (.not. (h%shift > 1.0_wp .and. h%shift/8.0_wp <= 1.0_wp)) error stop 'test_inertia FAILED: wrong shift'
    if (inertia%shift_last /= h%shift) error stop 'test_inertia FAILED: shift_last not set'
    shift = h%shift

    ! the same working set and shift again: no new factorization
    n_factor0 = kkt%solver%n_factor
    call inertia%correct(kkt, h, jac, [-1, 0, 0], changed, ok)
    if (.not. ok .or. changed .or. kkt%solver%n_factor /= n_factor0) error stop 'test_inertia FAILED: outcome not reused'

    ! new matrices: the search restarts from a third of the last shift
    ! (two factorizations: the zero shift, then that one, which is enough)
    call kkt%new_matrices()
    h%shift = 0.0_wp
    call inertia%correct(kkt, h, jac, [1, 0, 0], changed, ok)
    if (.not. (ok .and. changed) .or. kkt%solver%n_factor /= n_factor0 + 2) error stop 'test_inertia FAILED: restart'
    if (abs(h%shift - shift/3.0_wp) > 1.0e-12_wp) error stop 'test_inertia FAILED: restart shift'

    ! the second variable is at a bound of the working set: nothing is left
    ! of the null space, so no shift is needed
    call kkt%new_matrices()
    h%shift = 0.0_wp
    call inertia%correct(kkt, h, jac, [1, 0, 1], changed, ok, n_negative)
    if (.not. ok .or. changed .or. n_negative /= 0) error stop 'test_inertia FAILED: variable at a bound'

    ! an empty working set: the Hessian itself is indefinite
    call inertia%correct(kkt, h, jac, [0, 0, 0], changed, ok, n_negative)
    if (.not. (ok .and. changed) .or. n_negative /= 1 .or. .not. h%shift > 1.0_wp) then
        error stop 'test_inertia FAILED: empty working set'
    end if

    ! a shift that may not exceed shift_max can't correct it
    call kkt%new_matrices()
    h%shift     = 0.0_wp
    h%shift_max = 0.5_wp
    call inertia%correct(kkt, h, jac, [1, 0, 0], changed, ok)
    if (ok .or. .not. changed .or. h%shift /= 0.5_wp) error stop 'test_inertia FAILED: shift_max'
    h%shift_max = 1.0e10_wp

    ! zero curvature (a singular KKT matrix) is not negative curvature
    call kkt%new_matrices()
    call h%set_values([1.0_wp, 0.0_wp], decay=.false.)
    h%shift = 0.0_wp
    call inertia%correct(kkt, h, jac, [1, 0, 0], changed, ok, n_negative)
    if (.not. ok .or. changed .or. n_negative /= 0) error stop 'test_inertia FAILED: zero curvature'

    ! without the KKT matrix, a correction does nothing, and inertia control turns itself off
    call kkt%destroy()
    if (kkt%enabled .or. kkt%solver%n_factor /= 0) error stop 'test_inertia FAILED: destroy'
    call inertia%correct(kkt, h, jac, [1, 0, 0], changed, ok)
    if (ok .or. changed .or. inertia%enabled) error stop 'test_inertia FAILED: correction after destroy'

    print '(A)', 'test_inertia [correction] PASSED'

    end subroutine test_correction

    subroutine test_solve(label, qp_mode, use_trust_region)
    !! minimize \( -x_1x_2 \) subject to \( x_1^2 + x_2^2 \le 1 \), \( 0 \le x \le 2 \),
    !! whose Hessian is indefinite (the solution is \( x_1 = x_2 = 1/\sqrt{2} \), `f = -1/2`),
    !! with inertia control, twice (the second solve must repeat the first)

    character(len=*), intent(in) :: label            !! the configuration, for the output
    integer,          intent(in) :: qp_mode          !! `options%qp_solver_mode`
    logical,          intent(in) :: use_trust_region !! `trust_region%enabled`

    type(sqpopt_type)              :: solver
    type(sqpopt_problem_type)      :: problem
    type(sqpopt_options_type)      :: options
    type(sqpopt_trust_region_type) :: trust_region
    type(sqpopt_results_type)      :: r, r2
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp, 0.0_wp], x_ub=[2.0_wp, 2.0_wp], c_lb=[-big], c_ub=[1.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_hessian_sparsity(nnz=3, irow=[1,2,2], icol=[1,1,2])
    call problem%set_functions(fc=fc, gjac=gjac, hess=hess)

    options%hessian_mode    = sqpopt_hessian_exact
    options%inertia_control = .true.
    options%qp_solver_mode  = qp_mode
    trust_region%enabled    = use_trust_region

    call solver%initialize(problem=problem, options=options, trust_region=trust_region)
    call solver%solve([0.3_wp, 0.6_wp], istat)
    call solver%get_results(r)
    print '(3A,I0,A,F12.8,3(A,I0))', 'test_inertia [', label, ']: istat = ', istat, ', f = ', r%f, ', iterations = ', &
        r%iterations, ', factorizations = ', r%n_factorizations, ', shift increases = ', r%n_hessian_resets

    if (istat /= sqpopt_success .or. r%istat /= sqpopt_success) error stop 'test_inertia FAILED: '//label//' did not converge'
    if (abs(r%f + 0.5_wp) > 1.0e-6_wp .or. maxval(abs(r%x - sqrt(0.5_wp))) > 1.0e-4_wp) then
        error stop 'test_inertia FAILED: '//label//' wrong solution'
    end if
    ! (every iteration that solves a QP factors at least once, and the
    ! indefinite Hessian at the starting point must have been shifted)
    if (r%n_factorizations < r%iterations - 1) error stop 'test_inertia FAILED: '//label//' too few factorizations'
    if (r%n_hessian_resets < 1) error stop 'test_inertia FAILED: '//label//' the Hessian was never shifted'

    ! no state carries over to a second solve:
    call solver%solve([0.3_wp, 0.6_wp], istat)
    call solver%get_results(r2)
    if (istat /= r%istat .or. r2%iterations /= r%iterations .or. r2%n_factorizations /= r%n_factorizations .or. &
        r2%n_eval_fc /= r%n_eval_fc .or. any(r2%x /= r%x)) error stop 'test_inertia FAILED: '//label//' second solve differs'

    print '(A)', 'test_inertia ['//label//'] PASSED'

    end subroutine test_solve

    subroutine test_unavailable()
    !! without MUMPS, `options%inertia_control` is invalid input

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    integer :: istat

    call problem%set_problem_size(n=2, m=1)
    call problem%set_bounds(x_lb=[0.0_wp, 0.0_wp], x_ub=[2.0_wp, 2.0_wp], c_lb=[-big], c_ub=[1.0_wp])
    call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
    call problem%set_hessian_sparsity(nnz=3, irow=[1,2,2], icol=[1,1,2])
    call problem%set_functions(fc=fc, gjac=gjac, hess=hess)

    options%hessian_mode    = sqpopt_hessian_exact
    options%inertia_control = .true.
    call solver%initialize(problem=problem, options=options)
    call solver%solve([0.3_wp, 0.6_wp], istat)
    call solver%get_results(r)
    print '(A)', 'without MUMPS: '//solver%status_message()
    if (istat /= sqpopt_invalid_input .or. r%istat /= sqpopt_invalid_input) then
        error stop 'test_inertia FAILED: inertia_control was accepted without MUMPS'
    end if
    if (index(solver%status_message(), 'MUMPS') == 0) error stop 'test_inertia FAILED: the message does not say why'

    print '(A)', 'test_inertia [not built with MUMPS: option rejected] PASSED'

    end subroutine test_unavailable

    subroutine fc(x, f, c, status, data)
    !! the objective and the constraint
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f    = -x(1)*x(2)
    c(1) = x(1)**2 + x(2)**2
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective's gradient and the Jacobian's nonzeros
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` (in its sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! the accuracy asked for (the derivatives here are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    g       = [-x(2), -x(1)]
    jac_val = [2.0_wp*x(1), 2.0_wp*x(2)]
    end subroutine gjac

    subroutine hess(x, lambda, hess_val, status, data)
    !! the Hessian of the Lagrangian \( \nabla^2 f - \lambda_1 \nabla^2 c_1 \), lower triangle
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    ! (1,1) (2,1) (2,2):
    hess_val = [-2.0_wp*lambda(1), -1.0_wp, -2.0_wp*lambda(1)]
    end subroutine hess

end program test_inertia
