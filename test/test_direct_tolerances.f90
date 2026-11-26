program test_direct_tolerances

    !! Test of the direct QP method's opt-in tolerance settings
    !! (`qp_solver%direct_relative_tol` and `qp_solver%direct_accept_inexact`, the `relative` and
    !! `accept_inexact` arguments of [[direct_qp_step]]):
    !!
    !! * a QP whose step is very long: minimize \( p^T H p/2 + g^T p \) subject to
    !!   \( p_1 + 2 p_2 = 0.7 \), with \( H \) nearly singular along the row's null space \( (2,-1) \),
    !!   so the step is about \( 10^{10} \) and the row can only be met to roundoff relative to
    !!   that. With `relative`, the row is tested relative to \( |J_1| \, |p| \), and the face
    !!   is solved; the row's residual must be within the tolerance of that size.
    !! * the same QP with an impossible tolerance (`1e-30`): without `accept_inexact`, the face
    !!   is inconsistent and the method gives up (`sqpopt_direct_singular`); with it, the step
    !!   is returned with `sqpopt_direct_solved_inexact`, and it must be the long step.
    !! * Hock-Schittkowski problem 71, solved with `direct_qp` and both settings on, must
    !!   converge to its known solution (the settings don't change a well-scaled solve).

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_type
    use sqpopt_kkt_module,       only: sqpopt_kkt_type
    use sqpopt_qp_direct_module, only: direct_qp_step, sqpopt_direct_solved, sqpopt_direct_singular, &
                                       sqpopt_direct_solved_inexact
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type
    use sqpopt_types_module,     only: sqpopt_success, sqpopt_sparse_matrix
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    call test_long_step()
    call test_hs71()
    print '(A)', 'test_direct_tolerances PASSED'

    contains

    subroutine test_long_step()
    !! the direct method on the long-step QP, with each setting
    type(sqpopt_kkt_type)      :: kkt
    type(sqpopt_hessian_type)  :: h
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: p(2), lambda(1), size_terms, residual
    integer  :: status(3), n_changes, outcome
    logical  :: ok
    real(wp), parameter :: x(2) = 0.0_wp, c(1) = 0.0_wp, g(2) = [-2.3_wp, 1.1_wp], b = 0.7_wp
    real(wp), parameter :: tol = 1.0e-8_wp

    call h%initialize(2, 1)
    call h%set_exact([1, 2, 2], [1, 1, 2])
    call h%set_values([1.0_wp + 1.0e-10_wp, 2.0_wp, 4.0_wp + 1.0e-10_wp], decay=.false.)
    jac%nrows = 1
    jac%ncols = 2
    jac%nnz   = 2
    jac%irow  = [1, 1]
    jac%icol  = [1, 2]
    jac%val   = [1.0_wp, 2.0_wp]
    call kkt%initialize(2, 1, jac%irow, jac%icol, ok, hess_irow=h%h_irow, hess_icol=h%h_icol)
    if (.not. ok) error stop 'test_direct_tolerances FAILED: the KKT matrix could not be set up'

    ! the default test (absolute in the row's target): only reported
    status = [-1, 0, 0]
    call direct_qp_step(kkt, h, jac, x, g, c, [-big, -big], [big, big], [b], [b], &
                        10, tol, status, p, lambda, n_changes, outcome)
    print '(A,I0,A,2ES12.4,A,ES10.2)', 'long step, default test: outcome ', outcome, ', p = ', p, &
                                       ', row residual ', abs(p(1) + 2.0_wp*p(2) - b)

    ! relative to the terms: solved, with the row met to the tolerance of their size
    call kkt%new_matrices()
    status = [-1, 0, 0]
    call direct_qp_step(kkt, h, jac, x, g, c, [-big, -big], [big, big], [b], [b], &
                        10, tol, status, p, lambda, n_changes, outcome, relative=.true.)
    size_terms = abs(p(1)) + 2.0_wp*abs(p(2))
    residual   = abs(p(1) + 2.0_wp*p(2) - b)
    print '(A,I0,A,2ES12.4,A,ES10.2,A,ES10.2)', 'long step, relative test: outcome ', outcome, ', p = ', p, &
                                                ', row residual ', residual, ' of terms ', size_terms
    if (outcome /= sqpopt_direct_solved) error stop 'test_direct_tolerances FAILED: relative test, not solved'
    if (.not. size_terms > 1.0e9_wp) error stop 'test_direct_tolerances FAILED: the step is not long'
    if (residual > tol*size_terms) error stop 'test_direct_tolerances FAILED: the row is not met relative to its terms'

    ! an impossible tolerance: given up without accept_inexact ...
    call kkt%new_matrices()
    status = [-1, 0, 0]
    call direct_qp_step(kkt, h, jac, x, g, c, [-big, -big], [big, big], [b], [b], &
                        10, 1.0e-30_wp, status, p, lambda, n_changes, outcome)
    print '(A,I0)', 'impossible tolerance, default: outcome ', outcome
    if (outcome /= sqpopt_direct_singular) error stop 'test_direct_tolerances FAILED: an inconsistent face was accepted'
    if (any(status /= [-1, 0, 0])) error stop 'test_direct_tolerances FAILED: a failed step changed the working set'

    ! ... and accepted as inexact with it
    call kkt%new_matrices()
    status = [-1, 0, 0]
    call direct_qp_step(kkt, h, jac, x, g, c, [-big, -big], [big, big], [b], [b], &
                        10, 1.0e-30_wp, status, p, lambda, n_changes, outcome, accept_inexact=.true.)
    print '(A,I0,A,2ES12.4)', 'impossible tolerance, accept_inexact: outcome ', outcome, ', p = ', p
    if (outcome /= sqpopt_direct_solved_inexact) error stop 'test_direct_tolerances FAILED: not accepted as inexact'
    if (.not. abs(p(1)) + abs(p(2)) > 1.0e9_wp) error stop 'test_direct_tolerances FAILED: the inexact step is not the long step'
    if (any(status /= [-1, 0, 0])) error stop 'test_direct_tolerances FAILED: wrong working set'

    call kkt%destroy()

    end subroutine test_long_step

    subroutine test_hs71()
    !! HS71 with `direct_qp` and both settings on
    type(sqpopt_problem_type)   :: problem
    type(sqpopt_options_type)   :: options
    type(sqpopt_qp_solver_type) :: qp_solver
    type(sqpopt_type)           :: solver
    real(wp) :: x(4), lambda(2)
    integer  :: istat
    real(wp), parameter :: x_star(4) = [1.0_wp, 4.74299963_wp, 3.82114998_wp, 1.37940829_wp]

    call problem%set_problem_size(n=4, m=2)
    call problem%set_bounds(x_lb=[1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp], x_ub=[5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp], &
                            c_lb=[40.0_wp, 25.0_wp], c_ub=[40.0_wp, big])
    call problem%set_jacobian_sparsity(nnz=8, irow=[1,1,1,1,2,2,2,2], icol=[1,2,3,4,1,2,3,4])
    call problem%set_functions(fc=fc, gjac=gjac)
    options%direct_qp = .true.
    qp_solver%direct_relative_tol   = .true.
    qp_solver%direct_accept_inexact = .true.
    call solver%initialize(problem=problem, options=options, qp_solver=qp_solver)
    call solver%solve([1.0_wp, 5.0_wp, 5.0_wp, 1.0_wp], istat)
    call solver%get_solution(x, lambda)
    print '(A,I0,A,4F12.8)', 'hs71 with both settings: istat ', istat, ', x = ', x
    if (istat /= sqpopt_success) error stop 'test_direct_tolerances FAILED: hs71 did not converge'
    if (maxval(abs(x - x_star)) > 1.0e-6_wp) error stop 'test_direct_tolerances FAILED: hs71 wrong solution'

    end subroutine test_hs71

    subroutine fc(x, f, c, status, data)
    !! HS71's objective and constraints
    real(wp), dimension(:), intent(in)    :: x      !! the variables `dimension(4)`
    real(wp),               intent(out)   :: f      !! the objective
    real(wp), dimension(:), intent(out)   :: c      !! the constraints `dimension(2)`
    integer,                intent(inout) :: status !! (not used)
    class(*), optional,     intent(inout) :: data   !! (not used)
    f = x(1)*x(4)*(x(1)+x(2)+x(3)) + x(3)
    c = [sum(x**2), product(x)]
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! HS71's gradient and Jacobian
    real(wp), dimension(:), intent(in)    :: x        !! the variables `dimension(4)`
    real(wp), dimension(:), intent(out)   :: g        !! the gradient `dimension(4)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! the Jacobian's nonzeros `dimension(8)`
    integer,                intent(in)    :: accuracy !! (not used)
    integer,                intent(inout) :: status   !! (not used)
    class(*), optional,     intent(inout) :: data     !! (not used)
    g = [x(4)*(2.0_wp*x(1)+x(2)+x(3)), x(1)*x(4), x(1)*x(4) + 1.0_wp, x(1)*(x(1)+x(2)+x(3))]
    jac_val(1:4) = 2.0_wp*x
    jac_val(5:8) = [x(2)*x(3)*x(4), x(1)*x(3)*x(4), x(1)*x(2)*x(4), x(1)*x(2)*x(3)]
    end subroutine gjac

end program test_direct_tolerances
