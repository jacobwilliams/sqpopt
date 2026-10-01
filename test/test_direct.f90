program test_direct

    !! Test of the options that solve by sparse factorizations
    !! (`options%direct_qp` and `options%direct_least_squares`, with and
    !! without `options%inertia_control`), which need a library built with
    !! MUMPS (the `HAS_MUMPS` preprocessor directive):
    !!
    !!    fpm test test_direct --flag "-DHAS_MUMPS -I$CONDA_PREFIX/include" --link-flag "-ldmumps_seq"
    !!
    !! With MUMPS, two small problems are solved with each Hessian mode and
    !! several combinations of those options, with both QP solvers, with the
    !! trust region, and on two OpenMP threads
    !! (`options%factorization_threads`). Each run must converge to the known solution, solve
    !! QPs directly (if `direct_qp`), factor matrices, and give the same
    !! results when solved again (no state carries over):
    !!
    !! * Hock-Schittkowski problem 71;
    !! * the Maratos example (see `test_maratos`), whose second-order
    !!   corrections use the direct least-squares solver.
    !!
    !! It also gives [[direct_qp_step]] small QPs that take each of its special
    !! paths: a nonconvex face (which ends it without inertia control, and
    !! raises the Hessian's shift with it), a face whose KKT matrix is
    !! singular (which is regularized), an infeasible QP, and the limit on
    !! the changes of the working set.
    !!
    !! Without MUMPS, it checks that each option is rejected as invalid input.

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_hessian_module,   only: sqpopt_hessian_bfgs, sqpopt_hessian_sr1, sqpopt_hessian_exact, sqpopt_hessian_type
    use sqpopt_kkt_module,       only: sqpopt_kkt_type
    use sqpopt_inertia_module,   only: sqpopt_inertia_type
    use sqpopt_qp_direct_module, only: direct_qp_step, sqpopt_direct_solved, sqpopt_direct_nonconvex, &
                                       sqpopt_direct_singular, sqpopt_direct_max_changes
    use sqpopt_qp_solver_module, only: sqpopt_qp_auto, sqpopt_qp_reduced_hessian
    use sqpopt_trust_region_module, only: sqpopt_trust_region_type
    use sqpopt_symmetric_solver_module, only: sqpopt_has_mumps
    use sqpopt_types_module,     only: sqpopt_success, sqpopt_acceptable, sqpopt_stalled, sqpopt_invalid_input, &
                                       sqpopt_results_type, sqpopt_sparse_matrix
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    type :: config_type
        !! one combination of options
        character(len=32) :: name = ''                !! for the output
        integer :: hessian_mode = sqpopt_hessian_bfgs !! `options%hessian_mode`
        logical :: inertia   = .false.                !! `options%inertia_control`
        logical :: direct    = .false.                !! `options%direct_qp`
        logical :: direct_ls = .false.                !! `options%direct_least_squares`
        integer :: threads = 1                        !! `options%factorization_threads`
        logical :: trust_region = .false.             !! `trust_region%enabled`
        integer :: qp_mode = sqpopt_qp_auto           !! `options%qp_solver_mode`
    end type config_type

    type(config_type), parameter :: configs(11) = [ &
        config_type(name='exact, direct', hessian_mode=sqpopt_hessian_exact, direct=.true.), &
        config_type(name='exact, inertia, direct', hessian_mode=sqpopt_hessian_exact, inertia=.true., direct=.true.), &
        config_type(name='exact, inertia, direct, sparse', hessian_mode=sqpopt_hessian_exact, inertia=.true., &
                    direct=.true., qp_mode=sqpopt_qp_reduced_hessian), &
        config_type(name='exact, inertia, direct, TR', hessian_mode=sqpopt_hessian_exact, inertia=.true., &
                    direct=.true., trust_region=.true.), &
        config_type(name='L-BFGS, direct', direct=.true.), &
        config_type(name='L-BFGS, direct, TR', direct=.true., trust_region=.true.), &
        config_type(name='SR1, inertia', hessian_mode=sqpopt_hessian_sr1, inertia=.true.), &
        config_type(name='SR1, inertia, direct', hessian_mode=sqpopt_hessian_sr1, inertia=.true., direct=.true.), &
        config_type(name='L-BFGS, least squares', direct_ls=.true.), &
        config_type(name='exact, all three', hessian_mode=sqpopt_hessian_exact, inertia=.true., direct=.true., &
                    direct_ls=.true.), &
        config_type(name='exact, all three, 2 threads', hessian_mode=sqpopt_hessian_exact, inertia=.true., &
                    direct=.true., direct_ls=.true., threads=2) ]

    integer :: i

    write(*,*) '----------------------------'
    write(*,*) 'test_direct'
    write(*,*) '----------------------------'

    if (sqpopt_has_mumps) then
        call test_direct_step()
        do i = 1, size(configs)
            call run('hs71', configs(i))
            call run('maratos', configs(i))
        end do
    else
        call test_unavailable()
    end if

    print '(A)', 'test_direct PASSED'

    contains

    subroutine test_direct_step()
    !! the special paths of the direct QP method, on QPs with two variables
    !! (`x = 0`, so the step's bounds are the variables')
    type(sqpopt_kkt_type)      :: kkt
    type(sqpopt_inertia_type)  :: inertia
    type(sqpopt_hessian_type)  :: h
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: p(2), lambda(1)
    integer  :: status(3), n_changes, outcome
    logical  :: ok
    real(wp), parameter :: x(2) = 0.0_wp, c(1) = 0.0_wp, tol = 1.0e-8_wp

    call h%initialize(2, 1)
    call h%set_exact([1, 2], [1, 2])
    jac%nrows = 1
    jac%ncols = 2
    jac%nnz   = 2
    jac%irow  = [1, 1]
    jac%icol  = [1, 2]
    jac%val   = [1.0_wp, -1.0_wp]
    call kkt%initialize(2, 1, jac%irow, jac%icol, ok, hess_irow=h%h_irow, hess_icol=h%h_icol)
    if (.not. ok) error stop 'test_direct FAILED: the KKT matrix could not be set up'

    ! ---- a nonconvex face: H = diag(1,-1), no constraint in the working set ----
    call h%set_values([1.0_wp, -1.0_wp], decay=.false.)
    status = 0
    call direct_qp_step(kkt, h, jac, x, [0.0_wp, 0.5_wp], c, [-1.0_wp, -1.0_wp], [1.0_wp, 1.0_wp], [-big], [big], &
                        10, tol, status, p, lambda, n_changes, outcome)
    if (outcome /= sqpopt_direct_nonconvex) error stop 'test_direct FAILED: nonconvex face not reported'
    if (any(status /= 0) .or. h%shift /= 0.0_wp) error stop 'test_direct FAILED: a failed direct step changed its inputs'
    ! with inertia control, the shift is raised (beyond 1) and the QP solved: p2 = -0.5/(-1 + shift)
    inertia%enabled = .true.
    call direct_qp_step(kkt, h, jac, x, [0.0_wp, 0.5_wp], c, [-1.0_wp, -1.0_wp], [1.0_wp, 1.0_wp], [-big], [big], &
                        10, tol, status, p, lambda, n_changes, outcome, inertia=inertia)
    print '(A,ES10.2,A,2F10.6)', 'nonconvex face: shift = ', h%shift, ', p = ', p
    if (outcome /= sqpopt_direct_solved .or. .not. h%shift > 1.0_wp) error stop 'test_direct FAILED: shift not raised'
    if (abs(p(1)) > 1.0e-10_wp .or. abs(p(2) + 0.5_wp/(h%shift - 1.0_wp)) > 1.0e-10_wp) then
        error stop 'test_direct FAILED: wrong step on the shifted face'
    end if

    ! ---- a singular face: minimize |p|^2/2 - 10(p1+p2) with p1 - p2 = 0 and p <= 0.5. The step on
    !      the row's face is (10,10), so both bounds are added, and the row has no free variable
    !      left. The solution is (0.5,0.5), where the row holds. ----
    call kkt%new_matrices()
    h%shift = 0.0_wp
    call h%set_values([1.0_wp, 1.0_wp], decay=.false.)
    status = [-1, 0, 0]
    call direct_qp_step(kkt, h, jac, x, [-10.0_wp, -10.0_wp], c, [-1.0_wp, -1.0_wp], [0.5_wp, 0.5_wp], &
                        [0.0_wp], [0.0_wp], 10, tol, status, p, lambda, n_changes, outcome)
    print '(A,I0,A,2F10.6,A,3I3)', 'singular face: ', n_changes, ' change(s), p = ', p, ', working set ', status
    if (outcome /= sqpopt_direct_solved) error stop 'test_direct FAILED: singular face not solved'
    if (maxval(abs(p - 0.5_wp)) > 1.0e-10_wp .or. any(status(2:3) /= 1)) error stop 'test_direct FAILED: singular face'

    ! ---- the limit on the changes: the same QP needs one, and none is allowed ----
    call kkt%new_matrices()
    status = [-1, 0, 0]
    call direct_qp_step(kkt, h, jac, x, [-10.0_wp, -10.0_wp], c, [-1.0_wp, -1.0_wp], [0.5_wp, 0.5_wp], &
                        [0.0_wp], [0.0_wp], 0, tol, status, p, lambda, n_changes, outcome)
    if (outcome /= sqpopt_direct_max_changes .or. any(status /= [-1, 0, 0])) then
        error stop 'test_direct FAILED: the limit on the changes'
    end if

    ! ---- an infeasible QP: p1 - p2 = 2 can't hold with p1 <= 0.5 and p2 >= -1 ----
    call kkt%new_matrices()
    status = [-1, 0, 0]
    call direct_qp_step(kkt, h, jac, x, [-10.0_wp, 10.0_wp], c, [-1.0_wp, -1.0_wp], [0.5_wp, 0.5_wp], &
                        [2.0_wp], [2.0_wp], 10, tol, status, p, lambda, n_changes, outcome)
    if (outcome == sqpopt_direct_solved) error stop 'test_direct FAILED: an infeasible QP was reported as solved'
    print '(A,I0)', 'infeasible QP: outcome ', outcome
    if (outcome /= sqpopt_direct_singular .and. outcome /= sqpopt_direct_max_changes) then
        error stop 'test_direct FAILED: unexpected outcome for an infeasible QP'
    end if

    call kkt%destroy()
    print '(A)', 'test_direct [the direct method''s special paths] PASSED'

    end subroutine test_direct_step

    subroutine setup(name, problem, x0, x_star)
    !! define one of the test problems
    character(len=*),                    intent(in)  :: name    !! `hs71` or `maratos`
    type(sqpopt_problem_type),           intent(out) :: problem !! the problem definition
    real(wp), dimension(:), allocatable, intent(out) :: x0      !! the starting point
    real(wp), dimension(:), allocatable, intent(out) :: x_star  !! the solution
    if (name == 'hs71') then
        call problem%set_problem_size(n=4, m=2)
        call problem%set_bounds(x_lb=spread(1.0_wp, 1, 4), x_ub=spread(5.0_wp, 1, 4), &
                                c_lb=[40.0_wp, 25.0_wp], c_ub=[40.0_wp, big])
        call problem%set_jacobian_sparsity(nnz=8, irow=[1,1,1,1,2,2,2,2], icol=[1,2,3,4,1,2,3,4])
        call problem%set_hessian_sparsity(nnz=10, irow=[1,2,3,4,2,3,4,3,4,4], icol=[1,1,1,1,2,2,2,3,3,4])
        call problem%set_functions(fc=fc_hs71, gjac=gjac_hs71, hess=hess_hs71)
        x0     = [1.0_wp, 5.0_wp, 5.0_wp, 1.0_wp]
        x_star = [1.0_wp, 4.7429994_wp, 3.8211500_wp, 1.3794083_wp]
    else
        call problem%set_problem_size(n=2, m=1)
        call problem%set_bounds(x_lb=[-big, -big], x_ub=[big, big], c_lb=[1.0_wp], c_ub=[1.0_wp])
        call problem%set_jacobian_sparsity(nnz=2, irow=[1,1], icol=[1,2])
        call problem%set_hessian_sparsity(nnz=2, irow=[1,2], icol=[1,2])
        call problem%set_functions(fc=fc_maratos, gjac=gjac_maratos, hess=hess_maratos)
        x0     = [cos(0.3_wp), sin(0.3_wp)]
        x_star = [1.0_wp, 0.0_wp]
    end if
    end subroutine setup

    subroutine run(name, cfg)
    !! solve a problem with one combination of options, twice, and check the results
    character(len=*),  intent(in) :: name !! the problem (see `setup`)
    type(config_type), intent(in) :: cfg  !! the options

    type(sqpopt_type)              :: solver
    type(sqpopt_problem_type)      :: problem
    type(sqpopt_options_type)      :: options
    type(sqpopt_trust_region_type) :: trust_region
    type(sqpopt_results_type)      :: r, r2
    real(wp), dimension(:), allocatable :: x0, x_star
    character(len=:), allocatable :: label
    integer :: istat

    label = name//', '//trim(cfg%name)
    call setup(name, problem, x0, x_star)

    options%max_iter             = 500
    options%hessian_mode         = cfg%hessian_mode
    options%inertia_control      = cfg%inertia
    options%direct_qp            = cfg%direct
    options%direct_least_squares = cfg%direct_ls
    options%factorization_threads = cfg%threads
    options%qp_solver_mode       = cfg%qp_mode
    trust_region%enabled         = cfg%trust_region

    call solver%initialize(problem=problem, options=options, trust_region=trust_region)
    call solver%solve(x0, istat)
    call solver%get_results(r)
    print '(A42,A,I2,A,I4,A,I3,A,I3,A,I4,A,ES9.2)', label, ': istat=', istat, ' iter=', r%iterations, ' direct=', &
        r%n_direct_qp, '/', r%n_qp_solves, ' fact=', r%n_factorizations, ' error=', maxval(abs(r%x - x_star))

    if (istat /= r%istat) error stop 'test_direct FAILED: '//label//': istat and results%istat differ'
    if (all(istat /= [sqpopt_success, sqpopt_acceptable, sqpopt_stalled])) then
        error stop 'test_direct FAILED: '//label//': did not converge'
    end if
    if (maxval(abs(r%x - x_star)) > 1.0e-4_wp) error stop 'test_direct FAILED: '//label//': wrong solution'
    if (cfg%direct .and. r%n_direct_qp == 0) error stop 'test_direct FAILED: '//label//': no QP was solved directly'
    if (r%n_direct_qp > r%n_qp_solves) error stop 'test_direct FAILED: '//label//': more direct QPs than QPs'
    if (.not. cfg%direct .and. r%n_direct_qp /= 0) error stop 'test_direct FAILED: '//label//': direct QPs without direct_qp'
    if ((cfg%direct .or. cfg%inertia) .and. r%n_factorizations == 0) then
        error stop 'test_direct FAILED: '//label//': nothing was factored'
    end if
    ! (the Maratos problem needs second-order corrections, which factor)
    if (cfg%direct_ls .and. name == 'maratos' .and. .not. cfg%trust_region .and. r%n_soc > 0 .and. &
        r%n_factorizations == 0) error stop 'test_direct FAILED: '//label//': the corrections did not factor'
    if (r%n_factorizations > 0 .and. .not. r%time_factorization >= 0.0_wp) then
        error stop 'test_direct FAILED: '//label//': no factorization time'
    end if

    ! no state carries over to a second solve:
    call solver%solve(x0, istat)
    call solver%get_results(r2)
    if (istat /= r%istat .or. r2%iterations /= r%iterations .or. r2%n_factorizations /= r%n_factorizations .or. &
        r2%n_direct_qp /= r%n_direct_qp .or. r2%n_eval_fc /= r%n_eval_fc .or. any(r2%x /= r%x)) then
        error stop 'test_direct FAILED: '//label//': the second solve differs'
    end if

    end subroutine run

    subroutine test_unavailable()
    !! without MUMPS, the options are invalid input
    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    real(wp), dimension(:), allocatable :: x0, x_star
    integer :: istat, k

    call setup('maratos', problem, x0, x_star)
    do k = 1, 2
        options = sqpopt_options_type()
        options%direct_qp            = k == 1
        options%direct_least_squares = k == 2
        call solver%initialize(problem=problem, options=options)
        call solver%solve(x0, istat)
        call solver%get_results(r)
        if (istat /= sqpopt_invalid_input .or. r%istat /= sqpopt_invalid_input) then
            error stop 'test_direct FAILED: an option was accepted without MUMPS'
        end if
        if (index(solver%status_message(), 'MUMPS') == 0) error stop 'test_direct FAILED: the message does not say why'
    end do
    print '(A)', 'test_direct [not built with MUMPS: options rejected] PASSED'

    end subroutine test_unavailable

    !------------------------------------------------------------------------
    ! Hock-Schittkowski problem 71
    !------------------------------------------------------------------------

    subroutine fc_hs71(x, f, c, status, data)
    !! the objective and the constraints
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f    = x(1)*x(4)*(x(1)+x(2)+x(3)) + x(3)
    c(1) = x(1)**2 + x(2)**2 + x(3)**2 + x(4)**2
    c(2) = x(1)*x(2)*x(3)*x(4)
    end subroutine fc_hs71

    subroutine gjac_hs71(x, g, jac_val, accuracy, status, data)
    !! the objective's gradient and the Jacobian's nonzeros
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` (in its sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! the accuracy asked for (the derivatives here are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    g = [x(4)*(2.0_wp*x(1)+x(2)+x(3)), x(1)*x(4), x(1)*x(4) + 1.0_wp, x(1)*(x(1)+x(2)+x(3))]
    jac_val = [2.0_wp*x(1), 2.0_wp*x(2), 2.0_wp*x(3), 2.0_wp*x(4), &
               x(2)*x(3)*x(4), x(1)*x(3)*x(4), x(1)*x(2)*x(4), x(1)*x(2)*x(3)]
    end subroutine gjac_hs71

    subroutine hess_hs71(x, lambda, hess_val, status, data)
    !! the Hessian of the Lagrangian, lower triangle
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    ! (1,1) (2,1) (3,1) (4,1) (2,2) (3,2) (4,2) (3,3) (4,3) (4,4):
    hess_val = [2.0_wp*x(4), x(4), x(4), 2.0_wp*x(1)+x(2)+x(3), 0.0_wp, 0.0_wp, x(1), 0.0_wp, x(1), 0.0_wp] &
             - lambda(1)*[2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp] &
             - lambda(2)*[0.0_wp, x(3)*x(4), x(2)*x(4), x(2)*x(3), 0.0_wp, x(1)*x(4), x(1)*x(3), 0.0_wp, x(1)*x(2), 0.0_wp]
    end subroutine hess_hs71

    !------------------------------------------------------------------------
    ! the Maratos example
    !------------------------------------------------------------------------

    subroutine fc_maratos(x, f, c, status, data)
    !! the objective and the constraint
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    f    = 2.0_wp*(x(1)**2 + x(2)**2 - 1.0_wp) - x(1)
    c(1) = x(1)**2 + x(2)**2
    end subroutine fc_maratos

    subroutine gjac_maratos(x, g, jac_val, accuracy, status, data)
    !! the objective's gradient and the Jacobian's nonzeros
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the Jacobian at `x` (in its sparsity pattern's order)
    integer,                intent(in)    :: accuracy !! the accuracy asked for (the derivatives here are analytic)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    g       = [4.0_wp*x(1) - 1.0_wp, 4.0_wp*x(2)]
    jac_val = [2.0_wp*x(1), 2.0_wp*x(2)]
    end subroutine gjac_maratos

    subroutine hess_maratos(x, lambda, hess_val, status, data)
    !! the (diagonal) Hessian of the Lagrangian
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
    real(wp), dimension(:), intent(out)   :: hess_val !! nonzero values of the Hessian of the Lagrangian at `x` (in its sparsity pattern's order)
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    hess_val = 4.0_wp - 2.0_wp*lambda(1)
    end subroutine hess_maratos

end program test_direct
