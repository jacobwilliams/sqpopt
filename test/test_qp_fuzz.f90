program test_qp_fuzz

    !! Randomized test of the two active-set QP solvers
    !! (`sqpopt_dense_qp_type` and `sqpopt_reduced_hessian_qp_type`) on many
    !! small random QPs
    !!
    !!   minimize   0.5 p^T B p + g^T p
    !!   subject to c_lb <= J p <= c_ub,  x_lb <= p <= x_ub
    !!
    !! (with `x=0` and `c=0`, so the rows are just `J p`), including
    !! degenerate cases: dependent (duplicated) rows, rows parallel to a
    !! variable bound, equality rows, fixed variables, and more active
    !! constraints than variables. Each solution is checked against the
    !! first-order KKT conditions -- primal feasibility, stationarity with
    !! correctly-signed bound multipliers, and sign/complementarity of the
    !! returned constraint multipliers -- which, for a convex QP, certify
    !! optimality. Also:
    !!
    !! * nonconvex QPs (an indefinite SR1 `B`, with every variable boxed so
    !!   the QP is bounded) must still return a KKT point;
    !! * QPs with two contradictory rows must be reported as
    !!   `sqpopt_infeasible`.
    !!
    !! The random sequence is fixed (seeded), so failures are reproducible.

    use sqpopt_hessian_module,            only: sqpopt_hessian_type
    use sqpopt_qp_dense_module,           only: sqpopt_dense_qp_type
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_reduced_hessian_qp_type
    use sqpopt_types_module,              only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_infeasible, sqpopt_infinity
    use sqpopt_kinds,                     only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: n_trials = 400
    integer, parameter :: solver_dense = 1, solver_rh = 2

    integer :: trial, solver, kind, n_fail(2), n_run(2)
    logical :: trial_ok                          !! result of the current trial (set by `fail`)
    character(len=:), allocatable :: trial_why   !! why the current trial failed
    character(len=*), parameter :: solver_name(2) = ['dense          ', 'reduced-Hessian']

    write(*,*) '----------------------------'
    write(*,*) 'test_qp_fuzz'
    write(*,*) '----------------------------'

    call seed_rng()

    n_fail = 0
    n_run  = 0
    do trial = 1, n_trials
        ! trial kinds: 1-6 convex feasible (various degeneracies),
        ! 7 nonconvex (SR1) feasible, 8 infeasible
        kind = 1 + mod(trial-1, 8)
        do solver = solver_dense, solver_rh
            n_run(solver) = n_run(solver) + 1
            if (.not. run_trial(trial, kind, solver)) n_fail(solver) = n_fail(solver) + 1
        end do
    end do

    do solver = solver_dense, solver_rh
        print '(A,A,A,I0,A,I0)', 'solver ', trim(solver_name(solver)), ': failures = ', n_fail(solver), ' / ', n_run(solver)
    end do
    if (any(n_fail > 0)) error stop 'test_qp_fuzz FAILED'
    print '(A)', 'test_qp_fuzz PASSED'

    contains

    subroutine seed_rng()
    integer :: n, i
    integer, dimension(:), allocatable :: seed
    call random_seed(size=n)
    allocate(seed(n))
    seed = [(12345 + 7919*i, i=1,n)]
    call random_seed(put=seed)
    end subroutine seed_rng

    real(wp) function urand(a, b)
    !! uniform random number in `[a,b)`
    real(wp), intent(in) :: a, b
    real(wp) :: r
    call random_number(r)
    urand = a + (b-a)*r
    end function urand

    integer function irand(a, b)
    !! uniform random integer in `[a,b]`
    integer, intent(in) :: a, b
    irand = min(b, a + int(urand(0.0_wp, 1.0_wp)*(b-a+1)))
    end function irand

    logical function run_trial(trial, kind, solver) result(ok)
    integer, intent(in) :: trial, kind, solver

    integer :: n, m, i, j, k
    real(wp), dimension(:,:), allocatable :: jd, bd, a
    real(wp), dimension(:), allocatable :: g, x_lb, x_ub, c_lb, c_ub, p, lambda, pf, jp, r, s, zero_n, zero_m
    type(sqpopt_hessian_type) :: hess
    type(sqpopt_sparse_matrix) :: jac
    type(sqpopt_dense_qp_type) :: dense_qp
    type(sqpopt_reduced_hessian_qp_type) :: rh_qp
    integer :: istat
    real(wp) :: tol, scale
    logical :: nonconvex, infeasible

    nonconvex  = kind == 7
    infeasible = kind == 8

    n = irand(1, 6)
    m = irand(0, 6)
    if (infeasible) m = max(m, 2)

    ! ---- the Hessian: a limited-memory approximation built from pairs (s, A s) ----
    allocate(a(n,n), s(n))
    do i = 1, n
        do j = 1, n
            a(i,j) = urand(-1.0_wp, 1.0_wp)
        end do
    end do
    if (nonconvex) then
        a = 0.5_wp*(a + transpose(a))                       ! symmetric, generally indefinite
        call hess%initialize(n, 5, use_sr1=.true.)
    else
        a = matmul(transpose(a), a)
        do i = 1, n
            a(i,i) = a(i,i) + 0.1_wp                        ! symmetric positive definite
        end do
        call hess%initialize(n, 5)
    end if
    do k = 1, irand(0, 5)
        do i = 1, n
            s(i) = urand(-1.0_wp, 1.0_wp)
        end do
        if (nonconvex) then
            call hess%update_sr1(s, matmul(a, s))
        else
            call hess%update_bfgs(s, matmul(a, s))
        end if
    end do
    allocate(bd(n,n))
    block
        real(wp), dimension(n) :: e, be
        do j = 1, n
            e = 0.0_wp; e(j) = 1.0_wp
            call hess%hv_product(e, be)
            bd(:,j) = be
        end do
    end block

    allocate(g(n))
    do i = 1, n
        g(i) = urand(-5.0_wp, 5.0_wp)
    end do

    ! ---- a random feasible point pf, and bounds around it ----
    allocate(pf(n), x_lb(n), x_ub(n))
    do i = 1, n
        pf(i) = urand(-1.0_wp, 1.0_wp)
        x_lb(i) = pf(i) - urand(0.0_wp, 2.0_wp)
        x_ub(i) = pf(i) + urand(0.0_wp, 2.0_wp)
        if (.not. nonconvex) then  ! a nonconvex QP needs a bounded box
            if (urand(0.0_wp, 1.0_wp) < 0.3_wp) x_lb(i) = -sqpopt_infinity
            if (urand(0.0_wp, 1.0_wp) < 0.3_wp) x_ub(i) =  sqpopt_infinity
        end if
        if (kind == 3 .and. i == 1) then  ! a fixed variable
            x_lb(i) = pf(i); x_ub(i) = pf(i)
        end if
    end do
    ! the solvers assume x (here 0) is within the variable bounds:
    x_lb = min(x_lb, 0.0_wp)
    x_ub = max(x_ub, 0.0_wp)
    pf = min(max(pf, x_lb), x_ub)

    ! ---- random constraint rows (with degenerate variations) ----
    allocate(jd(m,n))
    do i = 1, m
        do j = 1, n
            jd(i,j) = urand(-2.0_wp, 2.0_wp)
            if (urand(0.0_wp, 1.0_wp) < 0.3_wp) jd(i,j) = 0.0_wp
        end do
    end do
    if (kind == 2 .and. m >= 2) jd(2,:) = 2.0_wp*jd(1,:)                ! dependent (parallel) rows
    if (kind == 4 .and. m >= 1) then                                   ! a row parallel to a variable bound
        jd(1,:) = 0.0_wp; jd(1,1) = 1.0_wp
    end if
    if (kind == 5 .and. m >= 1) then                                   ! rows duplicating one another
        do i = 2, m
            jd(i,:) = jd(1,:)
        end do
    end if
    allocate(jp(m), c_lb(m), c_ub(m))
    jp = matmul(jd, pf)
    do i = 1, m
        c_lb(i) = jp(i) - urand(0.0_wp, 1.0_wp)
        c_ub(i) = jp(i) + urand(0.0_wp, 1.0_wp)
        if (urand(0.0_wp, 1.0_wp) < 0.25_wp .or. kind == 6) then       ! equality rows (kind 6: all of them,
            c_lb(i) = jp(i); c_ub(i) = jp(i)                            ! often more than n)
        else if (urand(0.0_wp, 1.0_wp) < 0.2_wp) then
            c_lb(i) = -sqpopt_infinity
        else if (urand(0.0_wp, 1.0_wp) < 0.2_wp) then
            c_ub(i) = sqpopt_infinity
        end if
    end do
    if (kind == 5 .and. m >= 1) then                                    ! duplicates: same bounds (consistent)
        c_lb = c_lb(1); c_ub = c_ub(1)
    end if
    if (infeasible) then                                               ! row 2 = row 1, with disjoint bounds
        jd(1,1) = 1.0_wp   ! (make sure row 1 is nonzero)
        jd(2,:) = jd(1,:)
        jp = matmul(jd, pf)
        c_lb(1) = jp(1) - 1.0_wp; c_ub(1) = jp(1)
        c_lb(2) = jp(1) + 0.5_wp; c_ub(2) = jp(1) + 1.0_wp
    end if

    ! ---- sparse Jacobian ----
    jac%nrows = m
    jac%ncols = n
    jac%nnz   = count(jd /= 0.0_wp)
    allocate(jac%irow(jac%nnz), jac%icol(jac%nnz), jac%val(jac%nnz))
    k = 0
    do i = 1, m
        do j = 1, n
            if (jd(i,j) /= 0.0_wp) then
                k = k + 1
                jac%irow(k) = i; jac%icol(k) = j; jac%val(k) = jd(i,j)
            end if
        end do
    end do

    allocate(p(n), lambda(m), zero_n(n), zero_m(m))
    zero_n = 0.0_wp
    zero_m = 0.0_wp
    if (solver == solver_dense) then
        call dense_qp%solve(hess, jac, zero_n, g, zero_m, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)
        tol = 1.0e-6_wp
    else
        rh_qp%max_pcg_iter = 4*n
        call rh_qp%solve(hess, jac, zero_n, g, zero_m, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)
        tol = 1.0e-5_wp
    end if

    trial_ok  = .true.
    trial_why = ''
    if (infeasible) then
        if (istat /= sqpopt_infeasible) call fail('infeasible QP not reported as sqpopt_infeasible')
    else if (istat /= sqpopt_success) then
        call fail('feasible QP not solved')
    else
        ! ---- KKT checks ----
        scale = 1.0_wp + maxval(abs(g)) + maxval(abs(matmul(bd, p)))
        jp = matmul(jd, p)
        if (any(p < x_lb - tol) .or. any(p > x_ub + tol))  call fail('variable bound violated')
        if (any(jp < c_lb - tol) .or. any(jp > c_ub + tol)) call fail('constraint violated')
        allocate(r(n))
        r = matmul(bd, p) + g - matmul(transpose(jd), lambda)
        do j = 1, n
            if (x_ub(j) - x_lb(j) <= tol) cycle                          ! fixed: free multiplier
            if (p(j) - x_lb(j) <= tol) then
                if (r(j) < -tol*scale) call fail('wrong-sign bound multiplier (lower)')
            else if (x_ub(j) - p(j) <= tol) then
                if (r(j) > tol*scale) call fail('wrong-sign bound multiplier (upper)')
            else
                if (abs(r(j)) > tol*scale) call fail('not stationary')
            end if
        end do
        do i = 1, m
            if (c_ub(i) - c_lb(i) <= tol) cycle                          ! equality: free multiplier
            if (lambda(i) > tol*scale .and. abs(jp(i)-c_lb(i)) > tol) call fail('lambda>0 off the lower bound')
            if (lambda(i) < -tol*scale .and. abs(jp(i)-c_ub(i)) > tol) call fail('lambda<0 off the upper bound')
        end do
    end if

    ok = trial_ok
    if (.not. ok) print '(A,I0,A,I0,A,A,A,I0,A,I0,A,I0,2A)', 'trial ', trial, ' (kind ', kind, ', ', &
        trim(solver_name(solver)), ', n=', n, ', m=', m, '): istat=', istat, ': ', trial_why

    end function run_trial

    subroutine fail(msg)
    !! record that the current trial failed (the first reason is kept)
    character(len=*), intent(in) :: msg
    if (trial_ok) trial_why = msg
    trial_ok = .false.
    end subroutine fail

end program test_qp_fuzz
