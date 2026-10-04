program test_kkt

    !! Unit tests of the sparse factorizations and what is built on them, with
    !! each sparse solver of the build (QDLDL always, and MUMPS in a build with
    !! the `HAS_MUMPS` preprocessor directive; see
    !! [[sqpopt_symmetric_solver_module]]):
    !!
    !! * the sparse symmetric solver: the inertia of small matrices, solves
    !!   (checked by their residuals), a singular matrix, and refactoring;
    !! * the KKT matrix ([[sqpopt_kkt_module]]) against a dense reference: its
    !!   solves and its count of negative curvature, for the exact Hessian and
    !!   for the BFGS and SR1 Hessians (their low-rank form), with variables at
    !!   bounds and rows outside the working set; and that a factorization is
    !!   reused;
    !! * the direct least-squares solver ([[sqpopt_least_squares_module]]): the
    !!   minimum-norm solution, also for dependent rows.
    !!
    !! It also tests the small dense routines those use
    !! (`dense_symmetric_inertia`, `dense_lu_factor`, and `dense_lu_solve`),
    !! and without MUMPS, that MUMPS reports itself unavailable.

    use sqpopt_symmetric_solver_module, only: sqpopt_symmetric_solver_type, sqpopt_has_mumps, sqpopt_has_lapack, &
                                              sqpopt_linear_solver_mumps, sqpopt_linear_solver_qdldl, &
                                              sqpopt_linear_solver_dense, sqpopt_linear_solver_lapack, &
                                              sqpopt_linear_solver_name
    use sqpopt_kkt_module,              only: sqpopt_kkt_type
    use sqpopt_least_squares_module,    only: sqpopt_least_squares_type
    use sqpopt_hessian_module,          only: sqpopt_hessian_type
    use sqpopt_dense_linalg_module,     only: dense_symmetric_inertia, dense_lu_factor, dense_lu_solve
    use sqpopt_types_module,            only: sqpopt_sparse_matrix
    use sqpopt_kinds,                   only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: tol = 1.0e-9_wp !! tolerance on the residuals and on the differences from the references
    integer :: k
    integer, dimension(4), parameter :: solvers = [sqpopt_linear_solver_mumps, sqpopt_linear_solver_qdldl, &
                                                   sqpopt_linear_solver_dense, sqpopt_linear_solver_lapack]

    write(*,*) '----------------------------'
    write(*,*) 'test_kkt'
    write(*,*) '----------------------------'

    call seed_rng()
    call test_dense()
    do k = 1, size(solvers)
        if (solvers(k) == sqpopt_linear_solver_mumps .and. .not. sqpopt_has_mumps) then
            call test_unavailable()
            cycle
        end if
        if (solvers(k) == sqpopt_linear_solver_lapack .and. .not. sqpopt_has_lapack) then
            call test_lapack_unavailable()
            cycle
        end if
        print '(2A)', 'sparse solver: ', sqpopt_linear_solver_name(solvers(k))
        call test_solver(solvers(k))
        call test_kkt_matrix('exact', solvers(k))
        call test_kkt_matrix('bfgs', solvers(k))
        call test_kkt_matrix('sr1', solvers(k))
        call test_least_squares(solvers(k))
    end do
    call test_dense_limit()

    print '(A)', 'test_kkt PASSED'

    contains

    subroutine seed_rng()
    !! seed the random number generator with a fixed sequence (so failures are reproducible)
    integer :: n, i
    integer, dimension(:), allocatable :: seed
    call random_seed(size=n)
    allocate(seed(n))
    seed = [(4321 + 104729*i, i=1,n)]
    call random_seed(put=seed)
    end subroutine seed_rng

    subroutine check_inertia(label, a, expected)
    !! check the inertia that `dense_symmetric_inertia` finds for a matrix
    character(len=*),         intent(in) :: label    !! name of the matrix, for the messages
    real(wp), dimension(:,:), intent(in) :: a        !! the symmetric matrix `dimension(n,n)`
    integer, dimension(3),    intent(in) :: expected !! its numbers of positive, negative, and zero eigenvalues
    integer :: n_pos, n_neg, n_zero
    call dense_symmetric_inertia(a, n_pos, n_neg, n_zero)
    print '(3A,3(I0,1X))', 'inertia of ', label, ' (positive, negative, zero): ', n_pos, n_neg, n_zero
    if (any([n_pos, n_neg, n_zero] /= expected)) error stop 'test_kkt FAILED: dense_symmetric_inertia: '//label
    end subroutine check_inertia

    subroutine test_dense()
    !! the inertia of `Q D Q^T` for a known diagonal `D` and a reflection `Q`, of matrices with
    !! zeros on the diagonal, and an LU solve
    integer, parameter :: n = 7
    real(wp), parameter :: d(n) = [3.0_wp, -2.0_wp, 0.5_wp, 0.0_wp, -1.0_wp, 4.0_wp, -0.25_wp]
    real(wp) :: q(n,n), a(n,n), v(n), b(n), x(n), lu(n,n)
    integer :: piv(n), i, n_pos, n_neg, n_zero
    logical :: ok

    call random_number(v)
    v = v/norm2(v)
    q = 0.0_wp
    do i = 1, n
        q(i,i) = 1.0_wp
        q(:,i) = q(:,i) - 2.0_wp*v*v(i)
    end do
    do i = 1, n
        a(:,i) = d(i)*q(:,i)
    end do
    a = matmul(a, transpose(q))
    call dense_symmetric_inertia(a, n_pos, n_neg, n_zero)
    print '(A,3(I0,1X))', 'inertia of Q*D*Q^T (positive, negative, zero): ', n_pos, n_neg, n_zero
    if (n_pos /= 3 .or. n_neg /= 3 .or. n_zero /= 1) error stop 'test_kkt FAILED: dense_symmetric_inertia'

    ! matrices with zeros on the diagonal, whose zero pivots are not zero eigenvalues:
    call check_inertia('[0 1; 1 0]', reshape([0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp], [2,2]), [1, 1, 0])
    call check_inertia('zero 2 by 2', reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], [2,2]), [0, 0, 2])
    call check_inertia('diag(0, 2, -1)', reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, &
                                                  0.0_wp, 0.0_wp, -1.0_wp], [3,3]), [1, 1, 1])
    call check_inertia('[0 1 0; 1 0 0; 0 0 -3]', reshape([0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
                                                          0.0_wp, 0.0_wp, -3.0_wp], [3,3]), [1, 2, 0])
    call check_inertia('[0 2 0; 2 0 3; 0 3 5]', reshape([0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 3.0_wp, &
                                                         0.0_wp, 3.0_wp, 5.0_wp], [3,3]), [2, 1, 0])
    ! (tridiagonal, with zeros on the diagonal and ones beside it: the eigenvalues are 2*cos(k*pi/6), k = 1..5)
    a(1:5,1:5) = 0.0_wp
    do i = 1, 4
        a(i,i+1) = 1.0_wp
        a(i+1,i) = 1.0_wp
    end do
    call check_inertia('tridiagonal (1, 0, 1)', a(1:5,1:5), [2, 2, 1])
    ! (the same eigenvalues, in a full matrix)
    a(1:5,1:5) = matmul(matmul(q(1:5,1:5), a(1:5,1:5)), transpose(q(1:5,1:5)))
    call check_inertia('tridiagonal (1, 0, 1), rotated', a(1:5,1:5), [2, 2, 1])

    ! (a nonsingular matrix: shift the zero eigenvalue)
    do i = 1, n
        a(i,i) = a(i,i) + 0.1_wp
    end do
    call random_number(b)
    lu = a
    call dense_lu_factor(lu, piv, ok)
    x = b
    call dense_lu_solve(lu, piv, x)
    if (.not. ok .or. maxval(abs(matmul(a, x) - b)) > tol) error stop 'test_kkt FAILED: dense LU solve'
    a(:,2) = a(:,1)
    call dense_lu_factor(a, piv, ok)
    if (ok) error stop 'test_kkt FAILED: singular matrix not detected by dense_lu_factor'

    print '(A)', 'test_kkt [dense routines] PASSED'
    end subroutine test_dense

    subroutine test_solver(which)
    !! the sparse solver on `[2 0 1; 0 -1 0; 1 0 0]` (one triangle given, with a duplicate entry),
    !! then on a singular matrix with the same pattern. (The third row's zero diagonal is given the
    !! sign `-1`, as QDLDL needs; MUMPS ignores it.)
    integer, intent(in) :: which !! the sparse solver (`sqpopt_linear_solver_*`)
    type(sqpopt_symmetric_solver_type) :: solver
    real(wp) :: b(3), x(3), y(3)
    logical :: ok
    integer :: i

    ! (the element (1,1) is given twice: 1.5 + 0.5)
    call solver%initialize(3, [1, 2, 3, 3, 1], [1, 2, 1, 3, 1], ok, solver=which, signs=[0, 0, -1])
    if (.not. (ok .and. solver%ready)) error stop 'test_kkt FAILED: the solver could not be started'
    call solver%factor([1.5_wp, -1.0_wp, 1.0_wp, 0.0_wp, 0.5_wp], ok)
    if (.not. ok) error stop 'test_kkt FAILED: factorization'
    if (solver%n_negative /= 2 .or. solver%n_null /= 0) error stop 'test_kkt FAILED: the solver''s inertia'
    b = [1.0_wp, 2.0_wp, 3.0_wp]
    x = b
    call solver%solve(x, ok)
    call solver%multiply(x, y)
    if (.not. ok .or. maxval(abs(y - b)) > tol) error stop 'test_kkt FAILED: solve'
    if (maxval(abs(x - [3.0_wp, -2.0_wp, -5.0_wp])) > tol) error stop 'test_kkt FAILED: wrong solution'

    ! new values, same pattern: the second variable has a zero row
    call solver%factor([1.5_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.5_wp], ok)
    if (.not. ok .or. solver%n_null /= 1 .or. solver%n_negative /= 1) error stop 'test_kkt FAILED: singular matrix'
    if (solver%n_factor /= 2) error stop 'test_kkt FAILED: factorization count'

    call solver%destroy()
    if (solver%ready .or. solver%n_factor /= 0) error stop 'test_kkt FAILED: destroy'
    x = b
    call solver%solve(x, ok)
    if (ok) error stop 'test_kkt FAILED: solve after destroy'

    ! the same system on two OpenMP threads, and with the number left to the environment
    do i = 0, 2, 2
        call solver%initialize(3, [1, 2, 3, 3, 1], [1, 2, 1, 3, 1], ok, threads=i, solver=which, signs=[0, 0, -1])
        if (ok) call solver%factor([1.5_wp, -1.0_wp, 1.0_wp, 0.0_wp, 0.5_wp], ok)
        x = b
        if (ok) call solver%solve(x, ok)
        if (.not. ok .or. maxval(abs(x - [3.0_wp, -2.0_wp, -5.0_wp])) > tol) error stop 'test_kkt FAILED: threads'
        call solver%destroy()
    end do

    print '(A)', 'test_kkt [sparse solver] PASSED'
    end subroutine test_solver

    subroutine test_kkt_matrix(mode, which)
    !! the KKT matrix of random problems against a dense reference, for the
    !! Hessian mode `mode` (`exact`, `bfgs`, or `sr1`), and several working sets
    character(len=*), intent(in) :: mode  !! the Hessian mode
    integer,          intent(in) :: which !! the sparse solver (`sqpopt_linear_solver_*`)
    integer, parameter :: n = 8, m = 3, n_trials = 20
    type(sqpopt_kkt_type)      :: kkt
    type(sqpopt_hessian_type)  :: h
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: a(n,n), hd(n,n), jd(m,n), kd(n+m,n+m), lu(n+m,n+m), s(n), b(n+m), x(n+m), xref(n+m)
    integer  :: status(m+n), piv(n+m), irow(n*(n+1)/2), icol(n*(n+1)/2)
    integer  :: i, j, k, trial, n_pos, n_neg, n_zero, n_factor0, n_checked
    logical  :: ok, nonsingular

    ! the sparsity patterns: a dense Jacobian, and the lower triangle of the Hessian
    jac%nrows = m
    jac%ncols = n
    jac%nnz   = m*n
    jac%irow  = [((i, j=1,n), i=1,m)]
    jac%icol  = [((j, j=1,n), i=1,m)]
    allocate(jac%val(m*n))
    k = 0
    do i = 1, n
        do j = 1, i
            k = k + 1
            irow(k) = i
            icol(k) = j
        end do
    end do
    if (mode == 'exact') then
        call kkt%initialize(n, m, jac%irow, jac%icol, ok, hess_irow=irow, hess_icol=icol, solver=which)
    else
        call kkt%initialize(n, m, jac%irow, jac%icol, ok, solver=which)
    end if
    if (.not. (ok .and. kkt%enabled)) error stop 'test_kkt FAILED: the KKT matrix could not be set up'

    n_checked = 0
    do trial = 1, n_trials

        ! a symmetric indefinite matrix, as the Hessian or as the source of the quasi-Newton pairs
        call random_number(a)
        a = a + transpose(a) - 1.0_wp
        select case (mode)
        case ('exact')
            call h%initialize(n, 1)
            call h%set_exact(irow, icol)
            call h%set_values([(a(irow(k), icol(k)), k=1, size(irow))], decay=.false.)
        case ('bfgs')
            a = matmul(a, a)   ! (positive definite, so every pair is accepted)
            do i = 1, n
                a(i,i) = a(i,i) + 0.5_wp
            end do
            call h%initialize(n, 4)
        case default
            call h%initialize(n, 4, use_sr1=.true.)
        end select
        if (mode /= 'exact') then
            do k = 1, 4
                call random_number(s)
                s = s - 0.5_wp
                if (mode == 'bfgs') then
                    call h%update_bfgs(s, matmul(a, s))
                else
                    call h%update_sr1(s, matmul(a, s))
                end if
            end do
            if (h%low_rank_size() == 0) error stop 'test_kkt FAILED: no quasi-Newton pairs were stored'
        end if
        h%shift = merge(0.3_wp, 0.0_wp, mod(trial, 2) == 0)
        call h%dense(hd)

        call random_number(jd)
        jd = jd - 0.5_wp
        jac%val = [((jd(i,j), j=1,n), i=1,m)]

        ! a random working set: each row and each bound in it with probability 1/2 and 1/4
        call random_number(b)
        status = 0
        where (b(1:m) < 0.5_wp) status(1:m) = 1
        call random_number(b(1:n))
        where (b(1:n) < 0.25_wp) status(m+1:) = -1

        ! the dense reference (see the module documentation of sqpopt_kkt_module)
        kd = 0.0_wp
        do j = 1, n
            if (status(m+j) /= 0) then
                kd(j,j) = 1.0_wp
            else
                do i = 1, n
                    if (status(m+i) == 0) kd(i,j) = hd(i,j)
                end do
            end if
        end do
        do i = 1, m
            if (status(i) == 0) then
                kd(n+i,n+i) = -1.0_wp
            else
                do j = 1, n
                    if (status(m+j) == 0) then
                        kd(n+i,j) = jd(i,j)
                        kd(j,n+i) = jd(i,j)
                    end if
                end do
            end if
        end do

        call kkt%new_matrices()
        call kkt%factor(h, jac, status, ok)
        if (.not. ok) error stop 'test_kkt FAILED: KKT factorization'

        ! the count of negative curvature, from the reference's inertia
        call dense_symmetric_inertia(kd, n_pos, n_neg, n_zero)
        lu = kd
        call dense_lu_factor(lu, piv, nonsingular)
        if (n_zero > 0 .or. .not. nonsingular .or. kkt%singular) cycle   ! (a singular face: nothing to compare)
        n_checked = n_checked + 1
        if (kkt%n_negative /= max(n_neg - m, 0)) then
            print '(3A,I0,2(A,I0))', 'mode ', mode, ', trial ', trial, ': n_negative = ', kkt%n_negative, &
                ', reference ', max(n_neg - m, 0)
            error stop 'test_kkt FAILED: wrong count of negative curvature'
        end if

        ! a solve
        call random_number(b)
        xref = b
        call dense_lu_solve(lu, piv, xref)
        x = b
        call kkt%solve(x, ok, hessian=h)
        if (.not. ok .or. maxval(abs(x - xref)) > tol*max(1.0_wp, maxval(abs(xref)))) then
            print '(3A,I0,A,ES10.2)', 'mode ', mode, ', trial ', trial, ': error ', maxval(abs(x - xref))
            error stop 'test_kkt FAILED: wrong KKT solve'
        end if

        ! the same working set again: the factorization is reused
        n_factor0 = kkt%solver%n_factor
        call kkt%factor(h, jac, status, ok)
        if (.not. ok .or. kkt%solver%n_factor /= n_factor0) error stop 'test_kkt FAILED: factorization not reused'

    end do
    if (n_checked < n_trials/2) error stop 'test_kkt FAILED: too few nonsingular trials'

    call kkt%destroy()
    if (kkt%enabled) error stop 'test_kkt FAILED: KKT destroy'
    print '(3A,I0,A)', 'test_kkt [KKT matrix, ', mode, ' Hessian: ', n_checked, ' working sets] PASSED'
    end subroutine test_kkt_matrix

    subroutine test_least_squares(which)
    !! the minimum-norm solution of `J_S d = r`, for independent rows, for a duplicated row, and
    !! for a Jacobian with small elements in the rows in use and large ones elsewhere, and the
    !! least-squares multipliers of that Jacobian
    integer, intent(in) :: which !! the sparse solver (`sqpopt_linear_solver_*`)
    integer, parameter :: n = 5, m = 3
    real(wp), parameter :: small = 1.0e-6_wp              !! scale of the small Jacobian
    real(wp), parameter :: large = 1.0e6_wp               !! scale, relative to it, of its row and column not in use
    real(wp), parameter :: lambda0(2) = [2.0_wp, -3.0_wp] !! the multipliers to recover
    type(sqpopt_least_squares_type) :: ls
    type(sqpopt_sparse_matrix) :: jac
    real(wp) :: jd(m,n), r(m), d(n), jjt(2,2), y(2), dref(n), lambda(m)
    integer  :: piv(2), i, j
    logical  :: ok, nonsingular

    call random_number(jd)
    jd = jd - 0.5_wp
    jd(3,:) = jd(1,:)   ! (row 3 duplicates row 1)
    jac%nrows = m
    jac%ncols = n
    jac%nnz   = m*n
    jac%irow  = [((i, j=1,n), i=1,m)]
    jac%icol  = [((j, j=1,n), i=1,m)]
    jac%val   = [((jd(i,j), j=1,n), i=1,m)]
    r = [1.0_wp, -2.0_wp, 1.0_wp]

    call ls%initialize(n, m, jac%irow, jac%icol, ok, solver=which)
    if (.not. (ok .and. ls%enabled)) error stop 'test_kkt FAILED: the least-squares solver could not be started'

    ! rows 1 and 2: d = J^T (J J^T)^{-1} r
    jjt = matmul(jd(1:2,:), transpose(jd(1:2,:)))
    call dense_lu_factor(jjt, piv, nonsingular)
    y = r(1:2)
    call dense_lu_solve(jjt, piv, y)
    dref = matmul(transpose(jd(1:2,:)), y)
    call ls%min_norm(jac, [.true., .true., .false.], r, d, ok)
    print '(A,ES10.2)', 'least squares, independent rows: error ', maxval(abs(d - dref))
    if (.not. ok .or. maxval(abs(d - dref)) > 1.0e-6_wp) error stop 'test_kkt FAILED: minimum-norm solution'

    ! all three rows (two the same, with the same right-hand side): the same solution
    call ls%min_norm(jac, [.true., .true., .true.], r, d, ok)
    print '(A,ES10.2)', 'least squares, a duplicated row: error ', maxval(abs(d - dref))
    if (.not. ok .or. maxval(abs(d - dref)) > 1.0e-6_wp) error stop 'test_kkt FAILED: dependent rows'

    ! a Jacobian with small elements (and the right-hand side scaled with it), and a
    ! large row that isn't used: the same solution, since the regularization is
    ! relative to the rows in use
    jac%val = small*jac%val
    jac%val(2*n+1:3*n) = large*jac%val(2*n+1:3*n)
    call ls%new_matrices()
    call ls%min_norm(jac, [.true., .true., .false.], small*r, d, ok)
    print '(A,ES10.2)', 'least squares, a small Jacobian: error ', maxval(abs(d - dref))
    if (.not. ok .or. maxval(abs(d - dref)) > 1.0e-6_wp) error stop 'test_kkt FAILED: minimum-norm solution, small Jacobian'

    ! and its multipliers: for g = J^T lambda0, they are lambda0 (the last variable
    ! is not free, and its column is large too)
    jac%val([n, 2*n]) = large*jac%val([n, 2*n])
    call ls%new_matrices()
    lambda = 0.0_wp
    call ls%multipliers(jac, [.true., .true., .false.], [(i < n, i=1,n)], &
                        small*matmul(transpose(jd(1:2,:)), lambda0), lambda, ok)
    print '(A,ES10.2)', 'least-squares multipliers, a small Jacobian: error ', maxval(abs(lambda(1:2) - lambda0))
    if (.not. ok .or. maxval(abs(lambda(1:2) - lambda0)) > 1.0e-6_wp) &
        error stop 'test_kkt FAILED: least-squares multipliers, small Jacobian'

    call ls%destroy()
    call ls%min_norm(jac, [.true., .true., .false.], r, d, ok)
    if (ok .or. ls%enabled) error stop 'test_kkt FAILED: least squares after destroy'

    print '(A)', 'test_kkt [least squares] PASSED'
    end subroutine test_least_squares

    subroutine test_dense_limit()
    !! the dense solver refuses a matrix above its largest order (the options then
    !! fall back on the matrix-free methods)
    type(sqpopt_symmetric_solver_type) :: solver
    integer, dimension(:), allocatable :: idx
    logical :: ok
    integer :: i
    idx = [(i, i = 1, 2001)]
    call solver%initialize(2001, idx, idx, ok, solver=sqpopt_linear_solver_dense)
    if (ok .or. solver%ready) error stop 'test_kkt FAILED: the dense solver accepted an order above its limit'
    call solver%initialize(2000, idx(1:2000), idx(1:2000), ok, solver=sqpopt_linear_solver_dense)
    if (.not. ok) error stop 'test_kkt FAILED: the dense solver refused an order within its limit'
    call solver%destroy()
    print '(A)', 'test_kkt [dense solver: order limit] PASSED'
    end subroutine test_dense_limit

    subroutine test_lapack_unavailable()
    !! without LAPACK, the sparse solver can't be started with it
    type(sqpopt_symmetric_solver_type) :: solver
    logical :: ok
    call solver%initialize(2, [1, 2], [1, 2], ok, solver=sqpopt_linear_solver_lapack)
    if (ok .or. solver%ready) error stop 'test_kkt FAILED: the LAPACK solver started without LAPACK'
    print '(A)', 'test_kkt [not built with LAPACK: unavailable] PASSED'
    end subroutine test_lapack_unavailable

    subroutine test_unavailable()
    !! without MUMPS, nothing can be started with it
    type(sqpopt_symmetric_solver_type) :: solver
    type(sqpopt_kkt_type) :: kkt
    type(sqpopt_least_squares_type) :: ls
    logical :: ok
    call solver%initialize(2, [1, 2], [1, 2], ok, solver=sqpopt_linear_solver_mumps)
    if (ok .or. solver%ready) error stop 'test_kkt FAILED: the solver started without MUMPS'
    call kkt%initialize(2, 1, [1], [1], ok, solver=sqpopt_linear_solver_mumps)
    if (ok .or. kkt%enabled) error stop 'test_kkt FAILED: the KKT matrix was set up without MUMPS'
    call ls%initialize(2, 1, [1], [1], ok, solver=sqpopt_linear_solver_mumps)
    if (ok .or. ls%enabled) error stop 'test_kkt FAILED: the least-squares solver started without MUMPS'
    print '(A)', 'test_kkt [not built with MUMPS: unavailable] PASSED'
    end subroutine test_unavailable

end program test_kkt
