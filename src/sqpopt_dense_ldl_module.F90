!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The dense backends of [[sqpopt_symmetric_solver_module]], for small
!  matrices: the matrix is stored dense (order `n`, so \( n^2 \) elements),
!  and factored by the Bunch-Kaufman \( LDL^T \) factorization (symmetric
!  pivoting, with 1 by 1 and 2 by 2 pivots; about \( n^3/3 \) operations).
!  That one factorization gives both the solves and the exact inertia: the
!  signs of the eigenvalues of `D`'s blocks (Sylvester's law). So, unlike
!  QDLDL, it always knows the inertia, also of a Hessian with zeros on its
!  diagonal, and, unlike MUMPS, it has almost no overhead per call. But
!  each factorization costs \( O(n^3) \): it is meant for small matrices,
!  and refuses to start above `dense_max_order` (the factorization options
!  then fall back on the matrix-free methods). There are two:
!
!  * [[sqpopt_dense_ldl_type]] (`options%linear_solver =
!    sqpopt_linear_solver_dense`): SQPOPT's own factorization,
!    [[dense_ldl_factor]] and [[dense_ldl_solve]]. Always available.
!  * [[sqpopt_lapack_ldl_type]] (`options%linear_solver =
!    sqpopt_linear_solver_lapack`): the same factorization by LAPACK's
!    `DSYTRF` and `DSYTRS`, whose blocked, BLAS-level code is faster on the
!    larger matrices. It needs a library built with LAPACK and BLAS (the
!    `HAS_LAPACK` preprocessor directive, and `-llapack -lblas`), in double
!    precision; without it, `sqpopt_has_lapack` is false and the option
!    can't be chosen. This module is the only one that refers to LAPACK.
!
!  Both store the factors in the same form (LAPACK's), so they differ only
!  in the two procedures that factor and solve; the inertia, and the
!  handling of a singular matrix, are shared ([[count_inertia]]). A
!  singular matrix is not an error: its zero eigenvalues (those of `D`'s
!  blocks below `dense_zero_tol` times the largest element) are counted
!  (`n_null`), and replaced in `D` by a large value, which makes that
!  component of a solution negligible, so that a solve returns one of the
!  solutions of a consistent system (as MUMPS's null-pivot detection does).

    module sqpopt_dense_ldl_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_sparse_ldl_module, only: sqpopt_sparse_ldl_type
    use sqpopt_dense_linalg_module, only: dense_ldl_factor, dense_ldl_solve
    use sqpopt_types_module, only: sqpopt_all_finite

    implicit none

    private

#ifdef HAS_LAPACK
#if defined(REAL32) || defined(REAL128)
#error "HAS_LAPACK needs the default (double precision) real kind: don't define REAL32 or REAL128 with it"
#endif
    interface
        subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
            !! LAPACK: the Bunch-Kaufman factorization of a symmetric matrix
            import :: wp
            implicit none
            character, intent(in) :: uplo         !! which triangle (`'L'`)
            integer,   intent(in) :: n            !! order
            integer,   intent(in) :: lda          !! leading dimension of `a`
            real(wp),  intent(inout) :: a(lda,*)  !! the matrix, overwritten by its factors
            integer,   intent(out) :: ipiv(*)     !! the pivots (negative for a 2 by 2 block)
            real(wp),  intent(inout) :: work(*)   !! workspace
            integer,   intent(in) :: lwork        !! its size (`-1`: a query)
            integer,   intent(out) :: info        !! `0`, or `i > 0` if `D(i,i)` is exactly zero
        end subroutine dsytrf
        subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
            !! LAPACK: solves with the factors of `DSYTRF`
            import :: wp
            implicit none
            character, intent(in) :: uplo         !! which triangle (`'L'`)
            integer,   intent(in) :: n            !! order
            integer,   intent(in) :: nrhs         !! number of right-hand sides
            integer,   intent(in) :: lda          !! leading dimension of `a`
            real(wp),  intent(in) :: a(lda,*)     !! the factors
            integer,   intent(in) :: ipiv(*)      !! the pivots
            integer,   intent(in) :: ldb          !! leading dimension of `b`
            real(wp),  intent(inout) :: b(ldb,*)  !! the right-hand sides, overwritten by the solutions
            integer,   intent(out) :: info        !! `0` on success
        end subroutine dsytrs
    end interface
    logical, parameter, public :: sqpopt_has_lapack = .true.  !! whether the library was built with LAPACK (the
                                                              !! `HAS_LAPACK` preprocessor directive), which
                                                              !! `sqpopt_linear_solver_lapack` requires
#else
    logical, parameter, public :: sqpopt_has_lapack = .false. !! whether the library was built with LAPACK (the
                                                              !! `HAS_LAPACK` preprocessor directive), which
                                                              !! `sqpopt_linear_solver_lapack` requires
#endif

    real(wp), parameter :: dense_zero_tol = 1.0e-5_wp*epsilon(1.0_wp) !! an eigenvalue at most this times the largest
                                                                     !! element is a zero eigenvalue (the null-pivot
                                                                     !! threshold of QDLDL and MUMPS: well below the
                                                                     !! constraint block's eigenvalues of a KKT matrix
                                                                     !! with a large Hessian shift)
    integer, parameter, public :: dense_max_order = 2000 !! the largest order the dense backend accepts (its two
                                                         !! matrices then take 64 MB in double precision)

    type, extends(sqpopt_sparse_ldl_type), public :: sqpopt_dense_ldl_type
        !! the dense backend, with SQPOPT's own factorization (see the module
        !! documentation)
        private
        integer :: n = 0                                 !! order of the matrix
        integer, dimension(:), allocatable :: irow       !! row indices of the pattern's entries
        integer, dimension(:), allocatable :: icol       !! column indices of the pattern's entries
        real(wp), dimension(:,:), allocatable :: a       !! the matrix last factored `dimension(n,n)`
        real(wp), dimension(:,:), allocatable :: ldl     !! its factors `L` and `D` (as LAPACK's `DSYTRF` stores
                                                         !! them, in the lower triangle) `dimension(n,n)`
        integer,  dimension(:),   allocatable :: piv     !! the pivots of the factorization `dimension(n)`
        contains
        procedure :: start      => dense_start
        procedure :: refactor   => dense_refactor
        procedure :: back_solve => dense_back_solve
        procedure :: multiply   => dense_multiply
        procedure :: free       => dense_free
        procedure, private :: factor        => dense_factor        !! factor `ldl` in place
        procedure, private :: solve_factors => dense_solve_factors !! one solve with the factors
    end type sqpopt_dense_ldl_type

    type, extends(sqpopt_dense_ldl_type), public :: sqpopt_lapack_ldl_type
        !! the dense backend with LAPACK's `DSYTRF`/`DSYTRS` (see the module
        !! documentation); it can only be started in a library built with
        !! `HAS_LAPACK`
        private
        real(wp), dimension(:), allocatable :: work !! `DSYTRF`'s workspace
        contains
        procedure :: start         => lapack_start
        procedure, private :: factor        => lapack_factor
        procedure, private :: solve_factors => lapack_solve_factors
    end type sqpopt_lapack_ldl_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  keep the pattern, and allocate the dense matrices. `ok` is false if the
!  order exceeds `dense_max_order`, or the matrices can't be allocated.
!  `threads` and `signs` are not used (this backend pivots, and is
!  single-threaded).

    subroutine dense_start(me, n, irow, icol, ok, threads, signs)

    class(sqpopt_dense_ldl_type), intent(inout) :: me
    integer,               intent(in)  :: n       !! order of the matrix
    integer, dimension(:), intent(in)  :: irow    !! row indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)  :: icol    !! column indices of the entries `dimension(nnz)`
    logical,               intent(out) :: ok      !! whether the backend is ready
    integer, optional,     intent(in)  :: threads !! (not used)
    integer, dimension(:), optional, intent(in) :: signs !! (not used)

    integer :: alloc_stat

    ok = .false.
    if (present(threads) .or. present(signs)) continue   ! (not used)
    if (n > dense_max_order) return
    allocate(me%a(n,n), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%ldl(n,n), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%piv(n), stat=alloc_stat)
    if (alloc_stat /= 0) then
        me%out_of_memory = .true.
        call me%free()
        return
    end if
    me%irow = irow
    me%icol = icol
    me%n = n
    ok = .true.

    end subroutine dense_start
!*******************************************************************************


!*******************************************************************************
!>
!  form the dense matrix from the values `val`, factor it, and count its
!  inertia from the factors (see the module documentation).

    subroutine dense_refactor(me, val, ok)

    class(sqpopt_dense_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: val !! values of the entries `dimension(nnz)`
    logical,                intent(out) :: ok  !! whether the factorization succeeded

    integer :: k, i, j
    real(wp) :: amax

    ok = .false.
    me%inertia_known = .true.
    if (.not. sqpopt_all_finite(val)) return
    do j = 1, me%n
        do i = 1, me%n
            me%a(i,j) = 0.0_wp
        end do
    end do
    do k = 1, size(val)
        i = me%irow(k)
        j = me%icol(k)
        me%a(i,j) = me%a(i,j) + val(k)
        if (i /= j) me%a(j,i) = me%a(j,i) + val(k)
    end do
    amax = 0.0_wp
    do j = 1, me%n
        do i = 1, me%n
            me%ldl(i,j) = me%a(i,j)
            amax = max(amax, abs(me%a(i,j)))
        end do
    end do
    call me%factor(ok)
    if (.not. ok) return
    call count_inertia(me, amax)

    end subroutine dense_refactor
!*******************************************************************************

!*******************************************************************************
!>
!  the inertia of the matrix, from the eigenvalues of the blocks of the
!  factor `D` (Sylvester's law). An eigenvalue at most `dense_zero_tol`
!  times the matrix's largest element `amax` is a zero eigenvalue: it is
!  counted in `n_null`, and replaced in `D` by a large positive value, so
!  that the solves never divide by it (see the module documentation).

    subroutine count_inertia(me, amax)

    class(sqpopt_dense_ldl_type), intent(inout) :: me
    real(wp), intent(in) :: amax !! the largest element of the matrix, in magnitude

    integer :: k
    real(wp) :: tol, big, p, q, r, mean, radius, lam1, lam2, c, s, h
    logical :: null1, null2

    tol = dense_zero_tol*max(amax, tiny(1.0_wp))
    big = max(amax, 1.0_wp)/sqrt(epsilon(1.0_wp))
    me%n_negative = 0
    me%n_null     = 0
    k = 1
    do while (k <= me%n)
        if (me%piv(k) > 0) then
            ! a 1 by 1 block:
            p = me%ldl(k,k)
            if (abs(p) <= tol) then
                me%n_null = me%n_null + 1
                me%ldl(k,k) = big
            else if (p < 0.0_wp) then
                me%n_negative = me%n_negative + 1
            end if
            k = k + 1
        else
            ! a 2 by 2 block [p, r; r, q] (r is nonzero: the pivoting chose it as the largest element of its column):
            p = me%ldl(k,k)
            r = me%ldl(k+1,k)
            q = me%ldl(k+1,k+1)
            mean   = 0.5_wp*(p + q)
            radius = hypot(0.5_wp*(p - q), r)
            lam1 = mean + sign(radius, mean)   ! (the eigenvalue of the larger magnitude)
            lam2 = (p*q - r*r)/lam1            ! (the other, without cancellation)
            null1 = abs(lam1) <= tol
            null2 = abs(lam2) <= tol
            if (null1) me%n_null = me%n_null + 1
            if (null2) me%n_null = me%n_null + 1
            if (.not. null1 .and. lam1 < 0.0_wp) me%n_negative = me%n_negative + 1
            if (.not. null2 .and. lam2 < 0.0_wp) me%n_negative = me%n_negative + 1
            if (null1 .or. null2) then
                ! rebuild the block from its eigenvectors, with large eigenvalues in place
                ! of the zero ones (different, so that it stays a 2 by 2 block):
                h = hypot(r, lam1 - p)
                c = r/h
                s = (lam1 - p)/h               ! ((c, s) is lam1's eigenvector)
                if (null1) lam1 = big
                if (null2) lam2 = 2.0_wp*big
                me%ldl(k,k)     = lam1*c*c + lam2*s*s
                me%ldl(k+1,k)   = (lam1 - lam2)*c*s
                me%ldl(k+1,k+1) = lam1*s*s + lam2*c*c
            end if
            k = k + 2
        end if
    end do

    end subroutine count_inertia
!*******************************************************************************

!*******************************************************************************
!>
!  factor `ldl` in place with SQPOPT's Bunch-Kaufman factorization
!  ([[dense_ldl_factor]], which never fails).

    subroutine dense_factor(me, ok)

    class(sqpopt_dense_ldl_type), intent(inout) :: me
    logical, intent(out) :: ok !! whether the factorization succeeded

    call dense_ldl_factor(me%ldl, me%piv)
    ok = .true.

    end subroutine dense_factor
!*******************************************************************************

!*******************************************************************************
!>
!  one solve with the factors ([[dense_ldl_solve]]).

    subroutine dense_solve_factors(me, v, ok)

    class(sqpopt_dense_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: v  !! the right-hand side, overwritten by the solution
    logical,                intent(out)   :: ok !! whether it was solved

    call dense_ldl_solve(me%ldl, me%piv, v)
    ok = .true.

    end subroutine dense_solve_factors
!*******************************************************************************

!*******************************************************************************
!>
!  one solve with the factors, without refinement.

    subroutine dense_back_solve(me, v, ok)

    class(sqpopt_dense_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: v  !! the right-hand side, overwritten by the solution
    logical,                intent(out)   :: ok !! whether it was solved

    call me%solve_factors(v, ok)
    if (ok) ok = sqpopt_all_finite(v)

    end subroutine dense_back_solve
!*******************************************************************************

!*******************************************************************************
!>
!  the product \( y = A x \) with the values last factored.

    subroutine dense_multiply(me, x, y)

    class(sqpopt_dense_ldl_type), intent(in) :: me
    real(wp), dimension(:), intent(in)  :: x !! the vector `dimension(n)`
    real(wp), dimension(:), intent(out) :: y !! the product `dimension(n)`

    integer :: i, j

    do i = 1, me%n
        y(i) = 0.0_wp
    end do
    do j = 1, me%n
        do i = 1, me%n
            y(i) = y(i) + me%a(i,j)*x(j)
        end do
    end do

    end subroutine dense_multiply
!*******************************************************************************

!*******************************************************************************
!>
!  free everything.

    subroutine dense_free(me)

    class(sqpopt_dense_ldl_type), intent(inout) :: me

    if (allocated(me%a))    deallocate(me%a)
    if (allocated(me%ldl))  deallocate(me%ldl)
    if (allocated(me%piv))  deallocate(me%piv)
    if (allocated(me%irow)) deallocate(me%irow)
    if (allocated(me%icol)) deallocate(me%icol)
    me%n = 0

    end subroutine dense_free
!*******************************************************************************

!*******************************************************************************
!>
!  start the LAPACK backend: as [[dense_start]], and `DSYTRF`'s workspace.
!  `ok` is false in a library built without LAPACK.

    subroutine lapack_start(me, n, irow, icol, ok, threads, signs)

    class(sqpopt_lapack_ldl_type), intent(inout) :: me
    integer,               intent(in)  :: n       !! order of the matrix
    integer, dimension(:), intent(in)  :: irow    !! row indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)  :: icol    !! column indices of the entries `dimension(nnz)`
    logical,               intent(out) :: ok      !! whether the backend is ready
    integer, optional,     intent(in)  :: threads !! (not used)
    integer, dimension(:), optional, intent(in) :: signs !! (not used)

#ifdef HAS_LAPACK
    integer :: info, alloc_stat
    real(wp), dimension(1) :: query
#endif

    ok = .false.
    if (present(threads) .or. present(signs)) continue   ! (not used)
#ifdef HAS_LAPACK
    call dense_start(me, n, irow, icol, ok)
    if (.not. ok) return
    ! (the optimal workspace, from a query)
    call dsytrf('L', max(n, 1), me%ldl, max(n, 1), me%piv, query, -1, info)
    allocate(me%work(max(1, int(query(1)))), stat=alloc_stat)
    if (alloc_stat /= 0) then
        me%out_of_memory = .true.
        call me%free()
        ok = .false.
    end if
#endif

    end subroutine lapack_start
!*******************************************************************************

!*******************************************************************************
!>
!  factor `ldl` in place with LAPACK's `DSYTRF`. (An exactly zero pivot,
!  which it reports with `info > 0`, is not a failure: [[count_inertia]]
!  replaces it.)

    subroutine lapack_factor(me, ok)

    class(sqpopt_lapack_ldl_type), intent(inout) :: me
    logical, intent(out) :: ok !! whether the factorization succeeded

#ifdef HAS_LAPACK
    integer :: info
#endif

    ok = .false.
#ifdef HAS_LAPACK
    if (me%n > 0) then
        call dsytrf('L', me%n, me%ldl, me%n, me%piv, me%work, size(me%work), info)
        if (info < 0) return
    end if
    ok = .true.
#endif

    end subroutine lapack_factor
!*******************************************************************************

!*******************************************************************************
!>
!  one solve with the factors, by LAPACK's `DSYTRS`.

    subroutine lapack_solve_factors(me, v, ok)

    class(sqpopt_lapack_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: v  !! the right-hand side, overwritten by the solution
    logical,                intent(out)   :: ok !! whether it was solved

#ifdef HAS_LAPACK
    integer :: info
#endif

    ok = .false.
#ifdef HAS_LAPACK
    if (me%n > 0) then
        call dsytrs('L', me%n, 1, me%ldl, me%n, me%piv, v, me%n, info)
        if (info /= 0) return
    end if
    ok = .true.
#endif

    end subroutine lapack_solve_factors
!*******************************************************************************

    end module sqpopt_dense_ldl_module
!*******************************************************************************
