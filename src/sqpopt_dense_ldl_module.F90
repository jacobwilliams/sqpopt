!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The dense backend of [[sqpopt_symmetric_solver_module]]
!  (`options%linear_solver = sqpopt_linear_solver_dense`), for small
!  matrices: the matrix is stored dense (order `n`, so \( n^2 \) elements),
!  its inertia is the exact count of [[dense_symmetric_inertia]]
!  (Householder tridiagonalization, then the signs of a Sturm sequence, with
!  2 by 2 pivots for coupled zeros, as Bunch-Kaufman's), and solves use an
!  LU factorization with partial pivoting. So, unlike QDLDL, it always
!  knows the inertia, also of a Hessian with zeros on its diagonal, and,
!  unlike MUMPS, it has almost no overhead per call. But each factorization
!  costs \( O(n^3) \): it is meant for matrices of order up to a few hundred,
!  and refuses to start above `dense_max_order` (the factorization options
!  then fall back on the matrix-free methods).
!
!  A singular matrix is not an error: its zero eigenvalues are counted
!  (`n_null`), and in the LU factorization a null column gets a large
!  pivot, which makes that component of a solution zero, so that a solve
!  returns one of the solutions of a consistent system (as MUMPS's
!  null-pivot detection does).
!
!  **With LAPACK** (`options%linear_solver = sqpopt_linear_solver_lapack`,
!  [[sqpopt_lapack_ldl_type]]): the same dense matrix is factored by
!  LAPACK's `DSYTRF` (Bunch-Kaufman \( LDL^T \), with 1 by 1 and 2 by 2
!  pivots), and solved by `DSYTRS`. One factorization, of about \( n^3/3 \)
!  operations with blocked, BLAS-level code, gives both the solves and the
!  inertia (the signs of the eigenvalues of `D`'s blocks, by Sylvester's
!  law), in place of the tridiagonalization and the LU. A null 1 by 1 pivot
!  gets a large value, as in the LU; a null 2 by 2 block (rare) makes the
!  solves fall back on the LU. This needs a library built with LAPACK and
!  BLAS (the `HAS_LAPACK` preprocessor directive, and `-llapack -lblas`),
!  in double precision; without it, `sqpopt_has_lapack` is false and the
!  option can't be chosen. This module is the only one that refers to
!  LAPACK.

    module sqpopt_dense_ldl_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_sparse_ldl_module, only: sqpopt_sparse_ldl_type
    use sqpopt_dense_linalg_module, only: dense_symmetric_inertia
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
                                                                     !! element is a zero eigenvalue (well above the
                                                                     !! roundoff of the tridiagonal reduction, and
                                                                     !! well below the constraint block's eigenvalues
                                                                     !! of a KKT matrix with a large Hessian shift)
    integer, parameter, public :: dense_max_order = 2000 !! the largest order the dense backend accepts (its two
                                                         !! matrices then take 64 MB in double precision)

    type, extends(sqpopt_sparse_ldl_type), public :: sqpopt_dense_ldl_type
        !! the dense backend (see the module documentation)
        private
        integer :: n = 0                                 !! order of the matrix
        integer, dimension(:), allocatable :: irow       !! row indices of the pattern's entries
        integer, dimension(:), allocatable :: icol       !! column indices of the pattern's entries
        real(wp), dimension(:,:), allocatable :: a       !! the matrix last factored `dimension(n,n)`
        real(wp), dimension(:,:), allocatable :: lu      !! its LU factors `dimension(n,n)`
        integer,  dimension(:),   allocatable :: piv     !! the row pivots of the LU factorization `dimension(n)`
        contains
        procedure :: start      => dense_start
        procedure :: refactor   => dense_refactor
        procedure :: back_solve => dense_back_solve
        procedure :: multiply   => dense_multiply
        procedure :: free       => dense_free
    end type sqpopt_dense_ldl_type

    type, extends(sqpopt_dense_ldl_type), public :: sqpopt_lapack_ldl_type
        !! the dense backend with LAPACK's `DSYTRF`/`DSYTRS` (see the module
        !! documentation); it can only be started in a library built with
        !! `HAS_LAPACK`
        private
        real(wp), dimension(:), allocatable :: work !! `DSYTRF`'s workspace
        logical :: lu_solves = .false.              !! whether the last factorization had a null 2 by 2 block, so
                                                    !! the solves use the LU factors instead
        contains
        procedure :: start      => lapack_start
        procedure :: refactor   => lapack_refactor
        procedure :: back_solve => lapack_back_solve
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
    if (alloc_stat == 0) allocate(me%lu(n,n), stat=alloc_stat)
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
!  form the dense matrix from the values `val`, count its inertia, and
!  factor it (see the module documentation).

    subroutine dense_refactor(me, val, ok)

    class(sqpopt_dense_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: val !! values of the entries `dimension(nnz)`
    logical,                intent(out) :: ok  !! whether the factorization succeeded

    integer :: k, i, j, n_positive

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
    call dense_symmetric_inertia(me%a, n_positive, me%n_negative, me%n_null, zero_tol=dense_zero_tol)
    do j = 1, me%n
        do i = 1, me%n
            me%lu(i,j) = me%a(i,j)
        end do
    end do
    call lu_factor_null_safe(me%lu, me%piv)
    ok = .true.

    end subroutine dense_refactor
!*******************************************************************************

!*******************************************************************************
!>
!  in-place LU factorization with partial pivoting of the dense matrix
!  `a`, in which a negligible pivot (a null column, below the diagonal) is
!  replaced by a large one, so that the factorization never fails (see the
!  module documentation).

    pure subroutine lu_factor_null_safe(a, piv)

    real(wp), dimension(:,:), intent(inout) :: a   !! matrix, overwritten by its `L` (unit, below the diagonal)
                                                   !! and `U` factors
    integer,  dimension(:),   intent(out)   :: piv !! `piv(p)` is the row swapped with row `p` at step `p`

    integer :: i, j, p, n
    real(wp) :: amax, tol, big, t

    n = size(a,1)
    amax = 0.0_wp
    do j = 1, n
        do i = 1, n
            amax = max(amax, abs(a(i,j)))
        end do
    end do
    tol = 1.0e-14_wp*max(amax, tiny(1.0_wp))
    big = max(amax, 1.0_wp)/sqrt(epsilon(1.0_wp))
    do p = 1, n
        ! the pivot: the largest element of column p, on or below the diagonal
        piv(p) = p
        do i = p+1, n
            if (abs(a(i,p)) > abs(a(piv(p),p))) piv(p) = i
        end do
        if (piv(p) /= p) then
            do j = 1, n
                t = a(p,j)
                a(p,j) = a(piv(p),j)
                a(piv(p),j) = t
            end do
        end if
        if (abs(a(p,p)) <= tol) then
            a(p,p) = big   ! (a null column: its component of a solution becomes negligible)
        end if
        do i = p+1, n
            a(i,p) = a(i,p)/a(p,p)
            do j = p+1, n
                a(i,j) = a(i,j) - a(i,p)*a(p,j)
            end do
        end do
    end do

    end subroutine lu_factor_null_safe
!*******************************************************************************

!*******************************************************************************
!>
!  one solve with the LU factors, without refinement.

    subroutine dense_back_solve(me, v, ok)

    class(sqpopt_dense_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: v  !! the right-hand side, overwritten by the solution
    logical,                intent(out)   :: ok !! whether it was solved

    integer :: i, j
    real(wp) :: s, t

    do i = 1, me%n
        if (me%piv(i) /= i) then
            t = v(i)
            v(i) = v(me%piv(i))
            v(me%piv(i)) = t
        end if
    end do
    do i = 2, me%n
        s = v(i)
        do j = 1, i-1
            s = s - me%lu(i,j)*v(j)
        end do
        v(i) = s
    end do
    do i = me%n, 1, -1
        s = v(i)
        do j = i+1, me%n
            s = s - me%lu(i,j)*v(j)
        end do
        v(i) = s/me%lu(i,i)
    end do
    ok = sqpopt_all_finite(v)

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
    if (allocated(me%lu))   deallocate(me%lu)
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
    call dsytrf('L', max(n, 1), me%lu, max(n, 1), me%piv, query, -1, info)
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
!  form the dense matrix from the values `val`, and factor it with
!  `DSYTRF`; the inertia is counted from `D`'s blocks (see the module
!  documentation).

    subroutine lapack_refactor(me, val, ok)

    class(sqpopt_lapack_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: val !! values of the entries `dimension(nnz)`
    logical,                intent(out) :: ok  !! whether the factorization succeeded

#ifdef HAS_LAPACK
    integer :: k, i, j, info, n
    real(wp) :: amax, tol, big, a11, a21, a22, det
#endif

    ok = .false.
    me%inertia_known = .true.
#ifdef HAS_LAPACK
    if (.not. sqpopt_all_finite(val)) return
    n = me%n
    do j = 1, n
        do i = 1, n
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
    do j = 1, n
        do i = 1, n
            me%lu(i,j) = me%a(i,j)
            amax = max(amax, abs(me%a(i,j)))
        end do
    end do
    tol = dense_zero_tol*max(amax, tiny(1.0_wp))
    big = max(amax, 1.0_wp)/sqrt(epsilon(1.0_wp))
    me%n_negative = 0
    me%n_null     = 0
    me%lu_solves  = .false.
    if (n > 0) then
        call dsytrf('L', n, me%lu, n, me%piv, me%work, size(me%work), info)
        if (info < 0) return
    end if
    ! the inertia, from D's blocks (a 2 by 2 block has piv(k) = piv(k+1) < 0):
    k = 1
    do while (k <= n)
        if (me%piv(k) > 0) then
            a11 = me%lu(k,k)
            if (abs(a11) <= tol) then
                me%n_null = me%n_null + 1
                me%lu(k,k) = big   ! (a null pivot: its component of a solution becomes negligible)
            else if (a11 < 0.0_wp) then
                me%n_negative = me%n_negative + 1
            end if
            k = k + 1
        else
            a11 = me%lu(k,k)
            a21 = me%lu(k+1,k)
            a22 = me%lu(k+1,k+1)
            det = a11*a22 - a21*a21
            if (abs(det) <= tol*max(abs(a11), abs(a21), abs(a22))) then
                ! (a null 2 by 2 block: one zero eigenvalue, the other with the sign of
                ! the trace; DSYTRS can't solve with it, so the LU does)
                me%n_null = me%n_null + 1
                if (a11 + a22 < 0.0_wp) me%n_negative = me%n_negative + 1
                me%lu_solves = .true.
            else if (det < 0.0_wp) then
                me%n_negative = me%n_negative + 1         ! (one positive, one negative eigenvalue)
            else if (a11 + a22 < 0.0_wp) then
                me%n_negative = me%n_negative + 2
            end if
            k = k + 2
        end if
    end do
    if (me%lu_solves) then
        ! (the LU factors, for the solves; the inertia stays DSYTRF's)
        do j = 1, n
            do i = 1, n
                me%lu(i,j) = me%a(i,j)
            end do
        end do
        call lu_factor_null_safe(me%lu, me%piv)
    end if
    ok = .true.
#endif

    end subroutine lapack_refactor
!*******************************************************************************

!*******************************************************************************
!>
!  one solve with `DSYTRS` (or with the LU factors, after a null 2 by 2
!  block), without refinement.

    subroutine lapack_back_solve(me, v, ok)

    class(sqpopt_lapack_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: v  !! the right-hand side, overwritten by the solution
    logical,                intent(out)   :: ok !! whether it was solved

#ifdef HAS_LAPACK
    integer :: info
#endif

    ok = .false.
    if (me%lu_solves) then
        call dense_back_solve(me, v, ok)
        return
    end if
#ifdef HAS_LAPACK
    if (me%n > 0) then
        call dsytrs('L', me%n, 1, me%lu, me%n, me%piv, v, me%n, info)
        if (info /= 0) return
    end if
    ok = sqpopt_all_finite(v)
#endif

    end subroutine lapack_back_solve
!*******************************************************************************

    end module sqpopt_dense_ldl_module
!*******************************************************************************
