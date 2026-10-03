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

    module sqpopt_dense_ldl_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_sparse_ldl_module, only: sqpopt_sparse_ldl_type
    use sqpopt_dense_linalg_module, only: dense_symmetric_inertia
    use sqpopt_types_module, only: sqpopt_all_finite

    implicit none

    private

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

    end module sqpopt_dense_ldl_module
!*******************************************************************************
