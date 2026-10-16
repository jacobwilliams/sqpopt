!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The QDLDL backend of [[sqpopt_symmetric_solver_module]]
!  (`options%linear_solver = sqpopt_linear_solver_qdldl`, the default):
!  [QDLDL](https://github.com/jacobwilliams/qdldl-fortran), a Fortran port
!  of the \( LDL^T \) solver of the OSQP QP solver, an fpm dependency, so it
!  is always available, in any real kind. It factors in a fill-reducing
!  order (approximate minimum degree, with dense rows last), *without
!  pivoting*, on one thread, with very little overhead per call.
!
!  Without pivoting, it is exact for quasi-definite matrices, and for the
!  KKT matrices \( \begin{bmatrix} H & J^T \\ J & 0 \end{bmatrix} \) with
!  \( H \) positive definite once each constraint's row comes after its
!  variables (the `signs` of [[qdldl_start]], see [[delay_negative_rows]]):
!  the least-squares systems, and the matrices of the quasi-Newton
!  Hessians. For an indefinite \( H \) the inertia is still exact as long
!  as no pivot is zero (Sylvester's law).
!
!  **Null pivots.** A pivot of at most `qdldl_zero_tol` times the largest
!  element is null. The first one stops a first factorization, which is
!  then repeated with each null pivot replaced by a large one (which makes
!  its row and column of the factors negligible, so that a solve returns
!  one of the solutions of a consistent system, as with MUMPS). Whether
!  those factors are still the matrix's is then checked (see
!  [[qdldl_faithful]]): they are if the matrix is singular, and `n_null`
!  counts its zero eigenvalues. If they aren't, the zero pivot was an
!  artifact of the order (e.g. a Hessian with zeros on its diagonal, which
!  would need a 2x2 pivot), and the inertia is unknown:
!  `inertia_known = .false.` (what that means for the optimization is
!  decided by [[sqpopt_kkt_module]]).

    module sqpopt_qdldl_ldl_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_sparse_ldl_module, only: sqpopt_sparse_ldl_type
    use qdldl_module, only: qdldl_type, qdldl_wp, qdldl_success, qdldl_error_out_of_memory, qdldl_order_amd, &
                            qdldl_error_zero_pivot

    implicit none

    private

    integer, parameter :: kind_check = 1/merge(1, 0, qdldl_wp == wp) !! (a compile-time error if QDLDL was built
                                                                     !! with another real kind)
    real(wp), parameter :: qdldl_zero_tol = 1.0e-5_wp*epsilon(1.0_wp) !! a pivot at most this times the largest
                                                                      !! element is a null pivot

    public :: delay_negative_rows

    type, extends(sqpopt_sparse_ldl_type), public :: sqpopt_qdldl_ldl_type
        !! the QDLDL backend (see the module documentation)
        private
        integer :: n = 0          !! order of the matrix
        type(qdldl_type) :: ldl   !! the QDLDL solver (with the pattern, the values, and the factors)
        contains
        procedure :: start      => qdldl_start
        procedure :: refactor   => qdldl_refactor
        procedure :: back_solve => qdldl_back_solve
        procedure :: multiply   => qdldl_multiply
        procedure :: free       => qdldl_free
    end type sqpopt_qdldl_ldl_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  analyse the pattern: AMD, then (with `signs`) each row with the sign
!  `-1` moved after all of its neighbours (see [[delay_negative_rows]]).
!  `threads` is not used (QDLDL is single-threaded).

    subroutine qdldl_start(me, n, irow, icol, ok, threads, signs)

    class(sqpopt_qdldl_ldl_type), intent(inout) :: me
    integer,               intent(in)  :: n       !! order of the matrix
    integer, dimension(:), intent(in)  :: irow    !! row indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)  :: icol    !! column indices of the entries `dimension(nnz)`
    logical,               intent(out) :: ok      !! whether the backend is ready
    integer, optional,     intent(in)  :: threads !! (not used)
    integer, dimension(:), optional, intent(in) :: signs !! the expected sign of each row's pivot `dimension(n)`

    integer :: istat
    integer, dimension(:), allocatable :: perm

    ok = .false.
    if (present(threads)) continue   ! (QDLDL is single-threaded)
    call me%ldl%analyze(n, irow, icol, istat, ordering=qdldl_order_amd)
    if (istat == qdldl_success .and. present(signs)) then
        if (size(signs) == n) then
            if (any(signs == -1)) then
                call me%ldl%get_permutation(perm)
                call delay_negative_rows(n, irow, icol, signs, perm)
                call me%ldl%analyze(n, irow, icol, istat, perm=perm)
            end if
        end if
    end if
    if (istat /= qdldl_success) then
        if (istat == qdldl_error_out_of_memory) me%out_of_memory = .true.
        call me%ldl%destroy()
        return
    end if
    ! (null pivots: see the module documentation; a pivot of the wrong sign
    ! is kept, so the inertia is that of the matrix)
    me%ldl%reg_eps        = 0.0_wp
    me%ldl%zero_pivot_tol = qdldl_zero_tol
    me%ldl%max_refine     = 0   ! (the refinement is the symmetric solver's, for every backend)
    me%n = n
    ok = .true.

    end subroutine qdldl_start
!*******************************************************************************

!*******************************************************************************
!>
!  change the fill-reducing order `perm` (row `perm(k)` of the matrix is
!  the `k`-th to be eliminated) so that each row `i` with `signs(i) == -1`
!  comes after all of its neighbours (the rows `j` with an entry `(i,j)`),
!  and otherwise as early as it can. In a KKT matrix
!  \( \begin{bmatrix} H & J^T \\ J & 0 \end{bmatrix} \) with \( H \) positive
!  definite, the pivot of a constraint's row is then an element of
!  \( -J H^{-1} J^T \) (less the earlier constraints' parts), which is zero
!  only if the rows of \( J \) are dependent. Without the delay, a
!  constraint's row eliminated before its variables has a zero pivot even
!  when the matrix is nonsingular, which a solver without pivoting can't
!  get around.

    subroutine delay_negative_rows(n, irow, icol, signs, perm)

    integer,               intent(in)    :: n     !! order of the matrix
    integer, dimension(:), intent(in)    :: irow  !! row indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)    :: icol  !! column indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)    :: signs !! the expected sign of each row's pivot `dimension(n)`
    integer, dimension(:), intent(inout) :: perm  !! the order, changed in place `dimension(n)`

    integer, dimension(:), allocatable :: ptr, adj, waiting, new_perm, fill
    logical, dimension(:), allocatable :: done
    integer :: k, i, j, v, c, p, nnew

    ! the neighbours of each row (both directions), without the diagonal:
    allocate(ptr(n+1), fill(n), waiting(n), done(n), new_perm(n))
    ptr = 0
    do k = 1, size(irow)
        i = irow(k)
        j = icol(k)
        if (i == j) cycle
        ptr(i) = ptr(i) + 1
        ptr(j) = ptr(j) + 1
    end do
    p = 1
    do i = 1, n
        fill(i) = p
        p = p + ptr(i)
        ptr(i) = fill(i)
    end do
    ptr(n+1) = p
    allocate(adj(max(p-1, 0)))
    do k = 1, size(irow)
        i = irow(k)
        j = icol(k)
        if (i == j) cycle
        adj(fill(i)) = j
        fill(i) = fill(i) + 1
        adj(fill(j)) = i
        fill(j) = fill(j) + 1
    end do

    ! how many neighbours of each delayed row are still to come (duplicates
    ! counted, consistently, in both places):
    waiting = 0
    do i = 1, n
        if (signs(i) /= -1) cycle
        do p = ptr(i), ptr(i+1) - 1
            if (signs(adj(p)) /= -1) waiting(i) = waiting(i) + 1
        end do
    end do

    done = .false.
    nnew = 0
    do k = 1, n
        v = perm(k)
        if (done(v)) cycle
        if (signs(v) == -1) then
            if (waiting(v) == 0) call emit(v)   ! (else it is emitted after its last neighbour)
        else
            call emit(v)
            do p = ptr(v), ptr(v+1) - 1
                c = adj(p)
                if (signs(c) /= -1) cycle
                waiting(c) = waiting(c) - 1
                if (waiting(c) == 0 .and. .not. done(c)) call emit(c)
            end do
        end if
    end do
    ! (anything left: rows whose neighbours were all delayed rows themselves)
    do k = 1, n
        if (.not. done(perm(k))) call emit(perm(k))
    end do
    perm = new_perm

    contains

        subroutine emit(row)
        !! append `row` to the new order
        integer, intent(in) :: row !! the row
        nnew = nnew + 1
        new_perm(nnew) = row
        done(row) = .true.
        end subroutine emit

    end subroutine delay_negative_rows
!*******************************************************************************

!*******************************************************************************
!>
!  factor the matrix with the values `val`, and set its inertia (see the
!  module documentation for the null pivots, and `inertia_known`).

    subroutine qdldl_refactor(me, val, ok)

    class(sqpopt_qdldl_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: val !! values of the entries `dimension(nnz)`
    logical,                intent(out) :: ok  !! whether the factorization succeeded

    integer  :: istat
    integer  :: null_row !! the row of the first null pivot (0: none)
    real(wp) :: amax     !! the largest element

    ok = .false.
    me%inertia_known = .true.
    amax = 0.0_wp
    if (size(val) > 0) amax = maxval(abs(val))

    ! first without replacing null pivots: the first one stops the
    ! factorization, and tells its row
    me%ldl%regularize = .false.
    call me%ldl%factor(val, istat)
    null_row = 0
    if (istat == qdldl_error_zero_pivot) then
        null_row = me%ldl%zero_pivot_column
        ! again, with each null pivot replaced by a large one:
        me%ldl%regularize = .true.
        me%ldl%reg_delta  = max(amax, 1.0_wp)/sqrt(epsilon(1.0_wp))
        call me%ldl%factor(val, istat)
    end if
    if (istat /= qdldl_success) return
    me%n_negative = me%ldl%n_negative
    me%n_null     = me%ldl%n_zero
    if (null_row > 0) me%inertia_known = qdldl_faithful(me, amax)
    ok = .true.

    end subroutine qdldl_refactor
!*******************************************************************************

!*******************************************************************************
!>
!  after null pivots were replaced: whether the factors are still those of
!  the matrix (as they are if the matrix is singular, since the rest of a
!  null pivot's row is then zero too): for a fixed `x` and `b = A x`, the
!  solve with the factors satisfies `A x' = b` (`x'` may differ from `x`
!  along the null space).

    logical function qdldl_faithful(me, amax) result(faithful)

    class(sqpopt_qdldl_ldl_type), intent(inout) :: me
    real(wp), intent(in) :: amax !! the largest element of the matrix

    real(wp), parameter :: rtol = 1.0e-6_wp !! relative tolerance on the residual
    real(wp), dimension(:), allocatable :: x, b, ax
    integer :: i, istat

    allocate(x(me%n), b(me%n), ax(me%n))
    do i = 1, me%n
        x(i) = 1.0_wp + real(mod(i, 7), wp)/7.0_wp
    end do
    call me%ldl%multiply(x, b)
    x = b
    call me%ldl%solve(x, istat, refine=0)
    faithful = istat == qdldl_success
    if (.not. faithful) return
    call me%ldl%multiply(x, ax)
    faithful = maxval(abs(ax - b)) <= rtol*(amax*maxval(abs(x)) + maxval(abs(b)))

    end function qdldl_faithful
!*******************************************************************************

!*******************************************************************************
!>
!  one solve with the factors, without refinement.

    subroutine qdldl_back_solve(me, v, ok)

    class(sqpopt_qdldl_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: v  !! the right-hand side, overwritten by the solution
    logical,                intent(out)   :: ok !! whether it was solved

    integer :: istat

    call me%ldl%solve(v, istat, refine=0)
    ok = istat == qdldl_success

    end subroutine qdldl_back_solve
!*******************************************************************************

!*******************************************************************************
!>
!  the product \( y = A x \) with the values last factored.

    subroutine qdldl_multiply(me, x, y)

    class(sqpopt_qdldl_ldl_type), intent(in) :: me
    real(wp), dimension(:), intent(in)  :: x !! the vector `dimension(n)`
    real(wp), dimension(:), intent(out) :: y !! the product `dimension(n)`

    integer :: istat

    call me%ldl%multiply(x, y, istat)
    if (istat /= qdldl_success) y = 0.0_wp

    end subroutine qdldl_multiply
!*******************************************************************************

!*******************************************************************************
!>
!  free everything.

    subroutine qdldl_free(me)

    class(sqpopt_qdldl_ldl_type), intent(inout) :: me

    call me%ldl%destroy()
    me%n = 0

    end subroutine qdldl_free
!*******************************************************************************

    end module sqpopt_qdldl_ldl_module
!*******************************************************************************
