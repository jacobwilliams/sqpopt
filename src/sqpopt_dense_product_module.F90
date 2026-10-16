!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Dense matrix-vector and matrix-matrix products, in the caller's arrays:
!  the library's replacement for the `matmul` intrinsic, which it doesn't
!  use. An expression with `matmul` (nested, or with a `transpose`, or as an
!  actual argument) needs an array temporary for its result, and some
!  compilers put those on the stack (the Intel compilers by default,
!  gfortran with `-fstack-arrays` or `-Ofast`), which a moderately large
!  matrix then overflows (Windows's default stack is 1 MB, the size of a
!  350 by 350 matrix). These routines form nothing: the result goes into
!  an array the caller has allocated, and `op(A)` is never formed.
!
!  The products are explicit loops (column-oriented, the order of the
!  arrays' storage), in any real kind; in a library built with LAPACK (the
!  `HAS_LAPACK` preprocessor directive, in double precision, which links
!  BLAS too), BLAS's `DGEMV` and `DGEMM` instead, whenever the arrays are
!  contiguous (a section that isn't, e.g. the leading rows of a matrix,
!  uses the loops, since BLAS would need it copied). The two builds'
!  products then differ in rounding, and so can their iterates.

    module sqpopt_dense_product_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

#ifdef HAS_LAPACK
#if defined(REAL32) || defined(REAL128)
#error "HAS_LAPACK needs the default (double precision) real kind: don't define REAL32 or REAL128 with it"
#endif
    interface
        ! (declared `pure`: they change nothing but `y` and `c`)
        pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            !! BLAS: \( y \leftarrow \alpha \, op(A) x + \beta y \)
            import :: wp
            implicit none
            character, intent(in)    :: trans    !! `'N'`: \( op(A)=A \), `'T'`: \( op(A)=A^T \)
            integer,   intent(in)    :: m        !! number of rows of `a`
            integer,   intent(in)    :: n        !! number of columns of `a`
            real(wp),  intent(in)    :: alpha    !! the product's factor
            integer,   intent(in)    :: lda      !! leading dimension of `a`
            real(wp),  intent(in)    :: a(lda,*) !! the matrix
            real(wp),  intent(in)    :: x(*)     !! the vector
            integer,   intent(in)    :: incx     !! the stride of `x`
            real(wp),  intent(in)    :: beta     !! the factor of `y` (`0`: `y` isn't read)
            real(wp),  intent(inout) :: y(*)     !! the result
            integer,   intent(in)    :: incy     !! the stride of `y`
        end subroutine dgemv
        pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            !! BLAS: \( C \leftarrow \alpha \, op(A) \, op(B) + \beta C \)
            import :: wp
            implicit none
            character, intent(in)    :: transa   !! `'N'`: \( op(A)=A \), `'T'`: \( op(A)=A^T \)
            character, intent(in)    :: transb   !! the same for `b`
            integer,   intent(in)    :: m        !! number of rows of `c`
            integer,   intent(in)    :: n        !! number of columns of `c`
            integer,   intent(in)    :: k        !! the inner dimension
            real(wp),  intent(in)    :: alpha    !! the product's factor
            integer,   intent(in)    :: lda      !! leading dimension of `a`
            real(wp),  intent(in)    :: a(lda,*) !! the first matrix
            integer,   intent(in)    :: ldb      !! leading dimension of `b`
            real(wp),  intent(in)    :: b(ldb,*) !! the second matrix
            real(wp),  intent(in)    :: beta     !! the factor of `c` (`0`: `c` isn't read)
            integer,   intent(in)    :: ldc      !! leading dimension of `c`
            real(wp),  intent(inout) :: c(ldc,*) !! the result
        end subroutine dgemm
    end interface
#endif

    public :: dense_matvec
    public :: dense_matmul

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  \( y \leftarrow \alpha \, op(A) x + \beta y \), with \( op(A) = A \) or,
!  with `trans`, \( A^T \) (not formed). `y` must not overlap `a` or `x`.
!  With `beta = 0` (the default) `y` is only written, so it may hold
!  anything (even NaNs) on entry.

    pure subroutine dense_matvec(a, x, y, trans, alpha, beta)

    real(wp), dimension(:,:), intent(in)    :: a     !! the matrix `dimension(m,n)`
    real(wp), dimension(:),   intent(in)    :: x     !! the vector `dimension(n)` (`dimension(m)` with `trans`)
    real(wp), dimension(:),   intent(inout) :: y     !! the result `dimension(m)` (`dimension(n)` with `trans`)
    logical,  optional,       intent(in)    :: trans !! whether to multiply by \( A^T \) (default `.false.`)
    real(wp), optional,       intent(in)    :: alpha !! the product's factor (default `1`)
    real(wp), optional,       intent(in)    :: beta  !! the factor of `y` on entry (default `0`)

    real(wp) :: al, be
    logical :: tr
    integer :: j

    tr = .false.
    if (present(trans)) tr = trans
    al = 1.0_wp
    if (present(alpha)) al = alpha
    be = 0.0_wp
    if (present(beta)) be = beta

    if (size(y) == 0) return

#ifdef HAS_LAPACK
    if (size(a) > 0 .and. is_contiguous(a) .and. is_contiguous(x) .and. is_contiguous(y)) then
        call dgemv(merge('T', 'N', tr), size(a,1), size(a,2), al, a, size(a,1), x, 1, be, y, 1)
        return
    end if
#endif

    if (tr) then
        do j = 1, size(y)
            if (be == 0.0_wp) then
                y(j) = al*dot_product(a(:,j), x)
            else
                y(j) = al*dot_product(a(:,j), x) + be*y(j)
            end if
        end do
    else
        if (be == 0.0_wp) then
            y = 0.0_wp
        else if (be /= 1.0_wp) then
            y = be*y
        end if
        do j = 1, size(x)
            y = y + (al*x(j))*a(:,j)
        end do
    end if

    end subroutine dense_matvec
!*******************************************************************************

!*******************************************************************************
!>
!  \( C \leftarrow \alpha \, op(A) \, op(B) + \beta C \), with
!  \( op(A) = A \) or, with `trans_a`, \( A^T \) (not formed), and the same
!  for `b`. `c` must not overlap `a` or `b`. With `beta = 0` (the default)
!  `c` is only written, so it may hold anything (even NaNs) on entry.

    pure subroutine dense_matmul(a, b, c, trans_a, trans_b, alpha, beta)

    real(wp), dimension(:,:), intent(in)    :: a       !! the first matrix: `dimension(m,k)` (`dimension(k,m)`
                                                       !! with `trans_a`)
    real(wp), dimension(:,:), intent(in)    :: b       !! the second matrix: `dimension(k,n)` (`dimension(n,k)`
                                                       !! with `trans_b`)
    real(wp), dimension(:,:), intent(inout) :: c       !! the result `dimension(m,n)`
    logical,  optional,       intent(in)    :: trans_a !! whether to use \( A^T \) (default `.false.`)
    logical,  optional,       intent(in)    :: trans_b !! whether to use \( B^T \) (default `.false.`)
    real(wp), optional,       intent(in)    :: alpha   !! the product's factor (default `1`)
    real(wp), optional,       intent(in)    :: beta    !! the factor of `c` on entry (default `0`)

    real(wp) :: al, be
    logical :: ta, tb
    integer :: i, j, l, m, n, k

    ta = .false.
    if (present(trans_a)) ta = trans_a
    tb = .false.
    if (present(trans_b)) tb = trans_b
    al = 1.0_wp
    if (present(alpha)) al = alpha
    be = 0.0_wp
    if (present(beta)) be = beta

    m = size(c,1)
    n = size(c,2)
    k = merge(size(a,1), size(a,2), ta)
    if (m == 0 .or. n == 0) return

#ifdef HAS_LAPACK
    if (k > 0 .and. is_contiguous(a) .and. is_contiguous(b) .and. is_contiguous(c)) then
        call dgemm(merge('T', 'N', ta), merge('T', 'N', tb), m, n, k, al, a, size(a,1), b, size(b,1), &
                   be, c, size(c,1))
        return
    end if
#endif

    if (be == 0.0_wp) then
        c = 0.0_wp
    else if (be /= 1.0_wp) then
        c = be*c
    end if
    if (ta) then
        ! (each element a dot product of two columns, or of a column and a row)
        do j = 1, n
            do i = 1, m
                if (tb) then
                    c(i,j) = c(i,j) + al*dot_product(a(:,i), b(j,:))
                else
                    c(i,j) = c(i,j) + al*dot_product(a(:,i), b(:,j))
                end if
            end do
        end do
    else
        ! (each column a combination of the columns of `a`)
        do j = 1, n
            do l = 1, k
                if (tb) then
                    c(:,j) = c(:,j) + (al*b(j,l))*a(:,l)
                else
                    c(:,j) = c(:,j) + (al*b(l,j))*a(:,l)
                end if
            end do
        end do
    end if

    end subroutine dense_matmul
!*******************************************************************************

    end module sqpopt_dense_product_module
!*******************************************************************************
