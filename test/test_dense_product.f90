program test_dense_product

    !! Tests of the dense products ([[sqpopt_dense_product_module]]), the
    !! library's replacement for `matmul`, against `matmul` itself: in the
    !! default build they are explicit loops, and in the build with LAPACK
    !! BLAS's `DGEMV` and `DGEMM` (for contiguous arrays) or the loops (for
    !! sections that aren't).
    !!
    !! * `dense_matvec`, with and without the transpose, and `dense_matmul`
    !!   with each combination of transposes, on random rectangular matrices
    !!   of several sizes (one row or column, and empty ones, included),
    !!   with the default factors (`alpha = 1`, `beta = 0`) and others;
    !! * with `beta = 0` the result isn't read (it may hold NaNs);
    !! * the same on non-contiguous sections of larger arrays;
    !! * a 600 by 600 product and matrix-vector product (the matrices larger
    !!   than Windows's default 1 MB stack), checked against sums computed
    !!   here, without `matmul`.

    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use sqpopt_kinds,                only: wp => sqpopt_module_wp
    use sqpopt_dense_product_module, only: dense_matvec, dense_matmul

    implicit none

    integer, dimension(*), parameter :: sizes = [0, 1, 2, 7, 33]
    integer :: n_fail, i, j, l, it, ia, ib
    integer, allocatable :: seed(:)

    write(*,*) '----------------------------'
    write(*,*) 'test_dense_product'
    write(*,*) '----------------------------'

    call random_seed(size=i)
    allocate(seed(i))
    seed = [(17 + 3*j, j = 1, i)]
    call random_seed(put=seed)
    n_fail = 0

    do i = 1, size(sizes)
        do j = 1, size(sizes)
            do it = 0, 1
                call check_matvec(sizes(i), sizes(j), it == 1)
            end do
            do l = 1, size(sizes)
                do ia = 0, 1
                    do ib = 0, 1
                        call check_matmul(sizes(i), sizes(j), sizes(l), ia == 1, ib == 1)
                    end do
                end do
            end do
        end do
    end do
    call check_large()

    if (n_fail > 0) then
        write(*,'(I0,A)') n_fail, ' check(s) failed'
        error stop 'test_dense_product FAILED'
    end if
    print '(A)', 'test_dense_product PASSED'

    contains

    subroutine check_matvec(m, n, trans)
    !! `y = alpha*op(A)*x + beta*y` with an `m` by `n` matrix: contiguous,
    !! and as sections of larger arrays
    integer, intent(in) :: m     !! number of rows of the matrix
    integer, intent(in) :: n     !! number of columns of the matrix
    logical, intent(in) :: trans !! whether to multiply by the transpose
    real(wp), dimension(:,:), allocatable :: a, big
    real(wp), dimension(:), allocatable :: x, y, y0, expected, xbig, ybig
    integer :: nx, ny

    allocate(a(m,n), big(m+3,n+2))
    call random_number(a)
    a = a - 0.5_wp
    nx = merge(m, n, trans)
    ny = merge(n, m, trans)
    allocate(x(nx), y(ny), y0(ny), expected(ny), xbig(2*nx), ybig(2*ny))
    call random_number(x)
    call random_number(y0)
    if (trans) then
        expected = matmul(transpose(a), x)
    else
        expected = matmul(a, x)
    end if

    ! (the defaults, with NaNs in `y`, which `beta = 0` mustn't read)
    y = ieee_value(1.0_wp, ieee_quiet_nan)
    call dense_matvec(a, x, y, trans=trans)
    call compare(y, expected, max(m,n), 'matvec, defaults')

    y = y0
    call dense_matvec(a, x, y, trans=trans, alpha=-2.5_wp, beta=0.75_wp)
    call compare(y, -2.5_wp*expected + 0.75_wp*y0, max(m,n), 'matvec, alpha and beta')

    y = y0
    call dense_matvec(a, x, y, trans=trans, beta=1.0_wp)
    call compare(y, expected + y0, max(m,n), 'matvec, beta = 1')

    ! (non-contiguous sections: the leading rows of a matrix, every other element of the vectors)
    big = 0.0_wp
    big(1:m,1:n) = a
    xbig = 0.0_wp
    xbig(1::2) = x
    ybig = 0.0_wp
    ybig(1::2) = y0
    call dense_matvec(big(1:m,1:n), xbig(1::2), ybig(1::2), trans=trans, alpha=0.5_wp, beta=-1.0_wp)
    call compare(ybig(1::2), 0.5_wp*expected - y0, max(m,n), 'matvec, sections')
    if (any(ybig(2::2) /= 0.0_wp)) call fail('matvec, sections: an element outside the section changed')

    end subroutine check_matvec

    subroutine check_matmul(m, n, k, trans_a, trans_b)
    !! `C = alpha*op(A)*op(B) + beta*C` with an `m` by `n` result and the
    !! inner dimension `k`: contiguous, and as sections of larger arrays
    integer, intent(in) :: m       !! number of rows of the result
    integer, intent(in) :: n       !! number of columns of the result
    integer, intent(in) :: k       !! the inner dimension
    logical, intent(in) :: trans_a !! whether to use the transpose of `a`
    logical, intent(in) :: trans_b !! whether to use the transpose of `b`
    real(wp), dimension(:,:), allocatable :: a, b, c, c0, expected, abig, bbig, cbig

    if (trans_a) then
        allocate(a(k,m))
    else
        allocate(a(m,k))
    end if
    if (trans_b) then
        allocate(b(n,k))
    else
        allocate(b(k,n))
    end if
    call random_number(a)
    call random_number(b)
    a = a - 0.5_wp
    allocate(c(m,n), c0(m,n), expected(m,n))
    call random_number(c0)
    if (trans_a .and. trans_b) then
        expected = matmul(transpose(a), transpose(b))
    else if (trans_a) then
        expected = matmul(transpose(a), b)
    else if (trans_b) then
        expected = matmul(a, transpose(b))
    else
        expected = matmul(a, b)
    end if

    c = ieee_value(1.0_wp, ieee_quiet_nan)
    call dense_matmul(a, b, c, trans_a=trans_a, trans_b=trans_b)
    call compare(reshape(c, [m*n]), reshape(expected, [m*n]), k, 'matmul, defaults')

    c = c0
    call dense_matmul(a, b, c, trans_a=trans_a, trans_b=trans_b, alpha=3.0_wp, beta=-0.5_wp)
    call compare(reshape(c, [m*n]), reshape(3.0_wp*expected - 0.5_wp*c0, [m*n]), k, 'matmul, alpha and beta')

    ! (the leading rows of larger arrays)
    allocate(abig(size(a,1)+2, size(a,2)), bbig(size(b,1)+1, size(b,2)), cbig(m+4, n))
    abig = 0.0_wp
    abig(1:size(a,1),:) = a
    bbig = 0.0_wp
    bbig(1:size(b,1),:) = b
    cbig = 0.0_wp
    cbig(1:m,:) = c0
    call dense_matmul(abig(1:size(a,1),:), bbig(1:size(b,1),:), cbig(1:m,:), trans_a=trans_a, trans_b=trans_b, &
                      beta=1.0_wp)
    call compare(reshape(cbig(1:m,:), [m*n]), reshape(expected + c0, [m*n]), k, 'matmul, sections')
    if (any(cbig(m+1:,:) /= 0.0_wp)) call fail('matmul, sections: an element outside the section changed')

    end subroutine check_matmul

    subroutine check_large()
    !! a 600 by 600 product (2.9 MB) and matrix-vector product: their
    !! results against sums computed here, element by element
    integer, parameter :: n = 600
    real(wp), dimension(:,:), allocatable :: a, b, c
    real(wp), dimension(:), allocatable :: x, y
    integer :: ii, jj
    real(wp) :: s, err

    allocate(a(n,n), b(n,n), c(n,n), x(n), y(n))
    call random_number(a)
    call random_number(b)
    call random_number(x)
    call dense_matmul(a, b, c, trans_a=.true.)
    call dense_matvec(a, x, y)
    err = 0.0_wp
    do jj = 1, n, 37
        do ii = 1, n, 41
            s = dot_product(a(:,ii), b(:,jj))
            err = max(err, abs(c(ii,jj) - s)/s)
        end do
    end do
    do ii = 1, n
        s = 0.0_wp
        do jj = 1, n
            s = s + a(ii,jj)*x(jj)
        end do
        err = max(err, abs(y(ii) - s)/s)
    end do
    if (err > 1.0e-12_wp) call fail('the large products')

    end subroutine check_large

    subroutine compare(got, expected, k, what)
    !! check a result against the expected values, to a tolerance of the
    !! rounding error of sums of `k` products
    real(wp), dimension(:), intent(in) :: got      !! the result
    real(wp), dimension(:), intent(in) :: expected !! the expected values
    integer,                intent(in) :: k        !! the number of terms of each sum
    character(len=*),       intent(in) :: what     !! what is checked
    if (size(got) == 0) return
    if (.not. all(abs(got - expected) <= 1.0e-14_wp*(max(k,1) + 1)*max(1.0_wp, maxval(abs(expected))))) then
        call fail(what)
    end if
    end subroutine compare

    subroutine fail(what)
    !! report a failed check
    character(len=*), intent(in) :: what !! what failed
    n_fail = n_fail + 1
    print '(A)', 'FAILED: '//what
    end subroutine fail

end program test_dense_product
