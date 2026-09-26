program test_independent_columns

    !! Unit tests of [[independent_columns]] (the rank-revealing `LUSOL`
    !! factorization that picks the sparse QP's initial working set):
    !! parallel, duplicated, zero, and nearly dependent columns, more columns
    !! than rows, and steering which column of a dependent group is kept
    !! (by scaling, and by unit columns always being kept first).

    use sqpopt_linalg_module, only: independent_columns
    use sqpopt_kinds,         only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: rel_tol = 1.0e-8_wp
    real(wp), parameter :: abs_tol = 1.0e-14_wp

    write(*,*) '----------------------------'
    write(*,*) 'test_independent_columns'
    write(*,*) '----------------------------'

    ! parallel, zero, and sum columns: c1=(1,2,0), c2=2*c1, c3=e3, c4=0, c5=c1+c3 -> rank 2
    call check('parallel/zero/sum', reshape([ 1.0_wp, 2.0_wp, 0.0_wp, &
                                              2.0_wp, 4.0_wp, 0.0_wp, &
                                              0.0_wp, 0.0_wp, 1.0_wp, &
                                              0.0_wp, 0.0_wp, 0.0_wp, &
                                              1.0_wp, 2.0_wp, 1.0_wp ], [3,5]), expected_rank=2, &
               must_drop=[4])

    ! the same, with c1 scaled down: the parallel column c2 is kept instead
    call check('scaled preference', reshape([ 1.0e-4_wp, 2.0e-4_wp, 0.0_wp, &
                                              2.0_wp,    4.0_wp,    0.0_wp, &
                                              0.0_wp,    0.0_wp,    1.0_wp ], [3,3]), expected_rank=2, &
               must_keep=[2,3], must_drop=[1])

    ! a general column parallel to a unit column: the unit column is kept
    call check('unit column first', reshape([ 3.0_wp, 0.0_wp, &
                                              1.0_wp, 0.0_wp, &
                                              0.0_wp, 1.0_wp ], [2,3]), expected_rank=2, &
               must_keep=[2,3], must_drop=[1])

    ! duplicated columns
    call check('duplicates', reshape([ 1.0_wp, 1.0_wp, 1.0_wp, &
                                       1.0_wp, 1.0_wp, 1.0_wp, &
                                       1.0_wp, 1.0_wp, 1.0_wp, &
                                       1.0_wp,-1.0_wp, 0.0_wp ], [3,4]), expected_rank=2, &
               must_keep=[4])

    ! nearly dependent (1e-12 relative) is dependent; 1e-5 relative is not
    call check('nearly dependent', reshape([ 1.0_wp, 1.0_wp, 1.0_wp, &
                                             1.0_wp, 1.0_wp, 1.0_wp+1.0e-12_wp ], [3,2]), expected_rank=1)
    call check('barely independent', reshape([ 1.0_wp, 1.0_wp, 1.0_wp, &
                                               1.0_wp, 1.0_wp, 1.0_wp+1.0e-5_wp ], [3,2]), expected_rank=2)

    ! more columns than rows, and a full-rank square matrix
    call check('wide', reshape([ 1.0_wp, 2.0_wp, 3.0_wp, &
                                 0.5_wp,-1.0_wp, 2.0_wp, &
                                 4.0_wp, 0.0_wp, 1.0_wp, &
                                 1.0_wp, 1.0_wp, 1.0_wp, &
                                -2.0_wp, 3.0_wp, 0.5_wp ], [3,5]), expected_rank=3)
    call check('square', reshape([ 4.0_wp, 1.0_wp, 0.0_wp, &
                                   1.0_wp, 4.0_wp, 1.0_wp, &
                                   0.0_wp, 1.0_wp, 4.0_wp ], [3,3]), expected_rank=3)

    print '(A)', 'test_independent_columns PASSED'

    contains

    subroutine check(name, a, expected_rank, must_keep, must_drop)
    !! call [[independent_columns]] on the dense matrix `a` (as COO triplets),
    !! and check the number of columns kept, that the kept columns are
    !! independent and span the others, and any required choices
    character(len=*),         intent(in) :: name
    real(wp), dimension(:,:), intent(in) :: a
    integer,                  intent(in) :: expected_rank
    integer, dimension(:),    intent(in), optional :: must_keep, must_drop
    integer,  dimension(:), allocatable :: ir, ic
    real(wp), dimension(:), allocatable :: vv
    logical,  dimension(size(a,2)) :: indep
    integer :: i, j, istat
    ir = [integer ::]; ic = [integer ::]; vv = [real(wp) ::]
    do j = 1, size(a,2)
        do i = 1, size(a,1)
            if (a(i,j) /= 0.0_wp) then
                ir = [ir, i]; ic = [ic, j]; vv = [vv, a(i,j)]
            end if
        end do
    end do
    call independent_columns(size(a,1), size(a,2), ir, ic, vv, rel_tol, abs_tol, indep, istat)
    print '(A,T22,A,*(L2))', name, ': independent =', indep
    if (istat /= 0) call fail(name, 'istat /= 0')
    if (count(indep) /= expected_rank) call fail(name, 'wrong number of independent columns')
    if (rank_of(a(:, pack([(j, j=1,size(a,2))], indep))) /= expected_rank) call fail(name, 'kept columns are dependent')
    if (present(must_keep)) then
        if (.not. all(indep(must_keep))) call fail(name, 'a preferred column was dropped')
    end if
    if (present(must_drop)) then
        if (any(indep(must_drop))) call fail(name, 'a dependent column was kept')
    end if
    end subroutine check

    integer function rank_of(b)
    !! numerical rank of a small dense matrix (Gaussian elimination with complete pivoting)
    real(wp), dimension(:,:), intent(in) :: b
    real(wp), dimension(size(b,1),size(b,2)) :: w
    real(wp) :: tol
    integer :: k, loc(2), i
    w = b
    tol = 1.0e-10_wp*max(1.0_wp, maxval(abs(b)))
    rank_of = 0
    do k = 1, min(size(w,1), size(w,2))
        loc = maxloc(abs(w(k:,k:))) + k - 1
        if (abs(w(loc(1),loc(2))) <= tol) exit
        w([k,loc(1)],:) = w([loc(1),k],:)
        w(:,[k,loc(2)]) = w(:,[loc(2),k])
        do i = k+1, size(w,1)
            w(i,k:) = w(i,k:) - w(i,k)/w(k,k)*w(k,k:)
        end do
        rank_of = rank_of + 1
    end do
    end function rank_of

    subroutine fail(name, msg)
    character(len=*), intent(in) :: name, msg
    print '(A)', 'test_independent_columns FAILED: '//name//': '//msg
    error stop 1
    end subroutine fail

end program test_independent_columns
