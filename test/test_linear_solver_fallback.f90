program test_linear_solver_fallback

    !! Test of QDLDL's fallback on a solver that pivots
    !! (`options%linear_solver_fallback`, see [[sqpopt_symmetric_solver_module]]).
    !!
    !! QDLDL doesn't pivot, so it can't tell the inertia of
    !! \( \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix} \) (its first pivot is
    !! zero, but the matrix is nonsingular, with one negative eigenvalue). With
    !! the fallback, the matrix is factored again by MUMPS (in a build with
    !! it) or a dense solver, which then tells the inertia, solves, and is kept
    !! for the next factorizations. Without it, QDLDL stays. Above
    !! `sqpopt_fallback_max_order`, the fallback is MUMPS only.

    use sqpopt_symmetric_solver_module, only: sqpopt_symmetric_solver_type, sqpopt_has_mumps, &
                                              sqpopt_linear_solver_qdldl, sqpopt_linear_solver_mumps, &
                                              sqpopt_fallback_max_order
    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    write(*,*) '----------------------------'
    write(*,*) 'test_linear_solver_fallback'
    write(*,*) '----------------------------'

    call test_small(.false.)
    call test_small(.true.)
    call test_large()

    print '(A)', 'test_linear_solver_fallback PASSED'

    contains

    subroutine test_small(fallback)
    !! the 2x2 matrix of the program's documentation, with or without the fallback

    logical, intent(in) :: fallback !! whether QDLDL falls back on a solver that pivots

    type(sqpopt_symmetric_solver_type) :: solver
    real(wp), dimension(2) :: b
    logical :: ok

    call solver%initialize(2, [1, 2, 2], [1, 1, 2], ok, solver=sqpopt_linear_solver_qdldl, fallback=fallback)
    if (.not. ok) error stop 'test_linear_solver_fallback FAILED: QDLDL could not be started'
    call solver%factor([0.0_wp, 1.0_wp, 0.0_wp], ok)
    print '(A,L1,A,I0,A,L1,A,I0)', 'fallback = ', fallback, ': backend = ', solver%backend, &
        ', inertia known = ', solver%inertia_known, ', negative eigenvalues = ', solver%n_negative

    if (fallback) then
        if (.not. ok) error stop 'test_linear_solver_fallback FAILED: not factored'
        if (.not. solver%switched .or. solver%backend == sqpopt_linear_solver_qdldl) then
            error stop 'test_linear_solver_fallback FAILED: did not switch'
        end if
        if (.not. solver%inertia_known .or. solver%n_negative /= 1 .or. solver%n_null /= 0) then
            error stop 'test_linear_solver_fallback FAILED: wrong inertia after the switch'
        end if
        ! the solve: [0 1; 1 0] x = [2; 3] gives x = [3; 2]
        b = [2.0_wp, 3.0_wp]
        call solver%solve(b, ok)
        if (.not. ok .or. maxval(abs(b - [3.0_wp, 2.0_wp])) > 1.0e-12_wp) error stop 'test_linear_solver_fallback FAILED: solve'
        ! the next factorization uses the same solver, and is counted
        call solver%factor([2.0_wp, 1.0_wp, 2.0_wp], ok)
        if (.not. ok .or. solver%backend == sqpopt_linear_solver_qdldl .or. solver%n_negative /= 0 .or. &
            solver%n_factor /= 3) error stop 'test_linear_solver_fallback FAILED: refactorization'
    else
        if (solver%switched .or. solver%backend /= sqpopt_linear_solver_qdldl) then
            error stop 'test_linear_solver_fallback FAILED: switched without the fallback'
        end if
        if (ok .and. solver%inertia_known) error stop 'test_linear_solver_fallback FAILED: QDLDL told the inertia'
    end if
    call solver%destroy()
    if (solver%switched .or. solver%fallback) error stop 'test_linear_solver_fallback FAILED: destroy'

    print '(A,L1,A)', 'test_linear_solver_fallback [fallback = ', fallback, '] PASSED'

    end subroutine test_small

    subroutine test_large()
    !! the same 2x2 block, with the identity of order `sqpopt_fallback_max_order`
    !! after it: the fallback is MUMPS, or nothing without it

    type(sqpopt_symmetric_solver_type) :: solver
    integer,  dimension(:), allocatable :: irow, icol
    real(wp), dimension(:), allocatable :: val
    integer :: n, k
    logical :: ok

    n = sqpopt_fallback_max_order + 2
    allocate(irow(n + 1), icol(n + 1), val(n + 1))
    irow(1:3) = [1, 2, 2]
    icol(1:3) = [1, 1, 2]
    val(1:3)  = [0.0_wp, 1.0_wp, 0.0_wp]
    do k = 3, n
        irow(k + 1) = k
        icol(k + 1) = k
        val(k + 1)  = 1.0_wp
    end do

    call solver%initialize(n, irow, icol, ok, solver=sqpopt_linear_solver_qdldl, fallback=.true.)
    if (.not. ok) error stop 'test_linear_solver_fallback FAILED: QDLDL could not be started (large)'
    call solver%factor(val, ok)
    print '(A,I0,A,I0,A,L1)', 'order ', n, ': backend = ', solver%backend, ', inertia known = ', solver%inertia_known
    if (.not. solver%switched) error stop 'test_linear_solver_fallback FAILED: the switch was not tried (large)'
    if (sqpopt_has_mumps) then
        if (.not. ok .or. solver%backend /= sqpopt_linear_solver_mumps .or. solver%n_negative /= 1) then
            error stop 'test_linear_solver_fallback FAILED: no switch to MUMPS (large)'
        end if
    else
        if (solver%backend /= sqpopt_linear_solver_qdldl) then
            error stop 'test_linear_solver_fallback FAILED: a dense solver above its order (large)'
        end if
    end if
    call solver%destroy()

    print '(A)', 'test_linear_solver_fallback [large] PASSED'

    end subroutine test_large

end program test_linear_solver_fallback
