program test_convergence

    !! Unit tests of [[check_convergence]] on hand-constructed points:
    !!
    !! * a constraint at its lower bound needs `lambda>=0` (for the Lagrangian
    !!   `f - lambda^T c`); with the wrong sign the point is not a KKT point,
    !!   even though `g - J^T lambda = 0`;
    !! * an inactive constraint must have a zero multiplier (complementarity);
    !! * a feasible point where nothing is changing any more is reported as
    !!   `sqpopt_stalled`, not `sqpopt_success`;
    !! * an infeasible point that is stationary for the constraint violation
    !!   is reported as `sqpopt_infeasible`.

    use sqpopt_convergence_module, only: check_convergence
    use sqpopt_types_module,       only: sqpopt_sparse_matrix, sqpopt_success, sqpopt_stalled, sqpopt_infeasible
    use sqpopt_kinds,              only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: ktol = 1.0e-6_wp, ctol = 1.0e-8_wp
    real(wp), parameter :: x_lb(1) = [-10.0_wp], x_ub(1) = [10.0_wp]

    type(sqpopt_sparse_matrix) :: jac1, jac2
    logical :: converged
    integer :: istat

    write(*,*) '----------------------------'
    write(*,*) 'test_convergence'
    write(*,*) '----------------------------'

    ! one variable, one constraint c = x in [0, 5]:
    jac1 = sqpopt_sparse_matrix(nrows=1, ncols=1, nnz=1, irow=[1], icol=[1], val=[1.0_wp])

    ! min x at x=0 (c at its lower bound): g=1, lambda=1 is correct:
    call check_convergence([0.0_wp], [1.0_wp], jac1, [0.0_wp], x_lb, x_ub, [0.0_wp], [5.0_wp], &
                           [1.0_wp], ktol, ctol, converged, istat)
    call check('correct-sign multiplier at a lower bound', converged .and. istat == sqpopt_success)

    ! min -x at x=0: g=-1, lambda=-1 satisfies g = J^T lambda, but the sign is
    ! wrong at a lower bound (moving up decreases f), so not converged:
    call check_convergence([0.0_wp], [-1.0_wp], jac1, [0.0_wp], x_lb, x_ub, [0.0_wp], [5.0_wp], &
                           [-1.0_wp], ktol, ctol, converged, istat)
    call check('wrong-sign multiplier at a lower bound', .not. converged)

    ! min -x at x=5 (c at its upper bound): lambda=-1 is correct there:
    call check_convergence([5.0_wp], [-1.0_wp], jac1, [5.0_wp], x_lb, x_ub, [0.0_wp], [5.0_wp], &
                           [-1.0_wp], ktol, ctol, converged, istat)
    call check('correct-sign multiplier at an upper bound', converged .and. istat == sqpopt_success)

    ! inactive constraint (x=2 strictly inside [0,5]) with a nonzero multiplier:
    call check_convergence([2.0_wp], [1.0_wp], jac1, [2.0_wp], x_lb, x_ub, [0.0_wp], [5.0_wp], &
                           [1.0_wp], ktol, ctol, converged, istat)
    call check('nonzero multiplier on an inactive constraint', .not. converged)

    ! feasible, KKT not satisfied, but f and x unchanged from the previous iterate:
    call check_convergence([2.0_wp], [1.0_wp], jac1, [2.0_wp], x_lb, x_ub, [0.0_wp], [5.0_wp], &
                           [0.0_wp], ktol, ctol, converged, istat, &
                           f=2.0_wp, f_prev=2.0_wp, x_prev=[2.0_wp], ftol=1.0e-8_wp, xtol=1.0e-8_wp)
    call check('stalled progress', converged .and. istat == sqpopt_stalled)

    ! infeasible: c1 = x in [0,1] and c2 = x in [2,3]; at x=1.5 the (l2)
    ! violation is stationary, so no progress toward feasibility is possible:
    jac2 = sqpopt_sparse_matrix(nrows=2, ncols=1, nnz=2, irow=[1,2], icol=[1,1], val=[1.0_wp, 1.0_wp])
    call check_convergence([1.5_wp], [3.0_wp], jac2, [1.5_wp, 1.5_wp], x_lb, x_ub, [0.0_wp, 2.0_wp], [1.0_wp, 3.0_wp], &
                           [0.0_wp, 0.0_wp], ktol, ctol, converged, istat, &
                           f=2.25_wp, f_prev=2.0_wp, x_prev=[1.4_wp], ftol=1.0e-8_wp, xtol=1.0e-8_wp)
    call check('stationary infeasible point', converged .and. istat == sqpopt_infeasible)

    ! same constraints at x=1.2: still infeasible, but not stationary for the violation:
    call check_convergence([1.2_wp], [2.4_wp], jac2, [1.2_wp, 1.2_wp], x_lb, x_ub, [0.0_wp, 2.0_wp], [1.0_wp, 3.0_wp], &
                           [0.0_wp, 0.0_wp], ktol, ctol, converged, istat, &
                           f=1.44_wp, f_prev=1.0_wp, x_prev=[1.0_wp], ftol=1.0e-8_wp, xtol=1.0e-8_wp)
    call check('non-stationary infeasible point', .not. converged)

    print '(A)', 'test_convergence PASSED'

    contains

    subroutine check(label, ok)
    !! print the case, and stop with a failure message unless `ok`
    character(len=*), intent(in) :: label !! the case, for the message
    logical,          intent(in) :: ok    !! the condition that must hold
    print '(A,A,L1)', label, ': ', ok
    if (.not. ok) error stop 'test_convergence FAILED: '//label
    end subroutine check

end program test_convergence
