!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Inertia control for the Hessians that can be indefinite, the exact one
!  and SR1 (`options%inertia_control`): the shift \( \delta \) of
!  \( H + \delta I \) is found from a sparse factorization of the KKT
!  matrix of the QP subproblem's working set (see [[sqpopt_kkt_module]]),
!  by the sparse solver of `options%linear_solver` (MUMPS or QDLDL; with
!  the exact Hessian, MUMPS is the better one, see
!  [[sqpopt_symmetric_solver_module]]).
!
!  The QP subproblem's step is only a minimizer if the Hessian is positive
!  definite on the null space of the working set's constraints,
!  \( Z^THZ \succ 0 \). The factorization of the KKT matrix tells how many
!  directions of negative curvature that null space has, so the smallest
!  shift that leaves none is found by refactoring ([[inertia_correct]]),
!  much as in IPOPT: first the current shift, then `shift_min` times the
!  size of the Hessian (or a third of the last shift that was needed, if
!  that is larger), then 8 times as much at each attempt, up to
!  `shift_max`. (IPOPT multiplies by 100 until a shift has been needed
!  once. On the Hock-Schittkowski problems that was no better: 273 solved
!  and 3 failed, against 274 and 2, with 1% more function evaluations.
!  Factors of 2, 4, and 16 gave 271 or 272 solved.)
!
!  **The curvature test** (`options%inertia_test`, see
!  [[sqpopt_kkt_module]]). Instead of the inertia, the shift can be decided
!  by the inertia-free test of Chiang and Zavala (as HiOp's
!  `fact_acceptor = inertia_free`): the factorization is accepted, and the
!  QP solved, and the shift is raised only if the QP's step `p` has too
!  little curvature, \( p^T (H + \delta I) p < \kappa\, p^T p \) (see
!  [[enough_curvature]]); the QP is then re-solved. With
!  `sqpopt_inertia_test_curvature` that is always the test; with
!  `sqpopt_inertia_test_unknown`, only where the linear solver can't tell
!  the inertia. Before the QP solve, and in the direct QP method, there is
!  no step to test, and the face is taken as convex. On the
!  Hock-Schittkowski problems (QDLDL, release build), with the exact
!  Hessian: always testing the curvature solved the most problems (277, 1
!  failed, 12,599 evaluations, against the inertia's 273, 2, and 9,431,
!  with the fallback of [[sqpopt_symmetric_solver_module]]; the same with
!  MUMPS), but with the direct QP method it was worse (271 and 2 failed,
!  against 277 and 1), and with SR1 much worse (268 and 6 failed, against
!  275 and 1). Testing only where the inertia is unknown made no
!  difference with the fallback (the inertia is then always known), and
!  without it was 272 solved and 3 failed with 9,508 evaluations. So the
!  default is the inertia (`sqpopt_inertia_test_count`).
!
!  This only decides the shift. The QP subproblems are solved by
!  [[sqpopt_qp_solver_module]] with the shifted Hessian (see
!  [[sqpopt_iterate]]), which can reuse the factorization for a direct step
!  (see [[sqpopt_qp_direct_module]]).

    module sqpopt_inertia_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_kkt_module,     only: sqpopt_kkt_type

    implicit none

    private


    type, public :: sqpopt_inertia_type
        !! the inertia control of one solve: the state of the shift's search
        !! (see the module documentation). The KKT matrix it factors is a
        !! [[sqpopt_kkt_type]], which the caller owns.

        logical  :: enabled = .false.     !! whether inertia control is in use (set by the solver, and unset by
                                          !! [[inertia_correct]] if the KKT matrix can't be factored)
        real(wp) :: shift_last = 0.0_wp   !! the last shift that a correction ended with (`0` if none has
                                          !! been needed yet)

        contains

        procedure, public :: correct => inertia_correct
        procedure, public :: raise   => inertia_raise

    end type sqpopt_inertia_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  increase the Hessian's shift (`hessian%shift`), if necessary, until
!  \( H + \delta I \) has no negative curvature on the null space of the
!  working set `status` (see the module documentation for the shifts that
!  are tried). The shift is never decreased: the search starts from its
!  current value. On success, `kkt` holds the factorization for the final
!  shift.
!
!  `ok` is true if no negative curvature is left. It is false if inertia
!  control isn't enabled, if a factorization failed (inertia control is then
!  disabled for the rest of the solve), or if the largest shift,
!  `shift_max`, isn't enough.
!
!  Where the KKT matrix's inertia isn't used (`kkt%curvature_needed`, see
!  the module documentation), `step` is the QP's step for this working set
!  and the shift on entry: if it has too little curvature, the shift is
!  raised once (the caller then re-solves the QP, and calls this again
!  with the new step). Without `step`, the face is taken as convex.

    subroutine inertia_correct(me, kkt, hessian, jac, status, changed, ok, n_negative, step)

    class(sqpopt_inertia_type), intent(inout) :: me
    type(sqpopt_kkt_type),      intent(inout) :: kkt     !! the KKT matrix (see [[sqpopt_kkt_module]])
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! the Hessian (its `shift` is updated)
    type(sqpopt_sparse_matrix), intent(in)    :: jac     !! the constraint Jacobian `dimension(m,n)`
    integer, dimension(:),      intent(in)    :: status  !! the working set: nonzero for each general row, then each
                                                         !! variable bound, that is in it `dimension(m+n)`
    logical,                    intent(out)   :: changed !! whether the shift was increased
    logical,                    intent(out)   :: ok      !! whether no negative curvature is left (see above)
    integer, optional,          intent(out)   :: n_negative !! the number of directions of negative curvature
                                                            !! found with the shift on entry (`-1` if nothing
                                                            !! was factored)
    real(wp), dimension(:), optional, intent(in) :: step !! the QP's step for the working set `status` and the
                                                         !! shift on entry, for the curvature test (see above;
                                                         !! only used where `kkt%curvature_needed`)
                                                         !! `dimension(n)`

    real(wp) :: shift_max
    logical  :: first, factored

    changed = .false.
    ok      = .false.
    if (present(n_negative)) n_negative = -1
    if (.not. me%enabled) return

    shift_max = hessian%shift_max*hessian%magnitude()
    first = .true.
    do
        call kkt%factor(hessian, jac, status, factored)
        if (.not. factored) then
            ! (fall back on the matrix-free tests from here on)
            me%enabled = .false.
            return
        end if
        if (first .and. present(n_negative)) n_negative = kkt%n_negative
        if (kkt%n_negative == 0) then
            ! (where the inertia isn't used, the step's curvature decides, if
            ! there is a step for this shift: see the module documentation)
            ok = .true.
            if (.not. (kkt%curvature_needed .and. present(step) .and. first)) exit
            if (enough_curvature(hessian, step)) exit
            ok = .false.
            if (present(n_negative)) n_negative = 1
        end if
        if (hessian%shift >= shift_max) exit
        call me%raise(hessian)
        changed = .true.
        first   = .false.
    end do
    if (ok .and. changed) me%shift_last = hessian%shift

    end subroutine inertia_correct
!*******************************************************************************

!*******************************************************************************
!>
!  increase the Hessian's shift to the next one to try (see the module
!  documentation): 8 times the current one, and at least `shift_min` times
!  the size of the Hessian, or a third of the last shift that was needed;
!  at most `shift_max` times the size of the Hessian.

    subroutine inertia_raise(me, hessian)

    class(sqpopt_inertia_type), intent(in)    :: me
    type(sqpopt_hessian_type),  intent(inout) :: hessian !! the Hessian (its `shift` is updated)

    real(wp) :: shift

    shift = hessian%shift_min*hessian%magnitude()
    if (me%shift_last > 0.0_wp) shift = max(shift, me%shift_last/3.0_wp)
    shift = max(shift, 8.0_wp*hessian%shift)
    hessian%shift = min(shift, hessian%shift_max*hessian%magnitude())

    end subroutine inertia_raise
!*******************************************************************************

!*******************************************************************************
!>
!  the curvature test of Chiang and Zavala (as in HiOp): whether the step
!  `p` has enough curvature along the shifted Hessian,
!  \( p^T (H + \delta I) p \ge \kappa\, p^T p \), with \( \kappa \) =
!  `curvature_tol` times the size of the Hessian. A zero step passes.

    function enough_curvature(hessian, p) result(enough)

    type(sqpopt_hessian_type), intent(inout) :: hessian !! the Hessian, with its current shift
    real(wp), dimension(:),    intent(in)    :: p       !! the step `dimension(n)`
    logical :: enough

    real(wp), parameter :: curvature_tol = 1.0e-11_wp !! \( \kappa \), relative to the size of the Hessian (HiOp's
                                                      !! `neg_curv_test_fact`)
    real(wp), dimension(:), allocatable :: hp

    allocate(hp(size(p)))
    call hessian%hv_product(p, hp)
    enough = dot_product(p, hp) >= curvature_tol*hessian%magnitude()*dot_product(p, p)

    end function enough_curvature
!*******************************************************************************

    end module sqpopt_inertia_module
!*******************************************************************************
