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
!  This only decides the shift. The QP subproblems are solved by
!  [[sqpopt_qp_solver_module]] with the shifted Hessian (see
!  [[sqpopt_iterate]]), which can reuse the factorization for a direct step
!  (see [[sqpopt_qp_direct_module]]).

    module sqpopt_inertia_module

    use sqpopt_kinds,          only: wp => sqpopt_module_wp
    use sqpopt_types_module,   only: sqpopt_sparse_matrix
    use sqpopt_hessian_module, only: sqpopt_hessian_type
    use sqpopt_kkt_module,     only: sqpopt_kkt_type
    use sqpopt_symmetric_solver_module, only: sqpopt_has_mumps

    implicit none

    private

    public :: sqpopt_has_mumps ! (from [[sqpopt_symmetric_solver_module]])

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

    subroutine inertia_correct(me, kkt, hessian, jac, status, changed, ok, n_negative)

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
            ok = .true.
            exit
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

    end module sqpopt_inertia_module
!*******************************************************************************
