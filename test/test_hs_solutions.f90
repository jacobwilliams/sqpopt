program test_hs_solutions

    !! Validation of the Hock-Schittkowski collection's reference solutions
    !! (see [[hs_problems_module]]): for each of the 305 problems, the
    !! objective and constraints are evaluated at the recorded optimal point
    !! `x_star`, which must
    !!
    !! * reproduce the recorded optimal value `f_star`, to within `rel_tol`
    !!   (relative to \( \max(1,|f^*|) \), as in `test_hs_suite`), and
    !! * be feasible, up to the rounding of its recorded digits (many points
    !!   are recorded to only 5-8 significant digits): a variable bound must
    !!   hold to within `feas_tol` \( \max(1,|x_j|) \), and a constraint to
    !!   within `feas_tol` \( \max(1,|c_i|) + \sigma_i \), where
    !!   \( \sigma_i = \sum_j |c_i(x^* + h_j e_j) - c_i(x^*)| \),
    !!   \( h_j = \) `digits` \( |x^*_j| \), is how much the constraint can
    !!   change when each coordinate is rounded to 4 significant digits (the
    !!   coarsest recording in the collection, apart from TP377)
    !!   (by finite differences: some problems have no analytic derivatives).
    !!   A constraint made of large, nearly cancelling terms (e.g. TP359) is
    !!   very sensitive to that rounding. (The `rel.viol` column is the
    !!   violation relative to its tolerance, times `feas_tol`.)
    !!
    !! This guards the problem definitions (e.g. an edit to an objective or a
    !! constraint in `test/schittkowski_problems.f90`), and documents the
    !! collection's own inconsistencies, listed in `known_inconsistent` with
    !! the reason for each. The test fails if any other problem is
    !! inconsistent, and reports listed problems that no longer are (update
    !! the list).

    use hs_problems_module, only: hs_problem, hs_setup, hs_f, hs_c, hs_problem_ids, hs_n_problems
    use sqpopt_kinds, only: dp => sqpopt_module_wp

    implicit none

    real(dp), parameter :: rel_tol  = 1.0e-4_dp  !! objective tolerance (as `test_hs_suite`'s "solved" test)
    real(dp), parameter :: feas_tol = 1.0e-4_dp  !! relative feasibility tolerance
    real(dp), parameter :: digits   = 1.0e-4_dp  !! relative rounding of the recorded `x_star` (4 significant digits)

    !> problems whose recorded solution is inconsistent, in the original collection. (Two others, TP216 and
    !! TP340, had errors in their data, which are fixed in `test/schittkowski_problems.f90`.)
    !!
    !! * TP112: `x_star` gives f = -47.7076, not f* = -47.7611 (solvers reach f*: `x_star` is wrong)
    !! * TP368: several solutions (`NEX=3`), but the others are commented out; the remaining `x_star` gives
    !!   f = -1, lower than f* = -0.74998 (solvers reach f*, a local minimum)
    !! * TP377: `x_star` is recorded to only 2-3 digits (e.g. 10.0, 9.5), too coarse to be feasible
    !! * TP380: `x_star`'s digits reproduce f* only to 1.2e-4
    !! * TP381, TP382: `x_star` is infeasible, and doesn't reproduce f* (solvers reach f*: `x_star` is wrong)
    integer, dimension(*), parameter :: known_inconsistent = [112, 368, 377, 380, 381, 382]

    type(hs_problem) :: p
    integer  :: k, id, n_new, n_fixed
    real(dp) :: f, rel, viol
    real(dp), dimension(:), allocatable :: c, ch, sigma, xh
    integer  :: j
    logical  :: ok, listed

    write(*,*) '----------------------------'
    write(*,*) 'test_hs_solutions'
    write(*,*) '----------------------------'
    write(*,'(A)') '   TP  exact                f(x*)                   f*   rel.err  rel.viol'

    n_new   = 0
    n_fixed = 0
    do k = 1, hs_n_problems
        id = hs_problem_ids(k)
        call hs_setup(id, p)
        allocate(c(p%m), ch(p%m), sigma(p%m))
        call hs_f(id, p%x_star, f)
        sigma = 0.0_dp
        if (p%m > 0) then
            call hs_c(id, p%x_star, c)
            ! the constraints' sensitivity to rounding x_star:
            do j = 1, p%n
                if (p%x_star(j) == 0.0_dp) cycle
                xh = p%x_star
                xh(j) = xh(j) + digits*abs(xh(j))
                call hs_c(id, xh, ch)
                sigma = sigma + abs(ch - c)
            end do
        end if

        rel  = abs(f - p%f_star)/max(1.0_dp, abs(p%f_star))
        viol = maxval([0.0_dp, (max(p%x_lb - p%x_star, 0.0_dp) + max(p%x_star - p%x_ub, 0.0_dp)) &
                                / max(1.0_dp, abs(p%x_star))])
        if (p%m > 0) viol = max(viol, maxval((max(p%c_lb - c, 0.0_dp) + max(c - p%c_ub, 0.0_dp)) &
                                             / (max(1.0_dp, abs(c)) + sigma/feas_tol)))
        ok     = rel <= rel_tol .and. viol <= feas_tol
        listed = any(known_inconsistent == id)

        if (.not. ok .or. listed) then
            write(*,'(I5,L7,2ES21.12,2ES10.2,2X,A)') id, p%exact, f, p%f_star, rel, viol, trim(label(ok, listed))
        end if
        if (.not. ok .and. .not. listed) n_new = n_new + 1
        if (ok .and. listed) n_fixed = n_fixed + 1
        deallocate(c, ch, sigma)
    end do

    write(*,'(I0,A,I0,A)') hs_n_problems - size(known_inconsistent) + n_fixed - n_new, ' of ', hs_n_problems, &
        ' reference solutions are consistent'
    if (n_fixed > 0) write(*,'(I0,A)') n_fixed, ' problem(s) in known_inconsistent are now consistent: update the list'
    if (n_new > 0) then
        write(*,'(I0,A)') n_new, ' problem(s) not in known_inconsistent are inconsistent'
        error stop 'test_hs_solutions FAILED'
    end if
    write(*,*) 'test_hs_solutions PASSED'

contains

    pure function label(ok, listed) result(s)
    !! how a problem shown in the table compares with `known_inconsistent`
    logical, intent(in) :: ok     !! whether the problem's reference solution is consistent
    logical, intent(in) :: listed !! whether the problem is in `known_inconsistent`
    character(len=:), allocatable :: s
    if (ok) then
        s = '<-- now consistent'
    else if (listed) then
        s = '(known)'
    else
        s = '<-- inconsistent'
    end if
    end function label

end program test_hs_solutions
