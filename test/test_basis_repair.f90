program test_basis_repair

    !! Test of the sparse QP's basis repair and basis-selection settings
    !! (`qp_solver%sparse_qp%basis_repair`, `basis_rel_tol`, and `basis_eq_slack_weight`, the
    !! `repair`, `rel_tol`, and `eq_slack_weight` arguments of [[choose_basis]], and the
    !! `repair_rows` argument of [[independent_columns]]):
    !!
    !! * [[independent_columns]] on a matrix with a column whose pivot is too small for
    !!   `rel_tol` (though the matrix has full rank), and a row without any pivot: `repair_rows`
    !!   must be those two rows, and the kept columns with those rows' unit columns must make
    !!   a square nonsingular matrix.
    !! * a QP with two equality rows, \( p_1 + p_2 = 1 \) and \( p_1 + 1.001 p_2 = 1.001 \), whose
    !!   basis LUSOL rejects a column of with `basis_rel_tol = 1e-2`: without `basis_repair`, the
    !!   basis method gives up and the QP is solved by the LSQR method (`basis_used` is false);
    !!   with it, the basis is completed with a slack and the basis method solves it. Both must
    !!   give the QP's solution, as must the defaults (whose basis is square).
    !! * the same QP with a third, inactive, row (whose free slack is basic).
    !! * Hock-Schittkowski problem 71, solved with the sparse QP and the repair on, must converge
    !!   to its known solution (its bases are square, so the repair must not change the solve).
    !!
    !! (The other change, `independent_columns` retrying `lu1fac` up to 8 times with more
    !! workspace, is not exercised here: no small basis needs more than its first workspace.)

    use sqpopt_module,                    only: sqpopt_type
    use sqpopt_problem_module,            only: sqpopt_problem_type
    use sqpopt_options_module,            only: sqpopt_options_type
    use sqpopt_hessian_module,            only: sqpopt_hessian_type
    use sqpopt_qp_solver_module,          only: sqpopt_qp_solver_type, sqpopt_qp_reduced_hessian
    use sqpopt_qp_reduced_hessian_module, only: sqpopt_reduced_hessian_qp_type
    use sqpopt_linalg_module,             only: independent_columns
    use sqpopt_types_module,              only: sqpopt_success, sqpopt_sparse_matrix
    use sqpopt_kinds,                     only: wp => sqpopt_module_wp

    implicit none

    real(wp), parameter :: big = 1.0e20_wp !! sentinel value used for "unbounded" sides

    write(*,*) '----------------------------'
    write(*,*) 'test_basis_repair'
    write(*,*) '----------------------------'

    call test_repair_rows()
    call test_qp(with_free_row=.false.)
    call test_qp(with_free_row=.true.)
    call test_hs71()
    print '(A)', 'test_basis_repair PASSED'

    contains

    subroutine test_repair_rows()
    !! [[independent_columns]]'s `repair_rows`: the 3 x 3 matrix with columns (1,1,0), (1,1.001,0),
    !! and (0,0,0), with `rel_tol = 1e-2`: one of the first two columns is rejected (its pivot is about
    !! `1e-3` of its column of `U`), and row 3 has no pivot
    real(wp), parameter :: a(3,3) = reshape([1.0_wp, 1.0_wp, 0.0_wp, &
                                             1.0_wp, 1.001_wp, 0.0_wp, &
                                             0.0_wp, 0.0_wp, 0.0_wp], [3,3])
    logical  :: indep(3), rrows(3)
    real(wp) :: b(3,3)
    integer  :: istat, j, k, i

    call independent_columns(3, 3, [1,2,1,2], [1,1,2,2], [1.0_wp, 1.0_wp, 1.0_wp, 1.001_wp], &
                             1.0e-2_wp, 1.0e-14_wp, indep, istat, repair_rows=rrows)
    print '(A,3L2,A,3L2)', 'repair_rows: independent =', indep, ', repair_rows =', rrows
    if (istat /= 0) call fail('repair_rows', 'istat /= 0')
    if (count(indep) /= 1 .or. indep(3)) call fail('repair_rows', 'expected one of the first two columns')
    if (.not. rrows(3)) call fail('repair_rows', 'row 3 (no pivot) must be repaired')
    if (count(rrows) /= 3 - count(indep)) call fail('repair_rows', 'wrong number of rows to repair')

    ! the kept columns and the repaired rows' unit columns: a square nonsingular matrix
    b = 0.0_wp
    k = 0
    do j = 1, 3
        if (indep(j)) then
            k = k + 1
            b(:,k) = a(:,j)
        end if
    end do
    do i = 1, 3
        if (rrows(i)) then
            k = k + 1
            b(i,k) = 1.0_wp
        end if
    end do
    if (abs(det3(b)) < 1.0e-6_wp) call fail('repair_rows', 'the repaired basis is singular')
    print '(A)', 'test_repair_rows PASSED'

    end subroutine test_repair_rows

    subroutine test_qp(with_free_row)
    !! the QP: minimize \( (p_1^2+p_2^2+p_3^2)/2 - 2 p_3 \) subject to \( p_1 + p_2 = 1 \) and
    !! \( p_1 + 1.001 p_2 = 1.001 \) (and, with `with_free_row`, \( -10 \le p_3 \le 10 \) as a third,
    !! inactive, row). The rows fix \( p_1 = 0 \), \( p_2 = 1 \), and \( p_3 = 2 \).
    logical, intent(in) :: with_free_row !! add the inactive third row
    type(sqpopt_reduced_hessian_qp_type) :: qp
    character(len=:), allocatable :: name

    name = merge('qp + free row', 'qp           ', with_free_row)

    ! the defaults: the basis is square
    call run_qp(name//', defaults', qp, with_free_row, expect_basis=.true.)

    ! a tolerance that rejects a pivot: the basis method gives up without the repair
    qp%basis_rel_tol = 9.0e-1_wp
    call run_qp(name//', rel_tol 1e-2', qp, with_free_row, expect_basis=.false.)

    ! and the repair completes the basis with a slack
    qp%basis_repair = .true.
    call run_qp(name//', rel_tol 1e-2 + repair', qp, with_free_row, expect_basis=.true.)

    ! the equality slacks' weight doesn't change the solution
    qp%basis_eq_slack_weight = 1.0e-3_wp
    call run_qp(name//', + eq_slack_weight 1e-3', qp, with_free_row, expect_basis=.true.)

    end subroutine test_qp

    subroutine run_qp(label, qp, with_free_row, expect_basis)
    !! solve the QP of [[test_qp]] with `qp`'s settings, and check its solution and which method solved it
    character(len=*),                     intent(in)    :: label         !! the case, for the output
    type(sqpopt_reduced_hessian_qp_type), intent(inout) :: qp            !! the solver, with its settings
    logical,                              intent(in)    :: with_free_row !! add the inactive third row
    logical,                              intent(in)    :: expect_basis  !! whether the basis method must solve it
    real(wp), parameter :: pexpect(3) = [0.0_wp, 1.0_wp, 2.0_wp]
    type(sqpopt_hessian_type)  :: hessian
    type(sqpopt_sparse_matrix) :: jac
    real(wp), dimension(:), allocatable :: c, c_lb, c_ub, lambda
    real(wp) :: x(3), g(3), x_lb(3), x_ub(3), p(3)
    integer  :: istat, m

    m = merge(3, 2, with_free_row)
    call hessian%initialize(n=3, max_history=5)
    jac%nrows = m; jac%ncols = 3
    if (with_free_row) then
        jac%nnz = 5
        jac%irow = [1,1,2,2,3]; jac%icol = [1,2,1,2,3]; jac%val = [1.0_wp, 1.0_wp, 1.0_wp, 1.001_wp, 1.0_wp]
        c_lb = [1.0_wp, 1.001_wp, -10.0_wp]; c_ub = [1.0_wp, 1.001_wp, 10.0_wp]
    else
        jac%nnz = 4
        jac%irow = [1,1,2,2]; jac%icol = [1,2,1,2]; jac%val = [1.0_wp, 1.0_wp, 1.0_wp, 1.001_wp]
        c_lb = [1.0_wp, 1.001_wp]; c_ub = c_lb
    end if
    allocate(c(m), lambda(m))
    c = 0.0_wp
    x = 0.0_wp
    g = [0.0_wp, 0.0_wp, -2.0_wp]
    x_lb = -big
    x_ub = big
    if (allocated(qp%warm_status)) deallocate(qp%warm_status)

    call qp%solve(hessian, jac, x, g, c, x_lb, x_ub, c_lb, c_ub, p, lambda, istat)

    print '(A,T44,A,3F12.8,A,I0,A,L2)', label, ': p =', p, ', istat ', istat, ', basis_used', qp%basis_used
    if (istat /= sqpopt_success) call fail(label, 'istat /= sqpopt_success')
    if (maxval(abs(p - pexpect)) > 1.0e-6_wp) call fail(label, 'wrong p')
    if (qp%basis_used .neqv. expect_basis) call fail(label, 'solved by the wrong method')

    end subroutine run_qp

    subroutine test_hs71()
    !! HS71 with the sparse QP and the basis repair on
    type(sqpopt_problem_type)   :: problem
    type(sqpopt_options_type)   :: options
    type(sqpopt_qp_solver_type) :: qp_solver
    type(sqpopt_type)           :: solver
    real(wp) :: x(4), lambda(2)
    integer  :: istat
    real(wp), parameter :: x_star(4) = [1.0_wp, 4.74299963_wp, 3.82114998_wp, 1.37940829_wp]

    call problem%set_problem_size(n=4, m=2)
    call problem%set_bounds(x_lb=[1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp], x_ub=[5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp], &
                            c_lb=[40.0_wp, 25.0_wp], c_ub=[40.0_wp, big])
    call problem%set_jacobian_sparsity(nnz=8, irow=[1,1,1,1,2,2,2,2], icol=[1,2,3,4,1,2,3,4])
    call problem%set_functions(fc=fc, gjac=gjac)
    options%qp_solver_mode = sqpopt_qp_reduced_hessian
    qp_solver%sparse_qp%basis_repair  = .true.
    call solver%initialize(problem=problem, options=options, qp_solver=qp_solver)
    call solver%solve([1.0_wp, 5.0_wp, 5.0_wp, 1.0_wp], istat)
    call solver%get_solution(x, lambda)
    print '(A,I0,A,4F12.8)', 'hs71 with the basis repair: istat ', istat, ', x = ', x
    if (istat /= sqpopt_success) call fail('hs71', 'did not converge')
    if (maxval(abs(x - x_star)) > 1.0e-6_wp) call fail('hs71', 'wrong solution')

    end subroutine test_hs71

    pure real(wp) function det3(b)
    !! the determinant of a 3 x 3 matrix
    real(wp), intent(in) :: b(3,3) !! the matrix
    det3 = b(1,1)*(b(2,2)*b(3,3) - b(2,3)*b(3,2)) - b(1,2)*(b(2,1)*b(3,3) - b(2,3)*b(3,1)) &
         + b(1,3)*(b(2,1)*b(3,2) - b(2,2)*b(3,1))
    end function det3

    subroutine fc(x, f, c, status, data)
    !! HS71's objective and constraints
    real(wp), dimension(:), intent(in)    :: x      !! the variables `dimension(4)`
    real(wp),               intent(out)   :: f      !! the objective
    real(wp), dimension(:), intent(out)   :: c      !! the constraints `dimension(2)`
    integer,                intent(inout) :: status !! (not used)
    class(*), optional,     intent(inout) :: data   !! (not used)
    f = x(1)*x(4)*(x(1)+x(2)+x(3)) + x(3)
    c = [sum(x**2), product(x)]
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! HS71's gradient and Jacobian
    real(wp), dimension(:), intent(in)    :: x        !! the variables `dimension(4)`
    real(wp), dimension(:), intent(out)   :: g        !! the gradient `dimension(4)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! the Jacobian's nonzeros `dimension(8)`
    integer,                intent(in)    :: accuracy !! (not used)
    integer,                intent(inout) :: status   !! (not used)
    class(*), optional,     intent(inout) :: data     !! (not used)
    g = [x(4)*(2.0_wp*x(1)+x(2)+x(3)), x(1)*x(4), x(1)*x(4) + 1.0_wp, x(1)*(x(1)+x(2)+x(3))]
    jac_val(1:4) = 2.0_wp*x
    jac_val(5:8) = [x(2)*x(3)*x(4), x(1)*x(3)*x(4), x(1)*x(2)*x(4), x(1)*x(2)*x(3)]
    end subroutine gjac

    subroutine fail(name, msg)
    !! stop with a failure message
    character(len=*), intent(in) :: name !! the case
    character(len=*), intent(in) :: msg  !! what failed
    print '(A)', 'test_basis_repair FAILED: '//name//': '//msg
    error stop 1
    end subroutine fail

end program test_basis_repair
