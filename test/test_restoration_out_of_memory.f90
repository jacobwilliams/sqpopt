program test_restoration_out_of_memory

    !! A feasibility QP of the restoration phase that can't allocate its
    !! matrices must return `sqpopt_out_of_memory` (not an ordinary failed
    !! step), with the point unchanged, and set the `out_of_memory` flag of
    !! the phase's own QP solver (a copy of the main one), which
    !! [[sqpopt_iterate]] checks to stop the solve.
    !!
    !!   subject to  sum(x) >= 1,   at x = 0
    !!
    !! with `n = 6,000,000` variables and the dense QP solver (see
    !! `test_out_of_memory`: its `n x n` matrices need more than the whole
    !! 48-bit address space). The unconstrained step is turned off, though
    !! here it would not be feasible anyway (at the phase's start it is `0`).

    use sqpopt_restoration_module, only: sqpopt_restoration_type
    use sqpopt_problem_module,     only: sqpopt_problem_type
    use sqpopt_qp_solver_module,   only: sqpopt_qp_dense, sqpopt_qp_solver_type
    use sqpopt_linesearch_module,  only: l1_violation
    use sqpopt_types_module,       only: sqpopt_sparse_matrix, sqpopt_out_of_memory
    use sqpopt_kinds,              only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: n = 6000000

    type(sqpopt_problem_type)     :: problem
    type(sqpopt_qp_solver_type)   :: qp
    type(sqpopt_restoration_type) :: restoration
    type(sqpopt_sparse_matrix)    :: jac
    real(wp), dimension(:), allocatable :: x, x_new, big
    real(wp) :: c(1), alpha
    integer  :: i, istat

    write(*,*) '----------------------------'
    write(*,*) 'test_restoration_out_of_memory'
    write(*,*) '----------------------------'

    allocate(x(n), x_new(n), big(n))
    x   = 0.0_wp
    big = 1.0e20_wp

    call problem%set_problem_size(n=n, m=1)
    call problem%set_bounds(x_lb=-big, x_ub=big, c_lb=[1.0_wp], c_ub=[1.0e20_wp])
    deallocate(big)
    jac%nrows = 1
    jac%ncols = n
    jac%nnz   = n
    allocate(jac%irow(n), jac%icol(n), jac%val(n))
    jac%irow = 1
    jac%icol = [(i, i = 1, n)]
    jac%val  = 1.0_wp
    call problem%set_jacobian_sparsity(nnz=n, irow=jac%irow, icol=jac%icol)
    call problem%set_functions(fc=fc, gjac=gjac)
    call problem%reset_evaluations()
    call problem%c(x, c)

    qp%mode = sqpopt_qp_dense
    qp%unconstrained_step = .false.
    call restoration%enter(x, l1_violation(c, problem%c_lb, problem%c_ub), qp)
    call restoration%step(problem, jac, x, c, x_new, alpha, istat)

    print '(A,I0,A,L1)', 'istat=', istat, '  restoration%qp%out_of_memory=', restoration%qp%out_of_memory
    if (istat /= sqpopt_out_of_memory) error stop 'test_restoration_out_of_memory FAILED: status is not sqpopt_out_of_memory'
    if (.not. restoration%qp%out_of_memory) error stop 'test_restoration_out_of_memory FAILED: the flag is not set'
    if (any(x_new /= x) .or. alpha /= 0.0_wp) error stop 'test_restoration_out_of_memory FAILED: the point was changed'

    print '(A)', 'test_restoration_out_of_memory PASSED'

    contains

    subroutine fc(x, f, c, status, data)
    !! the objective (unused here) and the constraint
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    associate(unused => status); end associate
    if (present(data)) continue
    f = 0.0_wp
    c(1) = sum(x)
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient (unused here) and the constraint Jacobian
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the constraint Jacobian at `x`
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    associate(unused => accuracy, unused2 => status, unused3 => x); end associate
    if (present(data)) continue
    g = 0.0_wp
    jac_val = 1.0_wp
    end subroutine gjac

end program test_restoration_out_of_memory
