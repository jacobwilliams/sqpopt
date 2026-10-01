program test_out_of_memory

    !! A dense QP solve that can't allocate its matrices must end the solve
    !! with `sqpopt_out_of_memory`, not crash it.
    !!
    !!   minimize  sum (x_i - 1)**2      (no constraints)
    !!
    !! with `n = 6,000,000` variables and the dense QP solver forced
    !! (`sqpopt_qp_auto` would choose the sparse one): its `n x n` matrices
    !! need `8 n**2 = 2.9e14` bytes each (`real64`), more than the whole
    !! 48-bit address space, so the allocation fails on any machine. The
    !! starting point must be returned unchanged. (One L-BFGS pair, so the
    !! test's own memory is a few dozen vectors of length `n`.)

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_dense
    use sqpopt_types_module,     only: sqpopt_out_of_memory, sqpopt_results_type
    use sqpopt_kinds,            only: wp => sqpopt_module_wp

    implicit none

    integer, parameter :: n = 6000000

    type(sqpopt_type)         :: solver
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_results_type) :: r
    real(wp), dimension(:), allocatable :: x0, xsol, lam, big
    integer, dimension(0) :: none
    integer :: istat

    write(*,*) '----------------------------'
    write(*,*) 'test_out_of_memory'
    write(*,*) '----------------------------'

    allocate(x0(n), xsol(n), lam(0), big(n))
    x0  = 0.0_wp
    big = 1.0e20_wp

    call problem%set_problem_size(n=n, m=0)
    call problem%set_bounds(x_lb=-big, x_ub=big, c_lb=lam, c_ub=lam)
    call problem%set_jacobian_sparsity(nnz=0, irow=none, icol=none)
    call problem%set_functions(fc=fc, gjac=gjac)

    options%qp_solver_mode = sqpopt_qp_dense
    options%lbfgs_memory   = 1
    options%scaling        = .false.
    call solver%initialize(problem=problem, options=options)
    call solver%solve(x0, istat)
    call solver%get_solution(xsol, lam)
    call solver%get_results(r)

    print '(A,I0,2A)', 'istat=', istat, '  ', solver%status_message()
    if (istat /= sqpopt_out_of_memory) error stop 'test_out_of_memory FAILED: status is not sqpopt_out_of_memory'
    if (r%istat /= sqpopt_out_of_memory) error stop 'test_out_of_memory FAILED: results%istat'
    if (any(xsol /= x0)) error stop 'test_out_of_memory FAILED: the point was changed'

    print '(A)', 'test_out_of_memory PASSED'

    contains

    subroutine fc(x, f, c, status, data)
    !! the objective (there are no constraints)
    real(wp), dimension(:), intent(in)    :: x      !! point `dimension(n)`
    real(wp),               intent(out)   :: f      !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
    integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data   !! the user data passed to `set_functions` (if any)
    associate(unused => status, unused2 => c); end associate
    if (present(data)) continue
    f = sum((x - 1.0_wp)**2)
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! the objective gradient (the Jacobian is empty)
    real(wp), dimension(:), intent(in)    :: x        !! point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
    real(wp), dimension(:), intent(out)   :: jac_val  !! nonzero values of the constraint Jacobian at `x`
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop the solver
    class(*), optional,     intent(inout) :: data     !! the user data passed to `set_functions` (if any)
    associate(unused => accuracy, unused2 => status, unused3 => jac_val); end associate
    if (present(data)) continue
    g = 2.0_wp*(x - 1.0_wp)
    end subroutine gjac

end program test_out_of_memory
