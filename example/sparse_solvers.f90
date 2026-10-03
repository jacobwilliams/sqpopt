program sparse_solvers

    !! Compares the two sparse solvers of `options%linear_solver` (QDLDL and
    !! MUMPS, see [[sqpopt_symmetric_solver_module]]) on their own, on the
    !! kinds of matrices whose structure decides between them: grids in two
    !! and three dimensions. For each matrix it times the start (QDLDL's
    !! analysis), the first factorization (MUMPS's analysis and
    !! factorization), a second factorization with new values (what each
    !! major iteration costs), and a solve with iterative refinement, and
    !! checks the residual. Without MUMPS, only QDLDL is run.
    !!
    !!    fpm run --example sparse_solvers --profile release -- [--threads=T] [--max-n=N] [--matrix=NAME]
    !!
    !! or, with MUMPS (see the README):
    !!
    !!    pixi run run-mumps --example sparse_solvers --profile release
    !!
    !! `--threads=T` also runs MUMPS on `T` OpenMP threads (default: one
    !! thread only); `--max-n=N` skips matrices of order above `N`, and
    !! `--matrix=NAME` runs only the matrices of that kind (`grid2`, `grid3`,
    !! `kkt2`, or `kkt3`). Without `--max-n`, the largest 3-D matrices take
    !! QDLDL minutes.
    !!
    !! The matrices (`N` points per side):
    !!
    !! * `grid2`, `grid3`: the Laplacian of a 2-D (5-point) or 3-D (7-point)
    !!   grid plus the identity, positive definite;
    !! * `kkt2`, `kkt3`: the KKT matrix
    !!   \( \begin{bmatrix} L + I & J^T \\ J & 0 \end{bmatrix} \) of that
    !!   Hessian with one constraint for every second grid point, coupling it
    !!   with its neighbour in the first direction (`J` has the pattern of a
    !!   discretized derivative), with each constraint's row ordered after its
    !!   variables, as the KKT matrices of the solver are.

    use sqpopt_symmetric_solver_module, only: sqpopt_symmetric_solver_type, sqpopt_has_mumps, &
                                              sqpopt_linear_solver_qdldl, sqpopt_linear_solver_mumps, &
                                              sqpopt_linear_solver_name
    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use, intrinsic :: iso_fortran_env, only: int64

    implicit none

    integer :: threads, max_n, i, ios
    character(len=64) :: arg
    character(len=:), allocatable :: only !! `--matrix` (empty: all of them)

    threads = 1
    max_n = huge(1)
    only = ''
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        if (arg(1:10) == '--threads=') then
            read(arg(11:), *, iostat=ios) threads
            if (ios /= 0 .or. threads < 1) error stop 'sparse_solvers: bad --threads value'
        else if (arg(1:8) == '--max-n=') then
            read(arg(9:), *, iostat=ios) max_n
            if (ios /= 0 .or. max_n < 1) error stop 'sparse_solvers: bad --max-n value'
        else if (arg(1:9) == '--matrix=') then
            only = trim(arg(10:))
        else
            error stop 'sparse_solvers: unknown option (see the header of example/sparse_solvers.f90)'
        end if
    end do

    write(*,'(A8,A9,A10,2X,A14,4A10,A10)') 'matrix', 'order', 'nonzeros', 'solver', 'start', 'factor 1', &
        'factor 2', 'solve', 'residual'
    call run_all('grid2', 2, 300)
    call run_all('grid2', 2, 700)
    call run_all('kkt2',  2, 300)
    call run_all('kkt2',  2, 700)
    call run_all('grid3', 3, 30)
    call run_all('grid3', 3, 50)
    call run_all('grid3', 3, 60)
    call run_all('kkt3',  3, 30)
    call run_all('kkt3',  3, 50)

    contains

    subroutine run_all(name, dims, npts)
    !! one matrix, with each solver
    character(len=*), intent(in) :: name !! the kind of matrix (`grid2`, `grid3`, `kkt2`, or `kkt3`)
    integer,          intent(in) :: dims !! the number of dimensions of the grid (2 or 3)
    integer,          intent(in) :: npts !! the number of points per side
    integer, dimension(:), allocatable :: irow, icol, signs
    real(wp), dimension(:), allocatable :: val
    integer :: n
    if (only /= '' .and. only /= name) return
    call build(name(1:3) == 'kkt', dims, npts, n, irow, icol, val, signs)
    if (n > max_n) return
    call run(name, n, irow, icol, val, signs, sqpopt_linear_solver_qdldl, 1)
    if (sqpopt_has_mumps) then
        call run(name, n, irow, icol, val, signs, sqpopt_linear_solver_mumps, 1)
        if (threads > 1) call run(name, n, irow, icol, val, signs, sqpopt_linear_solver_mumps, threads)
    end if
    end subroutine run_all

    subroutine build(kkt, dims, npts, n, irow, icol, val, signs)
    !! the lower triangle of the matrix, in coordinate form, and the expected sign of each row's pivot
    logical, intent(in) :: kkt  !! the KKT matrix (else the grid's Hessian alone)
    integer, intent(in) :: dims !! the number of dimensions of the grid (2 or 3)
    integer, intent(in) :: npts !! the number of points per side
    integer, intent(out) :: n   !! the order of the matrix
    integer, dimension(:), allocatable, intent(out) :: irow  !! row indices
    integer, dimension(:), allocatable, intent(out) :: icol  !! column indices
    real(wp), dimension(:), allocatable, intent(out) :: val  !! values
    integer, dimension(:), allocatable, intent(out) :: signs !! expected sign of each row's pivot
    integer :: nv, m, k, i, j, l, p, c, nnz, stride(3)
    nv = npts**dims
    m = 0
    if (kkt) m = nv/2
    n = nv + m
    nnz = nv*(1 + dims) + 2*m
    allocate(irow(nnz), icol(nnz), val(nnz), signs(n))
    signs(1:nv) = 1
    signs(nv+1:) = -1
    stride = [1, npts, npts**2]
    p = 0
    do k = 1, nv
        p = p + 1
        irow(p) = k
        icol(p) = k
        val(p)  = 1.0_wp + 2.0_wp*dims
        ! (the neighbours with a smaller index, in each direction)
        l = k - 1
        do j = 1, dims
            i = mod(l/stride(j), npts)
            if (i > 0) then
                p = p + 1
                irow(p) = k
                icol(p) = k - stride(j)
                val(p)  = -1.0_wp
            end if
        end do
    end do
    ! the constraints: x(2c-1) - x(2c), for c = 1..m
    do c = 1, m
        p = p + 1
        irow(p) = nv + c
        icol(p) = 2*c - 1
        val(p)  = 1.0_wp
        p = p + 1
        irow(p) = nv + c
        icol(p) = 2*c
        val(p)  = -1.0_wp
    end do
    irow = irow(1:p)
    icol = icol(1:p)
    val  = val(1:p)
    end subroutine build

    subroutine run(name, n, irow, icol, val, signs, which, nthreads)
    !! time one solver on one matrix, and print a line
    character(len=*),       intent(in) :: name     !! the matrix's name
    integer,                intent(in) :: n        !! its order
    integer, dimension(:),  intent(in) :: irow     !! row indices
    integer, dimension(:),  intent(in) :: icol     !! column indices
    real(wp), dimension(:), intent(in) :: val      !! values
    integer, dimension(:),  intent(in) :: signs    !! expected sign of each row's pivot
    integer,                intent(in) :: which    !! the solver (`sqpopt_linear_solver_*`)
    integer,                intent(in) :: nthreads !! MUMPS's number of threads
    type(sqpopt_symmetric_solver_type) :: solver
    real(wp), dimension(:), allocatable :: b, x, ax, val2
    real(wp) :: t_start, t_f1, t_f2, t_solve, res
    character(len=14) :: label
    logical :: ok
    integer :: i
    allocate(b(n), x(n), ax(n))
    do i = 1, n
        b(i) = 1.0_wp + real(mod(i, 11), wp)/11.0_wp
    end do
    val2 = 1.01_wp*val   ! (new values, the same pattern)
    t_start = seconds()
    call solver%initialize(n, irow, icol, ok, threads=nthreads, solver=which, signs=signs)
    t_start = seconds() - t_start
    if (.not. ok) error stop 'sparse_solvers: the solver could not be started'
    t_f1 = seconds()
    call solver%factor(val, ok)
    t_f1 = seconds() - t_f1
    if (.not. ok) error stop 'sparse_solvers: factorization failed'
    t_f2 = seconds()
    call solver%factor(val2, ok)
    t_f2 = seconds() - t_f2
    if (.not. ok) error stop 'sparse_solvers: second factorization failed'
    x = b
    t_solve = seconds()
    call solver%solve(x, ok)
    t_solve = seconds() - t_solve
    if (.not. ok) error stop 'sparse_solvers: solve failed'
    call solver%multiply(x, ax)
    res = maxval(abs(ax - b))/maxval(abs(b))
    call solver%destroy()
    label = sqpopt_linear_solver_name(which)
    if (nthreads > 1) write(label, '(A,I0,A)') 'MUMPS, ', nthreads, ' thr.'
    write(*,'(A8,I9,I10,2X,A14,4F10.3,ES10.1)') name, n, size(irow), label, t_start, t_f1, t_f2, t_solve, res
    end subroutine run

    real(wp) function seconds()
    !! the wall-clock time, in seconds
    integer(int64) :: count, rate
    call system_clock(count, rate)
    seconds = real(count, wp)/real(rate, wp)
    end function seconds

end program sparse_solvers
