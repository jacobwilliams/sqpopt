!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The MUMPS backend of [[sqpopt_symmetric_solver_module]]
!  (`options%linear_solver = sqpopt_linear_solver_mumps`):
!  [MUMPS](https://mumps-solver.org), a multifrontal \( LDL^T \)
!  factorization with pivoting and OpenMP threads. Its inertia is always
!  known (`icntl(13) = 1` makes the count of negative pivots exact), and a
!  singular matrix's null pivots are detected (`icntl(24) = 1`).
!
!  This is the library's only interface to MUMPS, and the only source file
!  that depends on the `HAS_MUMPS` preprocessor directive: without it,
!  `sqpopt_has_mumps` is false and [[mumps_start]] always fails. With it,
!  the library must be linked with the sequential, double precision MUMPS
!  library, and compiled in double precision (the default kind, see
!  [[sqpopt_kinds]]).
!
!  **Threads.** MUMPS and the BLAS under it can use OpenMP threads, if they
!  were built with OpenMP (conda-forge's `mumps-seq` is): the `threads`
!  argument of [[mumps_start]], `1` by default, or `0` to leave it to the
!  OpenMP environment (`OMP_NUM_THREADS`, or every core). On small matrices
!  the threads cost far more time than they save (the 305
!  Hock-Schittkowski problems took 24.5 s with 10 threads, and 1.4 s with
!  one). Whether more threads help a large matrix depends on its factors:
!  the threads work inside the dense blocks of the factorization (its
!  frontal matrices), so they need those to be large. Measured (analysis,
!  two factorizations, and a solve; 1, 2, 4, and 8 threads): a 3-D grid
!  matrix of order 216,000 took 11.5, 7.0, 4.9, and 4.8 s; a 2-D grid
!  matrix of order 490,000 took 2.0 s with any number; and the banded KKT
!  matrices of `example/benchmark_large.f90` gained nothing either. MUMPS
!  sets the number of threads when it is called and restores it when it
!  returns, so the rest of the program is not affected.

    module sqpopt_mumps_ldl_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_sparse_ldl_module, only: sqpopt_sparse_ldl_type

    implicit none

    private

#ifdef HAS_MUMPS
#if defined(REAL32) || defined(REAL128)
#error "HAS_MUMPS needs the default (double precision) real kind: don't define REAL32 or REAL128 with it"
#endif
    include 'dmumps_struc.h'

    interface
        subroutine dmumps(id)
            !! the MUMPS driver (double precision): initializes, analyses, factors,
            !! solves, or frees, as `id%job` says
            import :: dmumps_struc
            implicit none
            type(dmumps_struc), intent(inout) :: id !! the MUMPS instance (its inputs and outputs)
        end subroutine dmumps
    end interface

    logical, parameter, public :: sqpopt_has_mumps = .true.  !! whether the library was built with MUMPS (the
                                                             !! `HAS_MUMPS` preprocessor directive), which
                                                             !! `sqpopt_linear_solver_mumps` requires
#else
    logical, parameter, public :: sqpopt_has_mumps = .false. !! whether the library was built with MUMPS (the
                                                             !! `HAS_MUMPS` preprocessor directive), which
                                                             !! `sqpopt_linear_solver_mumps` requires
#endif

    type, extends(sqpopt_sparse_ldl_type), public :: sqpopt_mumps_ldl_type
        !! the MUMPS backend (see the module documentation). It holds a
        !! pointer, so it must not be copied.
        private
        logical :: analysed = .false.               !! whether the sparsity pattern has been analysed
#ifdef HAS_MUMPS
        type(dmumps_struc), pointer :: id => null() !! the MUMPS instance (with the pattern and the values)
#endif
        contains
        procedure :: start      => mumps_start
        procedure :: refactor   => mumps_refactor
        procedure :: back_solve => mumps_back_solve
        procedure :: multiply   => mumps_multiply
        procedure :: free       => mumps_free
    end type sqpopt_mumps_ldl_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  start MUMPS for the pattern (it is analysed in the first factorization).
!  `ok` is false if the library was built without MUMPS, or MUMPS couldn't
!  be started. `signs` is not used (MUMPS pivots).

    subroutine mumps_start(me, n, irow, icol, ok, threads, signs)

    class(sqpopt_mumps_ldl_type), intent(inout) :: me
    integer,               intent(in)  :: n       !! order of the matrix
    integer, dimension(:), intent(in)  :: irow    !! row indices of the entries `dimension(nnz)`
    integer, dimension(:), intent(in)  :: icol    !! column indices of the entries `dimension(nnz)`
    logical,               intent(out) :: ok      !! whether the backend is ready
    integer, optional,     intent(in)  :: threads !! number of OpenMP threads (default `1`; `0`: as the OpenMP
                                                  !! environment says)
    integer, dimension(:), optional, intent(in) :: signs !! (not used)

#ifdef HAS_MUMPS
    integer :: nnz, alloc_stat
#endif

    ok = .false.
    if (present(signs)) continue   ! (MUMPS pivots)
#ifdef HAS_MUMPS
    nnz = size(irow)
    allocate(me%id, stat=alloc_stat)
    if (alloc_stat /= 0) then
        nullify(me%id)
        me%out_of_memory = .true.
        return
    end if
    me%id%comm = 0   ! (the sequential library's MPI is a stub: no communicator is needed)
    me%id%sym  = 2   ! general symmetric
    me%id%par  = 1
    me%id%job  = -1
    call dmumps(me%id)
    if (me%id%infog(1) < 0) then
        deallocate(me%id)
        return
    end if

    me%id%icntl(1:3) = -1   ! no printed output
    me%id%icntl(4)   = 0
    me%id%icntl(13)  = 1    ! (needed for the number of negative pivots to be exact)
    me%id%icntl(16)  = 1    ! the number of OpenMP threads (see the module documentation)
    if (present(threads)) me%id%icntl(16) = max(threads, 0)
    me%id%icntl(24)  = 1    ! detect null pivots, rather than fail on a singular matrix

    me%id%n   = n
    me%id%nnz = int(nnz, kind(me%id%nnz))
    allocate(me%id%irn(nnz), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%id%jcn(nnz), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%id%a(nnz), stat=alloc_stat)
    if (alloc_stat == 0) allocate(me%id%rhs(n), stat=alloc_stat)
    if (alloc_stat /= 0) then
        call me%free()
        me%out_of_memory = .true.
        return
    end if
    me%id%irn = irow
    me%id%jcn = icol
    ok = .true.
#endif

    end subroutine mumps_start
!*******************************************************************************

!*******************************************************************************
!>
!  factor the matrix with the values `val`. The pattern is analysed in the
!  first call, with that call's values. If MUMPS runs out of its estimated
!  workspace, the factorization is repeated with twice the allowance, up
!  to 20 times (as IPOPT does), and the larger allowance is kept for the
!  following factorizations.

    subroutine mumps_refactor(me, val, ok)

    class(sqpopt_mumps_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: val !! values of the entries `dimension(nnz)`
    logical,                intent(out) :: ok  !! whether the factorization succeeded

#ifdef HAS_MUMPS
    integer :: attempt
    integer, parameter :: max_attempts = 20 !! (each doubles the workspace increase; a very indefinite matrix,
                                            !! whose pivots are delayed, can need many times the estimate)
#endif

    ok = .false.
    me%inertia_known = .true.
#ifdef HAS_MUMPS
    me%id%a = val
    do attempt = 1, max_attempts
        me%id%job = merge(2, 4, me%analysed)   ! (4: analyse, then factor)
        call dmumps(me%id)
        if (me%id%infog(1) /= -8 .and. me%id%infog(1) /= -9) exit
        me%id%icntl(14) = 2*max(me%id%icntl(14), 20)   ! more workspace (a percentage of the estimate)
    end do
    if (me%id%infog(1) < 0) then
        ! (-13: an allocation failed; -5, -7: an allocation of the analysis failed;
        ! -8, -9: still not enough workspace)
        if (any(me%id%infog(1) == [-13, -5, -7, -8, -9])) me%out_of_memory = .true.
        return
    end if
    me%analysed   = .true.
    me%n_negative = me%id%infog(12)
    me%n_null     = me%id%infog(28)
    ok = .true.
#endif

    end subroutine mumps_refactor
!*******************************************************************************

!*******************************************************************************
!>
!  one solve with the factors, without refinement.

    subroutine mumps_back_solve(me, v, ok)

    class(sqpopt_mumps_ldl_type), intent(inout) :: me
    real(wp), dimension(:), intent(inout) :: v  !! the right-hand side, overwritten by the solution
    logical,                intent(out)   :: ok !! whether it was solved

    ok = .false.
#ifdef HAS_MUMPS
    me%id%rhs = v
    me%id%job = 3
    call dmumps(me%id)
    ok = me%id%infog(1) >= 0
    if (ok) then
        v = me%id%rhs
    else if (me%id%infog(1) == -13) then
        me%out_of_memory = .true.   ! (an allocation failed)
    end if
#endif

    end subroutine mumps_back_solve
!*******************************************************************************

!*******************************************************************************
!>
!  the product \( y = A x \) with the values last factored.

    subroutine mumps_multiply(me, x, y)

    class(sqpopt_mumps_ldl_type), intent(in) :: me
    real(wp), dimension(:), intent(in)  :: x !! the vector `dimension(n)`
    real(wp), dimension(:), intent(out) :: y !! the product `dimension(n)`

#ifdef HAS_MUMPS
    integer :: k, i, j
#endif

    y = 0.0_wp
#ifdef HAS_MUMPS
    if (.not. associated(me%id)) return
    do k = 1, size(me%id%a)
        i = me%id%irn(k)
        j = me%id%jcn(k)
        y(i) = y(i) + me%id%a(k)*x(j)
        if (i /= j) y(j) = y(j) + me%id%a(k)*x(i)
    end do
#endif

    end subroutine mumps_multiply
!*******************************************************************************

!*******************************************************************************
!>
!  free everything (it does nothing if MUMPS was never started).

    subroutine mumps_free(me)

    class(sqpopt_mumps_ldl_type), intent(inout) :: me

#ifdef HAS_MUMPS
    if (associated(me%id)) then
        if (associated(me%id%irn)) deallocate(me%id%irn)
        if (associated(me%id%jcn)) deallocate(me%id%jcn)
        if (associated(me%id%a))   deallocate(me%id%a)
        if (associated(me%id%rhs)) deallocate(me%id%rhs)
        me%id%job = -2
        call dmumps(me%id)
        deallocate(me%id)
    end if
#endif
    me%analysed = .false.

    end subroutine mumps_free
!*******************************************************************************

    end module sqpopt_mumps_ldl_module
!*******************************************************************************
