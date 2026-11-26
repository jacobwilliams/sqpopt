!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The interface of a sparse \( LDL^T \) backend: the abstract type that
!  each sparse solver extends ([[sqpopt_qdldl_ldl_type]],
!  [[sqpopt_mumps_ldl_type]]), and that [[sqpopt_symmetric_solver_type]]
!  holds and calls. A backend only factors and solves; the choice of the
!  backend, the iterative refinement, the counts, and the timing are
!  [[sqpopt_symmetric_solver_module]]'s, and what the inertia means for the
!  optimization (e.g. what to do when it is unknown) is
!  [[sqpopt_kkt_module]]'s.
!
!  A backend reports the inertia of the matrix it factored: `n_negative`
!  negative and `n_null` zero eigenvalues, if `inertia_known`. A backend
!  without pivoting may not always be able to tell (see
!  [[sqpopt_qdldl_ldl_module]]): it then says so, with
!  `inertia_known = .false.`, rather than guess.

    module sqpopt_sparse_ldl_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    type, abstract, public :: sqpopt_sparse_ldl_type
        !! a sparse symmetric \( LDL^T \) backend, for matrices with one
        !! sparsity pattern (see the module documentation)
        integer :: n_negative = 0           !! of the last factorization: the number of negative eigenvalues
        integer :: n_null = 0               !! of the last factorization: the number of zero eigenvalues
        logical :: inertia_known = .true.   !! of the last factorization: whether `n_negative` and `n_null` are
                                            !! the matrix's (else they are those of the factors, which are not
                                            !! the matrix's: see the backend's documentation)
        logical :: out_of_memory = .false.  !! whether an allocation of the backend failed (it stays set)
        contains
        procedure(start_interface),      deferred :: start      !! analyse the sparsity pattern
        procedure(refactor_interface),   deferred :: refactor   !! factor new values
        procedure(back_solve_interface), deferred :: back_solve !! one solve with the factors
        procedure(multiply_interface),   deferred :: multiply   !! the product with the last values
        procedure(free_interface),       deferred :: free       !! free everything
    end type sqpopt_sparse_ldl_type

    abstract interface

        subroutine start_interface(me, n, irow, icol, ok, threads, signs)
            !! start the backend for symmetric matrices of order `n` with the
            !! sparsity pattern `irow`/`icol` (each entry stands for itself and,
            !! off the diagonal, its mirror image; entries given more than once
            !! are added together). `ok` is false if it couldn't be started.
            import :: sqpopt_sparse_ldl_type
            implicit none
            class(sqpopt_sparse_ldl_type), intent(inout) :: me
            integer,               intent(in)  :: n       !! order of the matrix
            integer, dimension(:), intent(in)  :: irow    !! row indices of the entries `dimension(nnz)`
            integer, dimension(:), intent(in)  :: icol    !! column indices of the entries `dimension(nnz)`
            logical,               intent(out) :: ok      !! whether the backend is ready
            integer, optional,     intent(in)  :: threads !! number of OpenMP threads, if the backend uses them
            integer, dimension(:), optional, intent(in) :: signs !! the expected sign of each row's pivot (`-1`,
                                                                 !! `0`, or `+1`), if the backend uses it
                                                                 !! `dimension(n)`
        end subroutine start_interface

        subroutine refactor_interface(me, val, ok)
            !! factor the matrix with the values `val` (in the order of the
            !! pattern), and set `n_negative`, `n_null`, and `inertia_known`
            import :: sqpopt_sparse_ldl_type, wp
            implicit none
            class(sqpopt_sparse_ldl_type), intent(inout) :: me
            real(wp), dimension(:), intent(in)  :: val !! values of the entries `dimension(nnz)`
            logical,                intent(out) :: ok  !! whether the factorization succeeded
        end subroutine refactor_interface

        subroutine back_solve_interface(me, v, ok)
            !! one solve with the factors, without refinement: `v` is overwritten
            !! by \( A^{-1} v \) (if not `ok`, `v` is unchanged)
            import :: sqpopt_sparse_ldl_type, wp
            implicit none
            class(sqpopt_sparse_ldl_type), intent(inout) :: me
            real(wp), dimension(:), intent(inout) :: v  !! the right-hand side, overwritten by the solution
            logical,                intent(out)   :: ok !! whether it was solved
        end subroutine back_solve_interface

        subroutine multiply_interface(me, x, y)
            !! the product \( y = A x \) with the values last factored
            import :: sqpopt_sparse_ldl_type, wp
            implicit none
            class(sqpopt_sparse_ldl_type), intent(in) :: me
            real(wp), dimension(:), intent(in)  :: x !! the vector `dimension(n)`
            real(wp), dimension(:), intent(out) :: y !! the product `dimension(n)`
        end subroutine multiply_interface

        subroutine free_interface(me)
            !! free everything the backend holds
            import :: sqpopt_sparse_ldl_type
            implicit none
            class(sqpopt_sparse_ldl_type), intent(inout) :: me
        end subroutine free_interface

    end interface

    end module sqpopt_sparse_ldl_module
!*******************************************************************************
