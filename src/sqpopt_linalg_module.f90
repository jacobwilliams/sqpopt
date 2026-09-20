!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Linear algebra utility routines, including conversions between dense
!  and sparse matrix formats.

    module sqpopt_linalg_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_sparse_matrix

    implicit none

    private

    public :: dense_to_sparse
    public :: sparse_to_dense
    public :: solve_linear_system

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  convert a dense matrix to a sparse (COO) matrix, discarding zero elements.

    subroutine dense_to_sparse(dense, sparse)

    real(wp), dimension(:,:),      intent(in)  :: dense  !! dense matrix `dimension(nrows,ncols)`
    type(sqpopt_sparse_matrix),    intent(out) :: sparse !! equivalent sparse matrix

    ! TODO: implement

    end subroutine dense_to_sparse
!*******************************************************************************

!*******************************************************************************
!>
!  convert a sparse (COO) matrix to a dense matrix.

    subroutine sparse_to_dense(sparse, dense)

    type(sqpopt_sparse_matrix), intent(in)  :: sparse !! sparse matrix
    real(wp), dimension(:,:),   intent(out) :: dense  !! equivalent dense matrix `dimension(nrows,ncols)`

    ! TODO: implement

    end subroutine sparse_to_dense
!*******************************************************************************

!*******************************************************************************
!>
!  solve the dense linear system \( A x = b \).

    subroutine solve_linear_system(a, b, x, istat)

    real(wp), dimension(:,:), intent(in)  :: a      !! coefficient matrix `dimension(n,n)`
    real(wp), dimension(:),   intent(in)  :: b      !! right-hand-side vector `dimension(n)`
    real(wp), dimension(:),   intent(out) :: x      !! solution vector `dimension(n)`
    integer,                   intent(out) :: istat  !! status code (see [[sqpopt_types_module]])

    ! TODO: implement

    end subroutine solve_linear_system
!*******************************************************************************

    end module sqpopt_linalg_module
!*******************************************************************************
