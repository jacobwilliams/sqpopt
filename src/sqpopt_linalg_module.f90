!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Sparse linear algebra utility routines. All matrices in `sqpopt` are
!  stored in coordinate (COO) format (1-based `irow`/`icol`/`val` triplets,
!  see [[sqpopt_types_module(module):sqpopt_sparse_matrix(type)]]) -- dense
!  \( n \times n \) or \( m \times n \) arrays are never formed. This is
!  the same triplet convention used by the `lusol`, `LSQR`, and `LSMR`
!  dependencies, so that a [[sqpopt_sparse_matrix]] can be passed directly
!  to those solvers.

    module sqpopt_linalg_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_sparse_matrix

    implicit none

    private

    public :: sparse_matvec
    public :: sparse_matvec_transpose

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  compute the sparse matrix-vector product \( y = A x \), operating
!  directly on the COO triplets (never forms a dense matrix).

    subroutine sparse_matvec(a, x, y)

    type(sqpopt_sparse_matrix), intent(in)  :: a  !! sparse matrix, `dimension(nrows,ncols)`
    real(wp), dimension(:),     intent(in)  :: x  !! vector to multiply `dimension(ncols)`
    real(wp), dimension(:),     intent(out) :: y  !! result vector `dimension(nrows)`

    integer :: k  !! nonzero element counter

    y = 0.0_wp
    do k = 1, a%nnz
        y(a%irow(k)) = y(a%irow(k)) + a%val(k)*x(a%icol(k))
    end do

    end subroutine sparse_matvec
!*******************************************************************************

!*******************************************************************************
!>
!  compute the sparse transposed matrix-vector product \( x = A^T y \),
!  operating directly on the COO triplets (never forms a dense matrix).

    subroutine sparse_matvec_transpose(a, y, x)

    type(sqpopt_sparse_matrix), intent(in)  :: a  !! sparse matrix, `dimension(nrows,ncols)`
    real(wp), dimension(:),     intent(in)  :: y  !! vector to multiply `dimension(nrows)`
    real(wp), dimension(:),     intent(out) :: x  !! result vector `dimension(ncols)`

    integer :: k  !! nonzero element counter

    x = 0.0_wp
    do k = 1, a%nnz
        x(a%icol(k)) = x(a%icol(k)) + a%val(k)*y(a%irow(k))
    end do

    end subroutine sparse_matvec_transpose
!*******************************************************************************

    end module sqpopt_linalg_module
!*******************************************************************************
