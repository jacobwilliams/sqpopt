!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Sparse linear algebra utility routines. All matrices in `sqpopt` are
!  stored in coordinate (COO) format (1-based `irow`/`icol`/`val` triplets,
!  see [[sqpopt_types_module(module):sqpopt_sparse_matrix(type)]]) -- dense
!  \( n \times n \) or \( m \times n \) arrays are never formed. This is
!  the same triplet convention used by the `lusol`, `LSQR`, and `LSMR`
!  dependencies, and by the sparse mode of `nlesolver-fortran`, so that
!  a [[sqpopt_sparse_matrix]] can be passed directly to those solvers.

    module sqpopt_linalg_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_sparse_matrix

    implicit none

    private

    integer, parameter, public :: sqpopt_linsolve_lusol = 1  !! direct sparse solve via `lusol_ez_module` (LU factorization)
    integer, parameter, public :: sqpopt_linsolve_lsqr  = 2  !! iterative sparse solve via `lsqr_module`
    integer, parameter, public :: sqpopt_linsolve_lsmr  = 3  !! iterative sparse solve via `LSMRmodule`

    public :: sparse_matvec
    public :: sparse_matvec_transpose
    public :: solve_sparse_linear_system

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

    ! TODO: implement

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

    ! TODO: implement

    end subroutine sparse_matvec_transpose
!*******************************************************************************

!*******************************************************************************
!>
!  solve the sparse linear system \( A x = b \), dispatching to the
!  solver requested by `method` (one of `sqpopt_linsolve_lusol`,
!  `sqpopt_linsolve_lsqr`, or `sqpopt_linsolve_lsmr`).

    subroutine solve_sparse_linear_system(a, b, method, x, istat)

    type(sqpopt_sparse_matrix), intent(in)  :: a       !! coefficient matrix, `dimension(n,n)`
    real(wp), dimension(:),     intent(in)  :: b       !! right-hand-side vector `dimension(n)`
    integer,                    intent(in)  :: method  !! linear solver strategy to use
    real(wp), dimension(:),     intent(out) :: x       !! solution vector `dimension(n)`
    integer,                    intent(out) :: istat   !! status code (see [[sqpopt_types_module]])

    ! TODO: implement

    end subroutine solve_sparse_linear_system
!*******************************************************************************

    end module sqpopt_linalg_module
!*******************************************************************************
