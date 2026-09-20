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
    use sqpopt_types_module, only: sqpopt_sparse_matrix, sqpopt_error
    use lusol_ez_module, only: lusol_solve => solve
    use lsqr_module,     only: lsqr_solver_ez
    use lsmrModule,      only: lsmr_ez

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

!*******************************************************************************
!>
!  solve the sparse linear system \( A x = b \), dispatching to the
!  solver requested by `method` (one of `sqpopt_linsolve_lusol`,
!  `sqpopt_linsolve_lsqr`, or `sqpopt_linsolve_lsmr`). `lusol` performs a
!  direct sparse LU factorization; `lsqr`/`lsmr` are matrix-free iterative
!  least-squares solvers (also applicable to non-square/rank-deficient `A`).

    subroutine solve_sparse_linear_system(a, b, method, x, istat)

    type(sqpopt_sparse_matrix), intent(in)  :: a       !! coefficient matrix, `dimension(n,n)`
    real(wp), dimension(:),     intent(in)  :: b       !! right-hand-side vector `dimension(n)`
    integer,                    intent(in)  :: method  !! linear solver strategy to use
    real(wp), dimension(:),     intent(out) :: x       !! solution vector `dimension(n)`
    integer,                    intent(out) :: istat   !! status code (see [[sqpopt_types_module]])

    type(lsqr_solver_ez) :: lsqr
    integer :: lsqr_istop, lsmr_istop, lsmr_itn
    real(wp) :: lsmr_normA, lsmr_condA, lsmr_normr, lsmr_normAr, lsmr_normx

    select case (method)
    case (sqpopt_linsolve_lusol)
        call lusol_solve(a%ncols, a%nrows, a%nnz, a%irow, a%icol, a%val, b, x, istat)
    case (sqpopt_linsolve_lsqr)
        call lsqr%initialize(a%nrows, a%ncols, a%val, a%irow, a%icol)
        call lsqr%solve(b, 0.0_wp, x, lsqr_istop)
        istat = merge(0, sqpopt_error, lsqr_istop >= 0 .and. lsqr_istop <= 3)
    case (sqpopt_linsolve_lsmr)
        call lsmr_ez(a%nrows, a%ncols, a%irow, a%icol, a%val, b, 0.0_wp, &
                     1.0e-10_wp, 1.0e-10_wp, 1.0e8_wp, 4*max(a%nrows,a%ncols), 0, 0, &
                     x, lsmr_istop, lsmr_itn, lsmr_normA, lsmr_condA, lsmr_normr, lsmr_normAr, lsmr_normx)
        istat = merge(0, sqpopt_error, lsmr_istop >= 0 .and. lsmr_istop <= 3)
    case default
        istat = sqpopt_error
    end select

    end subroutine solve_sparse_linear_system
!*******************************************************************************

    end module sqpopt_linalg_module
!*******************************************************************************
