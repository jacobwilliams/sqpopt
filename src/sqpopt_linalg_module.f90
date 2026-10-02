!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Sparse linear algebra utility routines. All matrices in `sqpopt` are
!  stored in coordinate (COO) format (1-based `irow`/`icol`/`val` triplets,
!  see [[sqpopt_types_module(module):sqpopt_sparse_matrix(type)]]) -- dense
!  \( n \times n \) or \( m \times n \) arrays are never formed. This is
!  the same triplet convention used by the `lusol` and `LSQR`
!  dependencies, so that a [[sqpopt_sparse_matrix]] can be passed directly
!  to those solvers.

    module sqpopt_linalg_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_sparse_matrix
    use lusol,               only: lu1fac, lu6sol, lu8rpc
    use lusol_precision,     only: ip, rp

    implicit none

    private

    public :: sparse_matvec
    public :: sparse_matvec_transpose
    public :: independent_columns

    type, public :: sqpopt_lu_type
        !! sparse LU factors of a square matrix (`LUSOL`), for repeated solves
        !! with it and its transpose. See [[sqpopt_lu_type(type):factorize(bound)]].
        private
        integer(ip) :: n = 0
        integer(ip) :: lena = 0
        integer(ip) :: luparm(30) = 0
        real(rp)    :: parmlu(30) = 0.0_rp
        real(rp),    dimension(:), allocatable :: a
        integer(ip), dimension(:), allocatable :: indc, indr, p, q, lenc, lenr, locc, locr
        real(rp),    dimension(:), allocatable :: v, w !! work vectors of the solves and the column
                                                       !! replacements `dimension(n)` (kept here so that a
                                                       !! solve allocates nothing)
        contains
        procedure, public :: factorize => lu_factorize
        procedure, public :: solve     => lu_solve
        procedure, public :: replace_column => lu_replace_column
    end type sqpopt_lu_type

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
!  find a maximal linearly independent subset of the columns of the sparse
!  matrix `A` (`nrows x ncols`, COO triplets), with one rank-revealing
!  sparse LU factorization (`LUSOL`'s `lu1fac` with threshold complete
!  pivoting).
!
!  Column `j` is judged dependent if it gets no pivot, or if its diagonal
!  of `U` is at most `rel_tol` times the largest element of its column of
!  `U` (or at most `abs_tol`).
!
!  Which columns of a dependent group are kept is up to the pivoting, but
!  it can be steered. Unit columns (a single entry of exactly `+/-1`) are
!  always kept first (at most one per row), and complete pivoting prefers
!  large elements, so columns scaled up are preferred over columns scaled
!  down. (Scaling doesn't affect the relative test.)
!
!  `istat` is `0` on success, or else `lu1fac`'s `inform` code (it is never
!  `1`, which only means that some columns are dependent).

    subroutine independent_columns(nrows, ncols, irow, icol, val, rel_tol, abs_tol, independent, istat)

    integer,                intent(in)  :: nrows       !! number of rows of `A`
    integer,                intent(in)  :: ncols       !! number of columns of `A`
    integer,  dimension(:), intent(in)  :: irow        !! row indices of the nonzeros of `A`
    integer,  dimension(:), intent(in)  :: icol        !! column indices of the nonzeros of `A`
    real(wp), dimension(:), intent(in)  :: val         !! the nonzeros of `A` (no duplicate `(irow,icol)` pairs)
    real(wp),               intent(in)  :: rel_tol     !! relative singularity tolerance (e.g. `1e-8`)
    real(wp),               intent(in)  :: abs_tol     !! absolute singularity tolerance
    logical,  dimension(:), intent(out) :: independent !! `dimension(ncols)`: whether each column is in the subset
    integer,                intent(out) :: istat       !! status (0 = success)

    integer(ip) :: m, n, nelem, lena, inform, attempt
    integer(ip) :: luparm(30)
    real(rp)    :: parmlu(30)
    real(rp),    dimension(:), allocatable :: a, w
    integer(ip), dimension(:), allocatable :: indc, indr, p, q, lenc, lenr, locc, locr, &
                                              iploc, iqloc, ipinv, iqinv

    independent = .false.
    istat = 0
    if (ncols == 0) return
    if (nrows == 0 .or. size(val) == 0) return   ! (every column is zero)

    m     = nrows
    n     = ncols
    nelem = size(val)
    lena  = 1 + max(5*nelem, 10*m, 10*n, 10000_ip)
    allocate(p(m), q(n), lenc(n), lenr(m), locc(n), locr(m), &
             iploc(n), iqloc(m), ipinv(m), iqinv(n), w(n))

    do attempt = 1, 3   ! (enlarging the workspace if `lu1fac` asks for it)
        if (allocated(a)) deallocate(a, indc, indr)
        allocate(a(lena), indc(lena), indr(lena))
        a(1:nelem)    = real(val, rp)
        indc(1:nelem) = irow   ! (LUSOL: row indices in `indc`, column indices in `indr`)
        indr(1:nelem) = icol

        luparm = 0
        luparm(1) = 6      ! nout
        luparm(2) = -1     ! lprint: no output (not even the "singular" message)
        luparm(3) = 5      ! maxcol
        luparm(6) = 2      ! TCP: threshold complete pivoting (rank revealing)
        luparm(8) = 1      ! keepLU (needed for the relative singularity test)
        parmlu = 0.0_rp
        parmlu(1) = 5.0_rp                          ! Ltol1 (small, for a reliable rank)
        parmlu(2) = 5.0_rp                          ! Ltol2
        parmlu(3) = epsilon(1.0_rp)**0.8_rp         ! small: entries treated as zero
        parmlu(4) = real(abs_tol, rp)               ! Utol1
        parmlu(5) = real(rel_tol, rp)               ! Utol2
        parmlu(6) = 3.0_rp                          ! Uspace
        parmlu(7) = 0.3_rp                          ! dens1
        parmlu(8) = 0.5_rp                          ! dens2

        call lu1fac(m, n, nelem, lena, luparm, parmlu, a, indc, indr, p, q, &
                    lenc, lenr, locc, locr, iploc, iqloc, ipinv, iqinv, w, inform)
        if (inform /= 7) exit
        lena = max(2*lena, luparm(13) + 1)
    end do

    if (inform == 0 .or. inform == 1) then
        independent = w > 0.0_rp
    else
        istat = int(inform)
    end if

    end subroutine independent_columns
!*******************************************************************************

!*******************************************************************************
!>
!  factorize the square `n x n` sparse matrix `A` (COO triplets, no
!  duplicate `(irow,icol)` pairs) with `LUSOL`'s `lu1fac` (threshold rook
!  pivoting, for stability). `istat` is `0` on success, `1` if `A` appears
!  to be singular (see [[independent_columns]] for the test), or else
!  `lu1fac`'s `inform` code; the factors can only be used if it is `0`.

    subroutine lu_factorize(me, n, irow, icol, val, rel_tol, istat)

    class(sqpopt_lu_type),  intent(inout) :: me
    integer,                intent(in)    :: n       !! order of `A`
    integer,  dimension(:), intent(in)    :: irow    !! row indices of the nonzeros of `A`
    integer,  dimension(:), intent(in)    :: icol    !! column indices of the nonzeros of `A`
    real(wp), dimension(:), intent(in)    :: val     !! the nonzeros of `A`
    real(wp),               intent(in)    :: rel_tol !! relative singularity tolerance (e.g. `1e-10`)
    integer,                intent(out)   :: istat   !! status (0 = success)

    integer(ip) :: nelem, inform, attempt
    integer(ip), dimension(:), allocatable :: iploc, iqloc, ipinv, iqinv
    real(rp),    dimension(:), allocatable :: w

    me%n  = n
    nelem = size(val)
    istat = 0
    if (n == 0) return
    me%lena = 1 + max(10*nelem, 20*me%n, 10000_ip)   ! (with room for column replacements)
    if (allocated(me%p)) deallocate(me%p, me%q, me%lenc, me%lenr, me%locc, me%locr, me%v, me%w)
    allocate(me%p(n), me%q(n), me%lenc(n), me%lenr(n), me%locc(n), me%locr(n), me%v(n), me%w(n), &
             iploc(n), iqloc(n), ipinv(n), iqinv(n), w(n))

    do attempt = 1, 3   ! (enlarging the workspace if `lu1fac` asks for it)
        if (allocated(me%a)) deallocate(me%a, me%indc, me%indr)
        allocate(me%a(me%lena), me%indc(me%lena), me%indr(me%lena))
        me%a(1:nelem)    = real(val, rp)
        me%indc(1:nelem) = irow
        me%indr(1:nelem) = icol

        me%luparm = 0
        me%luparm(1) = 6      ! nout
        me%luparm(2) = -1     ! lprint: no output
        me%luparm(3) = 5      ! maxcol
        me%luparm(6) = 1      ! TRP: threshold rook pivoting
        me%luparm(8) = 1      ! keepLU
        me%parmlu = 0.0_rp
        me%parmlu(1) = 10.0_rp                      ! Ltol1
        me%parmlu(2) = 10.0_rp                      ! Ltol2
        me%parmlu(3) = epsilon(1.0_rp)**0.8_rp      ! small
        me%parmlu(4) = epsilon(1.0_rp)**0.67_rp     ! Utol1
        me%parmlu(5) = real(rel_tol, rp)            ! Utol2
        me%parmlu(6) = 3.0_rp                       ! Uspace
        me%parmlu(7) = 0.3_rp                       ! dens1
        me%parmlu(8) = 0.5_rp                       ! dens2

        call lu1fac(me%n, me%n, nelem, me%lena, me%luparm, me%parmlu, me%a, me%indc, me%indr, &
                    me%p, me%q, me%lenc, me%lenr, me%locc, me%locr, iploc, iqloc, ipinv, iqinv, w, inform)
        if (inform /= 7) exit
        me%lena = max(2*me%lena, me%luparm(13) + 1)
    end do
    istat = int(inform)

    end subroutine lu_factorize
!*******************************************************************************

!*******************************************************************************
!>
!  solve \( A x = b \) (or \( A^T x = b \) if `transpose`) with the factors
!  from [[sqpopt_lu_type(type):factorize(bound)]].

    subroutine lu_solve(me, b, x, transpose)

    class(sqpopt_lu_type),  intent(inout) :: me
    real(wp), dimension(:), intent(in)    :: b         !! right-hand side `dimension(n)`
    real(wp), dimension(:), intent(out)   :: x         !! solution `dimension(n)`
    logical,                intent(in)    :: transpose !! solve with \( A^T \) instead of \( A \)

    integer(ip) :: inform

    if (me%n == 0) return
    if (transpose) then
        me%w(1:me%n) = real(b, rp)   ! (mode 6: `v` solves `A'v = w`; `w` is destroyed)
        call lu6sol(6_ip, me%n, me%n, me%v, me%w, me%lena, me%luparm, me%parmlu, me%a, me%indc, me%indr, &
                    me%p, me%q, me%lenc, me%lenr, me%locc, me%locr, inform)
        x = real(me%v, wp)
    else
        me%v(1:me%n) = real(b, rp)   ! (mode 5: `w` solves `A w = v`; `v` is altered)
        call lu6sol(5_ip, me%n, me%n, me%v, me%w, me%lena, me%luparm, me%parmlu, me%a, me%indc, me%indr, &
                    me%p, me%q, me%lenc, me%lenr, me%locc, me%locr, inform)
        x = real(me%w, wp)
    end if

    end subroutine lu_solve
!*******************************************************************************

!*******************************************************************************
!>
!  update the factors when column `jrep` of the matrix is replaced by the
!  sparse column (`irow`, `val`), with `LUSOL`'s `lu8rpc` (a Bartels-Golub
!  update). `istat` is `0` on success; otherwise (the new matrix appears to
!  be singular, the update seemed unstable, or the factors ran out of
!  storage) the factors can no longer be used, and the matrix should be
!  factorized again.

    subroutine lu_replace_column(me, jrep, irow, val, istat)

    class(sqpopt_lu_type),  intent(inout) :: me
    integer,                intent(in)    :: jrep  !! the column to replace
    integer,  dimension(:), intent(in)    :: irow  !! row indices of the nonzeros of the new column
    real(wp), dimension(:), intent(in)    :: val   !! the nonzeros of the new column
    integer,                intent(out)   :: istat !! status (0 = success)

    real(rp)    :: diag, vnorm
    integer(ip) :: inform, nrank0
    integer :: k

    me%v(1:me%n) = 0.0_rp
    do k = 1, size(val)
        me%v(irow(k)) = me%v(irow(k)) + real(val(k), rp)
    end do
    nrank0 = me%luparm(16)
    call lu8rpc(1_ip, 1_ip, me%n, me%n, int(jrep, ip), me%v, me%w, me%lena, me%luparm, me%parmlu, &
                me%a, me%indc, me%indr, me%p, me%q, me%lenc, me%lenr, me%locc, me%locr, inform, diag, vnorm)
    if (inform == 0 .and. me%luparm(16) == nrank0) then
        istat = 0
    else
        istat = max(1, int(abs(inform)))
    end if

    end subroutine lu_replace_column
!*******************************************************************************

    end module sqpopt_linalg_module
!*******************************************************************************
