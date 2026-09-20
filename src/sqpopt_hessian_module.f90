!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Limited-memory quasi-Newton approximation to the Hessian of the
!  Lagrangian. Rather than storing a dense \( n \times n \) matrix, only
!  the last `max_history` step/gradient-change vector pairs \( (s,y) \)
!  are kept (`max_history` is a small constant, independent of `n`), and
!  Hessian(-inverse)-vector products are formed matrix-free using the
!  standard two-loop recursion. Supports (damped) BFGS and SR1 updates.

    module sqpopt_hessian_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    type, public :: sqpopt_hessian_type
        !! stores and updates a limited-memory approximation to the
        !! Hessian of the Lagrangian (never forms a dense `n x n` matrix).

        integer :: n           = 0  !! problem size
        integer :: max_history = 0  !! number of `(s,y)` pairs retained (independent of `n`)
        integer :: n_history   = 0  !! number of pairs currently stored (`<= max_history`)

        real(wp), dimension(:,:), allocatable :: s    !! stored step vectors `dimension(n,max_history)`
        real(wp), dimension(:,:), allocatable :: y    !! stored Lagrangian gradient-change vectors `dimension(n,max_history)`
        real(wp), dimension(:),   allocatable :: rho  !! `1/(y^T s)` for each stored pair `dimension(max_history)`
        real(wp) :: gamma = 1.0_wp  !! scaling of the initial Hessian \( H_0 = \gamma I \)

        contains

        procedure, public :: initialize             => hessian_initialize
        procedure, public :: update_bfgs             => hessian_update_bfgs
        procedure, public :: update_sr1              => hessian_update_sr1
        procedure, public :: hv_product              => hessian_vector_product
        procedure, public :: inverse_vector_product  => hessian_inverse_vector_product
        procedure, public :: reset                   => hessian_reset

    end type sqpopt_hessian_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the limited-memory Hessian approximation (equivalent to
!  \( H_0 = \gamma I \), with no `(s,y)` pairs stored).

    subroutine hessian_initialize(me, n, max_history)

    class(sqpopt_hessian_type), intent(inout) :: me
    integer, intent(in) :: n            !! problem size
    integer, intent(in) :: max_history  !! number of `(s,y)` pairs to retain

    ! TODO: implement

    end subroutine hessian_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  update the limited-memory Hessian approximation using the damped BFGS
!  update formula, given the step \( s = x_{k+1} - x_k \) and the change
!  in the Lagrangian gradient \( y = \nabla_x \mathcal{L}_{k+1} - \nabla_x \mathcal{L}_k \).
!  The oldest pair is discarded once `max_history` pairs are stored.

    subroutine hessian_update_bfgs(me, s, y)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: s  !! step vector `dimension(n)`
    real(wp), dimension(:), intent(in) :: y  !! Lagrangian gradient change `dimension(n)`

    ! TODO: implement

    end subroutine hessian_update_bfgs
!*******************************************************************************

!*******************************************************************************
!>
!  update the limited-memory Hessian approximation using the symmetric
!  rank-1 (SR1) update formula. The oldest pair is discarded once
!  `max_history` pairs are stored.

    subroutine hessian_update_sr1(me, s, y)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: s  !! step vector `dimension(n)`
    real(wp), dimension(:), intent(in) :: y  !! Lagrangian gradient change `dimension(n)`

    ! TODO: implement

    end subroutine hessian_update_sr1
!*******************************************************************************

!*******************************************************************************
!>
!  compute the matrix-free Hessian-vector product \( h_v = H v \), used
!  by the QP subproblem solver in place of an explicit dense matrix.

    subroutine hessian_vector_product(me, v, hv)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: v   !! input vector `dimension(n)`
    real(wp), dimension(:), intent(out) :: hv  !! result `dimension(n)`

    ! TODO: implement

    end subroutine hessian_vector_product
!*******************************************************************************

!*******************************************************************************
!>
!  compute the matrix-free inverse-Hessian-vector product \( d = H^{-1} v \)
!  using the standard two-loop recursion.

    subroutine hessian_inverse_vector_product(me, v, d)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in)  :: v  !! input vector `dimension(n)`
    real(wp), dimension(:), intent(out) :: d  !! result `dimension(n)`

    ! TODO: implement

    end subroutine hessian_inverse_vector_product
!*******************************************************************************

!*******************************************************************************
!>
!  reset the Hessian approximation, discarding all stored `(s,y)` pairs.

    subroutine hessian_reset(me)

    class(sqpopt_hessian_type), intent(inout) :: me

    ! TODO: implement

    end subroutine hessian_reset
!*******************************************************************************

    end module sqpopt_hessian_module
!*******************************************************************************
