!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  Quasi-Newton Hessian of the Lagrangian approximation.
!  Supports (damped) BFGS and SR1 updates.

    module sqpopt_hessian_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    type, public :: sqpopt_hessian_type
        !! stores and updates an approximation to the Hessian of the Lagrangian.

        integer :: n = 0  !! problem size (\( n \times n \) Hessian approximation)
        real(wp), dimension(:,:), allocatable :: b  !! current Hessian approximation `dimension(n,n)`

        contains

        procedure, public :: initialize   => hessian_initialize
        procedure, public :: update_bfgs   => hessian_update_bfgs
        procedure, public :: update_sr1    => hessian_update_sr1
        procedure, public :: reset        => hessian_reset

    end type sqpopt_hessian_type

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the Hessian approximation (typically to the identity matrix).

    subroutine hessian_initialize(me, n)

    class(sqpopt_hessian_type), intent(inout) :: me
    integer, intent(in) :: n  !! problem size

    ! TODO: implement

    end subroutine hessian_initialize
!*******************************************************************************

!*******************************************************************************
!>
!  update the Hessian approximation using the damped BFGS update formula,
!  given the step \( s = x_{k+1} - x_k \) and the change in the Lagrangian
!  gradient \( y = \nabla_x \mathcal{L}_{k+1} - \nabla_x \mathcal{L}_k \).

    subroutine hessian_update_bfgs(me, s, y)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: s  !! step vector `dimension(n)`
    real(wp), dimension(:), intent(in) :: y  !! Lagrangian gradient change `dimension(n)`

    ! TODO: implement

    end subroutine hessian_update_bfgs
!*******************************************************************************

!*******************************************************************************
!>
!  update the Hessian approximation using the symmetric rank-1 (SR1) update formula.

    subroutine hessian_update_sr1(me, s, y)

    class(sqpopt_hessian_type), intent(inout) :: me
    real(wp), dimension(:), intent(in) :: s  !! step vector `dimension(n)`
    real(wp), dimension(:), intent(in) :: y  !! Lagrangian gradient change `dimension(n)`

    ! TODO: implement

    end subroutine hessian_update_sr1
!*******************************************************************************

!*******************************************************************************
!>
!  reset the Hessian approximation back to its initial value.

    subroutine hessian_reset(me)

    class(sqpopt_hessian_type), intent(inout) :: me

    ! TODO: implement

    end subroutine hessian_reset
!*******************************************************************************

    end module sqpopt_hessian_module
!*******************************************************************************
