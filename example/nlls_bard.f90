!*******************************************************************************
!>
!  A least-squares fit with [[sqpopt_nlls_type]] (the guide's example of
!  it): Bard's problem, the model \( y = x_1 + u/(v x_2 + w x_3) \) fitted
!  to 15 points, with the variables kept positive.
!
!  Run it with `pixi run fpm run --example nlls_bard`.

    module nlls_bard_functions

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    integer, parameter, public :: n = 3  !! number of variables
    integer, parameter, public :: l = 15 !! number of residuals (data points)
    real(wp), dimension(l), parameter :: y = [0.14_wp, 0.18_wp, 0.22_wp, 0.25_wp, 0.29_wp, 0.32_wp, 0.35_wp, 0.39_wp, &
                                              0.37_wp, 0.58_wp, 0.73_wp, 0.96_wp, 1.34_wp, 2.10_wp, 4.39_wp] !! the data

    public :: residuals, jacobian

    contains

    subroutine residuals(x, r, c, status, data)
        !! the residuals: the data minus the model
        real(wp), dimension(:), intent(in)    :: x      !! variables `dimension(n)`
        real(wp), dimension(:), intent(out)   :: r      !! residuals `dimension(l)`
        real(wp), dimension(:), intent(out)   :: c      !! constraints (none here)
        integer,                intent(inout) :: status !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop
        class(*), optional,     intent(inout) :: data   !! user data (none here)
        integer :: i
        do i = 1, l
            r(i) = y(i) - (x(1) + i/((16 - i)*x(2) + min(i, 16 - i)*x(3)))
        end do
    end subroutine residuals

    subroutine jacobian(x, rjac_val, cjac_val, accuracy, status, data)
        !! the residuals' Jacobian, by rows (in the order of the pattern given to `initialize`)
        real(wp), dimension(:), intent(in)    :: x        !! variables `dimension(n)`
        real(wp), dimension(:), intent(out)   :: rjac_val !! nonzeros of the residuals' Jacobian `dimension(n*l)`
        real(wp), dimension(:), intent(out)   :: cjac_val !! nonzeros of the constraints' Jacobian (none here)
        integer,                intent(in)    :: accuracy !! requested accuracy (unused: the derivatives are analytic)
        integer,                intent(inout) :: status   !! `0` on entry; set `> 0` if `x` can't be evaluated, or `< 0` to stop
        class(*), optional,     intent(inout) :: data     !! user data (none here)
        integer :: i
        real(wp) :: d
        do i = 1, l
            d = ((16 - i)*x(2) + min(i, 16 - i)*x(3))**2
            rjac_val(3*i-2) = -1.0_wp
            rjac_val(3*i-1) = i*(16 - i)/d
            rjac_val(3*i)   = i*min(i, 16 - i)/d
        end do
    end subroutine jacobian

    end module nlls_bard_functions
!*******************************************************************************

    program nlls_bard

    use sqpopt_kinds,        only: wp => sqpopt_module_wp
    use sqpopt_nlls_module,  only: sqpopt_nlls_type
    use nlls_bard_functions, only: n, l, residuals, jacobian

    implicit none

    type(sqpopt_nlls_type) :: nlls
    real(wp), dimension(n) :: x
    real(wp), dimension(l) :: r
    real(wp) :: sum_of_squares
    integer, dimension(n*l) :: irow, icol
    integer :: i, j, istat

    ! the pattern of the residuals' Jacobian (dense here), by rows
    do i = 1, l
        do j = 1, n
            irow(n*(i-1)+j) = i
            icol(n*(i-1)+j) = j
        end do
    end do

    call nlls%initialize(n=n, n_residuals=l, residuals=residuals, jacobian=jacobian, &
                         rjac_irow=irow, rjac_icol=icol, x_lb=[0.0_wp, 0.0_wp, 0.0_wp])
    call nlls%solve([1.0_wp, 1.0_wp, 1.0_wp], istat)
    call nlls%get_solution(x, residuals=r, sum_of_squares=sum_of_squares)

    write(*,'(A,I0,2A)')   'status         = ', istat, ': ', nlls%status_message()
    write(*,'(A,3F12.7)')  'x              = ', x
    write(*,'(A,ES14.6)')  'sum of squares = ', sum_of_squares
    write(*,'(A,ES10.2)')  'largest residual = ', maxval(abs(r))

    call nlls%destroy()

    end program nlls_bard
!*******************************************************************************
