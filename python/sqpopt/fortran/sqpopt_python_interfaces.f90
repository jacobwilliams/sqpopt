!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The interfaces of the Python callbacks of the Python bindings (see
!  `sqpopt_python.f90`), shared by `sqpopt_python` (which PRIK wraps, and
!  reads this module for) and `sqpopt_python_core` (its implementation).

    module sqpopt_python_interfaces

    use, intrinsic :: iso_fortran_env, only: real64, int32

    implicit none

    private

    abstract interface
        subroutine py_fc_func(x, f, c, status)
            !! the objective and the constraints
            import :: real64, int32
            implicit none
            real(real64), dimension(:), intent(in)    :: x      !! point `dimension(n)`
            real(real64),               intent(out)   :: f      !! objective value at `x`
            real(real64), dimension(:), intent(out)   :: c      !! constraint values at `x` `dimension(m)`
            integer(int32),             intent(inout) :: status !! `0` on entry; `> 0`: can't evaluate at `x`; `< 0`: stop
        end subroutine py_fc_func
        subroutine py_gjac_func(x, g, jac_val, accuracy, status)
            !! the objective gradient and the nonzero values of the constraint Jacobian
            import :: real64, int32
            implicit none
            real(real64), dimension(:), intent(in)    :: x        !! point `dimension(n)`
            real(real64), dimension(:), intent(out)   :: g        !! objective gradient at `x` `dimension(n)`
            real(real64), dimension(:), intent(out)   :: jac_val  !! Jacobian values at `x`, in the pattern's order
            integer(int32),             intent(in)    :: accuracy !! requested accuracy (`1` fast, `2` accurate)
            integer(int32),             intent(inout) :: status   !! `0` on entry; `> 0`: can't evaluate at `x`; `< 0`: stop
        end subroutine py_gjac_func
        subroutine py_hess_func(x, lambda, hess_val, status)
            !! the nonzero values of the Hessian of the Lagrangian
            import :: real64, int32
            implicit none
            real(real64), dimension(:), intent(in)    :: x        !! point `dimension(n)`
            real(real64), dimension(:), intent(in)    :: lambda   !! constraint multipliers `dimension(m)`
            real(real64), dimension(:), intent(out)   :: hess_val !! Hessian values at `x`, in the pattern's order
            integer(int32),             intent(inout) :: status   !! `0` on entry; `> 0`: can't evaluate at `x`; `< 0`: stop
        end subroutine py_hess_func
        subroutine py_report_func(iter, x, f, c, lambda, stop)
            !! called once per major iteration
            import :: real64, int32
            implicit none
            integer(int32),             intent(in)    :: iter   !! major iteration number
            real(real64), dimension(:), intent(in)    :: x      !! point `dimension(n)`
            real(real64),               intent(in)    :: f      !! objective value at `x`
            real(real64), dimension(:), intent(in)    :: c      !! constraint values at `x` `dimension(m)`
            real(real64), dimension(:), intent(in)    :: lambda !! constraint multipliers `dimension(m)`
            integer(int32),             intent(inout) :: stop   !! `0` on entry; set nonzero to stop the solver
        end subroutine py_report_func
    end interface

    public :: py_fc_func, py_gjac_func, py_hess_func, py_report_func

    end module sqpopt_python_interfaces
!*******************************************************************************
