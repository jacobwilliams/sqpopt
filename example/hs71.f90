program hs71
!! Hock-Schittkowski problem 71 example using SQPOPT
!! This is the code from the user guide worked example.

use sqpopt_module,         only: sqpopt_type
use sqpopt_problem_module, only: sqpopt_problem_type
use sqpopt_options_module, only: sqpopt_options_type
use sqpopt_kinds,          only: wp => sqpopt_module_wp

implicit none

real(wp), parameter :: big = 1.0e20_wp  !! "no bound"
type(sqpopt_type)         :: solver
type(sqpopt_problem_type) :: problem
type(sqpopt_options_type) :: options
real(wp) :: x(4), lambda(2)
integer  :: istat

! equality constraint (index 1) first, then the inequality (index 2)
call problem%set_problem_size(n=4, m=2)
call problem%set_bounds(x_lb=[1.0_wp,1.0_wp,1.0_wp,1.0_wp], &
                        x_ub=[5.0_wp,5.0_wp,5.0_wp,5.0_wp], &
                        c_lb=[40.0_wp, 25.0_wp], c_ub=[40.0_wp, big])
call problem%set_jacobian_sparsity(nnz=8, irow=[1,1,1,1,2,2,2,2], &
                                          icol=[1,2,3,4,1,2,3,4])
call problem%set_functions(fc=fc, gjac=gjac)

options%print_level = 1 ! set verbosity level

call solver%initialize(problem=problem, options=options)
call solver%solve([1.0_wp, 5.0_wp, 5.0_wp, 1.0_wp], istat)
call solver%get_solution(x, lambda)
print *, solver%status_message()
print *, 'x = ', x

contains

    subroutine fc(x, f, c, status, data)
    !! evaluate the objective and constraints at the current point `x`
    real(wp), dimension(:), intent(in)    :: x !! current point `dimension(n)`
    real(wp),               intent(out)   :: f !! objective value at `x`
    real(wp), dimension(:), intent(out)   :: c !! constraint values at `x`
    integer,                intent(inout) :: status !! status flag (input/output)
    class(*), optional,     intent(inout) :: data !! user-provided data (optional)
    f = x(1)*x(4)*(x(1)+x(2)+x(3)) + x(3)
    c = [sum(x**2), product(x)]
    end subroutine fc

    subroutine gjac(x, g, jac_val, accuracy, status, data)
    !! evaluate the gradient of the objective and the Jacobian of the constraints at the current point `x`
    real(wp), dimension(:), intent(in)    :: x !! current point `dimension(n)`
    real(wp), dimension(:), intent(out)   :: g !! gradient of the objective at `x`
    real(wp), dimension(:), intent(out)   :: jac_val !! nonzero entries of the Jacobian at `x`
    integer,                intent(in)    :: accuracy !! requested accuracy: `sqpopt_derivatives_fast` or `sqpopt_derivatives_accurate`
    integer,                intent(inout) :: status !! status flag (input/output)
    class(*), optional,     intent(inout) :: data !! user-provided data (optional)
    g = [x(4)*(2.0_wp*x(1)+x(2)+x(3)), x(1)*x(4), x(1)*x(4) + &
         1.0_wp, x(1)*(x(1)+x(2)+x(3))]
    jac_val(1:4) = 2.0_wp*x                           ! row 1: d(sum x^2)/dx
    jac_val(5:8) = [x(2)*x(3)*x(4), x(1)*x(3)*x(4), & ! row 2: d(prod x)/dx
                    x(1)*x(2)*x(4), x(1)*x(2)*x(3)]
    end subroutine gjac

end program hs71