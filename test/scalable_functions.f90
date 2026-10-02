!*******************************************************************************
!>
!  Scalable test functions: objective functions of any number of variables
!  `n`, with bounds on the variables and no other constraints, for testing
!  the solver on large problems (see `test_scalable.f90`).
!
!  The functions, their bounds, and their starting points are those of the
!  Julia package NonlinearOptimizationTestFunctions.jl
!  (https://github.com/UweAlex/NonlinearOptimizationTestFunctions.jl, MIT
!  license), which collects them from the literature (mostly Jamil & Yang,
!  "A literature survey of benchmark functions for global optimization
!  problems", 2013). The gradients were converted from that package, and the
!  Hessians were derived here; `test_scalable` checks both by finite
!  differences. Only functions that are differentiable and defined for any
!  `n` were converted.
!
!  Differences from the package:
!
!  * `qing` starts at \( x = 1 \): the package's start, \( x = 0 \), is a
!    stationary point (a maximum in every variable), where a local method
!    stops at once.
!  * `schumer_steiglitz` starts at \( x_i = 5 \sin i \): the package's start
!    is random.
!  * `trid` is evaluated term by term, not as the difference of the two
!    sums of its definition. Near its minimizer those sums are each about
!    \( n^5/30 \) and their difference about \( n^3/6 \), so the
!    difference loses \( \log_{10}(n^2/5) \) digits. With 1000 variables
!    that noise (about 0.1 in an objective of \( -1.7 \times 10^8 \)) was
!    larger than the decrease of a quasi-Newton step near the solution, and
!    the line search failed there.
!
!  Each function has a `kind`, which says what a local method can be
!  expected to find:
!
!  * `scalable_unique`: every local minimizer is a global one, so the
!    solver must reach the known minimum `f_min`.
!  * `scalable_multimodal`: there are other local minimizers, so the solver
!    must only converge to a stationary point that is better than the start.
!
!  The functions with a sparse Hessian (diagonal, tridiagonal, or in 4 by 4
!  blocks) provide it (`has_hessian`), in the solver's format: one triangle,
!  with the pattern set by [[scalable_function_setup]]. The others have a dense
!  Hessian, and are for the quasi-Newton methods only.

    module scalable_functions_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp

    implicit none

    private

    integer, parameter, public :: scalable_unique     = 1 !! every local minimizer is a global one
    integer, parameter, public :: scalable_multimodal = 2 !! there are local minimizers that are not global

    integer, parameter, public :: n_scalable_functions = 17 !! number of functions

    ! the functions' identifiers (their order is the order of the tests):
    integer, parameter :: id_sphere = 1, id_ellipsoid = 2, id_schumer_steiglitz = 3, id_trid = 4, &
                          id_rosenbrock = 5, id_dixon_price = 6, id_powell_singular = 7, id_brown = 8, &
                          id_qing = 9, id_styblinski_tang = 10, id_rastrigin = 11, id_levy = 12, &
                          id_chung_reynolds = 13, id_schwefel12 = 14, id_zakharov = 15, id_griewank = 16, &
                          id_ackley = 17

    real(wp), parameter :: pi = acos(-1.0_wp)

    type, public :: scalable_function_type
        !! one of the functions, for a number of variables `n` (see [[scalable_function_setup]])
        integer :: id = 0                     !! which function (`1` to `n_scalable_functions`)
        integer :: n  = 0                     !! number of variables
        character(len=:), allocatable :: name !! the function's name
        integer :: kind = scalable_unique     !! `scalable_unique` or `scalable_multimodal`
        logical :: has_hessian = .false.      !! whether the (sparse) Hessian is available
        real(wp) :: f_min = 0.0_wp            !! the global minimum
        real(wp), dimension(:), allocatable :: x0   !! the starting point `dimension(n)`
        real(wp), dimension(:), allocatable :: x_lb !! lower bounds `dimension(n)`
        real(wp), dimension(:), allocatable :: x_ub !! upper bounds `dimension(n)`
        integer, dimension(:), allocatable :: hess_irow !! row indices of the Hessian's pattern (if `has_hessian`)
        integer, dimension(:), allocatable :: hess_icol !! column indices of the Hessian's pattern (if `has_hessian`)
    contains
        procedure :: f => scalable_f
        procedure :: g => scalable_g
        procedure :: h => scalable_h
    end type scalable_function_type

    public :: scalable_function_setup, scalable_function_name

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  the name of function `id` (blank if there is no such function).

    pure function scalable_function_name(id) result(name)

    integer, intent(in) :: id             !! which function
    character(len=:), allocatable :: name !! its name

    select case (id)
    case (id_sphere);            name = 'sphere'
    case (id_ellipsoid);         name = 'ellipsoid'
    case (id_schumer_steiglitz); name = 'schumer_steiglitz'
    case (id_trid);              name = 'trid'
    case (id_rosenbrock);        name = 'rosenbrock'
    case (id_dixon_price);       name = 'dixon_price'
    case (id_powell_singular);   name = 'powell_singular'
    case (id_brown);             name = 'brown'
    case (id_qing);              name = 'qing'
    case (id_styblinski_tang);   name = 'styblinski_tang'
    case (id_rastrigin);         name = 'rastrigin'
    case (id_levy);              name = 'levy'
    case (id_chung_reynolds);    name = 'chung_reynolds'
    case (id_schwefel12);        name = 'schwefel12'
    case (id_zakharov);          name = 'zakharov'
    case (id_griewank);          name = 'griewank'
    case (id_ackley);            name = 'ackley'
    case default;                name = ''
    end select

    end function scalable_function_name
!*******************************************************************************

!*******************************************************************************
!>
!  set up function `id` for `n` variables: its name, kind, minimum, starting
!  point, bounds, and the pattern of its Hessian (if it has one).
!  `powell_singular` needs `n` to be a multiple of 4, and every function
!  `n >= 2`.

    subroutine scalable_function_setup(id, n, fun)

    integer,                      intent(in)  :: id  !! which function (`1` to `n_scalable_functions`)
    integer,                      intent(in)  :: n   !! number of variables
    type(scalable_function_type), intent(out) :: fun !! the function

    integer :: i
    real(wp) :: lb, ub

    fun%id   = id
    fun%n    = n
    fun%name = scalable_function_name(id)
    fun%kind = scalable_unique
    fun%has_hessian = .true.
    fun%f_min = 0.0_wp
    allocate(fun%x0(n), fun%x_lb(n), fun%x_ub(n))

    select case (id)

    case (id_sphere)             ! sum x_i^2
        lb = -5.12_wp; ub = 5.12_wp
        fun%x0 = 3.0_wp
        call diagonal_pattern()

    case (id_ellipsoid)          ! the axis-parallel hyper-ellipsoid: sum i*x_i^2
        lb = -5.12_wp; ub = 5.12_wp
        fun%x0 = 1.0_wp
        call diagonal_pattern()

    case (id_schumer_steiglitz)  ! sum x_i^4
        lb = -10.0_wp; ub = 10.0_wp
        fun%x0 = [(5.0_wp*sin(real(i, wp)), i=1,n)]
        call diagonal_pattern()

    case (id_trid)               ! sum (x_i - 1)^2 - sum x_i*x_{i-1}
        lb = -real(n, wp)**2; ub = real(n, wp)**2
        fun%x0 = 0.0_wp
        fun%f_min = -real(n, wp)*real(n + 4, wp)*real(n - 1, wp)/6.0_wp
        call tridiagonal_pattern()

    case (id_rosenbrock)         ! sum 100*(x_{i+1} - x_i^2)^2 + (x_i - 1)^2
        lb = -30.0_wp; ub = 30.0_wp
        fun%x0 = -1.2_wp
        fun%x0(n) = 1.0_wp
        fun%kind = scalable_multimodal   ! (for n >= 4 there is a local minimizer near (-1, 1, ..., 1))
        call tridiagonal_pattern()

    case (id_dixon_price)        ! (x_1 - 1)^2 + sum i*(2*x_i^2 - x_{i-1})^2
        lb = -10.0_wp; ub = 10.0_wp
        fun%x0 = 1.0_wp
        fun%kind = scalable_multimodal   ! (it has a second local minimizer)
        call tridiagonal_pattern()

    case (id_powell_singular)    ! sum over blocks of 4 variables
        if (mod(n, 4) /= 0) error stop 'scalable_function_setup: powell_singular needs n to be a multiple of 4'
        lb = -4.0_wp; ub = 5.0_wp
        fun%x0 = [(merge(3.0_wp, merge(-1.0_wp, merge(0.0_wp, 1.0_wp, mod(i,4) == 3), mod(i,4) == 2), &
                         mod(i,4) == 1), i=1,n)]
        ! the lower triangle of each block, without its two zeros:
        allocate(fun%hess_irow(2*n), fun%hess_icol(2*n))
        do i = 1, n/4
            fun%hess_irow(8*i-7:8*i) = 4*(i-1) + [1, 2, 2, 3, 3, 4, 4, 4]
            fun%hess_icol(8*i-7:8*i) = 4*(i-1) + [1, 1, 2, 2, 3, 1, 3, 4]
        end do

    case (id_brown)              ! sum (x_i^2)^(x_{i+1}^2 + 1) + (x_{i+1}^2)^(x_i^2 + 1)
        lb = -1.0_wp; ub = 4.0_wp
        fun%x0 = 1.0_wp
        call tridiagonal_pattern()

    case (id_qing)               ! sum (x_i^2 - i)^2 (every minimizer, x_i = +-sqrt(i), is global)
        lb = -500.0_wp; ub = 500.0_wp
        fun%x0 = 1.0_wp
        call diagonal_pattern()

    case (id_styblinski_tang)    ! sum (x_i^4 - 16*x_i^2 + 5*x_i)/2
        lb = -5.0_wp; ub = 5.0_wp
        fun%x0 = 0.0_wp
        fun%kind = scalable_multimodal
        fun%f_min = -39.16616570377141_wp*real(n, wp)
        call diagonal_pattern()

    case (id_rastrigin)          ! 10*n + sum x_i^2 - 10*cos(2*pi*x_i)
        lb = -5.12_wp; ub = 5.12_wp
        fun%x0 = 2.0_wp
        fun%kind = scalable_multimodal
        call diagonal_pattern()

    case (id_levy)               ! (see `scalable_f`)
        lb = -10.0_wp; ub = 10.0_wp
        fun%x0 = 0.0_wp
        fun%kind = scalable_multimodal
        call diagonal_pattern()

    case (id_chung_reynolds)     ! (sum x_i^2)^2
        lb = -100.0_wp; ub = 100.0_wp
        fun%x0 = 1.0_wp
        fun%has_hessian = .false.

    case (id_schwefel12)         ! sum_i (sum_{j<=i} x_j)^2
        lb = -100.0_wp; ub = 100.0_wp
        fun%x0 = 50.0_wp
        fun%has_hessian = .false.

    case (id_zakharov)           ! sum x_i^2 + s^2 + s^4, s = sum i*x_i/2
        lb = -5.0_wp; ub = 10.0_wp
        fun%x0 = 1.0_wp
        fun%has_hessian = .false.

    case (id_griewank)           ! sum x_i^2/4000 - prod cos(x_i/sqrt(i)) + 1
        lb = -600.0_wp; ub = 600.0_wp
        fun%x0 = 300.0_wp
        fun%kind = scalable_multimodal
        fun%has_hessian = .false.

    case (id_ackley)             ! (see `scalable_f`)
        lb = -32.768_wp; ub = 32.768_wp
        fun%x0 = 16.0_wp
        fun%kind = scalable_multimodal
        fun%has_hessian = .false.

    case default
        error stop 'scalable_function_setup: no such function'

    end select

    fun%x_lb = lb
    fun%x_ub = ub
    if (.not. fun%has_hessian) allocate(fun%hess_irow(0), fun%hess_icol(0))

    contains

        subroutine diagonal_pattern()
        !! the pattern of a diagonal Hessian
        fun%hess_irow = [(i, i=1,n)]
        fun%hess_icol = [(i, i=1,n)]
        end subroutine diagonal_pattern

        subroutine tridiagonal_pattern()
        !! the pattern of a tridiagonal Hessian: the diagonal, then the subdiagonal
        fun%hess_irow = [(i, i=1,n), (i+1, i=1,n-1)]
        fun%hess_icol = [(i, i=1,n), (i,   i=1,n-1)]
        end subroutine tridiagonal_pattern

    end subroutine scalable_function_setup
!*******************************************************************************

!*******************************************************************************
!>
!  the value of the function at `x`.

    function scalable_f(me, x) result(f)

    class(scalable_function_type), intent(in) :: me
    real(wp), dimension(:),        intent(in) :: x !! the point `dimension(n)`
    real(wp) :: f                                  !! the function's value

    integer :: i, n
    real(wp) :: s, s2, w, p

    n = me%n
    select case (me%id)

    case (id_sphere)
        f = sum(x**2)

    case (id_ellipsoid)
        f = 0.0_wp
        do i = 1, n
            f = f + real(i, wp)*x(i)**2
        end do

    case (id_schumer_steiglitz)
        f = sum(x**4)

    case (id_trid)
        ! (term by term, as x_i*(x_i - x_{i-1}) - 2*x_i + 1: the two sums of the definition are
        ! each about n^5/30 near the minimizer, and their difference, about n^3/6, would lose
        ! log10(n^2/5) digits, which is enough to stop a line search before it has converged)
        f = x(1)*x(1) - 2.0_wp*x(1) + 1.0_wp
        do i = 2, n
            f = f + (x(i)*(x(i) - x(i-1)) - 2.0_wp*x(i) + 1.0_wp)
        end do

    case (id_rosenbrock)
        f = sum(100.0_wp*(x(2:n) - x(1:n-1)**2)**2 + (x(1:n-1) - 1.0_wp)**2)

    case (id_dixon_price)
        f = (x(1) - 1.0_wp)**2
        do i = 2, n
            f = f + real(i, wp)*(2.0_wp*x(i)**2 - x(i-1))**2
        end do

    case (id_powell_singular)
        f = 0.0_wp
        do i = 1, n, 4
            f = f + (x(i) + 10.0_wp*x(i+1))**2 + 5.0_wp*(x(i+2) - x(i+3))**2 + &
                    (x(i+1) - 2.0_wp*x(i+2))**4 + 10.0_wp*(x(i) - x(i+3))**4
        end do

    case (id_brown)
        f = 0.0_wp
        do i = 1, n-1
            f = f + (x(i)**2)**(x(i+1)**2 + 1.0_wp) + (x(i+1)**2)**(x(i)**2 + 1.0_wp)
        end do

    case (id_qing)
        f = 0.0_wp
        do i = 1, n
            f = f + (x(i)**2 - real(i, wp))**2
        end do

    case (id_styblinski_tang)
        f = 0.5_wp*sum(x**4 - 16.0_wp*x**2 + 5.0_wp*x)

    case (id_rastrigin)
        f = 10.0_wp*real(n, wp) + sum(x**2 - 10.0_wp*cos(2.0_wp*pi*x))

    case (id_levy)
        ! with w_i = 1 + (x_i - 1)/4:
        ! sin^2(pi*w_1) + sum_{i<n} (w_i - 1)^2 (1 + 10 sin^2(pi*w_i + 1)) + (w_n - 1)^2 (1 + sin^2(2*pi*w_n))
        w = 1.0_wp + (x(1) - 1.0_wp)/4.0_wp
        f = sin(pi*w)**2
        do i = 1, n-1
            w = 1.0_wp + (x(i) - 1.0_wp)/4.0_wp
            f = f + (w - 1.0_wp)**2*(1.0_wp + 10.0_wp*sin(pi*w + 1.0_wp)**2)
        end do
        w = 1.0_wp + (x(n) - 1.0_wp)/4.0_wp
        f = f + (w - 1.0_wp)**2*(1.0_wp + sin(2.0_wp*pi*w)**2)

    case (id_chung_reynolds)
        f = sum(x**2)**2

    case (id_schwefel12)
        f = 0.0_wp
        s = 0.0_wp
        do i = 1, n
            s = s + x(i)
            f = f + s**2
        end do

    case (id_zakharov)
        s = 0.0_wp
        do i = 1, n
            s = s + 0.5_wp*real(i, wp)*x(i)
        end do
        f = sum(x**2) + s**2 + s**4

    case (id_griewank)
        p = 1.0_wp
        do i = 1, n
            p = p*cos(x(i)/sqrt(real(i, wp)))
        end do
        f = sum(x**2)/4000.0_wp - p + 1.0_wp

    case (id_ackley)
        ! -20 exp(-0.2 sqrt(sum x_i^2/n)) - exp(sum cos(2*pi*x_i)/n) + 20 + e
        s  = sum(x**2)/real(n, wp)
        s2 = sum(cos(2.0_wp*pi*x))/real(n, wp)
        f = -20.0_wp*exp(-0.2_wp*sqrt(s)) - exp(s2) + 20.0_wp + exp(1.0_wp)

    case default
        error stop 'scalable_f: no such function'

    end select

    end function scalable_f
!*******************************************************************************

!*******************************************************************************
!>
!  the gradient of the function at `x`.

    subroutine scalable_g(me, x, g)

    class(scalable_function_type), intent(in)  :: me
    real(wp), dimension(:),        intent(in)  :: x !! the point `dimension(n)`
    real(wp), dimension(:),        intent(out) :: g !! the gradient `dimension(n)`

    integer :: i, n
    real(wp) :: s, s2, w, t, r, a, b
    real(wp), dimension(:), allocatable :: c, before, after

    n = me%n
    select case (me%id)

    case (id_sphere)
        g = 2.0_wp*x

    case (id_ellipsoid)
        g = [(2.0_wp*real(i, wp)*x(i), i=1,n)]

    case (id_schumer_steiglitz)
        g = 4.0_wp*x**3

    case (id_trid)
        g = 2.0_wp*(x - 1.0_wp)
        g(2:n)   = g(2:n)   - x(1:n-1)
        g(1:n-1) = g(1:n-1) - x(2:n)

    case (id_rosenbrock)
        g = 0.0_wp
        do i = 1, n-1
            t = x(i+1) - x(i)**2
            g(i)   = g(i) - 400.0_wp*t*x(i) + 2.0_wp*(x(i) - 1.0_wp)
            g(i+1) = g(i+1) + 200.0_wp*t
        end do

    case (id_dixon_price)
        g = 0.0_wp
        g(1) = 2.0_wp*(x(1) - 1.0_wp)
        do i = 2, n
            t = 2.0_wp*x(i)**2 - x(i-1)
            g(i)   = g(i) + 8.0_wp*real(i, wp)*t*x(i)
            g(i-1) = g(i-1) - 2.0_wp*real(i, wp)*t
        end do

    case (id_powell_singular)
        do i = 1, n, 4
            a = x(i+1) - 2.0_wp*x(i+2)
            b = x(i) - x(i+3)
            g(i)   =  2.0_wp*(x(i) + 10.0_wp*x(i+1)) + 40.0_wp*b**3
            g(i+1) = 20.0_wp*(x(i) + 10.0_wp*x(i+1)) + 4.0_wp*a**3
            g(i+2) = 10.0_wp*(x(i+2) - x(i+3)) - 8.0_wp*a**3
            g(i+3) = -10.0_wp*(x(i+2) - x(i+3)) - 40.0_wp*b**3
        end do

    case (id_brown)
        g = 0.0_wp
        do i = 1, n-1
            call brown_term(x(i), x(i+1), g(i), g(i+1))
            call brown_term(x(i+1), x(i), g(i+1), g(i))
        end do

    case (id_qing)
        g = [(4.0_wp*x(i)*(x(i)**2 - real(i, wp)), i=1,n)]

    case (id_styblinski_tang)
        g = 0.5_wp*(4.0_wp*x**3 - 32.0_wp*x + 5.0_wp)

    case (id_rastrigin)
        g = 2.0_wp*x + 20.0_wp*pi*sin(2.0_wp*pi*x)

    case (id_levy)
        g = 0.0_wp
        w = 1.0_wp + (x(1) - 1.0_wp)/4.0_wp
        g(1) = pi*sin(2.0_wp*pi*w)/4.0_wp
        do i = 1, n-1
            w = 1.0_wp + (x(i) - 1.0_wp)/4.0_wp
            t = pi*w + 1.0_wp
            g(i) = g(i) + (2.0_wp*(w - 1.0_wp)*(1.0_wp + 10.0_wp*sin(t)**2) + &
                           10.0_wp*pi*(w - 1.0_wp)**2*sin(2.0_wp*t))/4.0_wp
        end do
        w = 1.0_wp + (x(n) - 1.0_wp)/4.0_wp
        g(n) = g(n) + (2.0_wp*(w - 1.0_wp)*(1.0_wp + sin(2.0_wp*pi*w)**2) + &
                       2.0_wp*pi*(w - 1.0_wp)**2*sin(4.0_wp*pi*w))/4.0_wp

    case (id_chung_reynolds)
        g = 4.0_wp*sum(x**2)*x

    case (id_schwefel12)
        ! with the partial sums s_i: g_k = 2*sum_{i>=k} s_i
        allocate(c(n))
        s = 0.0_wp
        do i = 1, n
            s = s + x(i)
            c(i) = s
        end do
        s = 0.0_wp
        do i = n, 1, -1
            s = s + c(i)
            g(i) = 2.0_wp*s
        end do

    case (id_zakharov)
        s = 0.0_wp
        do i = 1, n
            s = s + 0.5_wp*real(i, wp)*x(i)
        end do
        g = [(2.0_wp*x(i) + (2.0_wp*s + 4.0_wp*s**3)*0.5_wp*real(i, wp), i=1,n)]

    case (id_griewank)
        ! (the product without factor i, from the products of the factors before and after it)
        allocate(c(n), before(n), after(n))
        c = [(cos(x(i)/sqrt(real(i, wp))), i=1,n)]
        before(1) = 1.0_wp
        do i = 2, n
            before(i) = before(i-1)*c(i-1)
        end do
        after(n) = 1.0_wp
        do i = n-1, 1, -1
            after(i) = after(i+1)*c(i+1)
        end do
        do i = 1, n
            r = sqrt(real(i, wp))
            g(i) = x(i)/2000.0_wp + before(i)*after(i)*sin(x(i)/r)/r
        end do

    case (id_ackley)
        s  = sum(x**2)/real(n, wp)
        s2 = sum(cos(2.0_wp*pi*x))/real(n, wp)
        r  = sqrt(s)
        g = exp(s2)*2.0_wp*pi*sin(2.0_wp*pi*x)/real(n, wp)
        ! (the first term is not differentiable at x = 0, its minimizer)
        if (r > 0.0_wp) g = g + 4.0_wp*exp(-0.2_wp*r)*x/(real(n, wp)*r)

    case default
        error stop 'scalable_g: no such function'

    end select

    contains

        subroutine brown_term(a, b, ga, gb)
        !! add the derivatives of \( (a^2)^{b^2+1} \) to `ga` and `gb`
        real(wp), intent(in)    :: a  !! the base's variable
        real(wp), intent(in)    :: b  !! the exponent's variable
        real(wp), intent(inout) :: ga !! the derivative with respect to `a`
        real(wp), intent(inout) :: gb !! the derivative with respect to `b`
        real(wp) :: p
        if (a == 0.0_wp) return   ! (the term and its derivatives are zero)
        p = b**2 + 1.0_wp
        ga = ga + 2.0_wp*p*a*(a**2)**(p - 1.0_wp)
        gb = gb + (a**2)**p*log(a**2)*2.0_wp*b
        end subroutine brown_term

    end subroutine scalable_g
!*******************************************************************************

!*******************************************************************************
!>
!  the Hessian of the function at `x`, in the pattern of
!  [[scalable_function_setup]] (only for the functions with `has_hessian`).

    subroutine scalable_h(me, x, h)

    class(scalable_function_type), intent(in)  :: me
    real(wp), dimension(:),        intent(in)  :: x !! the point `dimension(n)`
    real(wp), dimension(:),        intent(out) :: h !! the Hessian's values `dimension(size(hess_irow))`

    integer :: i, k, n
    real(wp) :: w, t, a, b, u

    n = me%n
    h = 0.0_wp
    select case (me%id)

    case (id_sphere)
        h = 2.0_wp

    case (id_ellipsoid)
        h = [(2.0_wp*real(i, wp), i=1,n)]

    case (id_schumer_steiglitz)
        h = 12.0_wp*x**2

    case (id_trid)
        h(1:n)  = 2.0_wp
        h(n+1:) = -1.0_wp

    case (id_rosenbrock)
        do i = 1, n-1
            h(i)   = h(i) + 1200.0_wp*x(i)**2 - 400.0_wp*x(i+1) + 2.0_wp
            h(i+1) = h(i+1) + 200.0_wp
            h(n+i) = -400.0_wp*x(i)
        end do

    case (id_dixon_price)
        h(1) = 2.0_wp
        do i = 2, n
            t = 2.0_wp*x(i)**2 - x(i-1)
            h(i)     = h(i) + 8.0_wp*real(i, wp)*(t + 4.0_wp*x(i)**2)
            h(i-1)   = h(i-1) + 2.0_wp*real(i, wp)
            h(n+i-1) = -8.0_wp*real(i, wp)*x(i)
        end do

    case (id_powell_singular)
        do i = 1, n, 4
            k = 2*(i - 1)
            a = 12.0_wp*(x(i+1) - 2.0_wp*x(i+2))**2
            b = 120.0_wp*(x(i) - x(i+3))**2
            h(k+1) = 2.0_wp + b          ! (1,1)
            h(k+2) = 20.0_wp             ! (2,1)
            h(k+3) = 200.0_wp + a        ! (2,2)
            h(k+4) = -2.0_wp*a           ! (3,2)
            h(k+5) = 10.0_wp + 4.0_wp*a  ! (3,3)
            h(k+6) = -b                  ! (4,1)
            h(k+7) = -10.0_wp            ! (4,3)
            h(k+8) = 10.0_wp + b         ! (4,4)
        end do

    case (id_brown)
        do i = 1, n-1
            call brown_term(x(i), x(i+1), h(i), h(i+1), h(n+i))
            call brown_term(x(i+1), x(i), h(i+1), h(i), h(n+i))
        end do

    case (id_qing)
        h = [(12.0_wp*x(i)**2 - 4.0_wp*real(i, wp), i=1,n)]

    case (id_styblinski_tang)
        h = 6.0_wp*x**2 - 16.0_wp

    case (id_rastrigin)
        h = 2.0_wp + 40.0_wp*pi**2*cos(2.0_wp*pi*x)

    case (id_levy)
        w = 1.0_wp + (x(1) - 1.0_wp)/4.0_wp
        h(1) = 2.0_wp*pi**2*cos(2.0_wp*pi*w)/16.0_wp
        do i = 1, n-1
            w = 1.0_wp + (x(i) - 1.0_wp)/4.0_wp
            t = pi*w + 1.0_wp
            u = w - 1.0_wp
            h(i) = h(i) + (2.0_wp*(1.0_wp + 10.0_wp*sin(t)**2) + 40.0_wp*pi*u*sin(2.0_wp*t) + &
                           20.0_wp*pi**2*u**2*cos(2.0_wp*t))/16.0_wp
        end do
        w = 1.0_wp + (x(n) - 1.0_wp)/4.0_wp
        u = w - 1.0_wp
        h(n) = h(n) + (2.0_wp*(1.0_wp + sin(2.0_wp*pi*w)**2) + 8.0_wp*pi*u*sin(4.0_wp*pi*w) + &
                       8.0_wp*pi**2*u**2*cos(4.0_wp*pi*w))/16.0_wp

    case default
        error stop 'scalable_h: the function has no Hessian'

    end select

    contains

        subroutine brown_term(a, b, haa, hbb, hab)
        !! add the second derivatives of \( (a^2)^{b^2+1} \) to `haa`, `hbb`, and `hab`
        real(wp), intent(in)    :: a   !! the base's variable
        real(wp), intent(in)    :: b   !! the exponent's variable
        real(wp), intent(inout) :: haa !! the second derivative with respect to `a`
        real(wp), intent(inout) :: hbb !! the second derivative with respect to `b`
        real(wp), intent(inout) :: hab !! the mixed second derivative
        real(wp) :: p, la, val
        p = b**2 + 1.0_wp
        if (a == 0.0_wp) then
            ! (only the term a^2, for b = 0, has a second derivative that isn't zero)
            if (b == 0.0_wp) haa = haa + 2.0_wp
            return
        end if
        la  = log(a**2)
        val = (a**2)**p
        haa = haa + 2.0_wp*p*(2.0_wp*p - 1.0_wp)*(a**2)**(p - 1.0_wp)
        hbb = hbb + val*(2.0_wp*la + (2.0_wp*b*la)**2)
        hab = hab + 4.0_wp*a*b*(a**2)**(p - 1.0_wp)*(1.0_wp + p*la)
        end subroutine brown_term

    end subroutine scalable_h
!*******************************************************************************

    end module scalable_functions_module
