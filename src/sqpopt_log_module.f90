!*******************************************************************************
!> author: Jacob Williams
!  license: MIT
!
!  The detailed log (`options%print_level >= 3`): the components that print
!  details of their work (the line search and the trust region) each hold a
!  [[sqpopt_log_type]], set from `options%print_level` and
!  `options%output_unit` at the start of each solve (so separate solver
!  objects log independently). Also the number formatting shared by the log
!  lines.

    module sqpopt_log_module

    use sqpopt_kinds, only: wp => sqpopt_module_wp
    use sqpopt_types_module, only: sqpopt_success, sqpopt_infeasible, sqpopt_qp_solve_failed, sqpopt_out_of_memory
    use, intrinsic :: iso_fortran_env, only: output_unit

    implicit none

    private

    integer, parameter, public :: sqpopt_log_detail = 3 !! the `print_level` of the detailed log

    type, public :: sqpopt_log_type
        !! where, and how much, to print (internal: set by `solve`)
        integer :: unit  = output_unit !! Fortran output unit
        integer :: level = 0           !! `options%print_level`
        contains
        procedure, public :: on  => log_on
        procedure, public :: put => log_put
    end type sqpopt_log_type

    public :: fmt_e, fmt_g, fmt_i, plural, qp_status_text

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  whether messages of the given `level` are printed.

    pure logical function log_on(me, level)

    class(sqpopt_log_type), intent(in) :: me
    integer,                intent(in) :: level !! the `print_level` a message needs

    log_on = me%level >= level

    end function log_on
!*******************************************************************************

!*******************************************************************************
!>
!  print a detail line (indented under the iteration it belongs to), if
!  messages of the given `level` are printed.

    subroutine log_put(me, level, line)

    class(sqpopt_log_type), intent(in) :: me
    integer,                intent(in) :: level !! the `print_level` the line needs
    character(len=*),       intent(in) :: line  !! the text (without the indentation, which is added)

    integer :: ios

    ! (a failed write, e.g. to a unit that isn't open, is ignored: the log
    ! must never stop the solver)
    if (me%level >= level) write(me%unit, '(A)', iostat=ios) '        . '//line

    end subroutine log_put
!*******************************************************************************

!*******************************************************************************
!>
!  a real number in the log's `ES10.3` style, without padding.

    pure function fmt_e(x) result(s)

    real(wp), intent(in) :: x !! the number
    character(len=:), allocatable :: s

    character(len=32) :: buf
    integer :: ios

    write(buf, '(ES10.3)', iostat=ios) x
    if (ios /= 0) then
        s = '****'
    else
        s = trim(adjustl(buf))
    end if

    end function fmt_e
!*******************************************************************************

!*******************************************************************************
!>
!  a real number with more digits (`ES16.9`), for objective and merit
!  values, whose changes are often in the later digits.

    pure function fmt_g(x) result(s)

    real(wp), intent(in) :: x !! the number
    character(len=:), allocatable :: s

    character(len=32) :: buf
    integer :: ios

    write(buf, '(ES16.9)', iostat=ios) x
    if (ios /= 0) then
        s = '****'
    else
        s = trim(adjustl(buf))
    end if

    end function fmt_g
!*******************************************************************************

!*******************************************************************************
!>
!  an integer, without padding.

    pure function fmt_i(i) result(s)

    integer, intent(in) :: i !! the number
    character(len=:), allocatable :: s

    character(len=32) :: buf
    integer :: ios

    write(buf, '(I0)', iostat=ios) i
    if (ios /= 0) then
        s = '****'
    else
        s = trim(adjustl(buf))
    end if

    end function fmt_i
!*******************************************************************************

!*******************************************************************************
!>
!  a short description of a QP solve's status, for the log.

    pure function qp_status_text(istat) result(s)

    integer, intent(in) :: istat !! the QP solve's status code (see [[sqpopt_types_module]])
    character(len=:), allocatable :: s

    select case (istat)
    case (sqpopt_success);         s = 'ok'
    case (sqpopt_infeasible);      s = 'inconsistent (elastic solution)'
    case (sqpopt_qp_solve_failed); s = 'failed'
    case (sqpopt_out_of_memory);   s = 'out of memory'
    case default;                  s = 'status '//fmt_i(istat)
    end select

    end function qp_status_text
!*******************************************************************************

!*******************************************************************************
!>
!  a count and a noun, e.g. `1 iteration`, `3 iterations`.

    pure function plural(n, one, many) result(s)

    integer,          intent(in) :: n    !! the count
    character(len=*), intent(in) :: one  !! the noun, singular
    character(len=*), intent(in) :: many !! the noun, plural
    character(len=:), allocatable :: s

    if (n == 1) then
        s = '1 '//one
    else
        s = fmt_i(n)//' '//many
    end if

    end function plural
!*******************************************************************************

    end module sqpopt_log_module
!*******************************************************************************
