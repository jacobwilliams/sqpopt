#!/usr/bin/env python3
"""
Convert K. Schittkowski's fixed-form Fortran 77 test problem collection
(PROB.FOR) to a free-form Fortran module (`schittkowski_problems.f90`), and
generate a modern wrapper module with the NLPQLP reference results from
TEST.DAT (`hs_problems.f90`). See the generated files' headers for what the
conversion does.

The original collection: https://klaus-schittkowski.de/test_probs_src.zip
(also at https://github.com/jacobwilliams/schittkowski-test-problems).

usage: python3 tools/convert_schittkowski.py <collection dir> test/
"""
import re, sys, os

src_dir, out_dir = sys.argv[1], sys.argv[2]

# ---------------------------------------------------------------- PROB.FOR
raw = open(os.path.join(src_dir, 'PROB.FOR'), encoding='latin-1').read().replace('\r', '')

# ---- 1. parse the fixed-form source into comments and statements.
# A statement is [label, parts], where `parts` holds the code of each of its
# lines (cols 7-72) and any comment lines interleaved between them.
items = []
for line in raw.split('\n'):
    line = line.replace('\t', ' ')
    if line.strip() == '':
        items.append(('blank',))
        continue
    if line[0] in 'Cc*!':
        items.append(('comment', ('!' + line[1:]).rstrip()))
        continue
    line = line[:72].ljust(72)
    label, cont, body = line[0:5].strip(), line[5], line[6:72]
    if cont not in (' ', '0'):
        # a continuation: attach it (and any comments since) to the last statement
        k = len(items) - 1
        while items[k][0] != 'stmt':
            k -= 1
        between = items[k+1:]
        del items[k+1:]
        items[k][2].extend([('c', b[1] if b[0] == 'comment' else '') for b in between])
        items[k][2].append(body)
    else:
        items.append(['stmt', label, [body]])

def code(st):
    """a statement's code, with the blanks removed (they are insignificant in fixed form)"""
    return ''.join(p for p in st[2] if isinstance(p, str)).replace(' ', '').upper()

def new_stmt(label, text):
    return ['stmt', label, [text], True]   # (the 4th element marks a generated statement)

# the names of all the procedures (which become module procedures):
proc_names = set(m.upper() for m in re.findall(
    r'^ {6}\s*(?:(?:DOUBLE\s*PRECISION|REAL\*8|INTEGER|LOGICAL|REAL)\s+)?(?:SUBROUTINE|FUNCTION|ENTRY)\s+(\w+)',
    raw, re.M | re.I))

# ---- 2. split into program units (from SUBROUTINE/FUNCTION to END), and modernize each
units, cur = [], []
for it in items:
    cur.append(it)
    if it[0] == 'stmt' and code(it) == 'END':
        units.append(cur)
        cur = []
trailer = cur

def top_level_split(s):
    """split `s` at top-level commas"""
    parts, depth, cur_ = [], 0, ''
    for ch in s:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        if ch == ',' and depth == 0:
            parts.append(cur_); cur_ = ''
        else:
            cur_ += ch
    parts.append(cur_)
    return parts

decl_re = re.compile(r'^(DOUBLEPRECISION|REAL\*8|INTEGER|LOGICAL|REAL)(?!FUNCTION)(.*)$')
arith_if_re = re.compile(r'^IF\((.*)\)(\d+),(\d+),(\d+)$')
do_re = re.compile(r'^DO(\d+),?([A-Z]\w*)=(.*)$')
goto_re = re.compile(r'GOTO(\d+)')
cgoto_re = re.compile(r'GOTO\(([\d,]+)\)')

stats = dict(arith_if=0, do=0, decl=0, labels=0, continues=0, entry=0)

entry_re = re.compile(r'^ENTRY(\w+)\(MODE\)$')
host_re = re.compile(r'^SUBROUTINE(\w+)\(MODE\)$')
spec_re = re.compile(r'^(SUBROUTINE|FUNCTION|DOUBLEPRECISION|REAL|INTEGER|LOGICAL|CHARACTER|IMPLICIT|COMMON|'
                     r'DIMENSION|EXTERNAL|INTRINSIC|SAVE|PARAMETER|EQUIVALENCE|DATA)')

def remove_entries(unit):
    """replace the ENTRY statements of a `TPnnn(MODE)` routine (which are
    obsolescent, and which gfortran compiles into a "master" function that
    breaks `--coverage`): the routine becomes `TPnnn_SHARED(MODE,IENTRY)`,
    each ENTRY a labeled CONTINUE that the start of its executable part
    jumps to (for `IENTRY>0`), and every entry point (including the
    original routine) a one-line wrapper around it. Returns the new unit
    and the wrappers."""
    stmts = [it for it in unit if it[0] == 'stmt']
    entries = [entry_re.match(code(it)).group(1) for it in stmts if entry_re.match(code(it))]
    if not entries:
        return unit, []
    host = host_re.match(code(stmts[0])).group(1)
    used = {it[1] for it in stmts if it[1]}
    labels = [str(l) for l in range(9001, 9100) if str(l) not in used][:len(entries)]
    shared = host + '_SHARED'
    out, k, dispatched = [], 0, False
    for it in unit:
        if it[0] == 'stmt':
            c = code(it)
            if it is stmts[0]:
                out.append(new_stmt(it[1], f'SUBROUTINE {shared}(MODE,IENTRY)'))
                continue
            if not dispatched and not spec_re.match(c):
                out.append(new_stmt('', 'INTEGER, INTENT(IN) :: IENTRY'))
                out.extend(new_stmt('', f'IF (IENTRY == {j+1}) GOTO {lab}') for j, lab in enumerate(labels))
                dispatched = True
            if entry_re.match(c):
                out.append(new_stmt(labels[k], 'CONTINUE'))
                k += 1
                stats['entry'] += 1
                continue
        out.append(it)
    wrappers = [('blank',)]
    for j, name in enumerate([host] + entries):
        wrappers += [new_stmt('', f'SUBROUTINE {name}(MODE)'), new_stmt('', f'CALL {shared}(MODE,{j})'),
                     new_stmt('', 'END'), ('blank',)]
    return out, wrappers

def modernize(unit):
    unit, wrappers = remove_entries(unit)
    out = []
    for it in unit:
        if it[0] != 'stmt':
            out.append(it); continue
        c = code(it)
        # (a) no local type declarations of the module's own functions
        #     (which would make them external procedures instead):
        m = decl_re.match(c)
        if m and '(' not in m.group(1):
            names = top_level_split(m.group(2))
            keep = [nm for nm in names if nm.split('(')[0] not in proc_names]
            if len(keep) != len(names):
                stats['decl'] += len(names) - len(keep)
                if keep:
                    out.append(new_stmt(it[1], m.group(1).replace('DOUBLEPRECISION', 'DOUBLE PRECISION') + ' ' + ','.join(keep)))
                continue
        # (b) arithmetic IF -> logical IFs (without comparing reals for equality):
        m = arith_if_re.match(c)
        if m and m.group(1).count('(') == m.group(1).count(')'):
            e, l1, l2, l3 = m.groups()
            stats['arith_if'] += 1
            if l1 == l2 == l3:
                new = ['GOTO ' + l1]
            elif l1 == l2:
                new = [f'IF (({e}) > 0) GOTO {l3}', 'GOTO ' + l1]
            elif l2 == l3:
                new = [f'IF (({e}) < 0) GOTO {l1}', 'GOTO ' + l2]
            elif l1 == l3:
                new = [f'IF (({e}) < 0 .OR. ({e}) > 0) GOTO {l1}', 'GOTO ' + l2]
            else:
                new = [f'IF (({e}) < 0) GOTO {l1}', f'IF (({e}) > 0) GOTO {l3}', 'GOTO ' + l2]
            out.append(new_stmt(it[1], new[0]))
            out.extend(new_stmt('', t) for t in new[1:])
            continue
        out.append(it)

    # labels referenced by a GOTO (plain or computed):
    refs = set()
    for it in out:
        if it[0] == 'stmt':
            c = code(it)
            refs.update(goto_re.findall(c))
            for lst in cgoto_re.findall(c):
                refs.update(lst.split(','))

    # (c) labeled DO loops -> DO ... END DO (one END DO for each loop sharing a terminal label):
    result, open_do = [], {}
    for it in out:
        if it[0] != 'stmt':
            result.append(it); continue
        c = code(it)
        m = do_re.match(c)
        if m:
            lab, var, rest = m.groups()
            open_do[lab] = open_do.get(lab, 0) + 1
            stats['do'] += 1
            result.append(new_stmt(it[1], f'DO {var}={rest}'))
            continue
        lab = it[1]
        if lab and open_do.get(lab):
            k = open_do.pop(lab)
            if not (c == 'CONTINUE' and lab not in refs):
                result.append(it)     # (still a real statement, or a GOTO target)
            result.extend(new_stmt('', 'END DO') for _ in range(k))
            continue
        result.append(it)
    assert not open_do, open_do

    # (d) drop labels that are no longer referenced, and the CONTINUEs they leave:
    final = []
    for it in result:
        if it[0] == 'stmt' and it[1] and it[1] not in refs:
            stats['labels'] += 1
            if code(it) == 'CONTINUE':
                stats['continues'] += 1
                continue
            it = [it[0], '', it[2]] + it[3:]
        final.append(it)
    return final + wrappers

def emit(items):
    lines = []
    for it in items:
        if it[0] == 'blank':
            lines.append('')
        elif it[0] == 'comment':
            lines.append(it[1])
        else:
            label = it[1]
            prefix = (label.ljust(5) + ' ') if label else '      '
            if len(it) > 3:   # generated: one line (split if very long)
                text, first = it[2][0], True
                while len(text) > 100:
                    lines.append((prefix if first else '    &') + text[:100] + '&')
                    text, first = text[100:], False
                lines.append((prefix if first else '    &') + text)
            else:             # original: keep its lines, joined with `&...&`
                first, last_code = True, None
                for p in it[2]:
                    if isinstance(p, tuple):
                        lines.append(p[1]); continue
                    if first:
                        lines.append((prefix + p).rstrip())
                        first = False
                    else:
                        lines[last_code] = lines[last_code].rstrip() + '&'
                        lines.append('    &' + p.strip())
                    last_code = len(lines) - 1
    return lines

out = []
for u in units:
    out.extend(emit(modernize(u)))
out.extend(emit(trailer))
print(stats)

header = """!*******************************************************************************
!>
!  The test problems of
!
!  * W. Hock, K. Schittkowski (1981): Test Examples for Nonlinear Programming
!    Codes, Lecture Notes in Economics and Mathematical Systems, Vol. 187,
!    Springer
!  * K. Schittkowski (1987): More Test Examples for Nonlinear Programming
!    Codes, Lecture Notes in Economics and Mathematical Systems, Vol. 282,
!    Springer
!  * K. Schittkowski (2010): An Updated Set of 306 Test Problems for
!    Nonlinear Programming with Validated Optimal Solutions, Report,
!    Department of Computer Science, University of Bayreuth
!
!  Author of the original Fortran 77 code: K. Schittkowski
!  (https://klaus-schittkowski.de/test_probs_src.zip, file `PROB.FOR`).
!
!  This file was generated from `PROB.FOR` by `tools/convert_schittkowski.py`
!  (do not edit it by hand), by a mechanical conversion:
!
!  * fixed-form source to free form (comments, continuation lines);
!  * every routine placed in this module (so they have explicit
!    interfaces), with local type declarations of the module's own
!    functions removed;
!  * labeled DO loops (including shared and non-CONTINUE terminations)
!    rewritten as `DO`/`END DO`, arithmetic IFs as logical IFs, and labels
!    that are no longer referenced removed;
!  * ENTRY statements removed: a routine with alternate entry points
!    becomes `TPnnn_SHARED(MODE,IENTRY)`, which jumps to the code at entry
!    point `IENTRY`, and each entry point a one-line wrapper around it.
!
!  The problems' code is otherwise unchanged. It still communicates through
!  the original COMMON blocks, and keeps its original implicit typing --
!  which is why this module has no `implicit none` and no module variables
!  (a module variable would be host-associated into routines that rely on
!  implicitly typed locals). Use the modern interface in
!  [[hs_problems_module]] instead of calling these routines directly.

    module schittkowski_problems_module

    contains
!*******************************************************************************
"""
footer = """
    end module schittkowski_problems_module
!*******************************************************************************
"""
with open(os.path.join(out_dir, 'schittkowski_problems.f90'), 'w') as f:
    f.write(header)
    f.write('\n'.join(out).rstrip() + '\n')
    f.write(footer)

# ------------------------------------------- the problem list, from TEST.DAT
# (TEST.DAT also has a "TP 200" line, but there is no problem 200: the
# original CONV.FOR dispatch falls through to TP1 for it, so it is skipped)
ref = []
for line in open(os.path.join(src_dir, 'TEST.DAT'), encoding='latin-1'):
    t = line.split()
    if len(t) >= 12 and t[0] == 'TP':
        ntp, n, me, m, ifail, nf, ndf, nef = map(int, t[1:9])
        if ntp == 200:
            continue
        ref.append((ntp, n, me, m, nf, ndf))

# every problem must have an entry point in PROB.FOR:
entries = set(int(k) for k in re.findall(r'^\s+(?:SUBROUTINE|ENTRY)\s+TP(\d+)\s*\(MODE\)', raw.replace('\r',''),
                                          re.M | re.I))
missing = [r[0] for r in ref if r[0] not in entries]
assert not missing, missing

ids = [r[0] for r in ref]
def fmt_list(vals, per=15):
    chunks = [', '.join(f'{v:4d}' for v in vals[i:i+per]) for i in range(0, len(vals), per)]
    return (', &\n' + ' '*12).join(chunks)

dispatch = '\n'.join(f'        case ({i:3d}); call tp{i}(mode)' for i in ids)

wrapper = f"""!*******************************************************************************
!>
!  A modern interface to the Schittkowski/Hock-Schittkowski test problem
!  collection (see [[schittkowski_problems_module]], generated from the
!  original `PROB.FOR`), for benchmarking `sqpopt`:
!
!      minimize f(x)  s.t.  c_lb <= c(x) <= c_ub,  x_lb <= x <= x_ub
!
!  [[hs_setup]] returns a problem's size, starting point, bounds, and
!  validated optimal solution; [[hs_f]], [[hs_g]], [[hs_c]], and [[hs_jac]] evaluate
!  it. The original problems' constraints are `g(x) >= 0` (inequalities,
!  first) and `g(x) = 0` (equalities, last); here they are returned in that
!  order, with bounds `[0,+inf]` and `[0,0]`. Absent variable bounds are
!  `-/+ hs_infinity`.
!
!  Also included, for comparison, are the numbers of function and gradient
!  evaluations used by Schittkowski's NLPQLP on each problem (from the
!  original `TEST.DAT`; its gradients were approximated by forward
!  differences).
!
!  This file was generated by `tools/convert_schittkowski.py` (do not edit
!  it by hand).

    module hs_problems_module

    use, intrinsic :: iso_fortran_env, only: real64
    use schittkowski_problems_module

    implicit none

    private

    integer, parameter :: dp = real64
    integer, parameter :: nmax = 101, mmax = 50  !! sizes of the original COMMON blocks

    real(dp), parameter, public :: hs_infinity = 1.0e20_dp  !! value used for absent bounds

    integer, parameter, public :: hs_n_problems = {len(ids)}  !! number of problems in the collection
    integer, dimension(hs_n_problems), parameter, public :: hs_problem_ids = [ &
            {fmt_list(ids)} ] !! the problem numbers

    ! NLPQLP's evaluation counts on each problem (from the original TEST.DAT):
    integer, dimension(hs_n_problems), parameter, public :: hs_nlpqlp_nf = [ &
            {fmt_list([r[4] for r in ref])} ] !! objective evaluations
    integer, dimension(hs_n_problems), parameter, public :: hs_nlpqlp_ndf = [ &
            {fmt_list([r[5] for r in ref])} ] !! gradient evaluations

    type, public :: hs_problem
        !! a test problem's definition and known solution
        integer  :: id = 0   !! problem number
        integer  :: n  = 0   !! number of variables
        integer  :: m  = 0   !! number of constraints
        integer  :: me = 0   !! number of equality constraints (the last `me` of the `m`)
        real(dp), dimension(:), allocatable :: x0    !! starting point `dimension(n)`
        real(dp), dimension(:), allocatable :: x_lb  !! variable lower bounds `dimension(n)`
        real(dp), dimension(:), allocatable :: x_ub  !! variable upper bounds `dimension(n)`
        real(dp), dimension(:), allocatable :: c_lb  !! constraint lower bounds `dimension(m)`
        real(dp), dimension(:), allocatable :: c_ub  !! constraint upper bounds `dimension(m)`
        real(dp) :: f_star = 0.0_dp                  !! the (validated) optimal objective value
        real(dp), dimension(:), allocatable :: x_star !! an optimal point `dimension(n)`
        logical  :: exact = .false.                  !! whether that solution is known analytically (else it
                                                     !! was computed numerically, and validated)
    end type hs_problem

    public :: hs_setup, hs_f, hs_g, hs_c, hs_jac, hs_index

    ! the original COMMON blocks (see `PROB.FOR`):
    integer  :: n_, nili, ninl, neli, nenl
    real(dp) :: x_(nmax), g_(mmax), gf_(nmax), gg_(mmax,nmax), fx_
    logical  :: index1(mmax), index2(mmax), lxl(nmax), lxu(nmax), lex
    real(dp) :: xl(nmax), xu(nmax), fex, xex(nmax)
    integer  :: nex
    common /l1/  n_, nili, ninl, neli, nenl
    common /l2/  x_
    common /l3/  g_
    common /l4/  gf_
    common /l5/  gg_
    common /l6/  fx_
    common /l9/  index1
    common /l10/ index2
    common /l11/ lxl
    common /l12/ lxu
    common /l13/ xl
    common /l14/ xu
    common /l20/ lex, nex, fex, xex

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  position of problem `id` in [[hs_problem_ids]] (0 if there is no such problem).

    pure integer function hs_index(id)

    integer, intent(in) :: id

    integer :: k

    hs_index = 0
    do k = 1, hs_n_problems
        if (hs_problem_ids(k) == id) then
            hs_index = k
            return
        end if
    end do

    end function hs_index
!*******************************************************************************

!*******************************************************************************
!>
!  call the original problem routine for problem `id` (the equivalent of
!  the original `CONV.FOR`).

    subroutine call_problem(id, mode)

    integer, intent(in) :: id, mode

    select case (id)
{dispatch}
    case default
        error stop 'hs_problems_module: no such problem'
    end select

    end subroutine call_problem
!*******************************************************************************

!*******************************************************************************
!>
!  the definition of problem `id`: size, starting point, bounds, and known
!  solution. This must be called before evaluating a problem (and again
!  when switching problems: the original code keeps constant derivative
!  elements, set here, in its COMMON blocks).

    subroutine hs_setup(id, prob)

    integer,          intent(in)  :: id
    type(hs_problem), intent(out) :: prob

    integer :: i

    if (hs_index(id) == 0) error stop 'hs_setup: no such problem'
    lxl = .false.; lxu = .false.; lex = .false.
    xl = 0.0_dp; xu = 0.0_dp; fex = 0.0_dp; xex = 0.0_dp
    ! (many problems set the constant elements of their derivatives here, in
    ! mode 1, and never again, so these must be cleared now -- and only now:)
    gf_ = 0.0_dp; gg_ = 0.0_dp
    call call_problem(id, 1)

    prob%id = id
    prob%n  = n_
    prob%me = neli + nenl
    prob%m  = nili + ninl + prob%me
    prob%x0 = x_(1:n_)
    allocate(prob%x_lb(n_), prob%x_ub(n_))
    do i = 1, n_
        prob%x_lb(i) = merge(xl(i), -hs_infinity, lxl(i))
        prob%x_ub(i) = merge(xu(i),  hs_infinity, lxu(i))
    end do
    allocate(prob%c_lb(prob%m), prob%c_ub(prob%m))
    prob%c_lb = 0.0_dp
    prob%c_ub = hs_infinity
    prob%c_ub(prob%m-prob%me+1:prob%m) = 0.0_dp
    prob%exact  = lex
    prob%f_star = fex
    prob%x_star = xex(1:n_)

    end subroutine hs_setup
!*******************************************************************************

!*******************************************************************************
!>
!  the objective of problem `id` at `x` (after [[hs_setup]] for that problem,
!  as for all the evaluation routines).

    subroutine hs_f(id, x, f)

    integer,                intent(in)  :: id
    real(dp), dimension(:), intent(in)  :: x
    real(dp),               intent(out) :: f

    x_(1:size(x)) = x
    call call_problem(id, 2)
    f = fx_

    end subroutine hs_f
!*******************************************************************************

!*******************************************************************************
!>
!  the objective gradient of problem `id` at `x`. (Some problems' original
!  code provides no gradients, or incorrect ones -- check them before use.)

    subroutine hs_g(id, x, g)

    integer,                intent(in)  :: id
    real(dp), dimension(:), intent(in)  :: x
    real(dp), dimension(:), intent(out) :: g

    x_(1:size(x)) = x
    call call_problem(id, 2)   ! (some problems compute shared terms in mode 2)
    call call_problem(id, 3)
    g = gf_(1:size(x))

    end subroutine hs_g
!*******************************************************************************

!*******************************************************************************
!>
!  the constraints of problem `id` at `x`: the inequalities first, then the
!  equalities (see [[hs_problem]]).

    subroutine hs_c(id, x, c)

    integer,                intent(in)  :: id
    real(dp), dimension(:), intent(in)  :: x
    real(dp), dimension(:), intent(out) :: c

    if (size(c) == 0) return
    x_(1:size(x)) = x
    index1 = .true.
    call call_problem(id, 4)
    c = g_(1:size(c))

    end subroutine hs_c
!*******************************************************************************

!*******************************************************************************
!>
!  the dense constraint Jacobian of problem `id` at `x` (`jac(i,j)` is the
!  derivative of constraint `i` with respect to `x(j)`). (As for [[hs_g]],
!  check before use.)

    subroutine hs_jac(id, x, jac)

    integer,                  intent(in)  :: id
    real(dp), dimension(:),   intent(in)  :: x
    real(dp), dimension(:,:), intent(out) :: jac

    if (size(jac,1) == 0) return
    x_(1:size(x)) = x
    index1 = .true.
    index2 = .true.
    call call_problem(id, 4)   ! (some problems compute shared terms in mode 4)
    call call_problem(id, 5)
    jac = gg_(1:size(jac,1), 1:size(jac,2))

    end subroutine hs_jac
!*******************************************************************************

    end module hs_problems_module
!*******************************************************************************
"""
with open(os.path.join(out_dir, 'hs_problems.f90'), 'w') as f:
    f.write(wrapper)

print(f'{len(out)} lines converted; {len(ids)} problems')
