"""SQPOPT on the CUTEst test problems, next to scipy's SLSQP and trust-constr, through PyCUTEst and the
Python bindings (see plan/OPENSQP_COMPARISON.md, whose paper ran such a comparison on 575 problems).

It is a developer's tool, not a test: it needs the CUTEst library, the SIF decoder, and the SIF problem
files, which `setup` puts in the directory `cutest/` (not in the repository; about 4 GB, most of it the
problem files). Run from the repository root, in the pixi environment:

    pixi run build-python     # the bindings' extension, if not built yet
    pixi run cutest setup     # once: download what is missing
    pixi run cutest list      # the problems selected (by --max-n and --max-m)
    pixi run cutest run       # solve them (continues an interrupted run)
    pixi run cutest report    # the table of the results

`run` solves each problem in a process of its own, with a time limit (a crash, or a solver that doesn't
return, then loses only that problem), and adds a line per problem and solver to the results file
(JSON lines), so that it can be interrupted and continued. Several problems are solved at once (`--jobs`),
which doesn't change the results, but makes the times recorded depend on the machine's load: use
`--jobs=1` for times that will be written down. Options (after the command):

    --max-n=N, --max-m=M     select the problems with at most N variables and M constraints (default 100, 100)
    --problems=A,B,...       only these problems
    --no-overdetermined      leave out the problems with more constraints than variables (144 of the default
                             664; most are fits that can't be feasible). For quick runs: the numbers of the
                             Performance page are of the whole set. `list` and `run` go by the sizes in
                             CUTEst's classification, which differ from the decoded problem's for a few
                             problems; `report`, which leaves them out of the table, by the sizes recorded.
    --jobs=J                 (run) the number of problems solved at once (default: the number of processors)
    --solvers=a,b,...        of sqpopt, slsqp, trust-constr (default: all three)
    --maxiter=K              iteration limit of every solver (default 250, as in the paper)
    --timeout=S              seconds allowed for a problem, all solvers together (default 120)
    --results=FILE           the results file (default cutest/results.jsonl)
    --sqpopt=name=value,...  SQPOPT options (e.g. ktol=1.22e-4,ctol=2e-6), under the solver name sqpopt[...]
    --redo                   solve again the problems that are already in the results file
    --failures[=SOLVER]      (report) also list the problems each solver (or that one) didn't solve

The CUTEst library is the directory `CUTEst_binaries*` in the repository root (a binary release of
https://github.com/ralna/CUTEst), or `$CUTEST`, or else the one `setup` downloads. PyCUTEst assumes
Homebrew on macOS (it runs `brew --prefix`, and looks for the Fortran run-time library under it), so a
stand-in `brew` is put on the path that points it at the pixi environment's library instead.
"""

from __future__ import annotations

import concurrent.futures
import json
import os
import pathlib
import platform
import subprocess
import sys
import tarfile
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
HOME = ROOT / 'cutest'
INFINITY = 1.0e20            # CUTEst's "no bound"
RELEASES = {'cutest': ('ralna/CUTEst', 'v2.7.1', 'CUTEst_binaries.v2.7.1'),
            'sifdecode': ('ralna/SIFDecode', 'v3.1.1', 'SIFDecode_binaries.v3.1.1')}
SIF_REPOSITORY = 'https://bitbucket.org/optrove/sif.git'
SOLVERS = ('sqpopt', 'slsqp', 'trust-constr')


def release_platform() -> str:
    """the platform's name in the names of the binary releases"""
    machine = {'arm64': 'aarch64', 'aarch64': 'aarch64', 'x86_64': 'x86_64'}[platform.machine()]
    system = {'Darwin': 'apple-darwin', 'Linux': 'linux-gnu'}[platform.system()]
    return f'{machine}-{system}-libgfortran5'


def cutest_directory() -> pathlib.Path | None:
    """where the CUTEst library is: `$CUTEST`, a `CUTEst_binaries*` in the repository root, or `cutest/cutest`"""
    candidates = [pathlib.Path(os.environ['CUTEST'])] if 'CUTEST' in os.environ else []
    candidates += sorted(ROOT.glob('CUTEst_binaries*')) + [HOME / 'cutest']
    for d in candidates:
        if (d / 'lib' / 'libcutest_double.a').is_file():
            return d
    return None


def setup_environment() -> list[str]:
    """set the environment PyCUTEst needs (before it is imported); returns what is missing"""
    prefix = pathlib.Path(os.environ.get('CONDA_PREFIX', sys.prefix))
    missing = []
    cutest = cutest_directory()
    if cutest is None:
        missing.append('the CUTEst library')
    else:
        os.environ['CUTEST'] = str(cutest)
    if not (HOME / 'sifdecode' / 'bin' / 'sifdecoder').is_file():
        missing.append('the SIF decoder')
    if not any((HOME / 'sif').glob('*.SIF')):
        missing.append('the SIF problem files')
    os.environ['SIFDECODE'] = str(HOME / 'sifdecode')
    os.environ['MASTSIF'] = str(HOME / 'sif')
    os.environ['PYCUTEST_CACHE'] = str(HOME / 'cache')
    (HOME / 'cache').mkdir(parents=True, exist_ok=True)
    # (the binary releases look for the Fortran run-time library on their run path)
    os.environ['DYLD_FALLBACK_LIBRARY_PATH'] = str(prefix / 'lib')
    os.environ['LD_LIBRARY_PATH'] = str(prefix / 'lib') + os.pathsep + os.environ.get('LD_LIBRARY_PATH', '')
    if platform.system() == 'Darwin':
        # the stand-in for Homebrew: `brew --prefix` gives a directory whose `Cellar/gcc/*/lib/gcc/*` is
        # the pixi environment's library directory
        shim = HOME / 'shim'
        gcc = shim / 'prefix' / 'Cellar' / 'gcc' / 'pixi' / 'lib' / 'gcc'
        gcc.mkdir(parents=True, exist_ok=True)
        link = gcc / 'current'
        if link.is_symlink() and link.resolve() != (prefix / 'lib').resolve():
            link.unlink()
        if not link.exists():
            link.symlink_to(prefix / 'lib')
        # (written only if it isn't there as it should be, and put in place whole: `run`'s child processes
        # come through here too, while others of them are running it)
        brew = shim / 'brew'
        text = f'#!/bin/sh\n# stands in for Homebrew: see tools/cutest_benchmark.py\necho "{shim / "prefix"}"\n'
        if not (brew.is_file() and brew.read_text() == text and os.access(brew, os.X_OK)):
            new = shim / f'brew.{os.getpid()}'
            new.write_text(text)
            new.chmod(0o755)
            new.replace(brew)
        os.environ['PATH'] = str(shim) + os.pathsep + os.environ['PATH']
    return missing


def setup():
    """download what is missing into `cutest/`"""
    HOME.mkdir(exist_ok=True)
    for what, target in (('cutest', HOME / 'cutest'), ('sifdecode', HOME / 'sifdecode')):
        if what == 'cutest' and cutest_directory() is not None:
            print(f'CUTEst library: {cutest_directory()}')
            continue
        if what == 'sifdecode' and (target / 'bin' / 'sifdecoder').is_file():
            print(f'SIF decoder: {target}')
            continue
        repository, tag, name = RELEASES[what]
        url = f'https://github.com/{repository}/releases/download/{tag}/{name}.{release_platform()}.tar.gz'
        print(f'downloading {url}')
        archive = HOME / f'{what}.tar.gz'
        urllib.request.urlretrieve(url, archive)
        target.mkdir(exist_ok=True)
        with tarfile.open(archive) as tar:
            tar.extractall(target)
        archive.unlink()
    if any((HOME / 'sif').glob('*.SIF')):
        print(f'SIF problem files: {HOME / "sif"}')
    else:
        print(f'cloning {SIF_REPOSITORY} (about 4 GB)')
        subprocess.run(['git', 'clone', '--depth', '1', SIF_REPOSITORY, str(HOME / 'sif')], check=True)
    missing = setup_environment()
    print('missing: ' + ', '.join(missing) if missing else 'ready')


def select(args: dict) -> list[str]:
    """the names of the problems to solve, sorted: those given, or those with fixed sizes within the limits
    (a problem whose size is a parameter of its SIF file is not selected: its default size isn't in
    CUTEst's classification), without those with more constraints than variables if asked to leave them
    out"""
    if args['problems']:
        return sorted(args['problems'].split(','))
    import pycutest
    names = []
    for name in pycutest.find_problems():
        properties = pycutest.problem_properties(name)
        n, m = properties['n'], properties['m']
        if isinstance(n, int) and isinstance(m, int) and 1 <= n <= int(args['max-n']) and m <= int(args['max-m']):
            if not (args['no-overdetermined'] and m > n):
                names.append(name)
    return sorted(names)


class Functions:
    """a CUTEst problem's functions, with counts of their evaluations"""

    def __init__(self, problem):
        self.p = problem
        self.n_f = self.n_g = self.n_c = self.n_jac = 0

    def f(self, x):
        self.n_f += 1
        return self.p.obj(x)

    def g(self, x):
        self.n_g += 1
        return self.p.obj(x, gradient=True)[1]

    def c(self, x):
        self.n_c += 1
        return self.p.cons(x)

    def jac(self, x):
        self.n_jac += 1
        return self.p.cons(x, gradient=True)[1]

    def violation(self, x) -> float:
        """the largest violation of a bound or a constraint at `x`"""
        import numpy as np
        v = max(0.0, float(np.max(self.p.bl - x)), float(np.max(x - self.p.bu)))
        if self.p.m > 0:
            c = self.p.cons(x)
            v = max(v, float(np.max(self.p.cl - c)), float(np.max(c - self.p.cu)))
        return v


def solve_one(name: str, solvers: list[str], maxiter: int, sqpopt_options: dict, label: str):
    """solve one problem with each solver, printing a JSON line for each (the child process of `run`)"""
    import numpy as np
    import pycutest
    import scipy.optimize
    sys.path.insert(0, str(ROOT / 'python'))
    import sqpopt

    p = pycutest.import_problem(name)
    lb = np.where(p.bl <= -INFINITY, -np.inf, p.bl)
    ub = np.where(p.bu >= INFINITY, np.inf, p.bu)
    if p.m > 0:
        cl = np.where(p.cl <= -INFINITY, -np.inf, p.cl)
        cu = np.where(p.cu >= INFINITY, np.inf, p.cu)
    for solver in solvers:
        fun = Functions(p)
        record = dict(problem=name, solver=label if solver == 'sqpopt' else solver, n=int(p.n), m=int(p.m))
        t0 = time.perf_counter()
        try:
            if solver == 'sqpopt':
                cons = [sqpopt.NonlinearConstraint(fun.c, cl, cu, jac=fun.jac)] if p.m > 0 else []
                r = sqpopt.minimize(fun.f, p.x0, jac=fun.g, bounds=list(zip(lb, ub)), constraints=cons,
                                    options=dict(max_iter=maxiter, **sqpopt_options))
            else:
                cons = [scipy.optimize.NonlinearConstraint(fun.c, cl, cu, jac=fun.jac)] if p.m > 0 else []
                bounds = scipy.optimize.Bounds(lb, ub)
                x0 = np.clip(p.x0, lb, ub)
                if solver == 'slsqp':
                    r = scipy.optimize.minimize(fun.f, x0, jac=fun.g, bounds=bounds, constraints=cons,
                                                method='SLSQP', options=dict(maxiter=maxiter, ftol=1.0e-6))
                else:
                    r = scipy.optimize.minimize(fun.f, x0, jac=fun.g, bounds=bounds, constraints=cons,
                                                method='trust-constr',
                                                options=dict(maxiter=maxiter, gtol=2.0e-5, xtol=2.0e-100))
            x = np.asarray(r.x, dtype=float)
            record.update(success=bool(r.success), status=int(r.status), message=str(r.message)[:80],
                          f=float(p.obj(x)), violation=fun.violation(x), iterations=int(r.nit))
        except Exception as e:   # (a solver that raises is a failure of that solver on this problem)
            record.update(success=False, status=-1, message=f'{type(e).__name__}: {e}'[:80])
        record.update(time=time.perf_counter() - t0, n_f=fun.n_f, n_g=fun.n_g, n_c=fun.n_c, n_jac=fun.n_jac)
        print(json.dumps(record), flush=True)


def solve_in_child(name: str, todo: list[str], args: dict) -> tuple[list[dict], str]:
    """solve one problem in a child process: the records its solvers printed, and a note if it crashed or
    reached the time limit"""
    command = [sys.executable, __file__, 'solve', f'--problems={name}', f'--solvers={",".join(todo)}',
               f'--maxiter={args["maxiter"]}', f'--sqpopt={args["sqpopt"]}']
    lines, note = [], ''
    try:
        out = subprocess.run(command, capture_output=True, text=True, timeout=float(args['timeout']))
        lines = out.stdout.splitlines()
        if out.returncode != 0:
            note = f'exit code {out.returncode}: ' + out.stderr.strip().splitlines()[-1][:80] if out.stderr.strip() \
                else f'exit code {out.returncode}'
    except subprocess.TimeoutExpired as e:
        lines = (e.stdout.decode() if isinstance(e.stdout, bytes) else e.stdout or '').splitlines()
        note = 'time limit'
    return [json.loads(line) for line in lines if line.startswith('{')], note


def run(args: dict):
    """solve the selected problems, each in a child process and `--jobs` of them at once, adding to the
    results file (in the order they finish)"""
    results = pathlib.Path(args['results'])
    solvers = args['solvers'].split(',')
    label = 'sqpopt' + (f'[{args["sqpopt"]}]' if args['sqpopt'] else '')
    labels = [label if s == 'sqpopt' else s for s in solvers]
    jobs = max(1, int(args['jobs']))
    done = set()
    if results.exists() and not args['redo']:
        done = {(r['problem'], r['solver']) for r in read_results(results)}
    names = select(args)
    print(f'{len(names)} problems, solvers {", ".join(labels)}, {jobs} at once, results in {results}')
    todo = {name: [s for s, lab in zip(solvers, labels) if (name, lab) not in done] for name in names}
    t0 = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        # (the threads only wait for their child processes; this one writes the results)
        futures = {pool.submit(solve_in_child, name, todo[name], args): name for name in names if todo[name]}
        try:
            for k, future in enumerate(concurrent.futures.as_completed(futures), 1):
                name = futures[future]
                records, note = future.result()
                solved = {r['solver'] for r in records}
                for s, lab in zip(solvers, labels):   # (the solvers that didn't finish: the crash or the time limit)
                    if s in todo[name] and lab not in solved:
                        records.append(dict(problem=name, solver=lab, success=False, status=-2,
                                            message=note or 'no result'))
                with results.open('a') as f:
                    for r in records:
                        f.write(json.dumps(r) + '\n')
                print(f'{k:4d}/{len(futures)} {name:10s} ' + '  '.join(
                    f'{r["solver"]}: {"ok" if r["success"] else "--"}' for r in records) +
                    (f'  ({note})' if note else ''), flush=True)
        except KeyboardInterrupt:
            pool.shutdown(wait=False, cancel_futures=True)   # (the problems being solved finish; no more start)
            raise
    print(f'{time.perf_counter() - t0:.0f} s')


def read_results(path: pathlib.Path) -> list[dict]:
    """the records of a results file (the last one of each problem and solver)"""
    records = {}
    for line in path.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            records[(r['problem'], r['solver'])] = r
    return list(records.values())


def report(args: dict):
    """the table of the results: for each solver, the problems it reported success on, how many of those
    are feasible with the best objective found, and its evaluations and time on the problems all solved"""
    records = read_results(pathlib.Path(args['results']))
    if args['no-overdetermined']:
        # (a record of a crash or the time limit has no sizes: the problem's other records do)
        over = {r['problem'] for r in records if r.get('m', 0) > r.get('n', 0)}
        records = [r for r in records if r['problem'] not in over]
    solvers = sorted({r['solver'] for r in records})
    problems = sorted({r['problem'] for r in records})
    by = {(r['problem'], r['solver']): r for r in records}
    feas_tol, rel_tol = 1.0e-5, 1.0e-4

    def good(r) -> bool:
        return r is not None and r.get('success', False) and r.get('violation', 1.0) <= feas_tol \
            and r.get('f') is not None and r['f'] == r['f']

    best = {}
    for name in problems:
        values = [by[(name, s)]['f'] for s in solvers if good(by.get((name, s)))]
        if values:
            best[name] = min(values)
    common = [name for name in problems if all(good(by.get((name, s))) for s in solvers)]
    print(f'{len(problems)} problems; on {len(common)} of them every solver reported success at a feasible point')
    print(f'{"solver":28s} {"success":>8s} {"feasible":>9s} {"best f":>7s} {"evaluations":>12s} {"time (s)":>9s}')
    for s in solvers:
        mine = [by.get((name, s)) for name in problems]
        n_success = sum(1 for r in mine if r is not None and r.get('success'))
        n_good = sum(1 for r in mine if good(r))
        n_best = sum(1 for name in problems if good(by.get((name, s))) and
                     by[(name, s)]['f'] <= best[name] + rel_tol*max(1.0, abs(best[name])))
        evaluations = sum(by[(name, s)][k] for name in common for k in ('n_f', 'n_g', 'n_c', 'n_jac'))
        seconds = sum(by[(name, s)]['time'] for name in common)
        print(f'{s:28s} {n_success:8d} {n_good:9d} {n_best:7d} {evaluations:12d} {seconds:9.2f}')
    print('success: the solver said so. feasible: and the point violates no bound or constraint by more than '
          f'{feas_tol:g}.')
    print(f'best f: and its objective is within {rel_tol:g} (relative) of the best such point of any solver.')
    print('evaluations and time: on the problems of the first line (objective, gradient, constraint, and '
          'Jacobian calls).')
    if args['failures']:
        for s in solvers:
            if args['failures'] not in ('all', s):
                continue
            print(f'\n{s}: not solved')
            for name in problems:
                r = by.get((name, s))
                if r is not None and not good(r):
                    print(f'  {name:10s} n={r.get("n", "?"):>3} m={r.get("m", "?"):>3} status {r.get("status")}'
                          f' viol {r.get("violation", float("nan")):.1e}  {r.get("message", "")}')


def main():
    command = sys.argv[1] if len(sys.argv) > 1 else ''
    args = {'max-n': '100', 'max-m': '100', 'problems': '', 'solvers': ','.join(SOLVERS), 'maxiter': '250',
            'timeout': '120', 'results': str(HOME / 'results.jsonl'), 'sqpopt': '', 'redo': False, 'failures': '',
            'no-overdetermined': False, 'jobs': str(os.cpu_count() or 1)}
    for a in sys.argv[2:]:
        key, _, value = a[2:].partition('=')
        if not a.startswith('--') or key not in args:
            sys.exit(f'unknown option {a}\n{__doc__}')
        if key in ('redo', 'no-overdetermined'):
            args[key] = True
        elif key == 'failures':
            args[key] = value or 'all'
        else:
            args[key] = value
    if command == 'setup':
        setup()
        return
    if command == 'report':
        report(args)
        return
    if command not in ('list', 'run', 'solve'):
        sys.exit(__doc__)
    missing = setup_environment()
    if missing:
        sys.exit('missing: ' + ', '.join(missing) + ' (run: pixi run cutest setup)')
    if command == 'list':
        names = select(args)
        print(f'{len(names)} problems with at most {args["max-n"]} variables and {args["max-m"]} constraints' +
              (', and no more constraints than variables' if args['no-overdetermined'] else '') + ':')
        print(' '.join(names))
    elif command == 'solve':
        options = {}
        for item in filter(None, args['sqpopt'].split(',')):
            key, _, value = item.partition('=')
            options[key] = float(value) if any(ch in value for ch in '.eE') else int(value)
        label = 'sqpopt' + (f'[{args["sqpopt"]}]' if args['sqpopt'] else '')
        for name in args['problems'].split(','):
            solve_one(name, args['solvers'].split(','), int(args['maxiter']), options, label)
    else:
        run(args)


if __name__ == '__main__':
    main()
