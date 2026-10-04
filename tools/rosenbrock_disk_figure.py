"""The figure of the first example of the guide's Performance page (web/performance.html): the iterates of
SQPOPT on the Rosenbrock function on the unit disk (test/test_rosenbrock_disk.f90), drawn over the contours
of the function, with the constraint's boundary. Two panels: the whole path, and the part near the solution.

The figure is inline SVG whose colors are the page's CSS variables (classes `rd-*` in web/css/style.css), so it
follows the light and dark themes. It is written into the page between the lines

    <!-- rosenbrock-disk-figure:begin -->  ...  <!-- rosenbrock-disk-figure:end -->

together with a table of the iterates (with the solver's KKT error at each). Run from the repository root:

    pixi run fpm test test_rosenbrock_disk -- --path=build/rosenbrock_disk_path.txt
    pixi run python tools/rosenbrock_disk_figure.py
"""

import re
import sys
from pathlib import Path

import contourpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PATH_FILE = ROOT / 'build' / 'rosenbrock_disk_path.txt'
PAGE = ROOT / 'web' / 'performance.html'
BEGIN, END = '<!-- rosenbrock-disk-figure:begin -->', '<!-- rosenbrock-disk-figure:end -->'

SIZE = 400                                  # width and height of a panel's plot area (px)
MARGIN = dict(left=52, right=12, top=12, bottom=36)


def objective(x, y):
    return 100*(y - x**2)**2 + (1 - x)**2


class Panel:
    """one square panel: the window [x0, x1] x [y0, y1] of the plane, with its contour levels"""

    def __init__(self, x0, x1, y0, y1, levels, ticks_x, ticks_y, labelled=()):
        self.x0, self.x1, self.y0, self.y1 = x0, x1, y0, y1
        self.levels, self.ticks_x, self.ticks_y, self.labelled = levels, ticks_x, ticks_y, labelled

    def px(self, x):
        return MARGIN['left'] + (np.asarray(x) - self.x0)/(self.x1 - self.x0)*SIZE

    def py(self, y):
        return MARGIN['top'] + (self.y1 - np.asarray(y))/(self.y1 - self.y0)*SIZE

    def inside(self, x, y):
        return self.x0 <= x <= self.x1 and self.y0 <= y <= self.y1

    def points(self, xs, ys):
        return ' '.join(f'{self.px(x):.1f},{self.py(y):.1f}' for x, y in zip(xs, ys))


def number(v, digits=4):
    """a number for a label or a tooltip, with a real minus sign"""
    return f'{v:.{digits}f}'.replace('-', '−')


def power(v):
    """a positive number in scientific notation, as m × 10^e (HTML)"""
    m, e = f'{v:.1e}'.split('e')
    return f'{m}&nbsp;&times;&nbsp;10<sup>{int(e)}</sup>'.replace('-', '−')


def panel_svg(p, path, title, zoom_of=None, zoom_box=None):
    """the SVG of one panel. `path` is the array of the iterates (iteration, x, y, f, c)"""
    w = MARGIN['left'] + SIZE + MARGIN['right']
    h = MARGIN['top'] + SIZE + MARGIN['bottom']
    left, top = MARGIN['left'], MARGIN['top']
    clip = f'rd-clip-{id(p)}'
    out = [f'<svg class="rd-figure" viewBox="0 0 {w} {h}" width="{w}" height="{h}" role="img" aria-label="{title}">',
           f'<title>{title}</title>',
           f'<defs><clipPath id="{clip}"><rect x="{left}" y="{top}" width="{SIZE}" height="{SIZE}"/></clipPath></defs>',
           f'<g clip-path="url(#{clip})">']

    # the infeasible region (outside the unit circle), and the constraint's boundary
    cx, cy = p.px(0.0), p.py(0.0)
    r = SIZE/(p.x1 - p.x0)
    out.append(f'<path class="rd-infeasible" fill-rule="evenodd" d="M{left},{top}h{SIZE}v{SIZE}h{-SIZE}z '
               f'M{cx - r:.1f},{cy:.1f}a{r:.1f},{r:.1f} 0 1,0 {2*r:.1f},0a{r:.1f},{r:.1f} 0 1,0 {-2*r:.1f},0z"/>')

    # the contours of the objective (and, for a labelled level, the point of it nearest to where
    # its label was asked for)
    xs = np.linspace(p.x0, p.x1, 401)
    ys = np.linspace(p.y0, p.y1, 401)
    gx, gy = np.meshgrid(xs, ys)
    generator = contourpy.contour_generator(x=gx, y=gy, z=objective(gx, gy))
    wanted = dict(p.labelled)
    label_at = {}
    for level in p.levels:
        lines = generator.lines(level)
        for line in lines:
            out.append(f'<polyline class="rd-contour" points="{p.points(line[:, 0], line[:, 1])}"/>')
        if level in wanted and lines:
            pts = np.vstack(lines)
            label_at[level] = pts[np.argmin(np.hypot(pts[:, 0] - wanted[level][0], pts[:, 1] - wanted[level][1]))]

    out.append(f'<circle class="rd-boundary" cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}"/>')

    # the window of the other panel
    if zoom_box is not None:
        out.append(f'<rect class="rd-zoom" x="{p.px(zoom_box.x0):.1f}" y="{p.py(zoom_box.y1):.1f}" '
                   f'width="{p.px(zoom_box.x1) - p.px(zoom_box.x0):.1f}" '
                   f'height="{p.py(zoom_box.y0) - p.py(zoom_box.y1):.1f}"/>')

    # the path, and a point for each iterate (with its values as a tooltip)
    out.append(f'<polyline class="rd-path" points="{p.points(path[:, 1], path[:, 2])}"/>')
    last = len(path) - 1
    for k, (it, x, y, f, c, kkt) in enumerate(path):
        if not p.inside(x, y):
            continue
        tip = f'iteration {int(it)}: x = {number(x)}, y = {number(y)}, f = {f:.4g}, x² + y² = {c:.4f}, KKT error = {kkt:.1e}'
        if k == last:
            s = 7.5   # (the solution: a diamond)
            out.append(f'<path class="rd-solution" d="M{p.px(x):.1f},{p.py(y) - s:.1f}l{s},{s}l{-s},{s}l{-s},{-s}z">'
                       f'<title>{tip} (the solution)</title></path>')
        else:
            out.append(f'<circle class="rd-point" cx="{p.px(x):.1f}" cy="{p.py(y):.1f}" r="4.5"><title>{tip}</title></circle>')
    out.append('</g>')

    # the frame, the ticks, and the axis labels
    out.append(f'<rect class="rd-frame" x="{left}" y="{top}" width="{SIZE}" height="{SIZE}"/>')
    for t in p.ticks_x:
        out.append(f'<text class="rd-tick" x="{p.px(t):.1f}" y="{top + SIZE + 16}" text-anchor="middle">{number(t, 2).rstrip("0").rstrip(".") or "0"}</text>')
    for t in p.ticks_y:
        out.append(f'<text class="rd-tick" x="{left - 7}" y="{p.py(t) + 4:.1f}" text-anchor="end">{number(t, 2).rstrip("0").rstrip(".") or "0"}</text>')
    out.append(f'<text class="rd-axis" x="{left + SIZE/2}" y="{top + SIZE + 32}" text-anchor="middle">x</text>')
    out.append(f'<text class="rd-axis" x="10" y="{top + SIZE/2 + 4}" text-anchor="middle">y</text>')

    # direct labels
    for level, (lx, ly) in label_at.items():
        out.append(f'<text class="rd-label rd-muted" x="{p.px(lx):.1f}" y="{p.py(ly) + 4:.1f}" text-anchor="middle">{level:g}</text>')
    for text, (lx, ly), anchor, cls in p.notes:
        out.append(f'<text class="rd-label {cls}" x="{p.px(lx):.1f}" y="{p.py(ly) + 4:.1f}" text-anchor="{anchor}">{text}</text>')
    out.append('</svg>')
    return '\n'.join(out)


def main():
    if not PATH_FILE.exists():
        sys.exit(f'{PATH_FILE} not found: run  pixi run fpm test test_rosenbrock_disk -- --path={PATH_FILE.relative_to(ROOT)}')
    path = np.loadtxt(PATH_FILE)
    x_end, y_end = path[-1, 1], path[-1, 2]

    whole = Panel(-1.6, 1.6, -1.45, 1.75, levels=[0.3, 1, 3, 10, 30, 100, 300],
                  ticks_x=[-1, 0, 1], ticks_y=[-1, 0, 1],
                  labelled=[(10, (0.0, -0.30)), (100, (-0.45, -0.79)), (300, (0.8, -1.09))])
    near = Panel(0.735, 0.805, 0.555, 0.625, levels=[0.047, 0.05, 0.055, 0.065, 0.08, 0.1, 0.13, 0.17],
                 ticks_x=[0.74, 0.76, 0.78, 0.80], ticks_y=[0.56, 0.58, 0.60, 0.62],
                 labelled=[(0.05, (0.790, 0.588)), (0.08, (0.790, 0.577)), (0.13, (0.790, 0.566))])
    whole.notes = [('start', (path[0, 1], path[0, 2] + 0.13), 'middle', ''),
                   ('1', (path[1, 1] + 0.09, path[1, 2]), 'start', ''),
                   ('2', (path[2, 1] + 0.02, path[2, 2] + 0.12), 'start', ''),
                   ('3', (path[3, 1] - 0.09, path[3, 2] + 0.02), 'end', ''),
                   ('solution', (x_end + 0.12, y_end - 0.04), 'start', ''),
                   ('x² + y² = 1', (-0.40, 0.56), 'middle', ''),
                   ('outside the disk', (-1.02, -1.25), 'middle', 'rd-muted')]
    near.notes = [(str(int(it)), (x - 0.0022, y + 0.0002), 'end', '') for it, x, y, f, c, kkt in path[5:9]]
    near.notes += [('9', (path[9, 1] + 0.0022, path[9, 2] + 0.0014), 'start', ''),
                   ('solution (10\u201312)', (x_end - 0.003, y_end - 0.0035), 'end', ''),
                   ('outside the disk', (0.7955, 0.6125), 'middle', 'rd-muted')]

    rows = '\n'.join(f'<tr><td class="num">{int(it)}</td><td class="num">{number(x, 6)}</td><td class="num">{number(y, 6)}</td>'
                     f'<td class="num">{f:.6f}</td><td class="num">{c:.6f}</td><td class="num">{power(kkt)}</td></tr>'
                     for it, x, y, f, c, kkt in path)
    fragment = f'''{BEGIN}
<div class="rd-panels">
<figure>
{panel_svg(whole, path, 'The iterates on the contours of the Rosenbrock function, with the unit circle', zoom_box=near)}
<figcaption>The whole path. Contours of the objective at 0.3, 1, 3, 10, 30, 100, and 300.</figcaption>
</figure>
<figure>
{panel_svg(near, path, 'The iterates near the solution')}
<figcaption>Near the solution: iterations 5 to 12. Contours at 0.047 to 0.17.</figcaption>
</figure>
</div>

<details class="rd-data">
<summary>The iterates as a table</summary>
<div class="table-wrap"><table class="bench compact">
<thead><tr><th class="num">iteration</th><th class="num">x</th><th class="num">y</th><th class="num">objective</th><th class="num">x&sup2; + y&sup2;</th><th class="num">KKT error</th></tr></thead>
<tbody>
{rows}
</tbody>
</table></div>
<p class="small">The KKT error is the one of the solver's <a href="index.html#optimality-test">optimality test</a>, for the scaled problem (the objective is scaled at the starting point so that its largest gradient element is 100, which is the first row's error): the solver stops when it is below <code>ktol</code> (10<sup>−6</sup>) and the constraint is satisfied to <code>ctol</code>.</p>
</details>
{END}'''

    page = PAGE.read_text()
    pattern = re.compile(re.escape(BEGIN) + '.*?' + re.escape(END), re.S)
    if not pattern.search(page):
        sys.exit(f'the markers {BEGIN} ... {END} are not in {PAGE}')
    PAGE.write_text(pattern.sub(lambda m: fragment, page))
    print(f'wrote the figure ({len(path)} iterates) into {PAGE.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
