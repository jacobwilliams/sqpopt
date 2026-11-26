// SQPOPT HS results page: a tab per comparison (NLPQLP, SLSQP), each with a
// summary, filters, a scatter chart, and a sortable table. The data are written
// by test/test_hs_suite.f90 (js/hs_results_data.js, with NLPQLP's published
// counts) and test/test_hs_slsqp.f90 (js/hs_slsqp_data.js), with --web-data.
(function () {
  'use strict';

  var SVGNS = 'http://www.w3.org/2000/svg';

  // (the theme toggle is main.js's, which the page also loads)

  var data = window.SQPOPT_HS_RESULTS;
  if (!data || !data.problems) {
    document.getElementById('hs-meta').textContent =
      'No results data found (js/hs_results_data.js): generate it with the test suite\'s --web-data option.';
    return;
  }
  var slsqp = window.SQPOPT_HS_SLSQP && window.SQPOPT_HS_SLSQP.problems ? window.SQPOPT_HS_SLSQP : null;

  // ---------- helpers ----------
  function el(tag, cls, text) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined) e.textContent = text;
    return e;
  }
  function svg(tag, attrs) {
    var e = document.createElementNS(SVGNS, tag);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    return e;
  }
  function clear(node) { while (node.firstChild) node.removeChild(node.firstChild); }
  function fmtInt(v) { return v === null || v === undefined ? '—' : v.toLocaleString('en-US'); }
  function fmtSci(v, d) { return v === null || v === undefined ? '—' : v.toExponential(d); }
  function fmtRatio(v) { return v === null || v === undefined ? '—' : v.toFixed(2); }
  function median(a) {
    if (!a.length) return null;
    var s = a.slice().sort(function (x, y) { return x - y; }), h = s.length >> 1;
    return s.length % 2 ? s[h] : 0.5 * (s[h - 1] + s[h]);
  }
  function sum(a, f) { return a.reduce(function (s, p) { return s + f(p); }, 0); }

  var OUTCOME_LABEL = { solved: 'Solved', local: 'Local', failed: 'Failed' };
  var probs = data.problems.slice().sort(function (a, b) { return a.id - b.id; });
  var total = probs.length;
  if (slsqp) {
    var byId = {};
    slsqp.problems.forEach(function (s) { byId[s.id] = s; });
    probs.forEach(function (p) { p.s = byId[p.id] || null; });
  }
  probs.forEach(function (p) {
    p.ratio_q = p.q_nf > 0 ? p.nf / p.q_nf : null;
    p.ratio_s = p.s && p.s.nf > 0 ? p.nf / p.s.nf : null;
    if (p.s) {
      var a = p.outcome === 'solved', b = p.s.outcome === 'solved';
      p.joint = a && b ? 'both' : a ? 'sqpopt' : b ? 'slsqp' : 'neither';
    }
  });

  // ---------- the comparisons ----------
  // the column descriptions shared by both tables
  function sqpoptCols() {
    return {
      id: { label: 'TP', type: 'int', get: function (p) { return p.id; },
            desc: 'Problem number in Schittkowski’s collection (TP1–TP395).' },
      n: { label: 'n', type: 'int', get: function (p) { return p.n; }, desc: 'Number of variables.' },
      m: { label: 'm', type: 'int', get: function (p) { return p.m; },
           desc: 'Number of general constraints, equality and inequality (not counting variable bounds).' },
      me: { label: 'me', type: 'int', get: function (p) { return p.me; }, desc: 'Number of equality constraints.' },
      iter: { label: 'iter', type: 'int', get: function (p) { return p.iter; }, desc: 'SQPOPT’s major (SQP) iterations.' },
      nf: { label: 'fc', type: 'int', get: function (p) { return p.nf; },
            desc: 'SQPOPT’s calls of the objective and constraint function fc, including the line search’s trial points.' },
      ng: { label: 'gjac', type: 'int', get: function (p) { return p.ng; },
            desc: 'SQPOPT’s calls of the gradient and Jacobian function gjac.' },
      f: { label: 'f', type: 'sci3', get: function (p) { return p.f; }, desc: 'The objective at SQPOPT’s final point.' },
      f_star: { label: 'f*', type: 'sci3', get: function (p) { return p.f_star; },
                desc: 'The validated optimal objective, from the collection.' },
      rel: { label: 'rel. err', type: 'sci1', get: function (p) { return p.rel; },
             desc: 'Relative error of SQPOPT’s objective, (f − f*) / max(1, |f*|). Negative: better than the published optimum.' },
      viol: { label: 'viol', type: 'sci1', get: function (p) { return p.viol; },
              desc: 'The largest violation of a constraint or variable bound at SQPOPT’s final point (of the original, unscaled problem).' },
      status: { label: 'status', type: 'text', get: function (p) { return p.status; }, desc: 'SQPOPT’s termination status.' },
      fd: { label: 'notes', type: 'fd', get: function (p) { return p.fd; },
            desc: 'fd: a problem without analytic derivatives in the collection, solved with central-difference derivatives (by both solvers).' }
    };
  }
  var C = sqpoptCols();
  var resultDesc = 'Solved: the final point is feasible (violation ≤ ' + data.feas_tol.toExponential(0) +
    ') and its objective is within ' + data.rel_tol.toExponential(0) + ' (relative) of the validated optimum, ' +
    'or better. Local: feasible and converged, but at a worse objective (a different local solution). Failed: anything else.';

  var MODES = {
    nlpqlp: {
      tab: 'vs NLPQLP',
      meta: 'SQPOPT options: ' + data.options + ' · generated ' + data.generated + ' · ' + data.compiler +
            ' · NLPQLP: the counts published with the collection',
      chartTitle: 'Function evaluations: SQPOPT vs NLPQLP',
      chartSub: 'One dot per problem: SQPOPT’s calls of fc against NLPQLP’s function evaluations (both log scales). ' +
                'Below the diagonal, SQPOPT needed fewer. NLPQLP’s counts exclude the extra evaluations of its ' +
                'finite-difference gradients, so they understate its cost. Click a dot to find the problem in the table.',
      tableSub: 'Hover over (or tab to) a column heading for its description, and click it to sort. fc and gjac are ' +
                'SQPOPT’s function calls; Q:nf and Q:ng are NLPQLP’s.',
      xLabel: 'NLPQLP function evaluations',
      x: function (p) { return p.q_nf; },
      cat: function (p) { return p.outcome; },
      cats: [['solved', 'Solved'], ['local', 'Local'], ['failed', 'Failed']],
      tip: function (p) { return [[fmtInt(p.nf), 'SQPOPT fc calls'], [fmtInt(p.q_nf), 'NLPQLP function evaluations']]; },
      tiles: function () {
        var solved = probs.filter(function (p) { return p.outcome === 'solved'; });
        var med = median(solved.map(function (p) { return p.ratio_q; }).filter(function (r) { return r !== null; }));
        function n(o) { return probs.filter(function (p) { return p.outcome === o; }).length; }
        return [
          [n('solved') + ' / ' + total, 'solved'],
          [String(n('local')), 'local solutions'],
          [String(n('failed')), 'failed'],
          [fmtInt(sum(solved, function (p) { return p.nf; })),
           'fc calls on the solved problems (NLPQLP: ' + fmtInt(sum(solved, function (p) { return p.q_nf; })) + ')'],
          [solved.filter(function (p) { return p.nf < p.q_nf; }).length + ' / ' + solved.length,
           'solved with fewer evaluations than NLPQLP'],
          [med === null ? '—' : med.toFixed(2) + '×', 'median ratio of SQPOPT’s to NLPQLP’s evaluations']
        ];
      },
      cols: [C.id, C.n, C.m, C.me,
        { label: 'result', type: 'outcome', get: function (p) { return p.outcome; }, desc: resultDesc },
        C.iter, C.nf, C.ng,
        { label: 'Q:nf', type: 'int', get: function (p) { return p.q_nf; },
          desc: 'NLPQLP’s function evaluations, as published with the collection. They exclude the evaluations ' +
                'for its finite-difference gradients (see Q:ng).' },
        { label: 'Q:ng', type: 'int', get: function (p) { return p.q_ng; },
          desc: 'NLPQLP’s gradient evaluations. They were finite-difference approximations, each costing about ' +
                'n more function evaluations (not included in Q:nf).' },
        { label: 'fc / Q:nf', type: 'ratio', get: function (p) { return p.ratio_q; },
          desc: 'SQPOPT’s fc calls divided by NLPQLP’s function evaluations: below 1, SQPOPT needed fewer.' },
        C.f, C.f_star, C.rel, C.viol, C.status, C.fd]
    }
  };
  if (slsqp) {
    MODES.slsqp = {
      tab: 'vs SLSQP',
      meta: 'SQPOPT options: ' + data.options + ' · SLSQP options: ' + slsqp.options + ' · generated ' +
            data.generated + ' and ' + slsqp.generated + ' · ' + data.compiler,
      chartTitle: 'Function evaluations: SQPOPT vs SLSQP',
      chartSub: 'One dot per problem: SQPOPT’s calls of fc against SLSQP’s calls of its function (the objective ' +
                'and constraints together, like fc), both log scales. Below the diagonal, SQPOPT needed fewer. Both ' +
                'solve each problem from the same starting point with the same derivatives. Click a dot to find the ' +
                'problem in the table.',
      tableSub: 'Hover over (or tab to) a column heading for its description, and click it to sort. The S: columns ' +
                'are SLSQP’s.',
      xLabel: 'SLSQP function calls',
      x: function (p) { return p.s ? p.s.nf : null; },
      cat: function (p) { return p.joint; },
      cats: [['both', 'Both solved'], ['sqpopt', 'Only SQPOPT solved'], ['slsqp', 'Only SLSQP solved'],
             ['neither', 'Neither solved']],
      tip: function (p) { return [[fmtInt(p.nf), 'SQPOPT fc calls'], [fmtInt(p.s.nf), 'SLSQP function calls']]; },
      tiles: function () {
        function n(c) { return probs.filter(function (p) { return p.joint === c; }).length; }
        var both = probs.filter(function (p) { return p.joint === 'both'; });
        var med = median(both.map(function (p) { return p.ratio_s; }).filter(function (r) { return r !== null; }));
        return [
          [probs.filter(function (p) { return p.outcome === 'solved'; }).length + ' / ' + total, 'solved by SQPOPT'],
          [probs.filter(function (p) { return p.s && p.s.outcome === 'solved'; }).length + ' / ' + total, 'solved by SLSQP'],
          [n('sqpopt') + ' vs ' + n('slsqp'), 'solved by only SQPOPT vs only SLSQP (' + n('neither') + ' by neither)'],
          [fmtInt(sum(both, function (p) { return p.nf; })),
           'SQPOPT fc calls on the ' + both.length + ' problems both solve (SLSQP: ' +
           fmtInt(sum(both, function (p) { return p.s.nf; })) + ')'],
          [both.filter(function (p) { return p.nf < p.s.nf; }).length + ' / ' + both.length,
           'of those, solved with fewer function calls than SLSQP'],
          [med === null ? '—' : med.toFixed(2) + '×', 'median ratio of SQPOPT’s to SLSQP’s function calls']
        ];
      },
      cols: [C.id, C.n, C.m, C.me,
        { label: 'SQPOPT', type: 'outcome', get: function (p) { return p.outcome; }, desc: 'SQPOPT’s result. ' + resultDesc },
        { label: 'SLSQP', type: 'outcome', get: function (p) { return p.s.outcome; },
          desc: 'SLSQP’s result, by the same criteria: solved, local (feasible and converged, at a worse objective), or failed.' },
        C.iter,
        { label: 'S:iter', type: 'int', get: function (p) { return p.s.iter; }, desc: 'SLSQP’s iterations.' },
        C.nf,
        { label: 'S:nf', type: 'int', get: function (p) { return p.s.nf; },
          desc: 'SLSQP’s calls of its function routine (the objective and the constraints together, like fc).' },
        C.ng,
        { label: 'S:ng', type: 'int', get: function (p) { return p.s.ng; },
          desc: 'SLSQP’s calls of its gradient routine (the objective gradient and the constraint Jacobian, like gjac).' },
        { label: 'fc / S:nf', type: 'ratio', get: function (p) { return p.ratio_s; },
          desc: 'SQPOPT’s fc calls divided by SLSQP’s function calls: below 1, SQPOPT needed fewer.' },
        C.f,
        { label: 'S:f', type: 'sci3', get: function (p) { return p.s.f; }, desc: 'The objective at SLSQP’s final point.' },
        C.f_star, C.viol,
        { label: 'S:viol', type: 'sci1', get: function (p) { return p.s.viol; },
          desc: 'The largest violation of a constraint or variable bound at SLSQP’s final point.' },
        C.status,
        { label: 'S:status', type: 'text', get: function (p) { return p.s.status; }, desc: 'SLSQP’s termination status.' },
        C.fd]
    };
  } else {
    document.getElementById('hs-tab-slsqp').hidden = true;
  }

  var state = { mode: 'nlpqlp', cat: 'all', terms: [], hideFd: false, sortCol: 0, sortDir: 1 };
  function M() { return MODES[state.mode]; }

  // ---------- filters ----------
  function matches(p) {
    if (state.cat !== 'all' && M().cat(p) !== state.cat) return false;
    if (state.hideFd && p.fd) return false;
    var text = (p.status + ' ' + p.outcome + ' ' + (p.s ? p.s.status + ' ' + p.s.outcome : '') + (p.fd ? ' fd' : '')).toLowerCase();
    return state.terms.every(function (t) {
      if (/^\d+$/.test(t)) return String(p.id) === t;
      return text.indexOf(t) !== -1;
    });
  }

  var chipBox = document.getElementById('hs-outcome');
  var chips = [];
  function buildChips() {
    clear(chipBox);
    chips = [];
    [['all', 'All']].concat(M().cats).forEach(function (c) {
      var n = c[0] === 'all' ? total : probs.filter(function (p) { return M().cat(p) === c[0]; }).length;
      if (c[0] !== 'all' && n === 0) return;
      var b = el('button', 'hs-chip');
      b.type = 'button';
      b.setAttribute('role', 'radio');
      b.dataset.cat = c[0];
      if (c[0] !== 'all') b.appendChild(el('span', 'hs-key ' + c[0]));
      b.appendChild(document.createTextNode(c[1] + ' '));
      b.appendChild(el('span', 'hs-chip-n', String(n)));
      b.addEventListener('click', function () { state.cat = c[0]; update(); });
      chipBox.appendChild(b);
      chips.push(b);
    });
  }
  document.getElementById('hs-search').addEventListener('input', function (e) {
    state.terms = e.target.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
    update();
  });
  document.getElementById('hs-hide-fd').addEventListener('change', function (e) {
    state.hideFd = e.target.checked;
    update();
  });

  // ---------- summary, legend, texts ----------
  function buildTiles() {
    var box = document.getElementById('hs-tiles');
    clear(box);
    M().tiles().forEach(function (it) {
      var t = el('div', 'hs-tile');
      t.appendChild(el('div', 'hs-tile-value', it[0]));
      t.appendChild(el('div', 'hs-tile-label', it[1]));
      box.appendChild(t);
    });
  }
  function buildLegend() {
    var box = document.getElementById('hs-legend');
    clear(box);
    M().cats.forEach(function (c) {
      if (!probs.some(function (p) { return M().cat(p) === c[0]; })) return;
      var item = el('span', 'hs-legend-item');
      item.appendChild(el('span', 'hs-key ' + c[0]));
      item.appendChild(document.createTextNode(c[1]));
      box.appendChild(item);
    });
    var d = el('span', 'hs-legend-item');
    d.appendChild(el('span', 'hs-key-line'));
    d.appendChild(document.createTextNode('equal evaluations'));
    box.appendChild(d);
  }
  function setTexts() {
    document.getElementById('hs-meta').textContent = M().meta;
    document.getElementById('hs-chart-title').textContent = M().chartTitle;
    document.getElementById('hs-chart-sub').textContent = M().chartSub;
    document.getElementById('hs-table-sub').textContent = M().tableSub;
    document.getElementById('hs-svg').setAttribute('aria-label', 'Scatter plot of ' + M().chartTitle.replace('Function evaluations: ', '') +
      ' function evaluations, one dot per problem. The same data is in the table below.');
  }

  // ---------- chart ----------
  var chartBox = document.getElementById('hs-chart');
  var svgEl = document.getElementById('hs-svg');
  var tip = document.getElementById('hs-tooltip');
  var plotted = [];   // [{p, x, y}] in pixel coordinates, for the hover search
  var hoverRing = null;

  function domain() {
    // a common log domain for both axes, from all the problems (so filtering doesn't move the axes)
    var v = [];
    probs.forEach(function (p) { var x = M().x(p); if (p.nf > 0) v.push(p.nf); if (x > 0) v.push(x); });
    return [Math.floor(Math.log10(Math.min.apply(null, v))), Math.ceil(Math.log10(Math.max.apply(null, v)))];
  }
  function tickLabel(e) { var v = Math.pow(10, e); return v >= 1000 ? (v / 1000) + 'k' : String(v); }

  function drawChart(visible) {
    var dom = domain(), lo = dom[0], hi = dom[1];
    var w = chartBox.clientWidth;
    var h = Math.max(300, Math.min(460, Math.round(w * 0.62)));
    var m = { l: 56, r: 16, t: 14, b: 46 };
    var pw = w - m.l - m.r, ph = h - m.t - m.b;
    svgEl.setAttribute('viewBox', '0 0 ' + w + ' ' + h);
    svgEl.setAttribute('width', w);
    svgEl.setAttribute('height', h);
    clear(svgEl);
    function sx(v) { return m.l + (Math.log10(v) - lo) / (hi - lo) * pw; }
    function sy(v) { return m.t + ph - (Math.log10(v) - lo) / (hi - lo) * ph; }

    var g = svg('g', { 'class': 'hs-axis' });
    for (var e = lo; e <= hi; e++) {
      var v = Math.pow(10, e);
      g.appendChild(svg('line', { x1: sx(v), x2: sx(v), y1: m.t, y2: m.t + ph, 'class': 'hs-grid' }));
      g.appendChild(svg('line', { x1: m.l, x2: m.l + pw, y1: sy(v), y2: sy(v), 'class': 'hs-grid' }));
      var tx = svg('text', { x: sx(v), y: m.t + ph + 18, 'text-anchor': 'middle' }); tx.textContent = tickLabel(e); g.appendChild(tx);
      var ty = svg('text', { x: m.l - 8, y: sy(v) + 4, 'text-anchor': 'end' }); ty.textContent = tickLabel(e); g.appendChild(ty);
    }
    var xl = svg('text', { x: m.l + pw / 2, y: h - 6, 'text-anchor': 'middle', 'class': 'hs-axis-title' });
    xl.textContent = M().xLabel;
    g.appendChild(xl);
    var yl = svg('text', { x: 14, y: m.t + ph / 2, 'text-anchor': 'middle', 'class': 'hs-axis-title',
                           transform: 'rotate(-90 14 ' + (m.t + ph / 2) + ')' });
    yl.textContent = 'SQPOPT fc calls';
    g.appendChild(yl);
    svgEl.appendChild(g);

    var a = Math.pow(10, lo), b = Math.pow(10, hi);
    svgEl.appendChild(svg('line', { x1: sx(a), y1: sy(a), x2: sx(b), y2: sy(b), 'class': 'hs-diag' }));

    // the dots, in the categories' order (so the rarer ones are drawn on top),
    // except that the neutral "neither" is drawn first, underneath
    plotted = [];
    var dots = svg('g', {});
    var order = M().cats.map(function (c) { return c[0]; });
    if (order.indexOf('neither') !== -1) order = ['neither'].concat(order.filter(function (c) { return c !== 'neither'; }));
    order.forEach(function (c) {
      visible.forEach(function (p) {
        var xv = M().x(p);
        if (M().cat(p) !== c || !(p.nf > 0 && xv > 0)) return;
        var x = sx(xv), y = sy(p.nf);
        dots.appendChild(svg('circle', { cx: x, cy: y, r: 4, 'class': 'hs-dot ' + c }));
        plotted.push({ p: p, x: x, y: y });
      });
    });
    svgEl.appendChild(dots);
    hoverRing = svg('circle', { r: 7, 'class': 'hs-hover-ring', visibility: 'hidden' });
    svgEl.appendChild(hoverRing);
  }

  function nearest(evt) {
    var r = svgEl.getBoundingClientRect();
    var px = evt.clientX - r.left, py = evt.clientY - r.top, best = null, bd = 24 * 24;
    plotted.forEach(function (q) {
      var d = (q.x - px) * (q.x - px) + (q.y - py) * (q.y - py);
      if (d < bd) { bd = d; best = q; }
    });
    return best;
  }

  function showTip(q) {
    if (!q) { tip.hidden = true; if (hoverRing) hoverRing.setAttribute('visibility', 'hidden'); return; }
    clear(tip);
    tip.appendChild(el('div', 'hs-tip-title', 'TP' + q.p.id + '  (n=' + q.p.n + ', m=' + q.p.m + ')'));
    M().tip(q.p).forEach(function (row) {
      var d = el('div', 'hs-tip-row');
      d.appendChild(el('strong', null, row[0]));
      d.appendChild(el('span', null, ' ' + row[1]));
      tip.appendChild(d);
    });
    var c = M().cat(q.p), label = M().cats.filter(function (x) { return x[0] === c; })[0][1];
    var o = el('div', 'hs-tip-row');
    o.appendChild(el('span', 'hs-key-line ' + c));
    o.appendChild(el('span', null, label + (q.p.fd ? ' (fd)' : '')));
    tip.appendChild(o);
    tip.hidden = false;
    var bw = chartBox.clientWidth, tw = tip.offsetWidth;
    tip.style.left = Math.min(Math.max(q.x + 14, 0), bw - tw - 4) + 'px';
    tip.style.top = Math.max(q.y - tip.offsetHeight - 10, 0) + 'px';
    hoverRing.setAttribute('cx', q.x);
    hoverRing.setAttribute('cy', q.y);
    hoverRing.setAttribute('visibility', 'visible');
  }

  svgEl.addEventListener('pointermove', function (e) { var q = nearest(e); showTip(q); svgEl.style.cursor = q ? 'pointer' : ''; });
  svgEl.addEventListener('pointerleave', function () { showTip(null); });
  svgEl.addEventListener('click', function (e) {
    var q = nearest(e);
    if (!q) return;
    var row = rows[q.p.id];
    row.scrollIntoView({ behavior: 'smooth', block: 'center' });
    row.classList.remove('hs-flash');
    void row.offsetWidth;   // restart the highlight animation
    row.classList.add('hs-flash');
  });

  // ---------- table ----------
  // the column descriptions: a tooltip on hover and keyboard focus of each heading
  // (fixed-positioned, so the table's scrolling box can't clip it), and the
  // heading button's accessible description
  var headTip = el('div', 'hs-tooltip hs-head-tip');
  headTip.hidden = true;
  headTip.setAttribute('role', 'tooltip');
  document.body.appendChild(headTip);
  function showHeadTip(btn, text) {
    headTip.textContent = text;
    headTip.hidden = false;
    var r = btn.getBoundingClientRect(), tw = headTip.offsetWidth;
    headTip.style.left = Math.max(8, Math.min(r.left, window.innerWidth - tw - 8)) + 'px';
    headTip.style.top = (r.bottom + 8) + 'px';
  }
  function hideHeadTip() { headTip.hidden = true; }
  window.addEventListener('scroll', hideHeadTip, { passive: true });

  var head = document.getElementById('hs-head');
  var body = document.getElementById('hs-body');
  var ths = [], rows = {};

  function buildTable() {
    var cols = M().cols;
    clear(head);
    clear(body);
    ths = cols.map(function (c, i) {
      var th = el('th');
      th.scope = 'col';
      if (c.type !== 'text' && c.type !== 'outcome' && c.type !== 'fd') th.className = 'num';
      var btn = el('button', 'hs-sort', c.label);
      btn.type = 'button';
      var desc = el('span', 'hs-visually-hidden', c.desc);
      desc.id = 'hs-desc-' + state.mode + '-' + i;
      th.appendChild(desc);
      btn.setAttribute('aria-describedby', desc.id);
      btn.addEventListener('pointerenter', function () { showHeadTip(btn, c.desc); });
      btn.addEventListener('pointerleave', hideHeadTip);
      btn.addEventListener('focus', function () { showHeadTip(btn, c.desc); });
      btn.addEventListener('blur', hideHeadTip);
      btn.addEventListener('keydown', function (e) { if (e.key === 'Escape') hideHeadTip(); });
      btn.addEventListener('click', function () {
        if (state.sortCol === i) state.sortDir = -state.sortDir;
        else { state.sortCol = i; state.sortDir = 1; }
        update();
      });
      th.appendChild(btn);
      head.appendChild(th);
      return th;
    });
    rows = {};
    probs.forEach(function (p) {
      var tr = el('tr');
      tr.id = 'tp-' + p.id;
      cols.forEach(function (c) {
        var td = el('td');
        var v = c.get(p);
        switch (c.type) {
          case 'int': td.className = 'num'; td.textContent = fmtInt(v); break;
          case 'ratio': td.className = 'num'; td.textContent = fmtRatio(v); break;
          case 'sci3': td.className = 'num'; td.textContent = fmtSci(v, 3); break;
          case 'sci1': td.className = 'num'; td.textContent = fmtSci(v, 1); break;
          case 'outcome':
            // (colored keys only where the colors are the outcomes', as in the
            // chart; in the SLSQP comparison the colors are the joint categories)
            if (state.mode === 'nlpqlp') td.appendChild(el('span', 'hs-key ' + v));
            td.appendChild(document.createTextNode(OUTCOME_LABEL[v]));
            break;
          case 'fd': td.textContent = v ? 'fd' : ''; break;
          default: td.textContent = v;
        }
        tr.appendChild(td);
      });
      rows[p.id] = tr;
      body.appendChild(tr);
    });
  }

  function compare(a, b) {
    var c = M().cols[state.sortCol], x = c.get(a), y = c.get(b);
    if (c.type === 'outcome') { x = ['solved', 'local', 'failed'].indexOf(x); y = ['solved', 'local', 'failed'].indexOf(y); }
    if (x === null || x === undefined) return 1;    // (missing values last, either way)
    if (y === null || y === undefined) return -1;
    var r = typeof x === 'string' ? x.localeCompare(y) : (x < y ? -1 : x > y ? 1 : 0);
    return (r || a.id - b.id) * state.sortDir;
  }

  // ---------- tabs ----------
  var tabs = Array.prototype.slice.call(document.querySelectorAll('.hs-tab'));
  function setMode(mode, focus, initial) {
    if (!MODES[mode]) mode = 'nlpqlp';
    state.mode = mode;
    state.cat = 'all';
    state.sortCol = 0;
    state.sortDir = 1;
    tabs.forEach(function (t) {
      var on = t.dataset.mode === mode;
      t.setAttribute('aria-selected', String(on));
      t.tabIndex = on ? 0 : -1;
      if (on && focus) t.focus();
    });
    // (the page has other sections: keep the reader's anchor on load, and name this section's on a switch)
    if (!initial) {
      try { history.replaceState(null, '', mode === 'nlpqlp' ? '#hs-results' : '#' + mode); } catch (e) {}
    }
    setTexts();
    buildTiles();
    buildChips();
    buildLegend();
    buildTable();
    update();
  }
  tabs.forEach(function (t, i) {
    t.addEventListener('click', function () { setMode(t.dataset.mode); });
    t.addEventListener('keydown', function (e) {   // arrow keys move between the tabs
      var visibleTabs = tabs.filter(function (x) { return !x.hidden; }), k = visibleTabs.indexOf(t);
      if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') {
        var next = visibleTabs[(k + (e.key === 'ArrowRight' ? 1 : visibleTabs.length - 1)) % visibleTabs.length];
        setMode(next.dataset.mode, true);
        e.preventDefault();
      }
    });
  });

  // ---------- update everything below the filters ----------
  function update() {
    chips.forEach(function (b) { b.setAttribute('aria-checked', String(b.dataset.cat === state.cat)); });
    var visible = probs.filter(matches);
    var sorted = visible.slice().sort(compare);
    probs.forEach(function (p) { rows[p.id].hidden = true; });
    sorted.forEach(function (p) { rows[p.id].hidden = false; body.appendChild(rows[p.id]); });
    ths.forEach(function (th, i) {
      th.setAttribute('aria-sort', i === state.sortCol ? (state.sortDir > 0 ? 'ascending' : 'descending') : 'none');
    });
    document.getElementById('hs-count').textContent = 'Showing ' + visible.length + ' of ' + total;
    drawChart(visible);
  }

  if (window.ResizeObserver) {
    var lastW = 0;
    new ResizeObserver(function () {
      if (chartBox.clientWidth !== lastW) { lastW = chartBox.clientWidth; drawChart(probs.filter(matches)); }
    }).observe(chartBox);
  } else {
    window.addEventListener('resize', function () { drawChart(probs.filter(matches)); });
  }
  setMode(location.hash === '#slsqp' ? 'slsqp' : 'nlpqlp', false, true);
})();
