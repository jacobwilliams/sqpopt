// SQPOPT HS results page: summary, filters, scatter chart, and sortable table,
// from the data written by test/test_hs_suite.f90 (--web-data) to js/hs_results_data.js.
(function () {
  'use strict';

  var root = document.documentElement;
  var SVGNS = 'http://www.w3.org/2000/svg';

  // ---------- theme toggle (as in the guide) ----------
  function currentTheme() {
    var t = root.getAttribute('data-theme');
    if (t) return t;
    return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  }
  document.querySelector('.theme-toggle').addEventListener('click', function () {
    var next = currentTheme() === 'dark' ? 'light' : 'dark';
    root.setAttribute('data-theme', next);
    try { localStorage.setItem('sqpopt-theme', next); } catch (e) {}
  });

  var data = window.SQPOPT_HS_RESULTS;
  if (!data || !data.problems) {
    document.getElementById('hs-meta').textContent =
      'No results data found (js/hs_results_data.js): generate it with the test suite\'s --web-data option.';
    return;
  }

  var OUTCOMES = ['solved', 'local', 'failed'];
  var LABEL = { solved: 'Solved', local: 'Local', failed: 'Failed' };
  var probs = data.problems.slice().sort(function (a, b) { return a.id - b.id; });
  probs.forEach(function (p) { p.ratio = p.q_nf > 0 ? p.nf / p.q_nf : null; });
  var total = probs.length;

  var state = { outcome: 'all', terms: [], hideFd: false, sortKey: 'id', sortDir: 1 };

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
  function fmtInt(v) { return v === null || v === undefined ? '—' : v.toLocaleString('en-US'); }
  function fmtSci(v, d) { return v === null || v === undefined ? '—' : v.toExponential(d); }
  function fmtRatio(v) { return v === null ? '—' : v.toFixed(2); }
  function median(a) {
    if (!a.length) return null;
    var s = a.slice().sort(function (x, y) { return x - y; }), h = s.length >> 1;
    return s.length % 2 ? s[h] : 0.5 * (s[h - 1] + s[h]);
  }
  function count(outcome) { return probs.filter(function (p) { return p.outcome === outcome; }).length; }

  // ---------- meta and summary ----------
  document.getElementById('hs-meta').textContent =
    'Options: ' + data.options + ' · generated ' + data.generated + ' · ' + data.compiler;

  (function tiles() {
    var solved = probs.filter(function (p) { return p.outcome === 'solved'; });
    var nf = solved.reduce(function (s, p) { return s + p.nf; }, 0);
    var qnf = solved.reduce(function (s, p) { return s + p.q_nf; }, 0);
    var fewer = solved.filter(function (p) { return p.nf < p.q_nf; }).length;
    var med = median(solved.map(function (p) { return p.ratio; }).filter(function (r) { return r !== null; }));
    var items = [
      [count('solved') + ' / ' + total, 'solved'],
      [String(count('local')), 'local solutions'],
      [String(count('failed')), 'failed'],
      [fmtInt(nf), 'fc calls on the solved problems (NLPQLP: ' + fmtInt(qnf) + ')'],
      [fewer + ' / ' + solved.length, 'solved with fewer evaluations than NLPQLP'],
      [med === null ? '—' : med.toFixed(2) + '×', 'median ratio of SQPOPT’s to NLPQLP’s evaluations']
    ];
    var box = document.getElementById('hs-tiles');
    items.forEach(function (it) {
      var t = el('div', 'hs-tile');
      t.appendChild(el('div', 'hs-tile-value', it[0]));
      t.appendChild(el('div', 'hs-tile-label', it[1]));
      box.appendChild(t);
    });
  })();

  // ---------- filters ----------
  function matches(p) {
    if (state.outcome !== 'all' && p.outcome !== state.outcome) return false;
    if (state.hideFd && p.fd) return false;
    return state.terms.every(function (t) {
      if (/^\d+$/.test(t)) return String(p.id) === t;
      return (p.status + ' ' + p.outcome + (p.fd ? ' fd' : '')).toLowerCase().indexOf(t) !== -1;
    });
  }

  var chipBox = document.getElementById('hs-outcome');
  var chips = [];
  ['all'].concat(OUTCOMES).forEach(function (o) {
    var n = o === 'all' ? total : count(o);
    if (o !== 'all' && n === 0) return;
    var b = el('button', 'hs-chip');
    b.type = 'button';
    b.setAttribute('role', 'radio');
    b.dataset.outcome = o;
    if (o !== 'all') b.appendChild(el('span', 'hs-key ' + o));
    b.appendChild(document.createTextNode((o === 'all' ? 'All' : LABEL[o]) + ' '));
    b.appendChild(el('span', 'hs-chip-n', String(n)));
    b.addEventListener('click', function () { state.outcome = o; update(); });
    chipBox.appendChild(b);
    chips.push(b);
  });
  document.getElementById('hs-search').addEventListener('input', function (e) {
    state.terms = e.target.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
    update();
  });
  document.getElementById('hs-hide-fd').addEventListener('change', function (e) {
    state.hideFd = e.target.checked;
    update();
  });

  // ---------- chart ----------
  var chartBox = document.getElementById('hs-chart');
  var svgEl = document.getElementById('hs-svg');
  var tip = document.getElementById('hs-tooltip');
  var plotted = [];   // [{p, x, y}] in pixel coordinates, for the hover search
  var hoverRing = null;

  // a common log domain for both axes, from all the problems (so filtering doesn't move the axes)
  var allCounts = [];
  probs.forEach(function (p) { if (p.nf > 0) allCounts.push(p.nf); if (p.q_nf > 0) allCounts.push(p.q_nf); });
  var lo = Math.floor(Math.log10(Math.min.apply(null, allCounts)));
  var hi = Math.ceil(Math.log10(Math.max.apply(null, allCounts)));

  (function legend() {
    var box = document.getElementById('hs-legend');
    OUTCOMES.forEach(function (o) {
      if (!count(o)) return;
      var item = el('span', 'hs-legend-item');
      item.appendChild(el('span', 'hs-key ' + o));
      item.appendChild(document.createTextNode(LABEL[o]));
      box.appendChild(item);
    });
    var d = el('span', 'hs-legend-item');
    d.appendChild(el('span', 'hs-key-line'));
    d.appendChild(document.createTextNode('equal evaluations'));
    box.appendChild(d);
  })();

  function tickLabel(e) { var v = Math.pow(10, e); return v >= 1000 ? (v / 1000) + 'k' : String(v); }

  function drawChart(visible) {
    var w = chartBox.clientWidth;
    var h = Math.max(300, Math.min(460, Math.round(w * 0.62)));
    var m = { l: 56, r: 16, t: 14, b: 46 };
    var pw = w - m.l - m.r, ph = h - m.t - m.b;
    svgEl.setAttribute('viewBox', '0 0 ' + w + ' ' + h);
    svgEl.setAttribute('width', w);
    svgEl.setAttribute('height', h);
    while (svgEl.firstChild) svgEl.removeChild(svgEl.firstChild);
    function sx(v) { return m.l + (Math.log10(v) - lo) / (hi - lo) * pw; }
    function sy(v) { return m.t + ph - (Math.log10(v) - lo) / (hi - lo) * ph; }

    // grid and axes (decades)
    var g = svg('g', { 'class': 'hs-axis' });
    for (var e = lo; e <= hi; e++) {
      var v = Math.pow(10, e);
      g.appendChild(svg('line', { x1: sx(v), x2: sx(v), y1: m.t, y2: m.t + ph, 'class': 'hs-grid' }));
      g.appendChild(svg('line', { x1: m.l, x2: m.l + pw, y1: sy(v), y2: sy(v), 'class': 'hs-grid' }));
      var tx = svg('text', { x: sx(v), y: m.t + ph + 18, 'text-anchor': 'middle' }); tx.textContent = tickLabel(e); g.appendChild(tx);
      var ty = svg('text', { x: m.l - 8, y: sy(v) + 4, 'text-anchor': 'end' }); ty.textContent = tickLabel(e); g.appendChild(ty);
    }
    var xl = svg('text', { x: m.l + pw / 2, y: h - 6, 'text-anchor': 'middle', 'class': 'hs-axis-title' });
    xl.textContent = 'NLPQLP function evaluations';
    g.appendChild(xl);
    var yl = svg('text', { x: 14, y: m.t + ph / 2, 'text-anchor': 'middle', 'class': 'hs-axis-title',
                           transform: 'rotate(-90 14 ' + (m.t + ph / 2) + ')' });
    yl.textContent = 'SQPOPT fc calls';
    g.appendChild(yl);
    svgEl.appendChild(g);

    // the diagonal (equal evaluations)
    var a = Math.pow(10, lo), b = Math.pow(10, hi);
    svgEl.appendChild(svg('line', { x1: sx(a), y1: sy(a), x2: sx(b), y2: sy(b), 'class': 'hs-diag' }));

    // the dots: solved first, so the rarer outcomes are drawn on top
    plotted = [];
    var dots = svg('g', {});
    OUTCOMES.forEach(function (o) {
      visible.forEach(function (p) {
        if (p.outcome !== o || !(p.nf > 0 && p.q_nf > 0)) return;
        var x = sx(p.q_nf), y = sy(p.nf);
        dots.appendChild(svg('circle', { cx: x, cy: y, r: 4, 'class': 'hs-dot ' + o }));
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
    if (!q) { tip.hidden = true; hoverRing.setAttribute('visibility', 'hidden'); return; }
    while (tip.firstChild) tip.removeChild(tip.firstChild);
    tip.appendChild(el('div', 'hs-tip-title', 'TP' + q.p.id + '  (n=' + q.p.n + ', m=' + q.p.m + ')'));
    [[fmtInt(q.p.nf), 'SQPOPT fc calls'], [fmtInt(q.p.q_nf), 'NLPQLP function evaluations']].forEach(function (row) {
      var d = el('div', 'hs-tip-row');
      d.appendChild(el('strong', null, row[0]));
      d.appendChild(el('span', null, ' ' + row[1]));
      tip.appendChild(d);
    });
    var o = el('div', 'hs-tip-row');
    o.appendChild(el('span', 'hs-key-line ' + q.p.outcome));
    o.appendChild(el('span', null, LABEL[q.p.outcome] + (q.p.fd ? ' (fd)' : '')));
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
  var COLS = [
    { key: 'id', label: 'TP', type: 'int' },
    { key: 'n', label: 'n', type: 'int' },
    { key: 'm', label: 'm', type: 'int' },
    { key: 'me', label: 'me', type: 'int' },
    { key: 'outcome', label: 'result', type: 'outcome' },
    { key: 'iter', label: 'iter', type: 'int' },
    { key: 'nf', label: 'fc', type: 'int' },
    { key: 'ng', label: 'gjac', type: 'int' },
    { key: 'q_nf', label: 'Q:nf', type: 'int' },
    { key: 'q_ng', label: 'Q:ng', type: 'int' },
    { key: 'ratio', label: 'fc / Q:nf', type: 'ratio' },
    { key: 'f', label: 'f', type: 'sci3' },
    { key: 'f_star', label: 'f*', type: 'sci3' },
    { key: 'rel', label: 'rel. err', type: 'sci1' },
    { key: 'viol', label: 'viol', type: 'sci1' },
    { key: 'status', label: 'status', type: 'text' },
    { key: 'fd', label: 'notes', type: 'fd' }
  ];

  var head = document.getElementById('hs-head');
  var ths = COLS.map(function (c) {
    var th = el('th');
    th.scope = 'col';
    if (c.type !== 'text' && c.type !== 'outcome' && c.type !== 'fd') th.className = 'num';
    var btn = el('button', 'hs-sort', c.label);
    btn.type = 'button';
    btn.addEventListener('click', function () {
      if (state.sortKey === c.key) state.sortDir = -state.sortDir;
      else { state.sortKey = c.key; state.sortDir = 1; }
      update();
    });
    th.appendChild(btn);
    head.appendChild(th);
    return th;
  });

  var body = document.getElementById('hs-body');
  var rows = {};
  probs.forEach(function (p) {
    var tr = el('tr');
    tr.id = 'tp-' + p.id;
    COLS.forEach(function (c) {
      var td = el('td');
      var v = p[c.key];
      switch (c.type) {
        case 'int': td.className = 'num'; td.textContent = fmtInt(v); break;
        case 'ratio': td.className = 'num'; td.textContent = fmtRatio(v); break;
        case 'sci3': td.className = 'num'; td.textContent = fmtSci(v, 3); break;
        case 'sci1': td.className = 'num'; td.textContent = fmtSci(v, 1); break;
        case 'outcome':
          td.appendChild(el('span', 'hs-key ' + v));
          td.appendChild(document.createTextNode(LABEL[v]));
          break;
        case 'fd': td.textContent = v ? 'fd' : ''; break;
        default: td.textContent = v;
      }
      tr.appendChild(td);
    });
    rows[p.id] = tr;
    body.appendChild(tr);
  });

  function compare(a, b) {
    var k = state.sortKey, x = a[k], y = b[k];
    if (k === 'outcome') { x = OUTCOMES.indexOf(x); y = OUTCOMES.indexOf(y); }
    if (x === null || x === undefined) return 1;    // (missing values last, either way)
    if (y === null || y === undefined) return -1;
    var r = typeof x === 'string' ? x.localeCompare(y) : (x < y ? -1 : x > y ? 1 : 0);
    return (r || a.id - b.id) * state.sortDir;
  }

  // ---------- update everything below the filters ----------
  function update() {
    chips.forEach(function (b) { b.setAttribute('aria-checked', String(b.dataset.outcome === state.outcome)); });
    var visible = probs.filter(matches);
    var sorted = visible.slice().sort(compare);
    probs.forEach(function (p) { rows[p.id].hidden = true; });
    sorted.forEach(function (p) { rows[p.id].hidden = false; body.appendChild(rows[p.id]); });
    ths.forEach(function (th, i) {
      th.setAttribute('aria-sort', COLS[i].key === state.sortKey ? (state.sortDir > 0 ? 'ascending' : 'descending') : 'none');
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
  update();
})();
