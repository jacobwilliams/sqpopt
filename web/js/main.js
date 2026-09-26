// SQPOPT user guide: theme, navigation, code blocks, math, and option filtering.
(function () {
  'use strict';

  var root = document.documentElement;

  // ---------- theme toggle ----------
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

  // ---------- mobile navigation ----------
  var navToggle = document.querySelector('.nav-toggle');
  function setNav(open) {
    document.body.classList.toggle('nav-open', open);
    navToggle.setAttribute('aria-expanded', String(open));
  }
  navToggle.addEventListener('click', function () {
    setNav(!document.body.classList.contains('nav-open'));
  });
  document.querySelectorAll('.toc a').forEach(function (a) {
    a.addEventListener('click', function () { setNav(false); });
  });

  // ---------- heading anchors ----------
  document.querySelectorAll('.content h2[id], .content h3[id]').forEach(function (h) {
    var a = document.createElement('a');
    a.className = 'anchor';
    a.href = '#' + h.id;
    a.setAttribute('aria-label', 'Link to this section');
    a.textContent = '#';
    h.appendChild(a);
  });

  // ---------- scrollspy ----------
  var tocLinks = Array.prototype.slice.call(document.querySelectorAll('.toc a[href^="#"]'));
  var targets = tocLinks
    .map(function (a) { return document.getElementById(a.getAttribute('href').slice(1)); })
    .filter(Boolean);
  function updateActive() {
    var offset = 90, current = targets[0];
    for (var i = 0; i < targets.length; i++) {
      if (targets[i].getBoundingClientRect().top - offset <= 0) current = targets[i];
    }
    // at the very bottom, activate the last section
    if (window.innerHeight + window.scrollY >= document.body.scrollHeight - 4) current = targets[targets.length - 1];
    tocLinks.forEach(function (a) {
      a.classList.toggle('active', current && a.getAttribute('href') === '#' + current.id);
    });
  }
  var ticking = false;
  window.addEventListener('scroll', function () {
    if (!ticking) { ticking = true; requestAnimationFrame(function () { updateActive(); ticking = false; }); }
  }, { passive: true });
  updateActive();

  // ---------- syntax highlighting & copy buttons ----------
  if (window.hljs) {
    document.querySelectorAll('pre code').forEach(function (el) {
      if (el.classList.contains('language-toml')) { el.classList.remove('language-toml'); el.classList.add('language-ini'); }
      if (el.classList.contains('language-bash')) { el.classList.add('nohighlight'); return; }
      try { hljs.highlightElement(el); } catch (e) {}
    });
  }
  document.querySelectorAll('pre').forEach(function (pre) {
    var btn = document.createElement('button');
    btn.className = 'copy-btn';
    btn.type = 'button';
    btn.textContent = 'Copy';
    btn.addEventListener('click', function () {
      var text = pre.querySelector('code').innerText;
      var done = function () {
        btn.textContent = 'Copied';
        btn.classList.add('copied');
        setTimeout(function () { btn.textContent = 'Copy'; btn.classList.remove('copied'); }, 1500);
      };
      if (navigator.clipboard) navigator.clipboard.writeText(text).then(done, function () {});
    });
    pre.appendChild(btn);
  });

  // ---------- math ----------
  if (window.renderMathInElement) {
    renderMathInElement(document.querySelector('.content'), {
      delimiters: [
        { left: '\\[', right: '\\]', display: true },
        { left: '\\(', right: '\\)', display: false }
      ],
      ignoredTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code'],
      throwOnError: false
    });
  }

  // ---------- option filter ----------
  var input = document.getElementById('opt-filter');
  var count = document.getElementById('opt-filter-count');
  var wraps = Array.prototype.slice.call(document.querySelectorAll('.opt-table')).map(function (t) {
    return t.closest('.table-wrap');
  });
  var rows = Array.prototype.slice.call(document.querySelectorAll('.opt-table tbody tr'));
  rows.forEach(function (r) { r._text = r.textContent.toLowerCase(); });

  function applyFilter() {
    var q = input.value.trim().toLowerCase();
    var terms = q.split(/\s+/).filter(Boolean);
    var shown = 0;
    rows.forEach(function (r) {
      var match = terms.every(function (t) { return r._text.indexOf(t) !== -1; });
      r.classList.toggle('filter-hidden', !match);
      if (match) shown++;
    });
    wraps.forEach(function (w) {
      var any = w.querySelector('tbody tr:not(.filter-hidden)');
      w.classList.toggle('filter-hidden', !any);
      // also hide a group heading / "set as" line that only introduces this table
      var prev = w.previousElementSibling;
      while (prev && (prev.classList.contains('opt-group') || prev.classList.contains('set-on'))) {
        if (prev.classList.contains('opt-group')) prev.classList.toggle('filter-hidden', !any);
        prev = prev.previousElementSibling;
      }
    });
    count.textContent = terms.length ? shown + (shown === 1 ? ' match' : ' matches') : rows.length + ' entries';
  }
  input.addEventListener('input', applyFilter);
  input.addEventListener('keydown', function (e) {
    if (e.key === 'Escape') { input.value = ''; applyFilter(); input.blur(); }
  });
  // "/" focuses the filter, as on many docs sites
  document.addEventListener('keydown', function (e) {
    if (e.key === '/' && document.activeElement !== input && !/INPUT|TEXTAREA/.test(document.activeElement.tagName)) {
      e.preventDefault();
      input.focus();
      input.scrollIntoView({ block: 'nearest' });
    }
  });
  applyFilter();
})();
