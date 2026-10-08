/* Shared helpers for every page: API calls with a loading bar and an error banner,
 * number formatting, sortable tables, CSV export, cancer-type selects and the gene
 * autocomplete. Everything is exposed on window.App. */
(function (global) {
  'use strict';

  // Relative base so the app also works when mounted under a path prefix.
  const API_BASE = 'api';

  const COMPARTMENT_LABEL = { malignant: 'Malignant', non_malignant: 'Non-malignant' };

  // Human labels for every column that can show up in a table (gene_summary and friends).
  const COLUMN_LABELS = {
    dataset: 'Dataset',
    cancer_type: 'Cancer',
    cancer_name: 'Cancer name',
    gene: 'Gene',
    rank: '#',
    value: 'Value',
    value_q: 'q (this metric)',
    patient: 'Patient',
    compartment: 'Compartment',
    cell_type: 'Cell type',
    n_cells: 'Cells',
    n_patients: 'Patients',
    mean: 'Mean (CP10K)',
    detection_rate: 'Detection rate',
    sd: 'SD',
    cv: 'CV',
    cv2: 'CV²',
    mad: 'MAD',
    fano: 'Fano',
    cv2_resid: 'Mean-corr. CV²',
    cv_ds: 'CV (downsampled)',
    cv2_resid_ds: 'Mean-corr. CV² (downsampled)',
    mean_malignant: 'Mean Mal',
    mean_non_malignant: 'Mean Non-mal',
    detection_rate_malignant: 'Detection Mal',
    detection_rate_non_malignant: 'Detection Non-mal',
    sd_malignant: 'SD Mal',
    sd_non_malignant: 'SD Non-mal',
    cv_malignant: 'CV Mal',
    cv_non_malignant: 'CV Non-mal',
    cv2_malignant: 'CV² Mal',
    cv2_non_malignant: 'CV² Non-mal',
    mad_malignant: 'MAD Mal',
    mad_non_malignant: 'MAD Non-mal',
    fano_malignant: 'Fano Mal',
    fano_non_malignant: 'Fano Non-mal',
    cv2_resid_malignant: 'Mean-corr. CV² Mal',
    cv2_resid_non_malignant: 'Mean-corr. CV² Non-mal',
    n_patients_malignant: 'Patients Mal',
    n_patients_non_malignant: 'Patients Non-mal',
    n_pairs: 'Paired patients',
    log2fc_cv_median: 'log2(CV Mal/Non-mal) median',
    log2fc_cv_mean: 'log2(CV Mal/Non-mal) mean',
    log2fc_cv_of_medians: 'log2(CV Mal/Non-mal) of medians',
    wilcoxon_stat: 'Wilcoxon W',
    p_value: 'p (CV)',
    q_value: 'q (CV)',
    rank_biserial: 'Rank-biserial',
    log2fc_mad_median: 'log2(MAD Mal/Non-mal)',
    p_value_mad: 'p (MAD)',
    q_value_mad: 'q (MAD)',
    log2fc_fano_median: 'log2(Fano Mal/Non-mal)',
    p_value_fano: 'p (Fano)',
    q_value_fano: 'q (Fano)',
    delta_cv2_resid: 'Δ mean-corr. CV²',
    expressed: 'Expressed',
    cv_pooled_malignant: 'CV pooled Mal',
    cv_pooled_non_malignant: 'CV pooled Non-mal',
    depth2_pseudobulk: 'DEPTH2 (pseudo-bulk)',
    depth_pseudobulk: 'DEPTH (pseudo-bulk)',
    depth_cell_median: 'DEPTH (cell median)',
    depth2_cell_median: 'DEPTH2 (cell median)',
    median: 'Median',
  };

  // Ranking metrics shared by the top-genes and heatmap pages (same vocabulary as the API).
  const RANK_METRICS = [
    { key: 'log2fc_cv_median', label: 'Differential noise log2(CV Mal/CV Non-mal)', hasQ: true, diverging: true },
    { key: 'log2fc_mad_median', label: 'log2(MAD Mal/MAD Non-mal)', hasQ: true, diverging: true },
    { key: 'log2fc_fano_median', label: 'log2(Fano Mal/Fano Non-mal)', hasQ: true, diverging: true },
    { key: 'delta_cv2_resid', label: 'Δ mean-corrected CV² (Mal - Non-mal)', hasQ: true, diverging: true },
    { key: 'cv_malignant', label: 'CV, Malignant', compartment: 'malignant' },
    { key: 'cv_non_malignant', label: 'CV, Non-malignant', compartment: 'non_malignant' },
    { key: 'mad_malignant', label: 'MAD, Malignant', compartment: 'malignant' },
    { key: 'mad_non_malignant', label: 'MAD, Non-malignant', compartment: 'non_malignant' },
    { key: 'fano_malignant', label: 'Fano, Malignant', compartment: 'malignant' },
    { key: 'fano_non_malignant', label: 'Fano, Non-malignant', compartment: 'non_malignant' },
    { key: 'cv2_resid_malignant', label: 'Mean-corrected CV² (LOWESS residual), Malignant', compartment: 'malignant', diverging: true },
    { key: 'cv2_resid_non_malignant', label: 'Mean-corrected CV² (LOWESS residual), Non-malignant', compartment: 'non_malignant', diverging: true },
  ];

  class ApiError extends Error {
    constructor(message, status) {
      super(message);
      this.name = 'ApiError';
      this.status = status;
    }
  }

  // ------------------------------------------------------------------ loading / errors
  let pending = 0;
  function setLoading(on) {
    pending = Math.max(0, pending + (on ? 1 : -1));
    const bar = document.getElementById('loading');
    if (bar) bar.hidden = pending === 0;
  }

  function showError(message) {
    const banner = document.getElementById('error-banner');
    const text = document.getElementById('error-text');
    if (!banner || !text) return;
    text.textContent = message;
    banner.hidden = false;
  }

  function clearError() {
    const banner = document.getElementById('error-banner');
    if (banner) banner.hidden = true;
  }

  function buildQuery(params) {
    const q = new URLSearchParams();
    Object.entries(params || {}).forEach(([k, v]) => {
      if (v !== undefined && v !== null && v !== '') q.set(k, String(v));
    });
    const s = q.toString();
    return s ? '?' + s : '';
  }

  function apiUrl(path, params) {
    return API_BASE + path + buildQuery(params);
  }

  function errorDetail(body, resp) {
    const detail = body && body.detail;
    if (typeof detail === 'string') return detail;
    if (Array.isArray(detail)) return detail.map((d) => (d.loc ? d.loc.slice(1).join('.') + ': ' : '') + d.msg).join('; ');
    return resp.status + ' ' + resp.statusText;
  }

  /* GET a JSON endpoint. Non-2xx responses raise ApiError with the server's detail
   * message; unless opts.silent the message also goes to the error banner. */
  async function apiGet(path, params, opts) {
    opts = opts || {};
    setLoading(true);
    let err = null;
    try {
      const resp = await fetch(apiUrl(path, params));
      let body = null;
      try {
        body = await resp.json();
      } catch (e) {
        body = null;
      }
      if (!resp.ok) throw new ApiError(errorDetail(body, resp), resp.status);
      return body;
    } catch (e) {
      err = e instanceof ApiError ? e : new ApiError('Request failed: ' + e.message, 0);
      throw err;
    } finally {
      setLoading(false);
      if (err && !opts.silent) showError(err.message);
    }
  }

  // ------------------------------------------------------------------ formatting
  function isNum(x) {
    return typeof x === 'number' && Number.isFinite(x);
  }

  /* Three significant digits; integers are shown in full; tiny/huge values in e-notation. */
  function fmt(x, digits) {
    digits = digits || 3;
    if (x === null || x === undefined) return '–';
    if (typeof x === 'boolean') return x ? 'yes' : 'no';
    if (typeof x !== 'number') return String(x);
    if (!Number.isFinite(x)) return '–';
    if (Number.isInteger(x) && Math.abs(x) < 1e9) return x.toLocaleString('en-US');
    const a = Math.abs(x);
    if (a !== 0 && (a < 1e-3 || a >= 1e6)) return x.toExponential(digits - 1);
    return String(parseFloat(x.toPrecision(digits)));
  }

  /* p- and q-values in scientific notation. */
  function fmtP(x) {
    if (!isNum(x)) return '–';
    return x.toExponential(2);
  }

  function fmtInt(x) {
    return isNum(x) ? Math.round(x).toLocaleString('en-US') : '–';
  }

  function isPColumn(key) {
    return /^(p_value|q_value)(_|$)/.test(key) || key === 'p' || key === 'q' || key === 'value_q';
  }

  function columnLabel(key) {
    return COLUMN_LABELS[key] || key;
  }

  /* Column definitions for a table from raw keys: label, formatter and a class per cell. */
  function columns(keys, overrides) {
    overrides = overrides || {};
    return keys.map((key) => {
      const col = { key: key, label: columnLabel(key) };
      if (isPColumn(key)) col.format = fmtP;
      if (key === 'compartment') col.format = (v) => COMPARTMENT_LABEL[v] || v;
      if (key === 'compartment') col.className = (v) => (v === 'malignant' ? 'mal' : v === 'non_malignant' ? 'non' : '');
      if (/^q_value/.test(key) || key === 'value_q') col.className = (v) => (isNum(v) && v < 0.05 ? 'sig' : '');
      return Object.assign(col, overrides[key] || {});
    });
  }

  // ------------------------------------------------------------------ DOM helpers
  function el(tag, attrs, children) {
    const node = document.createElement(tag);
    Object.entries(attrs || {}).forEach(([k, v]) => {
      if (k === 'class') node.className = v;
      else if (k === 'text') node.textContent = v;
      else if (k === 'html') node.innerHTML = v;
      else if (k.startsWith('on')) node.addEventListener(k.slice(2), v);
      else if (v !== null && v !== undefined && v !== false) node.setAttribute(k, v === true ? '' : v);
    });
    (children || []).forEach((c) => node.appendChild(typeof c === 'string' ? document.createTextNode(c) : c));
    return node;
  }

  function $(id) {
    return document.getElementById(id);
  }

  function show(node, visible) {
    if (typeof node === 'string') node = $(node);
    if (node) node.hidden = !visible;
  }

  function debounce(fn, ms) {
    let t = null;
    return function () {
      const args = arguments;
      clearTimeout(t);
      t = setTimeout(() => fn.apply(null, args), ms);
    };
  }

  function getParams() {
    return new URLSearchParams(window.location.search);
  }

  /* Keep the current analysis in the URL so a view can be shared or reloaded. */
  function setParams(params) {
    const q = buildQuery(params);
    const url = window.location.pathname + q;
    window.history.replaceState(null, '', url);
  }

  // ------------------------------------------------------------------ tables
  function compareValues(a, b) {
    const an = a === null || a === undefined || (typeof a === 'number' && !Number.isFinite(a));
    const bn = b === null || b === undefined || (typeof b === 'number' && !Number.isFinite(b));
    if (an && bn) return 0;
    if (an) return 1; // nulls always last
    if (bn) return -1;
    if (typeof a === 'number' && typeof b === 'number') return a - b;
    return String(a).localeCompare(String(b), undefined, { numeric: true });
  }

  /* Render a sortable table. cols: [{key, label, format?, className?}]. Clicking a header
   * sorts by that column (numbers numerically, nulls last); a second click flips the order. */
  function renderTable(container, cols, rows, options) {
    if (typeof container === 'string') container = $(container);
    options = options || {};
    let sortKey = options.sortKey || null;
    let sortDir = options.sortDir || 'desc';
    const maxRows = options.maxRows || Infinity;

    container.innerHTML = '';
    if (!rows || !rows.length) {
      container.appendChild(el('div', { class: 'placeholder', text: options.emptyText || 'No rows.' }));
      return;
    }
    const table = el('table', { class: 'data-table' });
    const thead = el('thead');
    const headRow = el('tr');
    cols.forEach((c) => {
      const th = el('th', { text: c.label, title: c.title || c.key, 'data-key': c.key });
      th.addEventListener('click', () => {
        if (sortKey === c.key) sortDir = sortDir === 'asc' ? 'desc' : 'asc';
        else {
          sortKey = c.key;
          sortDir = typeof rows[0][c.key] === 'string' ? 'asc' : 'desc';
        }
        body();
      });
      headRow.appendChild(th);
    });
    thead.appendChild(headRow);
    table.appendChild(thead);
    const tbody = el('tbody');
    table.appendChild(tbody);

    function body() {
      let view = rows.slice();
      if (sortKey) {
        const sign = sortDir === 'asc' ? 1 : -1;
        view.sort((a, b) => {
          const r = compareValues(a[sortKey], b[sortKey]);
          const an = a[sortKey] === null || a[sortKey] === undefined;
          const bn = b[sortKey] === null || b[sortKey] === undefined;
          return an || bn ? r : r * sign;
        });
      }
      headRow.querySelectorAll('th').forEach((th) => {
        th.classList.remove('sorted-asc', 'sorted-desc');
        if (th.dataset.key === sortKey) th.classList.add(sortDir === 'asc' ? 'sorted-asc' : 'sorted-desc');
      });
      tbody.innerHTML = '';
      view.slice(0, maxRows).forEach((row) => {
        const tr = el('tr');
        cols.forEach((c) => {
          const raw = row[c.key];
          const text = c.format ? c.format(raw, row) : fmt(raw);
          const td = el('td', { text: text });
          if (typeof raw === 'string') td.classList.add('text');
          const cls = c.className ? c.className(raw, row) : '';
          if (cls) td.classList.add(cls);
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      });
    }
    body();
    const wrap = el('div', { class: 'table-wrap' }, [table]);
    container.appendChild(wrap);
    if (rows.length > maxRows) container.appendChild(el('p', { class: 'table-caption', text: 'Showing the first ' + maxRows + ' of ' + rows.length + ' rows; download the CSV for all of them.' }));
  }

  // ------------------------------------------------------------------ CSV (client side)
  function csvEscape(v) {
    if (v === null || v === undefined) return '';
    const s = String(v);
    return /[",\n\r]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }

  function downloadCsv(filename, keys, rows) {
    const lines = [keys.map(csvEscape).join(',')];
    rows.forEach((r) => lines.push(keys.map((k) => csvEscape(r[k])).join(',')));
    const blob = new Blob([lines.join('\n') + '\n'], { type: 'text/csv;charset=utf-8' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  }

  // ------------------------------------------------------------------ cancer types
  let cancerTypesPromise = null;
  function loadCancerTypes() {
    if (!cancerTypesPromise) {
      cancerTypesPromise = apiGet('/cancer-types').catch((e) => {
        cancerTypesPromise = null;
        throw e;
      });
    }
    return cancerTypesPromise;
  }

  function cancerLabel(c) {
    return c.cancer_name ? c.cancer_type + ' – ' + c.cancer_name : c.cancer_type;
  }

  function fillCancerSelect(select, cancers, selected) {
    if (typeof select === 'string') select = $(select);
    select.innerHTML = '';
    cancers.forEach((c) => select.appendChild(el('option', { value: c.cancer_type, text: cancerLabel(c) })));
    if (selected && cancers.some((c) => c.cancer_type === selected)) select.value = selected;
  }

  function fillSelect(select, items, selected) {
    if (typeof select === 'string') select = $(select);
    select.innerHTML = '';
    items.forEach((it) => select.appendChild(el('option', { value: it.key, text: it.label })));
    if (selected) select.value = selected;
  }

  // ------------------------------------------------------------------ gene autocomplete
  /* Suggestions from /api/genes/search under a text input. onSelect(gene) fires when a
   * suggestion is chosen with the mouse or with Enter; arrow keys move the highlight. */
  function attachGeneAutocomplete(input, onSelect) {
    const list = el('ul', { class: 'autocomplete-list', role: 'listbox' });
    list.hidden = true;
    input.setAttribute('autocomplete', 'off');
    input.parentNode.insertBefore(list, input.nextSibling);
    let items = [];
    let active = -1;
    let seq = 0;

    function close() {
      list.hidden = true;
      list.innerHTML = '';
      items = [];
      active = -1;
    }
    function render() {
      list.innerHTML = '';
      items.forEach((it, i) => {
        const li = el('li', { 'data-index': i, role: 'option' }, [el('span', { text: it.gene }), el('span', { class: 'meta', text: it.n_cancer_types + (it.n_cancer_types === 1 ? ' cancer type' : ' cancer types') })]);
        if (i === active) li.classList.add('active');
        list.appendChild(li);
      });
      list.hidden = items.length === 0;
    }
    function choose(i) {
      if (i < 0 || i >= items.length) return;
      const gene = items[i].gene;
      input.value = gene;
      close();
      onSelect(gene);
    }
    const search = debounce(async () => {
      const q = input.value.trim();
      if (!q) {
        close();
        return;
      }
      const mine = ++seq;
      try {
        const res = await apiGet('/genes/search', { q: q, limit: 15 }, { silent: true });
        if (mine !== seq) return;
        items = res.results || [];
        active = items.length ? 0 : -1;
        render();
      } catch (e) {
        close();
      }
    }, 120);

    input.addEventListener('input', search);
    input.addEventListener('keydown', (e) => {
      if (list.hidden) return;
      if (e.key === 'ArrowDown') {
        e.preventDefault();
        active = Math.min(items.length - 1, active + 1);
        render();
      } else if (e.key === 'ArrowUp') {
        e.preventDefault();
        active = Math.max(0, active - 1);
        render();
      } else if (e.key === 'Enter') {
        if (active >= 0) {
          e.preventDefault();
          choose(active);
        }
      } else if (e.key === 'Escape') {
        close();
      }
    });
    input.addEventListener('blur', () => setTimeout(close, 150));
    list.addEventListener('mousedown', (e) => {
      const li = e.target.closest('li');
      if (li) {
        e.preventDefault();
        choose(Number(li.dataset.index));
      }
    });
    return { close: close };
  }

  // ------------------------------------------------------------------ page chrome
  function initChrome() {
    const dismiss = $('error-dismiss');
    if (dismiss) dismiss.addEventListener('click', clearError);
    const version = $('footer-version');
    if (version) {
      apiGet('/health', null, { silent: true })
        .then((h) => {
          const parts = [];
          if (h.build && h.build.sctcna_version) parts.push('scTCNA v' + h.build.sctcna_version);
          if (h.n_datasets) parts.push(h.n_datasets + (h.n_datasets === 1 ? ' dataset' : ' datasets'));
          if (h.n_genes) parts.push(fmtInt(h.n_genes) + ' genes');
          if (h.build && h.build.built) parts.push('built ' + h.build.built.slice(0, 10));
          version.textContent = parts.join(' · ');
          if (h.status !== 'ok') showError('The server has no data tables loaded yet (' + (h.detail || h.status) + ').');
        })
        .catch(() => {});
    }
  }
  document.addEventListener('DOMContentLoaded', initChrome);

  global.App = {
    ApiError: ApiError,
    COMPARTMENT_LABEL: COMPARTMENT_LABEL,
    COLUMN_LABELS: COLUMN_LABELS,
    RANK_METRICS: RANK_METRICS,
    apiGet: apiGet,
    apiUrl: apiUrl,
    buildQuery: buildQuery,
    showError: showError,
    clearError: clearError,
    fmt: fmt,
    fmtP: fmtP,
    fmtInt: fmtInt,
    isNum: isNum,
    columnLabel: columnLabel,
    columns: columns,
    el: el,
    $: $,
    show: show,
    debounce: debounce,
    getParams: getParams,
    setParams: setParams,
    renderTable: renderTable,
    downloadCsv: downloadCsv,
    loadCancerTypes: loadCancerTypes,
    cancerLabel: cancerLabel,
    fillCancerSelect: fillCancerSelect,
    fillSelect: fillSelect,
    attachGeneAutocomplete: attachGeneAutocomplete,
  };
})(window);
