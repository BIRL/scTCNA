/* Gene Noise Analysis page: pan-cancer and cancer-specific views with per-patient detail. */
(function () {
  'use strict';
  const $ = App.$;
  const fmt = App.fmt;
  const fmtP = App.fmtP;

  const MAX_GENES = 50;
  const COLLAPSE_ABOVE = 8; // with more genes, every metric gets its own collapsible card
  const MIN_PAIRS = 6; // patients with both cell groups needed for the paired test

  // gene_summary columns in display order; the CSV downloads keep this layout (one row per gene or cancer type)
  const SUMMARY_KEYS = [
    'cancer_type', 'cancer_name', 'dataset', 'gene',
    'mean_malignant', 'mean_non_malignant', 'sd_malignant', 'sd_non_malignant',
    'detection_rate_malignant', 'detection_rate_non_malignant',
    'cv_malignant', 'cv_non_malignant', 'cv2_malignant', 'cv2_non_malignant',
    'mad_malignant', 'mad_non_malignant', 'fano_malignant', 'fano_non_malignant',
    'cv2_resid_malignant', 'cv2_resid_non_malignant',
    'n_patients_malignant', 'n_patients_non_malignant', 'n_pairs',
    'log2fc_cv_median', 'log2fc_cv_mean', 'log2fc_cv_of_medians', 'wilcoxon_stat', 'p_value', 'q_value', 'rank_biserial',
    'log2fc_mad_median', 'p_value_mad', 'q_value_mad',
    'log2fc_fano_median', 'p_value_fano', 'q_value_fano',
    'delta_cv2_resid', 'expressed', 'cv_pooled_malignant', 'cv_pooled_non_malignant',
  ];
  const PATIENT_KEYS = ['patient', 'compartment', 'n_cells', 'mean', 'detection_rate', 'cv', 'mad', 'fano', 'cv2_resid', 'cv_ds', 'cv2_resid_ds'];
  const CELLTYPE_KEYS = ['cell_type', 'n_patients', 'mean', 'detection_rate', 'cv', 'mad', 'fano', 'cv2_resid'];
  const METRIC_LABELS = {
    cv: 'CV = SD / mean',
    mad: 'MAD',
    fano: 'Fano',
    cv2_resid: 'mean-corrected CV² (LOWESS residual)',
    mean: 'mean (CP10K)',
    detection_rate: 'detection rate',
  };
  // per-patient tables also hold the depth-downsampled version of these two metrics
  const DOWNSAMPLED = { cv: 'cv_ds', cv2_resid: 'cv2_resid_ds' };

  const CELL_GROUPS = [['malignant', 'Malignant cells'], ['non_malignant', 'Non-malignant cells']];
  const POOLED_TIP = 'Pooling patients inflates malignant CV; the row above is the median across patients.';
  const PAIRED_TIP = 'The paired Wilcoxon test needs the gene in at least 10% of cells of both cell groups (median over patients) and at least ' + MIN_PAIRS + ' patients with both cell groups.';
  // Transposed summary tables: metrics in rows (`metric` rows are split by cell group), genes or cancer types in columns.
  const TABLE_ROWS = [
    { group: 'Mean (CP10K)', metric: 'mean' },
    { group: 'SD', metric: 'sd' },
    { group: 'CV = SD / mean (median across patients)', metric: 'cv' },
    { group: 'CV, all cells pooled (ignores patient)', metric: 'cv_pooled', title: POOLED_TIP },
    { group: 'CV²', metric: 'cv2' },
    { group: 'MAD', metric: 'mad' },
    { group: 'Fano', metric: 'fano' },
    { group: 'Detection rate', metric: 'detection_rate' },
    { label: 'Paired test', value: pairedTest, title: PAIRED_TIP },
    { key: 'log2fc_cv_median', label: 'log2 CV ratio (malignant / non-malignant), median over patients' },
    { key: 'q_value' },
    { group: 'Mean-corrected CV² (LOWESS residual)', metric: 'cv2_resid' },
    { group: 'Patients', metric: 'n_patients' },
    { key: 'n_pairs' },
    { key: 'log2fc_cv_mean' },
    { key: 'log2fc_cv_of_medians' },
    { key: 'wilcoxon_stat' },
    { key: 'p_value' },
    { key: 'rank_biserial' },
    { key: 'log2fc_mad_median' },
    { key: 'p_value_mad' },
    { key: 'q_value_mad' },
    { key: 'log2fc_fano_median' },
    { key: 'p_value_fano' },
    { key: 'q_value_fano' },
    { key: 'delta_cv2_resid' },
    { key: 'expressed', label: 'Expressed (≥ 10% of cells in both cell groups)' },
  ];

  const state = {
    mode: 'pan',
    cancerTypes: [],
    panSelected: new Set(),
    gene: null,
    rows: [], // gene_summary rows of the current gene, one per cancer type
    detailCancer: null,
    patients: null,
    celltypes: null,
    csRows: [],
    csCancer: null,
    lookupSeq: 0,
  };

  function pick(rows, key) {
    return rows.map((r) => (App.isNum(r[key]) ? r[key] : null));
  }

  function compartmentSeries(rows, metric, withError) {
    return ['malignant', 'non_malignant'].map((comp) => ({
      key: comp,
      values: pick(rows, metric + '_' + comp),
      error: withError ? pick(rows, 'sd_' + comp) : null,
    }));
  }

  function parseGenes(text) {
    const seen = new Set();
    const out = [];
    (text || '').split(/[\s,;]+/).forEach((tok) => {
      const key = tok.trim().toUpperCase();
      if (key && !seen.has(key)) {
        seen.add(key);
        out.push(tok.trim());
      }
    });
    return out;
  }

  function cancerName(ct) {
    const c = state.cancerTypes.find((x) => x.cancer_type === ct);
    return c ? App.cancerLabel(c) : ct;
  }

  function sigColors(rows, qKey) {
    return rows.map((r) => (App.isNum(r[qKey]) && r[qKey] < 0.05 ? Plots.COLORS.significant : Plots.COLORS.neutral));
  }

  const SIG_LEGEND = [
    { name: 'q < 0.05', color: Plots.COLORS.significant },
    { name: 'not significant', color: Plots.COLORS.neutral },
  ];

  /* "not tested" over every category whose differential-noise value is missing. */
  function notTestedNotes(categories, values) {
    return categories
      .map((c, i) => (App.isNum(values[i]) ? null : { x: c, y: 0, text: 'not tested', textangle: -90, yanchor: 'bottom', showarrow: false, font: { size: 11, color: '#6b7a89', weight: 'bold' } }))
      .filter(Boolean);
  }

  function pairedTest(r) {
    if (!r.expressed) return 'no, detected in < 10% of cells';
    if (App.isNum(r.p_value)) return 'yes';
    return App.isNum(r.n_pairs) && r.n_pairs < MIN_PAIRS ? 'no, ' + r.n_pairs + ' paired patients (' + MIN_PAIRS + ' needed)' : 'no';
  }

  /* Metric rows of the transposed table for `rows` (one column per row). */
  function metricRows(rows) {
    const out = [];
    TABLE_ROWS.forEach((spec) => {
      if (spec.metric) {
        CELL_GROUPS.forEach(([comp, label]) => {
          const key = spec.metric + '_' + comp;
          const col = App.columns([key])[0];
          out.push({ group: spec.group, groupTitle: spec.title, label: label, title: spec.title || key, values: rows.map((r) => r[key]), format: col.format, className: col.className });
        });
      } else if (spec.value) {
        out.push({ label: spec.label, title: spec.title, values: rows.map(spec.value) });
      } else {
        const col = App.columns([spec.key])[0];
        out.push({ label: spec.label || col.label, title: spec.key, values: rows.map((r) => r[spec.key]), format: col.format, className: col.className });
      }
    });
    return out;
  }

  function tcnaLink(genes, cancerTypes, analysisType) {
    const href = App.tcnaGeneUrl(genes, cancerTypes, analysisType);
    return href ? { href: href, text: 'View in TCNA ↗' } : null;
  }

  function setLink(id, href) {
    const a = $(id);
    a.hidden = !href;
    if (href) a.href = href;
  }

  function diffHover(rows) {
    return rows.map(
      (r) =>
        (r.cancer_name || r.cancer_type) +
        '<br>log2(CV malignant / non-malignant), median: ' + fmt(r.log2fc_cv_median) +
        '<br>paired patients: ' + fmt(r.n_pairs) +
        '<br>p = ' + fmtP(r.p_value) + ', q = ' + fmtP(r.q_value) +
        '<br>rank-biserial: ' + fmt(r.rank_biserial)
    );
  }

  /* Inline notes about genes that are missing or have no paired test. */
  function geneNotes(where, groups) {
    const lines = [];
    if (groups.unknown && groups.unknown.length) lines.push('Not found in any dataset: ' + groups.unknown.join(', ') + '.');
    if (groups.notIn && groups.notIn.length) lines.push('Not measured in ' + where + ': ' + groups.notIn.join(', ') + '.');
    if (groups.lowDetection && groups.lowDetection.length) lines.push('Detected in < 10% of cells, so no paired test and no differential-noise bar: ' + groups.lowDetection.join(', ') + '.');
    if (groups.fewPairs && groups.fewPairs.length) lines.push('Too few patients with both cell groups for the paired test (' + MIN_PAIRS + ' needed): ' + groups.fewPairs.join(', ') + '.');
    if (groups.failed && groups.failed.length) lines.push('Could not be loaded: ' + groups.failed.join(', ') + '.');
    return lines;
  }

  function setMode(mode) {
    state.mode = mode;
    $('mode-pan').checked = mode === 'pan';
    $('mode-cancer').checked = mode === 'cancer';
    App.show('pan-form', mode === 'pan');
    App.show('cancer-form', mode === 'cancer');
    App.show('pan-results', mode === 'pan' && state.rows.length > 0);
    App.show('cancer-results', mode === 'cancer' && state.csRows.length > 0);
  }

  /* ---- pan-cancer: cancer-type checkboxes ---- */

  function allSelected() {
    return state.cancerTypes.every((c) => state.panSelected.has(c.cancer_type));
  }

  function panParams() {
    return { mode: 'pan', gene: state.gene, cancer_types: allSelected() ? null : state.cancerTypes.filter((c) => state.panSelected.has(c.cancer_type)).map((c) => c.cancer_type).join(',') || 'none' };
  }

  function updatePicker() {
    document.querySelectorAll('#pan-checks input').forEach((box) => {
      box.checked = state.panSelected.has(box.value);
    });
    $('pan-count').textContent = state.panSelected.size + ' of ' + state.cancerTypes.length + ' selected';
  }

  function onSelectionChange() {
    updatePicker();
    if (state.gene) {
      App.setParams(panParams());
      renderPan();
    }
  }

  function buildPicker(initial) {
    const box = $('pan-checks');
    box.innerHTML = '';
    state.cancerTypes.forEach((c) => {
      const input = App.el('input', { type: 'checkbox', value: c.cancer_type });
      input.addEventListener('change', () => {
        if (input.checked) state.panSelected.add(c.cancer_type);
        else state.panSelected.delete(c.cancer_type);
        onSelectionChange();
      });
      box.appendChild(App.el('label', { class: 'check check-chip', title: c.cancer_name || c.cancer_type }, [input, App.el('span', { text: c.cancer_type })]));
    });
    const wanted = initial === null ? null : initial.split(',');
    state.panSelected = new Set(state.cancerTypes.map((c) => c.cancer_type).filter((ct) => !wanted || wanted.includes(ct)));
    updatePicker();
  }

  /* ---- pan-cancer: one gene across cancer types ---- */

  async function runPan(gene) {
    gene = (gene || '').trim();
    if (!gene) return;
    App.clearError();
    let res;
    try {
      res = await App.apiGet('/genes/' + encodeURIComponent(gene));
    } catch (e) {
      return;
    }
    state.gene = res.gene;
    state.rows = res.cancers || [];
    $('gene-input').value = res.gene;
    App.setParams(panParams());
    renderPan();
  }

  function panRows() {
    return state.rows.filter((r) => state.panSelected.has(r.cancer_type));
  }

  const PAN_PLOTS = ['plot-mean', 'plot-cv', 'plot-mad', 'plot-fano', 'plot-cv2resid', 'plot-diff'];

  function renderPan() {
    const rows = panRows();
    const gene = state.gene;
    const selected = state.cancerTypes.map((c) => c.cancer_type).filter((ct) => state.panSelected.has(ct));
    const notIn = selected.filter((ct) => !state.rows.some((r) => r.cancer_type === ct));
    App.show('pan-results', true);
    $('pan-title').textContent = gene + ' across ' + rows.length + (rows.length === 1 ? ' cancer type' : ' cancer types');
    const lowDetection = rows.filter((r) => !r.expressed).map((r) => r.cancer_type);
    const fewPairs = rows.filter((r) => r.expressed && !App.isNum(r.p_value)).map((r) => r.cancer_type);
    App.setNote('pan-note', [
      selected.length ? '' : 'Select at least one cancer type above.',
      notIn.length ? gene + ' is not measured in ' + notIn.join(', ') + '.' : '',
      lowDetection.length ? gene + ' is detected in < 10% of cells in ' + lowDetection.join(', ') + ', so there is no paired test and no differential-noise bar there.' : '',
      fewPairs.length ? 'Too few patients with both cell groups for the paired test (' + MIN_PAIRS + ' needed) in ' + fewPairs.join(', ') + '.' : '',
    ]);
    setLink('pan-tcna', rows.length ? App.tcnaGeneUrl([gene], rows.map((r) => r.cancer_type), 'pan-cancer') : null);
    if (!rows.length) {
      PAN_PLOTS.forEach((id) => Plots.empty(id, selected.length ? gene + ' is not measured in the selected cancer types.' : 'No cancer type selected.'));
      App.renderMetricTable('pan-table', { columns: [], rows: [], emptyText: 'No cancer type to show.' });
      App.show('detail-card', false);
      return;
    }
    const cats = rows.map((r) => r.cancer_type);
    Plots.groupedBars('plot-mean', cats, compartmentSeries(rows, 'mean', true), { ytitle: 'mean (CP10K)', filename: gene + '_mean_sd' });
    Plots.groupedBars('plot-cv', cats, compartmentSeries(rows, 'cv'), { ytitle: 'CV', filename: gene + '_cv' });
    Plots.groupedBars('plot-mad', cats, compartmentSeries(rows, 'mad'), { ytitle: 'MAD', filename: gene + '_mad' });
    Plots.groupedBars('plot-fano', cats, compartmentSeries(rows, 'fano'), { ytitle: 'Fano', filename: gene + '_fano' });
    Plots.groupedBars('plot-cv2resid', cats, compartmentSeries(rows, 'cv2_resid'), { ytitle: 'residual log10 CV²', rangemode: 'normal', filename: gene + '_cv2_resid' });
    const diff = pick(rows, 'log2fc_cv_median');
    Plots.bars('plot-diff', cats, diff, {
      colors: sigColors(rows, 'q_value'),
      legend: SIG_LEGEND,
      annotations: notTestedNotes(cats, diff),
      text: rows.map((r) => (App.isNum(r.n_pairs) && App.isNum(r.log2fc_cv_median) ? 'n=' + r.n_pairs : '')),
      hover: diffHover(rows),
      ytitle: 'log2(CV malignant / CV non-malignant)',
      filename: gene + '_differential_noise',
    }).then(() => {
      const div = $('plot-diff');
      if (div.removeAllListeners) div.removeAllListeners('plotly_click');
      div.on('plotly_click', (ev) => {
        const pt = ev.points && ev.points[0];
        if (pt && pt.x) selectDetailCancer(String(pt.x), true);
      });
    });

    const columns = rows.map((r) => ({ label: r.cancer_type, title: r.cancer_name || r.cancer_type }));
    App.renderMetricTable('pan-table', {
      corner: gene,
      columns: columns,
      rows: [
        { label: 'Bulk TCNA', link: true, values: rows.map((r) => tcnaLink([gene], [r.cancer_type], 'cancer-specific')) },
        { label: 'Cancer name', title: 'cancer_name', values: rows.map((r) => r.cancer_name) },
        { label: 'Dataset', title: 'dataset', values: rows.map((r) => r.dataset) },
      ].concat(metricRows(rows)),
    });

    // per-patient detail: keep the chosen cancer when the new gene has it, else the first cancer where the gene
    // is expressed (per-patient tables only hold expressed genes), else the first one
    App.show('detail-card', true);
    const available = rows.map((r) => r.cancer_type);
    const expressedIn = rows.filter((r) => r.expressed).map((r) => r.cancer_type);
    const keep = available.includes(state.detailCancer) ? state.detailCancer : expressedIn[0] || available[0];
    App.fillCancerSelect('detail-cancer', rows.map((r) => ({ cancer_type: r.cancer_type, cancer_name: r.cancer_name })), keep);
    selectDetailCancer(keep);
  }

  function selectDetailCancer(cancerType, scroll) {
    if (!cancerType) return;
    state.detailCancer = cancerType;
    $('detail-cancer').value = cancerType;
    // scroll only when the user picked a cancer type in a plot, not on every new search
    if (scroll && state.mode === 'pan') $('detail-card').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    loadDetail();
  }

  async function loadDetail() {
    const gene = state.gene;
    const ct = state.detailCancer;
    $('detail-csv').href = App.apiUrl('/download/gene_patients.csv', { gene: gene, cancer_type: ct });
    const [pat, cts] = await Promise.allSettled([
      App.apiGet('/genes/' + encodeURIComponent(gene) + '/patients', { cancer_type: ct }, { silent: true }),
      App.apiGet('/genes/' + encodeURIComponent(gene) + '/celltypes', { cancer_type: ct }, { silent: true }),
    ]);
    if (gene !== state.gene || ct !== state.detailCancer) return; // a newer request took over
    state.patients = pat.status === 'fulfilled' ? pat.value.rows : null;
    state.patientsError = pat.status === 'rejected' ? pat.reason.message : null;
    state.celltypes = cts.status === 'fulfilled' ? cts.value.rows : null;
    state.celltypesError = cts.status === 'rejected' ? cts.reason.message : null;
    renderPaired();
    renderCelltypes();
  }

  function renderPaired() {
    const metric = $('detail-metric').value;
    const dsKey = DOWNSAMPLED[metric];
    $('detail-ds').disabled = !dsKey;
    $('paired-title').textContent = 'Paired per-patient ' + METRIC_LABELS[metric] + ' – ' + state.detailCancer;
    const rows = state.patients;
    $('detail-csv').hidden = !rows;
    if (!rows) {
      Plots.empty('plot-paired', (state.patientsError || 'No per-patient values.').replace('compartments', 'cell groups'));
      $('paired-stats').textContent = '';
      App.renderTable('patients-table', [], [], { emptyText: (state.patientsError || 'No per-patient values (per-patient tables only hold genes expressed in both cell groups).').replace('compartments', 'cell groups') });
      return;
    }
    const byPatient = {};
    rows.forEach((r) => {
      byPatient[r.patient] = byPatient[r.patient] || {};
      byPatient[r.patient][r.compartment] = r;
    });
    const patients = Object.keys(byPatient).sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
    const value = (p, comp, key) => (byPatient[p][comp] && App.isNum(byPatient[p][comp][key]) ? byPatient[p][comp][key] : null);
    const pairs = patients.map((p) => ({ label: p, a: value(p, 'malignant', metric), b: value(p, 'non_malignant', metric) }));
    let overlay = null;
    if (dsKey && $('detail-ds').checked) {
      overlay = { label: METRIC_LABELS[metric] + ' after downsampling', pairs: patients.map((p) => ({ label: p + ' (downsampled)', a: value(p, 'malignant', dsKey), b: value(p, 'non_malignant', dsKey) })) };
    }
    Plots.pairedLines('plot-paired', pairs, { ytitle: METRIC_LABELS[metric], overlay: overlay, name: 'per patient', filename: state.gene + '_' + state.detailCancer + '_paired_' + metric });

    const summary = state.rows.find((r) => r.cancer_type === state.detailCancer) || {};
    const stats = { cv: ['log2fc_cv_median', 'p_value', 'q_value'], mad: ['log2fc_mad_median', 'p_value_mad', 'q_value_mad'], fano: ['log2fc_fano_median', 'p_value_fano', 'q_value_fano'] }[metric];
    let text = patients.length + ' patients';
    if (stats) {
      text += ' · median log2(malignant / non-malignant) = ' + fmt(summary[stats[0]]) + ' · paired Wilcoxon p = ' + fmtP(summary[stats[1]]) + ', q = ' + fmtP(summary[stats[2]]);
      if (metric === 'cv') text += ' · rank-biserial = ' + fmt(summary.rank_biserial);
    } else if (metric === 'cv2_resid') {
      text += ' · median difference malignant − non-malignant = ' + fmt(summary.delta_cv2_resid);
    }
    $('paired-stats').textContent = text;
    App.renderTable('patients-table', App.columns(PATIENT_KEYS), rows, { sortKey: 'patient', sortDir: 'asc' });
  }

  function renderCelltypes() {
    const metric = $('detail-metric').value;
    const rows = state.celltypes;
    $('celltype-title').textContent = METRIC_LABELS[metric].replace(/^./, (c) => c.toUpperCase()) + ' by cell type – ' + state.detailCancer;
    if (!rows) {
      Plots.empty('plot-celltypes', state.celltypesError || 'No cell-type values.');
      App.renderTable('celltypes-table', [], [], { emptyText: state.celltypesError || 'No cell-type values.' });
      return;
    }
    const sorted = rows.slice().sort((a, b) => (a.cell_type === 'Malignant' ? -1 : b.cell_type === 'Malignant' ? 1 : a.cell_type.localeCompare(b.cell_type)));
    Plots.bars('plot-celltypes', sorted.map((r) => r.cell_type), pick(sorted, metric), {
      colors: sorted.map((r) => (r.cell_type === 'Malignant' ? Plots.COLORS.malignant : Plots.COLORS.non_malignant)),
      legend: [
        { name: 'Malignant cells', color: Plots.COLORS.malignant },
        { name: 'Non-malignant lineages', color: Plots.COLORS.non_malignant },
      ],
      hover: sorted.map((r) => r.cell_type + '<br>' + METRIC_LABELS[metric] + ': ' + fmt(r[metric]) + '<br>patients: ' + fmt(r.n_patients) + '<br>mean: ' + fmt(r.mean)),
      ytitle: METRIC_LABELS[metric],
      filename: state.gene + '_' + state.detailCancer + '_celltypes_' + metric,
    });
    App.renderTable('celltypes-table', App.columns(CELLTYPE_KEYS), sorted);
  }

  /* ---- cancer-specific: a gene list within one cancer type ---- */

  function hideCancerResults() {
    state.csRows = [];
    App.show('cancer-results', false);
  }

  async function runCancer(cancerType, genes) {
    App.clearError();
    state.lookupSeq++; // a pending list check would only repeat what this analysis reports
    if (!genes.length) {
      hideCancerResults();
      App.setNote('genes-note', ['Enter at least one gene symbol.']);
      return;
    }
    if (genes.length > MAX_GENES) {
      hideCancerResults();
      App.setNote('genes-note', ['At most ' + MAX_GENES + ' genes per query (the list has ' + genes.length + ').']);
      return;
    }
    const results = await Promise.allSettled(genes.map((g) => App.apiGet('/genes/' + encodeURIComponent(g), null, { silent: true })));
    const rows = [];
    const groups = { unknown: [], notIn: [], lowDetection: [], fewPairs: [], failed: [] };
    results.forEach((r, i) => {
      if (r.status !== 'fulfilled') {
        (r.reason && r.reason.status === 404 ? groups.unknown : groups.failed).push(genes[i]);
        return;
      }
      const row = (r.value.cancers || []).find((c) => c.cancer_type === cancerType);
      if (!row) {
        groups.notIn.push(r.value.gene);
        return;
      }
      rows.push(row);
      if (!row.expressed) groups.lowDetection.push(row.gene);
      else if (!App.isNum(row.p_value)) groups.fewPairs.push(row.gene);
    });
    const notes = geneNotes(cancerName(cancerType), groups);
    if (!rows.length) {
      hideCancerResults();
      App.setNote('genes-note', notes.concat(['None of the requested genes is measured in ' + cancerName(cancerType) + '.']));
      return;
    }
    App.setNote('genes-note', []);
    state.csRows = rows;
    state.csCancer = cancerType;
    App.setParams({ mode: 'cancer', cancer_type: cancerType, genes: rows.map((r) => r.gene).join(',') });
    renderCancer(notes);
  }

  /* One collapsible card per metric when many genes are shown (first card open), else all open in a grid. */
  function layoutCards(collapsible) {
    const grid = $('cs-grid');
    grid.classList.toggle('collapsible', collapsible);
    grid.querySelectorAll('details.metric-card').forEach((d, i) => {
      d.open = !collapsible || i === 0;
    });
  }

  function renderCancer(notes) {
    const rows = state.csRows;
    const ct = state.csCancer;
    const cancer = state.cancerTypes.find((c) => c.cancer_type === ct);
    $('cs-title').textContent = rows.length + (rows.length === 1 ? ' gene in ' : ' genes in ') + cancerName(ct);
    App.setNote('cs-missing', notes);
    App.show('cancer-results', true);
    const genes = rows.map((r) => r.gene);
    layoutCards(genes.length > COLLAPSE_ABOVE);
    setLink('cs-tcna', App.tcnaGeneUrl(genes, [ct], 'cancer-specific'));
    const base = ct + '_genes';
    Plots.groupedBars('plot-cs-mean', genes, compartmentSeries(rows, 'mean', true), { ytitle: 'mean (CP10K)', filename: base + '_mean_sd' });
    Plots.groupedBars('plot-cs-cv', genes, compartmentSeries(rows, 'cv'), { ytitle: 'CV', filename: base + '_cv' });
    Plots.groupedBars('plot-cs-mad', genes, compartmentSeries(rows, 'mad'), { ytitle: 'MAD', filename: base + '_mad' });
    Plots.groupedBars('plot-cs-fano', genes, compartmentSeries(rows, 'fano'), { ytitle: 'Fano', filename: base + '_fano' });
    Plots.groupedBars('plot-cs-cv2resid', genes, compartmentSeries(rows, 'cv2_resid'), { ytitle: 'residual log10 CV²', rangemode: 'normal', filename: base + '_cv2_resid' });
    const diff = pick(rows, 'log2fc_cv_median');
    Plots.bars('plot-cs-diff', genes, diff, {
      colors: sigColors(rows, 'q_value'),
      legend: SIG_LEGEND,
      annotations: notTestedNotes(genes, diff),
      hover: rows.map((r) => r.gene + '<br>log2(CV malignant / non-malignant), median: ' + fmt(r.log2fc_cv_median) + '<br>paired patients: ' + fmt(r.n_pairs) + '<br>p = ' + fmtP(r.p_value) + ', q = ' + fmtP(r.q_value)),
      ytitle: 'log2(CV malignant / CV non-malignant)',
      filename: base + '_differential_noise',
    });
    // cancer name and study once above the table instead of a column each
    $('cs-context').textContent = (cancer ? App.cancerLabel(cancer) : ct) + (cancer && cancer.study ? ' · ' + cancer.study : '') + ' · dataset ' + (rows[0].dataset || (cancer && cancer.dataset) || '–');
    App.renderMetricTable('cs-table', {
      corner: 'Metric',
      columns: genes.map((g) => ({ label: g })),
      rows: [{ label: 'Bulk TCNA', link: true, values: genes.map((g) => tcnaLink([g], [ct], 'cancer-specific')) }].concat(metricRows(rows)),
    });
  }

  /* Check the typed or uploaded list against the selected cancer type and list unknown symbols inline. */
  async function checkGeneList(intro) {
    const genes = parseGenes($('genes-text').value);
    const ct = $('cancer-select').value;
    const seq = ++state.lookupSeq;
    if (!genes.length) {
      App.setNote('genes-note', intro ? [intro] : []);
      return;
    }
    const res = await App.lookupGenes(genes, ct);
    if (seq !== state.lookupSeq) return;
    const lines = intro ? [intro] : [];
    if (genes.length > MAX_GENES) lines.push('The list has ' + genes.length + ' symbols; at most ' + MAX_GENES + ' can be analysed at once.');
    if (res) {
      if (genes.length > App.MAX_LOOKUP) lines.push('Only the first ' + App.MAX_LOOKUP + ' symbols were checked.');
      lines.push.apply(lines, geneNotes(cancerName(ct), { unknown: res.unknown, notIn: res.not_in_cancer_type }));
      if (intro && !res.unknown.length && !res.not_in_cancer_type.length) lines.push('All symbols are measured in ' + cancerName(ct) + '.');
    }
    App.setNote('genes-note', lines);
  }

  function renderSampleCounts(cancers) {
    const cats = cancers.map((c) => c.cancer_type);
    const names = cancers.map((c) => c.cancer_name || c.cancer_type);
    Plots.groupedBars(
      'plot-cells',
      cats,
      [
        { key: 'malignant', values: pick(cancers, 'n_cells_malignant'), hover: cancers.map((c, i) => names[i] + '<br>Malignant cells: ' + fmt(c.n_cells_malignant)) },
        { key: 'non_malignant', values: pick(cancers, 'n_cells_non_malignant'), hover: cancers.map((c, i) => names[i] + '<br>Non-malignant cells: ' + fmt(c.n_cells_non_malignant)) },
      ],
      { ytitle: 'cells after QC', filename: 'sample_counts_cells' }
    );
    Plots.groupedBars(
      'plot-patients',
      cats,
      [
        { name: 'Patients', values: pick(cancers, 'n_patients'), color: Plots.COLORS.accent, hover: cancers.map((c, i) => names[i] + '<br>patients: ' + fmt(c.n_patients)) },
        { name: 'Patients with both cell groups', values: pick(cancers, 'n_patients_paired'), color: '#8fb8e8', hover: cancers.map((c, i) => names[i] + '<br>paired patients: ' + fmt(c.n_patients_paired)) },
      ],
      { ytitle: 'patients', filename: 'sample_counts_patients' }
    );
  }

  async function init() {
    document.querySelectorAll('input[name="mode"]').forEach((r) => r.addEventListener('change', () => setMode(r.value)));

    App.attachGeneAutocomplete($('gene-input'), (gene) => runPan(gene));
    $('pan-form').addEventListener('submit', (e) => {
      e.preventDefault();
      runPan($('gene-input').value);
    });
    $('pan-all').addEventListener('click', () => {
      state.cancerTypes.forEach((c) => state.panSelected.add(c.cancer_type));
      onSelectionChange();
    });
    $('pan-none').addEventListener('click', () => {
      state.panSelected.clear();
      onSelectionChange();
    });

    App.attachGeneAutocomplete(
      $('gene-add'),
      (gene) => {
        const ta = $('genes-text');
        const present = parseGenes(ta.value).map((g) => g.toUpperCase());
        if (!present.includes(gene.toUpperCase())) ta.value = (ta.value.trim() ? ta.value.trim() + '\n' : '') + gene;
        $('gene-add').value = '';
      },
      { cancerType: () => $('cancer-select').value }
    );
    App.attachGeneUpload($('genes-file'), $('genes-text'), (genes, filename, added) => {
      checkGeneList(filename + ': ' + added + (added === 1 ? ' new symbol' : ' new symbols') + ' added; the list now has ' + genes.length + '.');
    });
    $('genes-text').addEventListener('change', () => checkGeneList());
    $('cancer-select').addEventListener('change', () => checkGeneList());
    $('cancer-form').addEventListener('submit', (e) => {
      e.preventDefault();
      runCancer($('cancer-select').value, parseGenes($('genes-text').value));
    });
    document.querySelectorAll('#cs-grid details.metric-card').forEach((d) => {
      d.addEventListener('toggle', () => {
        if (d.open) Plots.resizeWithin(d);
      });
      d.querySelector('summary').addEventListener('click', (e) => {
        if (!$('cs-grid').classList.contains('collapsible')) e.preventDefault(); // all cards stay open in the grid
      });
    });

    $('detail-cancer').addEventListener('change', (e) => selectDetailCancer(e.target.value));
    $('detail-metric').addEventListener('change', () => {
      renderPaired();
      renderCelltypes();
    });
    $('detail-ds').addEventListener('change', renderPaired);
    $('pan-download').addEventListener('click', () => App.downloadCsv(state.gene + '_gene_summary.csv', SUMMARY_KEYS, panRows()));
    $('cs-download').addEventListener('click', () => App.downloadCsv(state.csCancer + '_gene_summary.csv', ['gene'].concat(SUMMARY_KEYS.filter((k) => k !== 'gene')), state.csRows));

    const params = App.getParams();
    try {
      state.cancerTypes = await App.loadCancerTypes();
    } catch (e) {
      return;
    }
    App.fillCancerSelect('cancer-select', state.cancerTypes, params.get('cancer_type'));
    buildPicker(params.get('cancer_types'));
    renderSampleCounts(state.cancerTypes);

    if (params.get('mode') === 'cancer') {
      setMode('cancer');
      if (params.get('genes')) {
        $('genes-text').value = parseGenes(params.get('genes')).join('\n');
        runCancer($('cancer-select').value, parseGenes(params.get('genes')));
      }
    } else if (params.get('gene')) {
      $('gene-input').value = params.get('gene');
      runPan(params.get('gene'));
    }
  }

  document.addEventListener('DOMContentLoaded', init);
})();
