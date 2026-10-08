/* Gene Noise Analysis page: pan-cancer view of one gene, cancer-specific view of a gene
 * list, per-patient / per-cell-type detail and the sample-count overview. */
(function () {
  'use strict';
  const $ = App.$;
  const fmt = App.fmt;
  const fmtP = App.fmtP;

  const MAX_GENES = 50;

  // gene_summary columns in display order (cancer first, then the metric blocks)
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

  const state = {
    mode: 'pan',
    cancerTypes: [],
    gene: null,
    rows: [], // gene_summary rows of the current gene, one per cancer type
    detailCancer: null,
    patients: null,
    celltypes: null,
    csRows: [],
    csCancer: null,
  };

  // ------------------------------------------------------------------ helpers
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

  function sigColors(rows, qKey) {
    return rows.map((r) => (App.isNum(r[qKey]) && r[qKey] < 0.05 ? Plots.COLORS.significant : Plots.COLORS.neutral));
  }

  function diffHover(rows) {
    return rows.map(
      (r) =>
        (r.cancer_name || r.cancer_type) +
        '<br>log2(CV mal/non-mal) median: ' + fmt(r.log2fc_cv_median) +
        '<br>paired patients: ' + fmt(r.n_pairs) +
        '<br>p = ' + fmtP(r.p_value) + ', q = ' + fmtP(r.q_value) +
        '<br>rank-biserial: ' + fmt(r.rank_biserial)
    );
  }

  // ------------------------------------------------------------------ mode switch
  function setMode(mode) {
    state.mode = mode;
    $('mode-pan').checked = mode === 'pan';
    $('mode-cancer').checked = mode === 'cancer';
    App.show('pan-form', mode === 'pan');
    App.show('cancer-form', mode === 'cancer');
    App.show('pan-results', mode === 'pan' && state.rows.length > 0);
    App.show('cancer-results', mode === 'cancer' && state.csRows.length > 0);
  }

  // ------------------------------------------------------------------ pan-cancer
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
    App.setParams({ mode: 'pan', gene: res.gene });
    renderPan();
  }

  function renderPan() {
    const rows = state.rows;
    const cats = rows.map((r) => r.cancer_type);
    $('pan-title').textContent = state.gene + ' across ' + rows.length + (rows.length === 1 ? ' cancer type' : ' cancer types');
    App.show('pan-results', true);
    const base = state.gene;
    Plots.groupedBars('plot-mean', cats, compartmentSeries(rows, 'mean', true), { ytitle: 'mean (CP10K)', filename: base + '_mean_sd' });
    Plots.groupedBars('plot-cv', cats, compartmentSeries(rows, 'cv'), { ytitle: 'CV', filename: base + '_cv' });
    Plots.groupedBars('plot-mad', cats, compartmentSeries(rows, 'mad'), { ytitle: 'MAD', filename: base + '_mad' });
    Plots.groupedBars('plot-fano', cats, compartmentSeries(rows, 'fano'), { ytitle: 'Fano', filename: base + '_fano' });
    Plots.groupedBars('plot-cv2resid', cats, compartmentSeries(rows, 'cv2_resid'), { ytitle: 'residual log10 CV²', rangemode: 'normal', filename: base + '_cv2_resid' });
    Plots.bars('plot-diff', cats, pick(rows, 'log2fc_cv_median'), {
      colors: sigColors(rows, 'q_value'),
      text: rows.map((r) => (App.isNum(r.n_pairs) ? 'n=' + r.n_pairs : '')),
      hover: diffHover(rows),
      ytitle: 'log2(CV mal / CV non-mal)',
      filename: base + '_differential_noise',
    }).then(() => {
      const div = $('plot-diff');
      if (!div.dataset.clickBound) {
        div.dataset.clickBound = '1';
        div.on('plotly_click', (ev) => {
          const pt = ev.points && ev.points[0];
          if (pt && pt.x) selectDetailCancer(String(pt.x));
        });
      }
    });
    App.renderTable('pan-table', App.columns(SUMMARY_KEYS), rows, { sortKey: 'cancer_type', sortDir: 'asc' });

    // per-patient detail: keep the chosen cancer when the new gene has it, else the first
    const available = rows.map((r) => r.cancer_type);
    const keep = available.includes(state.detailCancer) ? state.detailCancer : available[0];
    App.fillCancerSelect('detail-cancer', rows.map((r) => ({ cancer_type: r.cancer_type, cancer_name: r.cancer_name })), keep);
    selectDetailCancer(keep);
  }

  function selectDetailCancer(cancerType) {
    if (!cancerType) return;
    state.detailCancer = cancerType;
    $('detail-cancer').value = cancerType;
    if (state.mode === 'pan') $('detail-card').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
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
    if (!rows) {
      Plots.empty('plot-paired', state.patientsError || 'No per-patient values.');
      $('paired-stats').textContent = '';
      App.renderTable('patients-table', [], [], { emptyText: state.patientsError || 'No per-patient values (per-patient tables only hold genes expressed in both compartments).' });
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
      text += ' · median log2(Mal/Non-mal) = ' + fmt(summary[stats[0]]) + ' · paired Wilcoxon p = ' + fmtP(summary[stats[1]]) + ', q = ' + fmtP(summary[stats[2]]);
      if (metric === 'cv') text += ' · rank-biserial = ' + fmt(summary.rank_biserial);
    } else if (metric === 'cv2_resid') {
      text += ' · median difference Mal − Non-mal = ' + fmt(summary.delta_cv2_resid);
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
      hover: sorted.map((r) => r.cell_type + '<br>' + METRIC_LABELS[metric] + ': ' + fmt(r[metric]) + '<br>patients: ' + fmt(r.n_patients) + '<br>mean: ' + fmt(r.mean)),
      ytitle: METRIC_LABELS[metric],
      filename: state.gene + '_' + state.detailCancer + '_celltypes_' + metric,
    });
    App.renderTable('celltypes-table', App.columns(CELLTYPE_KEYS), sorted);
  }

  // ------------------------------------------------------------------ cancer-specific
  async function runCancer(cancerType, genes) {
    App.clearError();
    if (!genes.length) {
      App.showError('Enter at least one gene symbol.');
      return;
    }
    if (genes.length > MAX_GENES) {
      App.showError('At most ' + MAX_GENES + ' genes per query (got ' + genes.length + ').');
      return;
    }
    const results = await Promise.allSettled(genes.map((g) => App.apiGet('/genes/' + encodeURIComponent(g), null, { silent: true })));
    const rows = [];
    const missing = [];
    results.forEach((r, i) => {
      if (r.status !== 'fulfilled') {
        missing.push(genes[i]);
        return;
      }
      const row = (r.value.cancers || []).find((c) => c.cancer_type === cancerType);
      if (row) rows.push(row);
      else missing.push(r.value.gene + ' (not in ' + cancerType + ')');
    });
    if (!rows.length) {
      App.showError('None of the requested genes is present in ' + cancerType + '.');
      return;
    }
    state.csRows = rows;
    state.csCancer = cancerType;
    App.setParams({ mode: 'cancer', cancer_type: cancerType, genes: rows.map((r) => r.gene).join(',') });
    renderCancer(missing);
  }

  function renderCancer(missing) {
    const rows = state.csRows;
    const cancer = state.cancerTypes.find((c) => c.cancer_type === state.csCancer);
    $('cs-title').textContent = rows.length + (rows.length === 1 ? ' gene in ' : ' genes in ') + (cancer ? App.cancerLabel(cancer) : state.csCancer);
    const note = $('cs-missing');
    note.hidden = missing.length === 0;
    note.textContent = missing.length ? 'Not shown: ' + missing.join(', ') : '';
    App.show('cancer-results', true);
    const genes = rows.map((r) => r.gene);
    const base = state.csCancer + '_genes';
    Plots.groupedBars('plot-cs-mean', genes, compartmentSeries(rows, 'mean', true), { ytitle: 'mean (CP10K)', filename: base + '_mean_sd' });
    Plots.groupedBars('plot-cs-cv', genes, compartmentSeries(rows, 'cv'), { ytitle: 'CV', filename: base + '_cv' });
    Plots.groupedBars('plot-cs-mad', genes, compartmentSeries(rows, 'mad'), { ytitle: 'MAD', filename: base + '_mad' });
    Plots.groupedBars('plot-cs-fano', genes, compartmentSeries(rows, 'fano'), { ytitle: 'Fano', filename: base + '_fano' });
    Plots.groupedBars('plot-cs-cv2resid', genes, compartmentSeries(rows, 'cv2_resid'), { ytitle: 'residual log10 CV²', rangemode: 'normal', filename: base + '_cv2_resid' });
    Plots.bars('plot-cs-diff', genes, pick(rows, 'log2fc_cv_median'), {
      colors: sigColors(rows, 'q_value'),
      hover: rows.map((r) => r.gene + '<br>log2(CV mal/non-mal) median: ' + fmt(r.log2fc_cv_median) + '<br>paired patients: ' + fmt(r.n_pairs) + '<br>p = ' + fmtP(r.p_value) + ', q = ' + fmtP(r.q_value)),
      ytitle: 'log2(CV mal / CV non-mal)',
      filename: base + '_differential_noise',
    });
    App.renderTable('cs-table', App.columns(['gene'].concat(SUMMARY_KEYS.filter((k) => k !== 'gene'))), rows, { sortKey: 'gene', sortDir: 'asc' });
  }

  // ------------------------------------------------------------------ sample counts
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
        { name: 'Patients with both compartments', values: pick(cancers, 'n_patients_paired'), color: '#8fb8e8', hover: cancers.map((c, i) => names[i] + '<br>paired patients: ' + fmt(c.n_patients_paired)) },
      ],
      { ytitle: 'patients', filename: 'sample_counts_patients' }
    );
  }

  // ------------------------------------------------------------------ wiring
  async function init() {
    document.querySelectorAll('input[name="mode"]').forEach((r) => r.addEventListener('change', () => setMode(r.value)));

    App.attachGeneAutocomplete($('gene-input'), (gene) => runPan(gene));
    $('pan-form').addEventListener('submit', (e) => {
      e.preventDefault();
      runPan($('gene-input').value);
    });

    App.attachGeneAutocomplete($('gene-add'), (gene) => {
      const ta = $('genes-text');
      const present = parseGenes(ta.value).map((g) => g.toUpperCase());
      if (!present.includes(gene.toUpperCase())) ta.value = (ta.value.trim() ? ta.value.trim() + '\n' : '') + gene;
      $('gene-add').value = '';
    });
    $('cancer-form').addEventListener('submit', (e) => {
      e.preventDefault();
      runCancer($('cancer-select').value, parseGenes($('genes-text').value));
    });

    $('detail-cancer').addEventListener('change', (e) => selectDetailCancer(e.target.value));
    $('detail-metric').addEventListener('change', () => {
      renderPaired();
      renderCelltypes();
    });
    $('detail-ds').addEventListener('change', renderPaired);
    $('pan-download').addEventListener('click', () => App.downloadCsv(state.gene + '_gene_summary.csv', SUMMARY_KEYS, state.rows));
    $('cs-download').addEventListener('click', () => App.downloadCsv(state.csCancer + '_gene_summary.csv', ['gene'].concat(SUMMARY_KEYS.filter((k) => k !== 'gene')), state.csRows));

    const params = App.getParams();
    try {
      state.cancerTypes = await App.loadCancerTypes();
    } catch (e) {
      return;
    }
    App.fillCancerSelect('cancer-select', state.cancerTypes, params.get('cancer_type'));
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
