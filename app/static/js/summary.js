/* Cross-cancer summary: patient-level noise contrasts of malignant cells against all
 * non-malignant cells and against single lineages (GET /api/noise-contrasts). */
(function () {
  'use strict';
  const $ = App.$;
  const fmt = App.fmt;

  const ALL = 'non_malignant';
  const MIN_PATIENTS_LINEAGE = 3;
  const MEASURES = [
    { key: 'log2fc_cv', summary: 'median_log2fc_cv', p: 'p_log2fc_cv', label: 'log2(CV mal / CV ref), raw', short: 'log2 CV ratio (raw)', color: '#9aa7b4' },
    { key: 'log2fc_cv_ds', summary: 'median_log2fc_cv_ds', p: 'p_log2fc_cv_ds', label: 'log2(CV mal / CV ref), depth-matched', short: 'log2 CV ratio (depth-matched)', color: '#5b8fd1' },
    { key: 'excess_cv2_ds', summary: 'median_excess_cv2_ds', p: 'p_excess_cv2_ds', label: 'excess noise, log10 CV² (depth-matched, mean-corrected)', short: 'excess noise (log10 CV²)', color: '#e4572e' },
  ];
  const state = { rows: [], cancers: [] };

  function refLabel(ref) {
    return ref === ALL ? 'All non-malignant' : ref.replace(/_/g, ' ');
  }

  function renderAll() {
    const rows = state.rows.filter((r) => r.reference === ALL).sort((a, b) => (a.cancer_type < b.cancer_type ? -1 : 1));
    const cats = rows.map((r) => r.cancer_type);
    Plots.groupedBars(
      'plot-all',
      cats,
      MEASURES.map((m) => ({
        key: m.key,
        name: m.label,
        color: m.color,
        values: rows.map((r) => r[m.summary]),
        hover: rows.map((r) => r.cancer_type + ' · ' + m.label + '<br>median over patients: ' + fmt(r[m.summary]) + '<br>p: ' + App.fmtP(r[m.p]) + ' · patients: ' + r.n_patients),
      })),
      { ytitle: 'median over patients', rangemode: 'normal', filename: 'sctcna_noise_contrast_all' }
    );
    const keys = ['cancer_type', 'n_patients', 'median_log2fc_cv', 'p_log2fc_cv', 'median_log2fc_cv_ds', 'p_log2fc_cv_ds', 'median_excess_cv2_ds', 'p_excess_cv2_ds', 'frac_patients_excess_positive'];
    const labels = {
      median_log2fc_cv: { label: 'log2 CV ratio (raw)' },
      p_log2fc_cv: { label: 'p (raw)', format: App.fmtP },
      median_log2fc_cv_ds: { label: 'log2 CV ratio (depth-matched)' },
      p_log2fc_cv_ds: { label: 'p (depth-matched)', format: App.fmtP },
      median_excess_cv2_ds: { label: 'Excess noise (log10 CV²)' },
      p_excess_cv2_ds: { label: 'p (excess)', format: App.fmtP },
      frac_patients_excess_positive: { label: 'Patients with excess > 0', format: (v) => (App.isNum(v) ? Math.round(v * 100) + '%' : '–') },
    };
    App.renderTable('table-all', App.columns(keys, labels), rows, { sortKey: 'cancer_type', sortDir: 'asc' });
    $('csv-all').href = App.apiUrl('/download/noise_contrasts.csv');
  }

  function renderLineages() {
    const rows = state.rows.filter((r) => r.reference !== ALL && r.n_patients >= MIN_PATIENTS_LINEAGE);
    if (!rows.length) {
      Plots.empty('plot-lineages', 'No lineage has enough cells in enough patients.');
      return;
    }
    const cancers = [...new Set(state.rows.map((r) => r.cancer_type))].sort();
    const lineages = [...new Set(rows.map((r) => r.reference))].sort();
    const values = lineages.map((ref) => cancers.map((ct) => {
      const r = rows.find((x) => x.reference === ref && x.cancer_type === ct);
      return r ? r.median_excess_cv2_ds : null;
    }));
    Plots.heatmap('plot-lineages', lineages.map(refLabel), cancers, values, { diverging: true, title: 'excess noise (log10 CV²)', colorbar: 'log10 CV² excess', filename: 'sctcna_noise_contrast_lineages' });
  }

  async function renderPatients() {
    const ct = $('cancer').value;
    const measure = MEASURES.find((m) => m.key === $('measure').value) || MEASURES[2];
    App.setParams({ cancer_type: ct, measure: measure.key });
    let data;
    try {
      data = await App.apiGet('/cancer-types/' + encodeURIComponent(ct) + '/noise-contrasts');
    } catch (e) {
      return;
    }
    const refs = [ALL, ...[...new Set(data.rows.map((r) => r.reference))].filter((r) => r !== ALL).sort()];
    const traces = [
      {
        type: 'box',
        x: data.rows.map((r) => refLabel(r.reference)),
        y: data.rows.map((r) => r[measure.key]),
        boxpoints: 'all',
        jitter: 0.4,
        pointpos: 0,
        fillcolor: 'rgba(0,0,0,0)',
        line: { color: '#9aa7b4' },
        marker: { color: measure.color, size: 7, opacity: 0.85 },
        text: data.rows.map((r) => r.patient + ' · ' + r.n_genes + ' genes'),
        hovertemplate: '%{text}<br>%{y:.3f}<extra></extra>',
        showlegend: false,
      },
    ];
    Plots.draw(
      'plot-patients',
      traces,
      {
        xaxis: { title: 'reference group', type: 'category', categoryorder: 'array', categoryarray: refs.map(refLabel), tickangle: refs.length > 6 ? -30 : 0 },
        yaxis: { title: measure.short, zeroline: true, zerolinecolor: '#55606b' },
        showlegend: false,
      },
      'sctcna_noise_contrast_' + ct
    );
    const keys = ['reference', 'n_patients', 'median_log2fc_cv', 'p_log2fc_cv', 'median_log2fc_cv_ds', 'p_log2fc_cv_ds', 'median_excess_cv2_ds', 'p_excess_cv2_ds', 'frac_patients_excess_positive'];
    const cols = App.columns(keys, {
      reference: { label: 'Reference', format: refLabel },
      median_log2fc_cv: { label: 'log2 CV ratio (raw)' },
      p_log2fc_cv: { label: 'p (raw)', format: App.fmtP },
      median_log2fc_cv_ds: { label: 'log2 CV ratio (depth-matched)' },
      p_log2fc_cv_ds: { label: 'p (depth-matched)', format: App.fmtP },
      median_excess_cv2_ds: { label: 'Excess noise (log10 CV²)' },
      p_excess_cv2_ds: { label: 'p (excess)', format: App.fmtP },
      frac_patients_excess_positive: { label: 'Patients with excess > 0', format: (v) => (App.isNum(v) ? Math.round(v * 100) + '%' : '–') },
    });
    App.renderTable('table-patients', cols, data.summary, { sortKey: 'n_patients', sortDir: 'desc' });
  }

  async function init() {
    App.fillSelect(
      $('measure'),
      MEASURES.map((m) => ({ key: m.key, label: m.label })),
      'excess_cv2_ds'
    );
    try {
      const [contrasts, cancers] = await Promise.all([App.apiGet('/noise-contrasts'), App.loadCancerTypes()]);
      state.rows = contrasts.rows;
      state.cancers = cancers;
    } catch (e) {
      return;
    }
    const params = App.getParams();
    App.fillCancerSelect($('cancer'), state.cancers, params.get('cancer_type'));
    if (params.get('measure')) $('measure').value = params.get('measure');
    renderAll();
    renderLineages();
    renderPatients();
    $('cancer').addEventListener('change', renderPatients);
    $('measure').addEventListener('change', renderPatients);
  }

  document.addEventListener('DOMContentLoaded', init);
})();
