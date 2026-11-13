/* Cross-cancer summary page: patient-level noise contrasts (GET /api/noise-contrasts). */
(function () {
  'use strict';
  const $ = App.$;
  const fmt = App.fmt;

  const ALL = 'non_malignant';
  const MIN_PATIENTS_LINEAGE = 3;
  const MEASURES = [
    { key: 'log2fc_cv', summary: 'median_log2fc_cv', p: 'p_log2fc_cv', label: 'raw log2(CV malignant / CV comparison group)', short: 'log2 CV ratio (raw)', color: '#9aa7b4' },
    { key: 'log2fc_cv_ds', summary: 'median_log2fc_cv_ds', p: 'p_log2fc_cv_ds', label: 'depth-matched log2(CV malignant / CV comparison group)', short: 'log2 CV ratio (depth-matched)', color: '#5b8fd1' },
    { key: 'excess_cv2_ds', summary: 'median_excess_cv2_ds', p: 'p_excess_cv2_ds', label: 'excess noise, log10 CV² (depth-matched, mean-corrected)', short: 'excess noise (log10 CV²)', color: '#e4572e' },
  ];
  const sig = (v) => (App.isNum(v) && v < 0.05 ? 'sig' : '');
  const pct = (v) => (App.isNum(v) ? Math.round(v * 100) + '%' : '–');
  // Rows of the transposed tables: measures in rows, cancer types or comparison groups in columns.
  const TABLE_ROWS = [
    { key: 'n_patients', label: 'Patients' },
    { key: 'median_log2fc_cv', label: 'log2 CV ratio, raw (median over patients)' },
    { key: 'q_log2fc_cv', label: 'q (raw)', format: App.fmtP, className: sig },
    { key: 'median_log2fc_cv_ds', label: 'log2 CV ratio, depth-matched (median over patients)' },
    { key: 'q_log2fc_cv_ds', label: 'q (depth-matched)', format: App.fmtP, className: sig },
    { key: 'median_excess_cv2_ds', label: 'Excess noise, log10 CV² (median over patients)' },
    { key: 'p_excess_cv2_ds', label: 'p (excess noise)', format: App.fmtP },
    { key: 'q_excess_cv2_ds', label: 'q (excess noise)', format: App.fmtP, className: sig },
    { key: 'frac_patients_excess_positive', label: 'Patients with excess noise > 0', format: pct },
  ];
  const state = { rows: [], cancers: [] };

  function refLabel(ref) {
    return ref === ALL ? 'All non-malignant cells' : ref.replace(/_/g, ' ');
  }

  function tableRows(rows) {
    return TABLE_ROWS.map((t) => ({ label: t.label, title: t.key, values: rows.map((r) => r[t.key]), format: t.format, className: t.className }));
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
        hover: rows.map((r) => r.cancer_type + ' · ' + m.label + '<br>median over patients: ' + fmt(r[m.summary]) + '<br>p: ' + App.fmtP(r[m.p]) + ' · q: ' + App.fmtP(r['q' + m.p.slice(1)]) + ' · patients: ' + r.n_patients),
      })),
      { ytitle: 'median over patients', rangemode: 'normal', filename: 'sctcna_noise_contrast_all' }
    );
    App.renderMetricTable('table-all', {
      corner: 'Measure',
      columns: rows.map((r) => {
        const c = state.cancers.find((x) => x.cancer_type === r.cancer_type);
        return { label: r.cancer_type, title: c && c.cancer_name ? c.cancer_name : r.cancer_type };
      }),
      rows: tableRows(rows),
    });
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
    // a later selection may have finished first; only the current one is drawn
    if (ct !== $('cancer').value || measure.key !== $('measure').value) return;
    const refs = [ALL, ...[...new Set(data.rows.map((r) => r.reference))].filter((r) => r !== ALL).sort()];
    const cancer = state.cancers.find((c) => c.cancer_type === ct);
    const traces = [
      {
        type: 'box',
        name: measure.short,
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
        showlegend: true,
      },
    ];
    Plots.draw(
      'plot-patients',
      traces,
      {
        title: { text: (cancer ? App.cancerLabel(cancer) : ct) + ': ' + measure.short + ' per patient', x: 0, xanchor: 'left', y: 0.98, yanchor: 'top', yref: 'container', font: { size: 14 } },
        margin: { l: 60, r: 16, t: 76, b: 56 },
        xaxis: { title: 'comparison group', type: 'category', categoryorder: 'array', categoryarray: refs.map(refLabel), tickangle: refs.length > 6 ? -30 : 0 },
        yaxis: { title: measure.short, zeroline: true, zerolinecolor: '#55606b' },
        showlegend: true,
      },
      'sctcna_noise_contrast_' + ct
    );
    // comparison groups as columns, most patients first
    const summary = data.summary;
    const order = summary.map((r, i) => i).sort((a, b) => App.isNum(summary[b].n_patients) - App.isNum(summary[a].n_patients) || (summary[b].n_patients || 0) - (summary[a].n_patients || 0));
    App.renderMetricTable('table-patients', {
      corner: 'Measure',
      columns: summary.map((r) => ({ label: refLabel(r.reference), title: r.reference })),
      rows: tableRows(summary),
      order: order,
      emptyText: 'No comparison groups for this cancer type.',
    });
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
