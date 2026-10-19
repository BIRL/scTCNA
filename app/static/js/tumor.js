/* Tumor Heterogeneity page: DEPTH / DEPTH2 per cell or per pseudo-bulk, plus the cross-cancer summary. */
(function () {
  'use strict';
  const $ = App.$;
  const fmt = App.fmt;
  const fmtP = App.fmtP;

  const SCORE_LABEL = { depth2: 'DEPTH2', depth: 'DEPTH' };
  const SAMPLE_KEYS = ['patient', 'compartment', 'n_cells', 'depth2_pseudobulk', 'depth_pseudobulk', 'depth2_cell_median', 'depth_cell_median'];
  const SUMMARY_KEYS = [
    'cancer_type', 'cancer_name', 'n_patients_paired', 'depth2_cell_median_malignant', 'depth2_cell_median_non_malignant',
    'depth2_cell_median_diff', 'depth2_cell_p', 'depth2_pseudobulk_median_diff', 'depth2_pseudobulk_p',
  ];
  const SUMMARY_LABELS = {
    depth2_cell_median_malignant: { label: 'DEPTH2 per cell, Malignant' },
    depth2_cell_median_non_malignant: { label: 'DEPTH2 per cell, Non-malignant' },
    depth2_cell_median_diff: { label: 'DEPTH2 per cell, median diff.' },
    depth2_cell_p: { label: 'p (DEPTH2 per cell)' },
    depth2_pseudobulk_median_diff: { label: 'DEPTH2 pseudo-bulk, median diff.' },
    depth2_pseudobulk_p: { label: 'p (DEPTH2 pseudo-bulk)' },
    depth_cell_median_diff: { label: 'DEPTH per cell, median diff.' },
    depth_cell_p: { label: 'p (DEPTH per cell)' },
  };

  const state = { cancerTypes: [] };

  function query() {
    return {
      cancer_type: $('cancer').value,
      score: document.querySelector('input[name="score"]:checked').value,
      level: document.querySelector('input[name="level"]:checked').value,
    };
  }

  function updateLinks() {
    const ct = $('cancer').value;
    $('csv-samples').href = App.apiUrl('/download/sample_scores.csv', { cancer_type: ct });
    $('csv-cells').href = App.apiUrl('/download/cell_scores.csv', { cancer_type: ct });
  }

  function statsText(test, what) {
    if (!test || !test.n_pairs) return 'No paired test: DEPTH is only defined for malignant cells (the non-malignant cells are its reference), or too few patients have both compartments.';
    return 'Paired Wilcoxon on ' + what + ': n = ' + test.n_pairs + ' patient pairs · median difference (Malignant − Non-malignant) = ' + fmt(test.median_diff) + ' · p = ' + fmtP(test.p_value) + ' · rank-biserial = ' + fmt(test.rank_biserial) + (App.isNum(test.p_value) ? '' : ' (test needs at least 6 pairs)');
  }

  async function run() {
    App.clearError();
    const q = query();
    updateLinks();
    App.setParams(q);
    let res;
    try {
      res = await App.apiGet('/cancer-types/' + encodeURIComponent(q.cancer_type) + '/tumor-heterogeneity', { score: q.score, level: q.level });
    } catch (e) {
      return;
    }
    render(res);
  }

  function render(res) {
    const label = SCORE_LABEL[res.score] || res.score;
    const cancer = state.cancerTypes.find((c) => c.cancer_type === res.cancer_type);
    const name = cancer ? App.cancerLabel(cancer) : res.cancer_type;
    App.show('results', true);
    const note = $('note');
    note.hidden = true;
    const filename = 'tumor_heterogeneity_' + res.cancer_type + '_' + res.score + '_' + res.level;

    if (res.level === 'cell') {
      const comp = res.compartments;
      $('title').textContent = label + ' per cell – ' + name;
      $('subtitle').textContent = 'Malignant: ' + App.fmtInt(comp.malignant.n) + ' cells, median ' + fmt(comp.malignant.median) + ' · Non-malignant: ' + App.fmtInt(comp.non_malignant.n) + ' cells, median ' + fmt(comp.non_malignant.median) + '. Violins show up to 5,000 cells per compartment; dots are the per-patient medians.';
      Plots.violins('plot-main', comp, res.patients, { ytitle: label + ' per cell (log2 CP10K + 1)', filename: filename });
      $('stats').textContent = statsText(res.paired_test, 'the per-patient medians');
      $('table-title').textContent = 'Per-patient medians of the per-cell ' + label;
      App.renderTable('table', App.columns(['patient', 'compartment', 'n_cells', 'median'], { median: { label: 'Median ' + label + ' per cell' } }), res.patients, { sortKey: 'patient', sortDir: 'asc' });
      return;
    }

    const key = res.score + '_pseudobulk';
    $('title').textContent = label + ' per pseudo-bulk profile – ' + name;
    $('subtitle').textContent = 'One profile per patient and compartment: summed counts, log2(CPM + 1).';
    const byPatient = {};
    res.rows.forEach((r) => {
      byPatient[r.patient] = byPatient[r.patient] || {};
      byPatient[r.patient][r.compartment] = r;
    });
    const patients = Object.keys(byPatient).sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
    const value = (p, comp) => (byPatient[p][comp] && App.isNum(byPatient[p][comp][key]) ? byPatient[p][comp][key] : null);
    const pairs = patients.map((p) => ({ label: p, a: value(p, 'malignant'), b: value(p, 'non_malignant') }));
    const hasNonMalignant = pairs.some((p) => App.isNum(p.b));
    if (hasNonMalignant) {
      Plots.pairedLines('plot-main', pairs, { ytitle: label + ' (pseudo-bulk)', filename: filename });
    } else {
      // DEPTH scores malignant profiles against the non-malignant mean, so only malignant values exist
      note.hidden = false;
      note.textContent = 'DEPTH scores each malignant pseudo-bulk against the mean of the non-malignant pseudo-bulks, so non-malignant profiles have no DEPTH value and no paired test is possible at this level. Use the cell level for a paired comparison of DEPTH.';
      const shown = pairs.filter((p) => App.isNum(p.a));
      Plots.bars('plot-main', shown.map((p) => p.label), shown.map((p) => p.a), { colors: Plots.COLORS.malignant, ytitle: label + ' (malignant pseudo-bulk)', xtitle: 'patient', filename: filename });
    }
    $('stats').textContent = hasNonMalignant ? statsText(res.paired_test, 'the pseudo-bulk scores') : '';
    $('table-title').textContent = 'Per-patient scores';
    App.renderTable('table', App.columns(SAMPLE_KEYS), res.rows, { sortKey: 'patient', sortDir: 'asc' });
  }

  async function loadSummary() {
    let res;
    try {
      res = await App.apiGet('/tumor-heterogeneity/summary');
    } catch (e) {
      return;
    }
    const rows = res.rows || [];
    const cats = rows.map((r) => r.cancer_type);
    Plots.groupedBars(
      'plot-summary',
      cats,
      [
        { key: 'malignant', values: rows.map((r) => (App.isNum(r.depth2_cell_median_malignant) ? r.depth2_cell_median_malignant : null)), hover: rows.map((r) => (r.cancer_name || r.cancer_type) + '<br>Malignant DEPTH2 per cell: ' + fmt(r.depth2_cell_median_malignant) + '<br>paired p = ' + fmtP(r.depth2_cell_p)) },
        { key: 'non_malignant', values: rows.map((r) => (App.isNum(r.depth2_cell_median_non_malignant) ? r.depth2_cell_median_non_malignant : null)), hover: rows.map((r) => (r.cancer_name || r.cancer_type) + '<br>Non-malignant DEPTH2 per cell: ' + fmt(r.depth2_cell_median_non_malignant)) },
      ],
      { ytitle: 'DEPTH2 per cell (median of patient medians)', filename: 'tumor_heterogeneity_summary' }
    );
    App.renderTable('summary-table', App.columns(SUMMARY_KEYS, SUMMARY_LABELS), rows, { sortKey: 'cancer_type', sortDir: 'asc' });
  }

  async function init() {
    const params = App.getParams();
    if (params.get('score') === 'depth') $('score-depth').checked = true;
    if (params.get('level') === 'sample') $('level-sample').checked = true;
    $('form').addEventListener('submit', (e) => {
      e.preventDefault();
      run();
    });
    $('form').addEventListener('change', run);
    loadSummary();
    try {
      state.cancerTypes = await App.loadCancerTypes();
    } catch (e) {
      return;
    }
    App.fillCancerSelect('cancer', state.cancerTypes, params.get('cancer_type'));
    if (state.cancerTypes.length) run();
  }

  document.addEventListener('DOMContentLoaded', init);
})();
