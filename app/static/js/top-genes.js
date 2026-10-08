/* Top genes page: ranked table and bar chart of one metric within one cancer type. */
(function () {
  'use strict';
  const $ = App.$;
  const fmt = App.fmt;
  const fmtP = App.fmtP;

  const MAX_BARS = 100;
  const TABLE_KEYS = [
    'rank', 'gene', 'value', 'value_q', 'mean_malignant', 'mean_non_malignant', 'detection_rate_malignant', 'detection_rate_non_malignant',
    'cv_malignant', 'cv_non_malignant', 'mad_malignant', 'mad_non_malignant', 'fano_malignant', 'fano_non_malignant',
    'cv2_resid_malignant', 'cv2_resid_non_malignant', 'log2fc_cv_median', 'delta_cv2_resid', 'p_value', 'q_value', 'rank_biserial', 'n_pairs',
  ];

  const state = { cancerTypes: [] };

  function metricInfo(key) {
    return App.RANK_METRICS.find((m) => m.key === key) || { key: key, label: key };
  }

  function query() {
    return {
      cancer_type: $('cancer').value,
      metric: $('metric').value,
      direction: document.querySelector('input[name="direction"]:checked').value,
      n: Math.min(500, Math.max(1, parseInt($('n').value, 10) || 50)),
      expressed_only: $('expressed').checked,
    };
  }

  function updateCsvLink() {
    $('csv-link').href = App.apiUrl('/download/top_genes.csv', query());
  }

  async function run() {
    App.clearError();
    const q = query();
    $('n').value = q.n;
    updateCsvLink();
    App.setParams(q);
    let res;
    try {
      res = await App.apiGet('/cancer-types/' + encodeURIComponent(q.cancer_type) + '/top-genes', { metric: q.metric, direction: q.direction, n: q.n, expressed_only: q.expressed_only });
    } catch (e) {
      return;
    }
    render(res, q);
  }

  function render(res, q) {
    const info = metricInfo(res.metric);
    const rows = res.rows.map((r, i) => Object.assign({ rank: i + 1 }, r));
    const cancer = state.cancerTypes.find((c) => c.cancer_type === res.cancer_type);
    $('title').textContent = (res.direction === 'noisy' ? 'Noisiest ' : 'Most stable ') + rows.length + ' genes in ' + (cancer ? App.cancerLabel(cancer) : res.cancer_type);
    $('subtitle').textContent = info.label + (q.expressed_only ? ' · expressed genes only' : ' · all genes');
    App.show('results', true);

    const shown = rows.slice(0, MAX_BARS);
    let colors;
    let chips = '';
    if (info.hasQ) {
      colors = shown.map((r) => (App.isNum(r.value_q) && r.value_q < 0.05 ? Plots.COLORS.significant : Plots.COLORS.neutral));
      chips = '<span class="chip chip-sig">q &lt; 0.05</span><span class="chip chip-ns">not significant</span>';
    } else {
      colors = Plots.COLORS[info.compartment] || Plots.COLORS.accent;
      chips = info.compartment === 'malignant' ? '<span class="chip chip-mal">Malignant</span>' : '<span class="chip chip-non">Non-malignant</span>';
    }
    $('chips').innerHTML = chips;
    Plots.bars('plot-top', shown.map((r) => r.gene), shown.map((r) => r.value), {
      horizontal: true,
      colors: colors,
      hover: shown.map(
        (r) =>
          r.gene + '<br>' + info.label + ': ' + fmt(r.value) +
          '<br>mean Mal / Non-mal: ' + fmt(r.mean_malignant) + ' / ' + fmt(r.mean_non_malignant) +
          '<br>CV Mal / Non-mal: ' + fmt(r.cv_malignant) + ' / ' + fmt(r.cv_non_malignant) +
          (info.hasQ ? '<br>q (this metric) = ' + fmtP(r.value_q) + ', paired patients = ' + fmt(r.n_pairs) : '')
      ),
      xtitle: info.label,
      height: Math.max(380, 18 * shown.length + 90),
      filename: 'top_genes_' + res.cancer_type + '_' + res.metric + '_' + res.direction,
    });
    $('plot-note').textContent = rows.length > MAX_BARS ? 'The chart shows the first ' + MAX_BARS + ' genes; the table and the CSV hold all ' + rows.length + '.' : '';
    App.renderTable('table', App.columns(TABLE_KEYS, { value: { label: info.label } }), rows);
  }

  async function init() {
    const params = App.getParams();
    App.fillSelect('metric', App.RANK_METRICS, params.get('metric') || 'log2fc_cv_median');
    if (params.get('direction') === 'stable') $('dir-stable').checked = true;
    if (params.get('n')) $('n').value = params.get('n');
    if (params.get('expressed_only') === 'false') $('expressed').checked = false;
    $('form').addEventListener('submit', (e) => {
      e.preventDefault();
      run();
    });
    $('form').addEventListener('change', updateCsvLink);
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
