/* Heatmap page: genes x cancer types for one ranking metric, optional row z-scores. */
(function () {
  'use strict';
  const $ = App.$;

  const MAX_GENES = 200;
  const state = { last: null };

  function metricInfo(key) {
    return App.RANK_METRICS.find((m) => m.key === key) || { key: key, label: key };
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

  function query() {
    return { genes: parseGenes($('genes').value).join(','), metric: $('metric').value, expressed_only: $('expressed').checked };
  }

  async function run() {
    App.clearError();
    const q = query();
    const genes = q.genes ? q.genes.split(',') : [];
    if (!genes.length) {
      App.showError('Enter at least one gene symbol.');
      return;
    }
    if (genes.length > MAX_GENES) {
      App.showError('At most ' + MAX_GENES + ' genes per heatmap (got ' + genes.length + ').');
      return;
    }
    let res;
    try {
      res = await App.apiGet('/heatmap', q);
    } catch (e) {
      return;
    }
    state.last = res;
    App.setParams({ genes: q.genes, metric: q.metric, expressed_only: q.expressed_only || null, zscore: $('zscore').checked || null });
    const link = $('csv-link');
    link.href = App.apiUrl('/download/heatmap.csv', q);
    link.hidden = false;
    render();
  }

  function render() {
    const res = state.last;
    if (!res) return;
    const info = metricInfo(res.metric);
    const zscore = $('zscore').checked;
    $('title').textContent = res.genes.length + (res.genes.length === 1 ? ' gene' : ' genes') + ' × ' + res.cancer_types.length + (res.cancer_types.length === 1 ? ' cancer type' : ' cancer types');
    $('subtitle').textContent = info.label + (zscore ? ' · row z-scores' : '') + ($('expressed').checked ? ' · expressed genes only' : '');
    const missing = $('missing');
    missing.hidden = !res.missing.length;
    missing.textContent = res.missing.length ? 'Not found in any dataset: ' + res.missing.join(', ') : '';
    App.show('results', true);
    Plots.heatmap('plot-heatmap', res.genes, res.cancer_types, res.values, {
      diverging: !!info.diverging,
      zscore: zscore,
      qValues: res.q_values,
      title: info.label,
      colorbar: info.label.length > 28 ? res.metric : info.label,
      filename: 'heatmap_' + res.metric + (zscore ? '_zscore' : ''),
    });
    const rows = res.genes.map((g, i) => {
      const row = { gene: g };
      res.cancer_types.forEach((ct, j) => {
        row[ct] = res.values[i][j];
      });
      return row;
    });
    App.renderTable('table', App.columns(['gene'].concat(res.cancer_types)), rows, { sortKey: null });
  }

  async function fillExample() {
    App.clearError();
    let cancers;
    try {
      cancers = await App.loadCancerTypes();
    } catch (e) {
      return;
    }
    if (!cancers.length) return;
    try {
      const res = await App.apiGet('/cancer-types/' + encodeURIComponent(cancers[0].cancer_type) + '/top-genes', { metric: 'log2fc_cv_median', direction: 'noisy', n: 25 });
      $('genes').value = res.rows.map((r) => r.gene).join('\n');
      $('metric').value = 'log2fc_cv_median';
      run();
    } catch (e) {
      /* the banner already shows the error */
    }
  }

  function init() {
    const params = App.getParams();
    App.fillSelect('metric', App.RANK_METRICS, params.get('metric') || 'log2fc_cv_median');
    if (params.get('zscore') === 'true') $('zscore').checked = true;
    if (params.get('expressed_only') === 'true') $('expressed').checked = true;
    $('form').addEventListener('submit', (e) => {
      e.preventDefault();
      run();
    });
    $('zscore').addEventListener('change', () => {
      App.setParams({ genes: query().genes, metric: query().metric, expressed_only: query().expressed_only || null, zscore: $('zscore').checked || null });
      render();
    });
    $('example').addEventListener('click', fillExample);
    if (params.get('genes')) {
      $('genes').value = parseGenes(params.get('genes')).join('\n');
      run();
    }
  }

  document.addEventListener('DOMContentLoaded', init);
})();
